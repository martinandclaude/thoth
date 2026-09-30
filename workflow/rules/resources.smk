# Resources are immutable, shared across projects, and live outside results/.
# Nothing here depends on a query wildcard — that separation is deliberate.


localrules:
    check_manifest,


# A login node without outbound network sets fetch_on_login_node: false in
# config/site.yaml; the fetch rules then run as cluster jobs.
if config.get("fetch_on_login_node", True):

    localrules:
        fetch_resource,
        fetch_index,


rule check_manifest:
    """HEAD every URL. Run this once before the first fetch — it turns a
    four-hour download failure into a ten-second report."""
    output:
        f"{OUT}/report/manifest_check.tsv",
    params:
        targets=fetched(MANIFEST),
    log:
        f"{OUT}/report/manifest_check.log",
    script:
        "../scripts/check_manifest.py"


rule fetch_resource:
    output:
        protected(RES + "/{dest}"),
    params:
        url=lambda w: DEST_URL[w.dest],
    log:
        RES + "/logs/{dest}.fetch.log",
    retries: 3
    shell:
        r"""
        mkdir -p $(dirname {output}) $(dirname {log})
        curl -fsSL --retry 3 --retry-delay 10 -o {output}.part '{params.url}' 2> {log}
        mv {output}.part {output}
        sha256sum {output} | awk '{{print $1}}' > {output}.sha256
        """


rule fetch_index:
    """GIAB ships .tbi next to the VCF; fetching beats re-indexing."""
    output:
        protected(RES + "/{dest}.tbi"),
    params:
        url=lambda w: DEST_URL[w.dest] + ".tbi",
    log:
        RES + "/logs/{dest}.tbi.fetch.log",
    retries: 3
    shell:
        r"""
        mkdir -p $(dirname {output}) $(dirname {log})
        curl -fsSL --retry 3 --retry-delay 10 -o {output}.part '{params.url}' 2> {log}
        mv {output}.part {output}
        """


rule extract_stratifications:
    input:
        RES + "/strat/{reference}/{version}/genome-stratifications-{reference}-{version}.tar.gz",
    output:
        directory(RES + "/strat/{reference}/{version}/{reference}@all"),
    log:
        RES + "/logs/extract_strat_{reference}_{version}.log",
    resources:
        runtime=120,
    shell:
        r"""
        # --strip-components=1 rather than extracting into the parent: the
        # tarball's own top-level directory is CHM13@all where this pipeline
        # calls the reference CHM13v2.0, and upstream naming has already been
        # wrong once (see the note in resources/manifest.tsv). Stripping it
        # makes the extracted path entirely ours. Safe because these tarballs
        # have exactly one top-level entry.
        mkdir -p {output}
        tar -xzf {input} -C {output} --strip-components=1 > {log} 2>&1
        """


rule absolutize_strat_tsv:
    """The stratification TSV with absolute BED paths, because the per-sample
    copy is written outside the extracted tree its relative paths point into.
    One per sample: only the query sample's own GenomeSpecific strata are kept
    (see the script)."""
    input:
        RES + "/strat/{reference}/{version}/{reference}@all",
    output:
        RES + "/strat/{reference}/{version}/{subset}.{sample}.abs.tsv",
    log:
        RES + "/logs/strat_{reference}_{version}_{subset}_{sample}.log",
    script:
        "../scripts/absolutize_strat_tsv.py"


REF_FASTAS = [r["fasta"] for r in REFERENCES.values() if r.get("fasta")]


rule faidx_reference:
    """A .fai beside every configured reference FASTA, for check_contigs. The
    wildcard is pinned to the configured FASTAs: an unconstrained {fasta}.fai
    rule would match every path in the workflow."""
    input:
        "{fasta}",
    output:
        "{fasta}.fai",
    log:
        "{fasta}.faidx.log",
    wildcard_constraints:
        fasta=_alternation(REF_FASTAS),
    container:
        container_for("samtools")
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        "samtools faidx {input} > {log} 2>&1"


rule rtg_format:
    """SDF for the vcfeval engine, one per reference. Built by the rtg bundled
    in the hap.py image, the same one vcfeval reads it with; it is not on the
    image PATH."""
    input:
        lambda w: REFERENCES[w.reference]["fasta"],
    output:
        directory(RES + "/sdf/{reference}.sdf"),
    params:
        # Unset, the rtg launcher sets -Xmx to 90 % of the HOST's RAM (it reads
        # the machine, not the cgroup) and the job is OOM-killed on a shared
        # node. 80 % of the request leaves room for the JVM's non-heap memory.
        rtg_mem=lambda w, resources: f"{int(resources.mem_mb * 0.8)}m",
    log:
        RES + "/logs/rtg_format_{reference}.log",
    container:
        container_for("happy")
    resources:
        mem_mb=16000,
        runtime=120,
    shell:
        "RTG_MEM={params.rtg_mem} /opt/hap.py/libexec/rtg-tools-install/rtg format -o {output} {input} > {log} 2>&1"


rule happy_reference:
    """The reference with its IUPAC ambiguity codes turned into N, for hap.py.

    GRCh38 carries 94 ambiguity bases (B, K, M, R, S, W, Y); GIAB's truth sets
    and callers write N there. hap.py rewrites the query's REF from its -r FASTA
    before vcfeval, so with the codes left in the query says TGTGB where the
    truth says TGTGN and vcfeval stops: "disagree on what the reference bases
    should be". Every other rule reads the original FASTA; the SDF does not
    matter (checked both ways).
    """
    input:
        lambda w: REFERENCES[w.reference]["fasta"],
    output:
        fasta=RES + "/reference/{reference}.iupacN.fa",
        fai=RES + "/reference/{reference}.iupacN.fa.fai",
    log:
        RES + "/logs/happy_reference_{reference}.log",
    container:
        container_for("samtools")
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        r"""
        sed '/^>/!s/[^ACGTNacgtn]/N/g' {input} > {output.fasta}.part 2> {log}
        mv {output.fasta}.part {output.fasta}
        samtools faidx {output.fasta} 2>> {log}
        """
