# The SV arm. Runs for every (query, caller) under a query's `vcf_sv`, against
# every truth set the manifest has an sv_vcf row for -- which, for GIAB, means
# HG002 only (v5.0q genome-wide, CMRG v1.00 in 273 genes). SV_ARMS in common.smk.
#
# SYMBOLIC ALLELES, read from the truvari 5.4.0 source: a symbolic <DEL>, <DUP>
# or <INV> is resolved against --reference (matching.py:186), but a symbolic
# <INS> cannot be -- its sequence is unknown. When either side is unresolved the
# pctseq check is SKIPPED, not failed (variant_record.py:822), so such a call is
# matched on size and position alone. Two callers where one writes resolved
# insertions and the other symbolic ones are therefore NOT held to the same
# standard in the same table.

# truvari 5.4.0 bench flags that config keys may name, long form, so a config
# key reads as the flag it sets.
_TRUVARI_BENCH_FLAGS = {
    "refdist", "pctseq", "pctsize", "pctovl", "sizemin", "sizefilt", "sizemax", "chunksize",
    "bnddist", "pick", "typeignore", "dup-to-ins", "no-ref", "max-resolve", "refine",
}


def truvari_settings(truth):
    """Global truvari settings with this benchmark's own `truvari:` block on top.

    Per benchmark because the right comparison depends on the benchmark: GIAB's
    v5.0q README recommends `--refine --pick ac -r 2000 -C 5000` because, unlike
    earlier releases, v5.0q does not exclude complex SVs, so one event can be
    written several valid ways. Settings that suit v5.0q need not suit CMRG.
    """
    merged = dict(config["structural_variants"]["truvari"])
    extra = [str(merged.pop("extra", "")).strip()]
    override = dict(BENCH[truth].get("truvari") or {})
    extra.append(str(override.pop("extra", "")).strip())
    merged.update(override)
    merged["extra"] = " ".join(e for e in extra if e)
    return merged


def truvari_args(truth):
    """truvari_settings(truth) as a flag string.

    Unknown keys fail at parse time (every truth is rendered once below) rather
    than as a truvari usage error in job 40 of 40. `extra` is passed through
    verbatim and is the escape hatch for anything not in the set above.
    """
    cfg = truvari_settings(truth)
    extra = cfg.pop("extra", "")
    parts = []
    for key, value in cfg.items():
        flag = key.replace("_", "-")
        if flag == "refine":  # run as its own command: see truvari_refine
            continue
        if flag not in _TRUVARI_BENCH_FLAGS:
            _fail(
                f"truvari setting '{key}' (for truth {truth}) is not a truvari "
                f"bench flag. Known: {sorted(_TRUVARI_BENCH_FLAGS)}. Put raw "
                f"flags in `extra`."
            )
        # Booleans are argparse store_true switches: present or absent, never
        # given a value. `--dup-to-ins True` is rejected by truvari's parser.
        if isinstance(value, bool):
            if value:
                parts.append(f"--{flag}")
            continue
        parts.append(f"--{flag} {value}")
    if extra:
        parts.append(extra)
    return " ".join(parts)


TRUVARI_ARGS = {t: truvari_args(t) for t in BENCH}
REFINE_TRUTHS = sorted(t for t in BENCH if truvari_settings(t).get("refine"))
_NO_REFINE_TRUTHS = sorted(t for t in BENCH if t not in REFINE_TRUTHS)


rule prep_sv_truth:
    """Drop ALT=* records from an SV truth VCF.

    GIAB's v5.0q README: "for truvari bench, variants with ALT=* are often
    incorrectly categorized, so users are encouraged to filter these variants
    from the benchmark VCF prior to benchmarking". They are spanning-deletion
    placeholders from decomposed multiallelics and carry no sequence, so
    refine's haplotype reconstruction gains nothing from them either. Applied to
    every SV truth set -- a no-op where there are none. Lives next to the truth
    file: derived once per truth set, shared by every query and caller.
    """
    input:
        RES + "/{dest}",
    output:
        vcf=RES + "/{dest}.noAltStar.vcf.gz",
        tbi=RES + "/{dest}.noAltStar.vcf.gz.tbi",
    log:
        RES + "/{dest}.noAltStar.log",
    container:
        container_for("bcftools")
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        r"""
        bcftools view -e 'ALT="*"' -Oz -o {output.vcf} {input} 2> {log}
        bcftools index -t {output.vcf} 2>> {log}
        echo "[thoth] kept $(bcftools index -n {output.vcf}) of $(bcftools view -H {input} | wc -l) records" >> {log}
        """


rule prep_sv_query:
    """Sort, recompress and index the caller's SV VCF into results/.

    Sorted because SV callers do not reliably emit sorted output and bcftools
    index refuses unsorted input; truvari needs an index on both sides. No
    filtering: `--passonly` in the truvari config does that, and doing it twice
    would make the params.json truvari writes disagree with what it was given.
    """
    input:
        vcf=lambda w: query_vcf(w.query, w.caller, "sv"),
        ok=f"{OUT}/{{query}}/{{caller}}/prep/sv.contigs.ok",
    output:
        vcf=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.sv.vcf.gz",
        tbi=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.sv.vcf.gz.tbi",
    log:
        f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.sv_prep.log",
    container:
        container_for("bcftools")
    threads: 2
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        r"""
        # -T is a mkdtemp TEMPLATE, not a directory: given "/tmp", bcftools
        # makes "/tmpXXXXXX" at the filesystem root. Same shape as its default.
        bcftools sort -T "${{TMPDIR:-/tmp}}/bcftools-sort.XXXXXX" -Oz -o {output.vcf} {input.vcf} 2> {log}
        bcftools index -t {output.vcf} 2>> {log}
        """


# THE OUTPUT DIRECTORY, from the truvari 5.4.0 source: bench exits if its -o
# directory already exists (bench.py:131) and creates it with os.mkdir
# (bench.py:321), so the parent must exist and the directory itself must not.
# Snakemake creates the parent of every declared output, which is the directory
# truvari refuses. So bench writes to a work directory that then replaces the
# output directory whole: merging into it would leave stale files from other
# settings, and mv cannot replace the phab_bench/ that refine leaves behind. The
# log lives beside that directory, not in it, because the directory is deleted.
_TRUVARI_DIR = f"{OUT}/{{query}}/{{caller}}/{{truth}}/truvari"


rule truvari:
    """truvari bench, one (query, caller, truth) arm, for truth sets without refine."""
    input:
        query=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.sv.vcf.gz",
        query_tbi=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.sv.vcf.gz.tbi",
        truth=lambda w: sv_truth_vcf(w.query, w.truth) + ".noAltStar.vcf.gz",
        truth_tbi=lambda w: sv_truth_vcf(w.query, w.truth) + ".noAltStar.vcf.gz.tbi",
        truth_bed=lambda w: sv_truth_bed(w.query, w.truth),
        ref=lambda w: ref_fasta(w.query),
    output:
        summary=f"{_TRUVARI_DIR}/summary.json",
        params=f"{_TRUVARI_DIR}/params.json",
    wildcard_constraints:
        truth=_alternation(_NO_REFINE_TRUTHS),
    params:
        outdir=lambda w, output: str(Path(output.summary).parent),
        workdir=lambda w, output: str(Path(output.summary).parent) + "_work",
        args=lambda w: TRUVARI_ARGS[w.truth],
        refine=lambda w: "true" if w.truth in REFINE_TRUTHS else "false",
    log:
        f"{OUT}/{{query}}/{{caller}}/{{truth}}/{{query}}.{{caller}}.{{truth}}.truvari.log",
    container:
        container_for("truvari")
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240,
    shell:
        r"""
        rm -rf {params.workdir}
        truvari bench \
            --base {input.truth} \
            --comp {input.query} \
            --output {params.workdir} \
            --reference {input.ref} \
            --includebed {input.truth_bed} \
            {params.args} \
            > {log} 2>&1
        if {params.refine}; then
            truvari refine --threads {threads} {params.workdir} >> {log} 2>&1
        fi
        rm -rf {params.outdir}
        mv {params.workdir} {params.outdir}
        """


# A second rule because Snakemake cannot make an output conditional: refine's
# summary is declared only for truth sets whose settings ask for --refine
# (v5.0q). Bench's own summary.json is kept too: the unrefined and refined
# numbers side by side show how much of the miss rate was representation rather
# than calling. GIAB's README runs `bench --refine`, which calls `truvari refine`
# with its defaults (bench.py:803): the POA aligner and 4 worker processes,
# whatever the job was given. Running refine as its own command is the same
# computation with the job's threads. The thread count cannot change a number:
# each region is aligned on its own and the output is sorted before the final
# bench (utils.compress_index_vcf); refine.variant_summary.json is byte-identical
# at 1, 4 and 8 threads on the test data.
#
# Measured on HG002 against v5.0q, genome-wide ONT calls (Sniffles, cuteSV):
# bench plus refine took 4.5 min on 8 CPUs and peaked at 59 and 77 GiB (sacct
# MaxRSS), almost all of it a few regions with large insertions, where abPOA's
# matrix grows with haplotype length. When a refine worker is killed for
# memory, truvari's pool waits forever for its region: the job idles until the
# time limit instead of failing. A refine job running far past minutes has hit
# the memory limit; raise mem_mb.
use rule truvari as truvari_refine with:
    output:
        summary=f"{_TRUVARI_DIR}/summary.json",
        params=f"{_TRUVARI_DIR}/params.json",
        refined=f"{_TRUVARI_DIR}/refine.variant_summary.json",
    wildcard_constraints:
        truth=_alternation(REFINE_TRUTHS),
    threads: 8
    resources:
        mem_mb=128000,
        runtime=240,
