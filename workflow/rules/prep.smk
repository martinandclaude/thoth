rule check_contigs:
    """Fail if a query VCF declares a contig the reference lacks: the cheapest
    guard against scoring an hg19 callset against hg38 truth and getting a
    plausible, meaningless number. Run per VCF, both variant classes."""
    input:
        vcf=lambda w: query_vcf(w.query, w.caller, w.kind),
        fai=lambda w: ref_fasta(w.query) + ".fai",
    output:
        f"{OUT}/{{query}}/{{caller}}/prep/{{kind}}.contigs.ok",
    log:
        f"{OUT}/{{query}}/{{caller}}/prep/{{kind}}.contigs.log",
    container:
        container_for("bcftools")
    shell:
        # No pipe feeds a decision: under Snakemake's pipefail, `! comm | grep -q .`
        # passes a large mismatch (grep exits early, comm gets SIGPIPE, 141 is
        # negated). Every stage logs a line: an empty log means no shell started.
        r"""
        exec 2>> {log}
        echo "[thoth] check_contigs {wildcards.query}/{wildcards.caller} ({wildcards.kind}) on $(hostname) at $(date -u +%FT%TZ)" >&2
        echo "[thoth]   query: {input.vcf}" >&2
        echo "[thoth]   fai:   {input.fai}" >&2

        bcftools view -h {input.vcf} | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' | sort -u > {output}.query
        echo "[thoth]   query header declares $(wc -l < {output}.query) contigs" >&2

        cut -f1 {input.fai} | sort -u > {output}.ref
        echo "[thoth]   reference .fai has $(wc -l < {output}.ref) contigs" >&2

        if [ ! -s {output}.query ]; then
            echo "[thoth] FAIL: no ##contig lines in the query header, so it cannot be checked" >&2
            exit 1
        fi

        comm -23 {output}.query {output}.ref > {output}.missing
        if [ -s {output}.missing ]; then
            echo "[thoth] FAIL: $(wc -l < {output}.missing) query contig(s) absent from the reference. First 10:" >&2
            head -n 10 {output}.missing >&2
            exit 1
        fi

        mv {output}.query {output}
        rm -f {output}.ref {output}.missing
        echo "[thoth] ok" >&2
        """


rule prep_query:
    """PASS-only, the configured regions, and genotypes that carry an ALT allele.

    `-i 'GT~"[1-9]"'` keeps every genotype with a non-reference allele index,
    half-calls such as ./1 included; 0/0, 0/. and ./. are dropped. Both
    `-i 'GT="alt"'` and the common `-e 'GT="ref" || GT="mis"'` treat a genotype
    with any missing allele as missing and drop ./1, although the ALT is
    asserted present. Joint callsets (GLnexus) are full of half-calls, and
    dropping them inflates FN.
    """
    input:
        vcf=lambda w: query_vcf(w.query, w.caller, "small"),
        ok=f"{OUT}/{{query}}/{{caller}}/prep/small.contigs.ok",
    output:
        vcf=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.prepped.vcf.gz",
        tbi=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.prepped.vcf.gz.tbi",
        counts=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.prep_counts.txt",
    params:
        pass_only=lambda w: "-f PASS,." if config["small_variants"]["pass_only"] else "",
        # -t (targets), NOT -r (regions): -r needs an index on the INPUT, and
        # callers hand over plain .vcf as often as indexed .vcf.gz. -t streams
        # the same contig filter with no index at all.
        regions=lambda w: (f"-t {regions_arg(w.query)}" if regions_arg(w.query) else ""),
    log:
        f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.prep.log",
    container:
        container_for("bcftools")
    threads: 2
    resources:
        mem_mb=8000,
        runtime=60,
    shell:
        r"""
        bcftools view {params.pass_only} {params.regions} -i 'GT~"[1-9]"' \
            --threads {threads} -Oz -o {output.vcf} {input.vcf} 2> {log}
        bcftools index -t {output.vcf} 2>> {log}
        {{
          printf 'stage\tvariants\n'
          printf 'input\t%s\n'   "$(bcftools index -n {input.vcf} 2>/dev/null || echo NA)"
          printf 'prepped\t%s\n' "$(bcftools index -n {output.vcf})"
        }} > {output.counts}
        """
