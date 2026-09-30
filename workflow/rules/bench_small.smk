# hap.py with --engine=vcfeval: vcfeval's haplotype-aware matching, and hap.py's
# counting per GIAB stratification in the same pass.


rule happy:
    input:
        query=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.prepped.vcf.gz",
        query_tbi=f"{OUT}/{{query}}/{{caller}}/prep/{{query}}.{{caller}}.prepped.vcf.gz.tbi",
        truth=lambda w: truth_vcf(w.query, w.truth),
        truth_tbi=lambda w: truth_vcf(w.query, w.truth) + ".tbi",
        truth_bed=lambda w: truth_bed(w.query, w.truth),
        strat=lambda w: strat_tsv(w.query, w.truth),
        ref=lambda w: happy_ref(w.query),
        sdf=lambda w: ref_sdf(w.query),
    output:
        summary=f"{OUT}/{{query}}/{{caller}}/{{truth}}/happy/{{query}}.{{caller}}.{{truth}}.summary.csv",
        extended=f"{OUT}/{{query}}/{{caller}}/{{truth}}/happy/{{query}}.{{caller}}.{{truth}}.extended.csv",
        roc=f"{OUT}/{{query}}/{{caller}}/{{truth}}/happy/{{query}}.{{caller}}.{{truth}}.roc.all.csv.gz",
        vcf=f"{OUT}/{{query}}/{{caller}}/{{truth}}/happy/{{query}}.{{caller}}.{{truth}}.vcf.gz",
        runinfo=f"{OUT}/{{query}}/{{caller}}/{{truth}}/happy/{{query}}.{{caller}}.{{truth}}.runinfo.json",
    params:
        prefix=lambda w: f"{OUT}/{w.query}/{w.caller}/{w.truth}/happy/{w.query}.{w.caller}.{w.truth}",
        engine=config["small_variants"]["engine"],
        score_field=lambda w: QUERIES[w.query]["score_field"],
        # Redundant with prep_query, deliberately: hap.py's own PASS handling is
        # what ends up in the runinfo.json a reviewer will read.
        pass_only=lambda w: "--pass-only" if config["small_variants"]["pass_only"] else "",
        locations=lambda w: (f"-l {regions_arg(w.query)}" if regions_arg(w.query) else ""),
        extra=config["small_variants"]["happy_extra"],
        # hap.py's bundled rtg.cfg sizes vcfeval's JVM at 90 % of the NODE's RAM
        # unless RTG_MEM is set; an exported RTG_MEM reaches the vcfeval
        # subprocess. Half the request; the rest is for hap.py's quantify.
        rtg_mem=lambda w, resources: f"{int(resources.mem_mb * 0.5)}m",
    log:
        f"{OUT}/{{query}}/{{caller}}/{{truth}}/happy/{{query}}.{{caller}}.{{truth}}.log",
    container:
        container_for("happy")
    threads: 16
    resources:
        # Set by the stratification set's size, not by the rule: see config.yaml.
        mem_mb=config["small_variants"]["happy_mem_mb"],
        runtime=480,
    shell:
        r"""
        # Full path: the image's PATH does not include /opt/hap.py/bin.
        export RTG_MEM={params.rtg_mem}
        /opt/hap.py/bin/hap.py \
            {input.truth} {input.query} \
            -f {input.truth_bed} \
            -r {input.ref} \
            -o {params.prefix} \
            --engine={params.engine} \
            --engine-vcfeval-template {input.sdf} \
            --stratification {input.strat} \
            --roc {params.score_field} \
            {params.pass_only} {params.locations} {params.extra} \
            --threads {threads} \
            > {log} 2>&1
        """
