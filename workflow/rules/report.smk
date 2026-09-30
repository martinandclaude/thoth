def _arm_path(q, c, t, tail):
    return f"{OUT}/{q}/{c}/{t}/{tail}"


def _happy_extended(wildcards=None):
    require_arms()
    return [_arm_path(q, c, t, f"happy/{q}.{c}.{t}.extended.csv") for q, c, t in SMVAR_ARMS]


rule tidy_metrics:
    """One row per (query, caller, truth, stratum, vartype, metric): long format,
    so a new arm adds rows, not columns."""
    input:
        happy=_happy_extended,
        truvari=[_arm_path(q, c, t, "truvari/summary.json") for q, c, t in SV_ARMS],
        truvari_refine=[
            _arm_path(q, c, t, "truvari/refine.variant_summary.json")
            for q, c, t in SV_ARMS
            if t in REFINE_TRUTHS
        ],
    output:
        tsv=f"{OUT}/report/{RUN_ID}.metrics.tsv",
    params:
        run_id=RUN_ID,
        queries=QUERIES,
        bench=BENCH,
        # From the flags each truth set actually ran with, global and its own.
        sv_filter={t: "PASS" if "--passonly" in TRUVARI_ARGS[t].split() else "ALL" for t in BENCH},
    log:
        f"{OUT}/report/{RUN_ID}.tidy.log",
    script:
        "../scripts/tidy_metrics.py"


rule runinfo:
    """Provenance record: code version, config hash, container digests, the
    checksums of the reference and truth files used, hap.py's command lines."""
    input:
        happy_runinfo=[_arm_path(q, c, t, f"happy/{q}.{c}.{t}.runinfo.json") for q, c, t in SMVAR_ARMS],
        prep_counts=[
            f"{OUT}/{q}/{c}/prep/{q}.{c}.prep_counts.txt"
            for q, c in sorted({(q, c) for q, c, _ in SMVAR_ARMS})
        ],
        fai=sorted({ref_fasta(q) + ".fai" for q in QUERIES}),
        metrics=f"{OUT}/report/{RUN_ID}.metrics.tsv",
    output:
        json=f"{OUT}/report/{RUN_ID}.runinfo.json",
    params:
        run_id=RUN_ID,
        repo=str(REPO),
        config=dict(config),
        queries=QUERIES,
        references={r: REFERENCES[r]["fasta"] for r in sorted({q_ref(q) for q in QUERIES})},
        containers={
            r["name"]: {k: r[k] for k in ("image", "digest", "built_from", "recipe")}
            for _, r in CONTAINERS.iterrows()
        },
        resources_dir=RES,
        resources=USED_RESOURCES,
    log:
        f"{OUT}/report/{RUN_ID}.runinfo.log",
    script:
        "../scripts/runinfo.py"


rule share_bundle:
    """The share gate: what may leave this machine, in results/share/<run_id>/.

    Rate-level metrics for every query. The annotated hap.py VCF, which holds
    every small-variant TP, FP and FN call, only for queries marked
    `share: public`; SV per-variant files are never staged. runinfo.json stays
    behind (it records user, host and every query's paths); MANIFEST.txt
    carries the code version, config hash and container digests instead.
    """
    input:
        metrics=f"{OUT}/report/{RUN_ID}.metrics.tsv",
        runinfo=f"{OUT}/report/{RUN_ID}.runinfo.json",
        vcfs=[
            _arm_path(q, c, t, f"happy/{q}.{c}.{t}.vcf.gz")
            for q, c, t in SMVAR_ARMS
            if QUERIES[q]["share"] == "public"
        ],
    output:
        manifest=f"{OUT}/share/{RUN_ID}/MANIFEST.txt",
    params:
        internal=sorted(q for q in QUERIES if QUERIES[q]["share"] == "internal"),
    log:
        f"{OUT}/report/{RUN_ID}.share.log",
    script:
        "../scripts/share_bundle.py"
