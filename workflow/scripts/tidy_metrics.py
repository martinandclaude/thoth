"""Collapse every hap.py and truvari result into one long table: one row per
(query, caller, truth, tool, type, subtype, subset, filter, metric). A new
stratification, truth set or arm adds rows, never columns."""

import json
from pathlib import Path

import pandas as pd
from thothlib import die, log_to

sm = snakemake  # noqa: F821
emit = log_to(sm.log[0])

# The hap.py columns worth scoring on. One that a hap.py version does not report
# is logged and skipped rather than filled with NA, so it shows as missing rows.
HAPPY_METRICS = [
    "METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score", "METRIC.Frac_NA",
    "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK",
    "FP.gt", "FP.al", "Subset.Size", "Subset.IS_CONF.Size",
]
ID_COLUMNS = ["Type", "Subtype", "Subset", "Filter", "Genotype", "QQ.Field"]
COLUMN_ORDER = [
    "run_id", "query", "caller", "sample", "reference", "origin", "share", "score_field",
    "truth", "strat", "strat_subset", "tool", "type", "subtype", "subset", "filter",
    "genotype", "qq_field", "metric", "value",
]


def arm(path):
    """(query, caller, truth) from results/<query>/<caller>/<truth>/<tool>/<file>."""
    return Path(path).parts[-5:-2]


def labels(query, caller, truth, tool):
    q = sm.params.queries[query]
    fields = ("sample", "reference", "origin", "share", "score_field")
    return {"run_id": sm.params.run_id, "query": query, "caller": caller, "truth": truth,
            "tool": tool, **{k: q[k] for k in fields}}


rows = []
for path in sm.input.happy:
    query, caller, truth = arm(path)
    table = pd.read_csv(path)
    absent = [c for c in HAPPY_METRICS if c not in table.columns]
    if absent:
        emit(f"{query}/{caller} x {truth}: hap.py did not report {absent}")
    ids = [c for c in ID_COLUMNS if c in table.columns]
    long = table.melt(id_vars=ids, value_vars=[c for c in HAPPY_METRICS if c not in absent],
                      var_name="metric", value_name="value")
    long = long[long["value"].notna()].rename(columns={c: c.lower().replace(".", "_") for c in ids})
    bench = sm.params.bench[truth]
    rows.append(long.assign(**labels(query, caller, truth, "happy"),
                            strat=bench["strat"], strat_subset=bench["strat_subset"]))
    emit(f"{query}/{caller} x {truth}: {len(long)} rows from {Path(path).name}")

# `truvari` is bench as-is, `truvari_refine` is after refine re-aligned the
# candidate regions. Both are kept: the gap between them is how much of the miss
# rate was representation, not calling.
for tool, paths in (("truvari", sm.input.truvari), ("truvari_refine", sm.input.truvari_refine)):
    for path in paths:
        query, caller, truth = arm(path)
        summary = json.loads(Path(path).read_text())
        scalars = {k: v for k, v in summary.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}
        skipped = [k for k in summary if k not in scalars]
        if skipped:
            emit(f"{query}/{caller} x {truth} ({tool}): non-scalar keys stay in the json: {skipped}")
        base = {**labels(query, caller, truth, tool), "type": "SV", "subtype": "*", "subset": "*",
                "filter": sm.params.sv_filter[truth], "genotype": "*"}
        if scalars:
            rows.append(pd.DataFrame([{**base, "metric": k, "value": v} for k, v in scalars.items()]))
        emit(f"{query}/{caller} x {truth}: {len(scalars)} rows from {tool} {Path(path).name}")

if not rows:
    die("no results to tidy -- every input list was empty")
metrics = pd.concat(rows, ignore_index=True).reindex(columns=COLUMN_ORDER).sort_values(
    ["query", "caller", "truth", "tool", "type", "subtype", "subset", "filter", "metric"], kind="stable"
)
# %.10g rather than the default repr: `value` holds both rates and counts, so it
# is one float64 column, and the default would write every count as 918.0. Ten
# significant digits is lossless for a proportion and a genome-sized base count.
metrics.to_csv(sm.output.tsv, sep="\t", index=False, float_format="%.10g")
emit(
    f"{len(metrics)} rows, {metrics['query'].nunique()} queries, {metrics['caller'].nunique()} callers, "
    f"{metrics['truth'].nunique()} truth sets -> {sm.output.tsv}"
)
