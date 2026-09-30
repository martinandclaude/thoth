"""Extract the numbers the Pages site shows from a thoth metrics.tsv.

    python docs/extract_evidence.py --query <query> [--caller <small-variant caller>]
        [--runinfo <run>.runinfo.json ...] <run>.metrics.tsv [<earlier run>.metrics.tsv]

Writes docs/data/evidence.json. Only aggregate rates and counts for the public
GIAB sample HG002 leave the metrics table, labelled by public caller names.
The seqsim baseline is copied from docs/data/seqsim_baseline.json (written by
docs/seqsim_baseline.py). Refine results come from the run's own table, or from
the optional second one: an earlier run whose refine results to show. Each
--runinfo, one per metrics file in the same order, adds checksums, the code
commit and container digests, and nothing else: no user, host, path, query or
caller name. Review the JSON before committing.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent / "data"
SAMPLE = "HG002"
CALLER_LABEL = {"deepvariant": "DeepVariant", "sniffles": "Sniffles", "cutesv": "cuteSV"}
TRUTH_LABEL = {"v4.2.1": "GIAB v4.2.1", "v5.0q": "GIAB v5.0q", "CMRG_v1.00": "GIAB CMRG v1.00"}
# index.html looks strata and truth sets up by these labels.
STRATA = [  # (subset in metrics.tsv, label on the site)
    ("*", "All benchmark regions"),
    ("notinalldifficultregions", "Easy: outside every difficult region"),
    ("alldifficultregions", "All difficult regions"),
    ("AllTandemRepeatsandHomopolymers_slop5", "Tandem repeats and homopolymers"),
    ("segdups", "Segmental duplications"),
    ("lowmappabilityall", "Low mappability"),
    ("gclt25orgt65_slop50", "Extreme GC (<25 % or >65 %)"),
    ("MHC", "MHC"),
]


def load(path, query):
    m = pd.read_csv(path, sep="\t", low_memory=False)
    m = m[m["query"] == query]
    if m.empty:
        sys.exit(f"no rows for query {query!r} in {path}")
    return m


def happy_rows(m, truth, subset, caller):
    h = m[(m.tool == "happy") & (m.caller == caller) & (m.truth == truth) & (m["filter"] == "PASS")
          & (m.subtype == "*") & (m.genotype == "*") & (m.subset == subset)]
    out = {}
    for vt in ("SNP", "INDEL"):
        x = h[h.type == vt].set_index("metric")["value"]
        if len(x):
            out[vt] = {"recall": round(float(x["METRIC.Recall"]), 5), "precision": round(float(x["METRIC.Precision"]), 5),
                       "truth_total": int(x["TRUTH.TOTAL"]), "tp": int(x["TRUTH.TP"]), "fn": int(x["TRUTH.FN"]),
                       "fp": int(x["QUERY.FP"])}
    return out


def sv_rows(m, tool):
    s = m[m.tool == tool]
    out = []
    for (caller, truth), g in s.groupby(["caller", "truth"]):
        x = g.set_index("metric")["value"]
        out.append({"caller": CALLER_LABEL.get(caller, caller), "truth": TRUTH_LABEL.get(truth, truth),
                    "recall": round(float(x["recall"]), 5), "precision": round(float(x["precision"]), 5),
                    "base": int(x["base cnt"]), "tp": int(x["TP-base"]), "fn": int(x["FN"]), "fp": int(x["FP"])})
    return out


def small_caller(m, wanted):
    callers = sorted(m.loc[m.tool == "happy", "caller"].unique())
    if wanted in callers or (wanted is None and len(callers) == 1):
        return wanted or callers[0]
    sys.exit(f"pick a small-variant caller with --caller; the table has: {', '.join(callers) or 'none'}")


def provenance(path):
    """A whitelist: runinfo.json also holds user, host, paths and the full config."""
    ri = json.loads(Path(path).read_text())
    git = ri.get("thoth", {}).get("git", {})
    rec = {"metrics_sha256": ri["outputs"]["metrics_sha256"], "config_sha256": ri["config_sha256"]}
    if git.get("available"):
        rec["commit"], rec["describe"] = git["commit"], git["describe"]
    # older runinfo files name the same SIF sha256 `sif_sha256`
    rec["containers"] = {name: c.get("digest") or c.get("sif_sha256") for name, c in ri["containers"].items()}
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", required=True, help="query name in metrics.tsv; read, never written out")
    ap.add_argument("--caller", help="small-variant caller; needed only when the table has more than one")
    ap.add_argument("--runinfo", action="append", default=[], help="<run>.runinfo.json, one per metrics file")
    ap.add_argument("metrics", nargs="+", help="<run>.metrics.tsv, then optionally an earlier run's refine results")
    args = ap.parse_args()
    if len(args.metrics) > 2 or len(args.runinfo) > len(args.metrics):
        ap.error("one or two metrics files, and at most one --runinfo per metrics file")
    baseline = DATA / "seqsim_baseline.json"
    if not baseline.exists():
        sys.exit(f"{baseline} missing: run docs/seqsim_baseline.py (see its docstring)")

    m = load(args.metrics[0], args.query)
    caller = small_caller(m, args.caller)
    label = CALLER_LABEL.get(caller, caller)
    ev = {
        "sample": SAMPLE,
        "technology": "Oxford Nanopore",
        "truth_sets": [{"truth": TRUTH_LABEL[t], "caller": label, **happy_rows(m, t, "*", caller)}
                       for t in ("v4.2.1", "v5.0q", "CMRG_v1.00")],
        "strata": [{"stratum": lab, "caller": label, "truth": TRUTH_LABEL["v4.2.1"], **happy_rows(m, "v4.2.1", s, caller)}
                   for s, lab in STRATA],
        "sv_bench": sv_rows(m, "truvari"),
    }
    missing = [r.get("stratum", r["truth"]) for r in ev["truth_sets"] + ev["strata"] if not {"SNP", "INDEL"} <= r.keys()]
    if missing:
        sys.exit(f"no SNP or INDEL rows for {caller} in: {', '.join(missing)}")
    src = load(args.metrics[1], args.query) if len(args.metrics) > 1 else m
    if (src.tool == "truvari_refine").any():
        ev["sv_refine_bench"] = sv_rows(src, "truvari")
        ev["sv_refine"] = sv_rows(src, "truvari_refine")
    ev["seqsim_baseline"] = json.loads(baseline.read_text())
    if args.runinfo:
        ev["provenance"] = [provenance(p) for p in args.runinfo]
    out = DATA / "evidence.json"
    out.write_text(json.dumps(ev, indent=1) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
