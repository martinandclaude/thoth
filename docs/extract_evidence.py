"""Extract the numbers the Pages site shows: ONT's public GIAB 2025.01 calls scored by thoth.

    python docs/extract_evidence.py [--runinfo <comparison>.runinfo.json] \
        --repro <reproduction dir> <comparison>.metrics.tsv

Writes docs/data/evidence.json with two keys. "ont": ONT's published numbers
(data/ont_giab_2025.01_published.tsv) beside the same public calls scored by
thoth (the metrics table) and ONT's SV recipe rerun (<reproduction dir>/bench/,
see data/README.md). "seqsim_baseline": copied from data/seqsim_baseline.json
(written by docs/seqsim_baseline.py). Only aggregate rates and counts for the
public GIAB samples leave the metrics table. --runinfo adds checksums and
container digests through a whitelist, and nothing else: no user, host, path,
query or caller name. Review the JSON before committing.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent / "data"
TRUTH_LABEL = {"v5.0q": "GIAB v5.0q", "V0.018": "GIAB draft V0.018"}


def happy_counts(m, sample, caller, truth, vtype):
    h = m[(m.tool == "happy") & (m["sample"] == sample) & (m.caller == caller) & (m.truth == truth)
          & (m["filter"] == "PASS") & (m.subtype == "*") & (m.genotype == "*") & (m.subset == "*")
          & (m.type == vtype)].set_index("metric")["value"]
    if h.empty:
        sys.exit(f"no hap.py {vtype} row for {sample} {caller} against {truth}")
    return {"recall": round(float(h["METRIC.Recall"]), 6), "precision": round(float(h["METRIC.Precision"]), 6),
            "truth_total": int(h["TRUTH.TOTAL"]), "tp": int(h["TRUTH.TP"]), "fn": int(h["TRUTH.FN"]),
            "fp": int(h["QUERY.FP"])}


def sv_counts(x):
    """A truvari summary (summary.json keys, or metrics.tsv metric -> value)."""
    return {"recall": round(float(x["recall"]), 6), "precision": round(float(x["precision"]), 6),
            "base": int(x["base cnt"]), "tp": int(x["TP-base"]), "fn": int(x["FN"]), "fp": int(x["FP"])}


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


def ont_section(metrics, repro):
    # rates to six digits: published and thoth's differ in the fifth
    pub = pd.read_csv(DATA / "ont_giab_2025.01_published.tsv", sep="\t", comment="#")
    m = pd.read_csv(metrics, sep="\t", low_memory=False)
    small = [{"sample": r.sample, "flowcell": r.flowcell, "type": r.type,
              "thoth": happy_counts(m, r.sample, "ont-" + r.flowcell.lower(), r.truth, r.type),
              "published": {"recall": round(r.recall, 6), "precision": round(r.precision, 6),
                            "truth_total": int(r.truth_total), "tp": int(r.tp), "fn": int(r.fn), "fp": int(r.fp)}}
             for r in pub[pub.tool == "happy"].itertuples()]
    b = Path(repro) / "bench"
    runs = {"bench": json.loads((b / "summary.json").read_text()),
            "refine": json.loads((b / "refine.variant_summary.json").read_text())}
    sv = []
    for r in pub[pub.type == "SV"].itertuples():
        method = "refine" if r.tool == "truvari_refine" else "bench"
        row = {"base cnt": r.truth_total, "TP-base": r.tp, "FN": r.fn, "FP": r.fp,
               "recall": r.recall, "precision": r.precision}
        for who, x in (("ONT, as published", row), ("ONT's recipe, rerun", runs[method])):
            sv.append({"scored_by": who, "truth": TRUTH_LABEL[r.truth], "method": method, **sv_counts(x)})
        caller = "ont-sniffles-" + r.flowcell.lower()   # the flow cell ONT published
        for truth in ("V0.018", "v5.0q"):
            x = m[(m.tool == r.tool) & (m.caller == caller) & (m.truth == truth)].set_index("metric")["value"]
            if x.empty:
                sys.exit(f"no {r.tool} rows for {caller} against {truth} in {metrics}")
            sv.append({"scored_by": "thoth", "truth": TRUTH_LABEL[truth], "method": method, **sv_counts(x)})
    # truvari 4.3.1 refine --recount: phab_bench's harmonised variants inside the
    # refined regions plus bench's own counts outside them.
    phab = json.loads((b / "phab_bench" / "summary.json").read_text())
    outside = runs["refine"]["base cnt"] - phab["base cnt"]
    split = {"regions_refined": int(pd.read_csv(b / "refine.regions.txt", sep="\t").refined.sum()),
             "truth_outside": outside, "bench_truth_inside": runs["bench"]["base cnt"] - outside,
             "refine_truth_inside": phab["base cnt"]}
    return {"small": small, "sv": sv, "recount_split": split}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repro", required=True, help="directory holding bench/ from ONT's SV recipe rerun")
    ap.add_argument("--runinfo", help="runinfo.json of the comparison run")
    ap.add_argument("metrics", help="metrics.tsv of ONT's public GIAB 2025.01 calls scored by thoth")
    args = ap.parse_args()
    baseline = DATA / "seqsim_baseline.json"
    if not baseline.exists():
        sys.exit(f"{baseline} missing: run docs/seqsim_baseline.py (see its docstring)")
    ont = ont_section(args.metrics, args.repro)
    if args.runinfo:
        ont["provenance"] = provenance(args.runinfo)
    out = DATA / "evidence.json"
    out.write_text(json.dumps({"ont": ont, "seqsim_baseline": json.loads(baseline.read_text())}, indent=1) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
