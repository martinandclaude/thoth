"""Assert the end-to-end run reproduced exactly what make_fixtures.py planted.

Every expectation here is a fact about the fixtures, derived in their
docstring -- not a tolerance. A rebuilt container, a dependency bump or a rule
change that alters any count fails this, which is the point.
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

W = Path(__file__).resolve().parent / "work"
R = W / "results"
failures = []


def check(label, got, want):
    ok = got == want
    print(f"  {'ok  ' if ok else 'FAIL'} {label}: {got!r}" + ("" if ok else f" (expected {want!r})"))
    if not ok:
        failures.append(label)


m = pd.read_csv(R / "report" / "e2e.metrics.tsv", sep="\t", low_memory=False)


def val(**sel):
    q = m
    for k, v in sel.items():
        q = q[q[k] == v]
    if len(q) != 1:
        return f"<{len(q)} rows>"
    v = q["value"].iloc[0]
    return int(v) if float(v).is_integer() else float(v)


print("small variants, hap.py (TESTSMALL, PASS, all regions)")
hs = dict(tool="happy", truth="TESTSMALL", caller="snvcaller", subset="*", subtype="*", filter="PASS", genotype="*")
for vt, tp, fn, fp in (("SNP", 9, 1, 2), ("INDEL", 2, 1, 0)):
    check(f"{vt} TRUTH.TP", val(**hs, type=vt, metric="TRUTH.TP"), tp)
    check(f"{vt} TRUTH.FN", val(**hs, type=vt, metric="TRUTH.FN"), fn)
    check(f"{vt} QUERY.FP", val(**hs, type=vt, metric="QUERY.FP"), fp)

print("stratification is per sample")
subsets = set(m[m.tool == "happy"].subset.dropna())
check("own GenomeSpecific stratum present", "HG002_v4.2.1_testregion" in subsets, True)
check("other sample's stratum absent", "HG007_v4.2.1_testregion" not in subsets, True)
check("SNP TP in stratum 'left'", val(**{**hs, "subset": "left"}, type="SNP", metric="TRUTH.TP"), 9)

print("an IUPAC code in the reference (B at 7004; truth and query write N)")
per_variant = R / "e2e-HG002" / "snvcaller" / "TESTSMALL" / "happy" / "e2e-HG002.snvcaller.TESTSMALL.vcf.gz"
rows = [line.rstrip("\n").split("\t") for line in gzip.open(per_variant, "rt") if not line.startswith("#")]


def bd(row, sample):
    return dict(zip(row[8].split(":"), row[sample].split(":")))["BD"]


check("deletion over it scored TP in truth and query (REF, truth, query)",
      [(r[3], bd(r, 9), bd(r, 10)) for r in rows if r[1] == "7000"], [("CTTTN", "TP", "TP")])

print("SVs with refine (TESTSV: one deletion written two ways)")
sv = dict(truth="TESTSV", caller="svcaller")
check("bench TP-base", val(**sv, tool="truvari", metric="TP-base"), 0)
check("bench FP", val(**sv, tool="truvari", metric="FP"), 2)
check("bench FN", val(**sv, tool="truvari", metric="FN"), 1)
check("refined TP-base", val(**sv, tool="truvari_refine", metric="TP-base"), 1)
check("refined TP-comp", val(**sv, tool="truvari_refine", metric="TP-comp"), 2)
check("refined FP", val(**sv, tool="truvari_refine", metric="FP"), 0)
check("refined FN", val(**sv, tool="truvari_refine", metric="FN"), 0)

print("SVs without refine (TESTSV2: symbolic DUP vs resolved INS, 30 kb symbolic DEL)")
sv2 = dict(truth="TESTSV2", caller="svcaller", tool="truvari")
check("TP-base", val(**sv2, metric="TP-base"), 2)
check("FP", val(**sv2, metric="FP"), 0)
check("FN", val(**sv2, metric="FN"), 0)
check("no refine rows for TESTSV2", len(m[(m.truth == "TESTSV2") & (m.tool == "truvari_refine")]), 0)
tp = R / "e2e-HG002" / "svcaller" / "TESTSV2" / "truvari" / "tp-comp.vcf.gz"
seqsim = sorted(f.split("=")[1] for line in gzip.open(tp, "rt") if not line.startswith("#")
                for f in line.split("\t")[7].split(";") if f.startswith("PctSeqSimilarity="))
check("both matched ON SEQUENCE (dup-to-ins, max-resolve)", seqsim, ["1", "1"])
check("SV rows labelled with the filter truvari ran with", sorted(set(m[m.type == "SV"]["filter"])), ["PASS"])

print("query preparation")
prepped = R / "e2e-HG002" / "snvcaller" / "prep" / "e2e-HG002.snvcaller.prepped.vcf.gz"
kept = [line.split("\t")[1] for line in gzip.open(prepped, "rt") if not line.startswith("#")]
check("./1 half-call kept", "59000" in kept, True)

print("truth preparation")
log = next((W / "resources" / "truth").glob("sv_refine.vcf.gz.noAltStar.log")).read_text()
check("ALT=* dropped from SV truth", "kept 1 of 2 records" in log, True)

print("provenance and sharing")
ri = json.loads((R / "report" / "e2e.runinfo.json").read_text())
check("every image recorded by digest", all(v["digest"].startswith("sha256:") for v in ri["containers"].values()), True)
check("callers recorded per query", sorted(ri["queries"]["e2e-HG002"]["vcf_sv"]), ["svcaller"])
check("config recorded as written, no defaults added", "share" not in ri["config"]["queries"][2], True)
check("reference identified by its .fai", ri["references"]["GRCh38"]["fai_sha256"],
      hashlib.sha256((W / "ref" / "chr1.fa.fai").read_bytes()).hexdigest())
check("truth and stratification files used, pinned or not",
      {k: v["pinned"] for k, v in ri["resources"].items()},
      {"truth/small.vcf.gz": False, "truth/small.bed": True, "truth/sv_refine.vcf.gz": False,
       "truth/sv_refine.bed": True, "truth/sv_plain.vcf.gz": False, "truth/sv_plain.bed": True,
       "strat/GRCh38/vT/genome-stratifications-GRCh38-vT.tar.gz": True})
check("each with its sha256", all(len(v["sha256"]) == 64 for v in ri["resources"].values()), True)
check("hap.py command line recorded", "--engine=vcfeval" in ri["tools"]["happy"]["e2e-HG002.snvcaller.TESTSMALL"], True)
shared = sorted(p.name for p in (R / "share" / "e2e").iterdir())
manifest = (R / "share" / "e2e" / "MANIFEST.txt").read_text()
check("public query's per-variant VCF shared", "e2e-HG002.snvcaller.TESTSMALL.vcf.gz" in shared, True)
check("internal query's files withheld", [f for f in shared if f.startswith("e2e-internal.")], [])
check("internal queries listed as withheld", "(share: internal): e2e-HG003 e2e-internal\n" in manifest, True)
check("runinfo.json not shared", "e2e.runinfo.json" not in shared, True)
check("MANIFEST carries the code hash", f"code_sha256 {ri['thoth']['code_sha256']}\n" in manifest, True)
check("MANIFEST carries the config hash", f"config_sha256: {ri['config_sha256']}\n" in manifest, True)
check("MANIFEST carries every container digest",
      all(f"{c['image']}@{c['digest']}" in manifest for c in ri["containers"].values()), True)

print()
if failures:
    print(f"{len(failures)} check(s) FAILED: {failures}")
    sys.exit(1)
print("all checks passed")
