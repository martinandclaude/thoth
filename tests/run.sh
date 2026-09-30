#!/usr/bin/env bash
# End-to-end test: the real pipeline, on planted data, in the real containers.
#
#   THOTH_CONTAINERS_DIR=tests/work/sifs workflow/scripts/fetch_containers.sh
#   tests/run.sh
#
# Needs apptainer and the driver environment (pip install --require-hashes -r
# requirements-driver.lock). Runs in CI on every push to main and every pull
# request (.github/workflows/tests.yml), and anywhere else the same way.
set -euo pipefail
cd "$(dirname "$0")/.."
W=tests/work

# Start from nothing but the fetched SIFs: a previous run's results/ and
# resources/ would take jobs out of the DAG and make the job-count check lie.
rm -rf -- "$W/results" "$W/resources" "$W/giab" "$W/calls"
python3 tests/make_fixtures.py

# The fixture truth VCFs must be BGZF with tabix indexes, as GIAB's are. Done in
# the pipeline's own hts image so the host needs no compression tools.
HTS=$(python3 -c 'import sys; sys.path.insert(0, "workflow/scripts"); import thothlib as t
r = t.containers().set_index("name").loc["bcftools"]; print("tests/work/sifs/" + t.sif_filename(r.image, r.digest))')
apptainer exec "$HTS" sh -c '
  set -e
  for v in tests/work/giab/truth/*.vcf; do bgzip -f "$v"; tabix -f -p vcf "$v.gz"; done'

echo "=== dry run: expected job counts ==="
snakemake --configfile "$W/test.yaml" -n > "$W/dry.txt" 2>&1 || { cat "$W/dry.txt"; exit 1; }
python3 - "$W/dry.txt" <<'PY'
import re, sys
text = open(sys.argv[1]).read()
counts = dict(re.findall(r"^(\w+)\s+(\d+)$", text, re.M))
want = {"happy": "2", "truvari": "1", "truvari_refine": "1", "prep_sv_truth": "2",
        "check_contigs": "3", "absolutize_strat_tsv": "1"}
bad = {k: (counts.get(k), v) for k, v in want.items() if counts.get(k) != v}
print("  job counts", "ok" if not bad else f"WRONG {bad}")
note = "[thoth] e2e-HG003/svcaller: skipped for TESTSV, TESTSV2 (no SV truth set for HG003 on GRCh38)"
print(f"  {'ok  ' if note in text else 'FAIL'} SV caller of a sample without SV truth skipped with a note")
sys.exit(1 if bad or note not in text else 0)
PY

echo "=== parse-time validation ==="
python3 - "$W" <<'PY'
import subprocess, sys, yaml
W = sys.argv[1]
base = yaml.safe_load(open(f"{W}/test.yaml"))
Q = {"name": "q", "sample": "HG002", "reference": "GRCh38"}
# label: (run-file changes, None removes a key; target; expected message, "" = must succeed)
cases = {
    "misspelt query key": ({"queries": [{**Q, "vcf_svs": {"c": "x.vcf"}}]}, "all", "unknown key(s) ['vcf_svs']"),
    "caller name with a dot": ({"queries": [{**Q, "vcf_small": {"deep.variant": "x.vcf"}}]}, "all",
                               "caller name 'deep.variant'"),
    "list instead of mapping": ({"queries": [{**Q, "vcf_small": ["x.vcf"]}]}, "all",
                                "must map caller names to VCF paths"),
    "sample not in the manifest": ({"queries": [{**Q, "sample": "NA24385", "vcf_small": {"c": "x.vcf"}}]}, "all",
                                   "sample 'NA24385' is not in"),
    "misspelt top-level key": ({"runid": "x"}, "all", "unknown key(s) ['runid'] in the config"),
    "misspelt benchmark key": ({"benchmarks": [{"truth": "TESTSMALL", "strat": "vT", "strat_subet": "all"}]},
                               "all", "unknown key(s) ['strat_subet'] in benchmark TESTSMALL"),
    "run_id that is a path": ({"run_id": "../x"}, "all", "run_id '../x' must be a plain name"),
    "no queries": ({"queries": []}, "all", "nothing to benchmark"),
    "no references, with queries": ({"references": None}, "all", "reference 'GRCh38' has no entry under `references`"),
    "no references, check_manifest": ({"references": None, "queries": []}, "check_manifest", ""),
}
bad = 0
for label, (change, target, expect) in cases.items():
    cfg = {k: v for k, v in {**base, **change}.items() if v is not None}
    path = f"{W}/bad.yaml"
    yaml.safe_dump(cfg, open(path, "w"))
    out = subprocess.run(["snakemake", target, "--configfile", path, "-n"], capture_output=True, text=True)
    ok = out.returncode == 0 if not expect else out.returncode != 0 and expect in out.stdout + out.stderr
    print(f"  {'ok  ' if ok else 'FAIL'} {label}")
    bad += not ok
sys.exit(1 if bad else 0)
PY

echo "=== full run ==="
snakemake --configfile "$W/test.yaml" --cores 2 --software-deployment-method apptainer --printshellcmds

echo "=== a fetched file that differs from its pinned sha256 ==="
python3 - "$W" <<'PY'
import subprocess, sys, yaml
W = sys.argv[1]
rows = [r.rsplit("\t", 1)[0] + "\t" + "0" * 64 if "\ttruth/small.bed\t" in r else r
        for r in open(f"{W}/manifest.tsv").read().splitlines()]
open(f"{W}/bad_manifest.tsv", "w").write("\n".join(rows) + "\n")
yaml.safe_dump(yaml.safe_load(open(f"{W}/test.yaml")) | {"manifest": f"{W}/bad_manifest.tsv"}, open(f"{W}/bad.yaml", "w"))
out = subprocess.run(["snakemake", "--configfile", f"{W}/bad.yaml", "-n"], capture_output=True, text=True)
ok = out.returncode != 0 and "truth/small.bed has sha256" in out.stdout + out.stderr
print(f"  {'ok  ' if ok else 'FAIL'} fails at parse time")
sys.exit(0 if ok else 1)
PY

echo "=== results ==="
python3 tests/check_results.py
