"""How similar unrelated random sequences look to truvari's seqsim (site section 4).

    apptainer exec <sif dir>/truvari-<digest>.sif python docs/seqsim_baseline.py

Runs inside the truvari image from workflow/containers.tsv, since truvari is
not in the driver environment. Deterministic (seed 2026). Writes
docs/data/seqsim_baseline.json, which extract_evidence.py merges into
docs/data/evidence.json.
"""
import json
import random
import statistics
from importlib.metadata import version
from pathlib import Path

from truvari.comparisons import seqsim  # the function bench uses for PctSeqSimilarity

PAIRS = 200
THRESHOLD = 0.7  # truvari's default pctseq
rng = random.Random(2026)


def random_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def quantile(v, p):
    return v[min(len(v) - 1, int(p * len(v)))]


rows = []
for n in (50, 100, 200, 300, 500, 1000, 2000, 5000):
    v = sorted(seqsim(random_seq(n), random_seq(n)) for _ in range(PAIRS))
    rows.append({"length": n, "pairs": len(v), "min": round(v[0], 4), "p05": round(quantile(v, 0.05), 4),
                 "median": round(statistics.median(v), 4), "p95": round(quantile(v, 0.95), 4),
                 "max": round(v[-1], 4),
                 "share_at_or_above_0_7": round(sum(x >= THRESHOLD for x in v) / len(v), 4)})

baseline = {
    "method": f"truvari.comparisons.seqsim on pairs of unrelated uniform-random sequences, {PAIRS} pairs per length, "
              f"truvari {version('truvari')} (docs/seqsim_baseline.py)",
    "threshold": THRESHOLD,
    "rows": rows,
}
out = Path(__file__).resolve().parent / "data" / "seqsim_baseline.json"
out.write_text(json.dumps(baseline, indent=1) + "\n")
print(f"wrote {out}")
