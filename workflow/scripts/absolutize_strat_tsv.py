"""Rewrite a GIAB stratification TSV to hold absolute BED paths, for one sample.

The per-sample TSV is written next to, not inside, the extracted <ref>@all
directory, so the shipped relative paths would no longer resolve from it.
Absolute paths make it independent of where it and hap.py run.

PER SAMPLE. The v3.6 set ships GenomeSpecific strata for all seven GIAB
genomes; only the query sample's own are kept. HG002's difficult regions say
nothing about HG007's calls, would read as results about the sample being
scored, and were 27 % of the intervals hap.py's quantify step holds in memory.
In v3.6 every GenomeSpecific BED names exactly one sample (90 files); one that
names none is KEPT and logged rather than dropped, so a future naming scheme
errs visibly.

Fails here, before any hap.py job is submitted, if a referenced BED is missing:
one stratification fewer coming out than going in is not something anyone spots
in a summary table.
"""

import re
from pathlib import Path

from thothlib import die, log_to

sm = snakemake  # noqa: F821
emit = log_to(sm.log[0])
strat_dir = Path(sm.input[0]).resolve()
out_path = Path(sm.output[0])
subset, sample = sm.wildcards.subset, sm.wildcards.sample
_SAMPLE_TOKEN = re.compile(r"(?<![A-Za-z0-9])(HG\d{3})(?![0-9])")

# strat_subset exists because the TSV naming has changed between stratification
# releases; on a miss, show what the tarball holds.
pattern = f"*{subset}-stratifications.tsv"
candidates = sorted(strat_dir.rglob(pattern))
if len(candidates) != 1:
    found = sorted(str(p.relative_to(strat_dir)) for p in strat_dir.rglob("*-stratifications.tsv"))
    die(
        f"{len(candidates)} files match {pattern!r} under {strat_dir}; set `strat_subset` to the part "
        f"before '-stratifications.tsv' of one of: {', '.join(found) or 'none (no stratification TSV in the tarball)'}"
    )
src = candidates[0]
emit(f"{src.relative_to(strat_dir)} -> {out_path}")

rows, missing, dropped, unnamed = [], [], [], []
for lineno, line in enumerate(src.read_text().splitlines(), start=1):
    if not line.strip() or line.startswith("#"):
        continue
    fields = line.split("\t")
    if len(fields) < 2:
        die(f"{src}:{lineno}: expected two tab-separated fields, got {len(fields)}")
    name, rel = fields[0], fields[1].strip()

    # Releases have used paths relative to the TSV and relative to the extracted
    # root; whichever exists wins, and if neither does the row is reported.
    for base in (src.parent, strat_dir):
        resolved = (base / rel).resolve()
        if resolved.exists():
            break
    else:
        missing.append((name, rel))
        continue
    if "GenomeSpecific" in Path(rel).parts:
        named = set(_SAMPLE_TOKEN.findall(Path(rel).name))
        if not named:
            unnamed.append(name)
        elif sample not in named:
            dropped.append(name)
            continue
    rows.append((name, str(resolved)))

if missing:
    listing = "\n".join(f"    {name}\t{rel}" for name, rel in missing[:20])
    more = f"\n    ... and {len(missing) - 20} more" if len(missing) > 20 else ""
    die(
        f"{len(missing)} stratification BED(s) referenced but not present:\n{listing}{more}\n"
        "the extracted tarball is incomplete -- delete it and re-fetch"
    )

out_path.write_text("".join(f"{name}\t{path}\n" for name, path in rows))

kept_gs = sum(1 for _, path in rows if "/GenomeSpecific/" in path)
emit(f"{sample}: kept {kept_gs} GenomeSpecific strata, dropped {len(dropped)} belonging to other samples")
if unnamed:
    emit(f"{len(unnamed)} GenomeSpecific strata name no sample and were KEPT: {unnamed}")
if kept_gs == 0:
    emit(f"note: no GenomeSpecific strata for {sample} in this release")
emit(f"wrote {len(rows)} absolute stratification paths")
