"""Stage the share bundle: metrics, public queries' hap.py VCFs, MANIFEST.txt.

The directory is replaced whole on every run, so a query that stops being
public, or leaves the run file, leaves nothing behind.
"""

import json
import shutil
from pathlib import Path

from thothlib import log_to

sm = snakemake  # noqa: F821
emit = log_to(sm.log[0])
manifest = Path(sm.output.manifest)
out = manifest.parent

shutil.rmtree(out, ignore_errors=True)
out.mkdir(parents=True)
for path in [sm.input.metrics, *sm.input.vcfs]:
    shutil.copy(path, out)
staged = sorted(p.name for p in out.iterdir())

ri = json.loads(Path(sm.input.runinfo).read_text())
git = ri["thoth"]["git"]
code = f"{git['describe']} (commit {git['commit']})" if git["available"] else "not a git checkout"
manifest.write_text(
    "\n".join(
        [
            f"run: {ri['run_id']}",
            f"thoth: {code}, code_sha256 {ri['thoth']['code_sha256']}",
            f"config_sha256: {ri['config_sha256']}",
            *(f"container {n}: {c['image']}@{c['digest']}" for n, c in ri["containers"].items()),
            "rate-level metrics: all queries",
            f"per-variant VCFs withheld (share: internal): {' '.join(sm.params.internal) or '-'}",
            "---",
            *staged,
        ]
    )
    + "\n"
)
emit(f"staged {len(staged)} file(s) in {out}")
