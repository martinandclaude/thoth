"""The provenance record: one JSON per run.

Everything in it is read from what ran -- the code, the container digests, the
checksums of the reference and truth files used, hap.py's own command lines --
rather than restated from documentation. It holds no interpretation (pass/fail,
thresholds): that belongs to whatever consumes it.
"""

import getpass
import hashlib
import importlib
import json
import os
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from thothlib import die, log_to, sha256_file

sm = snakemake  # noqa: F821
emit = log_to(sm.log[0])
repo = Path(sm.params.repo)
cfg = dict(sm.params.config)


def git(*args):
    return subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True, timeout=30
    ).stdout.strip()


def git_state():
    """Commit, branch and dirty flag of the thoth checkout, so a result can be
    traced to the code that made it. Recorded as {"available": false, "reason": ...}
    when thoth is not run from a git clone or git is missing."""
    try:
        return {
            "available": True,
            "commit": git("rev-parse", "HEAD"),
            "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(git("status", "--porcelain")),
            "describe": git("describe", "--always", "--dirty", "--tags"),
        }
    except FileNotFoundError:
        return {"available": False, "reason": "git is not installed"}
    except subprocess.CalledProcessError as exc:
        reason = (exc.stderr or "").strip().splitlines()
        return {"available": False, "reason": reason[0] if reason else "git failed"}
    except subprocess.TimeoutExpired:
        return {"available": False, "reason": "git timed out"}


def code_sha256():
    """Identifies the code with or without git: the sha256 of `sha256sum` output
    over workflow/, config/config.yaml and resources/manifest.tsv, sorted by path.
    From the repository root: find workflow config/config.yaml resources/manifest.tsv
    -type f ! -path '*/__pycache__/*' | LC_ALL=C sort | xargs sha256sum | sha256sum"""
    files = [p for p in (repo / "workflow").rglob("*") if p.is_file() and "__pycache__" not in p.parts]
    rel = sorted(p.relative_to(repo).as_posix() for p in files) + ["config/config.yaml", "resources/manifest.tsv"]
    listing = "".join(f"{sha256_file(repo / r)}  {r}\n" for r in sorted(rel))
    return hashlib.sha256(listing.encode()).hexdigest()


def resources():
    """Per manifest file this run read: its sha256, as fetch_resource recorded it
    (computed here for a file placed by hand), and whether the manifest pins it.
    A pinned file that does not match fails the run."""
    out = {}
    for path, pinned in sm.params.resources.items():
        sidecar = Path(f"{path}.sha256")
        got = sidecar.read_text().strip() if sidecar.exists() else sha256_file(path)
        if pinned and got != pinned:
            die(f"{path} has sha256 {got}, but the manifest pins {pinned}")
        if not pinned:
            emit(f"not pinned in the manifest, recorded as fetched: {path}")
        out[path.removeprefix(f"{sm.params.resources_dir}/")] = {"sha256": got, "pinned": bool(pinned)}
    return out


def happy_commandlines():
    """hap.py keeps argv under runInfo in its sidecar. Its own version string is
    empty in the pinned image, so the digest under `containers` identifies it."""
    out = {}
    for path in sm.input.happy_runinfo:
        info = json.loads(Path(path).read_text())
        out[Path(path).name.removesuffix(".runinfo.json")] = next(
            (r.get("value", "") for r in info.get("runInfo", []) if r.get("key") == "commandline"), ""
        )
    return out


def query_inputs():
    """Per query and caller: size and mtime, and sha256 only when
    runinfo.hash_query_vcfs is set, since hashing WGS VCFs on a shared
    filesystem takes minutes."""
    want_hash = cfg["runinfo"]["hash_query_vcfs"]
    out = {}
    for name, q in sm.params.queries.items():
        entry = {"sample": q["sample"], "reference": q["reference"], "vcf_small": {}, "vcf_sv": {}}
        for key in ("vcf_small", "vcf_sv"):
            for caller, path in q[key].items():
                rec = {"path": path}
                try:
                    stat = os.stat(path)
                except OSError as exc:
                    rec["stat_error"] = str(exc)
                else:
                    rec["bytes"] = stat.st_size
                    rec["mtime_utc"] = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat()
                    if want_hash:
                        rec["sha256"] = sha256_file(path)
                entry[key][caller] = rec
        out[name] = entry
    return out


def prep_counts():
    """Variant counts before and after prep_query's filters (PASS, contigs, ALT
    genotype), per <query>.<caller>. `input` is NA when the input VCF has no index."""
    out = {}
    for path in sm.input.prep_counts:
        lines = Path(path).read_text().splitlines()[1:]
        out[Path(path).name.removesuffix(".prep_counts.txt")] = dict(line.split("\t") for line in lines)
    return out


# The merged config exactly as Snakemake saw it: the validated, defaulted
# queries come in as a separate param, so config_sha256 hashes what was written.
canonical = json.dumps(cfg, sort_keys=True, default=str)
metrics = str(sm.input.metrics)
record = {
    "run_id": sm.params.run_id,
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "thoth": {"repo": str(repo), "git": git_state(), "code_sha256": code_sha256()},
    "environment": {
        "host": socket.gethostname(),
        "user": getpass.getuser(),
        "cwd": os.getcwd(),
        "python": sys.version.split()[0],
        # `snakemake` here is the injected job object, not the module.
        "snakemake": importlib.import_module("snakemake").__version__,
    },
    "config": cfg,
    "config_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
    "containers": sm.params.containers,
    "references": {
        ref: {"fasta": fasta, "fai_sha256": sha256_file(f"{fasta}.fai")}
        for ref, fasta in sm.params.references.items()
    },
    "resources": resources(),
    "queries": query_inputs(),
    "prep_counts": prep_counts(),
    "tools": {"happy": happy_commandlines()},
    "outputs": {"metrics_tsv": metrics, "metrics_sha256": sha256_file(metrics)},
}
with open(sm.output.json, "w") as fh:
    json.dump(record, fh, indent=2, default=str)
    fh.write("\n")
emit(f"provenance record -> {sm.output.json}")
