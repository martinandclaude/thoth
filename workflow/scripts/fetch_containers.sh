#!/usr/bin/env bash
# Download thoth's SIFs from GHCR, by digest, and verify them.
#
#   workflow/scripts/fetch_containers.sh           # fetch missing, verify present
#   workflow/scripts/fetch_containers.sh --test    # ...then run each test_cmd
#   DRY_RUN=1 workflow/scripts/fetch_containers.sh # say what would happen
#   THOTH_CONTAINERS_DIR=dir workflow/scripts/fetch_containers.sh
#       # test suite only; must equal the containers_dir Snakemake sees
#
# Run from the repository root, on a node with outbound HTTPS, in the
# environment you run snakemake from (python3 with PyYAML and pandas).
#
# WHY DOWNLOAD RATHER THAN PULL. The images are built and converted to SIF by
# this repo's containers workflow on GitHub, then published as ORAS artifacts.
# A SIF stored that way is a single registry blob whose digest is the file's
# own sha256, so fetching it is one anonymous HTTPS GET -- no apptainer, no
# conversion, no APPTAINER_TMPDIR -- and the digest recorded in
# workflow/containers.tsv is at once the pin and the integrity check. A file
# whose sha256 does not match is never moved into place.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case "${1:-}" in
    "") DO_TEST=0 ;;
    --test) DO_TEST=1 ;;
    *) echo "usage: $0 [--test]" >&2; exit 2 ;;
esac

exec python3 - "$REPO" "$DO_TEST" "${DRY_RUN:-0}" <<'PY'
import json, os, shutil, subprocess, sys, urllib.request
import yaml

repo, do_test, dry = sys.argv[1], sys.argv[2] == "1", sys.argv[3] not in ("", "0")
sys.path.insert(0, os.path.join(repo, "workflow", "scripts"))
from thothlib import containers, sha256_file, sif_filename

# Same resolution as common.smk: site.yaml over config.yaml; containers_dir,
# else <resources_dir>/containers.
cfg = yaml.safe_load(open(os.path.join(repo, "config", "config.yaml"))) or {}
site_path = os.path.join(repo, "config", "site.yaml")
if os.path.exists(site_path):
    cfg.update(yaml.safe_load(open(site_path)) or {})
cdir = os.environ.get("THOTH_CONTAINERS_DIR") or cfg.get("containers_dir")
if not cdir:
    if not cfg.get("resources_dir"):
        sys.exit("[thoth] no containers_dir or resources_dir: fill in config/site.yaml "
                 "(from config/site.example.yaml)")
    cdir = str(cfg["resources_dir"]).rstrip("/") + "/containers"
cdir = cdir.rstrip("/")

rows = containers()


def sha256(path):
    return "sha256:" + sha256_file(path)


def download(image, digest, dest):
    registry, name = image.split("/", 1)
    tok = json.load(urllib.request.urlopen(
        f"https://{registry}/token?scope=repository:{name}:pull", timeout=60))["token"]
    part = dest + ".part"
    # curl rather than urllib: the blob redirects to object storage, and curl
    # drops the Authorization header on a cross-host redirect, as it should.
    subprocess.run(["curl", "-fL", "--retry", "3", "--retry-delay", "10",
                    "-H", f"Authorization: Bearer {tok}",
                    "-o", part, f"https://{registry}/v2/{name}/blobs/{digest}"], check=True)
    got = sha256(part)
    if got != digest:
        os.remove(part)
        raise RuntimeError(f"downloaded {got}, expected {digest}")
    os.replace(part, dest)


print(f"[thoth] containers dir: {cdir}")
if not dry:
    os.makedirs(cdir, exist_ok=True)

bad = False
status = {}
for (image, digest), group in rows.groupby(["image", "digest"], sort=False):
    names = ", ".join(group["name"])
    if not digest.startswith("sha256:"):
        print(f"  {names:18s} {digest:9s} not yet published -- run the containers workflow, then record its digest")
        bad = True
        continue
    path = os.path.join(cdir, sif_filename(image, digest))
    if os.path.exists(path):
        if dry:
            print(f"  {names:18s} present  {path}")
            continue
        got = sha256(path)
        if got != digest:
            print(f"  {names:18s} SHA256 MISMATCH  expected {digest[:23]}..  file {got[:23]}..")
            bad = True
            continue
        print(f"  {names:18s} ok       {digest[:23]}..")
    else:
        if dry:
            print(f"  {names:18s} would download {image}@{digest[:23]}..")
            continue
        print(f"  {names:18s} downloading {image}@{digest[:23]}..")
        try:
            download(image, digest, path)
        except Exception as exc:  # noqa: BLE001
            print(f"  {names:18s} FAILED: {exc}")
            bad = True
            continue
        print(f"  {names:18s} ok       {digest[:23]}..  {os.path.getsize(path) / 1e6:.0f} MB")
    for n in group["name"]:
        status[n] = path

if do_test and not dry:
    runner = shutil.which("apptainer") or shutil.which("singularity")
    if not runner:
        sys.exit("[thoth] --test needs apptainer on PATH")
    print("[thoth] test_cmd inside each SIF")
    for _, r in rows.iterrows():
        if r["name"] not in status:
            continue
        ok = subprocess.run([runner, "exec", status[r["name"]], "sh", "-c", r["test_cmd"]],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
        print(f"  {r['name']:18s} {'OK' if ok else 'FAILED: ' + r['test_cmd']}")
        bad |= not ok

sys.exit(1 if bad else 0)
PY
