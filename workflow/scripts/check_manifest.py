"""HEAD every URL in the resource manifest and write a report.

Run this before the first fetch: it turns a four-hour download that dies on the
last file into a ten-second table.

Deliberately always exits 0. This is a report, not a gate: a non-zero exit would
make Snakemake delete the very table you asked for. Read the `ok` column.
"""

import os
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from thothlib import log_to

TIMEOUT = 30
WORKERS = 8
UA = "thoth-check-manifest/1.0"

sm = snakemake  # noqa: F821
emit = log_to(sm.log[0])


def _check_file_url(url):
    """file:// is supported: for a site that mirrors GIAB locally, or has no
    outbound network, the manifest can point at local copies -- and the test
    suite serves its fixtures this way."""
    path = urllib.parse.unquote(urllib.parse.urlparse(url).path)
    if not os.path.exists(path):
        return False, "MISSING", "", "no such file"
    if not os.access(path, os.R_OK):
        return False, "DENIED", "", "exists but is not readable"
    return True, "LOCAL", str(os.path.getsize(path)), ""


def _check_http_url(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
            size = resp.headers.get("Content-Length", "")
            return True, str(resp.status), size, ""
    except urllib.error.HTTPError as exc:
        # Some servers answer HEAD with 403 or 405 while serving GET fine.
        # A one-byte ranged GET distinguishes "not there" from "no HEAD here".
        if exc.code in (403, 405, 501):
            try:
                req = urllib.request.Request(
                    url, headers={"User-Agent": UA, "Range": "bytes=0-0"}
                )
                with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                    size = resp.headers.get("Content-Range", "").rsplit("/", 1)[-1]
                    return True, f"{resp.status} (GET)", size, f"HEAD gave {exc.code}"
            except Exception as inner:  # noqa: BLE001
                return False, str(exc.code), "", f"HEAD and ranged GET both failed: {inner}"
        return False, str(exc.code), "", exc.reason or ""
    except Exception as exc:  # noqa: BLE001
        return False, "ERROR", "", f"{type(exc).__name__}: {exc}"


def check(url):
    return _check_file_url(url) if url.startswith("file://") else _check_http_url(url)


targets = sm.params.targets
emit(f"checking {len(targets)} URLs")
with ThreadPoolExecutor(max_workers=WORKERS) as pool:
    results = list(pool.map(lambda t: check(t[2]), targets))

report = pd.DataFrame(
    [(kind, dest, "yes" if ok else "NO", status, size, note, url)
     for (kind, dest, url), (ok, status, size, note) in zip(targets, results)],
    columns=["kind", "dest", "ok", "http_status", "size_bytes", "note", "url"],
)
report.to_csv(sm.output[0], sep="\t", index=False)

bad = report[report["ok"] == "NO"]
total_gb = pd.to_numeric(report["size_bytes"], errors="coerce").sum() / 1e9
emit(f"{len(report) - len(bad)}/{len(report)} reachable, {total_gb:.1f} GB total, report at {sm.output[0]}")
if len(bad):
    emit("unreachable:\n" + "\n".join(f"    {s:>10}  {d}" for s, d in zip(bad["http_status"], bad["dest"])))
    emit("fix these in the manifest before running resources_all. Only http(s):// and file:// URLs can be fetched.")
