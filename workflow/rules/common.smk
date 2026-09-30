import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
import thothlib  # noqa: E402

REPO = Path(workflow.basedir).parent


def _fail(msg):
    raise WorkflowError(f"[thoth] {msg}")


def _alternation(names):
    """A regex matching exactly these names; '^$' when there are none."""
    return "|".join(re.escape(n) for n in names) if names else "^$"


def _check_keys(where, block, allowed):
    unknown = sorted(set(block) - set(allowed))
    if unknown:
        _fail(f"unknown key(s) {unknown} in {where}. Allowed: {sorted(allowed)}")


# --------------------------------------------------------------------------
# Config. A misspelt key that was silently ignored would run with the default
# instead, so unknown keys fail here.
# --------------------------------------------------------------------------
_check_keys(
    "the config",
    config,
    {
        "results_dir", "resources_dir", "containers_dir", "fetch_on_login_node",
        "manifest", "references", "benchmarks", "small_variants",
        "structural_variants", "runinfo", "run_id", "queries",
    },
)
_check_keys(
    "small_variants",
    config["small_variants"],
    {"engine", "pass_only", "regions", "default_score_field", "happy_extra", "happy_mem_mb"},
)
_check_keys("structural_variants", config["structural_variants"], {"truvari"})
_check_keys("runinfo", config["runinfo"], {"hash_query_vcfs"})

if not str(config.get("resources_dir") or "").strip():
    _fail("no resources_dir configured. Copy config/site.example.yaml to config/site.yaml and fill it in.")
RES = str(config["resources_dir"]).rstrip("/")
OUT = str(config["results_dir"]).rstrip("/")
RUN_ID = str(config["run_id"])
REFERENCES = config.get("references") or {}
# run_id becomes a directory that share_bundle deletes and recreates.
if not re.match(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", RUN_ID):
    _fail(f"run_id '{RUN_ID}' must be a plain name: letters, digits, '.', '_' and '-'")

MANIFEST_PATH = REPO / config["manifest"]
MANIFEST = pd.read_csv(MANIFEST_PATH, sep="\t", comment="#", dtype=str).fillna("")
DEST_URL = dict(zip(MANIFEST["dest"], MANIFEST["url"]))
DEST_SHA256 = dict(zip(MANIFEST["dest"], MANIFEST["sha256"]))
SAMPLES = sorted(set(MANIFEST["sample"]) - {"-"})
# Stratifications are fetched only for references that have truth rows.
USABLE = MANIFEST[MANIFEST["reference"].isin(MANIFEST.loc[MANIFEST["kind"] != "strat_tar", "reference"])]
STRAT_ROWS = USABLE[USABLE["kind"] == "strat_tar"]


def fetched(rows):
    """(kind, dest, url) of every file these manifest rows fetch: the row itself,
    and <dest>.tbi where index == yes."""
    out = []
    for kind, dest, url, index in zip(rows["kind"], rows["dest"], rows["url"], rows["index"]):
        out.append((kind, dest, url))
        if index == "yes":
            out.append((f"{kind}.tbi", f"{dest}.tbi", f"{url}.tbi"))
    return out


BENCH = {}
for b in config["benchmarks"]:
    _check_keys(f"benchmark {b.get('truth')}", b, {"truth", "strat", "strat_subset", "truvari"})
    if not b.get("truth") or not b.get("strat"):
        _fail(f"benchmark {b}: `truth` and `strat` are required")
    if b["truth"] in BENCH:
        _fail(f"duplicate benchmark truth: {b['truth']}")
    BENCH[b["truth"]] = {"strat_subset": "all", **b}


# --------------------------------------------------------------------------
# Containers: prebuilt SIFs on the shared filesystem, one per tool.
# --------------------------------------------------------------------------
CONTAINERS = thothlib.containers()
_IMAGES = {r["name"]: (r["image"], r["digest"]) for _, r in CONTAINERS.iterrows()}
CONTAINERS_DIR = str(config.get("containers_dir") or f"{RES}/containers").rstrip("/")


def container_for(name):
    """Local SIF path, not a registry reference, so Snakemake never pulls or
    converts an image. Not checked here (parse time); fetch_containers.sh
    reports missing SIFs."""
    if name not in _IMAGES:
        _fail(f"no container named '{name}' in {thothlib.CONTAINERS_TSV}")
    return f"{CONTAINERS_DIR}/{thothlib.sif_filename(*_IMAGES[name])}"


# --------------------------------------------------------------------------
# Queries: one sample on one reference, callers under vcf_small / vcf_sv (see
# runs/run.example.yaml). The unit of work is the arm (query, caller, truth).
# --------------------------------------------------------------------------
REQUIRED_KEYS = ("name", "sample", "reference")
OPTIONAL_KEYS = ("vcf_small", "vcf_sv", "origin", "share", "score_field", "notes")
VALID_SHARE = ("public", "internal")
# Names are path components and wildcard values: no dots, which would make
# `{query}.{caller}.{truth}.extended.csv` ambiguous, and no slashes.
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def _callers(qname, entry, key):
    block = entry.get(key) or {}
    if not isinstance(block, dict):
        _fail(
            f"query {qname}: `{key}` must map caller names to VCF paths, e.g.\n"
            f"    {key}:\n      deepvariant: /path/to/calls.vcf.gz\n"
            f"  got {type(block).__name__}: {block!r}"
        )
    for caller, path in block.items():
        if not _NAME_RE.match(str(caller)):
            _fail(
                f"query {qname}: caller name '{caller}' under {key} must match "
                f"{_NAME_RE.pattern} (it becomes a path component)"
            )
        if not isinstance(path, str) or not path:
            _fail(f"query {qname}: {key}.{caller} must be a VCF path, got {path!r}")
    return {str(c): str(v) for c, v in block.items()}


QUERIES = {}
for raw in config["queries"] or []:
    # A copy: defaults written back into config would change config_sha256.
    entry = dict(raw)
    missing = [k for k in REQUIRED_KEYS if not entry.get(k)]
    if missing:
        _fail(f"query entry is missing required key(s) {missing}: {entry}")
    name = entry["name"]
    _check_keys(f"query {name}", entry, REQUIRED_KEYS + OPTIONAL_KEYS)
    if not isinstance(name, str) or not _NAME_RE.match(name):
        _fail(f"query name '{name}' must match {_NAME_RE.pattern} (it becomes a path component)")
    if name in QUERIES:
        _fail(f"duplicate query name: {name}")
    if entry["reference"] not in REFERENCES:
        _fail(
            f"query {name}: reference '{entry['reference']}' has no entry under "
            f"`references` in config/site.yaml (have {sorted(REFERENCES)})"
        )
    if not REFERENCES[entry["reference"]].get("fasta"):
        _fail(f"query {name}: no FASTA configured for {entry['reference']}")
    if entry["sample"] not in SAMPLES:
        _fail(
            f"query {name}: sample '{entry['sample']}' is not in {MANIFEST_PATH}. "
            f"Use the GIAB ID, one of {SAMPLES}"
        )
    entry.setdefault("share", "internal")
    if entry["share"] not in VALID_SHARE:
        _fail(f"query {name}: share must be one of {VALID_SHARE}")
    entry.setdefault("origin", "unspecified")
    entry.setdefault("score_field", config["small_variants"]["default_score_field"])
    entry["vcf_small"] = _callers(name, entry, "vcf_small")
    entry["vcf_sv"] = _callers(name, entry, "vcf_sv")
    if not entry["vcf_small"] and not entry["vcf_sv"]:
        _fail(f"query {name}: needs at least one of vcf_small / vcf_sv")
    QUERIES[name] = entry

if QUERIES and RUN_ID == "unnamed-run":
    _fail("set run_id in the run file; it names the report and share directory")

CALLERS = sorted({c for q in QUERIES.values() for c in (*q["vcf_small"], *q["vcf_sv"])})


# --------------------------------------------------------------------------
# Manifest lookup
# --------------------------------------------------------------------------
def _row(sample, reference, release, kind):
    hit = MANIFEST[
        (MANIFEST["sample"] == sample)
        & (MANIFEST["reference"] == reference)
        & (MANIFEST["release"] == release)
        & (MANIFEST["kind"] == kind)
    ]
    return None if hit.empty else hit.iloc[0]


def _dest(sample, reference, release, kind):
    row = _row(sample, reference, release, kind)
    if row is None:
        _fail(f"no manifest row for sample={sample} reference={reference} release={release} kind={kind}")
    return f"{RES}/{row['dest']}"


# --------------------------------------------------------------------------
# The arms: each (query, caller) VCF against every truth set the manifest
# offers, in its variant class, for the query's sample and reference.
# --------------------------------------------------------------------------
def _build_arms(key, kind, label):
    """Arms of one variant class, noting per (query, caller) the truth sets that
    offer this class but not for the query's sample and reference."""
    offered = set(MANIFEST.loc[MANIFEST["kind"] == kind, "release"])
    arms = []
    for name, q in QUERIES.items():
        for caller in q[key]:
            hits = [t for t in BENCH if _row(q["sample"], q["reference"], t, kind) is not None]
            arms += [(name, caller, t) for t in hits]
            missed = [t for t in BENCH if t in offered and t not in hits]
            if missed:
                print(
                    f"[thoth] {name}/{caller}: skipped for {', '.join(missed)} "
                    f"(no {label} truth set for {q['sample']} on {q['reference']})",
                    file=sys.stderr,
                )
    return arms


SMVAR_ARMS = _build_arms("vcf_small", "smvar_vcf", "small-variant")
SV_ARMS = _build_arms("vcf_sv", "sv_vcf", "SV")
# Queries that resolve to no arm are a config error, decidable now (only the
# manifest is read). An EMPTY query list is not: the resource targets need no
# run file, and require_arms() guards the benchmark targets.
if QUERIES and not SMVAR_ARMS and not SV_ARMS:
    _fail(
        "no (query, caller, truth) arm resolved from "
        f"{len(QUERIES)} quer{'y' if len(QUERIES) == 1 else 'ies'}. Every one was "
        "skipped for want of a matching manifest row -- see the notes above."
    )


def require_arms(wildcards=None):
    """An input function, so it runs only when a benchmark target is requested."""
    if not SMVAR_ARMS and not SV_ARMS:
        _fail(
            "nothing to benchmark -- no queries. Pass a run file (targets go BEFORE --configfile):\n"
            "    snakemake --profile profiles/slurm --configfile runs/<run>.yaml\n"
            "  The resource targets (check_manifest, resources_all) need no run file."
        )
    return []


# --------------------------------------------------------------------------
# Accessors used by the benchmark rules
# --------------------------------------------------------------------------
def q_ref(query):
    return QUERIES[query]["reference"]


def q_sample(query):
    return QUERIES[query]["sample"]


def ref_fasta(query):
    return REFERENCES[q_ref(query)]["fasta"]


def ref_sdf(query):
    # Keyed by the reference name alone: a name must mean the same FASTA
    # everywhere this resources_dir is shared.
    return f"{RES}/sdf/{q_ref(query)}.sdf"


def happy_ref(query):
    """The reference hap.py reads: ambiguity codes as N (rule happy_reference).
    Keyed like ref_sdf."""
    return f"{RES}/reference/{q_ref(query)}.iupacN.fa"


def query_vcf(query, caller, kind):
    """The caller's VCF for this query. kind is `small` or `sv`."""
    return QUERIES[query][f"vcf_{kind}"][caller]


def truth_vcf(query, truth):
    return _dest(q_sample(query), q_ref(query), truth, "smvar_vcf")


def truth_bed(query, truth):
    return _dest(q_sample(query), q_ref(query), truth, "smvar_bed")


def sv_truth_vcf(query, truth):
    return _dest(q_sample(query), q_ref(query), truth, "sv_vcf")


def sv_truth_bed(query, truth):
    return _dest(q_sample(query), q_ref(query), truth, "sv_bed")


def strat_tar(query, truth):
    return _dest("-", q_ref(query), BENCH[truth]["strat"], "strat_tar")


def strat_tsv(query, truth):
    """Per sample: only the query sample's own GenomeSpecific strata are kept
    (see absolutize_strat_tsv.py)."""
    b = BENCH[truth]
    return f"{RES}/strat/{q_ref(query)}/{b['strat']}/{b['strat_subset']}.{q_sample(query)}.abs.tsv"


def regions_arg(query):
    """Contig list for hap.py -l. GRCh37 has no chr prefix."""
    spec = config["small_variants"]["regions"]
    if spec == "all":
        return ""
    if spec != "autosomes":
        return spec
    prefix = "" if q_ref(query).startswith("GRCh37") else "chr"
    return ",".join(f"{prefix}{i}" for i in range(1, 23))


# --------------------------------------------------------------------------
# The manifest files this run reads, each with the sha256 the manifest pins
# ("" = not pinned). A shared resources_dir may already hold a file another
# checkout fetched, so what fetch_resource recorded for it is checked here,
# before any job is submitted; runinfo checks again once everything is fetched.
# --------------------------------------------------------------------------
USED_RESOURCES = {
    p: DEST_SHA256[p.removeprefix(f"{RES}/")]
    for p in sorted(
        {p for q, _, t in SMVAR_ARMS for p in (truth_vcf(q, t), truth_bed(q, t), strat_tar(q, t))}
        | {p for q, _, t in SV_ARMS for p in (sv_truth_vcf(q, t), sv_truth_bed(q, t))}
    )
}
for _path, _pinned in USED_RESOURCES.items():
    _sidecar = Path(f"{_path}.sha256")
    _got = _sidecar.read_text().strip() if _sidecar.exists() else ""
    if _pinned and _got and _got != _pinned:
        _fail(
            f"{_path} has sha256 {_got}, but {MANIFEST_PATH} pins {_pinned}. "
            "Delete the file and its .sha256 to fetch it again, or correct the manifest."
        )


# --------------------------------------------------------------------------
# Wildcards. `dest` is pinned to exactly what the manifest declares, so the
# generic fetch rule can never collide with a derived resource (SDF, strat dir).
# --------------------------------------------------------------------------
wildcard_constraints:
    dest=_alternation(DEST_URL),
    query=_alternation(QUERIES),
    caller=_alternation(CALLERS),
    kind="small|sv",
    sample="[A-Za-z0-9]+",
    truth=_alternation(BENCH),
    reference="[A-Za-z0-9._]+",
    version="[A-Za-z0-9._]+",
    subset="[A-Za-z0-9._-]+",
