"""Helpers shared by the rules, the python scripts and fetch_containers.sh."""

import hashlib
import logging
import sys
from pathlib import Path

import pandas as pd

CONTAINERS_TSV = Path(__file__).resolve().parents[1] / "containers.tsv"


def containers():
    """workflow/containers.tsv: one row per tool name."""
    return pd.read_csv(CONTAINERS_TSV, sep="\t", comment="#", dtype=str).fillna("")


def sif_filename(image, digest):
    """<image basename>-<sha256 hex>.sif. Keyed by digest, so a rebuilt image is
    never mistaken for the one it replaces."""
    return f"{image.rsplit('/', 1)[-1]}-{digest.split(':', 1)[-1]}.sif"


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def log_to(path):
    """Log to the rule's log file and to stderr, and return the logging function.
    Snakemake does not redirect a python script's stdout/stderr, so each script
    opens its own log."""
    logging.basicConfig(
        level=logging.INFO,
        format="[thoth] %(message)s",
        handlers=[logging.FileHandler(path, "w"), logging.StreamHandler()],
        force=True,
    )
    return logging.info


def die(msg):
    logging.error(msg)
    sys.exit(1)
