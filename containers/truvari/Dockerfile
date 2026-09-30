# truvari 5.4.0 from PyPI, every Python dependency pinned with hashes in
# requirements.lock to the versions of the bioconda build -- see the note there.
# Debian packages come from the live archive, so the SIF digest in
# workflow/containers.tsv is the exact pin. Refine's default aligner is POA (the
# pyabpoa package); nothing external is called unless `--align mafft` is asked
# for, which thoth never does.
FROM python:3.11-slim-bookworm@sha256:a36c24f9cbdf4fd0f52d67f0823eeac19c2028c637cecc392d97f980d4fec56b AS build
RUN apt-get update && apt-get install -y --no-install-recommends gcc g++ make zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.lock build-constraints.txt /tmp/
# The pip in this base (24.0) applies PIP_CONSTRAINT inside isolated builds too.
ENV PIP_CONSTRAINT=/tmp/build-constraints.txt
RUN pip wheel --no-cache-dir --require-hashes -w /wheels -r /tmp/requirements.lock

# Install exactly the wheels built above: sdists built here would not match the
# PyPI hashes in the lock, so the lock is not read a second time.
FROM python:3.11-slim-bookworm@sha256:a36c24f9cbdf4fd0f52d67f0823eeac19c2028c637cecc392d97f980d4fec56b
COPY --from=build /wheels /wheels
RUN pip install --no-cache-dir --no-index --no-deps /wheels/*.whl && pip check && rm -rf /wheels
LABEL org.opencontainers.image.source=https://github.com/martinandclaude/thoth \
      org.opencontainers.image.description="truvari 5.4.0 for thoth"
