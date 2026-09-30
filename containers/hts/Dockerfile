# bcftools + samtools + htslib 1.24, built from the upstream release tarballs.
# Serves the `bcftools` and `samtools` names in workflow/containers.tsv.
# Pinned: base image by digest, release tarballs by sha256. Debian packages come
# from the live bookworm archive, so a rebuild is functionally equivalent
# (tests/run.sh), not byte-identical; the SIF digest in workflow/containers.tsv
# is the exact pin.
FROM debian:bookworm-slim@sha256:3783cc01769c7b2b1b83a5c5ad96c815348e28ed7da68e2e3687004faa906251 AS build
ARG HTS_VERSION=1.24
ARG BCFTOOLS_SHA256=8caddc22610ee2851666047c859bb91da0c1e32d0c2ec553db6f153ad130e46f
ARG SAMTOOLS_SHA256=89b2a440123eeaa400392ce1736e7d60ce9041843027d76819753c5a8246bfdd
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl bzip2 gcc make libc6-dev zlib1g-dev libbz2-dev \
        liblzma-dev libcurl4-openssl-dev libdeflate-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
RUN curl -fsSL -o bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/${HTS_VERSION}/bcftools-${HTS_VERSION}.tar.bz2 \
    && echo "${BCFTOOLS_SHA256}  bcftools.tar.bz2" | sha256sum -c - \
    && curl -fsSL -o samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/${HTS_VERSION}/samtools-${HTS_VERSION}.tar.bz2 \
    && echo "${SAMTOOLS_SHA256}  samtools.tar.bz2" | sha256sum -c - \
    && tar -xjf bcftools.tar.bz2 && tar -xjf samtools.tar.bz2
# htslib comes bundled in the bcftools tarball; installing it gives bgzip/tabix.
RUN cd bcftools-${HTS_VERSION}/htslib-${HTS_VERSION} && ./configure --prefix=/opt/hts --enable-libcurl && make -j"$(nproc)" && make install \
    && cd /src/bcftools-${HTS_VERSION} && ./configure --prefix=/opt/hts && make -j"$(nproc)" && make install \
    && cd /src/samtools-${HTS_VERSION} && ./configure --prefix=/opt/hts --without-curses && make -j"$(nproc)" && make install
# Their licences (MIT/Expat; modified BSD for htslib's cram/) require the
# notices to ship with the binaries.
RUN install -Dm644 /src/bcftools-${HTS_VERSION}/LICENSE /licenses/bcftools/LICENSE \
    && install -Dm644 /src/bcftools-${HTS_VERSION}/htslib-${HTS_VERSION}/LICENSE /licenses/htslib/LICENSE \
    && install -Dm644 /src/samtools-${HTS_VERSION}/LICENSE /licenses/samtools/LICENSE

FROM debian:bookworm-slim@sha256:3783cc01769c7b2b1b83a5c5ad96c815348e28ed7da68e2e3687004faa906251
RUN apt-get update && apt-get install -y --no-install-recommends \
        zlib1g libbz2-1.0 liblzma5 libcurl4 libdeflate0 ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/hts /opt/hts
COPY --from=build /licenses/ /usr/share/doc/
ENV PATH=/opt/hts/bin:$PATH
LABEL org.opencontainers.image.source=https://github.com/martinandclaude/thoth \
      org.opencontainers.image.description="bcftools, samtools and htslib 1.24 for thoth"
