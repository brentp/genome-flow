# Build from repository root: docker build -f docker/align.Dockerfile .
FROM rust:1-bookworm@sha256:82150a52ec202c1b14d7817e14516c392bb7f5cfebd88f1ed531cb37ebd39922 AS rustbuild
RUN cargo install --locked dupblaster --version 0.3.0 \
    && cargo install --locked fg-mako --version 0.1.7

FROM debian:trixie-slim@sha256:f6b08b09a05725f9431502d3b5714f5c4576e82c1a84fa1debe0c24a7e0c3ec4 AS cbuild
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential clang-19 libomp-19-dev cmake pkg-config \
    autoconf automake autoconf-archive libtool zlib1g-dev libdeflate-dev ca-certificates \
    && rm -rf /var/lib/apt/lists/*
# Release tarball with submodules bundled; checksum from the v0.16.0 release notes.
ADD https://github.com/fg-labs/bwa-mem3/releases/download/v0.16.0/Source_code_including_submodules.tar.gz /tmp/bwa-mem3.tar.gz
RUN echo "e138fbd5cc32538c534240b9223f418c4612bf211cf35650ed97a4ee9423e016  /tmp/bwa-mem3.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/bwa-mem3.tar.gz -C /tmp \
    && make -C /tmp/bwa-mem3-0.16.0 -j4 CXX=clang++-19 CC=clang-19 \
    && cp /tmp/bwa-mem3-0.16.0/bwa-mem3 /usr/local/bin/

FROM debian:trixie-slim@sha256:f6b08b09a05725f9431502d3b5714f5c4576e82c1a84fa1debe0c24a7e0c3ec4
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools=1.21-1 bcftools=1.21-1 python3=3.13.5-1 ca-certificates procps libgomp1 libomp5-19 libdeflate0 libjemalloc2 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=rustbuild /usr/local/cargo/bin/dupblaster /usr/local/bin/dupblaster
COPY --from=rustbuild /usr/local/cargo/bin/mako /usr/local/bin/mako
COPY --from=cbuild /usr/local/bin/bwa-mem3 /usr/local/bin/bwa-mem3
COPY scripts/stream_align.py /usr/local/bin/genome-flow-stream-align
COPY scripts/make_intervals.py /usr/local/bin/genome-flow-make-intervals
RUN bwa-mem3 version && dupblaster --version && mako --version \
    && samtools --version | head -n 1 && bcftools --version | head -n 1

# Upstream static x86-64 release; checksum makes a changed asset fail the build.
ADD https://github.com/brentp/mosdepth/releases/download/v0.3.14/mosdepth /usr/local/bin/mosdepth
RUN echo "c5182b74a8f1b66710efa16e122cbc8a197834874b103e7c5c0bd9a6265ae7b6  /usr/local/bin/mosdepth" | sha256sum -c - \
    && chmod 755 /usr/local/bin/mosdepth && mosdepth --version

# Upstream static Somalier release for identity QC in the shared QC task.
ADD https://github.com/brentp/somalier/releases/download/v0.3.6/somalier /usr/local/bin/somalier
RUN echo "6d25bbf7210f2f8c20783cc0b438a5a319dcb253a5d8c9d030954fb3eb1fbc1f  /usr/local/bin/somalier" | sha256sum -c - \
    && chmod 755 /usr/local/bin/somalier && somalier 2>&1 | head -n 3

# Reuse the already tested source build; no second GLnexus compilation here.
COPY --from=brentp/glnexus:1.4.5@sha256:92da22de417ce4ef2c72910928492e74a35e1b0e11a956f77d4fdcce11debfc4 /usr/local/bin/glnexus_cli /usr/local/bin/
COPY --from=brentp/glnexus:1.4.5@sha256:92da22de417ce4ef2c72910928492e74a35e1b0e11a956f77d4fdcce11debfc4 /usr/share/doc/glnexus/LICENSE /usr/share/doc/glnexus/LICENSE
RUN glnexus_cli --help
