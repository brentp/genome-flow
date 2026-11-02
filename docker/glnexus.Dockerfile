# Follow the v1.4.5 upstream build environment; no vendor source or wrapper scripts.
FROM ubuntu:18.04@sha256:152dc042452c496007f07ca9127571cb9c29697f42acbfad72324b2bb2e43c98 AS builder
ENV DEBIAN_FRONTEND=noninteractive LC_ALL=C.UTF-8 LANG=C.UTF-8
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    ca-certificates git g++ cmake autoconf make curl wget file \
    libjemalloc-dev libzip-dev libsnappy-dev libbz2-dev zlib1g-dev liblzma-dev libzstd-dev \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 --branch v1.4.5 https://github.com/dnanexus-rnd/GLnexus.git /GLnexus \
    && test "$(git -C /GLnexus rev-parse HEAD)" = fae491bcd3e19ac9d5095ccb8bb98643fa3f47bd
WORKDIR /GLnexus
# Original fcmm repository is gone. Six surviving forks preserve v1.0.1 at
# this commit; change only its download location and verify the archive hash.
RUN sed -i 's@URL https://github.com/giacomodrago/fcmm/archive/v1.0.1.zip@URL https://github.com/asmodehn/fcmm/archive/51935acee053208aed166ebd181417d4b6dd6a9f.zip\n    URL_HASH SHA256=df66c871a1184bc4e20cdb4d82bedaa3a2f8b3609d9e333fb8beccd22914e74d@' CMakeLists.txt
RUN cmake -DCMAKE_BUILD_TYPE=Release . && make -j4 \
    && apt-get update -qq && apt-get install -y --no-install-recommends python3-pyvcf \
    && ctest -R '^unit_tests$' --output-on-failure

FROM ubuntu:20.04@sha256:8feb4d8ca5354def3d8fce243717141ce31e2c428701f6682bd2fafe15388214
ENV DEBIAN_FRONTEND=noninteractive LC_ALL=C.UTF-8 LANG=C.UTF-8
RUN apt-get update -qq && apt-get install -y --no-install-recommends libjemalloc2 \
    && rm -rf /var/lib/apt/lists/*
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2
COPY --from=builder /GLnexus/glnexus_cli /usr/local/bin/
COPY --from=builder /GLnexus/LICENSE /usr/share/doc/glnexus/LICENSE
LABEL org.opencontainers.image.source="https://github.com/dnanexus-rnd/GLnexus" \
      org.opencontainers.image.revision="fae491bcd3e19ac9d5095ccb8bb98643fa3f47bd" \
      org.opencontainers.image.version="1.4.5"
CMD ["glnexus_cli"]
