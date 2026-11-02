#!/usr/bin/env bash
# Full workflow on a synthetic two-sample fixture, scattered over one and two
# DeepVariant shards; compare.py checks the outputs agree. Needs Docker and miniwdl.
# IMAGE overrides the tools image, e.g. a local build of docker/align.Dockerfile.
set -euo pipefail
cd "$(dirname "$0")/../.."
python3 testdata/boundary/make_fixture.py
fixture="$(pwd)/testdata/boundary/runs/fixture"
image="${IMAGE:-brentp/genome-flow-align:0.8.0}"
docker run --rm -v "$fixture:/work" -w /work "$image" \
  sh -c 'samtools faidx ref.fa && bwa-mem3 index ref.fa'
for count in 1 2; do
  miniwdl run wdl/genome_flow.wdl align_docker="$image" -i "$fixture/inputs.$count.json" \
    --dir "testdata/boundary/runs/scatter-$count"
done
docker run --rm -v "$(pwd):$(pwd):ro" -w "$(pwd)" "$image" \
  python3 testdata/boundary/compare.py \
  testdata/boundary/runs/scatter-1/_LAST/outputs.json \
  testdata/boundary/runs/scatter-2/_LAST/outputs.json
