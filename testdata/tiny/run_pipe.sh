#!/usr/bin/env bash
# Exercise streaming and SHM reuse in a fresh, disposable Docker IPC namespace.
# Do not add --ipc=host: the SHM test owns and cleans up this private namespace.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
IMAGE="${IMAGE:-brentp/genome-flow-align:0.8.0}"
docker run --rm -v "$ROOT:/src:ro" -w /src -e GENOME_FLOW_INTEGRATION=1 -e GENOME_FLOW_SHM_TEST=1 \
  "$IMAGE" python3 -m unittest discover -s tests -v
