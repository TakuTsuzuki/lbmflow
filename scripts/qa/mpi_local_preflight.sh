#!/usr/bin/env bash
# Local MPI preflight for the ME-3 cluster campaign.
#
# This intentionally uses tiny 3D weak-scaling cases (16^3 per rank by
# default) so the workstation sanity check finishes quickly. It validates the
# new D3Q19 bench_mpi path at n = 1, 2, 4 only; it is not a cluster scaling
# measurement.

set -euo pipefail

MPI_HOME="${MPI_HOME:-$HOME/.local/openmpi}"
export PATH="$MPI_HOME/bin:$PATH"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

LOCAL="${LOCAL:-16}"
STEPS="${STEPS:-5}"
WARMUP="${WARMUP:-1}"
RANKS=(${RANKS:-1 2 4})

if ! command -v mpirun >/dev/null || ! command -v mpicc >/dev/null; then
    echo "FAIL: mpirun/mpicc not found (MPI_HOME=$MPI_HOME)" >&2
    exit 1
fi
if [ "$(uname -sm)" = "Darwin arm64" ] && ! file "$(command -v mpirun)" | grep -q arm64; then
    echo "FAIL: $(command -v mpirun) is not an arm64 binary (PATH order?)" >&2
    file "$(command -v mpirun)" >&2
    exit 1
fi

echo "using $(command -v mpirun): $(mpirun --version | head -1)"
echo "== build (--features mpi) =="
cargo build -p lbm-core --release --features mpi --example bench_mpi

BIN="$ROOT/target/release/examples/bench_mpi"
LOG_DIR="${LOG_DIR:-$ROOT/target/mpi-preflight}"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/mpi_local_preflight.log"
: > "$LOG"

for n in "${RANKS[@]}"; do
    echo "== mpirun -n $n bench_mpi --mode weak3d --local-edge $LOCAL --steps $STEPS --warmup $WARMUP =="
    mpirun --oversubscribe -n "$n" "$BIN" \
        --mode weak3d --local-edge "$LOCAL" --steps "$STEPS" --warmup "$WARMUP" \
        | tee -a "$LOG"
done

echo "== aggregate =="
"$ROOT/scripts/qa/aggregate_mpi_results.sh" "$LOG" | tee "$LOG_DIR/mpi_local_preflight.csv"
echo "mpi_local_preflight: PASS (log: $LOG)"
