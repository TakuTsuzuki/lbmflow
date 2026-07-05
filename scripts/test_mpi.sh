#!/usr/bin/env bash
# T13-MPI verification driver (docs/MPI_GUIDE.md).
#
# Builds lbm-core with `--features mpi` and runs the distributed-vs-
# monolithic equivalence example under mpirun -n {1, 2, 4, 8}:
#   n = 1, 2, 4 : 2D TGV / cavity (lid over the seam) / cylinder + probe on
#                 the seam / Shan-Chen droplet (2x2 corner at n = 4)
#   n = 8       : 3D TGV, 2x2x2 decomposition (D3Q19)
# Prints PASS/FAIL per case and exits non-zero on any failure.
#
# Requires a native (arm64 on Apple silicon) MPI. Default: the source-built
# Open MPI under $HOME/.local/openmpi (override with MPI_HOME). Note that an
# x86_64 MPI earlier on PATH (e.g. /usr/local) breaks both rsmpi's build
# probe and mpirun, hence the explicit PATH prefix + arch check.

set -uo pipefail

MPI_HOME="${MPI_HOME:-$HOME/.local/openmpi}"
export PATH="$MPI_HOME/bin:$PATH"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

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
cargo build -p lbm-core --release --features mpi --example mpi_t13 || exit 1
BIN="$ROOT/target/release/examples/mpi_t13"

fail=0
for n in 1 2 4 8; do
    echo "== mpirun -n $n mpi_t13 =="
    if ! mpirun --oversubscribe -n "$n" "$BIN" "$@"; then
        fail=1
    fi
done

echo
if [ "$fail" -eq 0 ]; then
    echo "test_mpi.sh: ALL PASS"
else
    echo "test_mpi.sh: FAILURES DETECTED" >&2
    exit 1
fi
