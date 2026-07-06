#!/usr/bin/env bash
# T13-MPI verification driver (docs/MPI_GUIDE.md).
#
# Builds lbm-core with `--features mpi` and runs the distributed-vs-
# monolithic equivalence example under mpirun -n {1, 2, 3, 4, 5, 6, 8}.
# The example chooses a surface-minimising decomposition for each rank count;
# n = 3, 5, 6 cover non-dividing dimensions. The script also runs a 2-rank
# negative test that must reject mismatched viscosity with an explicit error.
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

echo "== build (--features mpi, debug dirty-guard check) =="
cargo build -p lbm-core --features mpi --example mpi_t13 || exit 1
DEBUG_BIN="$ROOT/target/debug/examples/mpi_t13"

echo "== build (--features mpi) =="
cargo build -p lbm-core --release --features mpi --example mpi_t13 || exit 1
BIN="$ROOT/target/release/examples/mpi_t13"

fail=0
for n in 1 2 3 4 5 6 8; do
    echo "== mpirun -n $n mpi_t13 =="
    if ! mpirun --oversubscribe -n "$n" "$BIN" "$@"; then
        fail=1
    fi
done

echo "== mpirun -n 2 mpi_t13 mismatch-nu (expected failure) =="
tmp="$(mktemp)"
if mpirun --oversubscribe -n 2 "$BIN" mismatch-nu >"$tmp" 2>&1; then
    echo "FAIL: mismatch-nu unexpectedly succeeded" >&2
    cat "$tmp" >&2
    fail=1
elif ! grep -q "MPI rank specification mismatch: nu differs across ranks" "$tmp"; then
    echo "FAIL: mismatch-nu did not name the offending item" >&2
    cat "$tmp" >&2
    fail=1
fi
rm -f "$tmp"

echo "== mpirun -n 2 mpi_t13 dirty-mismatch (debug expected failure) =="
tmp="$(mktemp)"
if mpirun --oversubscribe -n 2 "$DEBUG_BIN" dirty-mismatch >"$tmp" 2>&1; then
    echo "FAIL: dirty-mismatch unexpectedly succeeded in debug" >&2
    cat "$tmp" >&2
    fail=1
elif ! grep -q "MpiSolver host_dirty mismatch across ranks before step" "$tmp"; then
    echo "FAIL: dirty-mismatch did not trip the host_dirty debug assertion" >&2
    cat "$tmp" >&2
    fail=1
fi
rm -f "$tmp"

echo "== mpirun -n 2 mpi_t13 dirty-mismatch (release assert compiled out) =="
if ! mpirun --oversubscribe -n 2 "$BIN" dirty-mismatch; then
    echo "FAIL: dirty-mismatch failed in release; debug assertion should be compiled out" >&2
    fail=1
fi

echo
if [ "$fail" -eq 0 ]; then
    echo "test_mpi.sh: ALL PASS"
else
    echo "test_mpi.sh: FAILURES DETECTED" >&2
    exit 1
fi
