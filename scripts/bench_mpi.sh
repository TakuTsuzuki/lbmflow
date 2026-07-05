#!/usr/bin/env bash
# Weak-scaling table driver (docs/MPI_GUIDE.md): 512^2 per rank over
# ranks {1, 2, 4, 8}, aggregating examples/bench_mpi.rs RESULT lines into a
# table with parallel efficiency relative to n = 1.
#
# Single-node numbers measure shared-memory MPI only — see the caveats in
# examples/bench_mpi.rs and docs/MPI_GUIDE.md.

set -uo pipefail

MPI_HOME="${MPI_HOME:-$HOME/.local/openmpi}"
export PATH="$MPI_HOME/bin:$PATH"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

LOCAL="${LOCAL:-512}"
STEPS="${STEPS:-200}"
RANKS=(${RANKS:-1 2 4 8})

command -v mpirun >/dev/null || { echo "mpirun not found (MPI_HOME=$MPI_HOME)" >&2; exit 1; }

echo "== build (--features mpi) =="
cargo build -p lbm-core2 --release --features mpi --example bench_mpi || exit 1
BIN="$ROOT/target/release/examples/bench_mpi"

results=()
for n in "${RANKS[@]}"; do
    echo "== mpirun -n $n bench_mpi $LOCAL $STEPS =="
    line=$(mpirun --oversubscribe -n "$n" "$BIN" "$LOCAL" "$STEPS" | grep '^RESULT') || exit 1
    echo "$line"
    results+=("$line")
done

echo
echo "Weak scaling (${LOCAL}^2 per rank, $STEPS steps, f64 D2Q9 TGV, serial backend per rank):"
printf '%s\n' "${results[@]}" | awk '
    {
        for (i = 1; i <= NF; i++) {
            split($i, kv, "=");
            if (kv[1] == "ranks") n = kv[2];
            if (kv[1] == "mlups_total") m = kv[2];
            if (kv[1] == "time_s") t = kv[2];
        }
        if (base == 0) base = m / n;
        printf "  ranks=%-3s time=%7.3fs  MLUPS_total=%8.1f  MLUPS/rank=%7.1f  efficiency=%5.1f%%\n",
               n, t, m, m / n, 100.0 * (m / n) / base;
    }' base=0
