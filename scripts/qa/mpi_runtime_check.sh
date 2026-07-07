#!/usr/bin/env bash
# Ready-to-run MPI runtime V&V harness.
#
# This wrapper keeps scripts/test_mpi.sh as the source of truth for T13-MPI
# correctness, adds an MPI launcher preflight, and captures optional scaling
# artifacts under target/qa/mpi-runtime/.

set -uo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT" || exit 1

MPI_HOME="${MPI_HOME:-$HOME/.local/openmpi}"
export PATH="$MPI_HOME/bin:$PATH"

MPIRUN="${MPIRUN:-mpirun}"
MPICC="${MPICC:-mpicc}"
RANK_COUNTS="${RANK_COUNTS:-${MPI_RANKS:-1 2 4 8}}"
SCALING_RANKS="${SCALING_RANKS:-1 2 4 8}"
RUN_SCALING="${RUN_SCALING:-1}"
WEAK_LOCAL_EDGE="${WEAK_LOCAL_EDGE:-32}"
STRONG_GLOBAL_EDGE="${STRONG_GLOBAL_EDGE:-64}"
SCALING_STEPS="${SCALING_STEPS:-20}"
SCALING_WARMUP="${SCALING_WARMUP:-5}"
MPI_LAUNCH_EXTRA="${MPI_LAUNCH_EXTRA:---oversubscribe}"
DRY_RUN=0
LOG_DIR="${LOG_DIR:-$ROOT/target/qa/mpi-runtime/$(date -u +%Y%m%dT%H%M%SZ)}"

usage() {
    cat <<'EOF'
Usage: scripts/qa/mpi_runtime_check.sh [--dry-run] [--skip-scaling]
       [--rank-counts "1 2 4 8"] [--scaling-ranks "1 2 4 8"]
       [--log-dir DIR]

Environment overrides:
  MPI_HOME             MPI prefix prepended to PATH (default: $HOME/.local/openmpi)
  MPIRUN, MPICC        launcher/compiler wrapper names (default: mpirun/mpicc)
  RANK_COUNTS          T13-MPI rank counts (default: 1 2 4 8)
  RUN_SCALING          1 to run weak/strong scaling, 0 to skip (default: 1)
  SCALING_RANKS        scaling rank counts (default: 1 2 4 8)
  WEAK_LOCAL_EDGE      weak3d local edge per rank (default: 32)
  STRONG_GLOBAL_EDGE   strong3d global edge (default: 64)
  SCALING_STEPS        scaling measurement steps (default: 20)
  SCALING_WARMUP       scaling warmup steps (default: 5)
  MPI_LAUNCH_EXTRA     extra args for benchmark mpirun (default: --oversubscribe)
EOF
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --dry-run)
            DRY_RUN=1
            ;;
        --skip-scaling)
            RUN_SCALING=0
            ;;
        --with-scaling)
            RUN_SCALING=1
            ;;
        --rank-counts)
            shift
            RANK_COUNTS="${1:-}"
            ;;
        --scaling-ranks)
            shift
            SCALING_RANKS="${1:-}"
            ;;
        --log-dir)
            shift
            LOG_DIR="${1:-}"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "FAIL: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
    shift
done

print_setup_steps() {
    cat <<EOF

MPI runtime setup is required before this harness can execute:

  macOS arm64 / local Open MPI:
    mkdir -p ~/.local/src && cd ~/.local/src
    curl -sL https://download.open-mpi.org/release/open-mpi/v5.0/openmpi-5.0.9.tar.bz2 | tar xj
    cd openmpi-5.0.9
    ./configure --prefix=\$HOME/.local/openmpi CC=clang CXX=clang++ --disable-mpi-fortran
    make -j\$(sysctl -n hw.ncpu) && make install
    export PATH=\$HOME/.local/openmpi/bin:\$PATH

  Linux / AWS hpc7g:
    sudo apt-get update
    sudo apt-get install -y openmpi-bin libopenmpi-dev
    # or: sudo apt-get install -y mpich libmpich-dev

  Cluster module environment:
    module avail mpi
    module load openmpi     # or mpich / cray-mpich / site MPI module
    export MPIRUN=mpirun MPICC=mpicc

Then rerun:
  scripts/qa/mpi_runtime_check.sh
EOF
}

mpi_path() {
    command -v "$1" 2>/dev/null || true
}

print_environment_report() {
    local mpirun_path mpicc_path version
    mpirun_path="$(mpi_path "$MPIRUN")"
    mpicc_path="$(mpi_path "$MPICC")"

    echo "== MPI runtime environment =="
    echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "root=$ROOT"
    echo "uname=$(uname -a)"
    echo "MPI_HOME=$MPI_HOME"
    echo "MPIRUN=$MPIRUN"
    echo "MPICC=$MPICC"
    echo "mpirun_path=${mpirun_path:-not-found}"
    echo "mpicc_path=${mpicc_path:-not-found}"
    echo "rank_counts=$RANK_COUNTS"
    echo "run_scaling=$RUN_SCALING"
    echo "scaling_ranks=$SCALING_RANKS"
    echo "log_dir=$LOG_DIR"

    if [ -n "$mpirun_path" ]; then
        version="$("$MPIRUN" --version 2>&1 | head -5)"
        echo "mpirun_version:"
        printf '%s\n' "$version" | sed 's/^/  /'
        if printf '%s\n' "$version" | grep -qi 'open mpi'; then
            echo "mpi_vendor=openmpi"
        elif printf '%s\n' "$version" | grep -Eqi 'mpich|hydra'; then
            echo "mpi_vendor=mpich-family"
        else
            echo "mpi_vendor=unknown"
        fi
    fi

    if [ "$(uname -sm)" = "Darwin arm64" ] && [ -n "$mpirun_path" ]; then
        echo "mpirun_file:"
        file "$mpirun_path" | sed 's/^/  /'
    fi
}

missing_mpi=0
print_environment_report
if ! command -v "$MPIRUN" >/dev/null 2>&1 || ! command -v "$MPICC" >/dev/null 2>&1; then
    missing_mpi=1
    echo
    echo "MPI status: unavailable (mpirun and mpicc are both required)."
    print_setup_steps
fi

if [ "$DRY_RUN" -eq 1 ]; then
    echo
    echo "DRY-RUN: execution skipped."
    echo "Would verify launcher: $MPIRUN -n 2 echo hi"
    echo "Would run T13-MPI: MPI_RANKS=\"$RANK_COUNTS\" scripts/test_mpi.sh"
    if [ "$RUN_SCALING" -eq 1 ]; then
        echo "Would run weak scaling: ranks {$SCALING_RANKS}, weak3d local-edge=$WEAK_LOCAL_EDGE"
        echo "Would run strong scaling: ranks {$SCALING_RANKS}, strong3d global-edge=$STRONG_GLOBAL_EDGE"
    else
        echo "Would skip scaling: RUN_SCALING=0"
    fi
    exit 0
fi

if [ "$missing_mpi" -ne 0 ]; then
    exit 2
fi

mkdir -p "$LOG_DIR"

echo
echo "== mpirun hello preflight =="
HELLO_LOG="$LOG_DIR/hello.log"
"$MPIRUN" -n 2 echo hi >"$HELLO_LOG" 2>&1
hello_status=$?
cat "$HELLO_LOG"
if [ "$hello_status" -ne 0 ]; then
    if grep -qi 'Operation not permitted' "$HELLO_LOG"; then
        echo "FAIL: mpirun could not bind sockets in this environment."
        echo "Known issue: macOS sandboxed sessions can block MPI socket binding."
        echo "Run this harness off-sandbox or on the target cluster."
    else
        echo "FAIL: mpirun hello preflight failed."
    fi
    exit 1
fi

echo
echo "== T13-MPI runtime check =="
T13_LOG="$LOG_DIR/t13-mpi.log"
env MPI_RANKS="$RANK_COUNTS" "$ROOT/scripts/test_mpi.sh" >"$T13_LOG" 2>&1
t13_status=$?
cat "$T13_LOG"

validate_t13_log() {
    local log="$1"
    local fail=0
    local n expected total exact

    for n in $RANK_COUNTS; do
        if ! grep -Fq "mpi_t13 [n=$n]: ALL PASS" "$log"; then
            echo "FAIL: missing T13-MPI ALL PASS line for rank count $n" >&2
            fail=1
        fi
    done

    if grep -Eq 'FAILURES DETECTED|^FAIL:|^FAIL ' "$log"; then
        echo "FAIL: T13-MPI log contains a failure marker" >&2
        fail=1
    fi

    expected=0
    for n in $RANK_COUNTS; do
        expected=$((expected + 4))
        if [ "$n" = "6" ] || [ "$n" = "8" ]; then
            expected=$((expected + 1))
        fi
    done

    total="$(awk '/^PASS / && /max field/ { c++ } END { print c + 0 }' "$log")"
    exact="$(awk '/^PASS / && /max field/ && /max field .*0\.0e[+]?0/ { c++ } END { print c + 0 }' "$log")"
    if [ "$total" -lt "$expected" ] || [ "$exact" -ne "$total" ]; then
        echo "FAIL: T13 stored-reference field comparison was not bit-identical in every PASS line" >&2
        echo "      expected_at_least=$expected pass_lines=$total exact_field_zero=$exact" >&2
        fail=1
    else
        echo "T13 stored-reference field comparison: PASS ($exact/$total case lines report max field diff 0.0)"
    fi

    return "$fail"
}

validate_t13_log "$T13_LOG"
reference_status=$?
if [ "$t13_status" -ne 0 ] || [ "$reference_status" -ne 0 ]; then
    echo "mpi_runtime_check: FAIL (T13-MPI runtime validation)"
    exit 1
fi

scaling_status=0
if [ "$RUN_SCALING" -eq 1 ]; then
    echo
    echo "== MPI weak/strong scaling smoke =="
    cargo build -p lbm-core --release --features mpi --example bench_mpi >"$LOG_DIR/bench-build.log" 2>&1
    build_status=$?
    cat "$LOG_DIR/bench-build.log"
    if [ "$build_status" -ne 0 ]; then
        echo "FAIL: bench_mpi build failed" >&2
        exit 1
    fi

    BENCH_BIN="$ROOT/target/release/examples/bench_mpi"
    SCALING_LOG="$LOG_DIR/scaling.log"
    : > "$SCALING_LOG"

    run_scaling_case() {
        local mode="$1"
        local ranks="$2"
        local status=0
        local n
        for n in $ranks; do
            if [ "$mode" = "weak3d" ]; then
                echo "== $MPIRUN $MPI_LAUNCH_EXTRA -n $n bench_mpi --mode weak3d --local-edge $WEAK_LOCAL_EDGE --steps $SCALING_STEPS --warmup $SCALING_WARMUP =="
                # shellcheck disable=SC2086
                "$MPIRUN" $MPI_LAUNCH_EXTRA -n "$n" "$BENCH_BIN" \
                    --mode weak3d --local-edge "$WEAK_LOCAL_EDGE" \
                    --steps "$SCALING_STEPS" --warmup "$SCALING_WARMUP" \
                    >>"$SCALING_LOG" 2>&1 || status=1
            else
                echo "== $MPIRUN $MPI_LAUNCH_EXTRA -n $n bench_mpi --mode strong3d --global-edge $STRONG_GLOBAL_EDGE --steps $SCALING_STEPS --warmup $SCALING_WARMUP =="
                # shellcheck disable=SC2086
                "$MPIRUN" $MPI_LAUNCH_EXTRA -n "$n" "$BENCH_BIN" \
                    --mode strong3d --global-edge "$STRONG_GLOBAL_EDGE" \
                    --steps "$SCALING_STEPS" --warmup "$SCALING_WARMUP" \
                    >>"$SCALING_LOG" 2>&1 || status=1
            fi
        done
        return "$status"
    }

    run_scaling_case weak3d "$SCALING_RANKS" || scaling_status=1
    run_scaling_case strong3d "$SCALING_RANKS" || scaling_status=1
    cat "$SCALING_LOG"

    if [ "$scaling_status" -eq 0 ]; then
        echo "== scaling aggregate =="
        RANK_CSV="$LOG_DIR/scaling-ranks.csv" \
            "$ROOT/scripts/qa/aggregate_mpi_results.sh" "$SCALING_LOG" \
            | tee "$LOG_DIR/scaling-summary.csv"
    else
        echo "FAIL: one or more scaling cases failed; see $SCALING_LOG" >&2
    fi
else
    echo
    echo "== MPI weak/strong scaling smoke =="
    echo "SKIP: RUN_SCALING=0"
fi

if [ "$scaling_status" -ne 0 ]; then
    echo "mpi_runtime_check: FAIL (scaling smoke)"
    exit 1
fi

echo
echo "mpi_runtime_check: PASS"
echo "artifacts=$LOG_DIR"
