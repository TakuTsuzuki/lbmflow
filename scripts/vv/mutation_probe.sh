#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/vv/mutation_probe.sh --list
  bash scripts/vv/mutation_probe.sh <mutation-id>

The script applies one exact temporary physics mutation, runs the configured
sentinel command, and restores every touched engine file on exit. A failing
sentinel means the mutation was killed. A passing sentinel means the mutation
survived and is a false-assurance finding.
USAGE
}

TARGETS=()
BACKUP_DIR=""

restore_targets() {
  local status=$?
  if [[ -n "${BACKUP_DIR}" && -d "${BACKUP_DIR}" ]]; then
    for file in "${TARGETS[@]:-}"; do
      if [[ -f "${BACKUP_DIR}/${file}" ]]; then
        cp "${BACKUP_DIR}/${file}" "${file}"
      fi
    done
    rm -rf "${BACKUP_DIR}"
  fi
  if ((${#TARGETS[@]} > 0)); then
    git diff --quiet -- "${TARGETS[@]}" || {
      echo "ERROR: mutation restore left a diff in target files:" >&2
      git diff -- "${TARGETS[@]}" >&2
      exit 99
    }
  fi
  exit "$status"
}

trap restore_targets EXIT INT TERM

mutation_catalog() {
  cat <<'CATALOG'
guo-f2-velocity-removed
forcing-sign-flipped
trt-relaxation-swapped
d2q9-opposite-broken
d3q19-opposite-broken
d3q27-face-unknown-broken
halfway-wall-shifted
moving-wall-sign-flipped
zou-he-pressure-normal-sign-flipped
pressure-outlet-correction-removed
outflow-stale-slot-broken
mpi-halo-x-direction-swapped
probe-force-physicalization-removed
shan-chen-force-sign-flipped
contact-angle-wall-term-sign-flipped
f32-deviation-storage-disabled
CATALOG
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if [[ "${1:-}" == "--list" ]]; then
  mutation_catalog
  exit 0
fi

MUTATION="${1:-}"
if [[ -z "$MUTATION" ]]; then
  usage >&2
  exit 2
fi

replace_exact() {
  local file="$1"
  local find="$2"
  local repl="$3"
  FIND="$find" REPL="$repl" perl -0pi -e '
    BEGIN {
      our $find = $ENV{"FIND"};
      our $repl = $ENV{"REPL"};
      our $count = 0;
    }
    our ($find, $repl, $count);
    $count += s/\Q$find\E/$repl/g;
    END {
      die "exact mutation pattern not found\n" if $count == 0;
    }
  ' "$file"
}

replace_regex() {
  local file="$1"
  local find_regex="$2"
  local repl="$3"
  FIND_REGEX="$find_regex" REPL="$repl" perl -0pi -e '
    BEGIN {
      our $find = $ENV{"FIND_REGEX"};
      our $repl = $ENV{"REPL"};
      $repl =~ s/\\n/\n/g;
      our $count = 0;
    }
    our ($find, $repl, $count);
    $count += s/$find/$repl/g;
    END {
      die "regex mutation pattern not found\n" if $count == 0;
    }
  ' "$file"
}

set_mutation() {
  DESC=""
  CMD=()
  TARGETS=()
  case "$MUTATION" in
    guo-f2-velocity-removed)
      DESC="Remove Guo F/2 from physical velocity diagnostics."
      TARGETS=(crates/lbm-core/src/backend.rs crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test validation_conservation t6_periodic_uniform_force_adds_exact_momentum -- --exact)
      ;;
    forcing-sign-flipped)
      DESC="Flip the sign of the Guo source term."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test validation_conservation t6_periodic_uniform_force_adds_exact_momentum -- --exact)
      ;;
    trt-relaxation-swapped)
      DESC="Use omega+ for the TRT antisymmetric relaxation."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test validation_channel t2_trt_magic_poiseuille_is_exact_and_symmetric -- --exact)
      ;;
    d2q9-opposite-broken)
      DESC="Break one D2Q9 opposite-direction table entry."
      TARGETS=(crates/lbm-core/src/lattice.rs)
      CMD=(cargo test --release -p lbm-core lattice::tests::d2q9_invariants -- --exact)
      ;;
    d3q19-opposite-broken)
      DESC="Break one D3Q19 opposite-direction table entry."
      TARGETS=(crates/lbm-core/src/lattice.rs)
      CMD=(cargo test --release -p lbm-core lattice::tests::d3q19_invariants -- --exact)
      ;;
    d3q27-face-unknown-broken)
      DESC="Break a D3Q27 face-unknown table entry."
      TARGETS=(crates/lbm-core/src/lattice.rs)
      CMD=(cargo test --release -p lbm-core lattice::tests::d3q27_invariants -- --exact)
      ;;
    halfway-wall-shifted)
      DESC="Shift the T2 half-way wall force reference by one cell in the stream bounce-back path."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test validation_channel t2_trt_magic_poiseuille_is_exact_and_symmetric -- --exact)
      ;;
    moving-wall-sign-flipped)
      DESC="Flip the moving-wall bounce-back momentum-injection sign."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test validation_channel t3_top_wall_couette_exact_for_bgk_and_trt_all_taus -- --exact)
      ;;
    zou-he-pressure-normal-sign-flipped)
      DESC="Flip the pressure-Zou-He normal velocity closure sign."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test t15_3d t15_1c_zou_he_3d_enforces_prescribed_moments -- --exact)
      ;;
    pressure-outlet-correction-removed)
      DESC="Remove the tangential correction in Zou-He reconstruction."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test validation_open_bc t4_velocity_inlet_pressure_outlet_channel_all_four_orientations -- --exact)
      ;;
    outflow-stale-slot-broken)
      DESC="Let streaming overwrite open-face unknown slots instead of preserving stale slots."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test stream_contract cpu_stream_preserves_open_face_unknowns_d2q9 -- --exact)
      ;;
    mpi-halo-x-direction-swapped)
      DESC="Swap MPI x-direction halo buffer assignment."
      TARGETS=(crates/lbm-core/src/dist.rs)
      CMD=(cargo test --release -p lbm-core --features mpi dist::tests::mpi_choose_decomp_basic -- --exact)
      ;;
    probe-force-physicalization-removed)
      DESC="Remove +2w physical-population conversion from momentum-exchange probe."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test accuracy_audit_probe a2_steady_poiseuille_wall_friction_balance -- --exact)
      ;;
    shan-chen-force-sign-flipped)
      DESC="Flip the single-component Shan-Chen cohesion force sign."
      TARGETS=(crates/lbm-core/src/compat/multiphase.rs)
      CMD=(cargo test --release -p lbm-core --test validation_multiphase t11_laplace_single_radius_smoke -- --exact)
      ;;
    contact-angle-wall-term-sign-flipped)
      DESC="Flip the legacy Shan-Chen wall-adhesion term sign."
      TARGETS=(crates/lbm-core/src/compat/multiphase.rs)
      CMD=(cargo test --release -p lbm-core --test validation_contact_angle t11b_wall_adhesion_contact_angles_are_monotone_and_frozen -- --exact)
      ;;
    f32-deviation-storage-disabled)
      DESC="Remove the +1 rest-state term from deviation-storage density reconstruction."
      TARGETS=(crates/lbm-core/src/kernels.rs)
      CMD=(cargo test --release -p lbm-core --test validation_conservation t6_f32_mass_and_momentum_hold_with_tightened_tolerance -- --exact)
      ;;
    *)
      echo "Unknown mutation id: $MUTATION" >&2
      echo "Known ids:" >&2
      mutation_catalog >&2
      exit 2
      ;;
  esac
}

apply_mutation() {
  case "$MUTATION" in
    guo-f2-velocity-removed)
      replace_exact crates/lbm-core/src/backend.rs "acc += m + 0.5 * fa;" "acc += m;"
      replace_exact crates/lbm-core/src/kernels.rs "let half = T::r(0.5);" "let half = T::zero();"
      ;;
    forcing-sign-flipped)
      replace_exact crates/lbm-core/src/kernels.rs "src[q] = p.wr[q] * (three * (cf - uf) + nine * cu * cf);" "src[q] = -p.wr[q] * (three * (cf - uf) + nine * cu * cf);"
      ;;
    trt-relaxation-swapped)
      replace_exact crates/lbm-core/src/kernels.rs "let rm = p.omega_m * (fm - em);" "let rm = op * (fm - em);"
      ;;
    d2q9-opposite-broken)
      replace_exact crates/lbm-core/src/lattice.rs "const D2Q9_OPP: [usize; 9] = opp_table(&D2Q9_C);" "const D2Q9_OPP: [usize; 9] = [0, 3, 4, 1, 2, 8, 7, 6, 5];"
      ;;
    d3q19-opposite-broken)
      replace_exact crates/lbm-core/src/lattice.rs "const D3Q19_OPP: [usize; 19] = opp_table(&D3Q19_C);" "const D3Q19_OPP: [usize; 19] = [0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17];"
      ;;
    d3q27-face-unknown-broken)
      replace_exact crates/lbm-core/src/lattice.rs "const D3Q27_UNK_XN: [usize; 9] = face_unknowns(&D3Q27_C, [1, 0, 0]);" "const D3Q27_UNK_XN: [usize; 9] = [2, 8, 10, 14, 16, 20, 22, 24, 26];"
      ;;
    halfway-wall-shifted)
      replace_exact crates/lbm-core/src/kernels.rs "let fout = f[L::OPP[q] * np + i];" "let fout = f[L::OPP[q] * np + si];"
      ;;
    moving-wall-sign-flipped)
      replace_exact crates/lbm-core/src/kernels.rs "let fin = fout + six * p.wr[q] * rho_row[x] * cu;" "let fin = fout - six * p.wr[q] * rho_row[x] * cu;"
      ;;
    zou-he-pressure-normal-sign-flipped)
      replace_exact crates/lbm-core/src/kernels.rs "let un = T::one() - closure / rho_bc;" "let un = closure / rho_bc - T::one();"
      ;;
    pressure-outlet-correction-removed)
      replace_exact crates/lbm-core/src/kernels.rs "let tcorr = half * (r * ut - (f[q_t * np + i] - f[q_mt * np + i]));" "let tcorr = T::zero();"
      replace_exact crates/lbm-core/src/kernels.rs "let n1 = c13 * r * ut1 - half * qt1;" "let n1 = T::zero();"
      replace_exact crates/lbm-core/src/kernels.rs "let n2 = c13 * r * ut2 - half * qt2;" "let n2 = T::zero();"
      ;;
    outflow-stale-slot-broken)
      replace_regex crates/lbm-core/src/kernels.rs "if !halo\\[2 \\* a\\] \\{\\s*continue 'dirs;\\s*\\}" "if !halo[2 * a] {\n                        // MUTATION: fall through and read the halo slot.\n                    }"
      replace_regex crates/lbm-core/src/kernels.rs "\\} else if s\\[a\\] >= geom\\.core\\[a\\] as isize && !halo\\[2 \\* a \\+ 1\\] \\{\\s*continue 'dirs;\\s*\\}" "} else if s[a] >= geom.core[a] as isize && !halo[2 * a + 1] {\n                    // MUTATION: fall through and read the halo slot.\n                }"
      ;;
    mpi-halo-x-direction-swapped)
      replace_exact crates/lbm-core/src/dist.rs "const TAG_F: Tag = 100;" "const TAG_F: Tag = 101;"
      ;;
    probe-force-physicalization-removed)
      replace_exact crates/lbm-core/src/kernels.rs "let ftot = fout + fin + two * p.wr[q];" "let ftot = fout + fin;"
      ;;
    shan-chen-force-sign-flipped)
      replace_exact crates/lbm-core/src/compat/multiphase.rs "-psi_i * (self.g * sx + self.g_wall * ax)," "psi_i * (self.g * sx + self.g_wall * ax),"
      replace_exact crates/lbm-core/src/compat/multiphase.rs "-psi_i * (self.g * sy + self.g_wall * ay)," "psi_i * (self.g * sy + self.g_wall * ay),"
      ;;
    contact-angle-wall-term-sign-flipped)
      replace_exact crates/lbm-core/src/compat/multiphase.rs "self.g * sx + self.g_wall * ax" "self.g * sx - self.g_wall * ax"
      replace_exact crates/lbm-core/src/compat/multiphase.rs "self.g * sy + self.g_wall * ay" "self.g * sy - self.g_wall * ay"
      ;;
    f32-deviation-storage-disabled)
      replace_exact crates/lbm-core/src/kernels.rs "let r = T::one() + dr;" "let r = dr;"
      ;;
  esac
}

set_mutation

for file in "${TARGETS[@]}"; do
  if ! git diff --quiet -- "$file"; then
    echo "ERROR: refusing to mutate dirty target file: $file" >&2
    exit 3
  fi
done

BACKUP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lbmflow-mutation.XXXXXX")"
for file in "${TARGETS[@]}"; do
  mkdir -p "${BACKUP_DIR}/$(dirname "$file")"
  cp "$file" "${BACKUP_DIR}/${file}"
done

echo "MUTATION: $MUTATION"
echo "DESCRIPTION: $DESC"
echo "COMMAND: ${CMD[*]}"

apply_mutation

set +e
"${CMD[@]}"
status=$?
set -e

if [[ $status -eq 0 ]]; then
  echo "RESULT: SURVIVED (sentinel command exited 0)"
  exit 10
else
  echo "RESULT: KILLED (sentinel command exited $status)"
  exit 0
fi
