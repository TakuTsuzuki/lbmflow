#!/usr/bin/env bash
# Grep-based regression guard for the LBMFlow physics ban list.

set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

ROOT="$DEFAULT_ROOT"
WHITELIST=""
SELF_TEST=0
TARGET_OVERRIDE=0
TARGETS=()

usage() {
  cat <<'USAGE'
Usage: scripts/qa/banlist_check.sh [options]

Options:
  --self-test          Run the synthetic positive-control test and exit.
  --root PATH         Repository root. Defaults to the script's repo root.
  --whitelist PATH    Whitelist file. Defaults to scripts/qa/banlist_whitelist.txt.
  --target PATH       Scan only this path. May be repeated.
  --extra-path PATH   Add an extra scan path to the default targets.
  -h, --help          Show this help.

The script prints matches by check family and exits 1 if any unwhitelisted
HIGH-severity match remains.
USAGE
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --self-test)
      SELF_TEST=1
      shift
      ;;
    --root)
      ROOT="$2"
      shift 2
      ;;
    --whitelist)
      WHITELIST="$2"
      shift 2
      ;;
    --target)
      if [ "$TARGET_OVERRIDE" -eq 0 ]; then
        TARGETS=()
        TARGET_OVERRIDE=1
      fi
      TARGETS+=("$2")
      shift 2
      ;;
    --extra-path)
      TARGETS+=("$2")
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "banlist_check.sh: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

ROOT="$(cd "$ROOT" && pwd)"
if [ -z "$WHITELIST" ]; then
  WHITELIST="$ROOT/scripts/qa/banlist_whitelist.txt"
fi

if [ "$TARGET_OVERRIDE" -eq 0 ]; then
  DEFAULT_TARGETS=(
    "$ROOT/crates/lbm-core/src"
    "$ROOT/crates/lbm-cli/examples"
  )
  if [ "${#TARGETS[@]}" -gt 0 ]; then
    TARGETS=("${DEFAULT_TARGETS[@]}" "${TARGETS[@]}")
  else
    TARGETS=("${DEFAULT_TARGETS[@]}")
  fi
fi

require_rg() {
  if ! command -v rg >/dev/null 2>&1; then
    echo "banlist_check.sh: ripgrep (rg) is required" >&2
    exit 2
  fi
}

trim_left() {
  sed 's/^[[:space:]]*//'
}

source_text() {
  local record="$1"
  local rest="${record#*:}"
  printf '%s\n' "${rest#*:}"
}

is_full_line_comment() {
  local text
  text="$(source_text "$1" | trim_left)"
  case "$text" in
    ""|"//"*|"#"*)
      return 0
      ;;
  esac
  return 1
}

is_whitelisted() {
  local record="$1"
  local normalized
  local pattern

  normalized="$(printf '%s\n' "$record" | sed 's/^\(.*\):[0-9][0-9]*:/\1:/')"
  [ -f "$WHITELIST" ] || return 1
  while IFS= read -r pattern || [ -n "$pattern" ]; do
    case "$pattern" in
      ""|"#"*) continue ;;
    esac
    case "$record" in
      *"$pattern"*) return 0 ;;
    esac
    case "$normalized" in
      *"$pattern"*) return 0 ;;
    esac
  done < "$WHITELIST"

  return 1
}

run_rg() {
  local regex="$1"
  shift
  rg -n --no-heading --color never --glob '*.rs' -e "$regex" "$@" 2>/dev/null || true
}

HIGH_FLAGGED_TOTAL=0

run_check() {
  local name="$1"
  local severity="$2"
  local regex="$3"
  shift 3
  local targets=("$@")
  local raw flagged whitelisted record
  local raw_count=0
  local flagged_count=0
  local whitelisted_count=0

  raw="$(mktemp "${TMPDIR:-/tmp}/banlist-raw.XXXXXX")"
  flagged="$(mktemp "${TMPDIR:-/tmp}/banlist-flagged.XXXXXX")"
  whitelisted="$(mktemp "${TMPDIR:-/tmp}/banlist-whitelisted.XXXXXX")"

  run_rg "$regex" "${targets[@]}" > "$raw"

  while IFS= read -r record || [ -n "$record" ]; do
    [ -n "$record" ] || continue
    is_full_line_comment "$record" && continue
    raw_count=$((raw_count + 1))
    if is_whitelisted "$record"; then
      whitelisted_count=$((whitelisted_count + 1))
      printf '%s\n' "$record" >> "$whitelisted"
    else
      flagged_count=$((flagged_count + 1))
      printf '%s\n' "$record" >> "$flagged"
    fi
  done < "$raw"

  printf '\n== %s (%s) ==\n' "$name" "$severity"
  if [ "$flagged_count" -gt 0 ]; then
    sed "s/^/[$severity] /" "$flagged"
  else
    printf 'No unwhitelisted matches.\n'
  fi
  printf 'Count: %d flagged, %d whitelisted, %d raw code matches.\n' \
    "$flagged_count" "$whitelisted_count" "$raw_count"

  if [ "$severity" = "HIGH" ]; then
    HIGH_FLAGGED_TOTAL=$((HIGH_FLAGGED_TOTAL + flagged_count))
  fi

  rm -f "$raw" "$flagged" "$whitelisted"
}

existing_targets() {
  local t
  for t in "${TARGETS[@]}"; do
    [ -e "$t" ] && printf '%s\n' "$t"
  done
}

run_scan() {
  require_rg

  local scan_targets=()
  local t
  while IFS= read -r t; do
    scan_targets+=("$t")
  done < <(existing_targets)

  if [ "${#scan_targets[@]}" -eq 0 ]; then
    echo "banlist_check.sh: no scan targets exist" >&2
    return 2
  fi

  printf 'Ban-list grep sweep\n'
  printf 'Root: %s\n' "$ROOT"
  printf 'Whitelist: %s\n' "$WHITELIST"
  printf 'Targets:\n'
  printf '  %s\n' "${scan_targets[@]}"

  run_check \
    "CASE-IDENTITY" \
    "HIGH" \
    'sample_?name|harshness|case_?id|case_?name|protocol_?name|scenario_?name|preset_?name|fixture_?name|profile_?name|sample_?kind' \
    "${scan_targets[@]}"

  run_check \
    "CLAMPS ON TRANSPORTED QUANTITIES" \
    "HIGH" \
    '\.(clamp|min|max)[[:space:]]*\(' \
    "${scan_targets[@]}"

  run_check \
    "BARE CALIBRATED LITERALS" \
    "HIGH" \
    '(^|[^A-Za-z0-9_])(0\.0025|0\.15|0\.16|2\.5[eE]-5)([^A-Za-z0-9_]|$)' \
    "${scan_targets[@]}"

  run_check \
    "SILENT PHYSICAL DEFAULTS" \
    "HIGH" \
    '\.unwrap_or[[:space:]]*\(' \
    "${scan_targets[@]}"

  printf '\nSummary: %d unwhitelisted HIGH-severity matches.\n' "$HIGH_FLAGGED_TOTAL"
  [ "$HIGH_FLAGGED_TOTAL" -eq 0 ]
}

self_test() {
  local tmpdir fixture empty_whitelist output status
  tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/banlist-self-test.XXXXXX")"
  fixture="$tmpdir/known_bad.rs"
  empty_whitelist="$tmpdir/empty_whitelist.txt"
  output="$tmpdir/output.txt"

  cat > "$fixture" <<'EOF'
pub fn selected_by_case_id(case_id: &str, velocity: f64, rho: f64) -> f64 {
    if case_id == "gentle" {
        return velocity.max(0.16) + rho;
    }
    velocity
}
EOF
  : > "$empty_whitelist"

  "$0" --root "$ROOT" --target "$tmpdir" --whitelist "$empty_whitelist" > "$output" 2>&1
  status=$?

  if [ "$status" -eq 0 ]; then
    echo "banlist_check.sh self-test FAIL: synthetic banned pattern was not rejected" >&2
    cat "$output" >&2
    rm -rf "$tmpdir"
    return 1
  fi
  if ! grep -q 'CASE-IDENTITY' "$output" || ! grep -q 'CLAMPS ON TRANSPORTED QUANTITIES' "$output"; then
    echo "banlist_check.sh self-test FAIL: expected check families were not reported" >&2
    cat "$output" >&2
    rm -rf "$tmpdir"
    return 1
  fi

  rm -rf "$tmpdir"
  echo "banlist_check.sh self-test PASS"
}

if [ "$SELF_TEST" -eq 1 ]; then
  self_test
else
  run_scan
fi
