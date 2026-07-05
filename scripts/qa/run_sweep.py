#!/usr/bin/env python3
"""Physics-QA sweep driver: RUN -> COLLECT -> DETECT, one config at a time.

Usage:
  python3 scripts/qa/run_sweep.py --bin target/release/lbm --out out/qa-pass1
      [--only id1,id2] [--validate-only]

Writes per-config outputs under <out>/<id>/ (scenario.json, validate.json,
manifest.json plus field files) and a machine-readable <out>/results.json.
Exit code = number of failed checks (validate-only: number of build errors).
"""

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import matrix  # noqa: E402
import qa_checks as qc  # noqa: E402

# severity proposal on failure; final disposition is curated in the log
SEVERITY = {
    "finite/status": "S1",
    "|u| hard ceiling": "S1",
    "|u| scale sanity": "S1",
    "mass drift": "S1",
    "momentum growth": "S1",
    "field uniformity": "S1",
    "warnings audit": "S2",
    "Karman saturation guard": "S3",
}
DEFAULT_SEVERITY = "S2"  # band/accuracy misses


def sev_for(check_name):
    for key, s in SEVERITY.items():
        if check_name.startswith(key.split(" @")[0]):
            return s
    return DEFAULT_SEVERITY


def run_config(cfg, bin_path, out_root, validate_only):
    cdir = out_root / cfg["id"]
    cdir.mkdir(parents=True, exist_ok=True)
    sc_path = cdir / "scenario.json"
    sc_path.write_text(json.dumps(cfg["scenario"], indent=2) + "\n")

    val = subprocess.run([bin_path, "validate", str(sc_path)],
                         capture_output=True, text=True, timeout=300)
    validate_report = json.loads(val.stdout) if val.stdout.strip() else {
        "ok": False, "error": val.stderr}
    (cdir / "validate.json").write_text(json.dumps(validate_report, indent=2) + "\n")

    result = {"id": cfg["id"], "track": cfg.get("track"),
              "validate": validate_report, "findings": [], "manifest": None}
    if cfg.get("note"):
        result["note"] = cfg["note"]
    if validate_only or not validate_report.get("ok", False):
        return result

    proc = subprocess.run(
        [bin_path, "run", str(sc_path), "--out", str(cdir), "--json"],
        capture_output=True, text=True, timeout=3600)
    if proc.returncode != 0:
        result["findings"].append({
            "check": "runner exit", "ok": False, "severity": "S1",
            "expected": "exit 0", "observed": f"exit {proc.returncode}",
            "detail": (proc.stderr or proc.stdout)[-2000:]})
        return result
    manifest = json.loads(proc.stdout)
    result["manifest"] = {k: manifest[k] for k in
                          ("status", "stepsRun", "wallSeconds", "mlups",
                           "diagnostics", "warnings")}

    findings = [qc.check_finite_and_status(cdir, cfg),
                qc.check_speed_ceiling(cdir, cfg)]

    # warnings audit: expected advisories present? unexpected ones?
    got = sorted({w["field"] for w in manifest.get("warnings", [])})
    want = sorted(set(cfg.get("expect_warnings", [])))
    missing = [w for w in want if w not in got]
    extra = [w for w in got if w not in want]
    findings.append({
        "check": "warnings audit", "ok": not missing,
        "expected": f"validator advisories on {want or 'none'}",
        "observed": f"got {got or 'none'}",
        "detail": (f"missing={missing} " if missing else "")
                  + (f"unexpected={extra} (recorded, not a failure)" if extra else "")})

    diverged = manifest["status"] == "diverged"
    for chk in cfg["checks"]:
        if diverged:
            findings.append({"check": chk["name"], "ok": False,
                             "expected": chk.get("source", ""),
                             "observed": "skipped: run diverged", "detail": ""})
            continue
        fn = qc.CHECKS[chk["name"]]
        try:
            f = fn(cdir, cfg, chk.get("args"))
        except Exception as e:  # collection-surface failure is itself a finding
            f = {"check": chk["name"], "ok": False,
                 "expected": "check computable from collected outputs",
                 "observed": f"{type(e).__name__}: {e}", "detail": "harness/collection"}
        fs = f if isinstance(f, list) else [f]
        for one in fs:
            one["source"] = chk.get("source", "")
        findings.extend(fs)

    for f in findings:
        if not f["ok"] and "severity" not in f:
            f["severity"] = sev_for(f["check"])
    result["findings"] = findings
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default="target/release/lbm")
    ap.add_argument("--out", default="out/qa-pass1")
    ap.add_argument("--only", default="")
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()

    only = {s for s in args.only.split(",") if s}
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    results = []
    for cfg in matrix.CONFIGS:
        if only and cfg["id"] not in only:
            continue
        label = f"[{cfg['id']}]"
        print(f"{label} running...", flush=True)
        try:
            r = run_config(cfg, args.bin, out_root, args.validate_only)
        except subprocess.TimeoutExpired:
            r = {"id": cfg["id"], "findings": [{
                "check": "runner timeout", "ok": False, "severity": "S1",
                "expected": "completes within timeout", "observed": "timeout"}]}
        results.append(r)
        n_fail = sum(1 for f in r.get("findings", []) if not f["ok"])
        man = r.get("manifest")
        stat = man["status"] if man else "n/a"
        print(f"{label} status={stat} failed_checks={n_fail}", flush=True)

    # cross-config: BGK convergence order pairs
    by_id = {r["id"]: r for r in results}
    for cfg in matrix.CONFIGS:
        pid = cfg.get("pair_order_with")
        if not pid or cfg["id"] not in by_id or pid not in by_id:
            continue
        try:
            coarse = next(f["l2rel"] for f in by_id[pid]["findings"] if "l2rel" in f)
            fine = next(f["l2rel"] for f in by_id[cfg["id"]]["findings"] if "l2rel" in f)
            order = math.log2(coarse / fine)
            by_id[cfg["id"]]["findings"].append({
                "check": "BGK convergence order", "ok": order >= 1.7,
                "expected": "order >= 1.7 for H=8->16 (VALIDATION T2)",
                "observed": f"order = {order:.3f} (L2rel {coarse:.3g} -> {fine:.3g})",
                "severity": None if order >= 1.7 else "S2",
                "source": "VALIDATION T2"})
        except StopIteration:
            pass

    (out_root / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    failures = [(r["id"], f) for r in results
                for f in r.get("findings", []) if not f["ok"]]
    print(f"\n=== sweep done: {len(results)} configs, {len(failures)} failed checks ===")
    for cid, f in failures:
        print(f"  FAIL [{cid}] {f['check']}: expected {f['expected']}; "
              f"observed {f['observed']}")
    return len(failures)


if __name__ == "__main__":
    sys.exit(main())
