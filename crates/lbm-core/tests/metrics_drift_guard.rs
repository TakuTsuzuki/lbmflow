//! Drift guard: pins scripts/qa/metrics.py to crates/lbm-core/tests/common/metrics.rs
//!
//! Runs `python3 scripts/qa/metrics.py --drift-guard` on a fixed input vector
//! (SAME vector as `_drift_guard()` in metrics.py) and asserts every
//! returned metric agrees with the Rust computation to 1e-12. Any change
//! to metrics.rs that is not matched in metrics.py fails this test.

mod common;

use common::metrics::*;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::process::Command;

fn compute_rust() -> HashMap<&'static str, f64> {
    let reference = [1.0f64, 2.0, 3.0, 4.0, 5.0];
    let actual = [1.05f64, 1.9, 3.2, 4.1, 4.9];
    let h = [0.1f64, 0.05, 0.025, 0.0125];
    let err = [4e-2f64, 1.05e-2, 2.7e-3, 7e-4];
    let y: [f64; 5] = [0.5, 1.0, 2.0, 3.5, 5.0];
    let amp: Vec<f64> = y.iter().map(|&v| 0.7 * (-0.35 * v).exp()).collect();
    let omega = 2.0 * PI / 100.0;
    let t: Vec<f64> = (0..400).map(|i| i as f64).collect();
    let sig: Vec<f64> = t.iter().map(|&ti| 0.02 * (omega * ti + 0.6).sin()).collect();

    let of = order_fit(&h, &err);
    let ef = envelope_fit(&y, &amp);
    let (amp_fit, phase_fit_v) = phase_fit(&t, &sig, omega);
    let samples = [(1.0, 1.0), (2.0, 4.0), (3.0, 9.9), (4.0, 16.0)];
    let cg = curve_agreement(|x| x * x, &samples, 0.05, 0.0);

    HashMap::from([
        ("l2_rel", l2_rel(&actual, &reference)),
        ("linf_rel", linf_rel(&actual, &reference, 0.0)),
        ("order_slope", of.slope),
        ("order_r2", of.r2),
        ("env_slope", ef.slope),
        ("env_intercept", ef.intercept),
        ("phase_amp", amp_fit),
        ("phase_phase", phase_fit_v),
        ("monotonicity", monotonicity(&[5.0, 4.0, 3.0, 3.5, 1.0])),
        ("curve_worst_x", cg.worst_x),
        ("curve_max_rel_dev", cg.max_rel_dev),
    ])
}

#[test]
fn python_mirror_matches_rust_to_1e_minus_12() {
    let script = std::env::current_dir()
        .unwrap()
        .ancestors()
        .find(|p| p.join("scripts/qa/metrics.py").exists())
        .expect("scripts/qa/metrics.py not found from cwd")
        .join("scripts/qa/metrics.py");
    let out = match Command::new("python3")
        .arg(&script)
        .arg("--drift-guard")
        .output()
    {
        Ok(o) => o,
        Err(e) => {
            eprintln!("skipping metrics_drift_guard: python3 not runnable: {e}");
            return; // graceful skip when python3 is unavailable
        }
    };
    assert!(
        out.status.success(),
        "python3 metrics.py --drift-guard failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let mut py: HashMap<&str, f64> = HashMap::new();
    let rust = compute_rust();
    let rust_keys: Vec<&str> = rust.keys().copied().collect();
    for line in stdout.lines() {
        let (k, v) = line.split_once('=').unwrap_or_else(|| panic!("bad line: {line}"));
        let key = rust_keys
            .iter()
            .copied()
            .find(|&r| r == k)
            .unwrap_or_else(|| panic!("unexpected python key: {k}"));
        py.insert(key, v.trim().parse::<f64>().unwrap());
    }
    assert_eq!(py.len(), rust.len(), "python emitted {} keys, rust has {}", py.len(), rust.len());
    for (k, rv) in &rust {
        let pv = py[k];
        let denom = rv.abs().max(1.0);
        let rel = (pv - rv).abs() / denom;
        assert!(
            rel <= 1e-12,
            "drift on {k}: rust={rv:.17e}, python={pv:.17e}, rel={rel:e}"
        );
    }
}
