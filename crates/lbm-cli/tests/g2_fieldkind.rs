use serde_json::json;
use std::fs;
use std::path::Path;
use std::process::Command;

fn read_csv_field(path: &Path) -> Vec<f64> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    text.lines()
        .skip(1)
        .flat_map(|line| line.split(','))
        .map(|v| {
            v.parse::<f64>()
                .unwrap_or_else(|e| panic!("parse {v:?} in {}: {e}", path.display()))
        })
        .collect()
}

#[test]
fn g2_runner_fieldkind_dissipation_is_nu_times_shear_rate_squared() {
    let root = std::env::temp_dir().join(format!(
        "lbm_g2_fieldkind_{}_{}",
        std::process::id(),
        "dissipation"
    ));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let scenario_path = root.join("scenario.json");
    let out_dir = root.join("out");
    let nu = 0.1f64;
    let scenario = json!({
        "name": "g2-fieldkind",
        "grid": { "nx": 4, "ny": 10 },
        "physics": {
            "nu": nu,
            "collision": { "type": "trt" },
            "force": [1e-6, 0.0],
            "precision": "f64"
        },
        "edges": {
            "left": { "type": "periodic" },
            "right": { "type": "periodic" },
            "bottom": { "type": "bounceBack" },
            "top": { "type": "bounceBack" }
        },
        "run": { "steps": 2 },
        "outputs": [
            { "field": "shearRate", "format": "csv", "every": 0 },
            { "field": "dissipationRate", "format": "csv", "every": 0 }
        ]
    });
    fs::write(
        &scenario_path,
        serde_json::to_string_pretty(&scenario).unwrap(),
    )
    .unwrap();

    let status = Command::new(env!("CARGO_BIN_EXE_lbm"))
        .arg("run")
        .arg(&scenario_path)
        .arg("--out")
        .arg(&out_dir)
        .arg("--json")
        .status()
        .expect("run lbm scenario");
    assert!(status.success(), "lbm run failed with status {status}");

    let shear = read_csv_field(&out_dir.join("shearrate_2.csv"));
    let dissipation = read_csv_field(&out_dir.join("dissipationrate_2.csv"));
    assert_eq!(shear.len(), dissipation.len());
    let mut max_abs = 0.0f64;
    for (i, (gamma, eps)) in shear.iter().zip(&dissipation).enumerate() {
        assert!(
            gamma.is_finite() && eps.is_finite(),
            "G2 FieldKind output must be finite at flat cell {i}: gamma={gamma:e}, eps={eps:e}"
        );
        let expected = nu * gamma * gamma;
        max_abs = max_abs.max((eps - expected).abs());
    }
    eprintln!("G2 FieldKind DissipationRate consistency max_abs={max_abs:.3e}");
    assert!(
        max_abs <= 1.0e-18,
        "G2 FieldKind DissipationRate != nu*ShearRate^2 pointwise: max_abs={max_abs:.3e}"
    );
}
