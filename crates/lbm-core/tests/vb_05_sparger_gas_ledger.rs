// VB-05 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-05--sparger-gas-ledger

const VB05_IGNORE_REASON: &str = "VB-05: waits on BCFD-046/047";

const GAS_LEDGER_RELATIVE_TOLERANCE: f64 = 0.02;
const INITIAL_TRANSIENT_STEPS: usize = 500;
const GAS_VOLUME_FLOW_Q: f64 = 1.0e-6;
const GAS_LEDGER_RUN_TIME: f64 = 20.0;
const MIN_ALLOWED_PHI: f64 = 0.0;

#[derive(Clone, Debug)]
struct GasLedgerSample {
    step: usize,
    time: f64,
    integrated_gas_volume: f64,
}

#[derive(Clone, Debug)]
struct GasLedgerRun {
    volumetric_flow_q: f64,
    samples: Vec<GasLedgerSample>,
    min_phi: f64,
}

#[ignore = "VB-05: waits on BCFD-046/047"]
#[test]
fn ring_sparger_closed_lid_gas_ledger_balances_q_times_t_after_transient() {
    let run = pending_ring_sparger_closed_lid_run();

    assert_gas_ledger_balances_expected_q_times_t(&run);
}

#[ignore = "VB-05: waits on BCFD-046/047"]
#[test]
fn ring_sparger_rejects_negative_phi() {
    let run = pending_ring_sparger_closed_lid_run();

    assert_no_negative_phi(&run);
}

#[ignore = "VB-05: waits on BCFD-046/047"]
#[test]
fn liquid_injection_through_gas_sparger_is_rejected() {
    let rejection = pending_liquid_injection_rejection();

    assert_liquid_injection_rejected(rejection);
}

fn pending_ring_sparger_closed_lid_run() -> GasLedgerRun {
    panic!(
        "{VB05_IGNORE_REASON}: run real ring-sparger resolved phase-field case at \
         Q={GAS_VOLUME_FLOW_Q}, run_time={GAS_LEDGER_RUN_TIME}; no mocked solver data"
    )
}

fn pending_liquid_injection_rejection() -> Result<(), String> {
    panic!("{VB05_IGNORE_REASON}: validate real sparger scenario rejects inlet_phase=liquid")
}

fn assert_gas_ledger_balances_expected_q_times_t(run: &GasLedgerRun) {
    assert!(
        run.volumetric_flow_q > MIN_ALLOWED_PHI,
        "VB-05 gas volumetric flow must be positive; Q={}",
        run.volumetric_flow_q
    );
    for sample in run
        .samples
        .iter()
        .filter(|sample| sample.step >= INITIAL_TRANSIENT_STEPS)
    {
        let expected_volume = run.volumetric_flow_q * sample.time;
        let relative_error =
            (sample.integrated_gas_volume - expected_volume).abs() / expected_volume.abs();
        assert!(
            relative_error <= GAS_LEDGER_RELATIVE_TOLERANCE,
            "VB-05 gas ledger balance: step={}, time={}, measured_volume={}, expected_Qt={}, \
             relative_error={relative_error}, tolerance={GAS_LEDGER_RELATIVE_TOLERANCE}; \
             denominator is Q*t",
            sample.step,
            sample.time,
            sample.integrated_gas_volume,
            expected_volume
        );
    }
}

fn assert_no_negative_phi(run: &GasLedgerRun) {
    assert!(
        run.min_phi >= MIN_ALLOWED_PHI,
        "VB-05 phase field must reject negative phi: min_phi={}, allowed_min={MIN_ALLOWED_PHI}",
        run.min_phi
    );
}

fn assert_liquid_injection_rejected(rejection: Result<(), String>) {
    assert!(
        rejection.is_err(),
        "VB-05 liquid injection through a gas sparger must be rejected with a structured error"
    );
}
