// VB-01 adversarial validation skeleton.
// Source of truth: docs/VALIDATION_BIOPROCESS.md#vb-01--single-phase-stirred-tank-np

const VB01_IGNORE_REASON: &str = "VB-01: waits on BCFD-030/031";

const RUSHTON_REFERENCE_NP: f64 = 5.0;
const PBT45_REFERENCE_NP: f64 = 1.3;
const TURBULENT_RE_MIN: f64 = 10_000.0;
const NP_RELATIVE_TOLERANCE: f64 = 0.15;
const FINE_GRID_RELATIVE_DIFFERENCE_MAX: f64 = 0.05;
const GCI_MIN_OBSERVED_ORDER: f64 = 1.0;
const THREE_GRID_RESOLUTIONS: [usize; 3] = [64, 96, 128];

#[derive(Clone, Copy, Debug)]
enum ImpellerGeometry {
    Rushton,
    Pbt45,
}

#[derive(Clone, Copy, Debug)]
struct NpOperatingPoint {
    geometry: ImpellerGeometry,
    reynolds: f64,
    measured_np: f64,
}

#[derive(Clone, Copy, Debug)]
struct GridNp {
    resolution: usize,
    np: f64,
}

#[ignore = "VB-01: waits on BCFD-030/031"]
#[test]
fn rushton_np_matches_published_correlation() {
    let point = pending_np_operating_point(ImpellerGeometry::Rushton);

    assert_turbulent_reference_regime(point.reynolds);
    assert_np_within_15_percent_of_published(point.measured_np, point.geometry);
}

#[ignore = "VB-01: waits on BCFD-030/031"]
#[test]
fn pbt45_np_matches_published_correlation() {
    let point = pending_np_operating_point(ImpellerGeometry::Pbt45);

    assert_turbulent_reference_regime(point.reynolds);
    assert_np_within_15_percent_of_published(point.measured_np, point.geometry);
}

#[ignore = "VB-01: waits on BCFD-030/031"]
#[test]
fn np_three_grid_convergence_between_two_finest_grids() {
    let grid_np = pending_three_grid_np_series();

    assert_expected_three_grid_resolutions(&grid_np);
    assert_three_grid_np_convergence(&grid_np);
}

fn pending_np_operating_point(geometry: ImpellerGeometry) -> NpOperatingPoint {
    panic!(
        "{VB01_IGNORE_REASON}: run the real 3D stirred-tank solver and return \
         statistically stationary Np for {geometry:?}; no mocked fluid solver data"
    )
}

fn pending_three_grid_np_series() -> [GridNp; 3] {
    panic!(
        "{VB01_IGNORE_REASON}: run the real stirred-tank solver at \
         {THREE_GRID_RESOLUTIONS:?} and return time-averaged Np"
    )
}

fn assert_turbulent_reference_regime(reynolds: f64) {
    assert!(
        reynolds > TURBULENT_RE_MIN,
        "VB-01 reference correlations require Re > {TURBULENT_RE_MIN}; measured Re={reynolds}"
    );
}

fn assert_np_within_15_percent_of_published(measured_np: f64, geometry: ImpellerGeometry) {
    let reference_np = published_reference_np(geometry);
    let relative_error = relative_error(measured_np, reference_np);
    assert!(
        relative_error <= NP_RELATIVE_TOLERANCE,
        "VB-01 Np for {geometry:?}: measured={measured_np}, reference={reference_np}, \
         relative_error={relative_error}, tolerance={NP_RELATIVE_TOLERANCE}; \
         denominator is published reference Np"
    );
}

fn assert_expected_three_grid_resolutions(grid_np: &[GridNp; 3]) {
    let actual = [
        grid_np[0].resolution,
        grid_np[1].resolution,
        grid_np[2].resolution,
    ];
    assert_eq!(
        actual, THREE_GRID_RESOLUTIONS,
        "VB-01 requires the three-grid sequence {THREE_GRID_RESOLUTIONS:?}"
    );
}

fn assert_three_grid_np_convergence(grid_np: &[GridNp; 3]) {
    let coarse = grid_np[0].np;
    let medium = grid_np[1].np;
    let fine = grid_np[2].np;
    let fine_pair_relative_difference = relative_error(fine, medium);
    let observed_order = observed_convergence_order(coarse, medium, fine);

    assert!(
        fine_pair_relative_difference <= FINE_GRID_RELATIVE_DIFFERENCE_MAX,
        "VB-01 Np fine-grid convergence: Np_medium={medium}, Np_fine={fine}, \
         relative_difference={fine_pair_relative_difference}, \
         tolerance={FINE_GRID_RELATIVE_DIFFERENCE_MAX}; denominator is medium-grid Np"
    );
    assert!(
        observed_order >= GCI_MIN_OBSERVED_ORDER,
        "VB-01 Np three-grid GCI/order check: observed_order={observed_order}, \
         required_order_floor={GCI_MIN_OBSERVED_ORDER}, values=({coarse}, {medium}, {fine})"
    );
}

fn published_reference_np(geometry: ImpellerGeometry) -> f64 {
    match geometry {
        ImpellerGeometry::Rushton => RUSHTON_REFERENCE_NP,
        ImpellerGeometry::Pbt45 => PBT45_REFERENCE_NP,
    }
}

fn relative_error(measured: f64, reference: f64) -> f64 {
    (measured - reference).abs() / reference.abs()
}

fn observed_convergence_order(coarse: f64, medium: f64, fine: f64) -> f64 {
    let numerator = (coarse - medium).abs();
    let denominator = (medium - fine).abs();
    (numerator / denominator).log2()
}
