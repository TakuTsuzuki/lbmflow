# E1 Ablation A/B: Velocity-Correction Footprint Measurement

Date: 2026-07-07

## Setup

Temporary scratch test: `crates/lbm-core/tests/scratch_e1_ablation.rs` (deleted after measurement).

Command used for both passes:

```bash
cargo test -p lbm-core --release --test scratch_e1_ablation -- --nocapture
```

Measurement: D3Q19 CentralMoment, `N=32`, `nu in {0.02, 0.10}`,
`u0 in {0.02, 0.04, 0.08}`. The fitted quantity is

```text
c(nu) = d(nu_eff / nu - 1) / d(u0^2)
```

using the same 3D TGV decay-rate and `linear_fit` machinery as
`crates/lbm-core/tests/accuracy_audit_cumulant.rs`.

Visual artifact: `/private/tmp/e1-ablation-footprint.png`.

## Prediction

The landed velocity correction applies the local modifier

```text
domega / omega = -0.16 |u|^2
```

to the central-moment shear relaxation. With

```text
nu = (1 / omega - 1/2) / 3
```

a perturbation gives

```text
dnu / nu = -(domega / omega) * 2 / (2 - omega)
         =  0.16 |u|^2 * 2 / (2 - omega)
```

For the classic 3D TGV used here,

```text
u = u0 (sin x cos y cos z, -cos x sin y cos z, 0)
```

so

```text
<|u|^2> / u0^2
  = <cos^2 z> * (<sin^2 x cos^2 y> + <cos^2 x sin^2 y>)
  = (1/2) * (1/4 + 1/4)
  = 1/4.
```

Therefore the expected slope footprint is

```text
delta_c(nu) = 0.16 * W * 2 / (2 - omega), W = 1/4,
omega = 1 / (3 nu + 0.5).
```

## Results

| nu | omega | c ON | c OFF | ON-OFF | predicted delta_c | rel. error |
|---:|---:|---:|---:|---:|---:|---:|
| 0.02 | 1.7857142857 | 8.3345195863 | 8.0602422387 | 0.2742773476 | 0.3733333333 | 26.53% |
| 0.10 | 1.2500000000 | -0.1289900434 | -0.2091509745 | 0.0801609312 | 0.1066666667 | 24.85% |

Raw ON output:

```text
SCRATCH_E1_ABLATION N=32 nu=2e-2 omega=1.7857142857142856e0 u0=[0.02, 0.04, 0.08] defects=[0.02582456386848686, 0.03566806315725213, 0.07577903998125324] c=8.334519586281777e0 intercept=2.2420567494075092e-2 r2=9.99990765656135e-1
SCRATCH_E1_ABLATION N=32 nu=1e-1 omega=1.25e0 u0=[0.02, 0.04, 0.08] defects=[-0.00023147397557476967, -0.00038828387979628065, -0.0010060881864858429] c=-1.289900433637922e-1 intercept=-1.807765592003462e-4 r2=9.999936808680314e-1
```

Raw OFF output:

```text
SCRATCH_E1_ABLATION N=32 nu=2e-2 omega=1.7857142857142856e0 u0=[0.02, 0.04, 0.08] defects=[0.02571303857704521, 0.03522361669472107, 0.07402058781965515] c=8.06024223869964e0 intercept=2.241706942878148e-2 r2=9.999896471335668e-1
SCRATCH_E1_ABLATION N=32 nu=1e-1 omega=1.25e0 u0=[0.02, 0.04, 0.08] defects=[-0.0002637615563282347, -0.0005172741829633232, -0.0015195112226491503] c=-2.09150974543313e-1 intercept=-1.812262585922929e-4 r2=9.999962321509573e-1
```

## Verdict

The ON-OFF difference matches the derived `delta_c(nu)` within 30 percent at
both viscosities. Under the requested rule, the `-0.16 |u|^2` term is
classified as a real Galilean-class correction: candidate B.

## Behavior-Validity Review

Pattern: enabling the local velocity correction increases the fitted
`u0^2` viscosity-defect slope at both viscosities, with a larger shift at
lower `nu`.

Mechanism: the local `-0.16 |u|^2` relaxation-rate perturbation maps through
`nu(omega)` into a positive effective-viscosity footprint weighted by the TGV
spatial mean `<|u|^2> / u0^2 = 1/4`.

Resolved vs closure: the TGV decay-rate measurement is resolved by the LBM
operator; the tested non-resolved term is the central-moment velocity-dependent
shear-rate modifier.

Artifacts checked: periodic domain, no walls, no outlets, no clamps, no seams
(`LocalPeriodic` with `[1, 1, 1]`). The linear fits are near-perfect
(`r2 > 0.99998`) in all four rows, so the slope comparison is not dominated by
fit scatter.

Verdict: PHYSICAL for the measured footprint. The term remains tied to the
current CentralMoment TGV evidence; broader Galilean invariance claims still
depend on the separate holdout tests.

Routing: none.
