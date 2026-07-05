# T15.5 — 3D Lid-Driven Cavity (Cube) Reference Data Collection Record

- Created: 2026-07-05 (validation data curation, Claude)
- Purpose: A record of collecting the reference tables for VALIDATION.md **T15.5**
  (D3Q19 3D cavity validation), **with citations**, and performing an
  **ingestion-time smoothness check** (lesson from the Ghia typo incident, see
  PHYSICS.md 2026-07-05). codex should freeze the test constant tables based on
  this document.
- Every number in this document has been double-checked via "transcribe → tabulate
  → re-cross-check against the original text." No number without a citation is
  included. The smoothness-check arithmetic is fully disclosed in §4, so codex can
  recompute and cross-check it by script at ingestion time.

---

## 0. Summary (what could be frozen, what could not)

| Item | Result |
|---|---|
| **Re=1000 cube, symmetry-plane centerline profiles (u/U 17 points, w/U 17 points)** | **Collected** (redistributed data from the Albensoeder & Kuhlmann 2005 table. §3.1–3.2) |
| Cross-comparison of Re=1000 extrema (u_min, w_min, w_max and locations) across multiple sources | **Collected** (5 sources. §3.3, §5) |
| **Re=100 / Re=400 profile numeric tables** | **Numeric values not obtained**. The primary candidates (Ku et al. 1987, Wong & Baker 2002, Yang et al. 1998, Iwatsu et al. 1989) are all behind paywalls, and no public reprint or redistribution could be found either (§7). |
| Smoothness check | Performed on all 30 interior points of the 2 candidate frozen tables, **no points of concern** (§4) |

**Recommendation**: Freeze the first round of T15.5 with **Re=1000 only** (the
reference quality is highest here, and it aligns with the same Ghia-style
17-point format as 2D T7). Hold Re=100/400 until the actions in §7 (obtaining
the Ku et al. / Wong & Baker originals via a library) are completed, and note
explicitly in the test spec that "Re=100/400 reference is not yet frozen."

---

## 1. Problem Definition and Coordinate Convention (canonical form used in this document)

- Domain: unit cube [0,1]³. x: lid-driving direction, y: spanwise direction, z: vertical (bottom z=0, lid z=1).
- Boundary conditions: no-slip walls on all 6 faces. Only the top face z=1 moves with tangential velocity (U, 0, 0) (the spanwise faces y=0, 1 are also **rigid stationary walls**. Not spanwise-periodic).
- Re = U·L/ν (L = cavity edge length).
- Symmetry plane: y = 0.5 (the steady solution is mirror-symmetric about this plane for Re ≲ 1900. See §6.5).
- Centerlines:
  - **Vertical centerline**: **u/U** (lid-parallel component) at (x, y) = (0.5, 0.5), z ∈ [0,1]
  - **Horizontal centerline**: **w/U** (vertical component) at (y, z) = (0.5, 0.5), x ∈ [0,1]
- The table extraction coordinates are the same **17 points from Ghia et al. (1982)** used in 2D T7
  (the u-line uses Ghia's y-coordinates, the w-line uses Ghia's x-coordinates, identical values).

---

## 2. Source Survey Results (accessibility and presence/absence of numeric tables)

| Reference | Type | Access | Numeric table | Use in this document |
|---|---|---|---|---|
| Albensoeder & Kuhlmann (2005) JCP 206:536–558 | Chebyshev collocation (spectral) | Original behind paywall. **Numeric values obtained via 2 redistribution routes** | Yes (Table 5/6: Re=1000 centerline, Ghia-style 17 points) | **Primary reference** (§3.1–3.2) |
| Ben Beya & Lili (2008) C.R. Mécanique 336:863–872 | FVM 64³ multigrid | **Fully read from public PDF (academy website)** | Yes (Table 1: 6-source comparison of Re=1000 extrema) | Independent cross-check of extrema and A&K values (§3.3) |
| Zhang, Shi & Wang (2015) arXiv:1503.03337 | LBM D3Q14/18 MRT | **Fully read (public, arXiv)** | No profile table (numeric values for the u_min grid-refinement series are in the body text) | Reference point for LBM's grid dependence (§5) |
| Ku, Hirsh & Taylor (1987) JCP 70:439–462 | Pseudospectral (original benchmark) | Paywall (Elsevier) | Unknown (inaccessible) | **Numeric values not obtained** (§7) |
| Jiang, Lin & Povinelli (1993/1994) NASA TM-105904 | Least-squares FEM 50³ | **Fully read from public PDF (NTRS)** | **None (Fig.4–6 figures only)** | Confirmed figures only (§7) |
| Kato, Kawai & Tanahashi (1990) JSME Int. J. II 33(4):649–658 | GSMAC-FEM 40×40×20 (half domain) | **Fully read from public PDF (J-STAGE)** | **None (Fig.3/7/13 only)**. Contains a description of Ku's grid | Confirmed figures only, supplementary bibliographic info for Ku (§7) |
| Wong & Baker (2002) IJNMF 38:99–123 | FEM 48³ | Paywall (Wiley) | Unknown | **Numeric values not obtained** (§7) |
| Yang, Yang, Chen & Hsu (1998) JCP 146:464–487 | Implicit WENO | Paywall | Unknown (used as a reference point by Žunič et al.) | **Numeric values not obtained** (§7) |
| Lo, Murugesan & Young (2005) IJNMF 47:1469–1487 | Velocity–vorticity FDM 101³ | Paywall | Yes (extrema only, obtained via Ben Beya) | Extrema comparison (§3.3) |
| Ding, Shu, Yeo & Xu (2006) CMAME 195:516–533 | Local MQ-DQ 49³ | Paywall | Yes (extrema only, obtained via Ben Beya) | Extrema comparison (§3.3) |
| Žunič, Hriberšek, Škerget & Ravnik (2006, ECCOMAS CFD) | BEM-FEM | **Fully read from public PDF (ljll.fr mirror)** | None (figures only) | Status check only |
| Iwatsu et al. (1989) Fluid Dyn. Res. 5:173–189 | FDM 81³ stretched grid | Paywall | Unknown | **Numeric values not obtained** (§7) |
| Kuhlmann & Romanò (2019) chapter "The lid-driven cavity" | Review | **Public PDF (TU Wien)** | No profile table (within scope checked) | Background/stability context |
| Feldman & Gelfgat (2010) Phys. Fluids 22:093602 / Liberzon et al. arXiv:1107.0449 | FVM / PIV experiment | Fully read from arXiv version | No profile table | Basis for Re_cr ≈ 1914 (§6.5) |

---

## 3. Candidate Frozen Data (Re = 1000, cube, all walls no-slip)

### 3.1 Origin and Provenance

**Primary source**: S. Albensoeder, H.C. Kuhlmann,
"Accurate three-dimensional lid-driven cavity flow",
*Journal of Computational Physics* **206** (2005) 536–558.
doi:10.1016/j.jcp.2004.12.024
(Chebyshev collocation method. The grid that produced the extrema was
96×96×64 — per the citation in Ben Beya & Lili (2008) Table 1. The original's
Table 5 = centerline v distribution, Table 6 = centerline u distribution,
both tabulated at Ghia's 17-point coordinates)

**Route by which the numeric values were obtained (redistribution repository)**: GitHub `tum-pbs/PICT`
(TUM Physics-based Simulation group),
in `tests/validations.py`'s `lid_driven_cavity_3D()`,
arrays `reference_coords_x_T5` / `reference_vel_v` / `reference_coords_y_T6` /
`reference_vel_u`, key `(1000, 1, 1, False)` (= Re 1000, aspect ratio 1×1,
spanwise end walls are rigid walls).

- URL: <https://github.com/tum-pbs/PICT> —
  raw: <https://raw.githubusercontent.com/tum-pbs/PICT/main/tests/validations.py>
- Retrieval date: 2026-07-05. Last commit to the relevant file:
  `a95d7f9d0713262a1bff2bd9e2be5a203ee69208` (2025-05-22, "added PICT")
- In-file comment: "reference data from paper "Accurate three-dimensional
  lid-driven cavity flow" (Albensoeder, Kuhlmann; 2004)" "normalized with /Re"
  (A&K print their table with velocity non-dimensionalized by ν/L, so
  converting to u/U requires dividing by Re. The PICT array holds the
  already-divided values, i.e., in u/U units)
- Cross-checking the same file's 2D Ghia table (Re=100/1000/5000/10000) against
  the canonical Ghia et al. (1982) table showed **agreement at every point**,
  so the transcription quality of this redistributor is judged to be good.

**Independent cross-check (secondary source)**: B. Ben Beya, T. Lili, "Three-dimensional incompressible
flow in a two-sided non-facing lid-driven cubical cavity",
*C. R. Mécanique* **336** (2008) 863–872. doi:10.1016/j.crme.2008.10.004
Public PDF: <https://comptes-rendus.academie-sciences.fr/mecanique/item/10.1016/j.crme.2008.10.004.pdf>
(Table 1 reproduces A&K's extrema to 7 significant digits. Used for the
cross-check in §3.3. Additionally, the same paper's Fig. 2 plots A&K's table
points as symbols, which agree at the figure level with the profiles below)

**Coordinate/sign transformation (PICT array → this document's standard form §1)**:
A&K place the lid at the plane coordinate −0.5 and print their table over
[−0.5, 0.5]³. A single orthogonal-frame mapping moves it to this document's
standard form:

- Vertical centerline (Table 5): `z = 0.5 − x_AK`, **u/U = +(table value)**
  (check: at the table's endpoint x_AK=−0.5 → z=1, value 1.00000 = lid speed ✓; x_AK=+0.5 → z=0, value 0 ✓)
- Horizontal centerline (Table 6): `x = y_AK + 0.5`, **w/U = −(table value)**
  (the sign reflects the difference in lid-drive orientation. Check: after
  transformation w_min ≈ −0.434 @ x=0.9063 = downward flow on the downstream-wall
  side, w_max ≈ +0.244 @ x=0.0938 = upward flow on the upstream-wall side, which
  matches the A&K extrema in Ben Beya Table 1 (w_min −0.4350186 @ x=0.90957,
  w_max +0.2466511 @ x=0.10913) in position, sign, and magnitude ✓)

**Note on precision**: the PICT array has 5 decimal digits. The A&K original has
at least 7 digits (e.g., Ben Beya's citation of −0.2803833), so this table
should be treated as "the A&K table rounded to 5 digits (or /Re-normalized at
5 digits)," with a transcription resolution of ±1×10⁻⁵ assumed.
**A cross-check against the original Table 5/6 is recommended before final
freezing** (§7 action A-1).

### 3.2 Numeric Tables (candidate frozen data)

**Table A: u/U — vertical centerline (x=0.5, y=0.5), Re=1000, cube (rigid spanwise end walls)**
Source: Albensoeder & Kuhlmann (2005) Table 5 (96×96×64), via PICT (§3.1).
z is the height from the bottom. The endpoints z=0, 1 are exactly the boundary condition values.

| # | z/L | u/U |
|---|--------|----------|
| 0 | 0.0000 | 0.00000 |
| 1 | 0.0547 | -0.20623 |
| 2 | 0.0625 | -0.22283 |
| 3 | 0.0703 | -0.23696 |
| 4 | 0.1016 | -0.27293 |
| 5 | 0.1719 | -0.25160 |
| 6 | 0.2813 | -0.10999 |
| 7 | 0.4531 | -0.00612 |
| 8 | 0.5000 | 0.00802 |
| 9 | 0.6172 | 0.03905 |
| 10 | 0.7344 | 0.07334 |
| 11 | 0.8516 | 0.12183 |
| 12 | 0.9531 | 0.33171 |
| 13 | 0.9609 | 0.39821 |
| 14 | 0.9688 | 0.48443 |
| 15 | 0.9766 | 0.58964 |
| 16 | 1.0000 | 1.00000 |

**Table B: w/U — horizontal centerline (y=0.5, z=0.5), Re=1000, cube (rigid spanwise end walls)**
Source: same as above, A&K Table 6, via PICT, with sign transformation w/U = −(table value) (§3.1).
x=0 is the upstream wall (the side the lid moves away from), x=1 is the downstream wall.

| # | x/L | w/U |
|---|--------|----------|
| 0 | 0.0000 | 0.00000 |
| 1 | 0.0625 | 0.21738 |
| 2 | 0.0703 | 0.22746 |
| 3 | 0.0781 | 0.23503 |
| 4 | 0.0938 | 0.24407 |
| 5 | 0.1563 | 0.22924 |
| 6 | 0.2266 | 0.17580 |
| 7 | 0.2344 | 0.16987 |
| 8 | 0.5000 | 0.03674 |
| 9 | 0.8047 | -0.15223 |
| 10 | 0.8594 | -0.31117 |
| 11 | 0.9063 | -0.43423 |
| 12 | 0.9453 | -0.33511 |
| 13 | 0.9531 | -0.29032 |
| 14 | 0.9609 | -0.24095 |
| 15 | 0.9688 | -0.18864 |
| 16 | 1.0000 | 0.00000 |

Reference (raw data for transcription verification): the PICT arrays as-is —
Table 5 series: coordinates `[0.5000, 0.4453, 0.4375, 0.4297, 0.3984, 0.3281, 0.2187, 0.0469, 0.0000, -0.1172, -0.2344, -0.3516, -0.4531, -0.4609, -0.4688, -0.4766, -0.5000]`,
values `[0.00000, -0.20623, -0.22283, -0.23696, -0.27293, -0.25160, -0.10999, -0.00612, 0.00802, 0.03905, 0.07334, 0.12183, 0.33171, 0.39821, 0.48443, 0.58964, 1.00000]` (corresponding in the same order).
Table 6 series: coordinates `[-0.5000, -0.4375, -0.4297, -0.4219, -0.4062, -0.3437, -0.2734, -0.2656, 0.0000, 0.3047, 0.3594, 0.4063, 0.4453, 0.4531, 0.4609, 0.4688, 0.5000]`,
values `[0.00000, -0.21738, -0.22746, -0.23503, -0.24407, -0.22924, -0.17580, -0.16987, -0.03674, 0.15223, 0.31117, 0.43423, 0.33511, 0.29032, 0.24095, 0.18864, 0.00000]`.
(Table A keeps the sign as-is with coordinate z=0.5−x; Table B has the sign
flipped with coordinate x=y+0.5.
Verified point-by-point that the coordinates exactly match Ghia's 17 points.)

Note that PICT also contains, in the same format, A&K columns for "spanwise
aspect ratio 3" and "spanwise periodic boundary" (keys `(1000,1,3,False)`,
`(1000,1,1,True)`, etc.). **Do not reuse the spanwise-periodic system for
T15.5** — see the warning in §6.5.

### 3.3 Cross-Comparison Table of Re=1000 Extrema

Source: Ben Beya & Lili (2008) **Table 1** (visually transcribed from the
public PDF, full text). u_min is the extremum on the vertical centerline;
w_min/w_max are the extrema on the horizontal centerline. The values in ( )
are the relative difference from A&K as printed by the same paper — these
agree with this document's recomputation (§5), confirming mutual consistency
of the transcription.

| Source (grid) | u_min | z(u_min) | w_min | x(w_min) | w_max | x(w_max) |
|---|---|---|---|---|---|---|
| **Albensoeder & Kuhlmann 2005 (96×96×64, spectral)** | **-0.2803833** | **0.12419** | **-0.4350186** | **0.90957** | **0.2466511** | **0.10913** |
| Ben Beya & Lili 2008 (64³ FVM) | -0.2769995 (1.2%) | 0.1227 (1.2%) | -0.4295692 (1.25%) | 0.9041 (0.6%) | 0.2438928 (1.11%) | 0.1084 (0.6%) |
| Ben Beya & Lili (48³) | -0.2744732 | 0.1239 | -0.4294606 | 0.9100 | 0.2416644 | 0.1063 |
| Ben Beya & Lili (32³) | -0.2670297 | 0.1287 | -0.4133509 | 0.9214 | 0.2351729 | 0.1023 |
| Lo, Murugesan & Young 2005 (101³ FDM) | -0.26714 | 0.12 | -0.41534 | 0.92 | 0.23647 | 0.12 |
| Ding, Shu, Yeo & Xu 2006 (49³ MQ-DQ) | -0.258 | 0.12 | -0.414 | 0.92 | 0.225 | 0.12 |

LBM grid-refinement series (from the body text of Zhang, Shi & Wang 2015,
arXiv:1503.03337. D3Q14 MRT, z=0.5 symmetric half-domain, lid corner
regularization stated to be "the same setup as Ku"):

| Grid (half domain) | u_min (Re=1000) | Difference from A&K |
|---|---|---|
| 49×49×25 | -0.2619 | -6.6% |
| 65×65×33 | -0.2693 | -4.0% |
| 81×81×41 | -0.2730 | -2.6% |
| 97×97×49 | -0.2751 | -1.9% |

(Note that Table A's coordinates 0.1016 / 0.1719 straddle the true extremum
location z=0.12419, so it is correct behavior that the table's sample value
−0.27293 is shallower than the true extremum −0.28038. The test's extrema
comparison should compare the extremum obtained via a parabolic fit on the
grid (or similar) — not the "table point value" — against the A&K extrema
(the 7-digit values in this section).)

---

## 4. Smoothness Check (record of the mandatory ingestion-time inspection)

### 4.1 Method

Detect typo-type errors of the Ghia Re=400 v(0.9063)=−0.23827 kind (a single
point discontinuous with its neighbors). Since the coordinates are non-uniform,
use the **second-order divided difference**:

```
s_i  = (f_{i+1} − f_i) / (x_{i+1} − x_i)            (interval slope)
D2_i = (s_i − s_{i−1}) / (x_{i+1} − x_{i−1})        (second-order divided difference at interior point i)
```

D2 is equal to "the quadratic coefficient of the quadratic interpolant through
3 adjacent points." A single transcription error ε produces **alternating
large spikes of ± sign** (of magnitude ~ε/h²) at 3 consecutive points in D2,
so we confirm that the progression of sign and magnitude in the D2 sequence is
physically monotonic and smooth. For any suspicious point, we additionally
inspect the residual between the quadratic-interpolation predicted value Q(x_i)
from one-sided 3 points (i−1, i+1, i+2 or i−2, i−1, i+1) and the actual value
(a point is flagged as suspicious if the smaller of the two |residual| values
from either side substantially exceeds ~0.005·U).

> Implementation note: due to this session's constraints (code execution
> disallowed), this check was performed by **hand calculation (manual divided
> differences)**, and key values were double-checked. At test-data freeze time,
> codex should recompute the D2 columns below by script and confirm they agree
> with this document's values (to within rounding of the last digit).

### 4.2 Table A (u/U vertical centerline) results — all points pass

Interval slope s_i and second-order divided difference D2_i (rounded to 3–4 digits):

| i | z_i | u_i | s_i (i→i+1) | D2_i |
|---|--------|----------|---------|--------|
| 0 | 0.0000 | 0.00000 | -3.770 | — |
| 1 | 0.0547 | -0.20623 | -2.128 | +26.3 |
| 2 | 0.0625 | -0.22283 | -1.812 | +20.3 |
| 3 | 0.0703 | -0.23696 | -1.149 | +16.9 |
| 4 | 0.1016 | -0.27293 | +0.303 | +14.3 |
| 5 | 0.1719 | -0.25160 | +1.294 | +5.51 |
| 6 | 0.2813 | -0.10999 | +0.605 | -2.45 |
| 7 | 0.4531 | -0.00612 | +0.301 | -1.39 |
| 8 | 0.5000 | 0.00802 | +0.265 | -0.224 |
| 9 | 0.6172 | 0.03905 | +0.293 | +0.119 |
| 10 | 0.7344 | 0.07334 | +0.414 | +0.517 |
| 11 | 0.8516 | 0.12183 | +2.068 | +7.56 |
| 12 | 0.9531 | 0.33171 | +8.526 | +59.1 |
| 13 | 0.9609 | 0.39821 | +10.91 | +152 |
| 14 | 0.9688 | 0.48443 | +13.49 | +164 |
| 15 | 0.9766 | 0.58964 | +17.54 | +130 |
| 16 | 1.0000 | 1.00000 | — | — |

Judgment: D2 shows the physically natural progression "positive, monotonically
decaying near the bottom wall (+26→+14) → smooth sign transition through a
small value at the center → positive, sharply increasing in the lid boundary
layer (+7.6→+164)." **No alternating ± spikes. No points of concern.**

Additional inspection (quadratic-interpolation residuals, hand-calculation
examples for 2 representative points — for recomputation cross-check):

- i=4 (z=0.1016, near u_min): Newton quadratic interpolation via (z₂,z₃,z₅)
  Q(0.1016) = −0.22283 + (−1.81154)(0.0391) + 15.2418·(0.0391)(0.0313)
  = −0.27501. Actual value −0.27293 → residual +0.0021 (0.21% U) ✓
- i=15 (z=0.9766, lid boundary layer): via (z₁₃,z₁₄,z₁₆)
  Q(0.9766) = 0.39821 + 10.9139·(0.0157) + 143.498·(0.0157)(0.0078)
  = 0.58713. Actual value 0.58964 → residual +0.0025 (0.25% U) ✓
- i=12 (z=0.9531) appears to have a residual of −0.039 when predicted one-sidedly
  from the coarse-spacing side (i≤11), but is smooth at −0.0015 from the
  fine-spacing side (i13–i15). The degradation of the one-sided prediction at
  the point where the grid spacing switches from 0.1015→0.0078 is consistent
  with second-order truncation error (on the order of the product of curvature
  D2~59 and spacing 0.1). **Passes under the both-sides-minimum-residual criterion.**

### 4.3 Table B (w/U horizontal centerline) results — all points pass

| i | x_i | w_i | s_i (i→i+1) | D2_i |
|---|--------|----------|---------|--------|
| 0 | 0.0000 | 0.00000 | +3.478 | — |
| 1 | 0.0625 | 0.21738 | +1.292 | -31.1 |
| 2 | 0.0703 | 0.22746 | +0.971 | -20.6 |
| 3 | 0.0781 | 0.23503 | +0.576 | -16.8 |
| 4 | 0.0938 | 0.24407 | -0.237 | -10.4 |
| 5 | 0.1563 | 0.22924 | -0.760 | -3.94 |
| 6 | 0.2266 | 0.17580 | -0.760 | -0.001 |
| 7 | 0.2344 | 0.16987 | -0.501 | +0.947 |
| 8 | 0.5000 | 0.03674 | -0.620 | -0.209 |
| 9 | 0.8047 | -0.15223 | -2.906 | -6.36 |
| 10 | 0.8594 | -0.31117 | -2.624 | +2.77 |
| 11 | 0.9063 | -0.43423 | +2.542 | +60.1 |
| 12 | 0.9453 | -0.33511 | +5.742 | +68.4 |
| 13 | 0.9531 | -0.29032 | +6.329 | +37.6 |
| 14 | 0.9609 | -0.24095 | +6.622 | +18.6 |
| 15 | 0.9688 | -0.18864 | +6.046 | -14.7 |
| 16 | 1.0000 | 0.00000 | — | — |

Judgment: negative-curvature monotonic decay near the upstream wall (w_max
side), increasing positive curvature toward w_min (x≈0.91), an inflection
just before the wall — physically natural. **No alternating ± spikes. No
points of concern.**

**Case study — the point x=0.9063 (the same extraction coordinate as the Ghia typo incident)**:
w(0.9063) = −0.43423 also deviates from the quadratic prediction Q = −0.3584
(a deviation of 0.076) from the one-sided 3 points (x₉,x₁₀,x₁₂) that do not
straddle the extremum. Viewed in isolation, this looks like a similar
construction to the Ghia Re=400 typo (v(0.9063)=−0.23827, deviating ~0.14 from
its neighbor), but
- the prediction from the opposite side (x₁₀,x₁₂,x₁₃) is Q = −0.4418, consistent
  with a residual of +0.0076 (0.76% U)
- an independent source (the A&K extrema in Ben Beya Table 1) confirms that
  "the true minimum −0.4350186 is located at x=0.90957" (the table's sample
  point is nearly at the extremum)

For these reasons, this is **not a typo but a genuine sharp extremum**. "Do not
reject based solely on a one-sided prediction deviation. Judge using
both-sided prediction plus independent-source cross-check" — the same
two-stage judgment used for the Ghia-incident check should also be implemented
in the test ingestion script as an operational rule.

### 4.4 Endpoint Inspection

Table A: u(0)=0 (stationary wall), u(1)=1 (lid speed). Table B: w(0)=w(1)=0
(stationary wall). Both are exact boundary-condition values, correct
independently of the interior points. ✓

---

## 5. Cross-Source Comparison (order of the differences)

Relative comparison of the §3.3 extrema against the A&K (spectral) baseline:

| Source | \|u_min\| diff | \|w_min\| diff | w_max diff |
|---|---|---|---|
| Ben Beya 64³ | -1.2% | -1.25% | -1.1% |
| Ben Beya 48³ | -2.1% | -1.28% | -2.0% |
| Ben Beya 32³ | -4.8% | -5.0% | -4.7% |
| Lo 101³ | -4.7% | -4.5% | -4.1% |
| Ding 49³ | -8.0% | -4.8% | -8.8% |
| Zhang LBM 97×97×49 | -1.9% | — | — |
| Zhang LBM 65×65×33 | -4.0% | — | — |
| Zhang LBM 49×49×25 | -6.6% | — | — |

Observations:
1. Non-spectral solutions (N≈32–101) produce extrema that are **1–8% shallower**
   (systematically on the numerical-diffusion side). Refining the grid brings
   them monotonically closer to A&K (both in the Ben Beya series and the Zhang series).
2. **Grid resolution dominates** over differences between clever methods. ~1–4% at
   N≈64 (equivalent), ~2–5% at N≈100 (method-dependent), LBM is ~2% at N≈97 (half domain).
3. Agreement in location (z, x) is better than in value, within ±0.01–0.02L.
4. Where full-profile cross-checking was possible (A&K symbols in Ben Beya Fig.2
   vs. this document's Table A/B), agreement is within figure-reading precision (±0.02 U).

→ This scatter is the empirical basis for the acceptance bands in §6.

---

## 6. Recommendations to codex (guidance for building the tests)

### 6.1 Recommended Test Configuration

- **T15.5a (first frozen round, Re=1000 only)**: cube N³ (full-domain
  computation; do not use the symmetric-half-domain shortcut — symmetry itself
  is a check target of the test), top face MovingWall {[U,0,0]}, the other 5
  faces BounceBack, TRT Λ=3/16, U=0.1, Re=1000 → ν = U(N−2)/1000.
- Reference: §3.2 Table A (u/U, 17 points) and Table B (w/U, 17 points), §3.3 A&K extrema.
- Comparison follows the same convention as 2D T7: project from cell-center
  values onto the centerline (via 2-point linear interpolation to the
  reference coordinates if needed), and take the RMS over the 17 points
  (endpoints included). Estimate extrema via a local parabolic fit rather than
  the raw grid value, and compare against that.

### 6.2 Grid Resolution Guidance and Stability Constraints

- **Grid-Re constraint (empirical guidance from T10)**: Re/(N−2) = U/ν ≲ 15.
  At Re=1000, **N=64 gives 16.1, slightly over** (likely to still run with TRT
  3/16, but with no margin). **Recommend N ≥ 72 (Re/L=14.3)**; setting the
  standard at N=80 (12.8) is safe. N=96 (10.6) is for the #[ignore] full
  validation.
- Convergence: steady-state criterion ε=1e-8 (Δ=500 steps) or an upper bound of
  ~500k steps. The convective time L/U = (N−2)/U ≈ 780 steps (N=80, U=0.1), so
  convergence is slower than 2D T7 (the 3D vortex structure takes ~50–100 L/U to develop).

### 6.3 Recommended Acceptance Bands (initial bands, as a precursor to empirically frozen ones)

From the 2D Ghia empirical rule (T7: RMS ≤ 0.02·U (Re100/400) / 0.03·U (Re1000),
N=129, measured ~0.005·U) and the §5 cross-source scatter in 3D:

| Grid | u-line RMS | w-line RMS | Extrema |
|---|---|---|---|
| N=64–80 (standard suite) | ≤ 0.030·U | ≤ 0.035·U | \|u_min\|, \|w_min\| within ±6% of A&K, w_max ±8%, location ±0.03L |
| N=96–128 (#[ignore]) | ≤ 0.020·U | ≤ 0.025·U | ±3%, location ±0.02L |

Breakdown of the rationale (approximate): uncertainty in the reference itself
~1e-5 (negligible), LBM-side grid bias (from the §5 Zhang series: extrema
−3 to −4% at N≈64-equivalent, = profile RMS contribution ~0.005–0.010·U),
interpolation error to the reference points (D2~150 at the 4 lid-boundary-layer
points → ~0.003·U at h=1/78), O(Ma²) error (~0.003·U at U=0.1). The combined
expected RMS is on the order of 0.008–0.015·U, and the table above carries a
2–3x margin. This is the same margin ratio as the 2D track record (measured
0.005 against a band of 0.02).
**The reason the w-line band is set wider than the u-line band** is that
w_min (−0.435) is deeper than u_min (−0.280) so the grid bias has a larger
absolute magnitude, and the extraction point x=0.9063 sits nearly exactly at the extremum.

### 6.4 Additional Asserts Worth Adding (cheap, high-sensitivity)

1. **Symmetry**: the spanwise velocity on the symmetry plane of the steady
   solution, max|v(x, y=mid, z)| ≤ 1e-3·U (physically zero. A violation
   detects either non-convergence or a symmetry bug. If N is even, mid falls
   between cells, so evaluate as the average of the 2 adjacent planes).
2. **2D-confusion detection**: the RMS difference between Table A and the 2D
   Ghia Re=1000 u-table is approximately ~0.1·U (maximum point difference
   ~0.21, u_min: 3D −0.280 vs 2D −0.383). Adding an inverse assert ("fail if
   the 3D solution incorrectly matches the 2D table," i.e., RMS vs 2D Ghia
   ≥ 0.05·U) can detect a regression where the z-direction implementation is
   dead and the solution has effectively collapsed to 2D.
3. **Convergence trend**: |u_min(N=96) − A&K| < |u_min(N=64) − A&K|
   (same convention as T8).
4. **Mass conservation**: since it's an all-wall BB box, total mass drift
   ≤ 1e-11 relative per 10⁴ steps (same as T6).

### 6.5 Scope Warnings (physics)

- **Upper bound on steadiness**: for the cube (all walls), the
  steady-to-oscillatory transition is at Re_cr ≈ 1914 (Feldman & Gelfgat 2010;
  confirmed by the PIV experiment of Liberzon et al. arXiv:1107.0449, ω≈0.575 —
  confirmed by full-text reading. Note that the linear-stability value
  ≈1919.5 from Kuhlmann & Albensoeder 2014 is a secondhand citation via the
  nekStab repository's example README, and the original has not been
  verified). **Re=1000 is safely within the steady regime**, but adding an
  Re=2000 case out of 2D habit is not acceptable, as it would land just above
  the critical value. If Re is to be increased, go up to about 1500 at most.
- **Do not confuse with spanwise-periodic boundaries**: the spanwise-periodic
  cavity is already in the Taylor–Görtler-type instability regime at Re=1000,
  and A&K's periodic-boundary column (bundled in PICT) is neither a 2D
  solution nor a "z-invariant 3D solution" but a **cellular 3D steady
  solution**. Whether a D3Q19 run with z-periodicity, started from a
  symmetric initial condition, stays at the 2D solution or falls into TGL
  cells depends on rounding perturbations, making it indeterminate as a test.
  T15.5 uses **only the all-wall (rigid end wall) configuration**.
- **Difference in lid-corner regularization**: even among the references, the
  treatment of the lid-edge singularity differs (A&K treats the singularity
  analytically; Ku/Zhang apply a u=0.3 ramp at the lid-adjacent nodes).
  LBMFlow's half-way BB moving wall plus "faster wall wins" corner rule is yet
  another discrete regularization, and this difference concentrates in the 4
  points at z ≥ 0.95. During calibration, log the RMS both over "all 17
  points" and over "the 13 points excluding the 4 lid-boundary-layer points"
  (the spec's band applies to all 17 points).

### 6.6 Convention for Embedding the Data

- Embed Table A/B and the A&K extrema as constants in the test file, with a
  comment stating "Albensoeder & Kuhlmann, JCP 206 (2005) 536–558, Tables 5/6;
  transcribed via github.com/tum-pbs/PICT tests/validations.py @a95d7f9;
  smoothness-checked: docs/T15_5_CAVITY3D_REFERENCE.md §4" (the same convention
  as T7's Ghia table. A short numeric table is factual data and poses no
  copyright concern).
- Have the ingestion script recompute the §4 D2 columns and cross-check them
  against this document's values (to within the last digit) before freezing
  (this doubles as an independent verification of the hand calculation).

---

## 7. What Could Not Be Obtained, and Future Actions

**Numeric values not obtained (paywall, etc.)** — reference name only listed, no numeric values:

1. **Ku, Hirsh & Taylor (1987)**, "A pseudospectral method for solution of the
   three-dimensional incompressible Navier-Stokes equations",
   *J. Comput. Phys.* **70**:439–462.
   <https://www.sciencedirect.com/science/article/pii/0021999187901902>
   The original benchmark for Re=100/400/1000. **Numeric values not obtained**.
   Only surrounding information confirmed: the grid is Re=400: 25×25×13 /
   Re=1000: 31×31×16 (non-uniform, z=0.5 symmetric half domain) [per the
   description in Kato et al. 1990], velocity ramp regularization at the
   lid-adjacent nodes [per the description in Zhang et al. 2015]. Kato et al.
   (1990) Fig.3/7 reproduces the Re=400/1000 profiles as circular symbols, but
   this is figures only, no numeric table.
2. **Wong & Baker (2002)**, "A 3D incompressible Navier–Stokes velocity–vorticity
   weak form finite element algorithm", *Int. J. Numer. Meth. Fluids*
   **38**:99–123 (Re=100/400, 48³). **Numeric values not obtained** (Wiley
   paywall; no public mirror or reprint via a thesis, etc. could be found either).
3. **Yang, Yang, Chen & Hsu (1998)**, "Implicit weighted ENO schemes for the
   three-dimensional incompressible Navier–Stokes equations",
   *J. Comput. Phys.* **146**:464–487. **Numeric values not obtained**.
   (A public PDF of a 1999 IJNMF **31**:747–765 paper by authors of the same
   names was obtained, but that is a **2D paper** and has no 3D table — take
   care not to confuse the two.)
4. **Lo, Murugesan & Young (2005)** IJNMF **47**:1469–1487 (101³) and
   **Ding, Shu, Yeo & Xu (2006)** CMAME **195**:516–533 (49³):
   profile tables not obtained. **Only the Re=1000 extrema** were obtained via
   Ben Beya & Lili (2008) Table 1 (§3.3).
5. **Iwatsu, Ishii, Kawamura, Kuwahara & Hyun (1989)**, "Numerical simulation
   of three-dimensional flow structure in a driven cavity",
   *Fluid Dyn. Res.* **5**:173–189 (81³). **Numeric values not obtained**.
6. **Albensoeder & Kuhlmann (2005) original PDF**: not obtained (this
   document's numeric values are via the redistribution and independent
   cross-check described in §3.1).

**Accessed, but the numeric table did not exist** (this is "non-existent," not "not obtained"):

7. **Jiang, Lin & Povinelli**, NASA TM-105904 / ICOMP-93-06 (1993)
   (= the TM version of *Comput. Methods Appl. Mech. Engrg.* **114** (1994) 213–231),
   <https://ntrs.nasa.gov/citations/19930013636> — read in full, all 22 pages.
   The results for the 3D cavity (Re=100/400/1000/2000/3200, three grids at
   the 50³ level) are **only vector plots and vorticity contour maps in
   Fig.4–6, with no numeric table of velocity**. This TM also explicitly
   states that "a steady solution could not be obtained (with their LSFEM) for
   Re ≥ 1000," so it is unsuitable for adoption as an Re=1000 reference.
8. **Kato, Kawai & Tanahashi (1990)**, JSME Int. J. Ser. II **33**(4):649–658,
   published on J-STAGE <https://www.jstage.jst.go.jp/article/jsmeb1988/33/4/33_4_649/_article>
   — read in full. Re=400/1000/3200 (40×40×20 half domain). **Figures only
   (Fig.3, 7, 13–16); the only numeric table is the parameter table (Table 1)**.
9. **Žunič, Hriberšek, Škerget & Ravnik (2006)**, "3D lid driven cavity flow by
   mixed boundary and finite element method", ECCOMAS CFD 2006. Public mirror
   <https://www.ljll.fr/~frey/papers/Navier-Stokes/Zunic%20Z.,%203D%20Lid%20driven%20cavity%20flow%20by%20mixed%20boundary%20and%20finite%20element%20method.pdf>
   — read in full. Re=100/1000, **figures only** (their reference is the
   point data of Yang et al. 1998).
10. **Kuhlmann & Romanò (2019)**, "The lid-driven cavity" (review chapter,
    *Computational Modelling of Bifurcations and Instabilities in Fluid
    Dynamics*, Springer), public PDF:
    <https://www.fluid.tuwien.ac.at/HendrikKuhlmann?action=AttachFile&do=get&target=LidDrivenCavity.pdf>
    — no reproduction of a cube centerline-profile numeric table found within
    the scope checked (useful as background literature).

**Actions (in priority order)**:

- **A-1**: obtain the A&K (2005) original via a university library / ILL, and
  cross-check the 34 values in §3.2 against the printed values (7 digits) in
  Tables 5/6. If there is a discrepancy, treat the original as authoritative
  and revise this document.
  (Also resolve, using the original, the discrepancy between PICT's comment
  grid notation "96×96×96" and Ben Beya's citation "96×96×64." This is a
  bibliographic note that does not affect the values.)
- **A-2**: at the same time, obtain Ku et al. (1987) and Wong & Baker (2002),
  and confirm whether numeric tables for Re=100/400 exist. If tables exist,
  append them to §3 of this document, run them through the same smoothness
  check as §4, and then file T15.5b (Re=100/400).
- **A-3**: until then, freeze the T15.5 tests with Re=1000 only (§6).

---

## 8. Collection Log (all accessed 2026-07-05)

| Order | Source | Method | Result |
|---|---|---|---|
| 1 | NTRS 19930013636 (Jiang TM) PDF | WebFetch + close reading of page images | Confirmed figures only |
| 2 | Žunič et al. ECCOMAS 2006 (ljll.fr) PDF | Same as above | Confirmed figures only |
| 3 | arXiv:1503.03337 (Zhang et al.) PDF | Same as above | Collected the u_min grid-refinement series |
| 4 | arXiv:1704.08521 (Gelfgat) / arXiv:1107.0449 (Liberzon et al.) PDF | Same as above | No profile table / confirmed Re_cr≈1914, ω≈0.575 |
| 5 | C.R. Mécanique 336:863 (Ben Beya & Lili) PDF | Same as above (visually transcribed Table 1, double-verified against the paper's % annotations and by recomputation) | §3.3 extrema table |
| 6 | github.com/tum-pbs/PICT `tests/validations.py` @a95d7f9 | Fetched raw, close reading with line numbers (transcription quality verified against the 2D Ghia table) | §3.2 Table A/B |
| 7 | J-STAGE jsmeb1988 33-4-649 (Kato et al.) PDF | Fetched via direct URL guess, close reading | Confirmed figures only, Ku grid information |
| 8 | fluid.tuwien.ac.at (Kuhlmann & Romanò chapter) PDF | Fetched via curl, read the opening | No numeric table (review) |
| 9 | Elsevier (Ku 1987, A&K 2005), Wiley (Wong & Baker, Lo et al.), IOP (Iwatsu et al.) | — | Abandoned due to paywall (§7) |

(The Semantic Scholar API returned 429 and was not used. web.archive.org could
not be reached from this environment. The CiteSeerX / astro.yale.edu mirrors
were not reachable via WebFetch due to a TLS certificate issue; "weno_99.pdf"
from astro.yale.edu, obtained via curl -k, turned out to be the 2D paper (§7-3).)
