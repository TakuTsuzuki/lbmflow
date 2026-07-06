# Dispersed Seeding 3D Demo

This example is a small runnable 3D dispersed-particle withdraw, eject, and
deposition demo. It parses the protocol JSON, builds D3Q19 LBM domains through
the public `lbm_core` API, extracts a particle batch from a tall reservoir, then
injects that batch into a shallow target tray and writes deposition readouts.

Run:

```bash
cargo run --release -p lbm-cli --example dispersed_seeding -- \
  crates/lbm-cli/examples/dispersed_seeding/sample_gentle.json
cargo run --release -p lbm-cli --example dispersed_seeding -- \
  crates/lbm-cli/examples/dispersed_seeding/sample_harsh.json
```

Outputs are written under the JSON `output.dir`:

- `density.csv`: `M x N` deposition count and normalized density.
- `metrics.json`: CV, max/mean, empty-bin fraction, deposited, suspended and
  extracted counts, and regime numbers.
- `reservoir_velocity.vtk` and `tray_velocity.vtk`: ASCII VTK 3D velocity
  volumes.

## Model

The world is 3D with horizontal `x,y` and vertical `z`. Gravity is `-z`. The
domains are fluid-filled; there is no free surface. This is an MVP
simplification, not a claim about a gas-liquid interface.

Particles are one-way Lagrangian markers. Diameter is sampled from the requested
lognormal distribution. The tray integrator uses Stokes drag with a
semi-implicit exponential update, buoyancy-corrected gravity, and the
translational agitation acceleration from the protocol. A particle that crosses
the floor is frozen and counted into the partition bin containing its final
`x,y` position.

The reservoir phase uses an explicit 1D settling-column extraction model. The
current public core API exposes face boundary conditions, not an internal
localized suction outlet with a simultaneously open replenishment top. The
example therefore runs a closed 3D reservoir LBM visualization phase and
extracts particles from the settled column concentration at the requested
withdraw depth. The extracted diameter histogram is printed so the substitution
remains measurable.

The tray phase uses a D3Q19 domain with wall sides and floor, plus a localized
top velocity profile through `set_inlet_profile_with`. The top face is a
velocity face with zero velocity outside nozzle disks. The particle integrator
combines that sampled LBM field with an example-local impinging wall-jet and
near-wall dispersion closure because the public face boundary cannot represent
both localized injection and open top outflow on the same face.

## Non-goals

This example does not implement free-surface flow, gas-liquid interfaces,
two-way or four-way coupling, particle-particle collisions, an inverse-design
solver, GPU execution, or any core API changes.
