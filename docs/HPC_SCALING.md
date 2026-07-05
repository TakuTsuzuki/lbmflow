# HPC Scaling Design Note (2026-07-05)

Current answer to "does it scale to a supercomputer" (M-D update): **distributed-memory support is done
(lbm-core feature `mpi`; through single-node multi-rank validation and weak-scaling reference values).
Multi-node measurement awaits a cluster.** Staged plan 1–3 are implemented (of 3, communication overlap
and parallel I/O are deferred to M-E) — the procedure, measurements, and cluster-measurement list are in docs/MPI_GUIDE.md.
The following is the original design guidance (kept as a record).

## Current assets (design items already in place that help distribution)

| Asset | Meaning for distribution |
|---|---|
| pull scheme + double buffering | isomorphic to the halo-exchange pattern as-is |
| collision / Shan-Chen force are fully cell-local | closes within a subdomain (only adds 1 halo layer of ψ) |
| boundary conditions = edge spec + solid mask (data) | distribution to subdomains is trivial |
| deviation storage (f−w) | f32-izing halo communication to halve bandwidth is realistic |
| scenario JSON layer | the entry point is unchanged by swapping in a distributed runner |

## Blockers

1. No domain-decomposition abstraction (global nx×ny, periodicity is global modulo)
2. No halo/ghost layer (stream reads the global neighborhood directly)
3. No inter-process communication (rayon shared memory only)
4. Diagnostics are whole-domain serial reductions (total_mass / steady-state criterion / probes)
5. Output assumes whole-field bulk (PNG/CSV/manifest)
6. wgpu is oriented to single-node/Web. Multi-GPU / RDMA are in the CUDA/HIP/SYCL family

## Staged plan (co-designed with Phase 10)

1. **Subdomain abstraction**: local grid + width-1 halo + neighbor links + `HaloExchange` trait.
   Current = single-subdomain implementation. Do it at the same time as the 3D-ization index abstraction (avoids doing the work twice).
2. **In-process multi-subdomain**: a match test that partitioned execution ≡ monolithic execution
   (near-bit match can be required in f64; a candidate for the new T13 category of the codex adversarial suite).
3. **MPI backend** (rsmpi): implement the exchange trait + Allreduce diagnostics +
   parallel I/O (parallel VTK or HDF5). Overlap of communication and internal computation
   (the 2-pass structure that defers the boundary layer is compatible with the current stream's row partitioning).
4. **Per-rank accelerator**: based on the results of the wgpu evaluation (phase9-wgpu branch),
   if HPC is the real target, consider cudarc/HIP as an interchangeable backend.

## ROI note

2D LBM is enough on a single node (even now, 1024² at 380 MLUPS). The value of distribution
first appears in 3D (D3Q19, 10⁹-lattice class). **Include the Subdomain abstraction as a mandatory
requirement in the Phase 10 design phase** (building 3D standalone first and then distributing it means a rewrite).
