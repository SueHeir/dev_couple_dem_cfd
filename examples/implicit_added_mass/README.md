# implicit_added_mass

This validation is a minimal dense-regime / added-mass stress case for the existing
DEM-CFD seam. A one-particle SOIL bed and one FIELD sub-App exchange kinematics and
fluid load through `crates/dem_cfd/src/seam.rs`. The FIELD force model is the
linearized interface map `v_tilde = -3 v + 1`, so plain explicit partitioning has
spectral radius `3` and must oscillate away from the fixed point `v* = 0.25`.

The explicit run is the stock `couple_two_way` path. The implicit run uses the same
export, CFD tick, import, and SOIL tick systems, with
`grass_multi::converge_outer_iter` and `Relaxation::Aitken` driving the outer
interface velocity. This demonstrates shared-state access: the parent can read and
relax the interface state in-place across sub-Apps, then inject the relaxed value
back through the same seam without adding a new coupling primitive or copying
through a separate broker.

![explicit divergence and Aitken convergence](plots/implicit_added_mass.png)

Caption: explicit `couple_two_way` diverges on the strong added-mass map, while the
same seam phases converge to the analytic fixed point with Aitken relaxation. PASS
requires explicit residual growth of at least `1e3`, Aitken convergence to `v*`
within `1e-9`, and final residual below `1e-10`.

Reproduce:

```bash
cargo run --release --example implicit_added_mass -- examples/implicit_added_mass/config.toml
$BENCH_PYTHON examples/implicit_added_mass/plot.py
```
