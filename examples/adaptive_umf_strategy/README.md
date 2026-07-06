# Adaptive UMF Coupling Strategy

This example compares coupling schedules for a reduced unresolved DEM-CFD bed near
minimum fluidization. The force model is the same science target used by
`fluidized_bed_umf`: MacDonald packed-bed drag plus pressure-gradient buoyancy is
assembled through the live DEM-CFD seam, and the measured minimum-fluidization
velocity is checked against the Wen & Yu correlation.

The strategy comparison keeps the GRASS execution graph predeclared. The parent
schedule is the existing `dem_cfd::couple_two_way` sequence:

```text
export DEM kinematics -> tick CFD seam -> import fluid force -> tick SOIL bed
```

The strategies differ only in how many times that schedule is ticked inside one
macro coupling interval:

| strategy | behavior |
|---|---|
| `fixed_explicit` | one export/seam/import/SOIL tick per macro interval |
| `fixed_subcycle_4` | four fixed substeps per macro interval |
| `residual_gated` | tries 1, 2, 4, then 8 substeps and accepts the first with residual below the configured gate |

The residual is
`\|a_z(measured after the macro interval) - a_z(live seam force / M_bed)\| / g`.
That makes the plot show both the measured-vs-reference fluidization onset and the
coupling-strategy behavior.

![adaptive UMF strategy validation](plots/adaptive_umf_strategy.png)

The plotted run passes: live-seam `U_mf = 0.5138 m/s` versus Wen & Yu
`0.5380 m/s` (`4.51%`, tolerance `15%`), all dynamic zero crossings match the seam
within `5%`, and the residual-gated strategy reduces the worst residual from
`0.0699 g` to `0.0133 g` (`5.25x`) while satisfying the `0.015 g` gate. Negative
controls still fail the Wen & Yu gate: omitting pressure-gradient buoyancy shifts
`U_mf` by `+80.4%`, and the corrupted epsilon-power reduction shifts it by
`-53.2%`.

Regenerate with:

```bash
$BENCH_PYTHON examples/adaptive_umf_strategy/plot.py
```
