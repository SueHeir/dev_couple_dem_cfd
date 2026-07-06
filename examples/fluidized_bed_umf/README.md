# fluidized_bed_umf

Minimum fluidization is measured by bisection on the live DEM-CFD seam force for
an FCC SOIL bed in an imposed upward gas stream, then cross-checked against the
Wen & Yu (1966) minimum-fluidization correlation. The measured seam uses the
independent MacDonald packed-bed closure, so the nonzero difference from Wen & Yu
is an inter-correlation validation spread rather than a self-comparison.

![U_mf versus Wen & Yu](plots/umf_wen_yu.png)

Figure: `fluidized_bed_umf` PASS. The live seam gives `U_mf = 0.5138 m/s` versus
Wen & Yu `0.5380 m/s`, a `4.51%` relative error inside the fixed `15%` gate; the
dynamic zero crossing matches the seam value, while both independently run
negative controls miss that same band. The shaded band and control markers show
the actual acceptance comparison, not a fitted reference curve.

This is a closure-level coupling check, not an experimental fluidized-bed
validation: gas flow is imposed, the bed remains at the prescribed packing
voidage during the short onset measurement, and the result does not establish
predictive accuracy for evolving industrial beds. Wen & Yu (1966) is the
external reference and MacDonald et al. (1979) supplies the distinct measured
closure; full citations are in [data/references.md](data/references.md). The
repository-level AI-authorship and review disclaimer also applies to this
example.

Regenerate the figure from the existing demo:

```bash
$BENCH_PYTHON examples/fluidized_bed_umf/plot.py
```
