#!/usr/bin/env python3
"""Run fluidized_bed_umf and plot the measured U_mf validation gate.

    $BENCH_PYTHON examples/fluidized_bed_umf/plot.py

The plot is regenerated from the example's own output: the live seam bisection,
the dynamic a_z zero crossing, the Wen-Yu reference, and the two negative controls.
"""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CFG = "examples/fluidized_bed_umf/config.toml"


def run():
    out = subprocess.run(
        ["cargo", "run", "--release", "--example", "fluidized_bed_umf", "--", CFG],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    text = out.stdout + out.stderr
    if out.returncode != 0:
        sys.stderr.write(text)
        sys.exit(out.returncode)
    return text


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    text = run()
    u_meas = _grab(text, r"U_mf MEASURED .*: ([\d.]+) m/s")
    u_ref = _grab(text, r"U_mf Wen&Yu .*: +([\d.]+) m/s")
    rel_err = _grab(text, r"rel\.err ([\d.]+)%")
    tol = _grab(text, r"rel\.err [\d.]+% +\(tol ([\d.]+)%\)") / 100.0
    u_dyn = _grab(text, r"U_mf DYNAMIC .*: +([\d.]+) m/s", required=False)
    u_nopg = _grab(text, r"negative controls: omit-\S+ U_mf ([\d.]+)")
    u_epsbug = _grab(text, r"eps-power-bug U_mf ([\d.]+)")
    passed = "VALIDATION: PASS" in text

    labels = ["full seam\nbisection", "omit grad-P\ncontrol", "eps-power\ncontrol"]
    values = [u_meas, u_nopg, u_epsbug]
    colors = ["#1f77b4", "#d62728", "#d62728"]

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.axhspan(
        u_ref * (1.0 - tol),
        u_ref * (1.0 + tol),
        color="#2ca02c",
        alpha=0.16,
        label=f"pass band: Wen-Yu +/- {tol * 100:.0f}%",
    )
    ax.axhline(u_ref, color="#2ca02c", lw=1.6, ls="--", label=f"Wen-Yu reference {u_ref:.4f} m/s")
    bars = ax.bar(labels, values, color=colors, alpha=0.82)
    bars[0].set_label(f"live seam bisection {u_meas:.4f} m/s")
    if u_dyn is not None:
        ax.scatter([0], [u_dyn], marker="D", s=55, color="#ff7f0e", zorder=4, label=f"dynamic onset {u_dyn:.4f} m/s")

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            value + 0.025,
            f"{value:.4f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.set_ylabel("minimum fluidization velocity U_mf [m/s]")
    ax.set_title(
        "fluidized_bed_umf - live U_mf bisection vs Wen-Yu\n"
        + ("PASS" if passed else "FAIL")
        + f"  ({rel_err:.2f}% error, tolerance {tol * 100:.0f}%)"
    )
    ax.set_ylim(0.0, max(values + [u_ref * (1.0 + tol)]) * 1.18)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()

    outdir = HERE / "plots"
    outdir.mkdir(exist_ok=True)
    dst = outdir / "umf_validation.png"
    fig.savefig(dst, dpi=130)
    print(f"wrote {dst}  ({'PASS' if passed else 'FAIL'})")


def _grab(text, pat, required=True):
    m = re.search(pat, text)
    if not m:
        if required:
            sys.exit(f"could not parse {pat!r} from run output")
        return None
    return float(m.group(1))


if __name__ == "__main__":
    main()
