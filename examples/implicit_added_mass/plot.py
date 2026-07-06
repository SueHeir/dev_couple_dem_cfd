#!/usr/bin/env python3
"""Run the implicit_added_mass validation and plot explicit vs Aitken traces."""

import csv
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CFG = "examples/implicit_added_mass/config.toml"
CSV = HERE / "plots" / "trace.csv"


def run():
    out = subprocess.run(
        ["cargo", "run", "--release", "--example", "implicit_added_mass", "--", CFG],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    text = out.stdout + out.stderr
    if out.returncode != 0:
        print(text)
        sys.exit(out.returncode)
    print(text)
    return text


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    text = run()
    passed = "VALIDATION: PASS" in text
    growth = grab(text, r"explicit residual growth: .* = ([0-9.eE+-]+)")
    final_v = grab(text, r"implicit v_final: ([0-9.eE+-]+)")
    x_star = grab(text, r"analytic fixed point ([0-9.eE+-]+)")
    residual = grab(text, r"implicit residual: ([0-9.eE+-]+)")
    iters = int(grab(text, r"implicit residual: [0-9.eE+-]+; iters ([0-9]+)"))

    explicit = []
    with CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["kind"] == "explicit":
                explicit.append(
                    (
                        int(row["step"]),
                        abs(float(row["residual"])),
                        float(row["v_after"]),
                    )
                )

    fig, (ax_res, ax_v) = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=False)
    ax_res.semilogy(
        [r[0] for r in explicit],
        [r[1] for r in explicit],
        "o-",
        color="#b3261e",
        label=f"explicit couple_two_way residual ({growth:.1e}x growth)",
    )
    ax_res.axhline(residual, color="#1f7a3a", lw=1.6, label=f"Aitken final residual {residual:.1e}")
    ax_res.set_ylabel("|v_tilde - v|")
    ax_res.set_title(
        "Same DEM-CFD seam: explicit diverges, Aitken converges - "
        + ("PASS" if passed else "FAIL"),
        fontsize=11,
    )
    ax_res.grid(True, which="both", alpha=0.28)
    ax_res.legend(fontsize=8)

    ax_v.plot(
        [r[0] for r in explicit],
        [r[2] for r in explicit],
        "o-",
        color="#b3261e",
        label="explicit particle velocity",
    )
    ax_v.axhline(x_star, color="#444444", ls="--", lw=1.5, label=f"analytic fixed point {x_star:.4f}")
    ax_v.plot([iters], [final_v], "s", ms=8, color="#1f7a3a", label=f"Aitken converged in {iters} iters")
    ax_v.set_xlabel("step / outer iteration")
    ax_v.set_ylabel("interface velocity")
    ax_v.grid(True, alpha=0.28)
    ax_v.legend(fontsize=8)

    fig.tight_layout()
    outdir = HERE / "plots"
    outdir.mkdir(exist_ok=True)
    dst = outdir / "implicit_added_mass.png"
    fig.savefig(dst, dpi=150)
    print(f"wrote {dst} ({'PASS' if passed else 'FAIL'})")


def grab(text, pat):
    m = re.search(pat, text)
    if not m:
        sys.exit(f"could not parse pattern: {pat}")
    return float(m.group(1))


if __name__ == "__main__":
    main()
