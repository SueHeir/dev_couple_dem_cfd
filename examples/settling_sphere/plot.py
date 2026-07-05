#!/usr/bin/env python3
"""Run the settling_sphere DEM-CFD gate and plot its approach to terminal velocity.

    $BENCH_PYTHON examples/settling_sphere/plot.py

Runs the example, parses the `step  v_z` history and the header reference values,
and writes plots/terminal_velocity.png — the settling particle's vertical velocity
relaxing onto the Stokes (1851) / Schiller-Naumann terminal value through the live
DEM<->CFD drag seam. Committed so the figure renders in the README on Gitea.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CFG = "examples/settling_sphere/config.toml"


def run():
    out = subprocess.run(
        ["cargo", "run", "--release", "--example", "settling_sphere", "--", CFG],
        cwd=ROOT, capture_output=True, text=True,
    )
    return out.stdout + out.stderr


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    text = run()
    steps, vz = [], []
    for line in text.splitlines():
        m = re.match(r"\s*(\d+)\s+(-?\d+\.\d+)\s+", line)
        if m and "step" not in line:
            steps.append(int(m.group(1)))
            vz.append(-float(m.group(2)))  # plot settling speed (magnitude)
    v_stokes = _grab(text, r"v_stokes[^\d-]*([\d.]+)")
    v_sn = _grab(text, r"v_balance\(SN\)[^\d-]*([\d.]+)")
    passed = "VALIDATION: PASS" in text

    if not steps:
        sys.exit("could not parse a step/v_z history from the run output")

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(steps, vz, "o-", color="#1f77b4", lw=1.8, ms=4, label="settling speed |v_z| (DEM↔CFD seam)")
    if v_stokes:
        ax.axhline(v_stokes, color="#d62728", ls="--", lw=1.4, label=f"Stokes 1851  {v_stokes:.4f} m/s")
    if v_sn:
        ax.axhline(v_sn, color="#2ca02c", ls=":", lw=1.4, label=f"Schiller–Naumann  {v_sn:.4f} m/s")
    ax.set_xlabel("time step")
    ax.set_ylabel("settling speed  |v_z|  [m/s]")
    ax.set_title("settling_sphere — terminal velocity through the DEM↔CFD drag seam\n"
                 + ("PASS" if passed else "FAIL") + f"  (v_t={vz[-1]:.4f} vs Stokes {v_stokes:.4f} m/s)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    outdir = HERE / "plots"
    outdir.mkdir(exist_ok=True)
    dst = outdir / "terminal_velocity.png"
    fig.savefig(dst, dpi=130)
    print(f"wrote {dst}  ({'PASS' if passed else 'FAIL'})")


def _grab(text, pat):
    m = re.search(pat, text)
    return float(m.group(1)) if m else None


if __name__ == "__main__":
    main()
