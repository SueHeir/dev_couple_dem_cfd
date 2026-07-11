#!/usr/bin/env python3
"""Run fluidized_bed_umf and plot the U_mf validation gate.

    $BENCH_PYTHON examples/fluidized_bed_umf/plot.py

The figure is regenerated from the existing demo stdout: measured bed acceleration
versus superficial velocity, the live-seam U_mf, the Wen & Yu reference, the 15%
pass band, and the negative-control U_mf shifts.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CFG = "examples/fluidized_bed_umf/config.toml"


def run() -> str:
    proc = subprocess.run(
        ["cargo", "run", "--release", "--example", "fluidized_bed_umf", "--", CFG],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    text = proc.stdout + proc.stderr
    if proc.returncode != 0:
        print(text)
        sys.exit(proc.returncode)
    return text


def grab(text: str, pattern: str) -> float:
    match = re.search(pattern, text)
    if not match:
        sys.exit(f"could not parse pattern: {pattern}")
    return float(match.group(1))


def parse(text: str) -> dict[str, object]:
    row_re = re.compile(
        r"^\s*([0-9.]+)\s+([0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([0-9.]+)\s+",
        re.MULTILINE,
    )
    rows = [
        {
            "u": float(u),
            "re": float(re_p),
            "a_meas": float(a_meas),
            "a_force": float(a_force),
            "handoff": float(handoff),
        }
        for u, re_p, a_meas, a_force, handoff in row_re.findall(text)
    ]
    if not rows:
        sys.exit("could not parse U/a_z sweep rows from fluidized_bed_umf output")

    return {
        "passed": "VALIDATION: PASS" in text,
        "rows": rows,
        "u_wy": grab(text, r"U_mf Wen&Yu \(1966\) REFERENCE:\s+([0-9.]+)"),
        "u_seam": grab(text, r"U_mf MEASURED .*:\s+([0-9.]+)"),
        "u_dyn": grab(text, r"U_mf DYNAMIC .*:\s+([0-9.]+)"),
        "rel_err": grab(text, r"rel\.err\s+([0-9.]+)%"),
        "tol": grab(text, r"rel\.err\s+[0-9.]+%\s+\(tol\s+([0-9.]+)%\)") / 100.0,
        "u_nopg": grab(text, r"omit-∇P U_mf\s+([0-9.]+)"),
        "u_epsbug": grab(text, r"eps-power-bug U_mf\s+([0-9.]+)"),
        "handoff_worst": grab(text, r"worst \|a_z_meas .* =\s+([0-9.]+)"),
    }


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = parse(run())
    rows = data["rows"]
    u = [row["u"] for row in rows]
    a_meas = [row["a_meas"] for row in rows]
    a_force = [row["a_force"] for row in rows]
    u_wy = data["u_wy"]
    u_seam = data["u_seam"]
    u_dyn = data["u_dyn"]
    tol = data["tol"]

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.plot(u, a_meas, "o-", lw=2.0, ms=5, color="#1f77b4", label="integrated bed a_z")
    ax.plot(u, a_force, "s--", lw=1.5, ms=4, color="#2ca02c", label="net seam force / M_bed")
    ax.axhline(0.0, color="black", lw=1.0)
    ax.axvspan(
        u_wy * (1.0 - tol),
        u_wy * (1.0 + tol),
        color="0.88",
        label=f"Wen & Yu ±{100.0 * tol:.0f}% gate",
    )
    ax.axvline(u_wy, color="black", ls="--", lw=1.4, label="Wen & Yu reference")
    ax.axvline(u_seam, color="#d62728", lw=1.7, label="live seam U_mf")
    ax.plot([u_dyn], [0.0], marker="x", ms=9, mew=2, color="#d62728", label="dynamic zero crossing")
    ax.axvline(data["u_nopg"], color="#7f3c8d", ls=":", lw=1.4, label="omit-grad-P control")
    ax.axvline(data["u_epsbug"], color="#ff7f0e", ls=":", lw=1.4, label="epsilon-power control")
    ax.set_xlabel("superficial velocity U [m/s]")
    ax.set_ylabel("bed acceleration a_z [m/s^2]")
    status = "PASS" if data["passed"] else "FAIL"
    ax.set_title(
        "fluidized_bed_umf: live DEM-CFD U_mf vs Wen & Yu\n"
        f"{status}: seam {u_seam:.4f} m/s vs reference {u_wy:.4f} m/s "
        f"({data['rel_err']:.2f}% error)"
    )
    ax.grid(True, alpha=0.28)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()

    outdir = HERE / "plots"
    outdir.mkdir(exist_ok=True)
    dst = outdir / "umf_wen_yu.png"
    fig.savefig(dst, dpi=150)
    print(
        f"wrote {dst}  ({status}: U_mf {u_seam:.4f} vs Wen&Yu {u_wy:.4f}, "
        f"err {data['rel_err']:.2f}%, handoff {data['handoff_worst']:.4f})"
    )


if __name__ == "__main__":
    main()
