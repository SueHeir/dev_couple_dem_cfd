#!/usr/bin/env python3
"""Run the adaptive UMF strategy example and plot the validation metrics."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
EXAMPLE = Path(__file__).resolve().parent
CONFIG = EXAMPLE / "config.toml"
PLOTS = EXAMPLE / "plots"


def run_example() -> list[dict[str, str]]:
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "adaptive_umf_strategy",
        "--",
        str(CONFIG.relative_to(ROOT)),
    ]
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False)
    (PLOTS / "adaptive_umf_strategy.log").write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        raise SystemExit(proc.returncode)

    rows: list[dict[str, str]] = []
    for line in proc.stdout.splitlines():
        if not line.startswith("CSV,"):
            continue
        _, strategy, kind, u, re, a_meas, a_force, residual, substeps, accepted = next(
            csv.reader([line])
        )
        rows.append(
            {
                "strategy": strategy,
                "kind": kind,
                "U": u,
                "Re": re,
                "a_meas": a_meas,
                "a_force": a_force,
                "residual": residual,
                "substeps": substeps,
                "accepted": accepted,
            }
        )
    if not rows:
        raise SystemExit("example produced no CSV rows")
    return rows


def write_csv(rows: list[dict[str, str]]) -> Path:
    out = PLOTS / "adaptive_umf_strategy.csv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return out


def zero_crossing(xs: list[float], ys: list[float]) -> float | None:
    for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
        if y0 <= 0.0 < y1:
            return x0 + (x1 - x0) * (-y0) / (y1 - y0)
    return None


def plot(rows: list[dict[str, str]]) -> Path:
    PLOTS.mkdir(parents=True, exist_ok=True)
    by_strategy: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_strategy.setdefault(row["strategy"], []).append(row)

    u_ref = 0.53804
    tol = 0.15
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)

    for strategy, vals in by_strategy.items():
        vals = sorted(vals, key=lambda r: float(r["U"]))
        us = [float(r["U"]) for r in vals]
        ameas = [float(r["a_meas"]) for r in vals]
        residual = [float(r["residual"]) for r in vals]
        label = strategy.replace("_", " ")
        ax0.plot(us, ameas, marker="o", label=label)
        ax1.plot(us, residual, marker="o", label=label)
        u0 = zero_crossing(us, ameas)
        if u0 is not None:
            ax0.plot([u0], [0.0], marker="x", markersize=8, color=ax0.lines[-1].get_color())

    ax0.axhline(0.0, color="black", linewidth=1)
    ax0.axvline(u_ref, color="black", linestyle="--", linewidth=1, label="Wen & Yu U_mf")
    ax0.axvspan(u_ref * (1.0 - tol), u_ref * (1.0 + tol), color="0.85", label="15% U_mf band")
    ax0.set_ylabel("bed acceleration a_z [m/s^2]")
    ax0.set_title("Measured dynamic onset and coupling residual")
    ax0.legend(loc="upper left", fontsize=8)

    ax1.axhline(0.015, color="black", linestyle="--", linewidth=1, label="adaptive residual gate")
    ax1.set_xlabel("superficial velocity U [m/s]")
    ax1.set_ylabel("|a_meas - a_force| / g")
    ax1.legend(loc="upper left", fontsize=8)
    ax1.grid(True, alpha=0.25)
    ax0.grid(True, alpha=0.25)

    fig.tight_layout()
    out = PLOTS / "adaptive_umf_strategy.png"
    fig.savefig(out, dpi=180)
    return out


def main() -> None:
    PLOTS.mkdir(parents=True, exist_ok=True)
    rows = run_example()
    csv_path = write_csv(rows)
    png_path = plot(rows)
    print(f"wrote {csv_path}")
    print(f"wrote {png_path}")


if __name__ == "__main__":
    main()
