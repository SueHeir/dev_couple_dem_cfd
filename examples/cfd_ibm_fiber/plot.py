#!/usr/bin/env python3
"""Run cfd_ibm_fiber and plot the measured validation gates.

    $BENCH_PYTHON examples/cfd_ibm_fiber/plot.py

The figure is regenerated from the example's own stdout: Archimedes buoyancy
error from the exact cut-cell load check, and the measured drag anisotropy ratio
against Tirado-Garcia de la Torre with the tolerance band and isotropic null.
"""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CFG = "examples/cfd_ibm_fiber/config.toml"


def run():
    out = subprocess.run(
        ["cargo", "run", "--release", "--example", "cfd_ibm_fiber", "--", CFG],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    text = out.stdout + out.stderr
    if out.returncode != 0:
        print(text)
        sys.exit(out.returncode)
    return text


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    text = run()
    data = parse(text)

    fig, (ax_r, ax_b) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.0),
        gridspec_kw={"height_ratios": [1.35, 1.0]},
    )

    re_vals = [re for re, _ in data["ratios"]]
    r_vals = [r for _, r in data["ratios"]]
    ax_r.plot(re_vals, r_vals, "o-", lw=1.8, ms=5, color="#1f77b4", label="measured IBM fiber")
    ax_r.axhline(data["r_theory"], color="#2ca02c", lw=1.5, label=f"Tirado-GdlT {data['r_theory']:.4f}")
    ax_r.axhspan(
        data["r_theory"] * (1.0 - data["tol_ratio"]),
        data["r_theory"] * (1.0 + data["tol_ratio"]),
        color="#2ca02c",
        alpha=0.14,
        label=f"+/- {100.0 * data['tol_ratio']:.0f}% pass band",
    )
    ax_r.axhline(1.0, color="#d62728", ls="--", lw=1.4, label="isotropic control r=1")
    ax_r.set_xlabel("diameter Reynolds number")
    ax_r.set_ylabel("drag anisotropy ratio  F_perp / F_parallel")
    ax_r.set_title(
        "cfd_ibm_fiber drag anisotropy - "
        f"{'PASS' if data['passed'] else 'FAIL'} "
        f"(lowest Re err {100.0 * data['rel_low']:.2f}%)"
    )
    ax_r.grid(True, alpha=0.28)
    ax_r.legend(loc="best", fontsize=8)

    labels = ["fiber clump", "single-bead control"]
    vals = [100.0 * data["buoy_err"], 100.0 * data["single_err"]]
    colors = ["#1f77b4", "#d62728"]
    ax_b.bar(labels, vals, color=colors, alpha=0.86)
    ax_b.axhline(100.0 * data["tol_buoyancy"], color="#2ca02c", lw=1.5, label="Archimedes tolerance")
    ax_b.set_ylabel("buoyancy error vs Archimedes [%]")
    ax_b.set_title(
        "exact cut-cell load check - "
        f"uniform-pressure residual {data['uniform_ratio']:.1e}"
    )
    ax_b.grid(True, axis="y", alpha=0.28)
    ax_b.legend(loc="upper left", fontsize=8)
    ymax = max(vals + [100.0 * data["tol_buoyancy"]]) * 1.18
    ax_b.set_ylim(0, ymax)
    for i, v in enumerate(vals):
        ax_b.text(i, v + ymax * 0.025, f"{v:.2f}%", ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    outdir = HERE / "plots"
    outdir.mkdir(exist_ok=True)
    dst = outdir / "fiber_validation.png"
    fig.savefig(dst, dpi=140)
    print(f"wrote {dst}  ({'PASS' if data['passed'] else 'FAIL'})")


def parse(text):
    passed = "VALIDATION: PASS" in text
    ratios = [
        (float(re_s), float(r_s))
        for re_s, r_s in re.findall(r"#\s+r\(Re_d=([0-9.]+)\)\s+=\s+([0-9.]+)", text)
    ]
    # Keep the final result summary if the progress section also matched.
    if len(ratios) > 1:
        seen = {}
        for re_v, r_v in ratios:
            seen[re_v] = r_v
        ratios = sorted(seen.items())
    if not ratios:
        sys.exit("could not parse drag ratio summary from cfd_ibm_fiber output")

    r_theory = grab(text, r"r_theory .*?:\s+([0-9.]+)")
    rel_low = grab(text, r"lowest-Re rel\.err:\s+([0-9.]+)%") / 100.0
    tol_ratio = grab(text, r"lowest-Re rel\.err:.*\(tol\s+([0-9.]+)%\)") / 100.0
    buoy_err = grab(text, r"hydrostatic.*err\s+([0-9.]+)%\s+(?:<=|≤)") / 100.0
    tol_buoyancy = grab(text, r"hydrostatic.*err\s+[0-9.]+%\s+(?:<=|≤)\s+([0-9.]+)%\)") / 100.0
    uniform_ratio = grab(text, r"uniform pressure.*=\s+([0-9.eE+-]+)\s+\((?:<=|≤)")
    single_err = grab(text, r"SINGLE-bead integrator.*\(([0-9.]+)% off\)") / 100.0
    return {
        "passed": passed,
        "ratios": ratios,
        "r_theory": r_theory,
        "rel_low": rel_low,
        "tol_ratio": tol_ratio,
        "buoy_err": buoy_err,
        "tol_buoyancy": tol_buoyancy,
        "uniform_ratio": uniform_ratio,
        "single_err": single_err,
    }


def grab(text, pat):
    m = re.search(pat, text)
    if not m:
        sys.exit(f"could not parse pattern: {pat}")
    return float(m.group(1))


if __name__ == "__main__":
    main()
