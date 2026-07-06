#!/usr/bin/env python3
"""Run fixed_bed_ergun and plot the Ergun validation gate.

    $BENCH_PYTHON examples/fixed_bed_ergun/plot.py

The figure is regenerated from the example's stdout: MacDonald-seam pressure-drop
error versus the Ergun reference across the Reynolds sweep, the 25% pass line,
and the deliberately corrupted epsilon-power negative control.
"""
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CFG = "examples/fixed_bed_ergun/config.toml"


def run():
    out = subprocess.run(
        ["cargo", "run", "--release", "--example", "fixed_bed_ergun", "--", CFG],
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

    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    ax.semilogx(
        data["rep"],
        [100.0 * e for e in data["rel_err"]],
        "o-",
        lw=1.8,
        ms=5,
        color="#1f77b4",
        label="MacDonald seam vs Ergun reference",
    )
    ax.axhline(
        100.0 * data["tol"],
        color="#d62728",
        ls="--",
        lw=1.6,
        label=f"pass tolerance ({100.0 * data['tol']:.0f}%)",
    )
    ax.axhline(
        100.0 * data["neg_control"],
        color="#7f3c8d",
        ls=":",
        lw=1.8,
        label=f"corrupted seam worst ({100.0 * data['neg_control']:.1f}%)",
    )
    ax.plot(
        [max(data["rep"])],
        [100.0 * data["neg_control"]],
        "o",
        ms=6,
        color="#7f3c8d",
    )
    ax.annotate(
        "negative control fails",
        xy=(max(data["rep"]), 100.0 * data["neg_control"]),
        xytext=(data["rep"][-3], 170.0),
        arrowprops={"arrowstyle": "-", "color": "#7f3c8d"},
        color="#7f3c8d",
        fontsize=9,
    )
    ax.set_xlabel("modified particle Reynolds number  Re_p")
    ax.set_ylabel("pressure-drop relative error vs Ergun [%]")
    ax.set_title(
        "fixed_bed_ergun: DEM-CFD seam pressure drop vs Ergun\n"
        f"{'PASS' if data['passed'] else 'FAIL'}: measured error stays below tolerance; "
        "corrupted seam fails"
    )
    ax.grid(True, which="both", alpha=0.28)
    ax.legend(loc="center right", fontsize=9)
    ax.set_ylim(0.0, max(210.0, 100.0 * data["neg_control"] * 1.08))
    fig.tight_layout()

    outdir = HERE / "plots"
    outdir.mkdir(exist_ok=True)
    dst = outdir / "ergun_relative_error.png"
    fig.savefig(dst, dpi=140)
    print(f"wrote {dst}  ({'PASS' if data['passed'] else 'FAIL'})")


def parse(text):
    row_re = re.compile(
        r"^\s*([0-9.]+)\s+([0-9.]+)\s+"
        r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)%\s+",
        re.MULTILINE,
    )
    rows = row_re.findall(text)
    if not rows:
        sys.exit("could not parse Reynolds sweep rows from fixed_bed_ergun output")
    rep = [float(row[1]) for row in rows]
    rel_err = [float(row[4]) / 100.0 for row in rows]
    tol = grab(text, r"Ergun rel\.err spread.*\(tol\s+([0-9.]+)%\)") / 100.0
    neg_control = grab(text, r"negative control.*worst rel\.err\s+([0-9.]+)%") / 100.0
    return {
        "passed": "VALIDATION: PASS" in text,
        "rep": rep,
        "rel_err": rel_err,
        "tol": tol,
        "neg_control": neg_control,
    }


def grab(text, pat):
    m = re.search(pat, text)
    if not m:
        sys.exit(f"could not parse pattern: {pat}")
    return float(m.group(1))


if __name__ == "__main__":
    main()
