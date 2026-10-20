#!/usr/bin/env python3
"""Regenerate Supplementary Fig. S6 (thread scaling) from the Ice Lake sweep.

The original `figures/fig3/repeat_speedup.pdf` was produced in R 4.4.3 by a
script that is NOT in the paper repository (it arrived as a binary in 178a579,
"add Jarno's updates"). `plot_scaling.py` did not make it either -- that is a
matplotlib script whose outputs (`scaling.svg`/`scaling.png`) sit unused beside
it. This is therefore a re-implementation of the figure's STRUCTURE, from the
published caption, on measured data.

Data: figS6_phase_timings.tsv, the warm runs of thread-sweep job 38373939
(2026-10-12; $KS_WORK/ks_slurm/threadsweep_cold/, file
figS6_phase_timings.warm.tsv there) -- KaryoScope v2.3.0, hks 0.4.0, HG002 hap1
(2.95 Gb), `repeat` feature set, 1-32 threads, n=3, Intel Xeon Gold 6342
(Ice Lake, 2.80 GHz), the same hardware as Table S1. One job on one node with
no other tenant; the 36 points ran in a fixed random order (seed 20261012).
For each point the database, input, binaries and source were evicted from the
OS page cache (posix_fadvise DONTNEED, residency checked with fincore), then
the command ran twice: a cold run, and the warm run plotted here, with the
cache as the cold run left it. The cold runs (figS6_phase_timings.cold.tsv)
are not plotted; Supplementary Note 2 reports their extra assignment time.
This replaces the earlier sweep (jobs 35376001 / 35376043), whose cache state
was confounded with thread count and run order.

Two deliberate departures from the original, both forced by what is measurable:

  1. The original's bottom row splits peak memory by phase (assignment vs
     smoothing). Our logs record peak RSS for the whole `annotate` process only,
     so memory is shown ONCE, spanning both columns, and labelled as such.
     Inventing a per-phase split would be fabrication.
  2. The original describes the Rust smoother as single-threaded, drawn as a
     horizontal reference. In v2.3.0 it is threaded and does scale
     (42.2 s -> 13.8 s over 1 -> 8 threads, warm), so it is plotted as a curve.

Text is emitted as text, not outlined paths (`svg.fonttype = "none"`), so every
label stays editable in Illustrator.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

# Caption's colour convention: KMC red, HKS blue.
STYLE = {
    "kmc": dict(color="#d62728", marker="o", label="KaryoScope-KMC"),
    "hks": dict(color="#1f77b4", marker="s", label="KaryoScope-HKS"),
}
SMOOTH_LABEL = {"kmc": "KMC (Python smoother)", "hks": "HKS (Rust smoother)"}


def style_axes(ax, ticks):
    # Linear x-axis with ticks at the measured thread counts, matching the
    # original figure: positions are proportional to thread count, so 1-2-4 sit
    # close together and the gridlines fan out toward 32.
    ax.set_xscale("linear")
    ax.set_xlim(0, max(ticks) * 1.03)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.grid(True, which="major", alpha=0.25, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", type=Path,
                    default=Path(__file__).resolve().parent / "figS6_phase_timings.tsv")
    ap.add_argument("--out", type=Path,
                    default=Path("figS6_regenerated.svg"))
    ap.add_argument("--summary-tsv", type=Path,
                    default=Path("figS6_summary.regenerated.tsv"))
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df["smooth_s"] = pd.to_numeric(df["smooth_s"], errors="coerce")

    # 24 threads is dropped from the figure: it is the only non-power-of-two
    # point, and every difference between 16, 24 and 32 is inside replicate
    # scatter, so the extra point adds no resolvable information.
    df = df[df["threads"] != 24]

    g = (df.groupby(["backend", "threads"])
           .agg(n=("rep", "count"),
                assign_mb_s=("assign_Mb_s", "mean"), assign_sd=("assign_Mb_s", "std"),
                smooth_s=("smooth_s", "mean"),      smooth_sd=("smooth_s", "std"),
                rss=("peak_rss_GiB", "mean"),       rss_sd=("peak_rss_GiB", "std"))
           .reset_index().sort_values(["backend", "threads"]))
    args.summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(args.summary_tsv, sep="\t", index=False, float_format="%.3f")

    plt.rcParams.update({
        "svg.fonttype": "none",          # keep text editable in Illustrator
        # With svg.fonttype="none" the SVG only NAMES the font, so it must be one
        # the viewer has. matplotlib's default (DejaVu Sans) is bundled with
        # matplotlib but is not a system font, so Illustrator and browsers
        # substituted it and the labels came out wrong. Helvetica/Arial are
        # present on macOS and match the other manuscript figures.
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Helvetica Neue",
                            "Liberation Sans", "DejaVu Sans"],
        # Mathtext ($k$-mer in a title, 10^n tick labels) is typeset with the
        # mathtext font set, which ignores font.sans-serif and defaults to
        # DejaVu. Point it at the same family so no DejaVu name reaches the SVG.
        # Arial, not Helvetica: macOS ships Helvetica as a .ttc collection that
        # matplotlib cannot index by style, so "Helvetica:italic" silently
        # resolves to the regular face and the italic k in "$k$-mer" comes out
        # upright. Arial has a separate italic file and is metrically compatible.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
        "legend.fontsize": 9, "figure.dpi": 150,
    })

    fig = plt.figure(figsize=(9, 7.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.85], hspace=0.32, wspace=0.26)
    ax_assign = fig.add_subplot(gs[0, 0])
    ax_smooth = fig.add_subplot(gs[0, 1])
    ax_mem = fig.add_subplot(gs[1, :])
    ticks = sorted(df["threads"].unique().tolist())

    for b in ("kmc", "hks"):
        s = g[g["backend"] == b]
        st = STYLE[b]
        ax_assign.errorbar(s["threads"], s["assign_mb_s"], yerr=s["assign_sd"].fillna(0),
                           capsize=3, linewidth=1.5, markersize=6, **st)
        ax_smooth.errorbar(s["threads"], s["smooth_s"], yerr=s["smooth_sd"].fillna(0),
                           capsize=3, linewidth=1.5, markersize=6,
                           color=st["color"], marker=st["marker"], label=SMOOTH_LABEL[b])
        ax_mem.errorbar(s["threads"], s["rss"], yerr=s["rss_sd"].fillna(0),
                        capsize=3, linewidth=1.5, markersize=6, **st)

    ax_assign.set_title("$k$-mer feature assignment")
    ax_assign.set_ylabel("Throughput (Mb/s)")
    ax_assign.set_xlabel("Threads")
    ax_assign.legend(frameon=False, loc="upper left")
    style_axes(ax_assign, ticks)

    ax_smooth.set_title("Hierarchy-aware smoothing")
    ax_smooth.set_ylabel("Wall-clock (s)")
    ax_smooth.set_xlabel("Threads")
    ax_smooth.set_yscale("log", base=10)
    # Plain numeric tick labels rather than matplotlib's 10^n mathtext. Mathtext
    # is typeset in the mathtext font, which defaults to DejaVu regardless of
    # font.sans-serif, and with svg.fonttype="none" that name leaks into the SVG
    # and gets substituted by viewers that lack it.
    ax_smooth.set_yticks([20, 50, 100, 200, 500, 1000])
    ax_smooth.set_yticklabels(["20", "50", "100", "200", "500", "1000"])
    ax_smooth.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax_smooth.legend(frameon=False, loc="upper right")
    style_axes(ax_smooth, ticks)

    ax_mem.set_title("Peak resident memory (whole annotate process)")
    ax_mem.set_ylabel("Peak RSS (GiB)")
    ax_mem.set_xlabel("Threads")
    ax_mem.set_ylim(0, g["rss"].max() * 1.15)
    ax_mem.legend(frameon=False, loc="center right")
    style_axes(ax_mem, ticks)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight", format="svg")
    fig.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(args.out.with_suffix(".png"), bbox_inches="tight", dpi=300)
    print(f"[ok] {args.out}  (+ .pdf, .png)")
    print(f"[ok] {args.summary_tsv}")

    print("\nProse-ready:")
    for b in ("kmc", "hks"):
        s = g[g["backend"] == b].set_index("threads")
        print(f"  {b}: assignment {s.loc[1,'assign_mb_s']:.2f} -> "
              f"{s.loc[32,'assign_mb_s']:.2f} Mb/s (1->32 threads); "
              f"smoothing {s.loc[1,'smooth_s']:.0f} -> {s.loc[16,'smooth_s']:.0f} s (1->16); "
              f"peak RSS {s.loc[16,'rss']:.2f} GiB at 16")
    k16 = g[(g.backend=='kmc') & (g.threads==16)].iloc[0]
    h16 = g[(g.backend=='hks') & (g.threads==16)].iloc[0]
    print(f"  HKS/KMC at 16 threads: assignment {h16.assign_mb_s/k16.assign_mb_s:.2f}x throughput, "
          f"smoothing {k16.smooth_s/h16.smooth_s:.1f}x faster, "
          f"{k16.rss/h16.rss:.2f}x less memory")


main()
