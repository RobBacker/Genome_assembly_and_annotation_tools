#!/usr/bin/env python3
"""
# This is not production code, but included for publication posterity.

PAV figure builder for orthogroup presence/absence matrices.

Inputs
------
TSV with columns:
  Orthogroup, <sample1>, <sample2>, ..., [count], [category]
Sample columns are 0/1 integers. 'count' and 'category' are optional.

What this script does
---------------------
- Optional collapse of primary/alternate haplotypes.
- Honors existing 'category' unless --recompute-category is set.
- Enforces user sample order across panels (A, B, C, D).
- Panel A: custom UpSet-like plot with
    * exact selection: 1 Core, 15 Soft-core, 15 Shell, 15 Cloud (if available)
    * colored bars by category (your palette)
    * broken (split) y-axis auto-derived from data (to separate the huge Core bar)
    * dot-matrix underneath + group separators & labels
    * **percentage-driven grouping** (default): Core ≥95%, Soft-core ≥85%, Shell ≥20%, Cloud <20%
- Panel B: category composition pie (colored by palette)
- Panel C: rotated 90° heatmap (samples on Y), sorted for clean blocks; vertical category separators + colored header strip; Core/Soft-core compressed in width
- Panel D: stacked bars per accession, in your order but vertically flipped

Usage
-----
python upset_orthogroups.py \
  --input PAV_orthogroups.tsv \
  --outdir out_figs \
  --sample-order "West-Indian pure genome,Ashdot,Leola,Mike,Dusa,Choquette,Gottfried,Hass" \
  --collapse-haplotypes

Dependencies
------------
pip install pandas numpy matplotlib
"""

from __future__ import annotations
import argparse
import os
import re
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# -----------------------------
# Colors (mapped to PAV classes)
# -----------------------------

PAV_COLORS: Dict[str, str] = {
    "Core": "#1b9e77",
    "Soft-core": "#d95f02",
    "Shell": "#7570b3",
    "Cloud": "#e7298a",
}
CATEGORY_ORDER = ["Core", "Soft-core", "Shell", "Cloud"]

# Panel C width scaling (fraction of columns kept per category)
PANELC_WIDTH_SCALE = {"Core": 0.10, "Soft-core": 0.50, "Shell": 1.00, "Cloud": 1.00}

# Aliases help reconcile label variants in --sample-order
ALIASES = {
    "WI": "West-Indian pure genome",
    "West-Indian": "West-Indian pure genome",
    "FABI-033": "Ashdot",
    "Leola™": "Leola",
    "Dusa®": "Dusa",
}


# -----------------------------
# CLI
# -----------------------------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Build PAV figure panels from an orthogroup presence/absence matrix."
    )
    ap.add_argument("--input", required=True, help="Orthogroup PAV TSV")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--collapse-haplotypes", action="store_true",
                    help="Collapse *_primary/_alternate into a single accession (present if either is 1)")
    ap.add_argument("--ascii-labels", action="store_true",
                    help="Sanitize column labels to ASCII for plotting")
    ap.add_argument("--softcore-miss", type=int, default=1,
                    help="Soft-core (re)computation = present in N - {this many misses}. E.g., 1 => N-1")
    ap.add_argument("--cloud-max", type=int, default=1,
                    help="Cloud (re)computation if present in <= this many samples")
    ap.add_argument("--recompute-category", action="store_true",
                    help="Ignore existing 'category' and recompute after any collapsing")
    ap.add_argument("--dpi", type=int, default=600, help="Figure DPI")
    ap.add_argument(
        "--sample-order",
        type=str,
        default=None,
        help=('Comma-separated base names to enforce column order across panels, '
              'e.g. "West-Indian pure genome,Ashdot,Leola,Mike,Dusa,Choquette,Gottfried,Hass"')
    )
    ap.add_argument(
        "--haplotype-order",
        type=str,
        default=None,
        help=("If not collapsed, order within base (comma-separated), e.g. 'primary,alternate' "
              "or 'alternate,primary'. Default = keep file order per base.")
    )
    # --- Panel A percentage thresholds (ON by default) ---
    ap.add_argument("--panelA-use-percent", action="store_true", default=True,
                    help="Group Panel A intersections by percentage thresholds (default True)")
    ap.add_argument("--pct-core", type=float, default=0.95,
                    help="Core threshold as fraction (>=), e.g., 0.95")
    ap.add_argument("--pct-softcore", type=float, default=0.85,
                    help="Soft-core lower bound as fraction (>= and < Core), e.g., 0.85")
    ap.add_argument("--pct-shell", type=float, default=0.20,
                    help="Shell lower bound as fraction (>= and < Soft-core), e.g., 0.20; Cloud = < this")
    return ap.parse_args()


# -----------------------------
# Utilities
# -----------------------------
def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def sanitize_label(s: str) -> str:
    s = s.replace("®", "R").replace("™", "TM")
    s = re.sub(r"[_\s]+$", "", s)
    try:
        s = s.encode("ascii", "ignore").decode()
    except Exception:
        pass
    return s


def normalize_base(s: str) -> str:
    """Strip haplotype suffix and normalize aliases for matching."""
    s0 = re.sub(r"_(primary|alternate)$", "", s)
    s0 = s0.strip().replace("™", "").replace("®", "")
    return ALIASES.get(s0, s0)


def infer_sample_columns(df: pd.DataFrame) -> List[str]:
    excluded = {"Orthogroup", "count", "category"}
    # preserve file order
    return [c for c in df.columns if c not in excluded]


def collapse_haplotypes(
    df: pd.DataFrame, sample_cols: List[str]
) -> Tuple[pd.DataFrame, List[str], Dict[str, List[str]]]:
    """
    Collapse *_primary/_alternate => base name (present if either is 1).
    Returns (collapsed_df, new_cols, mapping_dict)
    """
    groups: Dict[str, List[str]] = {}
    for c in sample_cols:
        m = re.match(r"(.+?)_(primary|alternate)$", c)
        base = normalize_base(m.group(1) if m else c)
        groups.setdefault(base, []).append(c)

    new_df = df.copy()
    new_cols: List[str] = []
    for base, cols in groups.items():
        if len(cols) == 1:
            new_df[base] = new_df[cols[0]].astype(int)
        else:
            new_df[base] = (new_df[cols].max(axis=1)).astype(int)
        new_cols.append(base)

    keep_cols = ["Orthogroup"] + new_cols
    if "count" in new_df.columns:
        keep_cols.append("count")
    if "category" in new_df.columns:
        keep_cols.append("category")
    new_df = new_df[keep_cols]
    return new_df, new_cols, groups


def reorder_sample_cols(
    sample_cols: List[str],
    sample_order: Optional[str],
    haplotype_order: Optional[str],
) -> List[str]:
    """
    Reorder columns to follow a user-specified base order.
    If haplotypes exist and not collapsed, keep per-base haplotype order as in file,
    unless --haplotype-order is given.
    """
    if not sample_order:
        return sample_cols

    desired_bases = [normalize_base(x.strip()) for x in sample_order.split(",") if x.strip()]
    desired_hap = None
    if haplotype_order:
        desired_hap = [x.strip() for x in haplotype_order.split(",") if x.strip()]
        if set(desired_hap) - {"primary", "alternate"}:
            raise SystemExit("--haplotype-order must be a comma list of 'primary' and/or 'alternate'")

    groups: Dict[str, List[str]] = {}
    for c in sample_cols:
        m = re.match(r"(.+?)_(primary|alternate)$", c)
        base = normalize_base(m.group(1) if m else c)
        groups.setdefault(base, []).append(c)

    # within each base, apply desired hap order
    if desired_hap:
        for base, cols in groups.items():
            ordered = []
            for h in desired_hap:
                ordered.extend([c for c in cols if c.endswith("_" + h)])
            ordered.extend([c for c in cols if c not in ordered])
            groups[base] = ordered

    final_cols: List[str] = []
    seen = set()
    for base in desired_bases:
        if base in groups:
            final_cols.extend(groups[base])
            seen.add(base)
    # append leftovers in original relative order
    for c in sample_cols:
        b = normalize_base(re.sub(r"_(primary|alternate)$", "", c))
        if b not in seen and b in groups:
            final_cols.extend(groups[b])
            seen.add(b)

    # dedupe preserving order
    seen_cols = set()
    final_cols = [c for c in final_cols if not (c in seen_cols or seen_cols.add(c))]
    return final_cols


def compute_count_and_category(
    df: pd.DataFrame,
    sample_cols: List[str],
    softcore_miss: int = 1,
    cloud_max: int = 1,
    force_recompute_category: bool = False,
) -> pd.DataFrame:
    """
    Always recompute 'count' from the CURRENT matrix.
    Use existing 'category' unless --recompute-category is set or it's missing.
    """
    N = len(sample_cols)
    out = df.copy()
    out["count"] = out[sample_cols].sum(axis=1).astype(int)

    if "category" in out.columns and not force_recompute_category:
        out["category"] = (
            out["category"]
            .astype(str)
            .str.strip()
            .replace({"Dispensable": "Shell", "Specific": "Cloud"})
        )
        out["category"] = pd.Categorical(out["category"], categories=CATEGORY_ORDER, ordered=True)
        return out

    core_mask = out["count"] == N
    soft_mask = (out["count"] >= (N - softcore_miss)) & (out["count"] < N)
    cloud_mask = (out["count"] >= 1) & (out["count"] <= cloud_max)
    cat = np.where(core_mask, "Core",
          np.where(soft_mask, "Soft-core",
          np.where(cloud_mask, "Cloud", "Shell")))
    out["category"] = pd.Categorical(cat, categories=CATEGORY_ORDER, ordered=True)
    return out


# -----------------------------
# Intersections for Panel A
# -----------------------------
def pattern_counts_boolean(df_bool: pd.DataFrame) -> pd.DataFrame:
    """Return DataFrame with boolean columns for each sample, plus 'size' and 'degree'."""
    s = df_bool.groupby(list(df_bool.columns), observed=False).size()
    s.name = "size"
    df = s.reset_index()
    df["degree"] = df[df_bool.columns].astype(int).sum(axis=1)
    df = df.sort_values("size", ascending=False).reset_index(drop=True)
    return df


def select_panelA_intersections(
    pat_df: pd.DataFrame,
    sample_cols: List[str],
    softcore_miss: int,
    cloud_max: int,
    core_n: int = 1,
    soft_n: int = 15,
    shell_n: int = 15,
    cloud_n: int = 15,
    use_percent: bool = True,
    pct_core: float = 0.95,
    pct_softcore: float = 0.85,
    pct_shell: float = 0.20,
) -> Tuple[pd.DataFrame, List[int], List[str]]:
    """
    Select exact sets for Panel A, either by percentage thresholds (default) or absolute ranges.

    Percent mode (fractions of N):
      Core:       deg/N >= pct_core
      Soft-core:  pct_softcore <= deg/N < pct_core
      Shell:      pct_shell    <= deg/N < pct_softcore
      Cloud:      deg/N < pct_shell
    """
    N = len(sample_cols)
    df = pat_df.copy()

    if use_percent:
        core_min  = int(np.ceil(pct_core * N))
        soft_min  = int(np.ceil(pct_softcore * N))
        shell_min = int(np.ceil(pct_shell * N))

        core_mask  = df["degree"] >= core_min
        soft_mask  = (df["degree"] >= soft_min)  & (df["degree"] < core_min)
        shell_mask = (df["degree"] >= shell_min) & (df["degree"] < soft_min)
        cloud_mask = df["degree"] < shell_min
    else:
        # absolute-degree mode (misses / cloud_max)
        core_mask  = df["degree"] == N
        soft_mask  = (df["degree"] >= (N - softcore_miss)) & (df["degree"] < N)
        cloud_mask = (df["degree"] >= 1) & (df["degree"] <= cloud_max)
        shell_mask = (df["degree"] >= (cloud_max + 1)) & (df["degree"] <= (N - softcore_miss - 1))

    blocks, groups = [], []

    core_block = df[core_mask].sort_values("size", ascending=False).head(core_n).copy()
    if len(core_block): core_block["group"] = "Core";       blocks.append(core_block); groups.append("Core")

    soft_block = df[soft_mask].sort_values("size", ascending=False).head(soft_n).copy()
    if len(soft_block): soft_block["group"] = "Soft-core";  blocks.append(soft_block); groups.append("Soft-core")

    shell_block = df[shell_mask].sort_values("size", ascending=False).head(shell_n).copy()
    if len(shell_block): shell_block["group"] = "Shell";     blocks.append(shell_block); groups.append("Shell")

    cloud_block = df[cloud_mask].sort_values("size", ascending=False).head(cloud_n).copy()
    if len(cloud_block): cloud_block["group"] = "Cloud";     blocks.append(cloud_block); groups.append("Cloud")

    if not blocks:
        return df.head(0), [], []

    sel_df = pd.concat(blocks, axis=0, ignore_index=True)
    sel_df = sel_df[[*sample_cols, "size", "degree", "group"]]

    boundaries, start = [], 0
    for b in blocks:
        end = start + len(b)
        boundaries.append(end)
        start = end
    return sel_df, boundaries, groups


# -----------------------------
# Panel A: custom plot
# -----------------------------
def broken_ylim_from_sizes(sizes: np.ndarray) -> Optional[Tuple[Tuple[float, float], Tuple[float, float]]]:
    """Return ((ymin_top, ymax_top), (ymin_bot, ymax_bot)) or None if no break needed."""
    if len(sizes) == 0:
        return None
    smax = sizes.max()
    others = sizes[sizes < smax]
    if len(others) == 0:
        return None
    s2 = others.max()
    # Break only if the max is clearly larger than the rest
    if smax <= 1.8 * s2:
        return None
    # Bottom shows up to ~10% above second max; Top shows last 10% of the max bar
    bot_max = max(s2 * 1.10, s2 + 1)
    top_min = max(smax * 0.88, bot_max * 1.02)
    top_max = smax * 1.05
    return ((top_min, top_max), (0, bot_max))


def panel_A_custom(
    df: pd.DataFrame,
    sample_cols: List[str],
    softcore_miss: int,
    cloud_max: int,
    outbase: str,
    dpi: int,
    args: argparse.Namespace,
) -> pd.DataFrame:
    df_bool = df[sample_cols].astype(bool)
    pat = pattern_counts_boolean(df_bool)
    sel, boundaries, groups = select_panelA_intersections(
        pat, sample_cols, softcore_miss, cloud_max,
        core_n=1, soft_n=15, shell_n=15, cloud_n=15,
        use_percent=args.panelA_use_percent,
        pct_core=args.pct_core,
        pct_softcore=args.pct_softcore,
        pct_shell=args.pct_shell,
    )
    if sel.empty:
        print("Panel A: no intersections found to plot.")
        return sel

    sizes = sel["size"].to_numpy()
    colors = [PAV_COLORS[g] for g in sel["group"]]

    # Build layout: top bars (upper range), middle bars (lower range), bottom dot-matrix
    fig = plt.figure(figsize=(16, 8), constrained_layout=True)
    gs = fig.add_gridspec(nrows=3, ncols=1, height_ratios=[1.2, 1.6, 1.4], figure=fig)

    ax_top = fig.add_subplot(gs[0, 0])   # top of split y (shows only top slice of tallest bar)
    ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)  # bottom of split y (shows most bars)
    ax_dot = fig.add_subplot(gs[2, 0], sharex=ax_top)

    x = np.arange(len(sel))
    # Draw bars on both axes; y-lims will clip them appropriately
    for ax in (ax_top, ax_bot):
        ax.bar(x, sizes, color=colors, edgecolor="black", linewidth=0.6)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        ax.margins(x=0.01)

    # Configure split y-axis
    yl = broken_ylim_from_sizes(sizes)
    if yl is None:
        # No split needed: hide top, use bottom axis full range
        ax_top.set_visible(False)
        ax_bot.set_ylim(0, sizes.max() * 1.05)
        ax_bot.set_ylabel("Intersection size")
    else:
        (ymin_top, ymax_top), (ymin_bot, ymax_bot) = yl
        ax_top.set_ylim(ymin_top, ymax_top)
        ax_bot.set_ylim(ymin_bot, ymax_bot)
        ax_bot.set_ylabel("Intersection size")
        # Diagonal break marks
        d = .015
        kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False)
        ax_top.plot((-d, +d), (-d, +d), **kwargs)        # top-left
        ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right
        kwargs.update(transform=ax_bot.transAxes)
        ax_bot.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left
        ax_bot.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right

    # Group separators & labels (vertical lines & text at ax_bot)
    if boundaries:
        for b in boundaries[:-1]:
            ax_top.axvline(b - 0.5, color="k", lw=0.8, alpha=0.5)
            ax_bot.axvline(b - 0.5, color="k", lw=0.8, alpha=0.5)
            ax_dot.axvline(b - 0.5, color="k", lw=0.6, alpha=0.4)

        # group labels centered
        starts = [0] + boundaries[:-1]
        for g, s, e in zip(groups, starts, boundaries):
            mid = (s + e - 1) / 2.0
            ax_bot.text(mid, ax_bot.get_ylim()[1], g, ha="center", va="bottom", fontsize=10)

    # Dot-matrix (samples down Y)
    M = len(sample_cols)
    ax_dot.set_ylim(M - 0.5, -0.5)
    ax_dot.set_yticks(range(M))
    ax_dot.set_yticklabels(sample_cols)
    ax_dot.set_xlabel("Top Intersections (Core | Soft-core | Shell | Cloud)")
    ax_top.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax_bot.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

    mat = sel[sample_cols].astype(int).to_numpy().T  # shape: M x K
    for j in range(mat.shape[1]):
        col = mat[:, j]
        ys = np.where(col == 1)[0]
        ax_dot.scatter(np.full(M, j), range(M), s=18, facecolors="white", edgecolors="gray", linewidths=0.7)
        if len(ys) > 0:
            ax_dot.scatter(np.full(len(ys), j), ys, s=28, facecolors="black", edgecolors="black")
            ax_dot.plot([j, j], [ys.min(), ys.max()], color="black", lw=1.0)

    for ax in (ax_top, ax_bot, ax_dot):
        ax.set_xlim(-0.5, len(sel) - 0.5)

    for ext in ("png", "svg"):
        fig.savefig(f"{outbase}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    sel.to_csv(f"{outbase}.tsv", sep="\t", index=False)
    return sel


# -----------------------------
# Panel B (pie)
# -----------------------------
def panel_B_pie(df: pd.DataFrame, outbase: str, dpi: int) -> pd.DataFrame:
    counts = df["category"].value_counts().reindex(CATEGORY_ORDER).fillna(0).astype(int)
    fig = plt.figure(figsize=(5, 5))
    plt.pie(
        counts,
        labels=[f"{k}" for k in counts.index],
        autopct=lambda p: f"{p:.2f}%",
        startangle=90,
        colors=[PAV_COLORS[k] for k in counts.index],
    )
    plt.title("PAV composition")
    for ext in ("png", "svg"):
        fig.savefig(f"{outbase}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    counts.rename("count").to_csv(f"{outbase}.tsv", sep="\t")
    return counts.to_frame(name="count")


# -----------------------------
# Panel C (rotated 90°, compressed)
# -----------------------------
def _bitrank_frame(df_bin: pd.DataFrame) -> np.ndarray:
    """Rank rows by binary pattern using current column order (leftmost = most significant)."""
    weights = 2 ** np.arange(len(df_bin.columns) - 1, -1, -1, dtype=np.int64)
    return (df_bin.values * weights).sum(axis=1)


def panel_C_rotated(
    df: pd.DataFrame, sample_cols: List[str], outbase: str, dpi: int
) -> None:
    df2 = df.copy()
    # sort: category → count(desc) → bit-pattern(desc)
    df2["bitrank"] = _bitrank_frame(df2[sample_cols].astype(int))
    dsorted = df2.sort_values(by=["category", "count", "bitrank"], ascending=[True, False, False]).reset_index(drop=True)

    # compute per-category counts in this order
    cat_counts = dsorted["category"].value_counts().reindex(CATEGORY_ORDER).fillna(0).astype(int)

    # compress columns (orthogroups) per category for display width
    sel_idx = []
    start = 0
    scaled_counts = []
    for cat in CATEGORY_ORDER:
        n = int(cat_counts.get(cat, 0))
        if n <= 0:
            scaled_counts.append(0)
            continue
        keep = max(1, int(np.ceil(n * PANELC_WIDTH_SCALE.get(cat, 1.0))))
        idx_local = np.unique(np.round(np.linspace(0, n - 1, keep)).astype(int))
        sel_idx.append(np.arange(start, start + n)[idx_local])
        scaled_counts.append(len(idx_local))
        start += n
    if sel_idx:
        sel_idx = np.concatenate(sel_idx)
        dplot = dsorted.iloc[sel_idx]
    else:
        dplot = dsorted.iloc[[]]

    # build matrix (rows=OGs kept, cols=samples) and rotate → rows=samples, cols=OGs
    mat = dplot[sample_cols].to_numpy(dtype=float)
    matT = mat.T

    # figure size: keep reasonable width even for large OG counts
    width = min(14, max(8, 6 + 0.015 * matT.shape[1]))
    height = max(5.5, 0.5 + 0.35 * len(sample_cols))
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_subplot(111)
    im = ax.imshow(matT, aspect="auto", interpolation="nearest", vmin=0, vmax=1)

    ax.set_xticks([])
    ax.set_yticks(np.arange(len(sample_cols)))
    ax.set_yticklabels(sample_cols)
    ax.set_xlabel("Orthogroups (per category)")
    ax.set_ylabel("Accessions")

    # vertical separators & colored header strip using SCALED counts
    cum = np.cumsum(np.array(scaled_counts, dtype=int))
    for x in cum[:-1]:
        ax.axvline(x - 0.5, color="k", lw=0.8, alpha=0.5)

    x0 = 0
    for cat, n in zip(CATEGORY_ORDER, scaled_counts):
        if n <= 0:
            continue
        ax.add_patch(
            plt.Rectangle((x0 - 0.5, -1.2), width=n, height=0.6,
                          facecolor=PAV_COLORS[cat], edgecolor="none", clip_on=False)
        )
        x0 += n

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_ticks([0, 1])
    cbar.ax.set_yticklabels(["Absent", "Present"])

    for ext in ("png", "svg"):
        fig.savefig(f"{outbase}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Panel D (flipped vertically)
# -----------------------------
def panel_D_stacked_bars_flipped(
    df: pd.DataFrame, sample_cols: List[str], outbase: str, dpi: int
) -> pd.DataFrame:
    """Stacked bars per accession; flip order vertically (first in list at bottom)."""
    records = []
    for col in sample_cols:
        core = (df["category"] == "Core").sum()
        soft = ((df["category"] == "Soft-core") & (df[col] == 1)).sum()
        shell = ((df["category"] == "Shell") & (df[col] == 1)).sum()
        cloud = ((df["category"] == "Cloud") & (df[col] == 1)).sum()
        records.append((col, core, soft, shell, cloud))

    res = (
        pd.DataFrame(records, columns=["accession", "Core", "Soft-core", "Shell", "Cloud"])
        .set_index("accession")
    )

    # Flip order
    res = res.iloc[::-1]

    fig = plt.figure(figsize=(10, max(4, 0.4 * len(sample_cols))))
    ax = fig.add_subplot(111)
    left = np.zeros(len(res))
    for cat in CATEGORY_ORDER:
        ax.barh(res.index, res[cat], left=left, label=cat, color=PAV_COLORS[cat])
        left += res[cat].to_numpy()

    ax.set_xlabel("Number of orthogroups")
    ax.set_ylabel("")
    ax.legend(title="Category", bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set_title("Per-accession PAV composition")
    for ext in ("png", "svg"):
        fig.savefig(f"{outbase}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    res.to_csv(f"{outbase}.tsv", sep="\t")
    return res


# -----------------------------
# Main
# -----------------------------
def main():
    args = parse_args()
    ensure_outdir(args.outdir)

    df = pd.read_csv(args.input, sep="\t", dtype={})
    if "Orthogroup" not in df.columns:
        raise SystemExit("Input must have an 'Orthogroup' column.")

    # 1) sample columns in file order
    sample_cols = infer_sample_columns(df)

    # 2) optionally collapse haplotypes
    if args.collapse_haplotypes:
        df, sample_cols, _ = collapse_haplotypes(df, sample_cols)

    # 3) optionally sanitize labels
    if args.ascii_labels:
        newnames = {c: sanitize_label(c) for c in sample_cols}
        df = df.rename(columns=newnames)
        sample_cols = [newnames[c] for c in sample_cols]

    # 4) enforce sample order (and haplotype order if provided)
    sample_cols = reorder_sample_cols(sample_cols, args.sample_order, args.haplotype_order)
    df = df[["Orthogroup"] + sample_cols + [c for c in df.columns if c in ("count", "category")]]

    # 5) compute counts/categories once
    df = compute_count_and_category(
        df,
        sample_cols,
        softcore_miss=args.softcore_miss,
        cloud_max=args.cloud_max,
        force_recompute_category=args.recompute_category,
    )

    # ---------- Panel A (custom, colored bars, broken y, capped intersections, percentage thresholds) ----------
    panel_A_custom(
        df,
        sample_cols,
        softcore_miss=args.softcore_miss,
        cloud_max=args.cloud_max,
        outbase=os.path.join(args.outdir, "panelA_upset"),
        dpi=args.dpi,
        args=args,
    )

    # ---------- Panel B (pie) ----------
    panel_B_pie(df, os.path.join(args.outdir, "panelB_pie"), dpi=args.dpi)

    # ---------- Panel C (rotated + compressed) ----------
    panel_C_rotated(df, sample_cols, os.path.join(args.outdir, "panelC_heatmap"), dpi=args.dpi)

    # ---------- Panel D (flipped vertically) ----------
    panel_D_stacked_bars_flipped(df, sample_cols, os.path.join(args.outdir, "panelD_per_accession"), dpi=args.dpi)

    # Export the classified matrix (useful downstream)
    out_df = df[["Orthogroup"] + sample_cols + ["count", "category"]]
    out_df.to_csv(os.path.join(args.outdir, "PAV_classified.tsv"), sep="\t", index=False)
    print(f"Done. Figures and tables in: {args.outdir}")


if __name__ == "__main__":
    main()
