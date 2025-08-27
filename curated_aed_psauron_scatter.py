#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as mlines

# ─────────── 0) PARSE ARGUMENTS ───────────
parser = argparse.ArgumentParser(description="Regenerate AED plot from already curated GFF3")
parser.add_argument("-g", "--gff", required=True, help="Curated GFF3 with 'curation', 'aed_ev_tr', 'aed_ev_pr', and 'psauron' tags")
parser.add_argument("-o", "--output", default="aed_curation_categories.png", help="Output PNG filename")
args = parser.parse_args()

# ─────────── 1) THRESHOLDS & LABELS ───────────
TR_HIGH, PR_HIGH = 0.10, 0.05
TR_GOOD, PR_GOOD = 0.25, 0.10
TR_MOD,  PR_MOD  = 0.50, 0.20
AB_THRESHOLD     = 0.90

cat_order = ["cat1", "cat2", "cat3", "cat4", "cat5", "cat6"]
cat_labels = {
    "cat1": f"cat1: High (tr ≤ {TR_HIGH}, pr ≤ {PR_HIGH})",
    "cat2": f"cat2: Good (tr ≤ {TR_GOOD}, pr ≤ {PR_GOOD})",
    "cat3": f"cat3: Moderate (tr ≤ {TR_MOD}, pr ≤ {PR_MOD})",
    "cat4": f"cat4: One support only (tr ≤ {TR_MOD} or pr ≤ {PR_MOD})",
    "cat5": "cat5: Weak support",
    "cat6": f"cat6: Ab initio (AED ≥ {AB_THRESHOLD} for both)"
}
custom_palette = {
    "cat1": "#9467bd", "cat2": "#0492c2", "cat3": "#3bb143",
    "cat4": "#ff9913", "cat5": "#e3242b", "cat6": "#00CED1"
}

# ─────────── 2) LOAD CURATED GFF3 ───────────
records = []
with open(args.gff) as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.strip().split("\t")
        if cols[2] != "mRNA":
            continue
        attrs = dict(x.split("=", 1) for x in cols[8].split(";") if "=" in x)
        if not {"aed_ev_tr", "aed_ev_pr", "curation", "psauron"}.issubset(attrs):
            continue
        try:
            tr = float(attrs["aed_ev_tr"])
            pr = float(attrs["aed_ev_pr"])
            psauron = attrs["psauron"].upper() == "TRUE"
        except ValueError:
            continue
        tid = attrs.get("ID") or attrs.get("transcript_id")
        records.append({
            "transcript_id": tid,
            "aed_ev_tr": tr,
            "aed_ev_pr": pr,
            "curation": attrs["curation"],
            "psauron": psauron
        })

df = pd.DataFrame(records)
if df.empty:
    sys.exit("❌ No valid mRNA entries with all required tags found in GFF3.")

# ─────────── 3) PLOT ───────────
cat_counts = df["curation"].value_counts().to_dict()
ps_false_count = (~df["psauron"]).sum()

grid = sns.JointGrid(data=df, x="aed_ev_tr", y="aed_ev_pr", height=10, space=0.3)

# Colored points: psauron=TRUE
df_true = df[df["psauron"]]
sns.scatterplot(data=df_true, x="aed_ev_tr", y="aed_ev_pr",
                hue="curation", hue_order=cat_order, palette=custom_palette,
                s=5, alpha=0.7, edgecolor="none", linewidth=0,
                ax=grid.ax_joint, legend=False)

# Black points: psauron=FALSE
df_false = df[~df["psauron"]]
sns.scatterplot(data=df_false, x="aed_ev_tr", y="aed_ev_pr",
                color="black", s=5, alpha=0.7, edgecolor="none", linewidth=0,
                ax=grid.ax_joint, legend=False)

# Axis labels & thresholds
grid.ax_joint.set_xlabel("AED with transcript evidence", fontsize=12)
grid.ax_joint.set_ylabel("AED with protein evidence", fontsize=12)
grid.ax_marg_x.set_ylabel("# transcripts (transcript evidence)", fontsize=10)
grid.ax_marg_y.set_xlabel("# transcripts (protein evidence)", fontsize=10)
grid.ax_joint.axvline(x=TR_MOD, linestyle="--", color="black", linewidth=1)
grid.ax_joint.axhline(y=PR_MOD, linestyle="--", color="black", linewidth=1)

# Quadrant annotations
quad_counts = {
    "top left":    len(df[(df.aed_ev_tr <= TR_MOD) & (df.aed_ev_pr > PR_MOD)]),
    "top right":   len(df[(df.aed_ev_tr > TR_MOD) & (df.aed_ev_pr > PR_MOD)]),
    "bottom left": len(df[(df.aed_ev_tr <= TR_MOD) & (df.aed_ev_pr <= PR_MOD)]),
    "bottom right":len(df[(df.aed_ev_tr > TR_MOD) & (df.aed_ev_pr <= PR_MOD)])
}
positions = {"top left": (TR_MOD/2, 0.9), "top right": (0.75, 0.9),
             "bottom left": (TR_MOD/2, PR_MOD/2), "bottom right": (0.75, PR_MOD/2)}
for label, (x, y) in positions.items():
    grid.ax_joint.text(x, y, f"{quad_counts[label]}", fontsize=11, weight='bold', color="black",
                       ha="center", va="center",
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                 edgecolor="black", linewidth=0.5, alpha=0.65))

# Marginal histograms
sns.histplot(data=df_true, x="aed_ev_tr", hue="curation", hue_order=cat_order,
             palette=custom_palette, bins=40, ax=grid.ax_marg_x,
             element="step", stat="count", common_norm=False, legend=False)
sns.histplot(data=df_false, x="aed_ev_tr", color="black", bins=40,
             ax=grid.ax_marg_x, element="step", stat="count", alpha=0.3)
sns.histplot(data=df_true, y="aed_ev_pr", hue="curation", hue_order=cat_order,
             palette=custom_palette, bins=40, ax=grid.ax_marg_y,
             element="step", stat="count", common_norm=False, legend=False)
sns.histplot(data=df_false, y="aed_ev_pr", color="black", bins=40,
             ax=grid.ax_marg_y, element="step", stat="count", alpha=0.3)

# Custom legend with counts
handles = [mlines.Line2D([], [], color="black", marker='o', linestyle='None', markersize=6,
                         label=f"psauron=FALSE (n={ps_false_count})")]
for cat in cat_order:
    count = cat_counts.get(cat, 0)
    handles.append(mlines.Line2D([], [], color=custom_palette.get(cat, "gray"), marker='o', linestyle='None',
                                 markersize=6, label=f"{cat_labels.get(cat)} (n={count})"))
grid.ax_joint.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.27, -0.4),
                     ncol=1, frameon=False)

# Title and save
grid.fig.suptitle("AED (transcript vs protein) after auto-curation", y=1.02)
grid.fig.tight_layout()
grid.fig.savefig(args.output, dpi=600, bbox_inches="tight")
plt.close(grid.fig)

print(f"✅ AED plot regenerated and saved to: {args.output}")
