#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as mlines

# ─────────── 0) PARSE ARGUMENTS ───────────
parser = argparse.ArgumentParser(description="Categorize AED from GFF3 and mark PSAURON status")
parser.add_argument("-g", "--gff",      required=True, help="Input GFF3 with AED fields (aed_ev_tr/aed_ev_pr)")
parser.add_argument("-p", "--psauron", required=True, help="PSAURON CSV (skip 3 header rows)")
parser.add_argument("-w", "--workdir",  default=".",   help="Working directory for outputs")
args = parser.parse_args()

os.makedirs(args.workdir, exist_ok=True)
os.chdir(args.workdir)

# ─────────── 1) AED THRESHOLDS ───────────
TR_HIGH, PR_HIGH     = 0.10, 0.05
TR_GOOD, PR_GOOD     = 0.25, 0.10
TR_MOD,  PR_MOD      = 0.50, 0.20
AB_THRESHOLD         = 0.90  # AED threshold for ab initio (cat6)

# ─────────── 2) LOAD PSAURON STATUS ───────────
def load_psauron(fn):
    dfp = pd.read_csv(fn, skiprows=3)
    return dict(zip(dfp['description'], dfp['psauron_is_protein'].astype(bool)))

ps_map = load_psauron(args.psauron)

# ─────────── 3) READ GFF3 AND EXTRACT AED ───────────
records = []
with open(args.gff) as f:
    for line in f:
        if line.startswith("#"): continue
        cols = line.strip().split("\t")
        if cols[2] != "mRNA": continue
        attrs = dict(x.split("=",1) for x in cols[8].split(";") if "=" in x)
        if "aed_ev_tr" in attrs and "aed_ev_pr" in attrs:
            try:
                tr, pr = float(attrs["aed_ev_tr"]), float(attrs["aed_ev_pr"])
            except ValueError:
                continue
            tid = attrs.get("ID") or attrs.get("transcript_id")
            records.append({"transcript_id": tid, "aed_ev_tr": tr, "aed_ev_pr": pr})

df = pd.DataFrame(records)
if df.empty:
    sys.exit("❌ No mRNA with aed_ev_tr/aed_ev_pr found in your GFF3. Check attribute names!")

# ─────────── 4) CATEGORIZE ───────────
def classify(tr, pr):
    # Ab initio if both AEDs exceed threshold
    if tr >= AB_THRESHOLD and pr >= AB_THRESHOLD:
        return "cat6"
    if tr <= TR_HIGH and pr <= PR_HIGH:
        return "cat1"
    if tr <= TR_GOOD and pr <= PR_GOOD:
        return "cat2"
    if tr <= TR_MOD and pr <= PR_MOD:
        return "cat3"
    if (tr <= TR_MOD and pr > PR_MOD) or (pr <= PR_MOD and tr > TR_MOD):
        return "cat4"
    return "cat5"

df["curation"] = df.apply(lambda r: classify(r.aed_ev_tr, r.aed_ev_pr), axis=1)
df["psauron"] = df["transcript_id"].map(lambda tid: ps_map.get(tid, False))

# Compute category counts for legend
cat_counts = df['curation'].value_counts().to_dict()
ps_false_count = (~df['psauron']).sum()

# ─────────── 5) PLOT ───────────
cat_order = ["cat1","cat2","cat3","cat4","cat5","cat6"]
cat_labels = {
    "cat1": f"cat1: High (tr ≤ {TR_HIGH}, pr ≤ {PR_HIGH})",
    "cat2": f"cat2: Good (tr ≤ {TR_GOOD}, pr ≤ {PR_GOOD})",
    "cat3": f"cat3: Moderate (tr ≤ {TR_MOD}, pr ≤ {PR_MOD})",
    "cat4": f"cat4: One support only (tr ≤ {TR_MOD} or pr ≤ {PR_MOD})",
    "cat5": "cat5: Weak support",
    "cat6": f"cat6: Ab initio (AED ≥ {AB_THRESHOLD} for both)"
}
custom_palette = {
    "cat1":"#9467bd",   # purple
    "cat2":"#0492c2",   # blue
    "cat3":"#3bb143",   # green
    "cat4":"#ff9913",   # orange
    "cat5":"#e3242b",   # red
    "cat6":"#00CED1"    # DarkTurquoise
}

# prepare figure
grid = sns.JointGrid(data=df, x="aed_ev_tr", y="aed_ev_pr", height=10, space=0.3)

# colored dots for PSAURON=True
df_true = df[df["psauron"]]
sns.scatterplot(data=df_true, x="aed_ev_tr", y="aed_ev_pr",
                hue="curation", hue_order=cat_order, palette=custom_palette,
                s=5, alpha=0.7, edgecolor="none", linewidth=0,
                ax=grid.ax_joint, legend=False)

# black dots for PSAURON=False
df_false = df[~df["psauron"]]
sns.scatterplot(data=df_false, x="aed_ev_tr", y="aed_ev_pr",
                color="black", s=5, alpha=0.7,
                edgecolor="none", linewidth=0,
                ax=grid.ax_joint, legend=False)

# axis labels
grid.ax_joint.set_xlabel("AED with transcript evidence", fontsize=12)
grid.ax_joint.set_ylabel("AED with protein evidence", fontsize=12)
grid.ax_marg_x.set_ylabel("# transcripts (transcript evidence)", fontsize=10)
grid.ax_marg_y.set_xlabel("# transcripts (protein evidence)", fontsize=10)

# threshold lines
grid.ax_joint.axvline(x=TR_MOD, linestyle="--", color="black", linewidth=1)
grid.ax_joint.axhline(y=PR_MOD, linestyle="--", color="black", linewidth=1)

# quadrant counts and annotations
quad_counts = {
    "top left":    len(df[(df.aed_ev_tr<=TR_MOD)&(df.aed_ev_pr>PR_MOD)]),
    "top right":   len(df[(df.aed_ev_tr>TR_MOD)&(df.aed_ev_pr>PR_MOD)]),
    "bottom left": len(df[(df.aed_ev_tr<=TR_MOD)&(df.aed_ev_pr<=PR_MOD)]),
    "bottom right":len(df[(df.aed_ev_tr>TR_MOD)&(df.aed_ev_pr<=PR_MOD)])
}
positions = {"top left":(TR_MOD/2,0.9),"top right":(0.75,0.9),"bottom left":(TR_MOD/2,PR_MOD/2),"bottom right":(0.75,PR_MOD/2)}
for label,(x,y) in positions.items():
    grid.ax_joint.text(x,y,f"{quad_counts[label]}",fontsize=11,weight='bold',color="black",
                       ha="center",va="center",
                       bbox=dict(boxstyle="round,pad=0.3",facecolor="white",
                                 edgecolor="black",linewidth=0.5,alpha=0.65))

# marginal histograms without legends
sns.histplot(data=df_true, x="aed_ev_tr", hue="curation", hue_order=cat_order,
             palette=custom_palette, bins=40, ax=grid.ax_marg_x,
             element="step", stat="count", common_norm=False, legend=False)
sns.histplot(data=df_false, x="aed_ev_tr", color="black", bins=40,
             ax=grid.ax_marg_x, element="step", stat="count", common_norm=False,
             alpha=0.3, legend=False)
sns.histplot(data=df_true, y="aed_ev_pr", hue="curation", hue_order=cat_order,
             palette=custom_palette, bins=40, ax=grid.ax_marg_y,
             element="step", stat="count", common_norm=False, legend=False)
sns.histplot(data=df_false, y="aed_ev_pr", color="black", bins=40,
             ax=grid.ax_marg_y, element="step", stat="count", common_norm=False,
             alpha=0.3, legend=False)

# custom legend with counts
handles = [mlines.Line2D([],[],color="black",marker='o',linestyle='None',markersize=6,
                          label=f"psauron=FALSE (n={ps_false_count})")]
for cat in cat_order:
    count = cat_counts.get(cat,0)
    label = f"{cat_labels[cat]} (n={count})"
    handles.append(mlines.Line2D([],[],color=custom_palette[cat],marker='o',linestyle='None',markersize=6,
                                  label=label))

grid.ax_joint.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.27,-0.4), ncol=1, frameon=False)

# title & save
grid.fig.suptitle("AED (transcript vs protein) with custom curation categories", y=1.02)
grid.fig.tight_layout()
grid.fig.savefig("aed_curation_categories.png", dpi=600, bbox_inches="tight")
plt.close(grid.fig)

# ─────────── 6) OUTPUT TSV, KEEP & REMOVE LISTS ───────────
df[["transcript_id","aed_ev_tr","aed_ev_pr","curation","psauron"]].to_csv(
    "transcripts_aed_categorized.tsv", sep="\t", index=False)

keep = df[df.curation.isin(["cat1","cat2","cat3","cat4"]) | 
          ((df.curation.isin(["cat5","cat6"])) & df.psauron)]["transcript_id"].drop_duplicates()
keep.to_csv("curation_recommended_keep_transcripts.txt", index=False, header=False)

remove = df[df.curation.isin(["cat5","cat6"]) & (~df.psauron)]["transcript_id"].drop_duplicates()
remove.to_csv("curation_remove_transcripts.txt", index=False, header=False)

# ─── 7) ANNOTATE GFF3 ───
out_gff = os.path.basename(args.gff).replace(".gff3", ".uncurated.gff3")
header = f"""## AED curation categories (curation attribute)
## cat1: High confidence — transcript AED ≤ {TR_HIGH}, protein AED ≤ {PR_HIGH}
## cat2: Good support — transcript AED ≤ {TR_GOOD}, protein AED ≤ {PR_GOOD}
## cat3: Moderate support — transcript AED ≤ {TR_MOD}, protein AED ≤ {PR_MOD}
## cat4: One evidence type only — either transcript ≤ {TR_MOD} or protein ≤ {PR_MOD}
## cat5: Weak support — transcript & protein AED < {AB_THRESHOLD}
## cat6: Ab‑initio only — transcript & protein AED ≥ {AB_THRESHOLD}
## psauron=TRUE/FALSE indicates protein‑model status from PSAURON
"""
with open(args.gff) as fin, open(out_gff,"w") as fout:
    fout.write(header)
    for line in fin:
        if line.startswith("#"): fout.write(line); continue
        cols = line.rstrip("\n").split("\t")
        if cols[2] == "mRNA":
            attrs = dict(x.split("=",1) for x in cols[8].split(";") if "=" in x)
            tid = attrs.get("ID") or attrs.get("transcript_id")
            if tid in df.transcript_id.values:
                attrs["curation"] = df.loc[df.transcript_id==tid,"curation"].iat[0]
            attrs["psauron"] = str(ps_map.get(tid,False)).upper()
            cols[8] = ";".join(f"{k}={v}" for k,v in attrs.items())
        fout.write("\t".join(cols)+"\n")

print("✅ Written: aed_curation_categories.png, transcripts_aed_categorized.tsv, curation_recommended_keep_transcripts.txt, curation_remove_transcripts.txt, {out_gff}")

# ─────────── 8) WRITE README ───────────
readme_content = f"""
AED Curation Pipeline Results
=============================

**Overview**
This folder contains outputs and visualizations from the AED curation script. The model classifications and PSAURON support tags allow you to:
- Quickly assess model quality distribution
- Identify high‑confidence, moderately supported, and purely ab initio gene models
- Generate keep/remove lists for downstream filtering

**Threshold Definitions**
| Category | Description                                               |
|----------|-----------------------------------------------------------|
| cat1     | High confidence: tr ≤ {TR_HIGH}, pr ≤ {PR_HIGH}           |
| cat2     | Good support:    tr ≤ {TR_GOOD}, pr ≤ {PR_GOOD}           |
| cat3     | Moderate:        tr ≤ {TR_MOD}, pr ≤ {PR_MOD}             |
| cat4     | One evidence only: transcript ≤ {TR_MOD} or protein ≤ {PR_MOD} |
| cat5     | Weak support: transcripts/proteins above moderate, below AB threshold ({AB_THRESHOLD}) |
| cat6     | Ab initio:       tr ≥ {AB_THRESHOLD} and pr ≥ {AB_THRESHOLD} |

**Files and Descriptions**
1. **aed_curation_categories.png**
   - Joint scatter plot of `aed_ev_tr` vs. `aed_ev_pr` with marginal histograms.
   - Points colored by category; legend includes counts and PSAURON‑FALSE entries.

2. **transcripts_aed_categorized.tsv**
   - Columns: `transcript_id`, `aed_ev_tr`, `aed_ev_pr`, `curation`, `psauron`.
   - One row per mRNA model with classification and PSAURON support (TRUE/FALSE).

3. **curation_recommended_keep_transcripts.txt**
   - List of `transcript_id`s for categories cat1–4, plus any cat5/cat6 with PSAURON support.
   - Useful for extracting high‑confidence models.

4. **curation_remove_transcripts.txt**
   - List of `transcript_id`s for cat5/cat6 without PSAURON support.
   - Candidates for removal or manual review.

5. **{out_gff}**
   - Annotated GFF3 with two new attributes on each `mRNA` line:
     - `curation=<catX>`
     - `psauron=TRUE|FALSE`

**Usage**
```bash
python3 aed_curation.py \
  -g path/to/input.gff3 \
  -p path/to/psauron_output.csv \
  -w path/to/results_dir
```

- Ensure `psauron_output.csv` skips the first 3 header rows.
- Results directory will be created if it does not exist.

**Next Steps**
- Use the TSV or keep/remove lists to filter FASTA or GFF3 for high‑confidence models.
- Adjust AED threshold variables at the top of the script for custom stringency.

"""
readme_path = os.path.join(args.workdir, "README_AED_Curation_Results.txt")
with open(readme_path, "w") as fh:
    fh.write(readme_content)
print(f"📄 Detailed README written to {readme_path}")
