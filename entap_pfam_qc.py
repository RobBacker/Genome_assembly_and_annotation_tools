#!/usr/bin/env python3
"""
entap_pfam_qc.py

Parse EnTAP TSV outputs for multiple assemblies and produce:
 - PFAM count matrix (counts per 1k genes; isoforms collapsed)
 - Top-N stacked bar plot (SVG + PDF, dpi=600)  <-- uses colors + hatch patterns
 - Pairwise Pearson correlation (optional)
 - Optional: NLR & secondary-metabolism PFAM-proxy counts
 - AUDIT: write gene→PFAM mappings (sample and/or full)

Counting rule:
  • Aggregate all PFAMs across all rows per gene (after stripping isoform suffix)
  • Count each (gene, PFAM) pair ONCE (i.e., no double counting within a gene)

Author: Robert Backer
"""

from __future__ import annotations
import argparse, glob, os, re, sys
from collections import Counter, defaultdict
import pandas as pd, numpy as np

import matplotlib
matplotlib.use("Agg")
# Make hatches legible in vector outputs
matplotlib.rcParams.update({
    "hatch.linewidth": 0.4,
    "hatch.color": "k",
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import itertools

# ----------------------------
# Config
# ----------------------------
DEFAULT_TOP_PFAM = 30
DPI = 600
GENE_COL = "Query_Sequence"

PFAM_COLS = [
    "Database_UniProt_Protein_Domains",
    "Database_EggNOG_Protein_Domains",
    "Database_InterProScan_Protein_Description",
]

# Optional PFAM proxies (NLR & secondary metabolism)
PFAM_SETS = {
    "NB-ARC (NLR core)": {"PF00931"},
    "TIR": {"PF01582"},
    "RPW8": {"PF05659"},
    "LRR": {"PF00560","PF07723","PF12799","PF13306","PF13855"},
    "CYP450": {"PF00067"},
    "Terpene_synthase": {"PF01397"},
    "UGT (glycosyltransferase)": {"PF00201","PF13419"},
    "BAHD_acyltransferase": {"PF02458"},
    "O-methyltransferase": {"PF00891"},
    "PKS_family": {"PF00109","PF02801","PF00698"}
}

# ----------------------------
# Regex helpers
# ----------------------------
_isoform_re = re.compile(r'\.\d+$')
_pfam_re = re.compile(r'(PF\d{5})')

def strip_isoform(g: str) -> str:
    if g is None:
        return ""
    return _isoform_re.sub('', str(g).strip())

def safe_read_table(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", engine="python", dtype=str)
        if df.shape[1] == 1 and df.columns[0].count(",") > 0:
            df = pd.read_csv(path, sep=",", engine="python", dtype=str)
        return df.fillna('')
    except Exception as e:
        print(f"[WARN] failed to read {path}: {e}")
        return pd.DataFrame()

def detect_files(indir: str) -> list[str]:
    patterns = ["*entap*.tsv","*entap*.txt","*entap*.csv",
                "*EnTAP*.tsv","*EnTAP*.txt"]
    files = []
    for p in patterns:
        files.extend(glob.glob(os.path.join(indir, p)))
    return sorted(set(files))

def extract_pfams_from_cell(cell: str) -> list[str]:
    if not cell:
        return []
    return list(dict.fromkeys(_pfam_re.findall(str(cell))))

def collect_gene_to_pfams(df: pd.DataFrame,
                          pfam_cols: list[str],
                          gene_col_norm: str = "_gene_norm") -> dict[str, set]:
    """
    Aggregate PFAMs across *all rows* per gene; return {gene: set(PFxxxxx)}.
    """
    gene_to_pf = {}
    cols_present = [c for c in pfam_cols if c in df.columns]
    for g, gdf in df.groupby(gene_col_norm, sort=False):
        s = set()
        for col in cols_present:
            for cell in gdf[col].astype(str).tolist():
                s.update(extract_pfams_from_cell(cell))
        gene_to_pf[g] = s
    return gene_to_pf

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(description="QC EnTAP PFAM annotations across assemblies")
    ap.add_argument('--indir', default='.', help='Directory containing EnTAP files')
    ap.add_argument('--outdir', default='entap_pfam_results', help='Output directory')
    ap.add_argument('--top_pfam', type=int, default=DEFAULT_TOP_PFAM, help='Top N Pfam domains to plot')
    ap.add_argument('--no_corr', action='store_true', help='Skip pairwise Pearson correlations')
    ap.add_argument('--pfam_proxies', action='store_true', help='Also report PFAM proxy counts per assembly')
    ap.add_argument('--audit_full', action='store_true', help='Write FULL gene→PFAM mapping (can be large)')
    ap.add_argument('--audit_sample_n', type=int, default=500, help='Write SAMPLE (first N) gene→PFAM rows per assembly')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    files = detect_files(args.indir)
    if not files:
        sys.exit("[ERROR] No EnTAP files found in --indir.")

    print(f"[INFO] Found {len(files)} file(s). Parsing…")

    # Per-assembly storage
    asm_gene_to_pf: dict[str, dict[str, set]] = {}
    gene_counts: dict[str, int] = {}
    annotated_gene_counts: dict[str, int] = {}

    for fpath in files:
        fname = os.path.basename(fpath)
        asm = re.sub(r'\.(entap|EnTAP|annotation|annot).*', '', fname, flags=re.I)
        asm = re.sub(r'\.(tsv|txt|csv)$', '', asm, flags=re.I).strip()

        df = safe_read_table(fpath)
        if df.empty:
            print(f"[WARN] {fname} empty/unreadable; skipping.")
            continue

        if GENE_COL not in df.columns:
            print(f"[WARN] {fname} missing '{GENE_COL}'; skipping.")
            continue

        # Normalize gene IDs
        df[GENE_COL] = df[GENE_COL].astype(str).str.strip()
        df = df[df[GENE_COL] != '']
        if df.empty:
            print(f"[WARN] {fname} has no gene IDs after trimming; skipping.")
            continue
        df["_gene_norm"] = df[GENE_COL].apply(strip_isoform)

        # Unique genes in this assembly
        uniq_genes = df["_gene_norm"].unique().tolist()
        gene_counts[asm] = len(uniq_genes)

        # Build mapping gene to set(PFxxxxx)
        gene_to_pf = collect_gene_to_pfams(df, PFAM_COLS, gene_col_norm="_gene_norm")
        annotated_gene_counts[asm] = sum(1 for pfs in gene_to_pf.values() if pfs)
        asm_gene_to_pf[asm] = gene_to_pf

        print(f"[INFO] {asm}: unique genes={gene_counts[asm]}, with≥1 PFAM={annotated_gene_counts[asm]}")

        # ---- AUDIT: sample gene to PFAM mapping ----
        sample_n = max(0, int(args.audit_sample_n))
        if sample_n > 0:
            audit_path = os.path.join(args.outdir, f"{asm}_gene_to_pfam_SAMPLE.tsv")
            with open(audit_path, "w") as fh:
                print("Gene\tPFAMs", file=fh)
                for i, (g, s) in enumerate(gene_to_pf.items()):
                    if i >= sample_n:
                        break
                    print(f"{g}\t{';'.join(sorted(s))}", file=fh)
            print(f"[INFO] Wrote audit sample: {audit_path}")

        # ---- AUDIT: full mapping (optional) ----
        if args.audit_full:
            full_path = os.path.join(args.outdir, f"{asm}_gene_to_pfam_FULL.tsv")
            with open(full_path, "w") as fh:
                print("Gene\tPFAMs", file=fh)
                for g, s in gene_to_pf.items():
                    print(f"{g}\t{';'.join(sorted(s))}", file=fh)
            print(f"[INFO] Wrote FULL mapping: {full_path}")

    if not asm_gene_to_pf:
        sys.exit("[ERROR] No assemblies processed successfully.")

    # ----------------------------
    # Build PFAM count matrix (per 1k genes)
    # ----------------------------
    # Collect global PFAM universe
    all_pfams = set()
    for gene_to_pf in asm_gene_to_pf.values():
        for s in gene_to_pf.values():
            all_pfams.update(s)
    all_pfams = sorted(all_pfams)

    # (assembly x pfam) counts of UNIQUE (gene, PFAM) pairs
    mat_counts = pd.DataFrame(0, index=sorted(asm_gene_to_pf.keys()), columns=all_pfams, dtype=float)
    for asm, gene_to_pf in asm_gene_to_pf.items():
        c = Counter()
        for g, s in gene_to_pf.items():
            for pf in s:
                c[pf] += 1
        for pf, v in c.items():
            mat_counts.loc[asm, pf] = v

    # Normalize to per-1k genes
    mat_per1k = mat_counts.copy()
    for asm in mat_per1k.index:
        denom = float(max(gene_counts.get(asm, 0), 1))
        mat_per1k.loc[asm] = (mat_per1k.loc[asm] / denom) * 1000.0

    # Keep top N by mean abundance
    if DEFAULT_TOP_PFAM and mat_per1k.shape[1] > DEFAULT_TOP_PFAM:
        pass  # default is overridden by CLI below
    top_n = int(args.top_pfam) if args.top_pfam else DEFAULT_TOP_PFAM
    if mat_per1k.shape[1] > top_n:
        keep = mat_per1k.mean(axis=0).sort_values(ascending=False).index[:top_n]
        mat_top = mat_per1k.loc[:, keep]
    else:
        mat_top = mat_per1k

    outdir = args.outdir
    mat_per1k.to_csv(os.path.join(outdir, "PFAM_counts_per1k_all.tsv"), sep="\t")
    mat_top.to_csv(os.path.join(outdir, "PFAM_counts_per1k_TOP.tsv"), sep="\t")

    # ----------------------------
    # Pairwise correlation (optional)
    # ----------------------------
    if not args.no_corr and mat_top.shape[1] >= 2:
        corr = mat_top.T.corr(method="pearson")
        corr.to_csv(os.path.join(outdir, "PFAM_pairwise_pearson.tsv"), sep="\t")
        # Summary
        tril = corr.where(np.tril(np.ones(corr.shape), k=-1).astype(bool))
        vals = tril.values.flatten()
        vals = vals[~np.isnan(vals)]
        if len(vals):
            with open(os.path.join(outdir, "PFAM_corr_summary.txt"), "w") as fh:
                fh.write("PFAM Pearson r summary\n")
                fh.write(f"min\t{vals.min():.6f}\n")
                fh.write(f"median\t{np.median(vals):.6f}\n")
                fh.write(f"max\t{vals.max():.6f}\n")

    # ----------------------------
    # Plot stacked bar of Top-N PFAMs (colors + hatch patterns)
    # ----------------------------
    if not mat_top.empty:
        # color cycle (20 colorblind-safe colors), repeat as needed
        base_colors = list(plt.cm.tab20.colors)
        color_cycle = list(itertools.islice(itertools.cycle(base_colors), mat_top.shape[1]))

        # hatch cycle to augment colors when >20 series
        hatches = ["", "///", "\\\\\\", "xxx", "++", "...", "oo", "**", "//.", "\\\\.", "||", "--"]
        hatch_cycle = list(itertools.islice(itertools.cycle(hatches), mat_top.shape[1]))

        fig, ax = plt.subplots(figsize=(12, 6))
        bottom = np.zeros(mat_top.shape[0], dtype=float)

        # draw each PFAM layer with a unique (color, hatch)
        bar_containers = []
        for i, col in enumerate(mat_top.columns):
            vals = mat_top[col].values
            bars = ax.bar(
                mat_top.index, vals,
                bottom=bottom,
                color=color_cycle[i],
                edgecolor="black", linewidth=0.35,
                label=col
            )
            for b in bars:
                b.set_hatch(hatch_cycle[i])
            bottom += vals
            bar_containers.append(bars)

        ax.set_ylabel("Counts per 1k genes")
        ax.set_title(f"Top {mat_top.shape[1]} PFAM domains (per 1k genes)",
                     loc="left", fontsize=12, fontweight="bold")
        ax.tick_params(axis="x", rotation=45)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # legend showing both color + hatch
        legend_handles = [
            Patch(facecolor=color_cycle[i], edgecolor="black",
                  linewidth=0.35, hatch=hatch_cycle[i], label=col)
            for i, col in enumerate(mat_top.columns)
        ]
        ax.legend(
            handles=legend_handles,
            frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left",
            fontsize=7, title="PFAM", ncol=1
        )

        plt.tight_layout()
        for ext in ("svg", "pdf"):
            fig.savefig(os.path.join(outdir, f"PFAM_top_stacked.{ext}"),
                        dpi=DPI, bbox_inches="tight")
        plt.close(fig)

    # ----------------------------
    # Optional: PFAM proxy counts per assembly
    # ----------------------------
    if args.pfam_proxies:
        rows = []
        for asm, gene_to_pf in asm_gene_to_pf.items():
            total = gene_counts.get(asm, 0)
            rec = {"Assembly": asm, "TotalGenes": total}
            for name, pfset in PFAM_SETS.items():
                n = sum(1 for s in gene_to_pf.values() if s.intersection(pfset))
                rec[name] = n
            rows.append(rec)
        df_proxy = pd.DataFrame(rows).set_index("Assembly").fillna(0).astype(int)
        df_proxy.to_csv(os.path.join(outdir, "PFAM_proxies_counts.tsv"), sep="\t")
        # per-1k normalization
        norm = df_proxy.div(df_proxy["TotalGenes"], axis=0) * 1000.0
        if "TotalGenes" in norm.columns:
            norm = norm.drop(columns=["TotalGenes"])
        norm.to_csv(os.path.join(outdir, "PFAM_proxies_per1k.tsv"), sep="\t")

    # ----------------------------
    # Summary file
    # ----------------------------
    with open(os.path.join(outdir, "PFAM_qc_summary.txt"), "w") as fh:
        fh.write("Assembly\tUniqueGenes\tGenes_with_PFAM\tPct_with_PFAM\n")
        for asm in sorted(gene_counts.keys()):
            n = gene_counts[asm]
            a = annotated_gene_counts.get(asm, 0)
            pct = (a / n * 100.0) if n else 0.0
            fh.write(f"{asm}\t{n}\t{a}\t{pct:.1f}\n")

    print("[INFO] Done.")

if __name__ == "__main__":
    main()
