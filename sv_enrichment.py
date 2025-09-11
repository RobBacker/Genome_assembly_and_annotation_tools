#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NLR SV & SNV enrichment on pangene windows (WAVE-derived BED8s).

Inputs (produced by your Snippet 7 bash):
- HAL_OUT/sv/bed/wave.sv50.bed8              (SV BED8; type in col7, len in col8)
- HAL_OUT/pangene_loci.ref.bed6              (pangene locus spans on reference)
- HAL_OUT/sv/chrom.sizes                     (chrom\tlength)
Optional:
- HAL_OUT/sv/bed/wave.snp_mnp.bed8           (SNV BED8: SNP+MNP)

Also needs:
- HAL_OUT/gene_to_locus.tsv  (preferred) or HAL_OUT/pangene_members.tsv (fallback)
- A gene list file (one gene_id per line) to define “NLR windows”.

Outputs in --outdir:
- observed_counts.tsv            (per-window counts; includes per-type, incl. n_other)
- sv_density_summary.tsv         (global density summary; includes SNV if present)
- bootstrap_null.tsv             (SV null densities + per-type)
- bootstrap_null_snv.tsv         (SNV null densities + per-type, if SNV provided)
- bootstrap_summary.txt          (obs density, N, p-value, settings)
- figure_enrichment.svg/.pdf     (THREE panels: SV, SNV, per-type enrichment)
- type_enrichment.tsv            (Panel C data table)

Author: Robert Backer
"""

import argparse
import os
import sys
import tempfile
import random
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
from collections import defaultdict, Counter
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------- Helpers: IO, paths ----------

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)
    return p

def check_tool(cmd="bedtools"):
    try:
        subprocess.run([cmd, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
    except Exception as e:
        sys.exit(f"ERROR: required tool '{cmd}' not found in PATH: {e}")

def read_chrom_sizes(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom","length"], dtype={"chrom":str, "length":np.int64})
    return df

def read_bed6(path):
    cols = ["chrom","start","end","name","score","strand"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols,
                     dtype={"chrom":str, "start":np.int64, "end":np.int64, "name":str, "score":str, "strand":str})
    return df

def read_bed8(path):
    cols = ["chrom","start","end","name","score","strand","type","len"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols,
                     dtype={"chrom":str, "start":np.int64, "end":np.int64, "name":str, "score":str,
                            "strand":str, "type":str, "len":np.float64})
    df["type"] = df["type"].fillna("OTHER").str.upper()
    df["len"] = df["len"].fillna(0).astype(float)
    return df

def write_bed6(df, path):
    df = df.copy()
    df["score"] = df.get("score", 0)
    df["strand"] = df.get("strand", ".")
    df[["chrom","start","end","name","score","strand"]].to_csv(path, sep="\t", header=False, index=False)

def write_bed(df, path):
    df.to_csv(path, sep="\t", header=False, index=False)

def run_bedtools_intersect_count(a_bed, b_bed):
    """bedtools intersect -c -a A -b B -> rows (A fields + count)"""
    cmd = ["bedtools", "intersect", "-c", "-a", a_bed, "-b", b_bed]
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    lines = out.stdout.strip().splitlines()
    rows = []
    for line in lines:
        parts = line.split("\t")
        count = int(parts[-1])
        rows.append(parts[:-1] + [count])
    return rows

def run_bedtools_intersect_wb(a_bed, b_bed):
    """bedtools intersect -wa -wb -a A -b B -> list of [A..., B...] parts"""
    cmd = ["bedtools", "intersect", "-wa", "-wb", "-a", a_bed, "-b", b_bed]
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    return [l.split("\t") for l in out.stdout.strip().splitlines() if l]

# ---------- Type color palette ----------

def load_type_colors(path_or_none):
    """
    Load a type->color mapping from TSV/CSV/JSON file with columns/keys: type,color.
    If None or file missing/invalid, return a default palette.
    """
    default = {
        "SNP":      "#4C78A8",
        "MNP":      "#9ECAE1",
        "INS":      "#F58518",
        "DEL":      "#E45756",
        "INV":      "#72B7B2",
        "COMPLEX":  "#B279A2",
        "OTHER":    "#54A24B",
    }
    if not path_or_none or not os.path.exists(path_or_none):
        return default
    try:
        if path_or_none.lower().endswith(".json"):
            import json
            with open(path_or_none) as fh:
                data = json.load(fh)
            # normalize keys to upper
            return {str(k).upper(): str(v) for k, v in data.items()}
        else:
            df = pd.read_csv(path_or_none, sep=None, engine="python")
            if "type" in df.columns and "color" in df.columns:
                return {str(r["type"]).upper(): str(r["color"]) for _, r in df.iterrows()}
    except Exception as e:
        print(f"WARN: failed to load type colors from {path_or_none}: {e}", file=sys.stderr)
    return default

# ---------- Mapping: gene -> pangene_id ----------

def load_gene_to_pangene(hal_out):
    g2l = os.path.join(hal_out, "gene_to_locus.tsv")
    mem = os.path.join(hal_out, "pangene_members.tsv")
    if os.path.exists(g2l):
        df = pd.read_csv(g2l, sep="\t")
        if "gene_id" not in df or "pangene_id" not in df:
            raise RuntimeError(f"{g2l} missing required columns")
        mapping = df[["gene_id","pangene_id"]].dropna().drop_duplicates()
    elif os.path.exists(mem):
        df = pd.read_csv(mem, sep="\t", header=None, names=["pangene_id","region","gene_id","asm","ctx"])
        mapping = df[["gene_id","pangene_id"]].dropna().drop_duplicates()
    else:
        raise FileNotFoundError("Neither gene_to_locus.tsv nor pangene_members.tsv found in HAL_OUT")
    mapping = mapping.drop_duplicates(subset=["gene_id"], keep="first")
    return dict(zip(mapping["gene_id"], mapping["pangene_id"]))

def load_gene_list(path):
    with open(path) as fh:
        genes = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    return genes


# ---------- Build NLR windows from gene list + pangene loci ----------

def build_nlr_windows_from_genes(gene_list, gene2pangene, loci_bed6, chrom_sizes_df, flank_kb):
    loci_df = read_bed6(loci_bed6)
    loci_map = {r["name"]: (r["chrom"], int(r["start"]), int(r["end"]), r.get("strand",".")) for _, r in loci_df.iterrows()}

    chrom_max = dict(zip(chrom_sizes_df["chrom"], chrom_sizes_df["length"]))
    flank_bp = int(flank_kb * 1000)

    rows = []
    skipped = []
    seen_pids = set()

    for gid in gene_list:
        pid = gene2pangene.get(gid)
        if pid is None:
            skipped.append(gid)
            continue
        if pid in seen_pids:
            continue
        seen_pids.add(pid)
        if pid not in loci_map:
            skipped.append(gid)
            continue
        chr_, s, e, st = loci_map[pid]
        s0 = max(0, s - flank_bp)
        e0 = e + flank_bp
        L = chrom_max.get(chr_)
        if L is not None:
            e0 = min(e0, int(L))
        rows.append([chr_, s0, e0, pid, 0, "."])

    nlr_df = pd.DataFrame(rows, columns=["chrom","start","end","name","score","strand"])
    if len(skipped) > 0:
        print(f"NOTE: {len(skipped)} gene IDs had no pangene/locus match; ignoring (e.g. {skipped[:5]})", file=sys.stderr)
    return nlr_df


# ---------- Variant filtering ----------

def filter_bed8_by_args(bed8_df, types=None, min_len=None, max_len=None):
    df = bed8_df.copy()
    if types:
        allowed = {t.strip().upper() for t in types.split(",") if t.strip()}
        df = df[df["type"].isin(allowed)]
    if min_len is not None:
        df = df[df["len"].fillna(0) >= float(min_len)]
    if max_len is not None and str(max_len).lower() != "none":
        df = df[df["len"].fillna(0) <= float(max_len)]
    df = df.reset_index(drop=True)
    return df


# ---------- Counting overlaps (observed) ----------

def count_observed(nlr_df, var_bed_df, tmpdir):
    """
    Returns:
      per_win_df: per-window DataFrame with total and per-type counts
      obs_total_sv: int
      obs_bp: total bp across windows
      per_type_obs_density: dict type-> variants/kb
    """
    a_bed = os.path.join(tmpdir, "nlr_windows.bed")
    write_bed6(nlr_df, a_bed)

    b_bed = os.path.join(tmpdir, "variants.bed8")
    write_bed(var_bed_df[["chrom","start","end","name","score","strand","type","len"]], b_bed)

    rows = run_bedtools_intersect_count(a_bed, b_bed)
    per_win = pd.DataFrame(rows, columns=["chrom","start","end","name","score","strand","n_sv"])
    per_win[["start","end"]] = per_win[["start","end"]].astype(int)
    per_win["length_bp"] = per_win["end"] - per_win["start"]

    wb = run_bedtools_intersect_wb(a_bed, b_bed)
    type_counts = defaultdict(Counter)  # window name -> type -> count
    for rec in wb:
        win_name = rec[3]
        vtype = rec[12].upper() if len(rec) >= 13 else "OTHER"
        type_counts[win_name][vtype] += 1

    # attach per-type columns (include SV + SNV classes)
    wanted_types = ["DEL","INS","INV","COMPLEX","SNP","MNP","OTHER"]
    for t in wanted_types:
        per_win[f"n_{t.lower()}"] = per_win["name"].map(lambda w: type_counts[w][t])

    obs_total = int(per_win["n_sv"].sum())
    obs_bp = int(per_win["length_bp"].sum())

    per_type_obs_density = {}
    for t in wanted_types:
        n = per_win[f"n_{t.lower()}"].sum()
        per_type_obs_density[t] = (n / obs_bp * 1e3) if obs_bp > 0 else 0.0

    return per_win, obs_total, obs_bp, per_type_obs_density


# ---------- Null sampling helpers ----------

def sample_one_window(chrom, length, chrom_max):
    L = int(chrom_max[chrom])
    if L <= length:
        start = 0
    else:
        start = random.randint(0, int(L - length))
    return chrom, start, start + length

def generate_null_windows(observed_df, chrom_max, chrom_matched, exclude_bed_path=None, max_attempts=100):
    """
    Returns DataFrame with same number of windows, matched lengths.
    If exclude_bed_path is provided, reject samples overlapping it.
    """
    obs = observed_df[["chrom","start","end","name"]].copy()
    obs["length"] = obs["end"] - obs["start"]

    chroms = list(chrom_max.keys())
    lengths = [int(x) for x in obs["length"].values]

    rows = []
    for i, (orig_chr, Lw, wname) in enumerate(zip(obs["chrom"], lengths, obs["name"])):
        attempts = 0
        fixed_chr = orig_chr if chrom_matched else None
        while attempts < max_attempts:
            c = fixed_chr if fixed_chr is not None else random.choice(chroms)
            s_chrom, s, e = sample_one_window(c, Lw, chrom_max)
            rows.append([s_chrom, s, e, f"null_{i}", 0, "."])
            if exclude_bed_path:
                with tempfile.NamedTemporaryFile("w", delete=False) as tmpbed:
                    tmpbed.write(f"{s_chrom}\t{s}\t{e}\tX\t0\t.\n")
                    cand_path = tmpbed.name
                try:
                    cmd = ["bedtools", "intersect", "-u", "-a", cand_path, "-b", exclude_bed_path]
                    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    if out.returncode == 0 and out.stdout.strip() == "":
                        os.unlink(cand_path)
                        break
                finally:
                    if os.path.exists(cand_path):
                        os.unlink(cand_path)
                rows.pop()  # reject and retry
                attempts += 1
            else:
                break
        # if exceeded attempts, keep last try as-is
    null_df = pd.DataFrame(rows, columns=["chrom","start","end","name","score","strand"])
    return null_df

def count_windows_vs_variants(windows_df, var_bed_df, tmpdir, per_type=False):
    a_bed = os.path.join(tmpdir, "A_windows.bed")
    b_bed = os.path.join(tmpdir, "B_variants.bed8")
    write_bed6(windows_df, a_bed)
    write_bed(var_bed_df[["chrom","start","end","name","score","strand","type","len"]], b_bed)

    if not per_type:
        rows = run_bedtools_intersect_count(a_bed, b_bed)
        df = pd.DataFrame(rows, columns=["chrom","start","end","name","score","strand","n_sv"])
        total = int(df["n_sv"].astype(int).sum()) if len(df) else 0
        bp = int((df["end"].astype(int) - df["start"].astype(int)).sum()) if len(df) else 0
        return total, bp, None
    else:
        wb = run_bedtools_intersect_wb(a_bed, b_bed)
        tcount = Counter()
        for rec in wb:
            vtype = rec[12].upper() if len(rec) >= 13 else "OTHER"
            tcount[vtype] += 1
        bp = int((windows_df["end"] - windows_df["start"]).sum())
        return int(sum(tcount.values())), bp, tcount


# ---------- Top-level worker (pickle-safe) ----------

def one_rep_worker(seed, observed_windows_df, var_bed_df, chrom_max, chrom_matched, exclude_bed_path, wanted_types):
    random.seed(seed)
    with tempfile.TemporaryDirectory(prefix="nullrep_") as td:
        null_df = generate_null_windows(
            observed_windows_df, chrom_max, chrom_matched,
            exclude_bed_path=exclude_bed_path
        )
        total, bp, tcount = count_windows_vs_variants(null_df, var_bed_df, td, per_type=True)
        dens = (total / bp * 1e3) if bp > 0 else 0.0
        t_dens = {t: ((tcount.get(t, 0) / bp * 1e3) if bp > 0 else 0.0) for t in wanted_types}
    return dens, t_dens


# ---------- Bootstrap driver (seedable, pickle-safe) ----------

def bootstrap_null(n_reps, threads, observed_windows_df, var_bed_df,
                   chrom_max, chrom_matched, exclude_windows_df=None, seeds=None):
    exclude_bed_path = None
    excl_tmpdir = None
    if exclude_windows_df is not None and len(exclude_windows_df):
        excl_tmpdir = tempfile.mkdtemp(prefix="excl_")
        exclude_bed_path = os.path.join(excl_tmpdir, "exclude.bed")
        write_bed6(exclude_windows_df, exclude_bed_path)

    wanted_types = ["DEL","INS","INV","COMPLEX","SNP","MNP","OTHER"]
    if seeds is None:
        seeds = [random.randint(0, 2**31 - 1) for _ in range(n_reps)]
    else:
        assert len(seeds) == n_reps, "seeds length must equal n_reps"

    results = []
    try:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            futs = [
                ex.submit(
                    one_rep_worker, s,
                    observed_windows_df, var_bed_df,
                    chrom_max, chrom_matched,
                    exclude_bed_path, wanted_types
                )
                for s in seeds
            ]
            for fut in as_completed(futs):
                results.append(fut.result())
    finally:
        if excl_tmpdir and os.path.isdir(excl_tmpdir):
            shutil.rmtree(excl_tmpdir, ignore_errors=True)

    null_dens = [r[0] for r in results]
    per_type = defaultdict(list)
    for _, tdict in results:
        for t, v in tdict.items():
            per_type[t].append(v)

    return np.array(null_dens), {t: np.array(vals) for t, vals in per_type.items()}


# ---------- Panel C helpers (from existing outputs only) ----------

def load_type_panel_from_outputs(outdir):
    """Builds per-type log2(obs/null_median) using observed_counts.tsv and per-type nulls."""
    obs_path  = os.path.join(outdir, "observed_counts.tsv")
    sv_null   = os.path.join(outdir, "bootstrap_null.tsv")
    snv_null  = os.path.join(outdir, "bootstrap_null_snv.tsv")
    if not (os.path.exists(obs_path) and os.path.exists(sv_null)):
        return None

    oc = pd.read_csv(obs_path, sep="\t")
    nu_sv = pd.read_csv(sv_null, sep="\t")
    nu_snv = pd.read_csv(snv_null, sep="\t") if os.path.exists(snv_null) else None

    # total bp across pangene/NLR windows
    total_bp = float((oc["end"] - oc["start"]).sum())
    if total_bp <= 0:
        return None

    # observed per-kb by type (use columns if present)
    type_cols = [
        ("DEL", "n_del"),
        ("INS", "n_ins"),
        ("INV", "n_inv"),
        ("COMPLEX", "n_complex"),
        ("SNP", "n_snp"),
        ("MNP", "n_mnp"),
        ("OTHER","n_other"),
    ]
    obs_perkb = {}
    for t, col in type_cols:
        if col in oc.columns:
            obs_perkb[t] = (oc[col].fillna(0).sum() / total_bp) * 1e3

    # null per-kb arrays by type; prefer SNV null for SNP/MNP if available
    def _null_arr(t):
        col = f"{t.lower()}_per_kb_null"
        if t in {"SNP","MNP"} and (nu_snv is not None) and (col in nu_snv.columns):
            arr = nu_snv[col].dropna().astype(float).values
            return arr if arr.size else None
        if col in nu_sv.columns:
            arr = nu_sv[col].dropna().astype(float).values
            return arr if arr.size else None
        return None

    rows = []
    for t in [t for t,_ in type_cols]:
        if t not in obs_perkb:
            continue
        arr = _null_arr(t)
        if arr is None:
            continue
        obs = float(obs_perkb[t])
        med = float(np.median(arr))
        if med == 0:
            fold = float("inf") if obs > 0 else 0.0
        else:
            fold = obs / med
        if np.isfinite(fold) and fold > 0:
            log2_fold = np.log2(fold)
        elif fold == 0:
            log2_fold = float("-inf")
        else:
            log2_fold = float("inf")
        rows.append([t, obs, med, fold, log2_fold])

    if not rows:
        return None

    df = pd.DataFrame(rows, columns=["type","obs_perkb","null_median_perkb","fold","log2_fold"])
    order = ["SNP","MNP","INS","DEL","INV","COMPLEX","OTHER"]
    df["type"] = pd.Categorical(df["type"], categories=order + [t for t in df["type"] if t not in order], ordered=True)
    df = df.sort_values("type")
    df.to_csv(os.path.join(outdir, "type_enrichment.tsv"), sep="\t", index=False)
    return df

def plot_type_enrichment_panel(ax, type_df, type_colors=None):

    if type_df is None or type_df.empty:
        ax.text(0.5, 0.5, "No per-type data", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    # keep only non-empty types
    df = type_df[type_df["obs_perkb"] > 0].copy()
    if df.empty:
        ax.text(0.5, 0.5, "All types are zero", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    # values as-is
    df["log2_fold"] = pd.to_numeric(df["log2_fold"], errors="coerce")

    # prep axes + colors
    cmap   = (type_colors or {})
    xlabs  = df["type"].astype(str).tolist()
    y      = df["log2_fold"].astype(float).tolist()
    colors = [cmap.get(t, None) for t in xlabs]

    ax.bar(range(len(xlabs)), y, color=colors)
    ax.axhline(0, linestyle="--", linewidth=1, color="black")
    ax.set_xticks(range(len(xlabs)), xlabs, rotation=45, ha="right")
    ax.set_ylabel("log2(obs / null median) per type")
    ax.set_xlabel("Variant type")

    # data-driven y-limits (no symmetric scaling)
    finite = [v for v in y if np.isfinite(v)]
    if finite:
        ymin = min(0.0, min(finite))
        ymax = max(finite)
        pad  = max(0.1, (ymax - ymin) * 0.1)
        ax.set_ylim(ymin, ymax + pad)

    # labels (above for +, below for -)
    for i, v in enumerate(y):
        if np.isfinite(v) and v != 0:
            ax.text(i, v + (0.02*(ax.get_ylim()[1]-ax.get_ylim()[0]) if v > 0 else -0.02*(ax.get_ylim()[1]-ax.get_ylim()[0])),
                    f"{v:.2f}", ha="center", va="bottom" if v > 0 else "top", fontsize=8, rotation=90)

    ax.set_title("c", loc="left", fontweight="bold", fontsize=14)
    ax.yaxis.grid(True, linestyle=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# ---------- Plotting (A/B panels) ----------

def make_figure_two_panels(
    outdir,
    sv_null_dens,
    sv_obs_dens,
    sv_pval,
    snv_null_dens=None,
    snv_obs_dens=None,
    snv_pval=None,
    prefix="figure_enrichment",
):
    """(kept for compatibility; used internally by the 3-panel wrapper)."""
    import numpy as np
    import os
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # --- Panel A: SV ---
    ax1.hist(sv_null_dens, bins=40, edgecolor="black")
    ax1.axvline(sv_obs_dens, linestyle="--", linewidth=2, color="black", label="Observed SV (NLR)")
    ax1.set_xlabel("SV per kb (null)")
    ax1.set_ylabel("Frequency")
    ax1.set_title("a", loc="left", fontweight="bold", fontsize=14)
    ax1.yaxis.grid(True, linestyle=":", alpha=0.5)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.legend(loc="upper left", fontsize=9, frameon=True)
    ax1.text(
        0.98, 0.96, f"p = {sv_pval:.4g}",
        transform=ax1.transAxes, ha="right", va="top",
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="black"), fontsize=10
    )

    # --- Panel B: SNV (SNP+MNP) ---
    ax2.set_title("b", loc="left", fontweight="bold", fontsize=14)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.yaxis.grid(True, linestyle=":", alpha=0.5)

    if snv_null_dens is not None and snv_obs_dens is not None:
        ax2.hist(snv_null_dens, bins=40, edgecolor="black")
        ax2.axvline(snv_obs_dens, linestyle="--", linewidth=2, color="black", label="Observed substitution (SNP + MNP) (NLR)")
        ax2.set_xlabel("SNP + MNP per kb (null)")
        ax2.set_ylabel("Frequency")
        ax2.legend(loc="upper left", fontsize=9, frameon=True)
        ax2.text(
            0.98, 0.96, f"p = {snv_pval:.4g}" if snv_pval is not None else "p = n/a",
            transform=ax2.transAxes, ha="right", va="top",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="black"), fontsize=10
        )
    else:
        ax2.set_xlabel("SNP + MNP per kb (null)")
        ax2.set_ylabel("Frequency")
        ax2.text(0.5, 0.5, "No SNV data provided", ha="center", va="center", transform=ax2.transAxes)

    plt.tight_layout()
    svg = os.path.join(outdir, f"{prefix}.svg")
    pdf = os.path.join(outdir, f"{prefix}.pdf")
    plt.savefig(svg, dpi=600, bbox_inches="tight")
    plt.savefig(pdf, dpi=600, bbox_inches="tight")
    plt.close()
    return svg, pdf

def make_figure_three_panels_from_outputs(
    outdir, sv_null_dens, sv_obs_dens, sv_pval,
    snv_null_dens=None, snv_obs_dens=None, snv_pval=None,
    prefix="figure_enrichment",
    type_colors=None,
):
    import matplotlib.pyplot as plt
    plt.close('all')
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

    # Panel A
    ax1.hist(sv_null_dens, bins=40, edgecolor="black")
    ax1.axvline(sv_obs_dens, linestyle="--", linewidth=2, color="black", label="Observed SV (NLR)")
    ax1.set_xlabel("SV per kb (null)"); ax1.set_ylabel("Frequency")
    ax1.set_title("a", loc="left", fontweight="bold", fontsize=14)
    ax1.yaxis.grid(True, linestyle=":", alpha=0.5); ax1.spines["top"].set_visible(False); ax1.spines["right"].set_visible(False)
    ax1.legend(loc="upper left", fontsize=9, frameon=True)
    ax1.text(0.98, 0.96, f"p = {sv_pval:.4g}", transform=ax1.transAxes,
             ha="right", va="top", bbox=dict(boxstyle="round", facecolor="white", edgecolor="black"), fontsize=10)

    # Panel B
    ax2.set_title("b", loc="left", fontweight="bold", fontsize=14)
    ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
    ax2.yaxis.grid(True, linestyle=":", alpha=0.5)
    if snv_null_dens is not None and snv_obs_dens is not None:
        ax2.hist(snv_null_dens, bins=40, edgecolor="black")
        ax2.axvline(snv_obs_dens, linestyle="--", linewidth=2, color="black", label="Observed substitution (SNP+MNP)")
        ax2.set_xlabel("SNP+MNP per kb (null)"); ax2.set_ylabel("Frequency")
        ax2.legend(loc="upper left", fontsize=9, frameon=True)
        ax2.text(0.98, 0.96, f"p = {snv_pval:.4g}" if snv_pval is not None else "p = n/a",
                 transform=ax2.transAxes, ha="right", va="top",
                 bbox=dict(boxstyle="round", facecolor="white", edgecolor="black"), fontsize=10)
    else:
        ax2.set_xlabel("SNP+MNP per kb (null)"); ax2.set_ylabel("Frequency")
        ax2.text(0.5, 0.5, "No SNV data provided", ha="center", va="center", transform=ax2.transAxes)

    # Panel C
    type_df = load_type_panel_from_outputs(outdir)
    plot_type_enrichment_panel(ax3, type_df, type_colors=type_colors)

    plt.tight_layout()
    svg = os.path.join(outdir, f"{prefix}.svg")
    pdf = os.path.join(outdir, f"{prefix}.pdf")
    plt.savefig(svg, dpi=600, bbox_inches="tight")
    plt.savefig(pdf, dpi=600, bbox_inches="tight")
    plt.close()
    return svg, pdf

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description="NLR SV & SNV enrichment using WAVE BED8 and pangene loci.")
    ap.add_argument("--hal-out", required=True, help="HAL output dir (contains pangene_loci.ref.bed6 and sv/)")
    ap.add_argument("--gene-list", required=True, help="Text file with gene IDs (one per line)")
    ap.add_argument("--outdir", required=True, help="Output directory")

    # path flexibility / overrides
    ap.add_argument("--var-out", default=None, help="Path to VAR_OUT (default HAL_OUT/sv)")
    ap.add_argument("--sv-bed",  default=None, help="Explicit path to wave.sv50.bed8 (overrides autodetect)")
    ap.add_argument("--snv-bed", default=None, help="Explicit path to wave.snp_mnp.bed8 (overrides autodetect)")
    ap.add_argument("--type-colors", default=None,
                help="Path to TSV/CSV/JSON mapping variant type->color (columns/keys: type,color).")


    # resources / settings
    ap.add_argument("--chrom-sizes", default=None, help="Path to chrom.sizes (default: VAR_OUT/chrom.sizes)")
    ap.add_argument("--bootstrap", type=int, default=2000, help="Bootstrap iterations (default 2000)")
    ap.add_argument("--threads", type=int, default=32, help="Parallel workers (default 32)")
    ap.add_argument("--bootstrap-seed", type=int, default=None, help="Base seed to sync SV/SNV bootstraps")
    ap.add_argument("--flank-kb", type=int, default=5, help="Flank for provenance (windows already expanded in logic)")
    ap.add_argument("--types", default=None, help="Comma list of SV types to include (e.g., DEL,INS,INV,COMPLEX). Default: all in SV50 BED8")
    ap.add_argument("--min-sv-bp", type=float, default=50, help="Min SV length (bp) filter (default 50)")
    ap.add_argument("--max-sv-bp", default=None, help="Max SV length (bp) filter (default None)")
    ap.add_argument("--chrom-matched", action="store_true", help="Match chromosome in null sampling")
    ap.add_argument("--exclude-nlr", action="store_true", help="Exclude observed NLR windows from null sampling")
    args = ap.parse_args()

    # tools
    check_tool("bedtools")

    # resolve dirs
    HAL_OUT = os.path.abspath(args.hal_out)
    VAR_OUT = os.path.abspath(args.var_out) if args.var_out else os.path.join(HAL_OUT, "sv")
    OUTDIR  = ensure_dir(os.path.abspath(args.outdir))

    # pangene loci (try both HAL_OUT and VAR_OUT)
    CAND_LOCI = [
        os.path.join(HAL_OUT, "pangene_loci.ref.bed6"),
        os.path.join(VAR_OUT, "pangene_loci.ref.bed6"),
    ]
    loci_bed6 = next((p for p in CAND_LOCI if os.path.exists(p)), None)
    if not loci_bed6:
        sys.exit(f"ERROR: pangene loci not found; tried: {CAND_LOCI}")

    # SV BED8
    CAND_SV_BEDS = [
        args.sv_bed,
        os.path.join(VAR_OUT, "bed", "wave.sv50.bed8"),
        os.path.join(VAR_OUT, "wave.sv50.bed8"),
    ]
    wave_bed8 = next((p for p in CAND_SV_BEDS if p and os.path.exists(p)), None)
    if not wave_bed8:
        sys.exit(f"ERROR: SV BED8 not found; tried: {CAND_SV_BEDS}")

    # SNV BED8 (optional; SNP+MNP)
    CAND_SNV_BEDS = [
        args.snv_bed,
        os.path.join(VAR_OUT, "bed", "wave.snp_mnp.bed8"),
        os.path.join(VAR_OUT, "wave.snp_mnp.bed8"),
    ]
    snv_bed8 = next((p for p in CAND_SNV_BEDS if p and os.path.exists(p)), None)
    include_snv = snv_bed8 is not None

    # chrom.sizes
    chrom_sizes = args.chrom_sizes or os.path.join(VAR_OUT, "chrom.sizes")
    if not os.path.exists(chrom_sizes):
        sys.exit(f"ERROR: missing chrom.sizes: {chrom_sizes}")

    print(f"[paths] HAL_OUT   = {HAL_OUT}")
    print(f"[paths] VAR_OUT   = {VAR_OUT}")
    print(f"[paths] loci_bed6 = {loci_bed6}")
    print(f"[paths] wave_bed8 = {wave_bed8}")
    print(f"[paths] chrom_sizes = {chrom_sizes}")
    print(f"[paths] snv_bed8 = {snv_bed8 if include_snv else '<none found>'}")

    # load inputs
    chrom_df = read_chrom_sizes(chrom_sizes)
    chrom_max = dict(zip(chrom_df["chrom"], chrom_df["length"].astype(int)))

    bed8_df = read_bed8(wave_bed8)
    bed8_df = filter_bed8_by_args(
        bed8_df,
        types=args.types,
        min_len=args.min_sv_bp,
        max_len=args.max_sv_bp
    )

    genome_bp_bg = int(chrom_df["length"].sum())
    sv_total_genome = int(len(bed8_df))
    sv_per_kb_genome = (sv_total_genome / genome_bp_bg * 1e3) if genome_bp_bg > 0 else 0.0

    # NLR windows from gene list
    gene2pangene = load_gene_to_pangene(HAL_OUT)
    genes = load_gene_list(args.gene_list)
    nlr_df = build_nlr_windows_from_genes(genes, gene2pangene, loci_bed6, chrom_df, flank_kb=args.flank_kb)
    if len(nlr_df) == 0:
        sys.exit("ERROR: zero NLR windows built (check gene list and mappings)")

    # ----- Observed SV counts/densities -----
    with tempfile.TemporaryDirectory(prefix="obs_sv_") as td:
        per_win_df, obs_total_sv, obs_bp, _per_type_obs_sv = count_observed(nlr_df, bed8_df, td)
    obs_sv_per_kb = (obs_total_sv / obs_bp * 1e3) if obs_bp > 0 else 0.0

    # Prepare shared seeds (so SV/SNV nulls sample the same windows)
    shared_seeds = None
    if args.bootstrap_seed is not None:
        rng = random.Random(args.bootstrap_seed)
        shared_seeds = [rng.randint(0, 2**31 - 1) for _ in range(args.bootstrap)]

    # ----- SV bootstrap null -----
    null_dens, per_type_null = bootstrap_null(
        n_reps=args.bootstrap,
        threads=args.threads,
        observed_windows_df=nlr_df,
        var_bed_df=bed8_df,
        chrom_max=chrom_max,
        chrom_matched=args.chrom_matched,
        exclude_windows_df=(nlr_df if args.exclude_nlr else None),
        seeds=shared_seeds
    )
    # right-tail p-value (SV)
    sv_pval = (1.0 + float((null_dens >= obs_sv_per_kb).sum())) / (len(null_dens) + 1.0)

    # ----- Optional SNV (SNP+MNP) observed + null -----
    per_win_snv_df = None
    snv_obs_per_kb = None
    snv_null_dens = None
    per_type_null_snv = {}
    snv_total_genome = 0
    snv_per_kb_genome = None
    snv_pval = None

    if include_snv:
        snv_df = read_bed8(snv_bed8)
        snv_df = snv_df[snv_df["type"].isin(["SNP","MNP"])].reset_index(drop=True)
        snv_total_genome = int(len(snv_df))
        snv_per_kb_genome = (snv_total_genome / genome_bp_bg * 1e3) if genome_bp_bg > 0 else None

        with tempfile.TemporaryDirectory(prefix="obs_snv_") as td:
            per_win_snv_df, obs_total_snv, obs_bp_snv, _per_type_obs_snv = count_observed(nlr_df, snv_df, td)
        snv_obs_per_kb = (obs_total_snv / obs_bp_snv * 1e3) if obs_bp_snv > 0 else 0.0

        snv_null_dens, per_type_null_snv = bootstrap_null(
            n_reps=args.bootstrap,
            threads=args.threads,
            observed_windows_df=nlr_df,
            var_bed_df=snv_df,
            chrom_max=chrom_max,
            chrom_matched=args.chrom_matched,
            exclude_windows_df=(nlr_df if args.exclude_nlr else None),
            seeds=shared_seeds
        )
        # right-tail p-value (SNV)
        snv_pval = (1.0 + float((snv_null_dens >= snv_obs_per_kb).sum())) / (len(snv_null_dens) + 1.0)

    # ----- Combine observed per-window counts and write table -----
    obs_cols = ["chrom","start","end","n_sv","n_del","n_ins","n_inv","n_complex","n_other","length_bp","name"]
    per_win_obs = per_win_df[obs_cols].copy()

    if include_snv and per_win_snv_df is not None:
        snv_cols = per_win_snv_df[["name", "n_snp", "n_mnp"]].copy()
        per_win_obs = per_win_obs.merge(snv_cols, on="name", how="left")
        per_win_obs[["n_snp","n_mnp"]] = per_win_obs[["n_snp","n_mnp"]].fillna(0).astype(int)
    else:
        per_win_obs["n_snp"] = 0
        per_win_obs["n_mnp"] = 0

    out_obs = per_win_obs[[
        "chrom","start","end",
        "n_sv","n_del","n_ins","n_inv","n_complex","n_other",
        "n_snp","n_mnp",
        "length_bp","name"
    ]].rename(columns={"name":"pangene_id"})
    out_obs_path = os.path.join(OUTDIR, "observed_counts.tsv")
    out_obs.to_csv(out_obs_path, sep="\t", index=False)

    # ----- Save SV null -----
    null_df = pd.DataFrame({"sv_per_kb_null": null_dens})
    for t, arr in per_type_null.items():
        null_df[f"{t.lower()}_per_kb_null"] = arr
    null_path = os.path.join(OUTDIR, "bootstrap_null.tsv")
    null_df.to_csv(null_path, sep="\t", index=False)

    # ----- Save SNV null if present -----
    if include_snv:
        null_snv_df = pd.DataFrame({"snv_per_kb_null": snv_null_dens})
        for t, arr in per_type_null_snv.items():
            null_snv_df[f"{t.lower()}_per_kb_null"] = arr
        null_snv_path = os.path.join(OUTDIR, "bootstrap_null_snv.tsv")
        null_snv_df.to_csv(null_snv_path, sep="\t", index=False)

    # ----- Density summary -----
    dens_row = {
        "nlr_windows_bp": int(obs_bp),
        "sv_in_nlr_windows": int(obs_total_sv),
        "sv_per_kb_nlr": float(obs_sv_per_kb),
        "genome_bp_bg": int(genome_bp_bg),
        "sv_total_genome": int(sv_total_genome),
        "sv_per_kb_genome": float(sv_per_kb_genome),
        "types_filter": args.types or "ALL",
        "min_sv_bp": args.min_sv_bp,
        "max_sv_bp": args.max_sv_bp
    }
    if include_snv and per_win_snv_df is not None:
        dens_row.update({
            "snv_in_nlr_windows": int(per_win_snv_df["n_sv"].sum()),
            "snv_per_kb_nlr": float(snv_obs_per_kb),
            "snv_total_genome": int(snv_total_genome),
            "snv_per_kb_genome": float(snv_per_kb_genome) if snv_per_kb_genome is not None else np.nan,
        })
    dens_summary = pd.DataFrame([dens_row])
    dens_summary_path = os.path.join(OUTDIR, "sv_density_summary.tsv")
    dens_summary.to_csv(dens_summary_path, sep="\t", index=False)

    # ----- Bootstrap summary text -----
    with open(os.path.join(OUTDIR, "bootstrap_summary.txt"), "w") as fh:
        fh.write(
            "Observed SV/kb in NLR windows: {:.6f}\n"
            "Iterations: {}\n"
            "Threads: {}\n"
            "Right-tail p-value (SV): {:.6g}\n"
            "{}"
            "Median window bp: {}\n"
            "Chrom-matched: {}\n"
            "Exclude-NLR from null: {}\n"
            "Filters: types={}, min_sv_bp={}, max_sv_bp={}\n"
            "HAL_OUT: {}\n"
            "VAR_OUT: {}\n"
            "SV_BED : {}\n"
            "{}".format(
                obs_sv_per_kb, args.bootstrap, args.threads, sv_pval,
                (f"Right-tail p-value (SNV): {snv_pval:.6g}\n" if include_snv and snv_pval is not None else ""),
                int(np.median(out_obs["length_bp"])),
                args.chrom_matched, args.exclude_nlr,
                args.types or "ALL", args.min_sv_bp, args.max_sv_bp,
                HAL_OUT, VAR_OUT, wave_bed8,
                (f"SNV observed per kb: {snv_obs_per_kb:.6f}\nSNV bed: {snv_bed8}\n" if include_snv and snv_obs_per_kb is not None else "")
            )
        )

    # ----- Figure (THREE panels; SV, SNV, per-type) -----
    type_colors = load_type_colors(args.type_colors)
    svg, pdf = make_figure_three_panels_from_outputs(
        outdir=OUTDIR,
        sv_null_dens=null_dens,
        sv_obs_dens=obs_sv_per_kb,
        sv_pval=sv_pval,
        snv_null_dens=(snv_null_dens if include_snv else None),
        snv_obs_dens=(snv_obs_per_kb if include_snv else None),
        snv_pval=(snv_pval if include_snv else None),
        type_colors=type_colors,
    )

    print("\nDone.")
    print(f"- observed_counts.tsv      → {out_obs_path}")
    print(f"- sv_density_summary.tsv   → {dens_summary_path}")
    print(f"- bootstrap_null.tsv       → {null_path}")
    if include_snv:
        print(f"- bootstrap_null_snv.tsv   → {null_snv_path}")
    print(f"- bootstrap_summary.txt    → {os.path.join(OUTDIR, 'bootstrap_summary.txt')}")
    print(f"- figures                  → {svg}, {pdf}")
    return 0

def main_replot():
    """
    Replot-only entry point.
    Reads previously generated TSVs in --outdir and writes a fresh THREE-panel figure:
      - Panel A: SV null vs observed (NLR), with right-tail p-value
      - Panel B: SNV (SNP+MNP) null vs observed, with right-tail p-value (if available)
      - Panel C: Per-type log2(obs/null_median) from observed_counts + bootstrap_null(_snv)
    """
    import argparse
    ap = argparse.ArgumentParser(description="Replot SV/SNV enrichment from existing outputs (no recomputation).")
    ap.add_argument("--outdir", required=True, help="Directory containing observed_counts.tsv, sv_density_summary.tsv, bootstrap_null*.tsv")
    ap.add_argument("--figure-prefix", default="figure_enrichment", help="Basename for output figure files (default: figure_enrichment)")
    ap.add_argument("--type-colors", default=None, help="Path to TSV/CSV/JSON mapping variant type->color (columns/keys: type,color).")
    args = ap.parse_args()

    OUTDIR = os.path.abspath(args.outdir)
    prefix = args.figure_prefix

    # ---- Load density summary (observed per-kb values) ----
    dens_path = os.path.join(OUTDIR, "sv_density_summary.tsv")
    if not os.path.exists(dens_path):
        sys.exit(f"ERROR: missing {dens_path}")
    dens = pd.read_csv(dens_path, sep="\t")
    if dens.empty:
        sys.exit(f"ERROR: {dens_path} is empty")

    row = dens.iloc[0]
    # Observed SV density (required)
    if "sv_per_kb_nlr" not in dens.columns:
        sys.exit("ERROR: sv_per_kb_nlr not found in sv_density_summary.tsv")
    sv_obs_per_kb = float(row["sv_per_kb_nlr"])

    # Observed SNV density (optional; compute fallback if missing)
    snv_obs_per_kb = None
    if "snv_per_kb_nlr" in dens.columns and pd.notna(row["snv_per_kb_nlr"]):
        snv_obs_per_kb = float(row["snv_per_kb_nlr"])
    else:
        obs_counts_path = os.path.join(OUTDIR, "observed_counts.tsv")
        if os.path.exists(obs_counts_path):
            oc = pd.read_csv(obs_counts_path, sep="\t")
            if {"n_snp", "n_mnp", "start", "end"}.issubset(set(oc.columns)):
                total_snv = int(oc["n_snp"].fillna(0).sum() + oc["n_mnp"].fillna(0).sum())
                total_bp = int((oc["end"] - oc["start"]).sum())
                if total_bp > 0:
                    snv_obs_per_kb = total_snv / total_bp * 1e3

    # ---- Load SV null distribution ----
    sv_null_path = os.path.join(OUTDIR, "bootstrap_null.tsv")
    if not os.path.exists(sv_null_path):
        sys.exit(f"ERROR: missing {sv_null_path}")
    svnull = pd.read_csv(sv_null_path, sep="\t")

    sv_null_col = "sv_per_kb_null" if "sv_per_kb_null" in svnull.columns else next(
        (c for c in svnull.columns if c.endswith("_per_kb_null")), None
    )
    if sv_null_col is None:
        sys.exit(f"ERROR: sv_per_kb_null column not found in {sv_null_path}")

    sv_null_dens = svnull[sv_null_col].dropna().astype(float).values
    if len(sv_null_dens) == 0:
        sys.exit(f"ERROR: {sv_null_path} has no null values in {sv_null_col}")

    # Right-tail p-value for SV with add-one smoothing
    sv_pval = (1.0 + float((sv_null_dens >= sv_obs_per_kb).sum())) / (len(sv_null_dens) + 1.0)

    # ---- Load SNV null distribution (optional) ----
    snv_null_path = os.path.join(OUTDIR, "bootstrap_null_snv.tsv")
    snv_null_dens = None
    snv_pval = None
    if os.path.exists(snv_null_path):
        snvnull = pd.read_csv(snv_null_path, sep="\t")
        snv_null_col = "snv_per_kb_null" if "snv_per_kb_null" in snvnull.columns else next(
            (c for c in snvnull.columns if c.endswith("_per_kb_null")), None
        )
        if snv_null_col is not None and snv_obs_per_kb is not None:
            snv_null_arr = snvnull[snv_null_col].dropna().astype(float).values
            if len(snv_null_arr) > 0:
                snv_null_dens = snv_null_arr
                snv_pval = (1.0 + float((snv_null_dens >= snv_obs_per_kb).sum())) / (len(snv_null_dens) + 1.0)

    # ----- Figure (THREE panels; SV, SNV, per-type) -----
    type_colors = load_type_colors(args.type_colors)
    svg, pdf = make_figure_three_panels_from_outputs(
        outdir=OUTDIR,
        sv_null_dens=sv_null_dens,
        sv_obs_dens=sv_obs_per_kb,
        sv_pval=sv_pval,
        snv_null_dens=snv_null_dens,
        snv_obs_dens=snv_obs_per_kb,
        snv_pval=snv_pval,
        prefix=prefix,
        type_colors=type_colors,
    )

    print("\nReplot complete.")
    print(f"- figures → {svg}, {pdf}")
    print(f"- SV observed per kb:  {sv_obs_per_kb:.6f} (n={len(sv_null_dens)}, p={sv_pval:.6g})")
    if snv_obs_per_kb is not None and snv_null_dens is not None:
        print(f"- SNV observed per kb: {snv_obs_per_kb:.6f} (n={len(snv_null_dens)}, p={snv_pval:.6g})")
    elif snv_obs_per_kb is not None:
        print(f"- SNV observed per kb: {snv_obs_per_kb:.6f} (no null file found)")
    else:
        print("- SNV panel skipped (no observed SNV density available)")
    return 0

if __name__ == "__main__":
    if "--replot" in sys.argv:
        sys.argv.remove("--replot")
        sys.exit(main_replot())
    else:
        sys.exit(main())