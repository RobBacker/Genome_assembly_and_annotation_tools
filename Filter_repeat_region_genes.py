#!/usr/bin/env python3
"""
Filter_repeat_region_genes.py (fixed version)

Usage:
    Filter_repeat_region_genes.py \
        --gff helixer.agat.aed.gff3 \
        --repeat-gff repeatmasker.gff \
        --repeat-threshold 25 \
        --aed-threshold 0.75 \
        --out-keep keep_genes.txt \
        --out-drop drop_genes.txt \
        --out-stats gene_stats.tsv

Filters gene models by repeat overlap, using AED for rescue, with correct CDS merging
and de-duplicated repeat overlap calculations.
"""
import argparse
from collections import defaultdict
import sys


def parse_attrs(attrs):
    """Convert GFF attribute string into dict."""
    d = {}
    for part in attrs.strip().split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            d[k] = v
    return d


def load_gff3_data(gff3):
    """
    Parse GFF3 to collect CDS and AED per gene.
    Returns:
        gene_cds: dict gid -> list of (chrom, start0, end)
        gene_aed: dict gid -> float
    """
    gene_map = {}
    gene_cds = defaultdict(list)
    gene_aed = {}

    with open(gff3) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip().split('\t')
            if len(cols) < 9:
                continue
            feat = cols[2]
            attrs = parse_attrs(cols[8])

            if feat in ('mRNA', 'transcript'):
                tid = attrs.get('ID')
                gid = attrs.get('Parent')
                if tid and gid:
                    gene_map[tid] = gid
                    aed_str = attrs.get('aed_ev_tr')
                    if aed_str:
                        try:
                            gene_aed[gid] = float(aed_str)
                        except ValueError:
                            pass

            elif feat == 'CDS':
                tid = attrs.get('Parent')
                gid = gene_map.get(tid)
                if gid:
                    chrom = cols[0]
                    s = int(cols[3]) - 1
                    e = int(cols[4])
                    gene_cds[gid].append((chrom, s, e))

    return gene_cds, gene_aed


def load_repeats(gff_file):
    """
    Load repeat regions from a GFF file into a dictionary:
    chrom -> list of (start0, end)
    """
    rep_map = defaultdict(list)
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.rstrip().split('\t')
            if len(cols) < 5:
                continue

            feature = cols[2].lower()
            if feature not in {'dispersed_repeat', 'match', 'match_part', 'repeat_region', 'repeat'}:
                continue

            chrom = cols[0]
            try:
                start = int(cols[3]) - 1
                end = int(cols[4])
            except ValueError:
                continue
            rep_map[chrom].append((start, end))

    return rep_map


def merge_intervals(intervals):
    """Merge a list of (start,end) tuples and return merged list."""
    if not intervals:
        return []
    sorted_int = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_int[0]]
    for s, e in sorted_int[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))
    return merged


def compute_metrics(gene_cds, rep_map, gene_aed, repeat_threshold, aed_threshold):
    """
    For each gene, merge CDS, compute total_cds_bp, repeat_bp, pct_repeat, expressed flag.
    Returns dict of gene->metrics.
    """
    stats = {}
    for gid, exs in gene_cds.items():
        # group by chromosome
        chrom_map = defaultdict(list)
        for chrom, s, e in exs:
            chrom_map[chrom].append((s, e))

        total_bp = 0
        repeat_bp = 0
        for chrom, ivs in chrom_map.items():
            # merge CDS intervals
            merged_cds = merge_intervals(ivs)
            total_bp += sum(e - s for s, e in merged_cds)

            # collect overlaps with repeats
            overlaps = []
            reps = rep_map.get(chrom, [])
            for s, e in merged_cds:
                for rs, re in reps:
                    if re <= s or rs >= e:
                        continue
                    ov_s = max(s, rs)
                    ov_e = min(e, re)
                    overlaps.append((ov_s, ov_e))

            # merge repeat overlaps to avoid double-counting
            merged_repeats = merge_intervals(overlaps)
            repeat_bp += sum(end - start for start, end in merged_repeats)

        pct = (repeat_bp / total_bp * 100) if total_bp > 0 else 0
        aed = gene_aed.get(gid)
        expressed = 1 if (aed is not None and aed < aed_threshold) else 0
        keep = (expressed == 1 or pct <= repeat_threshold)
        stats[gid] = {
            'total_cds_bp': total_bp,
            'repeat_bp': repeat_bp,
            'pct_repeat': pct,
            'aed_ev_tr': aed,
            'expressed': expressed,
            'keep': keep
        }
    return stats


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', required=True, help='Input GFF3 with CDS and AED tags')
    parser.add_argument('--repeat-gff', required=True, help='Repeat GFF file')
    parser.add_argument('--repeat-threshold', type=float, default=25.0,
                        help='Maximum percent repeat overlap to drop a gene')
    parser.add_argument('--aed-threshold', type=float, default=0.75,
                        help='AED threshold below which genes are always kept')
    parser.add_argument('--out-keep', required=True, help='Output file for kept gene IDs')
    parser.add_argument('--out-drop', required=True, help='Output TSV for dropped genes')
    parser.add_argument('--out-stats', required=True, help='Output TSV of stats for all genes')
    args = parser.parse_args()

    sys.stderr.write('Loading GFF data...\n')
    gene_cds, gene_aed = load_gff3_data(args.gff)
    sys.stderr.write(f'  Loaded {len(gene_cds)} genes\n')

    sys.stderr.write('Loading repeats...\n')
    rep_map = load_repeats(args.repeat_gff)
    # optional: pre-merge repeat intervals per chromosome to speed overlap checks
    for chrom in list(rep_map):
        rep_map[chrom] = merge_intervals(rep_map[chrom])
    sys.stderr.write(f'  Loaded repeats on {len(rep_map)} chromosomes\n')

    sys.stderr.write('Computing metrics...\n')
    stats = compute_metrics(
        gene_cds, rep_map, gene_aed,
        args.repeat_threshold, args.aed_threshold
    )

    # Write stats
    with open(args.out_stats, 'w') as fs:
        fs.write('gene_id\ttotal_cds_bp\trepeat_bp\tpct_repeat\taed_ev_tr\texpressed\tkeep\n')
        for gid, m in stats.items():
            fs.write(
                f"{gid}\t{m['total_cds_bp']}\t{m['repeat_bp']}\t"
                f"{m['pct_repeat']:.2f}\t{m['aed_ev_tr']}\t{m['expressed']}\t{int(m['keep'])}\n"
            )

    # Write keep/drop
    with open(args.out_keep, 'w') as fk, open(args.out_drop, 'w') as fd:
        fd.write('gene_id\ttotal_cds_bp\trepeat_bp\tpct_repeat\taed_ev_tr\texpressed\n')
        for gid, m in stats.items():
            if m['keep']:
                fk.write(gid + '\n')
            else:
                fd.write(
                    f"{gid}\t{m['total_cds_bp']}\t{m['repeat_bp']}\t"
                    f"{m['pct_repeat']:.2f}\t{m['aed_ev_tr']}\t{m['expressed']}\n"
                )

    sys.stderr.write(
        f"Kept {sum(1 for m in stats.values() if m['keep'])} genes; "
        f"dropped {sum(1 for m in stats.values() if not m['keep'])} genes\n"
    )
