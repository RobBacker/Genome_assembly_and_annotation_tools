#!/usr/bin/env python3
"""
Pangenome orthogroup and gene support analysis using GFF3 curation tags.

Generates:
 1) PAV_orthogroups.tsv: orthogroup category and list of assemblies present
 2) Annotated <species>_gene_support_with_category.tsv: adds orthogroup category to each gene row with support status
 3) Stacked bar plots showing orthogroup and gene category distributions
 4) Optional: analysis of user-supplied gene list across assemblies

Usage:
    python3 pangenome_pav_support.py \
        --counts Orthogroups.GeneCount.tsv \
        --orthogroups Orthogroups.tsv \
        --input_dir path/to/search/_for_gff3_files \
        --output_dir results/ \
        [--gene_list my_genes.txt --target_outdir gene_list_results]
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = ['sans-serif']
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'Bitstream Vera Sans', 'Tahoma']
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import re
import argparse
import os

CATEGORIES_GENE = ['Core', 'Soft-core', 'Shell', 'Cloud-supported', 'Cloud-unsupported']
CATEGORIES_ORTHO = ['Core', 'Soft-core', 'Shell', 'Cloud']
COLORS = ['#2ca02c', '#ff7f0e', '#1f77b4', '#8F00FF', '#A9A9A9']

def categorize_orthogroup(count, num_species):
    proportion = count / num_species
    if proportion >= 0.95:   # Present in ≥95% of species
        return 'Core'
    elif proportion >= 0.85: # Present in 85%-95% of species
        return 'Soft-core'
    elif proportion >= 0.20: # Present in 20%-70% of species
        return 'Shell'
    else:                   # Present in <20% of species
        return 'Cloud'

def create_pav_table(counts_df, output_path, omit_set=None):
    if 'Total' in counts_df.columns:
        counts_df = counts_df.drop(columns=['Total'])
    if omit_set:
        drop_cols = [col for col in counts_df.columns if col.split('.')[0] in omit_set]
        counts_df = counts_df.drop(columns=drop_cols)
    new_cols = [counts_df.columns[0]] + [c.split('.')[0] for c in counts_df.columns[1:]]
    counts_df.columns = new_cols
    species_cols = new_cols[1:]
    num_species = len(species_cols)
    presence = counts_df[species_cols] >= 1
    counts_df['PresenceCount'] = presence.sum(axis=1)
    counts_df['Category'] = counts_df['PresenceCount'].apply(lambda x: categorize_orthogroup(x, num_species))
    counts_df['Assemblies'] = presence.apply(lambda row: ",".join([sp for sp, present in row.items() if present]), axis=1)
    pav_df = counts_df.rename(columns={counts_df.columns[0]: 'Orthogroup'})[['Orthogroup', 'Category', 'Assemblies']]
    pav_df = pav_df[pav_df['Assemblies'].str.strip() != '']
    pav_df.to_csv(output_path, sep='\t', index=False)
    print(f"Wrote PAV table to {output_path}")
    return pav_df

def strip_isoform_suffix(gene_id):
    """Remove isoform suffixes like .1, .2 at the end of gene IDs"""
    return re.sub(r'\.\d+$', '', gene_id)

def convert_gff3_filename_to_assembly_key(fname):
    """
    Convert GFF3 filename like 'Hass.primary.merged.aed.curated.isotag.gff3' 
    into standardized assembly key 'Hass_primary'
    """
    parts = fname.split('.')
    if len(parts) > 1:
        base = parts[0]
        suffix = parts[1]
        if suffix in ['primary', 'alternate']:
            return f"{base}_{suffix}"
    # fallback to base only if pattern unexpected
    return parts[0]

def parse_gff3_for_supported_mrna_ids(gff3_file):
    """
    Parse a GFF3 file to extract mRNA IDs with curation=cat1 to cat4 as supported.
    Returns a set of supported mRNA IDs including isoform suffixes.
    """
    supported_mrna_ids = set()
    curation_re = re.compile(r'curation=cat([1-4])')
    id_re = re.compile(r'ID=([^;]+)')
    with open(gff3_file, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type != 'mRNA':
                continue
            attributes = parts[8]
            curation_match = curation_re.search(attributes)
            id_match = id_re.search(attributes)
            if curation_match and id_match:
                supported_mrna_ids.add(id_match.group(1))
    return supported_mrna_ids

def find_supported_genes_per_assembly_from_gff3(input_dir, omit_set):
    """
    Walk the input directory looking for GFF3 files with names matching '*.merged.aed.curated.isotag.gff3'.
    Extract supported mRNA IDs per assembly keyed by the converted assembly key.
    """
    supported_per_assembly = {}
    for root, dirs, files in os.walk(input_dir):
        for fname in files:
            if fname.endswith(".merged.aed.curated.isotag.gff3"):
                assembly_name = convert_gff3_filename_to_assembly_key(fname)
                if assembly_name in omit_set:
                    continue
                filepath = os.path.join(root, fname)
                supported_ids = parse_gff3_for_supported_mrna_ids(filepath)
                supported_per_assembly[assembly_name] = supported_ids
                print(f"Assembly '{assembly_name}': found {len(supported_ids)} supported mRNA IDs")
    return supported_per_assembly

def create_species_gene_support(orthogroups_file, input_dir, output_dir, pav_df, omit_set):
    """
    Create per-species gene support annotation by checking genes against supported sets parsed from GFF3.
    Keep separate support info for each assembly column like 'Hass_primary' and 'Hass_alternate'.
    """
    os.makedirs(output_dir, exist_ok=True)
    supported_per_assembly = find_supported_genes_per_assembly_from_gff3(input_dir, omit_set)
    og_df = pd.read_csv(orthogroups_file, sep='\t')
    header = og_df.columns.tolist()
    new_header = [header[0]] + [h.split('.')[0] for h in header[1:]]
    og_df.columns = new_header
    species_cols = [sp for sp in new_header[1:] if sp not in omit_set]
    og_to_category = pav_df.set_index('Orthogroup')['Category'].to_dict()
    category_matrices = {sp: {cat: 0 for cat in CATEGORIES_GENE} for sp in species_cols}

    for species in species_cols:
        records = []
        # Look up supported set directly by assembly key == species name (e.g. Hass_primary)
        supported_set = supported_per_assembly.get(species, set())
        supported_normalized = set(strip_isoform_suffix(g) for g in supported_set)

        for _, row in og_df.iterrows():
            og = row['Orthogroup']
            cell = row.get(species, '')
            if pd.isna(cell) or str(cell).strip() == '':
                continue
            genes = [g.strip() for g in re.split(r',\s*', str(cell)) if g.strip()]
            for gene in genes:
                gene_norm = strip_isoform_suffix(gene)
                is_supported = gene_norm in supported_normalized
                base_category = og_to_category.get(og, 'NA')
                if base_category == 'Cloud':
                    full_category = 'Cloud-supported' if is_supported else 'Cloud-unsupported'
                else:
                    full_category = base_category
                support_status = 'supported' if is_supported else 'unsupported'
                records.append((gene, og, support_status, base_category))
                if full_category in category_matrices[species]:
                    category_matrices[species][full_category] += 1

        df_out = pd.DataFrame(records, columns=['Gene', 'Orthogroup', 'Support', 'Category'])
        out_file = os.path.join(output_dir, f"{species}_gene_support_with_category.tsv")
        df_out.to_csv(out_file, sep='\t', index=False)
        print(f"Wrote gene support with category for {species} to {out_file}")

    matrix_df = pd.DataFrame(category_matrices).fillna(0).astype(int).loc[CATEGORIES_GENE]
    plot_stacked(
        matrix_df,
        'Gene Count Distribution by Support Category',
        'gene',
        save_path=os.path.join(output_dir, 'gene_support_category_distribution.png'),
        categories=CATEGORIES_GENE
    )
    matrix_df.T.to_csv(os.path.join(output_dir, 'gene_support_category_counts.tsv'), sep='\t')

def create_orthogroup_category_support(pav_df, output_dir):
    pav_df_clean = pav_df[pav_df['Assemblies'].notna() & (pav_df['Assemblies'].str.strip() != '')]
    all_assemblies = set()
    for asm_list in pav_df_clean['Assemblies']:
        all_assemblies.update(asm_list.split(','))
    assemblies = sorted(all_assemblies)
    category_counts = {asm: {cat: 0 for cat in CATEGORIES_ORTHO} for asm in assemblies}
    for _, row in pav_df_clean.iterrows():
        cat = row['Category']
        present_in = row['Assemblies'].split(',')
        for asm in present_in:
            if cat in category_counts.get(asm, {}):
                category_counts[asm][cat] += 1
    matrix_df = pd.DataFrame(category_counts).fillna(0).astype(int).loc[CATEGORIES_ORTHO]
    matrix_df.T.to_csv(os.path.join(output_dir, 'orthogroup_category_counts.tsv'), sep='\t')
    plot_stacked(
        matrix_df,
        'Orthogroup Category Distribution',
        'orthogroup',
        save_path=os.path.join(output_dir, 'orthogroup_category_distribution.png'),
        categories=CATEGORIES_ORTHO
    )
    return matrix_df

def plot_stacked(mat, title, analysis_type, save_path=None, categories=None):
    if categories is None:
        categories = mat.index.tolist()
    fig, ax = plt.subplots(figsize=(12, 8))
    x = np.arange(len(mat.columns))
    bottom = np.zeros_like(x, dtype=float)
    for cat, color in zip(categories, COLORS):
        if cat not in mat.index:
            continue
        vals = mat.loc[cat].values
        ax.bar(x, vals, bottom=bottom, color=color, edgecolor='none', linewidth=0)
        bottom += vals
    ax.set_ylabel('Gene Count' if analysis_type=='gene' else 'Number of Orthogroups', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(mat.columns, rotation=45, ha='right', fontsize=10)
    handles, labels = [], []
    for cat, color in zip(categories, COLORS):
        if cat not in mat.index:
            continue
        lines = [rf"$\mathbf{{{cat}}}$"]
        for sp in mat.columns:
            lines.append(f"    {sp} [{mat.loc[cat, sp]:,}]")
        label = "\n".join(lines)
        handles.append(Patch(facecolor=color, edgecolor='white'))
        labels.append(label)
    ax.legend(
        handles, labels,
        title='Category Breakdown',
        bbox_to_anchor=(1.01, 0.5),
        loc='center left',
        fontsize=7,
        title_fontsize=9,
        borderaxespad=0.5,
        frameon=True
    )
    plt.subplots_adjust(bottom=0.15, top=0.9)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax*1.1)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
        base, _ = os.path.splitext(save_path)
        svg_path = base + ".svg"
        plt.savefig(svg_path, format='svg', bbox_inches='tight')
        print(f"Saved figure to {svg_path}")

    plt.close(fig)

def prompt_omit_assemblies(species_list):
    print("\nAssemblies/species found:")
    for i, sp in enumerate(species_list, 1):
        print(f"  {i}. {sp}")
    omit = input("\nEnter comma-separated numbers (or names) of assemblies to omit, or press Enter to keep all: ").strip()
    omit_set = set()
    if omit:
        for val in omit.split(','):
            val = val.strip()
            if val.isdigit():
                idx = int(val) - 1
                if 0 <= idx < len(species_list):
                    omit_set.add(species_list[idx])
            elif val in species_list:
                omit_set.add(val)
    return omit_set

def find_orthogroups_for_genes_ignore_isoform(gene_list, orthogroups_file):
    og_df = pd.read_csv(orthogroups_file, sep='\t')
    header = og_df.columns.tolist()
    new_header = [header[0]] + [h.split('.')[0] for h in header[1:]]
    og_df.columns = new_header
    gene_list_no_iso = set([strip_isoform_suffix(g) for g in gene_list])
    gene_to_og = {}
    assembly_for_gene = {}
    for assembly in new_header[1:]:
        for _, row in og_df.iterrows():
            genes = [g.strip() for g in str(row[assembly]).split(',') if g.strip()]
            genes_no_iso = [strip_isoform_suffix(g) for g in genes]
            for gene, gene_no_iso in zip(genes, genes_no_iso):
                if gene_no_iso in gene_list_no_iso:
                    gene_to_og[gene_no_iso] = row[new_header[0]]
                    assembly_for_gene[gene_no_iso] = assembly
    if not gene_to_og:
        return None, []
    assemblies = set(assembly_for_gene.values())
    if len(assemblies) > 1:
        print(f"Warning: Input genes found in multiple assemblies: {', '.join(assemblies)}")
    source_assembly = list(assemblies)[0]
    orthogroups = list(set(gene_to_og.values()))
    return source_assembly, orthogroups

def analyze_orthogroup_presence(orthogroups, orthogroups_file, output_dir, omit_set, gene_list_filename=None):
    og_df = pd.read_csv(orthogroups_file, sep='\t')
    header = og_df.columns.tolist()
    new_header = [header[0]] + [h.split('.')[0] for h in header[1:]]
    og_df.columns = new_header
    assemblies = [col for col in new_header[1:] if col not in omit_set]

    presence_data = []
    for og in orthogroups:
        og_row = og_df[og_df[new_header[0]] == og]
        if og_row.empty:
            continue
        og_row = og_row.iloc[0]
        counts = {asm: len([g for g in str(og_row[asm]).split(',') if g.strip()]) for asm in assemblies}
        counts['Orthogroup'] = og
        presence_data.append(counts)

    presence_df = pd.DataFrame(presence_data).set_index('Orthogroup')
    presence_df['Total'] = presence_df.sum(axis=1)
    presence_df_sorted = presence_df.sort_values(['Total', 'Orthogroup'], ascending=[False, True]).drop(columns=['Total'])

    COLORS = plt.cm.tab20.colors
    PATTERNS = ['', '///', 'xxx', '..........']
    n_ogs = presence_df_sorted.shape[0]

    color_cycle = (COLORS * ((n_ogs // len(COLORS)) + 1))[:n_ogs]
    hatch_cycle = []
    for i in range(n_ogs):
        if i >= len(COLORS):
            pattern = PATTERNS[(i // len(COLORS)) % (len(PATTERNS)-1) + 1]
            hatch_cycle.append(pattern)
        else:
            hatch_cycle.append(PATTERNS[0])

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = presence_df_sorted.T.plot(kind='bar', stacked=True, ax=ax, color=color_cycle, edgecolor='black', linewidth=0.2)

    for idx, patch_list in enumerate(ax.containers):
        hatch = hatch_cycle[idx]
        for patch in patch_list:
            patch.set_hatch(hatch)
            patch.set_edgecolor('black')
            patch.set_linewidth(0.2)

    legend_elements = []
    for idx, og in enumerate(presence_df_sorted.index):
        color = color_cycle[idx]
        hatch = hatch_cycle[idx]
        legend_elements.append(Patch(facecolor=color, edgecolor='black', hatch=hatch, label=og, linewidth=0.2))

    plot_title = 'Gene Counts in Target Orthogroups per Assembly'
    if gene_list_filename:
        plot_title += f' ({os.path.basename(gene_list_filename)})'
    plt.title(plot_title)
    plt.ylabel('Number of Genes')
    plt.xlabel('Assembly')
    plt.xticks(rotation=45)

    ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.01, 0.5),
        loc='center left',
        fontsize=7,
        title='Orthogroup',
        title_fontsize=9,
        borderaxespad=0.5,
        frameon=True,
        ncol=2
    )

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'orthogroup_gene_counts.png'), dpi=300, bbox_inches='tight')
    plt.close()
    presence_df_sorted.to_csv(os.path.join(output_dir, 'orthogroup_presence_counts.tsv'), sep='\t')
    return presence_df_sorted

def create_binary_presence_matrix(counts_df, output_path, omit_set=None):
    """
    Create a binary presence/absence matrix for orthogroups across assemblies.
    Writes a table with columns: Orthogroup, <assemblies...>, count, category.
    """
    if 'Total' in counts_df.columns:
        counts_df = counts_df.drop(columns=['Total'])
    if omit_set:
        drop_cols = [col for col in counts_df.columns if col.split('.')[0] in omit_set]
        counts_df = counts_df.drop(columns=drop_cols)

    new_cols = [counts_df.columns[0]] + [c.split('.')[0] for c in counts_df.columns[1:]]
    counts_df.columns = new_cols
    species_cols = new_cols[1:]
    num_species = len(species_cols)

    # binary presence/absence
    presence = (counts_df[species_cols] >= 1).astype(int)
    presence.insert(0, 'Orthogroup', counts_df[new_cols[0]])
    presence['count'] = presence[species_cols].sum(axis=1)
    presence['category'] = presence['count'].apply(lambda x: categorize_orthogroup(x, num_species))

    presence.to_csv(output_path, sep='\t', index=False)
    print(f"Wrote binary PAV matrix to {output_path}")
    return presence

def main():
    parser = argparse.ArgumentParser(description="Generate PAV and per-species gene support tables from GFF3 curation tags")
    parser.add_argument('--counts',      required=True, help='Path to Orthogroups.GeneCount.tsv')
    parser.add_argument('--orthogroups', required=True, help='Path to Orthogroups.tsv')
    parser.add_argument('--input_dir',   required=True, help='Root directory to search for GFF3 files')
    parser.add_argument('--output_dir',  default='.', help='Directory to save outputs')
    parser.add_argument('--gene_list', help='File containing list of target genes (one per line)')
    parser.add_argument('--target_outdir', default='orthogroup_analysis', help='Output directory for orthogroup analysis')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    counts_df = pd.read_csv(args.counts, sep='\t')
    species_cols = [c.split('.')[0] for c in counts_df.columns if c not in ['Orthogroup', 'Total']]
    omit_set = prompt_omit_assemblies(species_cols)
    if omit_set:
        print(f"Omitting assemblies: {', '.join(omit_set)}")

    pav_all_path = os.path.join(args.output_dir, 'PAV_all.tsv')
    pav_omit_path = os.path.join(args.output_dir, 'PAV_orthogroups.tsv')
    pav_df_all = create_pav_table(counts_df.copy(), pav_all_path, omit_set=None)
    pav_df = create_pav_table(counts_df.copy(), pav_omit_path, omit_set)
    binary_matrix_path = os.path.join(args.output_dir, 'PAV_binary_matrix.tsv')
    create_binary_presence_matrix(counts_df.copy(), binary_matrix_path, omit_set)


    if args.gene_list:
        with open(args.gene_list) as f:
            gene_list = [line.strip() for line in f if line.strip()]
        source_asm, orthogroups = find_orthogroups_for_genes_ignore_isoform(gene_list, args.orthogroups)
        if not orthogroups:
            print("Warning: No genes from the input list were found in any assembly!")
        else:
            print(f"Found {len(orthogroups)} orthogroups containing input genes (source assembly: {source_asm})")
            os.makedirs(args.target_outdir, exist_ok=True)
            analyze_orthogroup_presence(orthogroups, args.orthogroups, args.target_outdir, omit_set, gene_list_filename=args.gene_list)

    create_orthogroup_category_support(pav_df, args.output_dir)
    create_species_gene_support(args.orthogroups, args.input_dir, args.output_dir, pav_df, omit_set)

if __name__ == '__main__':
    main()



