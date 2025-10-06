#!/usr/bin/env python3
"""
# This is not production code, but included for publication posterity.
Chromosome Gene Visualization Script - Proportional Scale

This script creates a visual representation of genes on a chromosome with their connections.
All chromosomes are displayed proportionally to a standard 1,000,000 bp reference scale,
so shorter chromosomes appear smaller and longer chromosomes appear larger.

Input file format:
Line 1: chromosome_name start_position end_position
Following lines: gene_name start_position end_position connected_gene1,connected_gene2,...

Usage: python chromosome_viz.py input_file.txt output_file.png
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, ConnectionPatch
import sys
import argparse
import numpy as np

# Standard scale reference (1 million base pairs)
STANDARD_SCALE = 1000000
# Visual length for standard scale (in plot units)
STANDARD_VISUAL_LENGTH = 10

def parse_input_file(filename):
    """Parse the input file and return chromosome info and gene data."""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    # Parse chromosome info (first line)
    chr_info = lines[0].split()
    chromosome = {
        'name': chr_info[0],
        'start': int(chr_info[1]),
        'end': int(chr_info[2])
    }
    
    # Parse gene info
    genes = {}
    connections = {}
    
    for line in lines[1:]:
        parts = line.split()
        gene_name = parts[0]
        start_pos = int(parts[1])
        end_pos = int(parts[2])
        
        genes[gene_name] = {
            'start': start_pos,
            'end': end_pos,
            'center': (start_pos + end_pos) / 2
        }
        
        # Parse connections if they exist
        if len(parts) > 3:
            connected_genes = [g.strip() for g in parts[3].split(',') if g.strip()]
            connections[gene_name] = connected_genes
    
    return chromosome, genes, connections

def scale_to_reference(chromosome, genes):
    """Scale chromosome and gene positions to proportional reference scale."""
    # Calculate actual chromosome length
    actual_length = chromosome['end'] - chromosome['start']
    
    # Calculate proportional visual length based on standard scale
    visual_scale_factor = STANDARD_VISUAL_LENGTH / STANDARD_SCALE
    chromosome_visual_length = actual_length * visual_scale_factor
    
    # Scale chromosome positions to start at 0
    scaled_chromosome = {
        'name': chromosome['name'],
        'start': 0,
        'end': chromosome_visual_length,
        'actual_start': chromosome['start'],
        'actual_end': chromosome['end'],
        'actual_length': actual_length,
        'visual_length': chromosome_visual_length
    }
    
    # Scale gene positions proportionally
    scaled_genes = {}
    for gene_name, gene_info in genes.items():
        # Convert to relative positions within chromosome
        relative_start = gene_info['start'] - chromosome['start']
        relative_end = gene_info['end'] - chromosome['start']
        
        # Scale proportionally to visual scale
        scaled_start = relative_start * visual_scale_factor
        scaled_end = relative_end * visual_scale_factor
        
        scaled_genes[gene_name] = {
            'start': scaled_start,
            'end': scaled_end,
            'center': (scaled_start + scaled_end) / 2,
            'actual_start': gene_info['start'],
            'actual_end': gene_info['end']
        }
    
    return scaled_chromosome, scaled_genes

def create_chromosome_plot(chromosome, genes, connections, output_file):
    """Create the chromosome visualization plot with proportional scaling."""
    # Scale positions proportionally to reference scale
    scaled_chromosome, scaled_genes = scale_to_reference(chromosome, genes)
    
    fig, ax = plt.subplots(1, 1, figsize=(16, 8))
    
    # Set up the plot area to accommodate full reference scale
    margin = STANDARD_VISUAL_LENGTH * 0.1
    ax.set_xlim(-margin, STANDARD_VISUAL_LENGTH + margin)
    ax.set_ylim(-3, 4)
    
    # Draw chromosome line (black) - only as long as actual chromosome
    chromosome_y = 1
    ax.plot([scaled_chromosome['start'], scaled_chromosome['end']], [chromosome_y, chromosome_y], 
            'k-', linewidth=4, solid_capstyle='round')
    
    # Draw gene blocks (blue)
    gene_height = 0.3
    gene_positions = {}  # Store gene positions for connections
    
    for gene_name, gene_info in scaled_genes.items():
        # Create gene block starting exactly at gene start position
        gene_width = gene_info['end'] - gene_info['start']
        
        # Position gene block at the front of the chromosome line
        gene_y = chromosome_y + 0.05  # Slightly above chromosome line
        
        # Use FancyBboxPatch for rounded rectangles
        gene_block = FancyBboxPatch(
            (gene_info['start'], gene_y - gene_height/2),
            gene_width, gene_height,
            boxstyle="round,pad=0.02",
            facecolor='steelblue',
            edgecolor='darkblue',
            linewidth=1.5
        )
        ax.add_patch(gene_block)
        
        # Add gene name on top of block
        ax.text(gene_info['center'], gene_y + gene_height/2 + 0.2, 
                gene_name, ha='center', va='bottom', fontsize=6, 
                fontweight='bold', color='black', rotation=90)
        
        # Store position for connections
        gene_positions[gene_name] = {
            'center': gene_info['center'],
            'top': gene_y + gene_height/2,
            'bottom': gene_y - gene_height/2
        }
    
    # Draw connections (curved gray lines)
    for gene_name, connected_genes in connections.items():
        if gene_name not in gene_positions:
            continue
            
        for connected_gene in connected_genes:
            if connected_gene not in gene_positions:
                continue
                
            # Create curved connection from bottom of gene blocks
            start_x = gene_positions[gene_name]['center']
            end_x = gene_positions[connected_gene]['center']
            
            # Use ConnectionPatch for curved lines downward
            con = ConnectionPatch(
                (start_x, gene_positions[gene_name]['bottom']),
                (end_x, gene_positions[connected_gene]['bottom']),
                "data", "data",
                arrowstyle="-", shrinkA=0, shrinkB=0,
                connectionstyle="arc3,rad=0.15",  # Positive rad with bottom connection creates downward curve
                color='gray', linewidth=2, alpha=0.7
            )
            ax.add_patch(con)
    
    # Draw intron connections (dashed gray lines)
    intron_pairs = find_intron_pairs(scaled_genes)
    
    for gene1, gene2 in intron_pairs:
        if gene1 in gene_positions and gene2 in gene_positions:
            start_x = gene_positions[gene1]['center']
            end_x = gene_positions[gene2]['center']
            
            # Draw straight line above chromosome
            intron_y = chromosome_y + 0.3
            ax.plot([start_x, end_x], [intron_y, intron_y], 
                   color='gray', linewidth=1.2, linestyle='--', alpha=0.8)
            
            # Add small vertical connectors
            ax.plot([start_x, start_x], 
                   [gene_positions[gene1]['top'], intron_y], 
                   color='gray', linewidth=1.2, linestyle='--', alpha=0.8)
            ax.plot([end_x, end_x], 
                   [gene_positions[gene2]['top'], intron_y], 
                   color='gray', linewidth=1.2, linestyle='--', alpha=0.8)

    # Add reference scale bar and chromosome-specific scale
    add_reference_scale(ax, scaled_chromosome, chromosome_y - 2.5)
    
    # Set labels and title
    actual_length_kb = scaled_chromosome['actual_length'] / 1000
    percentage_of_reference = (scaled_chromosome['actual_length'] / STANDARD_SCALE) * 100
    
    ax.set_title(f'Gene Map - Chromosome {scaled_chromosome["name"]} '
                f'(Length: {actual_length_kb:.1f} kb, {percentage_of_reference:.1f}% of 1Mb reference)', 
                fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Position (proportional to 1Mb reference scale)', fontsize=12, fontweight='bold')
    
    # Remove y-axis ticks and labels
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], color='black', linewidth=6, label='Chromosome'),
        patches.Rectangle((0, 0), 1, 1, facecolor='steelblue', label='Genes'),
        plt.Line2D([0], [0], color='gray', linewidth=2, label='Gene Connections'),
        plt.Line2D([0], [0], color='gray', linewidth=1.2, linestyle='--', label='Intron Connections'),
        plt.Line2D([0], [0], color='lightgray', linewidth=2, linestyle=':', label='1Mb Reference Scale')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=500, bbox_inches='tight', format='png', transparent=True)
    print(f"Chromosome visualization saved to: {output_file}")
    print(f"Chromosome length: {actual_length_kb:.1f} kb ({percentage_of_reference:.1f}% of 1Mb reference)")
    print(f"Visual length: {scaled_chromosome['visual_length']:.2f} units (reference = {STANDARD_VISUAL_LENGTH} units)")

def find_intron_pairs(genes):
    """Find gene pairs that represent introns (e.g., Gene1.1 and Gene1.2)."""
    intron_pairs = []
    gene_names = list(genes.keys())
    
    for gene in gene_names:
        if '.' in gene:
            base_name = gene.rsplit('.', 1)[0]
            suffix = gene.rsplit('.', 1)[1]
            
            # Look for corresponding intron pair
            if suffix == '1':
                pair_name = f"{base_name}.2"
                if pair_name in genes:
                    intron_pairs.append((gene, pair_name))
    
    return intron_pairs

def add_reference_scale(ax, scaled_chromosome, y_position):
    """Add reference scale showing 1Mb scale and actual chromosome scale."""
    
    # Draw full reference scale line (light gray, dashed)
    reference_y = y_position - 0.8
    ax.plot([0, STANDARD_VISUAL_LENGTH], [reference_y, reference_y], 
           color='lightgray', linewidth=2, linestyle=':', alpha=0.8)
    ax.text(-0.5, reference_y + 0.1, '1Mb Reference Scale', ha='right', va='bottom', 
           fontsize=10, fontweight='bold', color='gray', alpha=0.8)
    
    # Add reference scale ticks (every 100kb)
    ref_tick_interval = STANDARD_VISUAL_LENGTH / 10  # 10 ticks for 1Mb = 100kb each
    for i in range(11):  # 0 to 10 inclusive
        tick_pos = i * ref_tick_interval
        ax.plot([tick_pos, tick_pos], [reference_y - 0.05, reference_y + 0.05], 
               color='lightgray', linewidth=1, alpha=0.8)
        if i % 2 == 0:  # Label every 200kb
            label = f'{i * 100}k'
            ax.text(tick_pos, reference_y - 0.15, label, ha='center', va='top', 
                   fontsize=8, color='gray', alpha=0.8)
    
    # Draw actual chromosome scale (darker)
    chromosome_y = y_position
    actual_length = scaled_chromosome['actual_length']
    
    # Calculate appropriate tick interval for chromosome
    if actual_length <= 10000:
        tick_interval_bp = 1000  # 1kb
    elif actual_length <= 100000:
        tick_interval_bp = 10000  # 10kb
    elif actual_length <= 1000000:
        tick_interval_bp = 100000  # 100kb
    else:
        tick_interval_bp = 1000000  # 1Mb
    
    # Convert tick interval to visual units
    visual_scale_factor = STANDARD_VISUAL_LENGTH / STANDARD_SCALE
    tick_interval_visual = tick_interval_bp * visual_scale_factor
    
    # Draw chromosome scale bar
    ax.plot([0, scaled_chromosome['end']], [chromosome_y, chromosome_y], 
           color='black', linewidth=2)
    ax.text(-0.5, chromosome_y + 0.1, f'Chromosome Scale (actual: {actual_length/1000:.0f}kb)', 
           ha='right', va='bottom', fontsize=10, fontweight='bold', color='black')
    
    # Add chromosome scale ticks
    current_pos = 0
    tick_count = 0
    while current_pos <= scaled_chromosome['end'] and tick_count < 20:  # Limit ticks
        ax.plot([current_pos, current_pos], [chromosome_y - 0.05, chromosome_y + 0.05], 
               color='black', linewidth=1.5)
        
        # Calculate actual position for label
        actual_pos = scaled_chromosome['actual_start'] + (current_pos / visual_scale_factor)
        if tick_interval_bp >= 1000:
            label = f'{actual_pos/1000:.0f}k'
        else:
            label = f'{actual_pos:.0f}'
            
        ax.text(current_pos, chromosome_y + 0.15, label, ha='center', va='bottom', 
               fontsize=8, color='black', rotation=45)
        
        current_pos += tick_interval_visual
        tick_count += 1

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(
        description='Create proportionally-scaled chromosome gene visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Example input file format:
Chr1 1000 10000
Gene1.1 1500 2000 Gene2,Gene3
Gene1.2 2500 3000
Gene2 4000 4500 Gene1.1
Gene3 6000 6800 Gene1.1,Gene2
Gene4 8000 8300

All chromosomes will be displayed proportionally to a {STANDARD_SCALE:,} bp reference scale.
- A 70,000 bp chromosome will appear as 7% the length of the full reference scale
- A 500,000 bp chromosome will appear as 50% the length of the full reference scale
- A 2,000,000 bp chromosome will appear as 200% the length (extending beyond reference)

Features:
- Proportional sizing relative to 1Mb reference
- Reference scale bar showing full 1Mb scale
- Chromosome-specific scale bar with actual positions
- Percentage of reference scale shown in title
- Curved connections between connected genes
- Dashed intron connections between gene pairs
        """
    )
    
    parser.add_argument('input_file', help='Input file containing chromosome and gene data')
    parser.add_argument('output_file', help='Output PNG file name')
    
    args = parser.parse_args()
    
    try:
        # Parse input file
        chromosome, genes, connections = parse_input_file(args.input_file)
        
        # Create visualization
        create_chromosome_plot(chromosome, genes, connections, args.output_file)
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()