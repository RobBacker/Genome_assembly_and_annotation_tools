#!/usr/bin/env python3
"""
# This is not production code, but included for publication posterity.
GFF3 NLR and TE Density Analyzer

This script reads three GFF3 files (main annotations, NLRs, TEs) and creates:
1. Chromosome visualization with NLR density heatmap (left) and TE density heatmap (right)
2. Correlation plot between NLR and TE densities

Usage: python gff3_density_analyzer.py main.gff3 nlrs.gff3 tes.gff3 output_prefix
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from collections import defaultdict
import argparse
import sys
import re
from scipy.stats import pearsonr
from scipy.ndimage import gaussian_filter1d

def parse_gff3_file(filename):
    """Parse GFF3 file and return features as a list of dictionaries."""
    features = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            if len(fields) < 9:
                continue
                
            feature = {
                'seqid': fields[0],
                'source': fields[1],
                'type': fields[2],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'score': fields[5],
                'strand': fields[6],
                'phase': fields[7],
                'attributes': fields[8]
            }
            features.append(feature)
    
    return features

def get_chromosome_boundaries(main_features):
    """Determine chromosome boundaries based on first and last gene positions."""
    chromosomes = {}
    
    # Group features by chromosome/contig
    chr_features = defaultdict(list)
    for feature in main_features:
        if feature['type'].lower() in ['gene', 'mrna', 'cds']:  # Consider gene-like features
            chr_features[feature['seqid']].append(feature)
    
    # Calculate boundaries for each chromosome
    for chr_name, features in chr_features.items():
        if not features:
            continue
            
        starts = [f['start'] for f in features]
        ends = [f['end'] for f in features]
        
        chromosomes[chr_name] = {
            'start': min(starts),
            'end': max(ends),
            'length': max(ends) - min(starts)
        }
    
    return chromosomes

def calculate_density(features, chromosomes, window_size=10000):
    """Calculate feature density across chromosomes using sliding windows."""
    densities = {}
    
    for chr_name, chr_info in chromosomes.items():
        chr_start = chr_info['start']
        chr_end = chr_info['end']
        chr_length = chr_end - chr_start
        
        # Create windows
        num_windows = max(1, chr_length // window_size)
        window_centers = []
        window_densities = []
        
        for i in range(num_windows):
            window_start = chr_start + (i * window_size)
            window_end = min(chr_start + ((i + 1) * window_size), chr_end)
            window_center = (window_start + window_end) / 2
            
            # Count features in this window
            count = 0
            for feature in features:
                if (feature['seqid'] == chr_name and 
                    feature['start'] <= window_end and 
                    feature['end'] >= window_start):
                    count += 1
            
            # Density per kb
            window_length_kb = (window_end - window_start) / 1000
            density = count / window_length_kb if window_length_kb > 0 else 0
            
            window_centers.append(window_center)
            window_densities.append(density)
        
        # Smooth the density using gaussian filter
        if len(window_densities) > 3:
            smoothed_densities = gaussian_filter1d(window_densities, sigma=1)
        else:
            smoothed_densities = window_densities
            
        densities[chr_name] = {
            'positions': window_centers,
            'densities': smoothed_densities,
            'raw_densities': window_densities
        }
    
    return densities

def create_density_visualization(chromosomes, nlr_densities, te_densities, output_file):
    """Create chromosome visualization with density heatmaps."""
    
    # Sort chromosomes by length (largest first)
    sorted_chrs = sorted(chromosomes.items(), key=lambda x: x[1]['length'], reverse=True)
    
    fig, ax = plt.subplots(1, 1, figsize=(16, 2 + len(sorted_chrs) * 0.8))
    
    # Color maps
    nlr_cmap = LinearSegmentedColormap.from_list("nlr", ["white", "red", "darkred"])
    te_cmap = LinearSegmentedColormap.from_list("te", ["white", "blue", "darkblue"])
    
    # Calculate max densities for normalization
    all_nlr_densities = []
    all_te_densities = []
    for chr_name in chromosomes:
        if chr_name in nlr_densities:
            all_nlr_densities.extend(nlr_densities[chr_name]['densities'])
        if chr_name in te_densities:
            all_te_densities.extend(te_densities[chr_name]['densities'])
    
    max_nlr_density = max(all_nlr_densities) if all_nlr_densities else 1
    max_te_density = max(all_te_densities) if all_te_densities else 1
    
    # Set up plot dimensions
    chr_height = 0.6
    chr_spacing = 1.0
    heatmap_width = 0.3
    
    y_positions = []
    chr_labels = []
    
    for i, (chr_name, chr_info) in enumerate(sorted_chrs):
        y_pos = len(sorted_chrs) - i - 1
        y_positions.append(y_pos)
        chr_labels.append(chr_name)
        
        chr_start = chr_info['start']
        chr_end = chr_info['end']
        chr_length = chr_end - chr_start
        
        # Scale chromosome length for display (relative to longest chromosome)
        max_length = sorted_chrs[0][1]['length']
        display_length = 10 * (chr_length / max_length)  # Max 10 units wide
        
        # Draw chromosome backbone
        chr_y = y_pos * chr_spacing
        ax.plot([0, display_length], [chr_y, chr_y], 'k-', linewidth=6, solid_capstyle='round')
        
        # Draw NLR density heatmap (left side)
        if chr_name in nlr_densities:
            positions = nlr_densities[chr_name]['positions']
            densities = nlr_densities[chr_name]['densities']
            
            for j, (pos, density) in enumerate(zip(positions, densities)):
                # Convert position to display coordinates
                rel_pos = (pos - chr_start) / chr_length * display_length
                
                # Normalize density for color mapping
                norm_density = density / max_nlr_density if max_nlr_density > 0 else 0
                color = nlr_cmap(norm_density)
                
                # Draw heatmap rectangle
                rect = patches.Rectangle(
                    (rel_pos - 0.05, chr_y + 0.1), 
                    0.1, heatmap_width,
                    facecolor=color, edgecolor='none', alpha=0.8
                )
                ax.add_patch(rect)
        
        # Draw TE density heatmap (right side)
        if chr_name in te_densities:
            positions = te_densities[chr_name]['positions']
            densities = te_densities[chr_name]['densities']
            
            for j, (pos, density) in enumerate(zip(positions, densities)):
                # Convert position to display coordinates
                rel_pos = (pos - chr_start) / chr_length * display_length
                
                # Normalize density for color mapping
                norm_density = density / max_te_density if max_te_density > 0 else 0
                color = te_cmap(norm_density)
                
                # Draw heatmap rectangle
                rect = patches.Rectangle(
                    (rel_pos - 0.05, chr_y - 0.1 - heatmap_width), 
                    0.1, heatmap_width,
                    facecolor=color, edgecolor='none', alpha=0.8
                )
                ax.add_patch(rect)
        
        # Add chromosome label
        ax.text(-0.5, chr_y, chr_name, ha='right', va='center', 
               fontsize=10, fontweight='bold')
        
        # Add length annotation
        ax.text(display_length + 0.2, chr_y, f'{chr_length/1000:.0f}kb', 
               ha='left', va='center', fontsize=8, color='gray')
    
    # Set plot limits and labels
    ax.set_xlim(-2, 12)
    ax.set_ylim(-0.5, len(sorted_chrs) * chr_spacing)
    
    ax.set_title('Chromosome Map with NLR and TE Density', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Relative Position', fontsize=12)
    
    # Remove axes
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    
    # Add legends
    # NLR legend
    nlr_legend_y = len(sorted_chrs) * chr_spacing - 0.5
    ax.text(-1.8, nlr_legend_y, 'NLR Density', ha='center', va='bottom', 
           fontsize=10, fontweight='bold', color='darkred', rotation=90)
    
    # Create NLR colorbar
    for i in range(10):
        intensity = i / 9
        color = nlr_cmap(intensity)
        rect = patches.Rectangle((-2.1, nlr_legend_y - 2 + i*0.2), 0.2, 0.2,
                               facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(rect)
    
    ax.text(-2.3, nlr_legend_y - 2, '0', ha='right', va='center', fontsize=8)
    ax.text(-2.3, nlr_legend_y - 0.2, f'{max_nlr_density:.1f}', ha='right', va='center', fontsize=8)
    
    # TE legend
    ax.text(-1.4, nlr_legend_y, 'TE Density', ha='center', va='bottom', 
           fontsize=10, fontweight='bold', color='darkblue', rotation=90)
    
    # Create TE colorbar
    for i in range(10):
        intensity = i / 9
        color = te_cmap(intensity)
        rect = patches.Rectangle((-1.7, nlr_legend_y - 2 + i*0.2), 0.2, 0.2,
                               facecolor=color, edgecolor='black', linewidth=0.5)
        ax.add_patch(rect)
    
    ax.text(-1.9, nlr_legend_y - 2, '0', ha='right', va='center', fontsize=8)
    ax.text(-1.9, nlr_legend_y - 0.2, f'{max_te_density:.1f}', ha='right', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='png')
    print(f"Chromosome density visualization saved to: {output_file}")

def create_correlation_plot(nlr_densities, te_densities, output_file):
    """Create correlation plot between NLR and TE densities."""
    
    # Collect paired density values
    nlr_values = []
    te_values = []
    
    for chr_name in nlr_densities:
        if chr_name in te_densities:
            nlr_chr = nlr_densities[chr_name]['densities']
            te_chr = te_densities[chr_name]['densities']
            
            # Match positions (assuming same window size)
            min_len = min(len(nlr_chr), len(te_chr))
            nlr_values.extend(nlr_chr[:min_len])
            te_values.extend(te_chr[:min_len])
    
    if not nlr_values or not te_values:
        print("Warning: No overlapping density data found for correlation analysis")
        return
    
    # Calculate correlation
    correlation, p_value = pearsonr(nlr_values, te_values)
    
    # Create scatter plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Create scatter plot with some transparency
    ax.scatter(nlr_values, te_values, alpha=0.6, c='blue', s=30, edgecolors='black', linewidth=0.5)
    
    # Add trend line
    z = np.polyfit(nlr_values, te_values, 1)
    p = np.poly1d(z)
    x_trend = np.linspace(min(nlr_values), max(nlr_values), 100)
    ax.plot(x_trend, p(x_trend), "r--", linewidth=2, alpha=0.8)
    
    # Set labels and title
    ax.set_xlabel('NLR Density (per kb)', fontsize=12, fontweight='bold')
    ax.set_ylabel('TE Density (per kb)', fontsize=12, fontweight='bold')
    ax.set_title(f'NLR vs TE Density Correlation\n'
                f'r = {correlation:.3f}, p = {p_value:.2e}', 
                fontsize=14, fontweight='bold')
    
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='png')
    print(f"Correlation plot saved to: {output_file}")
    print(f"Correlation coefficient: {correlation:.3f} (p-value: {p_value:.2e})")

def main():
    """Main function to run the analysis."""
    parser = argparse.ArgumentParser(
        description='Analyze NLR and TE density from GFF3 files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script creates two visualizations:
1. Chromosome map with NLR density (red, left) and TE density (blue, right) heatmaps
2. Scatter plot showing correlation between NLR and TE densities

Input files:
- main_gff3: Main annotation file (used to determine chromosome boundaries)
- nlr_gff3: NLR-specific annotations
- te_gff3: TE-specific annotations

Output files:
- [output_prefix]_chromosome_density.png: Chromosome visualization
- [output_prefix]_correlation.png: Correlation plot

Example usage:
python gff3_density_analyzer.py annotations.gff3 nlrs.gff3 tes.gff3 analysis
        """
    )
    
    parser.add_argument('main_gff3', help='Main GFF3 annotation file')
    parser.add_argument('nlr_gff3', help='NLR-specific GFF3 file')
    parser.add_argument('te_gff3', help='TE-specific GFF3 file')
    parser.add_argument('output_prefix', help='Prefix for output PNG files')
    parser.add_argument('--window-size', type=int, default=50000,
                       help='Window size for density calculation (default: 50000)')
    
    args = parser.parse_args()
    
    try:
        print("Reading GFF3 files...")
        
        # Parse GFF3 files
        main_features = parse_gff3_file(args.main_gff3)
        nlr_features = parse_gff3_file(args.nlr_gff3)
        te_features = parse_gff3_file(args.te_gff3)
        
        print(f"Loaded {len(main_features)} main features")
        print(f"Loaded {len(nlr_features)} NLR features")
        print(f"Loaded {len(te_features)} TE features")
        
        # Determine chromosome boundaries
        print("Calculating chromosome boundaries...")
        chromosomes = get_chromosome_boundaries(main_features)
        print(f"Found {len(chromosomes)} chromosomes/contigs")
        
        for chr_name, chr_info in sorted(chromosomes.items()):
            print(f"  {chr_name}: {chr_info['start']}-{chr_info['end']} ({chr_info['length']/1000:.0f}kb)")
        
        # Calculate densities
        print("Calculating NLR densities...")
        nlr_densities = calculate_density(nlr_features, chromosomes, args.window_size)
        
        print("Calculating TE densities...")
        te_densities = calculate_density(te_features, chromosomes, args.window_size)
        
        # Create visualizations
        print("Creating chromosome density visualization...")
        density_output = f"{args.output_prefix}_chromosome_density.png"
        create_density_visualization(chromosomes, nlr_densities, te_densities, density_output)
        
        print("Creating correlation plot...")
        correlation_output = f"{args.output_prefix}_correlation.png"
        create_correlation_plot(nlr_densities, te_densities, correlation_output)
        
        print("Analysis complete!")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()