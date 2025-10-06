#!/usr/bin/env python3
"""
# This is not production code, but included for publication posterity.
OrthoFinder Species Phylogenetic Tree Pipeline

This script automates the generation of species phylogenetic trees using OrthoFinder.
It handles the complete workflow from protein FASTA files to final tree visualization.

Requirements:
- OrthoFinder (orthofinder command in PATH)
- Python packages: biopython, ete3, matplotlib, pandas
- Optional: FastTree, IQ-TREE, or RAxML for custom tree inference

Usage: python3 orthofinder_phylogeny.py <protein_dir> [options]

Author: Phylogenetic Analysis Pipeline
"""

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import re
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set Arial font for plots
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']

def check_dependencies():
    """Check if required software is available."""
    
    print("Checking dependencies...")
    
    required_commands = ['orthofinder']
    optional_commands = ['fasttree', 'iqtree', 'raxmlHPC']
    
    missing_required = []
    available_optional = []
    
    for cmd in required_commands:
        if not shutil.which(cmd):
            missing_required.append(cmd)
        else:
            print(f"✓ {cmd} found")
    
    for cmd in optional_commands:
        if shutil.which(cmd):
            available_optional.append(cmd)
            print(f"✓ {cmd} found (optional)")
    
    if missing_required:
        print(f"Error: Missing required software: {', '.join(missing_required)}")
        print("Please install OrthoFinder: https://github.com/davidemms/OrthoFinder")
        sys.exit(1)
    
    print(f"Optional tree builders available: {available_optional if available_optional else 'None (will use default)'}")
    return available_optional

def validate_protein_files(protein_dir):
    """Validate input protein FASTA files."""
    
    protein_path = Path(protein_dir)
    if not protein_path.exists():
        print(f"Error: Directory {protein_dir} does not exist")
        sys.exit(1)
    
    # Find FASTA files
    fasta_extensions = ['.fasta', '.fa', '.fas', '.faa', '.pep']
    fasta_files = []
    
    for ext in fasta_extensions:
        fasta_files.extend(protein_path.glob(f'*{ext}'))
    
    if not fasta_files:
        print(f"Error: No protein FASTA files found in {protein_dir}")
        print(f"Expected extensions: {', '.join(fasta_extensions)}")
        sys.exit(1)
    
    print(f"Found {len(fasta_files)} protein files:")
    
    # Validate each file
    valid_files = []
    species_names = []
    
    for fasta_file in fasta_files:
        try:
            # Count sequences
            seq_count = 0
            with open(fasta_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    seq_count += 1
                    if seq_count > 10000:  # Stop counting after 10k for speed
                        break
            
            if seq_count == 0:
                print(f"Warning: {fasta_file.name} appears to be empty")
                continue
            
            # Extract species name from filename
            species_name = fasta_file.stem
            species_names.append(species_name)
            valid_files.append(fasta_file)
            
            print(f"  {fasta_file.name}: {seq_count}{'+' if seq_count >= 10000 else ''} proteins")
            
        except Exception as e:
            print(f"Warning: Could not read {fasta_file.name}: {e}")
    
    if len(valid_files) < 3:
        print("Error: Need at least 3 valid protein files for phylogenetic analysis")
        sys.exit(1)
    
    print(f"✓ {len(valid_files)} valid protein files ready for analysis")
    return valid_files, species_names

def run_orthofinder(protein_dir, output_dir, threads, tree_method, additional_args):
    """Run OrthoFinder analysis."""
    
    print(f"\nRunning OrthoFinder analysis...")
    print(f"Input: {protein_dir}")
    print(f"Output: {output_dir}")
    print(f"Threads: {threads}")
    print(f"Tree method: {tree_method}")
    
    # Build OrthoFinder command
    cmd = [
        'orthofinder',
        '-f', str(protein_dir),
        '-o', str(output_dir),
        '-t', str(threads),
        '-a', str(threads)  # Number of parallel alignment tasks
    ]
    
    # Add tree method if specified
    if tree_method:
        if tree_method.lower() == 'fasttree':
            cmd.extend(['-M', 'fasttree'])
        elif tree_method.lower() == 'iqtree':
            cmd.extend(['-M', 'iqtree'])
        elif tree_method.lower() == 'raxml':
            cmd.extend(['-M', 'raxml'])
    
    # Add any additional arguments
    if additional_args:
        cmd.extend(additional_args.split())
    
    print(f"Command: {' '.join(cmd)}")
    
    try:
        # Run OrthoFinder
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1
        )
        
        # Print output in real-time
        print("\n" + "="*60)
        print("ORTHOFINDER OUTPUT:")
        print("="*60)
        
        for line in process.stdout:
            print(line.rstrip())
        
        process.wait()
        
        if process.returncode != 0:
            print(f"Error: OrthoFinder failed with return code {process.returncode}")
            sys.exit(1)
        
        print("="*60)
        print("✓ OrthoFinder completed successfully")
        
    except FileNotFoundError:
        print("Error: OrthoFinder not found in PATH")
        sys.exit(1)
    except Exception as e:
        print(f"Error running OrthoFinder: {e}")
        sys.exit(1)

def find_orthofinder_results(output_dir):
    """Find OrthoFinder results directory."""
    
    output_path = Path(output_dir)
    
    # OrthoFinder creates timestamped directories
    result_dirs = list(output_path.glob('Results_*'))
    
    if not result_dirs:
        print(f"Error: No OrthoFinder results found in {output_dir}")
        sys.exit(1)
    
    # Get the most recent results
    latest_results = max(result_dirs, key=lambda x: x.stat().st_mtime)
    print(f"Using results from: {latest_results.name}")
    
    return latest_results

def extract_species_tree(results_dir):
    """Extract and analyze the species tree from OrthoFinder results."""
    
    print("\nExtracting species tree...")
    
    # Look for species tree file
    species_tree_dir = results_dir / "Species_Tree"
    
    if not species_tree_dir.exists():
        print("Error: Species_Tree directory not found in results")
        return None
    
    # Find tree files
    tree_files = list(species_tree_dir.glob("SpeciesTree_rooted*.txt"))
    
    if not tree_files:
        # Try unrooted tree
        tree_files = list(species_tree_dir.glob("SpeciesTree*.txt"))
    
    if not tree_files:
        print("Error: No species tree files found")
        return None
    
    # Use the first tree file found
    tree_file = tree_files[0]
    print(f"Found species tree: {tree_file.name}")
    
    # Read tree
    with open(tree_file, 'r') as f:
        tree_string = f.read().strip()
    
    return tree_file, tree_string

def analyze_orthofinder_stats(results_dir):
    """Analyze OrthoFinder statistics."""
    
    print("\nAnalyzing OrthoFinder statistics...")
    
    stats = {}
    
    # Statistics from Statistics_Overall directory
    stats_dir = results_dir / "Comparative_Genomics_Statistics"
    
    if stats_dir.exists():
        # Look for statistics files
        stats_files = {
            'orthogroups': list(stats_dir.glob("*Orthogroups*Statistics*.tsv")),
            'species_stats': list(stats_dir.glob("*Statistics_PerSpecies*.tsv")),
            'duplications': list(stats_dir.glob("*Duplications*.tsv"))
        }
        
        for stat_type, files in stats_files.items():
            if files:
                try:
                    df = pd.read_csv(files[0], sep='\t')
                    stats[stat_type] = df
                    print(f"✓ Loaded {stat_type} statistics: {df.shape}")
                except Exception as e:
                    print(f"Warning: Could not load {stat_type}: {e}")
    
    # Gene count statistics
    orthogroups_file = results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    if orthogroups_file.exists():
        try:
            df = pd.read_csv(orthogroups_file, sep='\t', index_col=0)
            stats['gene_counts'] = df
            print(f"✓ Loaded gene count matrix: {df.shape}")
        except Exception as e:
            print(f"Warning: Could not load gene counts: {e}")
    
    return stats

def visualize_species_tree(tree_file, tree_string, output_prefix):
    """Visualize species tree using matplotlib (simple method)."""
    
    print(f"Visualizing species tree...")
    
    try:
        # Try to import ete3 for better tree visualization
        from ete3 import Tree, TreeStyle, NodeStyle
        
        # Parse tree
        tree = Tree(tree_string, format=1)
        
        # Create tree style
        ts = TreeStyle()
        ts.show_leaf_names = True
        ts.show_branch_length = True
        ts.show_branch_support = True
        ts.branch_vertical_margin = 10
        ts.show_scale = True
        
        # Node styling
        for node in tree.traverse():
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "black"
            nstyle["size"] = 8
            if node.is_leaf():
                nstyle["fgcolor"] = "blue"
            node.set_style(nstyle)
        
        # Save tree images
        png_file = f"{output_prefix}_species_tree.png"
        pdf_file = f"{output_prefix}_species_tree.pdf"
        svg_file = f"{output_prefix}_species_tree.svg"
        
        tree.render(png_file, tree_style=ts, dpi=300)
        tree.render(pdf_file, tree_style=ts)
        tree.render(svg_file, tree_style=ts)
        
        print(f"✓ Tree visualizations saved: {png_file}, {pdf_file}, {svg_file}")
        
        # Also save the tree in Newick format
        tree_out = f"{output_prefix}_species_tree.newick"
        with open(tree_out, 'w') as f:
            f.write(tree_string)
        print(f"✓ Tree saved in Newick format: {tree_out}")
        
        return True
        
    except ImportError:
        print("ete3 not available, using simple text output")
        
        # Simple text-based tree representation
        tree_out = f"{output_prefix}_species_tree.newick"
        with open(tree_out, 'w') as f:
            f.write(tree_string)
        
        # Create a simple species list
        species_list = re.findall(r'([^(),:\s]+):', tree_string)
        
        txt_file = f"{output_prefix}_species_list.txt"
        with open(txt_file, 'w') as f:
            f.write("Species in phylogenetic tree:\n")
            f.write("=" * 30 + "\n")
            for i, species in enumerate(sorted(set(species_list)), 1):
                f.write(f"{i:2d}. {species}\n")
        
        print(f"✓ Tree saved: {tree_out}")
        print(f"✓ Species list saved: {txt_file}")
        return False

def create_summary_report(results_dir, stats, tree_file, output_prefix, args):
    """Create a comprehensive summary report."""
    
    print("\nGenerating summary report...")
    
    report_file = f"{output_prefix}_analysis_summary.txt"
    
    with open(report_file, 'w') as f:
        f.write("ORTHOFINDER SPECIES PHYLOGENETIC TREE ANALYSIS SUMMARY\n")
        f.write("=" * 65 + "\n\n")
        
        # Analysis parameters
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input Directory: {args.protein_dir}\n")
        f.write(f"Output Directory: {args.output}\n")
        f.write(f"Threads Used: {args.threads}\n")
        f.write(f"Tree Method: {args.tree_method if args.tree_method else 'Default'}\n")
        f.write(f"Results Directory: {results_dir.name}\n\n")
        
        # Species information
        if 'species_stats' in stats:
            f.write("SPECIES STATISTICS:\n")
            f.write("-" * 20 + "\n")
            df = stats['species_stats']
            for _, row in df.iterrows():
                f.write(f"{row.iloc[0]}: {row.iloc[1]} genes in {row.iloc[2]} orthogroups\n")
            f.write("\n")
        
        # Orthogroup statistics
        if 'orthogroups' in stats:
            f.write("ORTHOGROUP STATISTICS:\n")
            f.write("-" * 25 + "\n")
            df = stats['orthogroups']
            for _, row in df.iterrows():
                f.write(f"{row.iloc[0]}: {row.iloc[1]}\n")
            f.write("\n")
        
        # Gene count summary
        if 'gene_counts' in stats:
            df = stats['gene_counts']
            f.write("GENE COUNT SUMMARY:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total Orthogroups: {len(df)}\n")
            f.write(f"Species Analyzed: {len(df.columns) - 1}\n")  # -1 for Total column
            f.write(f"Average Genes per Species: {df.iloc[:, :-1].sum().mean():.0f}\n\n")
        
        # Tree information
        f.write("SPECIES TREE:\n")
        f.write("-" * 15 + "\n")
        f.write(f"Tree File: {tree_file.name if tree_file else 'Not found'}\n")
        f.write("Tree can be viewed in tree visualization software (FigTree, iTOL, etc.)\n\n")
        
        # File outputs
        f.write("OUTPUT FILES:\n")
        f.write("-" * 15 + "\n")
        f.write(f"• {report_file} - This summary report\n")
        f.write(f"• {output_prefix}_species_tree.newick - Species tree in Newick format\n")
        f.write(f"• {output_prefix}_species_tree.png - Tree visualization (if ete3 available)\n")
        f.write(f"• {output_prefix}_species_tree.pdf - Tree visualization (if ete3 available)\n")
        f.write(f"• {output_prefix}_species_tree.svg - Tree visualization (if ete3 available)\n")
        f.write(f"• {results_dir}/ - Complete OrthoFinder results\n\n")
        
        f.write("For detailed orthogroup analysis, check the OrthoFinder results directory.\n")
    
    print(f"✓ Summary report saved: {report_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate species phylogenetic trees using OrthoFinder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 orthofinder_phylogeny.py protein_sequences/
  python3 orthofinder_phylogeny.py proteins/ --threads 8 --tree-method fasttree
  python3 orthofinder_phylogeny.py proteins/ --output my_phylogeny --tree-method iqtree

Requirements:
  - OrthoFinder installed and in PATH
  - Protein FASTA files (.fasta, .fa, .fas, .faa, .pep) in input directory
  - Each file represents one species/genome
  - At least 3 species required for phylogenetic analysis

Optional for tree visualization:
  - ete3: conda install ete3
  - biopython: pip install biopython
        """
    )
    
    parser.add_argument("protein_dir",
                       help="Directory containing protein FASTA files (one per species)")
    parser.add_argument("--output", "-o", default="orthofinder_phylogeny",
                       help="Output directory prefix (default: orthofinder_phylogeny)")
    parser.add_argument("--threads", "-t", type=int, default=4,
                       help="Number of threads to use (default: 4)")
    parser.add_argument("--tree-method", "-m", choices=['fasttree', 'iqtree', 'raxml'],
                       help="Tree inference method (default: OrthoFinder default)")
    parser.add_argument("--additional-args", "-a",
                       help="Additional arguments to pass to OrthoFinder")
    parser.add_argument("--keep-temp", action="store_true",
                       help="Keep temporary OrthoFinder files")
    
    args = parser.parse_args()
    
    print("ORTHOFINDER SPECIES PHYLOGENETIC TREE PIPELINE")
    print("=" * 55)
    
    # Check dependencies
    available_methods = check_dependencies()
    
    # Validate tree method
    if args.tree_method:
        method_commands = {
            'fasttree': 'fasttree',
            'iqtree': 'iqtree',
            'raxml': 'raxmlHPC'
        }
        
        if method_commands[args.tree_method] not in available_methods:
            print(f"Warning: {args.tree_method} not found, using OrthoFinder default")
            args.tree_method = None
    
    # Validate input
    protein_files, species_names = validate_protein_files(args.protein_dir)
    
    # Create output directory
    output_dir = Path(f"{args.output}_orthofinder_results")
    output_dir.mkdir(exist_ok=True)
    
    # Run OrthoFinder
    run_orthofinder(args.protein_dir, output_dir, args.threads, 
                    args.tree_method, args.additional_args)
    
    # Find results
    results_dir = find_orthofinder_results(output_dir)
    
    # Extract species tree
    tree_file, tree_string = extract_species_tree(results_dir)
    
    if tree_file and tree_string:
        # Analyze statistics
        stats = analyze_orthofinder_stats(results_dir)
        
        # Visualize tree
        visualize_species_tree(tree_file, tree_string, args.output)
        
        # Create summary
        create_summary_report(results_dir, stats, tree_file, args.output, args)
        
        print(f"\n✓ Analysis complete! Results saved with prefix '{args.output}'")
        print(f"✓ OrthoFinder results: {results_dir}")
        
    else:
        print("Error: Could not extract species tree from results")
        sys.exit(1)

if __name__ == "__main__":
    main()