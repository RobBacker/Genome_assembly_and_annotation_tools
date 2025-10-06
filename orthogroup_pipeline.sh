#!/bin/bash


# This is not production code, but included for publication posterity.
# Orthogroup Analysis Pipeline
# This script processes orthogroup data to extract sequences, align them, and analyze with Shannon entropy

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <orthogroup_file> <fasta_directory> [output_directory]"
    echo "Example: $0 orthogroups.txt ./fasta_files/ ./output/"
    exit 1
fi

ORTHOGROUP_FILE="$1"
FASTA_DIR="$2"
OUTPUT_DIR="${3:-./orthogroup_analysis}"

# Create output directories
mkdir -p "$OUTPUT_DIR/raw_sequences"
mkdir -p "$OUTPUT_DIR/alignments"
mkdir -p "$OUTPUT_DIR/shannon_results"

# Check if required tools are available
command -v muscle >/dev/null 2>&1 || { echo "Error: muscle is required but not installed. Please install muscle for sequence alignment." >&2; exit 1; }
command -v python >/dev/null 2>&1 || { echo "Error: python is required but not installed." >&2; exit 1; }

if [ ! -f "shannonent.py" ]; then
    echo "Warning: shannonent.py not found in current directory. Please ensure it's available."
fi

echo "Starting orthogroup analysis pipeline..."
echo "Input file: $ORTHOGROUP_FILE"
echo "FASTA directory: $FASTA_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Function to extract sequence from FASTA file
extract_sequence() {
    local seq_id="$1"
    local fasta_file="$2"
    local output_file="$3"
    
    # Use awk to extract the sequence
    awk -v id="$seq_id" '
    BEGIN { found=0; printing=0 }
    /^>/ { 
        if ($0 ~ ">" id) {
            found=1
            printing=1
            print $0
        } else if (printing==1) {
            printing=0
        }
    }
    !/^>/ && printing==1 { print $0 }
    ' "$fasta_file" >> "$output_file"
}

# Function to find FASTA file for a given sequence ID
find_fasta_file() {
    local seq_id="$1"
    # Extract genome name from sequence ID (everything before #)
    local genome_name=$(echo "$seq_id" | cut -d'#' -f1)
    
    # Look for matching FASTA files
    for fasta_file in "$FASTA_DIR"/*.fasta "$FASTA_DIR"/*.fa "$FASTA_DIR"/*.fas; do
        if [ -f "$fasta_file" ]; then
            local basename=$(basename "$fasta_file")
            # Check if genome name is in the filename
            if [[ "$basename" == *"$genome_name"* ]]; then
                echo "$fasta_file"
                return 0
            fi
        fi
    done
    
    # If no match found, try a broader search
    for fasta_file in "$FASTA_DIR"/*.fasta "$FASTA_DIR"/*.fa "$FASTA_DIR"/*.fas; do
        if [ -f "$fasta_file" ]; then
            # Check if the sequence exists in this file
            if grep -q ">.*$seq_id" "$fasta_file" 2>/dev/null; then
                echo "$fasta_file"
                return 0
            fi
        fi
    done
    
    return 1
}

# First pass: collect all sequences for each orthogroup
declare -A orthogroup_sequences
declare -A orthogroup_counts

echo "First pass: collecting all sequences for each orthogroup..."

# Skip header and process each line to collect sequences
tail -n +2 "$ORTHOGROUP_FILE" | while IFS=
    
    # Only proceed with alignment if we have sequences
    if [ -s "$ortho_fasta" ] && [ $sequence_count -gt 1 ]; then
        echo "  Found $sequence_count sequences for $orthogroup"
        
        # Create alignment
        alignment_file="$OUTPUT_DIR/alignments/${orthogroup}_aligned.fasta"
        echo "  Creating alignment..."
        
        muscle -in "$ortho_fasta" -out "$alignment_file" -quiet 2>/dev/null
        
        if [ -f "$alignment_file" ] && [ -s "$alignment_file" ]; then
            echo "  Alignment created successfully"
            
            # Run Shannon entropy analysis if shannonent.py exists
            if [ -f "shannonent.py" ]; then
                echo "  Running Shannon entropy analysis..."
                shannon_output="$OUTPUT_DIR/shannon_results/${orthogroup}_shannon.txt"
                python shannonent.py "$alignment_file" > "$shannon_output" 2>&1
                
                if [ $? -eq 0 ]; then
                    echo "  Shannon entropy analysis completed"
                else
                    echo "  Warning: Shannon entropy analysis failed for $orthogroup"
                fi
            fi
        else
            echo "  Warning: Alignment failed for $orthogroup"
        fi
    else
        echo "  Warning: Insufficient sequences found for $orthogroup (found: $sequence_count)"
        # Remove empty files
        [ -f "$ortho_fasta" ] && [ ! -s "$ortho_fasta" ] && rm "$ortho_fasta"
    fi
    
    echo ""
done < "$unique_orthogroups"

# Clean up temporary files
rm -f "$temp_sequences" "$unique_orthogroups"

echo "Pipeline completed!"
echo ""
echo "Results:"
echo "- Raw sequences: $OUTPUT_DIR/raw_sequences/"
echo "- Alignments: $OUTPUT_DIR/alignments/"
echo "- Shannon analysis: $OUTPUT_DIR/shannon_results/"
echo ""
echo "Summary:"
echo "Raw sequence files: $(ls "$OUTPUT_DIR/raw_sequences/"*.fasta 2>/dev/null | wc -l)"
echo "Alignment files: $(ls "$OUTPUT_DIR/alignments/"*.fasta 2>/dev/null | wc -l)"
echo "Shannon result files: $(ls "$OUTPUT_DIR/shannon_results/"*.txt 2>/dev/null | wc -l)"\t' read -r line; do
    # Convert line to array
    IFS=
    
    # Only proceed with alignment if we have sequences
    if [ -s "$ortho_fasta" ] && [ $sequence_count -gt 1 ]; then
        echo "  Found $sequence_count sequences for $orthogroup"
        
        # Create alignment
        alignment_file="$OUTPUT_DIR/alignments/${orthogroup}_aligned.fasta"
        echo "  Creating alignment..."
        
        muscle -in "$ortho_fasta" -out "$alignment_file" -quiet 2>/dev/null
        
        if [ -f "$alignment_file" ] && [ -s "$alignment_file" ]; then
            echo "  Alignment created successfully"
            
            # Run Shannon entropy analysis if shannonent.py exists
            if [ -f "shannonent.py" ]; then
                echo "  Running Shannon entropy analysis..."
                shannon_output="$OUTPUT_DIR/shannon_results/${orthogroup}_shannon.txt"
                python shannonent.py "$alignment_file" > "$shannon_output" 2>&1
                
                if [ $? -eq 0 ]; then
                    echo "  Shannon entropy analysis completed"
                else
                    echo "  Warning: Shannon entropy analysis failed for $orthogroup"
                fi
            fi
        else
            echo "  Warning: Alignment failed for $orthogroup"
        fi
    else
        echo "  Warning: Insufficient sequences found for $orthogroup (found: $sequence_count)"
        # Remove empty files
        [ -f "$ortho_fasta" ] && [ ! -s "$ortho_fasta" ] && rm "$ortho_fasta"
    fi
    
    echo ""
done

echo "Pipeline completed!"
echo ""
echo "Results:"
echo "- Raw sequences: $OUTPUT_DIR/raw_sequences/"
echo "- Alignments: $OUTPUT_DIR/alignments/"
echo "- Shannon analysis: $OUTPUT_DIR/shannon_results/"
echo ""
echo "Summary:"
echo "Raw sequence files: $(ls "$OUTPUT_DIR/raw_sequences/"*.fasta 2>/dev/null | wc -l)"
echo "Alignment files: $(ls "$OUTPUT_DIR/alignments/"*.fasta 2>/dev/null | wc -l)"
echo "Shannon result files: $(ls "$OUTPUT_DIR/shannon_results/"*.txt 2>/dev/null | wc -l)"\t' read -ra fields <<< "$line"
    
    # Process each set of 3 columns (Gene, Orthogroup, Category)
    for ((i=0; i<${#fields[@]}; i+=3)); do
        if [ $((i+2)) -ge ${#fields[@]} ]; then
            break
        fi
        
        gene_id="${fields[i]}"
        gene_orthogroup="${fields[i+1]}"
        category="${fields[i+2]}"
        
        # Skip empty fields
        if [ -z "$gene_id" ] || [ -z "$gene_orthogroup" ] || [ -z "$category" ]; then
            continue
        fi
        
        # Only collect if category is "Core"
        if [ "$category" = "Core" ]; then
            echo "$gene_id|$gene_orthogroup" >> "$OUTPUT_DIR/temp_sequences.txt"
        fi
    done
done

# Process collected sequences by orthogroup
echo "Second pass: processing sequences by orthogroup..."

# Get unique orthogroups
sort "$OUTPUT_DIR/temp_sequences.txt" | cut -d'|' -f2 | sort | uniq > "$OUTPUT_DIR/unique_orthogroups.txt"

while read -r orthogroup; do
    echo "Processing orthogroup: $orthogroup"
    
    # Create output file for this orthogroup
    ortho_fasta="$OUTPUT_DIR/raw_sequences/${orthogroup}.fasta"
    > "$ortho_fasta"  # Clear the file
    
    sequence_count=0
    
    # Get all sequences for this orthogroup
    grep "|$orthogroup$" "$OUTPUT_DIR/temp_sequences.txt" | cut -d'|' -f1 | while read -r gene_id; do
        echo "  Extracting sequence: $gene_id"
        
        # Find the appropriate FASTA file
        fasta_file=$(find_fasta_file "$gene_id")
        if [ $? -eq 0 ] && [ -f "$fasta_file" ]; then
            extract_sequence "$gene_id" "$fasta_file" "$ortho_fasta"
            echo $((sequence_count + 1)) > "$OUTPUT_DIR/temp_count_${orthogroup}.txt"
        else
            echo "    Warning: Could not find FASTA file or sequence for $gene_id"
        fi
    done
    
    # Read final count
    if [ -f "$OUTPUT_DIR/temp_count_${orthogroup}.txt" ]; then
        sequence_count=$(cat "$OUTPUT_DIR/temp_count_${orthogroup}.txt")
    else
        sequence_count=0
    fi
    
    # Only proceed with alignment if we have sequences
    if [ -s "$ortho_fasta" ] && [ $sequence_count -gt 1 ]; then
        echo "  Found $sequence_count sequences for $orthogroup"
        
        # Create alignment
        alignment_file="$OUTPUT_DIR/alignments/${orthogroup}_aligned.fasta"
        echo "  Creating alignment..."
        
        muscle -in "$ortho_fasta" -out "$alignment_file" -quiet 2>/dev/null
        
        if [ -f "$alignment_file" ] && [ -s "$alignment_file" ]; then
            echo "  Alignment created successfully"
            
            # Run Shannon entropy analysis if shannonent.py exists
            if [ -f "shannonent.py" ]; then
                echo "  Running Shannon entropy analysis..."
                shannon_output="$OUTPUT_DIR/shannon_results/${orthogroup}_shannon.txt"
                python shannonent.py "$alignment_file" > "$shannon_output" 2>&1
                
                if [ $? -eq 0 ]; then
                    echo "  Shannon entropy analysis completed"
                else
                    echo "  Warning: Shannon entropy analysis failed for $orthogroup"
                fi
            fi
        else
            echo "  Warning: Alignment failed for $orthogroup"
        fi
    else
        echo "  Warning: Insufficient sequences found for $orthogroup (found: $sequence_count)"
        # Remove empty files
        [ -f "$ortho_fasta" ] && [ ! -s "$ortho_fasta" ] && rm "$ortho_fasta"
    fi
    
    echo ""
done

echo "Pipeline completed!"
echo ""
echo "Results:"
echo "- Raw sequences: $OUTPUT_DIR/raw_sequences/"
echo "- Alignments: $OUTPUT_DIR/alignments/"
echo "- Shannon analysis: $OUTPUT_DIR/shannon_results/"
echo ""
echo "Summary:"
echo "Raw sequence files: $(ls "$OUTPUT_DIR/raw_sequences/"*.fasta 2>/dev/null | wc -l)"
echo "Alignment files: $(ls "$OUTPUT_DIR/alignments/"*.fasta 2>/dev/null | wc -l)"
echo "Shannon result files: $(ls "$OUTPUT_DIR/shannon_results/"*.txt 2>/dev/null | wc -l)"