#!/bin/bash

# Check if the required arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path_to_csv> <upstream_distance> <downstream_distance>"
    exit 1
fi

CSV_FILE=$1  # Capture the first argument (path to CSV)
UPSTREAM_DISTANCE=$2  # Capture the second argument (upstream distance)
DOWNSTREAM_DISTANCE=$3  # Capture the third argument (downstream distance)

# Run the Python code with the CSV file and distance arguments
python << EOF
import pyensembl
import csv

# Load the pyensembl data for human (GRCh38, release 101)
data = pyensembl.EnsemblRelease(release=101, species="homo_sapiens")

# CSV file and distances provided as arguments
csv_file = "$CSV_FILE"
upstream_distance = int("$UPSTREAM_DISTANCE")
downstream_distance = int("$DOWNSTREAM_DISTANCE")

# Read the genes from the CSV file, skipping the first row (header)
genes_of_interest = []
with open(csv_file, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header row
    for row in reader:
        genes_of_interest.append(row[0].strip())  # Assuming each row has one gene name in the first column

# Function to get the canonical transcript and TSS of a gene
def get_tss(gene_name):
    try:
        # Get a list of genes by name
        genes = data.genes_by_name(gene_name)
        if genes:
            gene = genes[0]  # Take the first gene if there are multiple

            # Search for canonical transcript or the best-supported protein-coding transcript
            canonical_transcript = None
            for transcript in gene.transcripts:
                if transcript.biotype == 'protein_coding':
                    # Select transcript if it's the canonical or best supported one
                    if canonical_transcript is None or transcript.support_level == 1:
                        canonical_transcript = transcript
            
            if canonical_transcript:
                # Correct handling of TSS based on the strand of the canonical transcript
                if canonical_transcript.strand == '+':
                    tss = canonical_transcript.start  # TSS is at the start for positive strand
                else:
                    tss = canonical_transcript.end  # TSS is at the end for negative strand

                # Calculate upstream and downstream based on user input
                if canonical_transcript.strand == '+':
                    upstream = tss - upstream_distance
                    downstream = tss + downstream_distance
                else:
                    upstream = tss + upstream_distance
                    downstream = tss - downstream_distance
                
                return gene_name, gene.contig, tss, canonical_transcript.strand, upstream, downstream
            else:
                return gene_name, None, None, None, None, None
        else:
            return gene_name, None, None, None, None, None
    except ValueError:
        return gene_name, None, None, None, None, None

# Collect TSS locations for all genes in the list
tss_locations = [get_tss(gene) for gene in genes_of_interest]

# Display results with specified upstream and downstream distances
print(f"Gene Name\tChromosome\tTSS\tStrand\tUpstream (-{upstream_distance}bp)\tDownstream (+{downstream_distance}bp)")
for gene_name, contig, tss, strand, upstream, downstream in tss_locations:
    if contig:
        strand_symbol = '+' if strand == '+' else '-'
        print(f"{gene_name}\t{contig}\t{tss}\t{strand_symbol}\t{upstream}\t{downstream}")
    else:
        print(f"{gene_name} not found in the Ensembl database")
EOF
