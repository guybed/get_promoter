# get_promoter
A tool that pulls a range of genomic coordinates surrounding a TSS for a given list of genes. Useful for bigwig analysis.

This script processes a list of gene names from a CSV file, retrieves the canonical transcription start sites (TSS) from the Ensembl database, and calculates user-specified upstream and downstream regions relative to the TSS. The script is flexible, allowing users to specify how many base pairs upstream and downstream they wish to include.

Features
Retrieves the canonical or best-supported protein-coding transcript for each gene.
Calculates TSS based on the gene's strand (positive or negative).
Allows user-defined upstream and downstream distances.
Outputs a table with the following columns:
Gene Name
Chromosome
TSS
Strand
Upstream (-N bp from TSS)
Downstream (+M bp from TSS)

Requirements
Python 3.x
PyEnsembl: Python package for accessing Ensembl gene annotations.
set up ensenbl data for human as follows
import pyensembl
data = pyensembl.EnsemblRelease(release=101, species="homo_sapiens")
data.download()
data.index()
