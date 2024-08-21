#!/bin/bash

set -euo pipefail

# Configuration variables
CELLRANGER_IMAGE="path_to_cellranger_image.sif"  # Set your Cellranger image path here
OUTDIR="results/"

# Cellranger arguments (should be set before running the script)
TRANSCRIPTOME="path_to_transcriptome_dir"
FASTQS="path_to_fastq_dir"
SAMPLE=""

# Create necessary directories
mkdir -p "$OUTDIR"

# Change to output directory
cd "$OUTDIR"

# Run Cellranger count with Singularity
singularity run "$CELLRANGER_IMAGE" \
 cellranger count \
  --transcriptome "$TRANSCRIPTOME" \
  --fastqs "$FASTQS" \
  --sample "$SAMPLE" \
  --id "$SAMPLE" \
  --include-introns true \
  --chemistry "auto"

echo "Cellranger count completed for sample $SAMPLE"

# Set include introns to True when processin sn-RNAseq datasets
# Manually set chemistry for multiomic datasets