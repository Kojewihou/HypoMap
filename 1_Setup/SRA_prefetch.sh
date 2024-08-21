#!/bin/bash

set -euo pipefail

# Configuration variables
SRA_TOOLS_MODULE="sra-tools/"

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <accession> <output_directory>"
    exit 1
fi

ACCESSION=$1
OUTPUT_DIR=$2

# Load SRA tools module
module load $SRA_TOOLS_MODULE

# Create necessary directory
mkdir -p "$OUTPUT_DIR"

# Fetch SRA file if not already present
if [ ! -f "$OUTPUT_DIR/${ACCESSION}/${ACCESSION}.sra" ]; then
    prefetch "$ACCESSION" -L 0 -X u -O "$OUTPUT_DIR" -r yes
else
    echo "$ACCESSION already prefetched from SRA"
fi
