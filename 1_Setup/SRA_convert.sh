#!/bin/bash

set -euo pipefail

# Configuration variables
SRA_TOOLS_MODULE="sra-tools/"
TEMP_DIR="/tmp"

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <accession> <output_directory>"
    exit 1
fi

ACCESSION=$1
OUTPUT_DIR=$2

# Load SRA tools module
module load "$SRA_TOOLS_MODULE"


# Run fasterq-dump and compress the output FASTQ files
fasterq-dump "$ACCESSION" -p --split-files --include-technical --temp "$TEMP_DIR" -O "$OUTDIR"
find "$OUTDIR" -type f -iname "${ACCESSION}*.fastq" -exec pigz -v {} +
