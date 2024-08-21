#!/bin/bash

set -euo pipefail

BAM='path_to_bam_file'
OUTDIR='results/ID'

cellranger bamtofastq \
 --nthreads=8 \
 $BAM \
 $OUTDIR