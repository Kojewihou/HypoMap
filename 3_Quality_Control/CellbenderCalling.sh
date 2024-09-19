#!/bin/bash

set -euo pipefile

MTX_DIR=/path/to/mtx/directory
N_EXP_CELLS=10000

cellbender remove-background \
 --input "$MTX_DIR" \
 --output output.h5 \
 --expected-cells "$N_EXP_CELLS" \
 --model full \
 --fpr 0 \
 --cuda || exit 1

python -u mtx_to_h5ad.py \
 "$MTX_DIR" \
 -b output_cell_barcodes.csv \
 -O output.h5ad