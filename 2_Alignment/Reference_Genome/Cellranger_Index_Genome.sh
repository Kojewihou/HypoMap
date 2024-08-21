#!/bin/bash

set -euo pipefail

cellranger mkref \
 --genome GRCm39 \
 --fasta ./Mus_musculus.GRCm39.dna.toplevel.fa \
 --genes Mus_musculus.GRCm39.110.woGm28040.gtf