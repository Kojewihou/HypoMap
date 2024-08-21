#!/bin/bash

set -euo pipefail

STAR \
 --runThreadN 56 \
 --runMode genomeGenerate \
 --genomeDir GRCm39/star-v2.7.10a/ \
 --genomeFastaFiles Mus_musculus.GRCm39.dna.toplevel.fa \
 --sjdbGTFfile Mus_musculus.GRCm39.110.woGm28040.gtf