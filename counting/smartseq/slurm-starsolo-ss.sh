#!/bin/bash

#SBATCH -D .
#SBATCH -o ../logs/smartseqSTARsolo-%A.log
#SBATCH -c 28
#SBATCH -t 10:00:00
#SBATCH --mail-type=FAIL,REQUEUE

set -e

star_image_path=/rds/project/rds-O11U8YqSuCk/ehv20/singularity_images/STAR_image.sif
reference_genome_path=/rds/project/rds-O11U8YqSuCk/ehv20/reference_genome/GRCm39/star-v2.7.10a
gtf_file=/rds/project/rds-O11U8YqSuCk/ehv20/reference_genome/Mus_musculus.GRCm39.110.woGm28040.gtf
output_path=/rds/project/rds-O11U8YqSuCk/ehv20/counting/outs

refid=$1

echo $refid

manifest_file=$(realpath ../temp/${refid}_manifest.tsv)
ls -l $manifest_file

mkdir -p ../results/$refid
singularity run $star_image_path STAR \
--runMode=alignReads \
--runThreadN=56 \
--quantMode GeneCounts \
--soloType SmartSeq --soloStrand Unstranded --soloUMIdedup NoDedup \
--outFilterScoreMin 30 \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ../results/$refid/ \
--readFilesCommand zcat \
--genomeDir "$reference_genome_path" \
--sjdbGTFfile "$gtf_file" \
--readFilesManifest "$manifest_file" || exit 1

# Translate market matrix to h5
python matrix_h5_format.py \
../results/${refid}/Solo.out/Gene/raw \
../results/${refid}/Solo.out/Gene/Unique_matrix.h5 \
-m "matrix.mtx"
python matrix_h5_format.py \
../results/${refid}/Solo.out/Gene/filtered \
../results/${refid}/Solo.out/Gene/Filtered_matrix.h5 \
-m "matrix.mtx"

# Link up to count_matrices
ln -srvf ../results/${refid}/Solo.out/Gene/Unique_matrix.h5 $output_path/${refid}_unique_matrix.h5
ln -srvf ../results/${refid}/Solo.out/Gene/Filtered_matrix.h5 $output_path/${refid}_filtered_matrix.h5

echo "$0 finished!"
exit 0
