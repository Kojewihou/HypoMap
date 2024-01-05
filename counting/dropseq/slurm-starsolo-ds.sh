#!/bin/bash

#SBATCH -D .
#SBATCH -o ../logs/dropseqSTARsolo-%A_%a.log
#SBATCH -c 28
#SBATCH -t 04:00:00
#SBATCH --mail-type=FAIL,REQUEUE

set -e

# Define paths
star_image_path=/rds/project/rds-O11U8YqSuCk/ehv20/singularity_images/STAR_image.sif
reference_genome_path=/rds/project/rds-O11U8YqSuCk/ehv20/reference_genome/GRCm39/star-v2.7.10a
gtf_file=/rds/project/rds-O11U8YqSuCk/ehv20/reference_genome/Mus_musculus.GRCm39.110.woGm28040.gtf
output_path=/rds/project/rds-O11U8YqSuCk/ehv20/counting/outs

# Define barcode and UMI lengths
CB_start=1
CB_len=12
UMI_start=13
UMI_len=8

# Extract sample ID
sample_id=$(head -n $SLURM_ARRAY_TASK_ID ../temp/sample_list.txt | tail -1)
run_name="$sample_id"

#TODO check if file already exists

# Define manifests files
mkdir -p ../temp/manifests
manifest_path=../temp/manifest.tsv
temp_manifest_path="../temp/manifests/$SLURM_ARRAY_TASK_ID.tsv"

# Generate temporary manifest file
grep "$sample_id" "$manifest_path" > $temp_manifest_path

# Create results directory
mkdir -p ../results/$run_name

# Run Drop-seq STARsolo pipeline
singularity run $star_image_path STAR \
    --runMode=alignReads \
    --soloType CB_UMI_Simple --soloStrand Forward \
    --soloCBstart "$CB_start" --soloCBlen "$CB_len" --soloUMIstart "$UMI_start" --soloUMIlen "$UMI_len" --soloCBwhitelist None \
    --soloBarcodeReadLength 0 \
    --soloCellFilter EmptyDrops_CR \
    --soloUMIdedup 1MM_CR \
    --soloMultiMappers EM \
    --quantMode GeneCounts \
    --runThreadN 56 \
    --limitOutSJcollapsed 20000000 \
    --outSAMtype BAM Unsorted \
    --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
    --outFileNamePrefix ../results/$run_name/ \
    --readFilesCommand zcat \
    --genomeDir "$reference_genome_path" \
    --sjdbGTFfile "$gtf_file" \
    --readFilesManifest $temp_manifest_path || exit 1

# Translate market matrix to h5
python matrix_h5_format.py \
    ../results/${run_name}/Solo.out/Gene/raw \
    ../results/${run_name}/Solo.out/Gene/Unique_matrix.h5 \
    -m "matrix.mtx"
python matrix_h5_format.py \
    ../results/${run_name}/Solo.out/Gene/raw \
    ../results/${run_name}/Solo.out/Gene/UniqueAndMult_matrix.h5 \
    -m "UniqueAndMult-EM.mtx"
python matrix_h5_format.py \
    ../results/${run_name}/Solo.out/Gene/filtered \
    ../results/${run_name}/Solo.out/Gene/Filtered_matrix.h5 \
    -m "matrix.mtx"

# Link up to count_matrices
ln -srvf ../results/${run_name}/Solo.out/Gene/Unique_matrix.h5 $output_path/${run_name}_unique_matrix.h5
ln -srvf ../results/${run_name}/Solo.out/Gene/UniqueAndMult_matrix.h5 $output_path/${run_name}_multimap_matrix.h5
ln -srvf ../results/${run_name}/Solo.out/Gene/Filtered_matrix.h5 $output_path/${run_name}_filtered_matrix.h5

# Check if sam file generated, if one exists, delete it
if [ -f ../results/$run_name/Aligned.out.sam ]; then
    rm ../results/$run_name/Aligned.out.sam
    echo "Removed SAM File"
fi

# Sort and compress bam output
samtools sort -@ 56 -o ../results/$run_name/Aligned.out.bam ../results/$run_name/Aligned.out.bam && \
mv ../results/$run_name/Aligned.out.bam ../results/$run_name/Sorted.out.bam
samtools view -C -@ 56 -o ../results/$run_name/Sorted_bam.cram ../results/$run_name/Sorted.out.bam && \
rm ../results/$run_name/Sorted.out.bam

mv ../logs/dropseqSTARsolo-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log ../logs/${run_name}.log

echo "$0 finished!"
exit 0
