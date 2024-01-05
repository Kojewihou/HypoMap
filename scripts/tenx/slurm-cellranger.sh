#!/bin/bash

#SBATCH -D .
#SBATCH -o ../logs/slurm-cellranger-%A_%a.log
#SBATCH -c 32
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL,REQUEUE

echo $0

manifest=~/genehunter/raw/manifest.csv
transcriptome=~/genehunter/reference_genome/GRCm39/
cellranger_image_path=~/genehunter/singularity_images/cellranger_image.sif

# Extract sample ids from sample_list.txt
sample_name=$(head -n $SLURM_ARRAY_TASK_ID ../temp/sample_list.txt | tail -n 1)

# Extract the reference id
refid=$(echo "$sample_name" | grep -o "v[0-9]-[0-9]*")

# Reference the  manifest to identify if cell or nuclear RNA sequencing was used
type=$(cat ${manifest} | grep "$refid" | cut -d',' -f 3)

# Find all valids sample_ids matching the sample name
valid_sample_ids=$(cat ../temp/all_sample_ids.txt | grep "$sample_name" | paste -sd,)

echo "Input Samples for ${sample_name}($type): ${valid_sample_ids}"
echo -e "${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}:\t${sample_name}\t<----\t${valid_sample_ids}" >> ../temp/sample_to_job_ids.txt

# Set introns flag dependant on experiment type
if [[ $type = "scRNAseq" ]]; then
    introns_option="--include-introns=false"
elif [[ $type = "snRNAseq" ]]; then
    introns_option="--include-introns=true"
else
    echo "Unknown type"
    exit 1
fi

# Set multiomic chemsitry otherwise all autodetection
if grep "$refid" "$manifest" | grep -i "multiome"; then
    set_chemistry="--chemistry=ARC-v1"
else
    set_chemistry="--chemistry=auto"
fi

# Run cellranger count on fastqs per sample
cd ../results
singularity run "${cellranger_image_path}" \
 cellranger count --transcriptome $transcriptome --fastqs ../data \
 --sample $valid_sample_ids --id $sample_name \
 ${introns_option} \
 ${set_chemistry}\
 --localmem 140

# Check the exit status and exit if it's not 0
if [ $? -ne 0 ]; then
    exit 1  # Exit the script
fi

# Compress bam output
samtools view -C -@ 56 -o ../results/$sample_name/outs/possorted_genome_bam.cram  ../results/$sample_name/outs/possorted_genome_bam.bam  && \
rm ../results/$sample_name/outs/possorted_genome_bam.bam

# Link up to count table merger
ln -srf "${sample_name}/outs/filtered_feature_bc_matrix.h5" "../../outs/${sample_name}_filtered_matrix.h5"
ln -srf "${sample_name}/outs/raw_feature_bc_matrix.h5" "../../outs/${sample_name}_raw_matrix.h5"

# Remove excess cellranger debug files
rm -r "${sample_name}/SC_RNA_COUNTER_CS"

exit 0
