#!/bin/bash

#SBATCH -D .
#SBATCH -o ../logs/STARsolo-ss-%A.log
#SBATCH -p sapphire
#SBATCH -c 28
#SBATCH -t 12:00:00
#SBATCH --mail-type=FAIL,REQUEUE

refid=$1

echo "Post-processing $refid associated bam files"

# Check if sam file generated, if one exists, delete it
if [ -f ../results/$refid/Aligned.out.sam ]; then
        rm ../results/$refid/Aligned.out.sam
        echo "Removed SAM File"
fi

# Sort outputted bam file
echo "Soring $refid associated bam files"
samtools sort -@ 56 ../results/$refid/Aligned.out.bam -o ../results/$refid/Aligned_Sorted.out.bam
rm ../results/$refid/Aligned.out.bam

# Compress sorted bam file
echo "Compressing sorted $refid associated bam files"
samtools view -b -C -@ 56 -o  ../results/$refid/possorted_genome_bam.cram ../results/$refid/Aligned_Sorted.out.bam
rm ../results/$refid/Aligned_Sorted.out.bam

echo "$0 finished successfully"
exit 0

