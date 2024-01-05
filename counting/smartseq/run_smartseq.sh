#!/bin/bash

set -e

slurm_account_name="LAM-SL3-CPU"
slurm_partition_name="sapphire"
#raw_data_path="/rds/project/rds-O11U8YqSuCk/ehv20/raw/merged_data" #Uncomment to link data

cd ./scripts

#exec >"../logs/$0.log" 2>&1

# Create soft links to smart-seq data files in ../data
if [ ! -z "$raw_data_path" ];then
	bash link_data.sh "$raw_data_path" "smartseq|smart-seq"
fi

for refid in $(find ../data/ -maxdepth 1 ! -type d | grep -o "v[0-9]-[0-9]*" | sort | uniq); do
	echo $refid

	# Ignore dataset if already processed
	if [ -d ../results/$refid/Solo.out/Gene ]; then
		echo "$refid already processed - skipping"
		continue
	fi

	echo "Generating manifest for ${refid}'s associated files."
	bash generate_manifest-ss.sh $refid && echo "Manifest generated in ./temp/"

	echo "Submitting job for $refid"
	job_id=$(sbatch --parsable \
		-A "$slurm_account_name" -p "$slurm_partition_name" \
		slurm-starsolo-ss.sh $refid)
	sbatch --dependency=afterok:"$job_id" \
		-A "$slurm_account_name" -p "$slurm_partition_name" \
		slurm-bam_processing-ss.sh $refid
done

echo "$0 finished!"
exit 0
