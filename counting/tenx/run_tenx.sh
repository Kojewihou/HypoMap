#!/bin/bash

set -e

slurm_account_name="LAM-SL3-CPU"
slurm_partition_name="sapphire"
#raw_data_path="/rds/project/rds-O11U8YqSuCk/ehv20/raw/merged_data" #Uncomment to link data

cd ./scripts

#exec >"../logs/$0.log" 2>&1

# Create soft links to smart-seq data files in ../data
if [ ! -z "$raw_data_path" ];then
	bash link_data.sh "$raw_data_path" "10x|tenx"
fi

# Convert SRA formatted file names to bcl2 format for cellranger
bash rename_10x_files.sh && echo "Files renamed"

# Create a list of samples to be processed
bash generate_sample_list.sh && echo "Temp files created"

# Count number of samples to processed
count=$(cat ../temp/sample_list.txt | wc -l)
echo "$count samples found!"

# Submit job array only if there samples are to be processed
if [ -n "$count" ];then
	sbatch \
		-A "$slurm_account_name" -p "$slurm_partition_name" \
		--array 1-$count \
		slurm-cellranger.sh
fi

exit 0

echo "$0 finished!"
exit 0
