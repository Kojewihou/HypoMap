#!/bin/bash

set -e

slurm_account_name="LAM-SL3-CPU"
slurm_partition_name="sapphire"
raw_data_path="/rds/project/rds-O11U8YqSuCk/ehv20/raw/merged_data" #Uncomment to link data

cd ./scripts

# Create soft links to smart-seq data files in ../data
if [ ! -z "$raw_data_path" ];then
	bash link_data.sh "$raw_data_path" "dropseq|drop-seq"
fi

echo -ne "\033[KCreating sample list\r"
bash generate_sample_list.sh
echo "Sample list created in temp folder"

echo -ne "\033[KCreating manifest\r"
bash generate_manifest.sh
echo "Manifest created in temp folder"

sample_list="../temp/sample_list.txt"

echo "Submitting jobs"
count=$(cat $sample_list | wc -l)
sbatch  -a 1-$count \
	-A "$slurm_account_name" \
	-p "$slurm_partition_name" \
	slurm-starsolo-ds.sh

exit 0
