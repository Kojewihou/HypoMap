#!/bin/bash

#set -e
manifest_path="/rds/project/rds-O11U8YqSuCk/ehv20/raw/manifest.csv"

raw_reads_path=$(realpath $1)
technology_pattern=$2

# Purge old links from data file
find ../data/ -type l -o -xtype l | xargs -P 10 -I {} rm {} 2>/dev/null && echo "Purged links in ../data"

# Identify all refids associated with Drop-Seq technologies
refids=$(cat $manifest_path | grep -i -E "$technology_pattern" | cut -d',' -f 1)

# Get files asssociated with each refid
for refid in $refids;do
	file_count=0
	total_count=$(find "$raw_reads_path" -name "$refid*.fastq.gz" | wc -l)
	for file in $(find "$raw_reads_path" -name "$refid*.fastq.gz");do
		filename=$(basename $file)
		ln -sr "$file" "../data/$filename"
		echo -ne "\033[2K$refid:\t IN-PROGRESS\t(${file_count}/${total_count})\t../data$filename\r"
		((file_count ++))
	done
	echo -e "\033[2K$refid:\tDATASET LINKED\t(${file_count}/${total_count})\t"
done

echo "$0 finished!"
exit 0
