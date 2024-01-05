#!/bin/bash

#set -e

refid=$1

if [ -z "$refid" ];then
	echo "ERROR: No refid given"
	exit 1
fi


# Get all associated files
files=$(find ../data -name "$refid*.fastq.gz")

# Process the files using awk to extract cell IDs, read1, and read2
cell_ids=$(printf "%s\n" "$files" | awk -F_ '{print substr($0, 1, length($0)-length($NF)-1)}' | sort -u)
read1=$(printf "%s\n" "$files" | tr ' ' '\n' | grep '_\(1\|R1\)' | sort)
read2=$(printf "%s\n" "$files" | tr ' ' '\n' | grep '_\(2\|R2\)' | sort)

# Output all files into tsv manifest of 3 columns: 	read2	read1	sample_id
paste -d'\t' <(echo "${read2}") <(echo "${read1}") <(echo "$cell_ids") > ../temp/${refid}_manifest.tsv

exit 0
