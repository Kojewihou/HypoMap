#!/bin/bash

# Initialize a file to store all samples to be processed
>../temp/sample_list.txt
>../temp/all_sample_ids.txt

# Identify all sample names from files in fastqs folder
sample_ids=$(ls ../data/ | grep -o ".*_S[0-9]" | awk -F '_S[0-9]' '{print $1}' | sort | uniq)
echo -e "$sample_ids" > ../temp/all_sample_ids.txt
sample_ids=$(echo -e "$sample_ids" | cut -d'_' -f 1 | sort | uniq)

# Store samples to be processed only if not already processed i.e. not results directory found
for sample in $sample_ids; do
	if [[ ! -d "../results/$sample" ]];then
		#echo $sample
		echo $sample >> ../temp/sample_list.txt
	elif [[ -f "../results/__$sample.mro" ]];then
		# If job previously failed midway - pickup where it left off
		#echo $sample
		rm "../results/$sample/_lock" 2>/dev/null && echo "lock removed"; echo $sample >> ../temp/sample_list.txt
	else
		echo "$sample already processed"
	fi
done

exit 0