#!/bin/bash

mkdir -p data scripts temp logs results

repo="https://raw.githubusercontent.com/ehvr20/HypoMap/main/"

function get_repo_file () {
	filepath=$1
	output_dir=$2
	filename=$(basename $filepath)
	if [ ! -e "$output_dir/$filename" ];then
  	      wget -P $output_dir "$repo/$filepath"
	else
	      echo "$output_dir/$filename already retrieved"
	fi
}

get_repo_file "scripts/utils/matrix_h5_format.py" "scripts/"
get_repo_file "scripts/utils/link_data.sh" "scripts/"
get_repo_file "scripts/utils/extract_starsolo_results.sh" "scripts/"
get_repo_file "scripts/smartseq/generate_manifest-ss.sh" "scripts/"
get_repo_file "scripts/smartseq/slurm-starsolo-ss.sh" "scripts/"
get_repo_file "scripts/smartseq/slurm-bam_processing-ss.sh" "scripts/"
get_repo_file "scripts/smartseq/run_smartseq.sh" "."


echo "$0 finished!"
exit 0
