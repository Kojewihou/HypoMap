#!/bin/bash

mkdir -p data scripts temp logs results

repo="https://raw.githubusercontent.com/ehvr20/HypoMap/main/"

function get_repo_file () {
	filepath=$1
	output_dir=$2
	filename=$(basename $filepath)
	if [ ! -e "$output_dir/$filename" ];then
  	      wget -nv -P $output_dir "$repo/$filepath"
	else
	      echo "$output_dir/$filename already retrieved"
	fi
}

get_repo_file "counting/utils/matrix_h5_format.py" "scripts/"
get_repo_file "counting/utils/link_data.sh" "scripts/"
get_repo_file "counting/utils/extract_starsolo_results.sh" "scripts/"
get_repo_file "counting/smartseq/generate_manifest-ss.sh" "scripts/"
get_repo_file "counting/smartseq/slurm-starsolo-ss.sh" "scripts/"
get_repo_file "counting/smartseq/slurm-bam_processing-ss.sh" "scripts/"
get_repo_file "counting/smartseq/run_smartseq.sh" "."


echo "$0 finished!"
exit 0
