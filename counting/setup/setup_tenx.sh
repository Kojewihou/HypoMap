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

get_repo_file "scripts/utils/link_data.sh" "scripts/"
get_repo_file "scripts/utils/generate_sample_list.sh" "scripts/"

get_repo_file "scripts/tenx/rename_10x_files.sh" "scripts/"
get_repo_file "scripts/tenx/extract_10x_results.sh" "scripts/"
get_repo_file "scripts/tenx/slurm-cellranger.sh" "scripts/"

get_repo_file "scripts/tenx/map.txt" "temp/"

get_repo_file "scripts/tenx/run_tenx.sh" "."



echo "$0 finished!"
exit 0
