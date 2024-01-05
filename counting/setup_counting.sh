#!/bin/bash

mkdir -p tenx smartseq dropseq outs qc

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

get_repo_file "scripts/setup/setup_smartseq.sh" "smartseq/"
get_repo_file "scripts/setup/setup_dropseq.sh" "dropseq/"
get_repo_file "scripts/setup/setup_tenx.sh" "tenx/"

echo "$0 finished!"
exit 0
