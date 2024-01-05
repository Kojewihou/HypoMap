#!/bin/bash

function extract_outs () {
	dir_path=$1
	dir_name=$(echo "$dir_path" | cut -d'/' -f3)

	# Check the directory is finished
	if [[ -f "../results/__$dir_name.mro" ]];then
		echo "$dir_name only partially resolved -> please ensure successfuly cellranger run"; return
	fi

	ln -srf "$dir_path/outs/filtered_feature_bc_matrix.h5" "../../outs/${dir_name}_filtered_matrix.h5"
	ln -srf "$dir_path/outs/raw_feature_bc_matrix.h5" "../../outs/${dir_name}_unique_matrix.h5"

	return
}

export -f extract_outs

find ../results -maxdepth 1 -type d | sort | tail -n +2 | xargs -P 6 -I {} bash -c "extract_outs {}"

# Remove deadlinks
find ../../outs -xtype l -exec rm {} \;

exit 0
