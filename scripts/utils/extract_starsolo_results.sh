#!/bin/bash

output_path=$(realpath ../../outs)

for unique_matrix in $(find ../results -type f -name "Unique_matrix.h5");do
	refid=$(echo $unique_matrix | cut -d'/' -f 3)
	echo $refid
	echo $unique_matrix
	ln -srvf $unique_matrix $output_path/${refid}_unique_matrix.h5
done

for multimap_matrix in $(find ../results -type f -name "UniqueAndMult_matrix.h5");do
        refid=$(echo $multimap_matrix | cut -d'/' -f 3)
        echo $refid
        echo $multimap_matrix
        ln -srvf $multimap_matrix $output_path/${refid}_multimap_matrix.h5
done

for filtered_matrix in $(find ../results -type f -name "Filtered_matrix.h5");do
        refid=$(echo $unique_matrix | cut -d'/' -f 3)
        echo $refid
        echo $filtered_matrix
        ln -srvf $filtered_matrix $output_path/${refid}_filtered_matrix.h5
done

echo "$0 finished!"
exit 0
