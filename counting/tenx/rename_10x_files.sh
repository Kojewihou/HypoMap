#!/bin/bash
set -e

function rename () {
	file_path=$1
	refid=$(echo $file_path | grep -o "v[0-9]-[0-9]*")
	refid_map=$(cat ../temp/map.txt | grep "$refid")
	id=$(basename $file_path | grep -o "_[1-3]")
	column_number=$(echo "$refid_map" | awk -v val="$id" '{ for(i=1; i<=NF; i++) if($i==val) print i }')
	tag=$(awk -v col="$column_number" '{print $col}' <<< "$(head -n 1 ../temp/map.txt)")
	new_name=$(echo "$file_path" | sed "s|$id|_S1_L001_${tag}_001|")
	mv -v $file_path $new_name
}

export -f rename

echo "Starting $0"
find ../data/ -regex '.*/*[S,E]RR[0-9]*_[0-9].fastq.gz' | xargs -P 20 -I {} bash -c "rename {}"

exit 0
