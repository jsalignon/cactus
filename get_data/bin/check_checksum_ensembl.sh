#!/bin/bash

function check_checksum_and_gunzip () {
	downloaded_file=$1
  downloaded_file_name=$2
	checksums_file=$3
	local_checksum=$(sum $downloaded_file | tr -s ' ' | awk '{$1=$1;print}')
	ensembl_checksum=$(cat $checksums_file | grep $downloaded_file_name | tr -s ' ' | cut -f1,2 -d ' ' -)
  
  checksum_report="checksum__${downloaded_file}.txt"
  echo "downloaded file: $downloaded_file, checksum: $local_checksum" > $checksum_report
  echo "online file: $downloaded_file_name, checksum: $ensembl_checksum" >> $checksum_report
  
	if [[ "$local_checksum" != "$ensembl_checksum" ]] ; then
		echo "checksum is wrong for file ${downloaded_file}"
  else
    gunzip $downloaded_file
	fi
}

