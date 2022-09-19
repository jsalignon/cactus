#!/bin/bash

# cd $cactus
# docs/util/get_all_parameters.sh

# grep "**_params." docs/4_Prepro/ATAC_peaks.md

IFS=$'\n'
doc_files=($(find docs -name "*.md" -not -path "docs/3_Inputs/config.md"))
unset IFS
# echo ${doc_files[@]} | grep config

IFS=$'\n'
doc_files_sorted=($(sort <<<"${doc_files[*]}"))
unset IFS

all_config_file="docs/3_Inputs/all_config_entries.txt"

[[ -f $all_config_file ]] && rm $all_config_file
touch $all_config_file

for doc_file in ${doc_files_sorted[@]}
do
  nr=$(grep "**_params." $doc_file | wc -l)
  if [[ $nr -gt 0 ]]
  then
    echo -e "\n" >> $all_config_file
    header="${doc_file/docs\//}"
    header="${header/4_Prepro\//1. Preprocessing: }"
    header="${header/5_DA\//2. Differential Abundance: }"
    header="${header/6_Enrich\//3. Enrichment: }"
    header="${header/.md/}"
    echo "## $header" >> $all_config_file
    echo -e "" >> $all_config_file
    grep "**_params." $doc_file >> $all_config_file
  fi
done


