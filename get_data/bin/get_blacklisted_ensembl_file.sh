#!/usr/bin/env bash

function get_blacklisted_ensembl_file (){
  specie_short_name=$1
  SPECIE_CODE=$2
  NCBI_CODE=$3
  PREPRO_BL_DIR=$4
  BLACKLIST_FILE="${SPECIE_CODE}-blacklist.v2.bed.gz"
  MAPPING_FILE="${NCBI_CODE}_UCSC2ensembl.txt"
  MAPPING_PATH="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master"
  BLACKLIST_PATH="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists"
  wget -O blacklist_NCBI.bed.gz $BLACKLIST_PATH/$BLACKLIST_FILE
  gunzip blacklist_NCBI.bed.gz
  wget -O NCBI_to_Ensembl.txt $MAPPING_PATH/$MAPPING_FILE
  singularity run $singularity_dir/cvbio:3.0.0--hdfd78af_1 cvbio UpdateContigNames -i blacklist_NCBI.bed -o blacklist_Ensembl.bed -m NCBI_to_Ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
  mv blacklist_NCBI.bed $PREPRO_BL_DIR/${SPECIE_CODE}_blacklist_NCBI.bed
  mv NCBI_to_Ensembl.txt $PREPRO_BL_DIR/${SPECIE_CODE}_NCBI_to_Ensembl.txt
  mv blacklist_Ensembl.bed $specie_short_name/blacklisted_regions
}
