
gsm_to_srr (){
  species=$1
  samples_ids_dir=$2
  singularity_dir=$3
  my_gsm_ids="$(cat $samples_ids_dir/gsm_accession/gsm_$species.txt | awk '{print}' ORS=' OR ' RS='\r\n' )"
  singularity exec $singularity_dir/entrez-direct:16.2--he881be0_1 esearch -db sra -query "$my_gsm_ids" | singularity exec $singularity_dir/entrez-direct:16.2--he881be0_1 efetch -format runinfo | cut -d ',' -f 1 - | grep -v 'Run' - > "$samples_ids_dir/srr_accession/srr_$species.txt"
}

# $1 = ${prepro_dir}   => head -3 ${prepro_dir}/samplesheet/samplesheet.csv 
make_samples_info_file (){
  samples_info_file=${1}/samplesheet/samples_info.tsv
  samplesheet=${1}/samplesheet/samplesheet.csv
  cat $samplesheet | sed 's/\",\"/\|/g' | cut -d"|" -f5,6,15,17,20,28 | sed 's/,\| /_/g' | sed 's/__/_/g' | awk 'BEGIN { FS = "|"; OFS = "\t"} {print $2, $1, $3, $4, $5, $6}'| column -t > $samples_info_file
  cat $samples_info_file
}

make_fastq_info_file (){
  species=$1
  n_reads_atac=$2
  n_reads_mrna=$3
  
  prepro_dir="preprocessing/${species}"
  fastq_info_file="${prepro_dir}/samplesheet/fastq_info.tsv"
  design_dir="${species}/design"
  atac_fastq_file="${design_dir}/atac_fastq.tsv"
  mrna_fastq_file="${design_dir}/mrna_fastq.tsv"
  
  awk 'BEGIN {OFS=""} { \
    library_strategy = tolower($5) ; \
    gsub(/-seq/, "", library_strategy) ; \
    gsub(/rna/, "mrna", library_strategy) ; \
    library_layout = $4
    gsub(/SINGLE/, "", library_layout) ; \
    gsub(/PAIRED/, "_R1", library_layout) ; \
    if (NR != 1) print $7, " data/", library_strategy, "/sample_100K_reads_", library_strategy, "_", $1, "_", $2, library_layout, ".fastq.gz" \
  }' ${prepro_dir}/samplesheet/samples_info_1.tsv > ${fastq_info_file}
  
  sed -i -e "s/_100K_reads_atac/_${n_reads_atac}K_reads_atac/g" ${fastq_info_file}
  sed -i -e "s/_100K_reads_mrna/_${n_reads_mrna}K_reads_mrna/g" ${fastq_info_file}
  
  grep atac ${fastq_info_file} > $atac_fastq_file
  grep mrna ${fastq_info_file} > $mrna_fastq_file
  
  echo $atac_fastq_file ; cat $atac_fastq_file
  echo $mrna_fastq_file ; cat $mrna_fastq_file

}

# note : in awk, $7 corresponds to the sample_id column that indicate the sample name used throughout the pipeline, and that is made manually


replace_spaces_by_tabs_in_the_design_tsv_files (){
  species=$1
  tsv_files=$(ls ${species}/design/*.tsv)
  for tsv_file in ${tsv_files[@]}
  do
    awk -i inplace -v OFS="\t" '$1=$1' $tsv_file
  done
}

