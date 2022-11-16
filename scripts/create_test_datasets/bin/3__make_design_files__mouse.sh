#!/bin/bash

species="mouse"
n_reads_atac=$1
n_reads_mrna=$2
prepro_dir="preprocessing/${species}"
source $create_test_datasets_bin_dir/create_test_datasets_functions.sh


# full_test.yml
cat > ${species}/parameters/full_test.yml << EOL
res_dir                   : 'results/full_test'
species                   : 'mouse'
chromatin_state           : 'ENCFF809HLK'
split__threshold_type     : 'rank' 
split__threshold_values   : [ 200, 1000 ]
EOL

# ENCFF809HLK	mm10	ChromHMM 18 state model for kidney (postnatal 0 days), mesoderm,	excretory system,	mouse

# atac_fastq.tsv and mrna_fastq.tsv
make_samples_info_file ${prepro_dir}

awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(rna|atac|rep)/, "", sample_id); \
  gsub(/-seq_/, "", sample_id); \
  gsub(/-\/-/, "", sample_id); \
  gsub(/_control/, "", sample_id); \
  gsub(/glcstarv/, "Starv", sample_id); \
  gsub(/phf20/, "Phf", sample_id); \
  gsub(/wt/, "Wt", sample_id); \
  gsub(/_Starv/, "Starv", sample_id); \
  if (NR == 1) sample_id = "sample_id"; \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
make_fastq_info_file $species $n_reads_atac $n_reads_mrna

# comparisons.tsv
cat > ${species}/design/comparisons.tsv <<EOL
WtStarv Wt
PhfStarv Phf
PhfStarv WtStarv
Phf20 Wt
EOL

# groups.tsv
cat > ${species}/design/groups.tsv << EOL
all WtStarv_vs_Wt PhfStarv_vs_Phf PhfStarv_vs_WtStarv Phf_vs_Wt
Starv WtStarv_vs_Wt PhfStarv_vs_Phf
Phf Phf_vs_Wt PhfStarv_vs_WtStarv
EOL

# regions_to_remove.tsv
touch ${species}/design/regions_to_remove.tsv

# genes_to_remove.tsv
touch ${species}/design/genes_to_remove.tsv

# no_enrich.yml
yml_file="${species}/parameters/no_enrich.yml"
cp ${species}/parameters/full_test.yml $yml_file
sed -i 's/full_test/no_enrich/g' $yml_file
cat >> $yml_file << EOL
disable_all_enrichments   : true
EOL

replace_spaces_by_tabs_in_the_design_tsv_files $species

