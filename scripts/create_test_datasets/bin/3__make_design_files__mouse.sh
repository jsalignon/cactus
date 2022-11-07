#!/bin/bash

species="mouse"
n_reads_atac=$1
n_reads_mrna=$2
prepro_dir="preprocessing/${species}"
source $create_test_datasets_bin_dir/create_test_datasets_functions.sh


# run.yml
cat > ${species}/parameters/run.yml << EOL
res_dir                   : 'results/test_mouse'
species                    : 'mouse'
chromatin_state           : 'ENCFF809HLK'
split__threshold_type     : 'rank' 
split__threshold_values   : [ 200, 1000 ]
EOL

# ENCFF809HLK	mm10	ChromHMM 18 state model for kidney (postnatal 0 days), mesoderm,	excretory system,	mouse

# atac_fastq.tsv and mrna_fastq.tsv
make_samples_info_file ${prepro_dir}
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(mrna_|atac_|rep)/, "", sample_id); \
  gsub(/old_/, "Old", sample_id); \
  gsub(/young_/, "Yng", sample_id); \
  gsub(/kidney/, "Kid", sample_id); \
  gsub(/liver/, "Liv", sample_id); \
  if (NR == 1) sample_id = "sample_id"; \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
make_fastq_info_file $species $n_reads_atac $n_reads_mrna

# comparisons.tsv
cat > ${species}/design/comparisons.tsv <<EOL
YngKid OldKid
YngLiv OldLiv
YngKid YngLiv
OldKid OldLiv
EOL

# groups.tsv
cat > ${species}/design/groups.tsv << EOL
all YngKid_vs_OldKid YngLiv_vs_OldLiv YngKid_vs_YngLiv OldKid_vs_OldLiv
age YngKid_vs_OldKid YngLiv_vs_OldLiv 
tissue YngKid_vs_YngLiv OldKid_vs_OldLiv
EOL

# regions_to_remove.tsv
touch ${species}/design/regions_to_remove.tsv

# genes_to_remove.tsv
touch ${species}/design/genes_to_remove.tsv

# run__no_enrich.yml
yml_file="${species}/parameters/run__no_enrich.yml"
cp ${species}/parameters/run.yml $yml_file
sed -i 's/test_mouse/test_mouse__no_enrich/g' $yml_file
cat >> $yml_file << EOL
disable_all_enrichments   : true
EOL

replace_spaces_by_tabs_in_the_design_tsv_files $species

