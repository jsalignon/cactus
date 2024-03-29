#!/bin/bash

species="fly"
n_reads_atac=$1
n_reads_mrna=$2
prepro_dir="preprocessing/${species}"


# full_test.yml
cat > ${species}/parameters/full_test.yml << EOL
res_dir                   : 'results/full_test'
species                    : 'fly'
chromatin_state           : 'iHMM.M1K16.fly_L3'
split__threshold_type     : 'rank' 
split__threshold_values   : [ 200, 1000 ]
EOL

# atac_fastq.tsv and mrna_fastq.tsv
make_samples_info_file ${prepro_dir}
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(atac|rna)seq_rep/, "", sample_id); \
  gsub(/bap170/, "b170", sample_id); \
  gsub(/nurf301/, "n301", sample_id); \
  gsub(/lacz/, "ctl", sample_id); \
  if (NR == 1) sample_id = "sample_id"; \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
make_fastq_info_file $species $n_reads_atac $n_reads_mrna

# comparisons.tsv
cat > ${species}/design/comparisons.tsv <<EOL
gaf ctl
b170 ctl
n301 ctl
n301b170 ctl
b170 n301b170
n301 n301b170
EOL

# groups.tsv
cat > ${species}/design/groups.tsv << EOL
all gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
ctl gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl
n301b170 n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
EOL

# regions_to_remove.tsv
cat > ${species}/design/regions_to_remove.tsv << EOL
gaf gaf->3L:14,747,929-14,761,049
b170 bap170->2R:6,636,512-6,642,358
n301 nurf301->3L:233,926-246,912
n301b170 bap170->2R:6,636,512-6,642,358
n301b170 nurf301->3L:233,926-246,912
EOL

# genes_to_remove.tsv
cat > ${species}/design/genes_to_remove.tsv << EOL
gaf Trl
b170 Bap170
n301 E(bx)
n301b170 Bap170
n301b170 E(bx)
EOL

# genes_to_remove_empty.tsv
touch ${species}/design/genes_to_remove_empty.tsv

# no_enrich.yml
yml_file="${species}/parameters/no_enrich.yml"
cp ${species}/parameters/full_test.yml $yml_file
sed -i 's/full_test/no_enrich/g' $yml_file
cat >> $yml_file << EOL
disable_all_enrichments   : true
EOL

# no_enrich__no_gtr.yml
yml_file="${species}/parameters/no_enrich__no_gtr.yml"
cp ${species}/parameters/no_enrich.yml $yml_file
sed -i 's/no_enrich/no_enrich__no_gtr/g' $yml_file
cat >> $yml_file << EOL
design__genes_to_remove   : 'design/genes_to_remove_empty.tsv'
EOL



replace_spaces_by_tabs_in_the_design_tsv_files $species

