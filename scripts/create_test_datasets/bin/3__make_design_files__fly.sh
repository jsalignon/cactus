#!/bin/bash

specie="fly"
n_reads_atac=$1
n_reads_mrna=$2
run_config_file=$specie/conf/run.config
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_bin_dir/get_test_datasets_functions.sh


# run.config
cat > ${specie}/conf/run.config <<EOL
params {
  specie            = 'fly'
  chromatin_state   = 'iHMM.M1K16.fly_L3'
  threshold_type_for_splitting_subsets   = 'rank' 
  threshold_values_for_splitting_subsets = [ 200, 1000 ]
}
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
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

# comparisons.tsv
cat > ${specie}/design/comparisons.tsv <<EOL
gaf ctl
b170 ctl
n301 ctl
n301b170 ctl
b170 n301b170
n301 n301b170
EOL

# groups.tsv
cat > ${specie}/design/groups.tsv << EOL
all gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
ctl gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl
n301b170 n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
EOL

# regions_to_remove.tsv
cat > ${specie}/design/regions_to_remove.tsv << EOL
gaf gaf->3L:14,747,929-14,761,049
b170 bap170->2R:6,636,512-6,642,358
n301 nurf301->3L:233,926-246,912
n301b170 bap170->2R:6,636,512-6,642,358
n301b170 nurf301->3L:233,926-246,912
EOL


# replacing spaces by tabs in all tsv files
tsv_files=$(ls ${specie}/design/*.tsv)
for tsv_file in ${tsv_files[@]}
do
  awk -i inplace -v OFS="\t" '$1=$1' $tsv_file
done


cat > ${specie}/yml/run.yml << EOL
specie                                 : 'worm'
use_input_control                      : false
save_bed_type                          : 'all'
chip_ontology                          : 'all'
chromatin_state                        : 'iHMM.M1K16.worm_L3'
threshold_type_for_splitting_subsets   : 'rank' 
threshold_values_for_splitting_subsets : [ 200, 1000 ]
design__mrna_fastq                     : 'design/mrna_fastq.tsv'
design__atac_fastq                     : 'design/atac_fastq.tsv'
design__comparisons                    : 'design/comparisons.tsv'
design__regions_to_remove              : 'design/regions_to_remove.tsv'
design__groups                         : 'design/groups.tsv'
EOL


replace_spaces_by_tabs_in_the_design_tsv_files $specie

