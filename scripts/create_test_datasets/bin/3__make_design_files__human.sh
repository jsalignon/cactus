#!/bin/bash

species="human"
n_reads_atac=$1
n_reads_mrna=$2
prepro_dir="preprocessing/${species}"
source $create_test_datasets_bin_dir/create_test_datasets_functions.sh


# full_test.yml
cat > ${species}/parameters/full_test.yml << EOL
res_dir                   : 'results/full_test'
species                    : 'human'
chromatin_state           : 'ENCFF941SVR'
chip_ontology             : 'cell_type.fibroblast'
split__threshold_type     : 'rank' 
split__threshold_values   : [ 200, 1000 ]
EOL

## details on the cell line:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98758
 # we measured the chromatin accessibility landscape using ATAC-seq following mock treatment, SSRP1 knockdown, or SUPT16H knockdown in human fibroblasts 
# cell type: human secondary fibroblasts
# genotype/variation: hiF-T cells carrying DOX-inducible, polycistronic human OCT4/KLF4/c-MYC/SOX2 (OKMS) cassette
# passages/stage: 13-18
# hiF-T cells -> derive from hBJ fibroblasts (= cell line established from skin taken from the normal foreskin of a neonatal male) https://www.sciencedirect.com/science/article/pii/S009286741500700X

## details on the chromatin state file:
# ENCFF941SVR: ChromHMM 18-state model of BSS00066: AG09309 from donor(s) ENCDO002AAA, cell line,	fibroblast,	skin of body, connective tissue

## alternative parameters if the cell are considered more like stem cell:
# chip_ontology = 'cell_type.stem_cell'
# chromatin_state = 'ENCFF676VUR'
# details on the chromatin state file: ENCFF676VUR, ChromHMM 18-state model of BSS00735: iPS-11a male adult (36 years) from donor(s) ENCDO632AGT,iPS-11a	cell line,	stem cell, induced pluripotent stem cell, skin of body

# atac_fastq.tsv and mrna_fastq.tsv
make_samples_info_file ${prepro_dir}
awk 'BEGIN {OFS = "\t"} { \
  condition = tolower($6) ; \
  replicate = $6 ; \
  gsub(/.*rep/, "", replicate); \
  gsub(/.*mock.*/, "ctl", condition); \
  gsub(/.*ssrp1.*/, "ssrp1", condition); \
  gsub(/.*supt16h.*/, "supt16h", condition); \
  sample_id = condition "_" replicate ; \
  if (NR == 1) sample_id = "sample_id"; \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
make_fastq_info_file $species $n_reads_atac $n_reads_mrna

# comparisons.tsv
cat > ${species}/design/comparisons.tsv <<EOL
ssrp1 ctl
supt16h ctl
ssrp1 supt16h
EOL

# groups.tsv
cat > ${species}/design/groups.tsv << EOL
all ssrp1_vs_ctl supt16h_vs_ctl ssrp1_vs_supt16h
ctl ssrp1_vs_ctl supt16h_vs_ctl
supt16h supt16h_vs_ctl ssrp1_vs_supt16h
EOL

# regions_to_remove.tsv
cat > ${species}/design/regions_to_remove.tsv << EOL
ssrp1 ssrp1->chr11:57,325,986-57,335,892
supt16h supt16h->chr14:21,351,476-21,384,019
EOL

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



