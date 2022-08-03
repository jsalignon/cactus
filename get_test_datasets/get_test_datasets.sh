

cactus_dir=~/workspace/cactus
singularity_dir=~/workspace/singularity_containers

get_test_datasets_dir=$cactus_dir/software/get_test_datasets
samples_ids_dir=$get_test_datasets_dir/samples_ids
test_datasets_dir=$cactus_dir/test_datasets

export NXF_SINGULARITY_CACHEDIR=${singularity_dir}


##############################################
### Converting GSM ids to SRR ids
##############################################

# this was done since there was an issue with SRA that changed its data downloading procedures making the fetchngs pipeline not working anymore when downloading via SRA-tools. 
# a workaround that consist in using ENA instead of SRA is described here: https://github.com/nf-core/fetchngs/issues/98
# however this works only with SRA ids not GEO ids (because there is only a function to convert SRA its to ENA ids not from GEO ids in fetchngs) therefore we need to convert our GEO ids.
# this conversion was done using the entrez-direct API via a BioContainer


cd $singularity_dir
singularity pull https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1

cd $samples_ids_dir
source $get_test_datasets_dir/get_test_datasets_functions.sh
gsm_to_srr worm
gsm_to_srr human
gsm_to_srr mouse
gsm_to_srr fly


## making the run configuration file

cd $test_datasets_dir

cat > run.config << EOL

params {

  use_input_control = false
  
  save_bed_type = 'all'

  fdr_for_splitting_subsets = [ 0.2, 1.3 ] // setting -log10 FDR thresholds 

}

EOL




cd $test_datasets_dir


##############################################
### Human (GSE98758)
##############################################

specie="human"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_dir/get_test_datasets_functions.sh

mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'human'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'ENCFF941SVR'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'cell_type.fibroblast'\n/" $specie/conf/run.config
cat $specie/conf/run.config

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
# sed -i "5s/^/\n  chip_ontology = 'cell_type.stem_cell'\n/" $specie/conf/run.config
# sed -i "5s/^/\n  chromatin_state = 'ENCFF676VUR'\n/" $specie/conf/run.config
# details on the chromatin state file: ENCFF676VUR, ChromHMM 18-state model of BSS00735: iPS-11a male adult (36 years) from donor(s) ENCDO632AGT,iPS-11a	cell line,	stem cell, induced pluripotent stem cell, skin of body


# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 
# --force_sratools_download # => using this options results in files of the format "SRR7101009_R1.fastq.gz" instead of "SRX2794538_SRR5521297_R1.fastq.gz" which crash my parsing script

# creating the sample_info file
make_samples_info_file ${prepro_dir}

# renaming files
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX2794*
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX4029*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*

# subsampling reads
n_reads_atac=2000
n_reads_mrna=100
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  condition = tolower($6) ; \
  replicate = $6 ; \
  gsub(/.*rep/, "", replicate); \
  gsub(/.*mock.*/, "ctl", condition); \
  gsub(/.*ssrp1.*/, "ssrp1", condition); \
  gsub(/.*supt16h.*/, "supt16h", condition); \
  sample_id = condition "_" replicate ; \
  gsub(/sample_title_sample_title/, "sample_id", sample_id)
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
cat ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
ssrp1 ctl
supt16h ctl
ssrp1 supt16h
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all ssrp1_vs_ctl supt16h_vs_ctl ssrp1_vs_supt16h
ctl ssrp1_vs_ctl supt16h_vs_ctl
supt16h supt16h_vs_ctl ssrp1_vs_supt16h
EOL

# regions to remove
cat > ${specie}/design/regions_to_remove.tsv << EOL
ssrp1 ssrp1->chr11:57,325,986-57,335,892
supt16h supt16h->chr14:21,351,476-21,384,019
EOL


##############################################
### Worm (GSE98758)
##############################################

specie="worm"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_dir/get_test_datasets_functions.sh

mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'worm'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'iHMM.M1K16.worm_L3'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'all'\n/" $specie/conf/run.config
cat $specie/conf/run.config

# downloading our fastq samples of interest and subsampling them
# nextflow run nf-core/fetchngs --input samples_id/sra_accession/sra_acc_worm.txt --outdir preprocessing/worm -profile singularity -r 1.6 -resume
# nextflow run nf-core/fetchngs --input samples_ids/gsm_accession/gsm_worm.txt --outdir preprocessing/worm -profile singularity -r 1.6  --force_sratools_download 
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX30291{12..20}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX30291{24..35}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX2333004*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*

# subsampling reads
n_reads_atac=200
n_reads_mrna=100
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(rep|-|_control)/, "", sample_id); \
  gsub(/rluc/, "ctl", sample_id); \
  gsub(/ctrl/, "ctl", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

cp ${specie}/design/atac_fastq.tsv ${specie}/design/atac_fastq__with_input.tsv
grep -v input ${specie}/design/atac_fastq__with_input.tsv > ${specie}/design/atac_fastq__without_input.tsv

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
hmg4 ctl
spt16 ctl
hmg4 spt16
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all hmg4_vs_ctl spt16_vs_ctl hmg4_vs_spt16
ctl hmg4_vs_ctl spt16_vs_ctl
spt16 spt16_vs_ctl hmg4_vs_spt16
EOL

# regions to remove
cat > ${specie}/design/regions_to_remove.tsv << EOL
hmg4 Hmg4->chrIII:7,379,143-7,381,596
spt16 Spt16->chrI:10,789,130-10,793,152
EOL


##############################################
### mouse (GSE181797)
##############################################

specie="mouse"
prepro_dir="preprocessing/${specie}"
source $get_test_datasets_dir/get_test_datasets_functions.sh

mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'mouse'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'ENCFF809HLK'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'all'\n/" $specie/conf/run.config
cat $specie/conf/run.config

# ENCFF809HLK	mm10	ChromHMM 18 state model for kidney (postnatal 0 days), mesoderm,	excretory system,	mouse

# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}

# renaming files
cat ${prepro_dir}/samplesheet/samples_info.csv
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX117086{63..78}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX117086{79..90}*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*
ls ${prepro_dir}/fastq

# subsampling reads
n_reads_atac=4000
n_reads_mrna=100
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(mrna_|atac_|rep)/, "", sample_id); \
  gsub(/old_/, "Old", sample_id); \
  gsub(/young_/, "Yng", sample_id); \
  gsub(/kidney/, "Kid", sample_id); \
  gsub(/liver/, "Liv", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv
cat ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
YngKid OldKid
YngLiv OldLiv
YngKid YngLiv
OldKid OldLiv
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all YngKid_vs_OldKid YngLiv_vs_OldLiv YngKid_vs_YngLiv OldKid_vs_OldLiv
age YngKid_vs_OldKid YngLiv_vs_OldLiv 
tissue YngKid_vs_YngLiv OldKid_vs_OldLiv
EOL

# regions to remove
touch ${specie}/design/regions_to_remove.tsv


##############################################
### Fly (GSE149339)
##############################################

specie="fly"
prepro_dir="preprocessing/${specie}"

mkdir -p $specie/data/mrna $specie/data/atac $specie/conf $specie/design

# creating the run.config file
cp run.config $specie/conf
sed -i "3s/^/\n  specie = 'fly'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chromatin_state = 'iHMM.M1K16.fly_L3'\n/" $specie/conf/run.config
sed -i "5s/^/\n  chip_ontology = 'all'\n/" $specie/conf/run.config
cat $specie/conf/run.config

# downloading our fastq samples of interest and subsampling them
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_${specie}.txt" --outdir ${prepro_dir} -profile singularity -r 1.6 -resume

# creating a simple reference file
make_samples_info_file ${prepro_dir}
cat ${prepro_dir}/samplesheet/samples_info.tsv

# renaming files
cat ${prepro_dir}/samplesheet/samples_info.csv
rename -v 's/SRX/mrna_SRX/' ${prepro_dir}/fastq/SRX81740{44..53}*
rename -v 's/SRX/atac_SRX/' ${prepro_dir}/fastq/SRX81740{34..43}*
rename -v 's/_1.fastq.gz/_R1.fastq.gz/' ${prepro_dir}/fastq/*
rename -v 's/_2.fastq.gz/_R2.fastq.gz/' ${prepro_dir}/fastq/*
ls ${prepro_dir}/fastq

# subsampling reads
n_reads_atac=200
n_reads_mrna=100
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_atac --experiment atac
nextflow $get_test_datasets_dir/subsample_fastq.nf --specie $specie --thousand_reads $n_reads_mrna --experiment mrna

# adding the sample_id column (to edit manually)
awk 'BEGIN {OFS = "\t"} { \
  sample_id = tolower($6) ; \
  gsub(/(atac|rna)seq_rep/, "", sample_id); \
  gsub(/bap170/, "b170", sample_id); \
  gsub(/nurf301/, "n301", sample_id); \
  gsub(/lacz/, "ctl", sample_id); \
  print $1, $2, $3, $4, $5, $6, sample_id \
}' ${prepro_dir}/samplesheet/samples_info.tsv | column -t > ${prepro_dir}/samplesheet/samples_info_1.tsv

# making the fastq design files
make_fastq_info_file $specie $n_reads_atac $n_reads_mrna

# making the comparison design file
cat > ${specie}/design/comparisons.tsv <<EOL
gaf ctl
b170 ctl
n301 ctl
n301b170 ctl
b170 n301b170
n301 n301b170
EOL

# making the groups design file
cat > ${specie}/design/groups.tsv << EOL
all gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
ctl gaf_vs_ctl b170_vs_ctl n301_vs_ctl n301b170_vs_ctl
n301b170 n301b170_vs_ctl b170_vs_n301b170 n301_vs_n301b170
EOL

# regions to remove
cat > ${specie}/design/regions_to_remove.tsv << EOL
gaf gaf->3L:14,747,929-14,761,049
b170 bap170->2R:6,636,512-6,642,358
n301 nurf301->3L:233,926-246,912
n301b170 bap170->2R:6,636,512-6,642,358
n301b170 nurf301->3L:233,926-246,912
EOL




##############################################
### The end
##############################################

# saving compressed objects and their sizes
du -h -d1 > test_datasets_sizes.txt
# tar --use-compress-program="pigz -p 15 -k -r" -cf worm.tar.gz worm


