

cactus_dir=~/workspace/cactus
create_test_ds_dir=$cactus_dir/scripts/create_test_datasets
create_test_datasets_bin_dir=$create_test_ds_dir/bin
samples_ids_dir=$create_test_ds_dir/samples_ids
test_datasets_dir=$cactus_dir/test_datasets
singularity_dir=/home/jersal/workspace/singularity_containers/
cd $test_datasets_dir

source $create_test_datasets_bin_dir/create_test_datasets_functions.sh


##############################################
### Initialization
##############################################

# clean install: rm -r preprocessing worm fly mouse human work

source $create_test_datasets_bin_dir/0__initialization.sh
# note: one needs to press enter for the script to finish

## Genome and Transcriptome sizes (can help to (still empirically) determine how many reads should be sampled by species)
# species  genome   transcriptome
# worm    100.286  53.1519
# fly     143.726  88.6633
# mouse   2730.87  157.86
# human   3099.75  261.067


##############################################
### Time estimates
##############################################

## Approximate time needed
# species download sample_atac sample_mrna
#  worm   10m33s     1m33s       ?
#   fly   06m03s     1m20s       ?
#  mouse  33m38s    14m47s       ?
#  human  29m59s    15m45s       ?

## CPU hours
# species download sample_atac sample_mrna
#  worm    4.6        0.4        ?
#   fly    1.2          ?        ?
#  mouse  33m38s    14m47s       ?
#  human   55.1    15m45s       ?




##############################################
### Getting tests datasets
##############################################

# n_reads_atac=$1
# n_reads_mrna=$2



##############################################
### Worm (GSE98758)
##############################################

species="worm" ; n_reads_atac=200 ; n_reads_mrna=50

gsm_to_srr $species $samples_ids_dir $singularity_dir

source $create_test_datasets_bin_dir/1__get_fastq.sh $species
rename -v 's/SRX/mrna_SRX/' ${fastq_dir}/SRX30291{12..20}*
rename -v 's/SRX/atac_SRX/' ${fastq_dir}/SRX30291{24..35}*
rename -v 's/SRX/atac_SRX/' ${fastq_dir}/SRX2333004*
ls "preprocessing/${species}/fastq"

source $create_test_datasets_bin_dir/2__subsample_reads.sh $species $n_reads_atac $n_reads_mrna
source $create_test_datasets_bin_dir/3__make_design_files__worm.sh $n_reads_atac $n_reads_mrna


##############################################
### Fly (GSE149339)
##############################################

species="fly" ; n_reads_atac=300 ; n_reads_mrna=100

gsm_to_srr $species $samples_ids_dir $singularity_dir

source $create_test_datasets_bin_dir/1__get_fastq.sh $species
rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX81740{44..53}*
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX81740{34..43}*
ls "preprocessing/${species}/fastq"

source $create_test_datasets_bin_dir/2__subsample_reads.sh $species $n_reads_atac $n_reads_mrna
source $create_test_datasets_bin_dir/3__make_design_files__fly.sh $n_reads_atac $n_reads_mrna


##############################################
### mouse (GSE193392)
##############################################

species="mouse" ; n_reads_atac=6000 ; n_reads_mrna=150

gsm_to_srr $species $samples_ids_dir $singularity_dir

source $create_test_datasets_bin_dir/1__get_fastq.sh $species
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX136541{74..81}*
rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX137050{91..98}*
ls "preprocessing/${species}/fastq"

source $create_test_datasets_bin_dir/2__subsample_reads.sh $species $n_reads_atac $n_reads_mrna
source $create_test_datasets_bin_dir/3__make_design_files__mouse.sh $n_reads_atac $n_reads_mrna


##############################################
### Human (GSE98758)
##############################################

species="human" ; n_reads_atac=15000 ; n_reads_mrna=250

gsm_to_srr $species $samples_ids_dir $singularity_dir

source $create_test_datasets_bin_dir/1__get_fastq.sh $species
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX2794*
rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX4029*
ls "preprocessing/${species}/fastq"

source $create_test_datasets_bin_dir/2__subsample_reads.sh $species $n_reads_atac $n_reads_mrna
source $create_test_datasets_bin_dir/3__make_design_files__human.sh $n_reads_atac $n_reads_mrna


##############################################
### Application note: Human and Worms (GSE98758)
##############################################

source $create_test_datasets_bin_dir/make_application_note_dirs.sh



##############################################
### The end
##############################################

# saving compressed objects and their sizes
du -h -d1 **/data/ > preprocessing/report/test_datasets_sizes.txt
du -h -d1 preprocessing/**/fastq >> preprocessing/report/test_datasets_sizes.txt
grep -v md5 preprocessing/report/test_datasets_sizes.txt



##############################################
### If needed: updating all design files only
##############################################

species="worm" ; n_reads_atac=200 ; n_reads_mrna=50
source $create_test_datasets_bin_dir/3__make_design_files__${species}.sh $n_reads_atac $n_reads_mrna

species="fly" ; n_reads_atac=300 ; n_reads_mrna=100
source $create_test_datasets_bin_dir/3__make_design_files__${species}.sh $n_reads_atac $n_reads_mrna

species="mouse" ; n_reads_atac=6000 ; n_reads_mrna=150
source $create_test_datasets_bin_dir/3__make_design_files__${species}.sh $n_reads_atac $n_reads_mrna

species="human" ; n_reads_atac=15000 ; n_reads_mrna=250
source $create_test_datasets_bin_dir/3__make_design_files__${species}.sh $n_reads_atac $n_reads_mrna




##############################################
### If needed: testing nf-core/fetchngs
##############################################

cd $samples_ids_dir
gsm_to_srr test

cd $test_datasets_dir
nextflow run nf-core/fetchngs --input "$samples_ids_dir/srr_accession/srr_test.txt" --outdir test -profile singularity -r 1.6 -resume



