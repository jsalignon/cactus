

cactus_dir=~/workspace/cactus
get_test_datasets_bin_dir=$cactus_dir/software/get_test_datasets/bin
test_datasets_dir=$cactus_dir/test_datasets
cd $test_datasets_dir


##############################################
### Initialization
##############################################

# clean install: rm -r preprocessing worm fly mouse human work

source $get_test_datasets_bin_dir/0__initialization.sh
# note: one needs to press enter for the script to finish

## Genome and Transcriptome sizes (can help to (still empirically) determine how many reads should be sampled by specie)
# specie  genome   transcriptome
# worm    100.286  53.1519
# fly     143.726  88.6633
# mouse   2730.87  157.86
# human   3099.75  261.067


##############################################
### Time estimates
##############################################

## Approximate time needed
# specie download sample_atac sample_mrna
#  worm   10m33s     1m33s       ?
#   fly   06m03s     1m20s       ?
#  mouse  33m38s    14m47s       ?
#  human  29m59s    15m45s       ?

## CPU hours
# specie download sample_atac sample_mrna
#  worm    4.6        0.4        ?
#   fly    1.2          ?        ?
#  mouse  33m38s    14m47s       ?
#  human  29m59s    15m45s       ?




##############################################
### Getting tests datasets
##############################################

# n_reads_atac=$1
# n_reads_mrna=$2


##############################################
### Worm (GSE98758)
##############################################

specie="worm"
n_reads_atac=1000
n_reads_mrna=50

source $get_test_datasets_bin_dir/1__get_fastq.sh $specie
rename -v 's/SRX/mrna_SRX/' ${fastq_dir}/SRX30291{12..20}*
rename -v 's/SRX/atac_SRX/' ${fastq_dir}/SRX30291{24..35}*
rename -v 's/SRX/atac_SRX/' ${fastq_dir}/SRX2333004*
ls "preprocessing/${specie}/fastq"

source $get_test_datasets_bin_dir/2__subsample_reads.sh $specie $n_reads_atac $n_reads_mrna
source $get_test_datasets_bin_dir/3__make_design_files__worm.sh $n_reads_atac $n_reads_mrna


##############################################
### Fly (GSE149339)
##############################################

specie="fly"
n_reads_atac=300
n_reads_mrna=100

source $get_test_datasets_bin_dir/1__get_fastq.sh $specie
rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX81740{44..53}*
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX81740{34..43}*
ls $fastq_dir

source $get_test_datasets_bin_dir/2__subsample_reads.sh $specie $n_reads_atac $n_reads_mrna
source $get_test_datasets_bin_dir/3__make_design_files__fly.sh $n_reads_atac $n_reads_mrna


##############################################
### mouse (GSE181797)
##############################################

specie="mouse"
n_reads_atac=6000
n_reads_mrna=150

source $get_test_datasets_bin_dir/1__get_fastq.sh $specie
rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX117086{63..78}*
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX117086{79..90}*
ls $fastq_dir

source $get_test_datasets_bin_dir/2__subsample_reads.sh $specie $n_reads_atac $n_reads_mrna
source $get_test_datasets_bin_dir/3__make_design_files__mouse.sh $n_reads_atac $n_reads_mrna


##############################################
### Human (GSE98758)
##############################################

specie="human"
n_reads_atac=6000
n_reads_mrna=250

source $get_test_datasets_bin_dir/1__get_fastq.sh $specie
rename -v 's/SRX/atac_SRX/' $fastq_dir/SRX2794*
rename -v 's/SRX/mrna_SRX/' $fastq_dir/SRX4029*
ls $fastq_dir

source $get_test_datasets_bin_dir/2__subsample_reads.sh $specie $n_reads_atac $n_reads_mrna
source $get_test_datasets_bin_dir/3__make_design_files__human.sh $n_reads_atac $n_reads_mrna



##############################################
### The end
##############################################

# saving compressed objects and their sizes
du -h -d1 > preprocessing/report/test_datasets_sizes.txt
# tar --use-compress-program="pigz -p 15 -k -r" -cf worm.tar.gz worm


