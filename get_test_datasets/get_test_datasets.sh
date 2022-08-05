

cactus_dir=~/workspace/cactus
get_test_datasets_bin_dir=$cactus_dir/software/get_test_datasets/bin
test_datasets_dir=$cactus_dir/test_datasets
cd $test_datasets_dir


##############################################
### Initialization
##############################################

source $get_test_datasets_bin_dir/initialization.sh

# specie  genome   transcriptome
# worm    100.286  53.1519
# fly     143.726  88.6633
# mouse   2730.87  157.86
# human   3099.75  261.067


##############################################
### Getting tests datasets
##############################################

# n_reads_atac=$1
# n_reads_mrna=$2


##############################################
### Worm (GSE98758)
##############################################

source $get_test_datasets_bin_dir/worm__get_fastq.sh
# Duration    : 2m 22s
# CPU hours   : 5.2 (100% cached)

source $get_test_datasets_bin_dir/worm__subsample_reads.sh 200 50
# atac:
# Duration    : 1m 20s
# CPU hours   : 0.3

##############################################
### Fly (GSE149339)
##############################################

source $get_test_datasets_bin_dir/fly__get_fastq.sh
# Duration    : 13m 28s
# CPU hours   : 2.3 (17.5% failed)

source $get_test_datasets_bin_dir/fly__subsample_reads.sh 300 100
# atac:

# mrna:


##############################################
### mouse (GSE181797)
##############################################

source $get_test_datasets_bin_dir/mouse__get_fastq.sh
# Duration    : 33m 38s
# CPU hours   : 14.4

source $get_test_datasets_bin_dir/mouse__subsample_reads.sh 6000 150
# atac:
# Duration    : 14m 47s
# CPU hours   : 5.4

##############################################
### Human (GSE98758)
##############################################

source $get_test_datasets_bin_dir/human__get_fastq.sh
# Duration    : 29m 59s
# CPU hours   : 5.8

source $get_test_datasets_bin_dir/human__subsample_reads.sh 6000 250
# atac:
# Duration    : 15m 45s
# CPU hours   : 2.5



##############################################
### The end
##############################################

# saving compressed objects and their sizes
du -h -d1 > preprocessing/report/test_datasets_sizes.txt
# tar --use-compress-program="pigz -p 15 -k -r" -cf worm.tar.gz worm


