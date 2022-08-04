

cactus_dir=~/workspace/cactus
get_test_datasets_bin_dir=$cactus_dir/software/get_test_datasets/bin
test_datasets_dir=$cactus_dir/test_datasets
cd $test_datasets_dir


##############################################
### Initialization
##############################################

source $get_test_datasets_bin_dir/initialization.sh

cd /home/jersal/workspace/cactus/data

# get_transcriptome_size_in_Mb (){ cat $1 | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM / 1000000}' }

get_transcriptome_size_in_Mb (){ cat $1/genome/annotation/bed_regions/exons.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM / 1000000}' ; }

get_transcriptome_size_in_Mb human
get_transcriptome_size_in_Mb mouse
get_transcriptome_size_in_Mb worm
get_transcriptome_size_in_Mb fly


genome_size=$(cut -f2 fly/genome/annotation/filtered/chromosomes_sizes.txt | paste -sd+ | bc)
echo "scale=2 ; $genome_size / 10^6" | bc

genome_size=$(cut -f2 human/genome/annotation/filtered/chromosomes_sizes.txt | paste -sd+ | bc)
echo "scale=2 ; $genome_size / 10^6" | bc

 grep -v ">" human/genome/sequence/genome.fa | wc | awk '{print $3-$1}' 

 grep -v ">" human/genome/sequence/genome.fa | wc | awk '{print $3-$1}' 


grep -v ">" human/genome/sequence/genome.fa | wc | head

grep -v ">" human/genome/sequence/genome.fa | head | wc | awk '{print ($3-$1) / 10^6}'


##############################################
### Worm (GSE98758)
##############################################

source $get_test_datasets_bin_dir/worm__get_fastq.sh

source $get_test_datasets_bin_dir/worm__subsample_reads.sh 300 100


##############################################
### Fly (GSE149339)
##############################################

source $get_test_datasets_bin_dir/worm__get_fastq.sh

source $get_test_datasets_bin_dir/worm__subsample_reads.sh 500 100


##############################################
### mouse (GSE181797)
##############################################

source $get_test_datasets_bin_dir/worm__get_fastq.sh

source $get_test_datasets_bin_dir/worm__subsample_reads.sh 4000 100


##############################################
### Human (GSE98758)
##############################################

source $get_test_datasets_bin_dir/human__get_fastq.sh

source $get_test_datasets_bin_dir/human__subsample_reads.sh 2000 100



##############################################
### The end
##############################################

# saving compressed objects and their sizes
du -h -d1 > test_datasets_sizes.txt
# tar --use-compress-program="pigz -p 15 -k -r" -cf worm.tar.gz worm


