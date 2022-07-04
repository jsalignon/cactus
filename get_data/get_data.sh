

cactus_dir=~/workspace/cactus
singularity_dir=~/workspace/singularity_containers

get_data_dir=$cactus_dir/software/get_data
data_dir=$cactus_dir/data

export NXF_SINGULARITY_CACHEDIR=${singularity_dir}

tree -d data_archive/worm/




##############################################
### Creating all directories
##############################################

cd $data_dir ; mkdir worm ; cd worm
mkdir -p blacklisted_regions bowtie2_indexes_conta chromatin_states encode_CHIP genome/annotation genome/sequence homer_genome motifs_PWMS org_db refGene_UCSC
cd $data_dir ; cp -r worm fly ; cp -r worm mouse ; cp -r worm human
tree -d


##############################################
### Genome Fasta and Annotation files
##############################################

cd $data_dir 

sed -i 's/\r//g' $get_data_dir/bin/get_fasta_and_gff.sh # this command is needed to replace windows carriage return to unix
source $get_data_dir/bin/get_fasta_and_gff.sh

get_fasta_and_gff 100 dna_sm caenorhabditis_elegans WBcel235 toplevel worm
get_fasta_and_gff 100 dna_sm drosophila_melanogaster BDGP6.28 toplevel fly
get_fasta_and_gff 100 dna_sm mus_musculus GRCm38 primary_assembly mouse
get_fasta_and_gff 100 dna_sm homo_sapiens GRCh38 primary_assembly human


## in R
cd ~/workspace/cactus/data
R
cur_seq_info = rtracklayer::SeqinfoForUCSCGenome('ce11')
cur_seq_info@seqnames %<>% gsub('chr', '', .)
saveRDS(cur_seq_info, '~/workspace/cactus/data/worm/genome/sequence/cur_seqinfo.rds')


