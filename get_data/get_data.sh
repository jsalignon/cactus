

cactus_dir=~/workspace/cactus
singularity_dir=~/workspace/singularity_containers

get_data_dir=$cactus_dir/software/get_test_datasets
samples_ids_dir=$get_data_dir/samples_ids
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


function get_fasta_and_gff {
  
  release=$1
  dna_type=$2
  specie=$3
  genome=$4
  assembly=$5
  specie_short_name=$6
  
  URL=ftp://ftp.ensembl.org/pub/release-$release/
  
  gff3_file=${specie^}.$genome.$release.gff3.gz
  fasta_file=${specie^}.$genome.$dna_type.$assembly.fa.gz
  
  wget -O - $URL/gff3/$specie/$gff3_file | gunzip -c > $specie_short_name/annotation/anno.gff3
  wget -O - $URL/fasta/$specie/dna/$fasta_file | gunzip -c > $specie_short_name/sequence/genome.fa
  echo "gff3 file : $gff3_file" > $genome/README
  echo "fasta file: $fasta_file" >> $genome/README
  cd $GENF
  
}

get_fasta_and_gff 100 dna_sm caenorhabditis_elegans WBcel235 toplevel worm
get_fasta_and_gff 100 dna_sm drosophila_melanogaster BDGP6.28 toplevel fly
get_fasta_and_gff 100 dna_sm mus_musculus GRCm38 primary_assembly mice
get_fasta_and_gff 100 dna_sm homo_sapiens GRCh38 primary_assembly human


