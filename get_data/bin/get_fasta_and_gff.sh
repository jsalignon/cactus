#!/usr/bin/env bash

function get_fasta_and_gff() {
  
  release=$1
  dna_type=$2
  specie=$3
  genome=$4
  assembly=$5
  specie_short_name=$6
  
  out_folder=$specie_short_name/genome
  
  URL=ftp://ftp.ensembl.org/pub/release-$release/
  
  gff3_file=${specie^}.$genome.$release.gff3.gz
  fasta_file=${specie^}.$genome.$dna_type.$assembly.fa.gz
  
  wget -O annotation.gff3.gz $URL/gff3/$specie/$gff3_file 
  gunzip annotation.gff3.gz 
  mv annotation.gff3 $out_folder/annotation/annotation.gff3
  
  wget -O genome.fa.gz $URL/fasta/$specie/dna/$fasta_file 
  gunzip genome.fa.gz > $out_folder/sequence/genome.fa
  mv genome.fa $out_folder/sequence/genome.fa
  
  echo "gff3 file : $gff3_file" > $out_folder/README.txt
  echo "fasta file: $fasta_file" >> $out_folder/README.txt
  
}

