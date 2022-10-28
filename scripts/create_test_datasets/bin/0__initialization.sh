

##############################################
### Converting GSM ids to SRR ids
##############################################

# this was done since there was an issue with SRA that changed its data downloading procedures making the fetchngs pipeline not working anymore when downloading via SRA-tools. 
# a workaround that consist in using ENA instead of SRA is described here: https://github.com/nf-core/fetchngs/issues/98
# however this works only with SRA ids not GEO ids (because there is only a function to convert SRA its to ENA ids not from GEO ids in fetchngs) therefore we need to convert our GEO ids.
# this conversion was done using the entrez-direct API via a BioContainer


## commands to download the container if needed:
# cd $singularity_dir
# singularity pull https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1

cd $get_test_datasets_bin_dir/../samples_ids
source $get_test_datasets_bin_dir/get_test_datasets_functions.sh
gsm_to_srr worm
gsm_to_srr human
gsm_to_srr mouse
gsm_to_srr fly




##############################################
### Get genomes and transcriptome sizes
##############################################


get_genome_and_transcriptome_size () {

  get_fasta_size_in_Mb (){ grep -v ">" "/home/jersal/workspace/cactus/references/${1}/genome/sequence/${2}.fa" | wc | awk '{print ($3-$1) / 10^6}' ; }

  get_genome_and_transcriptome () { echo "$1 $(get_fasta_size_in_Mb $1 genome) $(get_fasta_size_in_Mb $1 transcriptome)" ; }

  get_genome_and_transcriptome_file () {
    cat | column -t > $1 << EOL
specie genome transcriptome
$(get_genome_and_transcriptome worm)
$(get_genome_and_transcriptome fly)
$(get_genome_and_transcriptome mouse)
$(get_genome_and_transcriptome human)
EOL
  }

  report_dir=preprocessing/report
  mkdir -p $report_dir
  report_file=$report_dir/genome_and_transcriptome_size.txt
  get_genome_and_transcriptome_file $report_file
  
  cat $report_file

}

get_genome_and_transcriptome_size




