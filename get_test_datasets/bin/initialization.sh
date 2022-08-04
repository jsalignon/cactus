

##############################################
### Converting GSM ids to SRR ids
##############################################

# this was done since there was an issue with SRA that changed its data downloading procedures making the fetchngs pipeline not working anymore when downloading via SRA-tools. 
# a workaround that consist in using ENA instead of SRA is described here: https://github.com/nf-core/fetchngs/issues/98
# however this works only with SRA ids not GEO ids (because there is only a function to convert SRA its to ENA ids not from GEO ids in fetchngs) therefore we need to convert our GEO ids.
# this conversion was done using the entrez-direct API via a BioContainer


cd $singularity_dir
singularity pull https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1

cd $get_test_datasets_bin_dir/../samples_ids
source $get_test_datasets_bin_dir/get_test_datasets_functions.sh
gsm_to_srr worm
gsm_to_srr human
gsm_to_srr mouse
gsm_to_srr fly


##############################################
### Making the run configuration file
##############################################

cd $test_datasets_dir

cat > run.config << EOL

params {

  use_input_control = false
  
  save_bed_type = 'all'

  fdr_for_splitting_subsets = [ 0.2, 1.3 ] // setting -log10 FDR thresholds 

}

EOL


