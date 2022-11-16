
homedir=~
eval homedir=$homedir

cactus_dir=$homedir/workspace/cactus
test_ds_dir=$cactus_dir/test_datasets
app_not_dir=$test_ds_dir/application_note
CACTUS=$homedir/workspace/cactus/main.nf


####################################################################
## create references

cd $cactus_dir/references

rm -r work worm fly mouse human util .cache .nextflow*

nextflow run $cactus_dir/scripts/create_references/create_references.nf -profile singularity


####################################################################
## test datasets

# worm
cd $test_ds_dir/worm
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity -resume
nextflow run ${CACTUS} -params-file parameters/no_rtr.yml -profile singularity -resume
nextflow run ${CACTUS} -params-file parameters/enrich_only_genes_self.yml -profile singularity -resume
nextflow run jsalignon/cactus -r e3e546b7e0f937019a6d0041923109c2207693dc -params-file parameters/no_enrich_fdr.yml -profile singularity -resume

# fly
cd $test_ds_dir/fly
nextflow run ${CACTUS} -params-file parameters/no_gtr.yml -profile singularity -resume
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity -resume

# mouse
cd $test_ds_dir/mouse
nextflow run ${CACTUS} -params-file parameters/no_enrich.yml -profile singularity -resume
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity -resume

# human
cd $test_ds_dir/human
nextflow run ${CACTUS} -params-file parameters/no_enrich.yml -profile singularity -resume
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity -resume


####################################################################
## application notes

# worm
cd $app_not_dir/worm
nextflow run jsalignon/cactus -params-file parameters/base.yml -profile singularity -resume -r 49913b378da7386fada4f23bebfdc00eba404533
nextflow run jsalignon/cactus -params-file parameters/vary_rank.yml -profile singularity -resume -r 49913b378da7386fada4f23bebfdc00eba404533
nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -resume -r 49913b378da7386fada4f23bebfdc00eba404533

# human
cd $app_not_dir/human
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity -resume -bg > nf_log.txt
nextflow run jsalignon/cactus -params-file parameters/vary_rank.yml -profile singularity -resume -r 49913b378da7386fada4f23bebfdc00eba404533
nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -resume -r 49913b378da7386fada4f23bebfdc00eba404533

nextflow run ${CACTUS} -params-file parameters/vary_FDR.yml -profile singularity -resume -bg


nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -resume -r 49913b378da7386fada4f23bebfdc00eba404533 -bg > nf_log.txt



####################################################################
## run times

# # application notes: worm
# 
# Completed at: 02-Nov-2022 08:52:47
# Duration    : 15h 45m 29s
# CPU hours   : 200.9
# Succeeded   : 4'508

# Completed at: 05-Nov-2022 00:29:04
# Duration    : 2d 15h 36m 8s
# CPU hours   : 769.8
# Succeeded   : 3'626

# nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -resume -r 49913b378da7386fada4f23bebfdc00eba404533
# Completed at: 12-Nov-2022 09:50:12
# Duration    : 4h 2m 13s
# CPU hours   : 224.4 (78.4% cached)
# Succeeded   : 3'979
# Cached      : 306
