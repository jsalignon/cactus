
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

# running all tests from scratch
for species in worm fly human mouse
do
  cd $test_ds_dir/$species
  rm -r work results
  nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity --executor_local_cpus 47 --executor_local_memory '250G' --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200]
done


for species in worm fly human mouse
do
  cd $test_ds_dir/$species
  nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity --executor_local_cpus 47 --executor_local_memory '250G' --res_dir 'results/full_test'  --split__peak_assignment ['all','prom','distNC'] --split__threshold_values [200,1000]
done

cd $cactus_dir/testing
nextflow run jsalignon/cactus/scripts/download/download.nf  -profile singularity --references --test_datasets --species worm -r main -latest --references_dir .


# running one by one specific files

# worm
cd $test_ds_dir/worm
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity
nextflow run ${CACTUS} -params-file parameters/no_rtr.yml -profile singularity
nextflow run ${CACTUS} -params-file parameters/enrich_only_genes_self.yml -profile singularity
nextflow run jsalignon/cactus -r e3e546b7e0f937019a6d0041923109c2207693dc -params-file parameters/no_enrich_fdr.yml -profile singularity

# fly
cd $test_ds_dir/fly
nextflow run ${CACTUS} -params-file parameters/no_gtr.yml -profile singularity
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity

# mouse
cd $test_ds_dir/mouse
nextflow run ${CACTUS} -params-file parameters/no_enrich.yml -profile singularity
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity

# human
cd $test_ds_dir/human
nextflow run ${CACTUS} -params-file parameters/no_enrich.yml -profile singularity
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity


# evaluation of the time needed to run Cactus on a typical laptop setup
cd $test_ds_dir/worm
rm -r work results
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity --executor_local_cpus 8 --executor_local_memory '16G' --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200]

cd $test_ds_dir/fly
rm -r work results
nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity --executor_local_cpus 8 --executor_local_memory '16G' --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200]

nextflow run jsalignon/cactus -params-file parameters/full_test.yml -profile singularity --executor_local_cpus 8 --executor_local_memory '16G' --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200] -r ea3ae1a882732e5ecdcaf70faf911e8679c68ff0 -resume


## worm
# Completed at: 17-Nov-2022 15:29:02
# Duration    : 26m 35s
# CPU hours   : 2.9
# Succeeded   : 874
# Execution status: Succeeded



nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest --test_datasets --references --references_dir refs -profile conda --species worm 

# running all tests from scratch


homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/testing

cpu_nb=47
memory_size='250G'

cd $test_dir 
# rm -r *

# for tools_manager in singularity conda
for tools_manager in mamba
do 
# for species in worm
for species in worm fly human mouse
  do
    cd $test_dir
    mkdir -p $tools_manager/$species
    cd $test_dir/$tools_manager/$species
    nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest --test_datasets --references --references_dir refs -profile $tools_manager --species $species -resume
    nextflow run jsalignon/cactus -r main -latest -params-file parameters/full_test.yml --references_dir $test_dir/$tools_manager/$species/refs -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200] -resume
  done
done

species=worm
tools_manager=mamba



homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/testing

# species=worm
# tools_manager=conda

for species in worm
# for species in worm fly human mouse
do 
# for tools_manager in singularity docker conda mamba 
for tools_manager in docker
  do
    cd $test_dir
    mkdir -p $tools_manager/$species
    cd $test_dir/$tools_manager/$species
    nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest --test_datasets --references --references_dir refs -profile $tools_manager --species $species -resume
    nextflow run jsalignon/cactus -r main -latest -params-file parameters/full_test.yml --references_dir $test_dir/$tools_manager/$species/refs -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200] -resume
  done
done


# worm conda (already made containers)
# Duration    : 10m 27s
# CPU hours   : 4.1




nextflow run /home/jersal/workspace/cactus/scripts/download/download.nf --test_datasets --references --species fly --references_dir refs --cactus_dir /home/jersal/workspace/cactus 

rm -r data parameters design refs
nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest --test_datasets --references --species fly --references_dir refs -profile singularity 

nextflow run jsalignon/cactus/scripts/download/download.nf -c https://github.com/jsalignon/cactus/scripts/download/nextflow.config -r main -latest --test_datasets --references --species fly --references_dir refs -profile singularity 
workflow 


nextflow run jsalignon/cactus/scripts/download/download.nf  -profile singularity --test_datasets --species fly -r ea3ae1a882732e5ecdcaf70faf911e8679c68ff0 --references_dir test


####################################################################
## application notes

# worm
cd $app_not_dir/worm
nextflow run jsalignon/cactus -params-file parameters/base.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533
nextflow run jsalignon/cactus -params-file parameters/vary_rank.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533
nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533

# human
cd $app_not_dir/human
nextflow run ${CACTUS} -params-file parameters/vary_FDR_no_enrich.yml -profile singularity -bg 



rm -r work

nextflow run jsalignon/cactus -profile singularity -r 8e43dd8c7045bd8bc8f3f426e3902c2ee52f93f5 -resume --executor_local_cpus 47 --executor_local_memory '250G' --species 'human' --res_dir 'results/2022_11_25__no_enrich' --chromatin_state 'ENCFF321DGG' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3,3,5,10] --split__peak_assignment ['all','distNC','prom'] --disable_all_enrichments true
mv .nextflow.log results/2022_11_25__no_enrich/Run_Info/nf.log


nextflow run jsalignon/cactus -profile singularity -r 8e43dd8c7045bd8bc8f3f426e3902c2ee52f93f5 -resume --executor_local_cpus 47 --executor_local_memory '250G' --species 'human' --res_dir 'results/2022_11_25__enrich_self' --chromatin_state 'ENCFF321DGG' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3,3,5,10] --split__peak_assignment ['all','distNC','prom'] --do_genes_self_enrichment true --do_peaks_self_enrichment true --do_gene_set_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true
mv .nextflow.log results/2022_11_25__no_enrich/Run_Info/nf.log

nextflow run jsalignon/cactus -profile singularity -r 8e43dd8c7045bd8bc8f3f426e3902c2ee52f93f5 -resume --executor_local_cpus 47 --executor_local_memory '250G' --species 'human' --res_dir 'results/2022_11_25__enrich' --chromatin_state 'ENCFF321DGG' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all']
mv .nextflow.log results/2022_11_25__enrich/Run_Info/nf.log

nextflow run jsalignon/cactus -profile singularity -r 8e43dd8c7045bd8bc8f3f426e3902c2ee52f93f5 -resume --executor_local_cpus 47 --executor_local_memory '250G' --species 'human' --res_dir 'results/2022_11_25__enrich_hihmm_chrom_state_1' --chromatin_state 'iHMM.M1K16.human_GM' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all'] --do_gene_set_enrichment false  --do_genes_self_enrichment false --do_peaks_self_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true
mv .nextflow.log results/2022_11_25__enrich_hihmm_chrom_state_1/Run_Info/nf.log

nextflow run jsalignon/cactus -profile singularity -r 8e43dd8c7045bd8bc8f3f426e3902c2ee52f93f5 -resume --executor_local_cpus 47 --executor_local_memory '250G' --species 'human' --res_dir 'results/2022_11_25__enrich_hihmm_chrom_state_2' --chromatin_state 'iHMM.M1K16.human_H1' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all'] --do_gene_set_enrichment false --do_genes_self_enrichment false --do_peaks_self_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true
mv .nextflow.log results/2022_11_25__enrich_hihmm_chrom_state_2/Run_Info/nf.log



cd $app_not_dir/worm

nextflow run jsalignon/cactus -profile singularity -r 8e43dd8c7045bd8bc8f3f426e3902c2ee52f93f5 -resume --executor_local_cpus 47 --executor_local_memory '250G' --species 'worm' --res_dir 'results/2022_11_25__enrich_self' --chromatin_state 'iHMM.M1K16.worm_L3' --chip_ontology 'all' --split__threshold_type 'FDR' --split__threshold_values [1.3,3,5,10] --split__peak_assignment ['all','distNC','prom'] --do_genes_self_enrichment true --do_peaks_self_enrichment true --do_gene_set_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true

nextflow run jsalignon/cactus -profile singularity -r 8e43dd8c7045bd8bc8f3f426e3902c2ee52f93f5 -resume --executor_local_cpus 47 --executor_local_memory '250G' --species 'worm' --res_dir 'results/2022_11_25__enrich' --chromatin_state 'iHMM.M1K16.worm_L3' --chip_ontology 'all' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all']






nextflow run ${CACTUS} -params-file parameters/full_test.yml -profile singularity -bg > nf_log.txt
nextflow run ${CACTUS} -params-file parameters/vary_FDR.yml -profile singularity -bg
nextflow run jsalignon/cactus -params-file parameters/vary_rank.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533
nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533
nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533 -bg > nf_log.txt



nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -r 874b7cd0280cecdf6cfc87407a00526d823e695c -resume

nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533 -resume
nextflow run jsalignon/cactus -params-file parameters/base.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533 -resume

nextflow run ${CACTUS} -params-file parameters/vary_FDR.yml -profile singularity -resume




####################################################################
## testing running Cactus from scratch

homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus

cd $cactus_dir/testing
nextflow drop jsalignon/cactus

nextflow run jsalignon/cactus/scripts/download/download.nf  -profile singularity --references --test_datasets --species fly -r main -latest --references_dir .



nextflow run jsalignon/cactus/scripts/download/download.nf  -profile singularity --references --test_datasets --species worm -r main -latest --references_dir .

nextflow run jsalignon/cactus -profile singularity -params-file parameters/full_test.yml -r main -latest -resume --references_dir . --


####################################################################
## run times


## building references
# Completed at: 16-Nov-2022 18:32:19
# Duration    : 31m 16s
# CPU hours   : 13.0 (89.3% cached)
# Succeeded   : 19
# Cached      : 3'666


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

# nextflow run jsalignon/cactus -params-file parameters/vary_FDR.yml -profile singularity -r 49913b378da7386fada4f23bebfdc00eba404533
# Completed at: 12-Nov-2022 09:50:12
# Duration    : 4h 2m 13s
# CPU hours   : 224.4 (78.4% cached)
# Succeeded   : 3'979
# Cached      : 306
