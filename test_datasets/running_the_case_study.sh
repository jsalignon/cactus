
# setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_ds_dir=$cactus_dir/test_datasets
case_study_dir=$test_ds_dir/application_note

# setting up variables
cactus_version=0.9.0
cpu_nb=47
memory_size='250G'


# running the worm case study
cd $case_study_dir/worm
# rm -r work results

nextflow run jsalignon/cactus -profile singularity -r $cactus_version -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'worm' --res_dir 'results/2022_11_25__enrich_self' --chromatin_state 'iHMM.M1K16.worm_L3' --chip_ontology 'all' --split__threshold_type 'FDR' --split__threshold_values [1.3,3,5,10] --split__peak_assignment ['all','distNC','prom'] --do_genes_self_enrichment true --do_peaks_self_enrichment true --do_gene_set_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true
mv .nextflow.log results/2022_11_25__enrich_self/Run_Info/nf.log

nextflow run jsalignon/cactus -profile singularity -r $cactus_version -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'worm' --res_dir 'results/2022_11_25__enrich' --chromatin_state 'iHMM.M1K16.worm_L3' --chip_ontology 'all' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all']
mv .nextflow.log results/2022_11_25__enrich/Run_Info/nf.log


# running the human case study
cd $case_study_dir/human
# rm -r work results

nextflow run jsalignon/cactus -profile singularity -r $cactus_version -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'human' --res_dir 'results/2022_11_25__enrich_self' --chromatin_state 'ENCFF321DGG' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3,3,5,10] --split__peak_assignment ['all','distNC','prom'] --do_genes_self_enrichment true --do_peaks_self_enrichment true --do_gene_set_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true
mv .nextflow.log results/2022_11_25__enrich_self/Run_Info/nf.log

nextflow run jsalignon/cactus -profile singularity -r $cactus_version -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'human' --res_dir 'results/2022_11_25__enrich' --chromatin_state 'ENCFF321DGG' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all']
mv .nextflow.log results/2022_11_25__enrich/Run_Info/nf.log

nextflow run jsalignon/cactus -profile singularity -r $cactus_version -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'human' --res_dir 'results/2022_11_25__enrich_hihmm_chrom_state_1' --chromatin_state 'iHMM.M1K16.human_GM' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all'] --do_gene_set_enrichment false  --do_genes_self_enrichment false --do_peaks_self_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true
mv .nextflow.log results/2022_11_25__enrich_hihmm_chrom_state_1/Run_Info/nf.log

nextflow run jsalignon/cactus -profile singularity -r $cactus_version -resume --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --species 'human' --res_dir 'results/2022_11_25__enrich_hihmm_chrom_state_2' --chromatin_state 'iHMM.M1K16.human_H1' --chip_ontology 'cell_type.fibroblast' --split__threshold_type 'FDR' --split__threshold_values [1.3] --split__peak_assignment ['all'] --do_gene_set_enrichment false --do_genes_self_enrichment false --do_peaks_self_enrichment false --do_motif_enrichment false --do_chip_enrichment false --do_chrom_state_enrichment true
mv .nextflow.log results/2022_11_25__enrich_hihmm_chrom_state_2/Run_Info/nf.log

