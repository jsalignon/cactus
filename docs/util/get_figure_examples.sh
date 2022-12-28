
testing_dir_name=testing2


####### running the test datasets

homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/$testing_dir_name
cpu_nb=47 ; memory_size='300G'
figshare_version=v4
tools_manager=singularity

species=worm ; cd $test_dir/$tools_manager/$species
nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version --res_dir 'results/almost_full_test' --species $species --chromatin_state 'iHMM.M1K16.worm_L3' --split__threshold_type 'rank' --split__threshold_values [200] 
nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version --res_dir 'results/no_enrich__no_rtr' --species $species --chromatin_state 'iHMM.M1K16.worm_L3' --split__threshold_type 'rank' --split__threshold_values [200] --disable_all_enrichments --design__regions_to_remove 'design/regions_to_remove_empty.tsv'

species=fly ; cd $test_dir/$tools_manager/$species
nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version --res_dir 'results/no_enrich' --species $species --chromatin_state 'iHMM.M1K16.fly_L3' --split__threshold_type 'rank' --split__threshold_values [200] --disable_all_enrichments
nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version --res_dir 'results/no_enrich__no_gtr' --species $species --chromatin_state 'iHMM.M1K16.fly_L3' --split__threshold_type 'rank' --split__threshold_values [200] --disable_all_enrichments --design__genes_to_remove 'design/genes_to_remove_empty.tsv'



####### getting example figures

homedir=~
eval homedir=$homedir
cactus_dir="$homedir/workspace/cactus"

cd $cactus_dir

examples_dir_pdf="docs/examples/pdf"
examples_dir_png="docs/examples/png"
examples_dir_html="docs/examples/html"
examples_dir_xlsx="docs/examples/xlsx"
# rm $examples_dir_png/* ; rm $examples_dir_pdf/* ; rm $examples_dir_html/* ; rm $examples_dir_xlsx/*

test_dir_2="$cactus_dir/$testing_dir_name/singularity"
res_worm_dir="$test_dir_2/worm/results/"
res_fly_dir="$test_dir_2/fly/results/"

run_worm_dir="$res_worm_dir/almost_full_test"
run_worm_figure_dir="$run_worm_dir/Figures_Individual"
run_worm_tables_dir="$run_worm_dir/Tables_Merged"

run_worm_no_rtr_dir="$res_worm_dir/no_enrich__no_rtr/Figures_Individual"
run_fly_wt_gtr_dir="$res_fly_dir/no_enrich/Figures_Individual"
run_fly_no_gtr_dir="$res_fly_dir/no_enrich__no_gtr/Figures_Individual"


# Volcano plots
cp $run_worm_no_rtr_dir/2_Differential_Analysis/ATAC__volcano/hmg4_vs_spt16__ATAC_volcano.pdf $examples_dir_pdf/hmg4_vs_spt16__ATAC_volcano__no_rtr.pdf
find $run_worm_tables_dir -name "*.xlsx" -exec cp "{}" $examples_dir_xlsx \;
cp $run_fly_no_gtr_dir/2_Differential_Analysis/mRNA__volcano/b170_vs_n301b170__mRNA_volcano.pdf $examples_dir_pdf/mRNA_volcano__no_gtr.pdf
cp $run_fly_wt_gtr_dir/2_Differential_Analysis/mRNA__volcano/b170_vs_n301b170__mRNA_volcano.pdf $examples_dir_pdf/mRNA_volcano__with_gtr.pdf

# Other plots
find $run_worm_figure_dir -name "*ctl_1*.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "*hmg4_vs_ctl_*_volcano.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "*hmg4_vs_ctl_*_PCA_*.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "*hmg4_vs_ctl_*other_plots.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "spearman_correlation_heatmap_without_outliers_without_control_cor.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__peaks__grouped__annotation_barplot.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__peaks__grouped__average_profile.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__peaks__grouped__distance_to_TSS.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "pca_top5000_without_control_pca.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__multiQC.html" -exec cp "{}" $examples_dir_html \;
find $run_worm_figure_dir -name "mRNA__multiQC.html" -exec cp "{}" $examples_dir_html \;
find $run_worm_figure_dir -name "all__down__200__hmg4_vs_ctl__venn_up_or_down.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "all__200__hmg4_vs_ctl__venn_up_and_down.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__all__down__200__hmg4_vs_ctl__*.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__all__down__200__hmg4_vs_spt16__motifs__barplot.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__all__200__all__*.pdf" -exec cp "{}" $examples_dir_pdf \;



cd $examples_dir_pdf

# for FILE in $(ls ATAC__all__down__200__hmg4_vs_spt16__motifs__barplot.pdf) 
# for FILE in $(ls mRNA_volcano*.pdf) 
for FILE in $(ls ATAC__all__200__all__*.pdf) 
# for FILE in $(ls *.pdf) 
do
  echo $FILE
  file_name=$(basename $FILE .pdf)
  pdftoppm -png -rx 300 -ry 300 $FILE > ${file_name}.png
done

mv *.png ../png
cd ../../..



## multiple files

cd $examples_dir_pdf

# for FILE in $(ls *.pdf) 
for FILE in $(ls *other_plots*.pdf) 
do
  echo $FILE
  file_name=$(basename $FILE .pdf)
  pdftoppm -png -rx 300 -ry 300 $FILE ${file_name}
done

mv *.png ../png
cd ../../..


FILE="ctl_1__reads_coverage.pdf"
file_name=$(basename $FILE .pdf)
pdftoppm -png -rx 300 -ry 300 $FILE > ${file_name}.png

mv *.png ../png
