
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
# nextflow run  jsalignon/cactus -r hotfix/0.8.1  -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version -params-file parameters/full_test.yml --split__threshold_values [1000] --res_dir 'results/almost_full_test'
nextflow run  jsalignon/cactus -r v0.8.6  -profile $tools_manager \
  --executor_local_cpus $cpu_nb --executor_local_memory $memory_size \
  --references_dir $cactus_dir/references/$figshare_version \
  -params-file parameters/full_test.yml --split__threshold_values [1000] \
  --res_dir 'results/v0.8.6'
# nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version -params-file parameters/full_test.yml --split__threshold_values [1000] --res_dir 'results/almost_full_test'
# nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version -params-file parameters/no_enrich__no_rtr.yml

# species=fly ; cd $test_dir/$tools_manager/$species
# nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version -params-file parameters/no_enrich.yml
# nextflow run  jsalignon/cactus -r main -latest -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --references_dir $cactus_dir/references/$figshare_version -params-file parameters/no_enrich__no_gtr.yml



####### getting example figures

homedir=~
eval homedir=$homedir
cactus_dir="$homedir/workspace/cactus"

cd $cactus_dir

## Paths where examples are stored
examples_dir_pdf="$cactus_dir/docs/examples/pdf"
examples_dir_png="$cactus_dir/docs/examples/png"
examples_dir_html="$cactus_dir/docs/examples/html"
examples_dir_xlsx="$cactus_dir/docs/examples/xlsx"
# rm $examples_dir_png/* ; rm $examples_dir_pdf/* ; rm $examples_dir_html/* ; rm $examples_dir_xlsx/*

## Paths where examples are copied from
test_dir_2="$cactus_dir/$testing_dir_name/singularity"
res_worm_dir="$test_dir_2/worm/results/"
res_fly_dir="$test_dir_2/fly/results/"

# run_worm_dir="$res_worm_dir/almost_full_test"
run_worm_dir="$res_worm_dir/v0.8.6"
run_worm_figure_dir="$run_worm_dir/Figures_Individual"
run_worm_tables_dir="$run_worm_dir/Tables_Merged"

run_worm_no_rtr_dir="$res_worm_dir/no_enrich__no_rtr/Figures_Individual"
run_fly_wt_gtr_dir="$res_fly_dir/no_enrich/Figures_Individual"
run_fly_no_gtr_dir="$res_fly_dir/no_enrich__no_gtr/Figures_Individual"

# Volcano plots
cp $run_worm_no_rtr_dir/2_Differential_Analysis/ATAC__volcano/hmg4_vs_spt16__ATAC_volcano.pdf $examples_dir_pdf/hmg4_vs_spt16__ATAC_volcano__no_rtr.pdf
find $run_worm_figure_dir/2_Differential_Analysis -name "hmg4_vs_spt16__ATAC_volcano.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_tables_dir -name "*.xlsx" -exec cp "{}" $examples_dir_xlsx \;
cp $run_fly_no_gtr_dir/2_Differential_Analysis/mRNA__volcano/b170_vs_n301b170__mRNA_volcano.pdf $examples_dir_pdf/mRNA_volcano__no_gtr.pdf
cp $run_fly_wt_gtr_dir/2_Differential_Analysis/mRNA__volcano/b170_vs_n301b170__mRNA_volcano.pdf $examples_dir_pdf/mRNA_volcano__with_gtr.pdf

# Venn diagrams
find $run_worm_figure_dir -name "all__down__1000__hmg4_vs_ctl__venn_up_or_down.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "all__1000__hmg4_vs_ctl__venn_up_and_down.pdf" -exec cp "{}" $examples_dir_pdf \;

# Barplots
find $run_worm_figure_dir -name "ATAC__all__down__1000__hmg4_vs_ctl__*.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "ATAC__all__down__1000__hmg4_vs_spt16__motifs__barplot.pdf" -exec cp "{}" $examples_dir_pdf \;

# Heatmaps
find $run_worm_figure_dir -name "ATAC__all__1000__all__*.pdf" -exec cp "{}" $examples_dir_pdf \;

# Other plots
find $run_worm_figure_dir -name "*ctl_1*.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir/2_Differential_Analysis -name "*hmg4_vs_ctl*.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir/1_Preprocessing/ATAC__peaks__grouped_plots -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "spearman_correlation_heatmap_without_outliers_without_control_cor.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -name "pca_top5000_without_control_pca.pdf" -exec cp "{}" $examples_dir_pdf \;
find $run_worm_figure_dir -regex ".*\(ATAC\|mRNA\)__multiQC.html" -exec cp "{}" $examples_dir_html \;


find $run_worm_figure_dir -name "hmg4_*.bw"

find $run_worm_figure_dir -name "ATAC__all__1000__all__*.pdf" -exec cp "{}" $examples_dir_pdf \;


## Making png images from PDF
cd $examples_dir_pdf

# for FILE in $(ls ATAC__all__down__200__hmg4_vs_spt16__motifs__barplot.pdf) 
# for FILE in $(ls mRNA_volcano*.pdf) 
for FILE in $(ls *heatmap*.pdf) 
# for FILE in $(ls hmg4_vs_spt16__ATAC_volcano.pdf) 
# for FILE in $(ls *.pdf) 
do
  echo $FILE
  file_name=$(basename $FILE .pdf)
  pdftoppm -png -rx 300 -ry 300 $FILE > $examples_dir_png/${file_name}.png
done

cd $cactus_dir


## Same but splitting the pages in the files *__{ATAC,mRNA}_other_plots.pdf in multiple png images

cd $examples_dir_pdf

# for FILE in $(ls *.pdf) 
for FILE in $(ls *other_plots*.pdf) 
do
  echo $FILE
  file_name=$(basename $FILE .pdf)
  pdftoppm -png -rx 300 -ry 300 $FILE ${file_name}
done

mv *.png ../png
cd $cactus_dir



# FILE="ctl_1__reads_coverage.pdf"
# file_name=$(basename $FILE .pdf)
# pdftoppm -png -rx 300 -ry 300 $FILE > ${file_name}.png
# 
# mv *.png ../png
