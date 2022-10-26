
cd $cactus


examples_dir_pdf="docs/examples/pdf"
examples_dir_png="docs/examples/png"

# volcano plot without removing specific regions
run_path="test_datasets/worm/results/test_worm__no_rtr/Figures_Individual"
cp $run_path/2_Differential_Abundance/ATAC__volcano/hmg4_vs_spt16__ATAC_volcano.pdf $examples_dir_pdf/hmg4_vs_spt16__ATAC_volcano__no_rtr.pdf 

run_path="test_datasets/worm/results/test_worm/Figures_Individual"
find $run_path -type f -name "*ctl_1*.pdf" -exec cp "{}" $examples_dir_pdf \;

find $run_path -name "*hmg4_vs_ctl_*_volcano.pdf" 
find $run_path -name "*hmg4_vs_ctl_*_PCA_*.pdf" 

find $run_path -name "*hmg4_vs_ctl_*_{volcano|PCA}.pdf" 


find $run_path -name "*ctl_1*.pdf" -exec cp -t docs/examples {} +

find $run_path -name "*hmg4_vs_ctl_*.pdf" 


find $run_path -type f -exec grep -q "*ctl*.pdf" {} \; -exec cp -t $examples_dir_pdf {} +

find $run_path -type f -exec grep -q "*ctl*.pdf" {} \; -exec cp -t $examples_dir_pdf {} +


find /home/shantanu/processed/ -name '*2011*.xml' -exec cp "{}" /home/shantanu/tosend  \;



find $run_path -type f -name "*ctl_1*.pdf"

find $run_path -type f -exec grep -q '*ctl_1*.pdf' {} \; -exec cp -t $examples_dir_pdf {} +


find $run_path -type f -name "*volcano*.pdf"


find $run_path -type f -exec grep -q "*volcano*.pdf" {} \; 

find . -type f -exec grep -q '^beginString' {} \; -exec cp -t /home/user/DestinationFolder {} +



find $run_path -type f -exec grep -q "*ctl_1*.pdf" {} \; 


find $run_path -type f -exec grep -q "*ctl_1*.pdf" {} \; -exec cp -t docs/examples {} +


find $run_path -type f -exec grep -q "*ctl_1*.pdf"

convert $examples_dir_pdf/hmg4_vs_spt16__ATAC_volcano__no_rtr.pdf $examples_dir_png/hmg4_vs_spt16__ATAC_volcano__no_rtr.jpg


cd $examples_dir_png


run_path="test_datasets/worm/results/test_worm/Figures_Individual"



/home/jersal/workspace/cactus/test_datasets/worm/results/test_worm/Figures_Individual/2_Differential_Abundance/ATAC__volcano


test_datasets/worm/results/test_worm/Figures_Individual/2_Differential_Abundance/ATAC__volcano

