
cd $cactus

run_path="test_datasets/worm/results/Cactus_v15.08.22/Figures_Individual"

find $run_path -name "*ctl_1*.pdf" -exec cp -t docs/examples {} +


find $run_path -name "*hmg4_vs_ctl_*.pdf" 

| grep "*ctl_1*.pdf" 


find $run_path -type f -exec grep -q "*ctl_1*.pdf" {} \; -exec cp -t docs/examples {} +


find . -type f -exec grep -q '^beginString' {} \; -exec cp -t /home/user/DestinationFolder {} +