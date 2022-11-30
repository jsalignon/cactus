
# setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/testing

# setting up variables
cactus_version=0.9.0
cpu_nb=47
memory_size='250G'

# running Cactus on all tests datasets and with all tools manager 
cd $test_dir 
# rm -r *

for tools_manager in singularity docker conda mamba 
do 
for species in worm fly human mouse
  do
    cd $test_dir
    mkdir -p $tools_manager/$species
    cd $test_dir/$tools_manager/$species
    nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest --test_datasets --references --references_dir refs -profile $tools_manager --species $species -resume
    nextflow run jsalignon/cactus -r main -latest -params-file parameters/full_test.yml --references_dir $test_dir/$tools_manager/$species/refs -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200] -resume
  done
done

