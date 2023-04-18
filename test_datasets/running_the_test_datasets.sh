

##################################################################
## Downloading references and running Cactus on all tests datasets 
## and with all tools manager 


## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/testing

## setting up variables

cpu_nb=47
memory_size='300G'
cactus_version=v0.8.0 ; latest_flag=''

## running Cactus on all tests datasets and with all tools manager 
cd $test_dir 
# rm -r *

for tools_manager in singularity docker conda mamba 
do 
for species in worm fly human mouse
  do
    cd $test_dir
    mkdir -p $tools_manager/$species
    cd $test_dir/$tools_manager/$species
    nextflow run jsalignon/cactus/scripts/download/download.nf -r $cactus_version $latest_flag --test_datasets --references --references_dir refs -profile $tools_manager --species $species
    nextflow run jsalignon/cactus -r $cactus_version $latest_flag -params-file parameters/full_test.yml --references_dir $test_dir/$tools_manager/$species/refs -profile $tools_manager --executor_local_cpus $cpu_nb --executor_local_memory $memory_size --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [1000]
  done
done


## resources needed for the analysis with these parameters:
# cpu_nb=47
# memory_size='300G'

## Worm analysis 
# Duration    : 10m 43s
# CPU hours   : 4.2
# Succeeded   : 874

## Fly analysis 
# Duration    : 21m 36s
# CPU hours   : 9.8
# Succeeded   : 1'496

## Mouse analysis 
# Duration    : 2h 11m 40s
# CPU hours   : 51.5
# Succeeded   : 981

## Human analysis 
# Duration    : 3h 27m 19s
# CPU hours   : 77.5
# Succeeded   : 801


## Downloading references and running Cactus on all tests datasets 
## and with all tools manager 
##################################################################




###############################################################
## Just running the worm dataset without downloading references

# Note: runs for the example figures on the docs are made in this scritp: 
# docs/util/get_figure_examples.sh

## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/testing2

## setting up variables

cpu_nb=47
memory_size='300G'
# cactus_version=v0.8.0 ; latest_flag=''
cactus_version=hotfix/0.8.1 ; latest_flag=''
# cactus_version=1fd6d2054839fa4978d8c8e2f765c90611f36b02 ; latest_flag=''
# cactus_version=main ; latest_flag='-latest'
# tools_manager=singularity
# species=worm

cd $test_dir 
# rm -r *

for tools_manager in singularity 
do 
for species in worm
  do
    cd $test_dir
    mkdir -p $tools_manager/$species
    cd $test_dir/$tools_manager/$species
    nextflow run jsalignon/cactus -r $cactus_version $latest_flag \
      -params-file parameters/full_test.yml -profile $tools_manager \
      --executor_local_cpus $cpu_nb --executor_local_memory $memory_size \
      --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] \
      --split__threshold_values [1000]
  done
done

## Just running the worm dataset without downloading references
###############################################################

