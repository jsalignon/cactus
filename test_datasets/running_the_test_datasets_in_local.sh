


 # => it is important to add the --cactus_dir param otherwise the 
 #    default ${HOME}/.nextflow/assets/jsalignon/cactus is used!




## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/testing

## setting up variables
cpu_nb=30
memory_size='200G'
cactus_version=main ; latest_flag='-latest'
cactus_version=v0.8.5 ; latest_flag=''
cactus_version=updating_the_bioconductor_container ; latest_flag=''


###########################
## Looping across executors

## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
# test_dir=$cactus_dir/test_datasets
# test_dir=$cactus_dir/testing

cpu_nb=30
memory_size='200G'

## running Cactus on all tests datasets and with all tools manager 
cd $test_dir 
# rm -r *

for species in fly
#for species in fly worm
#for species in human mouse
#for species in fly worm human mouse
do 
  for tools_manager in conda singularity docker mamba 
  #for tools_manager in docker mamba
  #for tools_manager in conda singularity
  #for tools_manager in conda singularity mamba
  #for tools_manager in singularity 
  #for tools_manager in conda
  #for tools_manager in mamba
  #for tools_manager in docker
  do
    # species=fly
    # tools_manager=mamba
    cd $test_dir
    mkdir -p $tools_manager/$species
    cur_dir="$test_dir/$tools_manager/$species"
    cd $cur_dir
    ref_dir="$cur_dir/refs"
    pwd ; ls
    nextflow run /home/jersal/workspace/cactus/scripts/download/download.nf \
      --references --references_dir $ref_dir \
      --test_datasets \
      -profile $tools_manager --species $species \
      --cactus_dir $cactus_dir
     

    nextflow run /home/jersal/workspace/cactus/main.nf \
      -params-file parameters/full_test.yml --references_dir $ref_dir \
      -profile $tools_manager --executor_local_cpus $cpu_nb \
      --executor_local_memory $memory_size \
      --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] \
      --split__threshold_values [1000] --cactus_dir $cactus_dir

  done
done


## Looping across executors
###########################



#########################
## Looping across species

## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
#test_dir=$cactus_dir/testing
#test_dir=/home/jersal/test/test_cactus
# test_dir=$cactus_dir/testing

cpu_nb=30
memory_size='200G'

## running Cactus on all tests datasets and with all tools manager 
cd $test_dir 
# rm -r *

#for species in fly
#for species in fly worm
#for species in human mouse
for species in fly worm human mouse
do 
  #for tools_manager in conda singularity docker mamba 
  #for tools_manager in docker mamba
  #for tools_manager in conda singularity
  #for tools_manager in conda singularity mamba
  for tools_manager in singularity 
  #for tools_manager in conda
  #for tools_manager in mamba
  #for tools_manager in docker
  do
    # species=fly
    # tools_manager=mamba
    cd $test_dir
    mkdir -p $tools_manager/$species
    cur_dir="$test_dir/$tools_manager/$species"
    cd $cur_dir
    ref_dir="$cur_dir/refs"
    pwd ; ls
    nextflow run /home/jersal/workspace/cactus/scripts/download/download.nf \
      --references --references_dir $ref_dir \
      --test_datasets \
      -profile $tools_manager --species $species \
      --cactus_dir /home/jersal/workspace/cactus
     
    nextflow run /home/jersal/workspace/cactus/main.nf \
      -params-file parameters/full_test.yml --references_dir $ref_dir \
      -profile $tools_manager --executor_local_cpus $cpu_nb \
      --executor_local_memory $memory_size \
      --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] \
      --split__threshold_values [1000] --cactus_dir /home/jersal/workspace/cactus

  done
done


## Looping across species
#########################



