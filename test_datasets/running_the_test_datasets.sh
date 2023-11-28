


#################
## Initialization

# please change paths and variables if needed

## setting up paths
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/testing

## setting up variables
cpu_nb=30
memory_size='300G'
cactus_version=main ; latest_flag='-latest'

cactus_version=updating_the_bioconductor_container ; latest_flag=''
# nextflow drop jsalignon/cactus # if needed

## Initialization
#################




###########################
## Looping across executors

# cd $test_dir ; rm -r *  # if needed

for species in worm
do 
  for tools_manager in conda singularity docker mamba 
  do
    # species=worm        # if needed
    # tools_manager=mamba # if needed

    cd $test_dir
    mkdir -p $tools_manager/$species
    cur_dir="$test_dir/$tools_manager/$species"
    cd $cur_dir
    ref_dir="$cur_dir/refs"
    pwd ; ls

    nextflow run jsalignon/cactus/scripts/download/download.nf \
      -r $cactus_version $latest_flag \
      -profile $tools_manager \
      --references_dir $ref_dir \
      --references \
      --test_datasets \
      --species $species

    nextflow run jsalignon/cactus/main.nf \
      -r $cactus_version $latest_flag \
      -profile $tools_manager \
      --references_dir $ref_dir \
      --res_dir 'results/almost_full_test'  \
      -params-file parameters/full_test.yml \
      --executor_local_cpus $cpu_nb \
      --executor_local_memory $memory_size \
      --split__peak_assignment ['all'] \
      --split__threshold_values [1000]

  done
done

## Looping across executors
###########################



#########################
## Looping across species

# cd $test_dir ; rm -r *  # if needed

for species in worm fly mouse human
do 
  for tools_manager in singularity 
  do
    # species=fly         # if needed
    # tools_manager=mamba # if needed

    cd $test_dir
    mkdir -p $tools_manager/$species
    cur_dir="$test_dir/$tools_manager/$species"
    cd $cur_dir
    ref_dir="$cur_dir/refs"
    pwd ; ls

    nextflow run jsalignon/cactus/scripts/download/download.nf \
      -r $cactus_version $latest_flag \
      -profile $tools_manager \
      --references_dir $ref_dir \
      --references \
      --test_datasets \
      --species $species

    nextflow run jsalignon/cactus/main.nf \
      -r $cactus_version $latest_flag \
      -profile $tools_manager \
      --references_dir $ref_dir \
      --res_dir 'results/almost_full_test'  \
      -params-file parameters/full_test.yml \
      --executor_local_cpus $cpu_nb \
      --executor_local_memory $memory_size \
      --split__peak_assignment ['all'] \
      --split__threshold_values [1000]

  done
done

## Looping across species
#########################


## resources needed for the analysis with these parameters:
# cpu_nb=30
# memory_size='300G'

## Worm analysis 
# Duration    : 21m
# CPU hours   : 7.9h
# Succeeded   : 1,212

## Fly analysis 
# Duration    : 48m
# CPU hours   : 20.1h
# Succeeded   : 1,793

## Mouse analysis 
# Duration    : 4h13
# CPU hours   : 85
# Succeeded   : 1,164

## Human analysis 
# Duration    : 5h54
# CPU hours   : 117
# Succeeded   : 1,115

