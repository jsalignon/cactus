

## Please note that these commands are for singularity but can easily be adjusted
# to other dependencies. For instance, to use conda one would only need to set 
# the CONDA_CACHE_DIR variable and to update the PROFILE variable

## Setting up variables.
# Please udate the variables needed
REFERENCES_DIR="${HOME}/workspace/cactus/references"
SINGULARITY_CACHE_DIR="${HOME}/workspace/singularity_containers"
TOWER_TOKEN="*"
ENABLE_TOWER="false" # set to true if a tower token is provided
TEST_DS_DIR="${HOME}/workspace/cactus/tests" ; mkdir -p $TEST_DS_DIR
SPECIES="worm"
PROFILE="singularity"

## Creating the singularity cache directory
mkdir -p $SINGULARITY_CACHE_DIR

## Setting up the general cactus parameters
cat <<EOT >> ${HOME}/.cactus.config
params.references_dir        = "$REFERENCES_DIR"
params.singularity_cache_dir = "$SINGULARITY_CACHE_DIR"
params.tower_token           = "$TOWER_TOKEN"
params.enable_tower          = $ENABLE_TOWER
EOT

## Exporting Nextflow variables to set the Nextflow version and the java runtime
export NXF_VER=22.10.8
export NXF_OPTS="-Xms1g -Xmx4g"

## Going to the test datasets directory
cd $TEST_DS_DIR

## Downloading a given test dataset and its associated references
nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest \
	--test_datasets --references -profile $PROFILE --species $SPECIES

## Running the test dataset
nextflow run jsalignon/cactus -profile $PROFILE \
	-params-file parameters/full_test.yml -r main -latest

## Same command but this time specifying that we use 8 cores, 16 gb of memory 
# and we restrict the peak assignments filters to "all" peaks and we keep the 
# top 200 results for each differential analysis results.
nextflow run jsalignon/cactus -profile $PROFILE \
	-params-file parameters/full_test.yml -r main -latest \
	--executor_local_cpus 8 --executor_local_memory '16G' \
	--res_dir 'results/almost_full_test'  \
	--split__peak_assignment ['all'] \
	--split__threshold_values [200] \
	--split__threshold_type "rank"

## Running the test dataset using Cactus version 0.8.6
nextflow run jsalignon/cactus -r v0.8.6 -profile $PROFILE \
	-params-file parameters/full_test.yml

