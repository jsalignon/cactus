
# Creating the Cactus configuration file (if needed)
cat >> ${HOME}/.cactus.config << EOL
params.references_dir        = "${HOME}/workspace/cactus/references"
params.singularity_cache_dir = "${HOME}/workspace/singularity_containers"
params.tower_token           = "*"
params.enable_tower          = false
EOL

# Selecting the Nextflow version
export NXF_VER=22.10.8

# Selecting the Cactus version 
cactus_version=0.9.0

# Downloading the worm references and test dataset
nextflow run jsalignon/cactus/scripts/download/download.nf -r $cactus_version \
	--test_datasets --references -profile singularity --species worm

# change folder and check design
cd worm
cat design/*.yml

