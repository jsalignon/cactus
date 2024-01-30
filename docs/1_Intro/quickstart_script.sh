
##########################
## Installing dependencies

# Please observe that the following commands are provided for convenience, but they may not be up to data and it is recommended to follow the installation guidelines from the respective tools (links are provided in the dependencies section of Cactus documentation).

# Variables to modify. Please udate the variables needed and run only these lines.
NEXTFLOW_DIR="$HOME/install/nextflow" ; mkdir -p $NEXTFLOW_DIR
DOCKER_DIR="$HOME/install/docker" ; mkdir -p $DOCKER_DIR
CONDA_DIR="$HOME/install/conda" ; mkdir -p $CONDA_DIR
SINGULARITY_DIR="$HOME/install/singularity" ; mkdir -p $SINGULARITY_DIR
SINGULARITY_VERSION="3.11.3-focal"

# Installing Nextflow version 22.10.8. (Alternatively, one can append the "NXF_VER=22.10.8 " prefix before each call to Cactus.
cd $NEXTFLOW_DIR
wget -qO- https://github.com/nextflow-io/nextflow/releases/download/v22.10.8/nextflow-22.10.8-all | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin

# Installing Singularity
cd $SINGULARITY_DIR
sudo apt-get update
sudo apt-get install -y \
   autoconf \
   automake \
   cryptsetup \
   fuse \
   fuse2fs \
   git \
   libfuse-dev \
   libglib2.0-dev \
   libseccomp-dev \
   libtool \
   pkg-config \
   runc \
   squashfs-tools \
   squashfs-tools-ng \
   uidmap \
   wget \
   zlib1g-dev
export VERSION=1.21.6 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc
export VERSION=$SINGULARITY_VERSION &&  \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
    tar -xzf singularity-ce-${VERSION}.tar.gz && \
    cd singularity-ce-${VERSION}

# Installing Docker
cd $DOCKER_DIR
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world

# installing conda and Mamba
cd $CONDA_DIR
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
yes
cd $CONDA_DIR/miniforge3
conda config --set auto_activate_base false
eval "$(${CONDA_DIR}/bin/conda shell.bash hook 2> /dev/null)"
conda --version # check version
mamba --version # check version

## Installing dependencies
##########################


#########################
## Running a test dataset

# Variables to modify. Please udate the variables needed.
REFERENCES_DIR="${HOME}/workspace/cactus/references"
SINGULARITY_CACHE_DIR="${HOME}/workspace/singularity_containers"
mkdir -p $SINGULARITY_CACHE_DIR
TOWER_TOKEN="*"
ENABLE_TOWER="false" # set to true if a tower token is provided
TEST_DS_DIR="${HOME}/workspace/cactus/tests" ; mkdir -p $TEST_DS_DIR
SPECIES="worm"

# setting up the general cactus parameters
cat <<EOT >> ${HOME}/.cactus.config
params.references_dir        = "$REFERENCES_DIR"
params.singularity_cache_dir = "$SINGULARITY_CACHE_DIR"
params.tower_token           = "$TOWER_TOKEN"
params.enable_tower          = $ENABLE_TOWER
EOT

# updating variables in the bashrc file
echo "export NXF_VER=22.10.8" >> ${HOME}/.bashrc 
export NXF_OPTS='-Xms1g -Xmx4g' >> ${HOME}/.bashrc 

# going to the test datasets directory
cd $TEST_DS_DIR

# downloading a given test dataset and its associated references
nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest \
	--test_datasets --references -profile singularity --species $SPECIES

# running the test dataset using singularity
nextflow run jsalignon/cactus -profile singularity \
	-params-file parameters/full_test.yml -r main -latest

# same command but this time specifying that we use 8 cores, 16 gb of memory 
# and we restrict the peak assignments filters to "all" peaks and we keep the 
#  top 200 results for each differential analysis results.
nextflow run jsalignon/cactus -profile singularity \
	-params-file parameters/full_test.yml -r main -latest \
	--executor_local_cpus 8 --executor_local_memory '16G' \
	--res_dir 'results/almost_full_test'  \
	--split__peak_assignment ['all'] \
	--split__threshold_values [200] \
	--split__threshold_type "rank"

# running the test dataset using conda and specifiying the Cactus version 0.8.6
nextflow run jsalignon/cactus -r v0.8.6 -profile conda \
	-params-file parameters/full_test.yml

## Running a test dataset
#########################


##################
## Troubleshooting

# Killing all nextflow processes running in the current directory and subdirectories
kill_nextflow_processes() {
  kill -9 `ps -aux | grep $(whoami) | grep "${PWD}/work" | awk '{print $2}'`
}
kill_nextflow_processes

# loading a singularity container
WORK_DIR_WITH_BUG=work/59/8a6fb9elkfl39j # to adjust
load_singularity_container() {
  container=$(grep SINGULARITY .command.run | sed 's/.*\/dev\/shm //g' | \
  	sed 's/.img.*/.img/g')
  singularity shell -B /home/jersal --containall --cleanenv --home $PWD \
  --workdir /dev/shm $container
}
cd $WORK_DIR_WITH_BUG
load_singularity_container
cat .command.sh

## Troubleshooting
##################
