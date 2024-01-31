

# Please observe that the following commands are provided for convenience, but 
# they may not be up to date; and it is recommended to follow the installation 
# guidelines from the respective tools (links are provided in the dependencies 
# section of the Cactus documentation).

#######################
## Setting up variables

# Please udate the variables needed and run only these lines.
NEXTFLOW_DIR="$HOME/install/nextflow" ; mkdir -p $NEXTFLOW_DIR
DOCKER_DIR="$HOME/install/docker" ; mkdir -p $DOCKER_DIR
CONDA_DIR="$HOME/install/conda" ; mkdir -p $CONDA_DIR
SINGULARITY_DIR="$HOME/install/singularity" ; mkdir -p $SINGULARITY_DIR
SINGULARITY_VERSION="3.11.3-focal"

## Setting up variables
#######################


######################
## Installing Nextflow

# Here we install Nextflow version 22.10.8, but one can also use an exisiting
# and more recent version of Nextflow and append the "NXF_VER=22.10.8 " prefix 
# before each call to Cactus, or export this variable in the ${HOME}/.bashrc 
# file).

cd $NEXTFLOW_DIR

wget -qO- https://github.com/nextflow-io/nextflow/releases/download/v22.10.8/nextflow-22.10.8-all | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin

## Installing Nextflow
######################


#########################
## Installing Singularity

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

## Installing Singularity
#########################


####################
## Installing Docker 

cd $DOCKER_DIR

sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings

sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg \
	-o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io \
	docker-buildx-plugin docker-compose-plugin
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world

## Installing Docker 
####################


#############################
## Installing conda and Mamba

cd $CONDA_DIR

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
yes

cd $CONDA_DIR/miniforge3
conda config --set auto_activate_base false
eval "$(${CONDA_DIR}/bin/conda shell.bash hook 2> /dev/null)"
conda --version # check version
mamba --version # check version

## Installing conda and Mamba
#############################
