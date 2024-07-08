
#### testing the cactus module
# https://docs.uppmax.uu.se/software/software-table/
# https://docs.uppmax.uu.se/cluster_guides/slurm_on_rackham/
# https://www.uu.se/centrum/uppmax/digitalAssets/560/c_560271-l_1-k_uppmax-tutorial_lecture.pdf

## Logging in to UPPMAX
ssh user@rackham.uppmax.uu.se

## Creating a .cactus.config file
cat > .cactus.config <<EOL
params.references_dir        = "/sw/bioinfo/cactus_atac/1.0.0/rackham/refrences"
params.singularity_cache_dir = "/sw/bioinfo/cactus_atac/1.0.0/rackham/singularity_containers"
params.tower_token           = "*"
params.enable_tower          = false
EOL

## Downloading the test dataset
module load bioinfo-tools
module load cactus_atac
nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest \
    --test_datasets -profile singularity --species worm

## making the run script
cat > run_cactus_atac.sh <<EOL
#!/bin/bash
#SBATCH -A "project_id"                # e.g., naiss2024-XX-XXX
#SBATCH -p core                        # Partition name
#SBATCH --job-name=cactus_atac         # Job name
#SBATCH --output=output_%j.txt         # Standard output log
#SBATCH --error=error_%j.txt           # Standard error log
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=6              # Number of CPU cores per task
#SBATCH --mem=16G                      # Memory per node
#SBATCH --time=04:00:00                # Time limit hrs:min:sec

# Load necessary modules
module load bioinfo-tools
module load cactus_atac

nextflow run jsalignon/cactus -profile singularity \
    -params-file parameters/full_test.yml -r main -latest
EOL

## running the test dataset
sbatch run_cactus_atac.sh

# monitoring the run in local
squeue -u user
  #    JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
  #    XXX   core      cactus_a     user  R       4:14      1 current_node
ssh current_node
htop
