

# Killing all nextflow processes running in the current directory and 
# subdirectories
kill_nextflow_processes() {
  kill -9 `ps -aux | grep $(whoami) | grep "${PWD}/work" | awk '{print $2}'`
}
kill_nextflow_processes

# loading a singularity container
WORK_DIR_WITH_BUG=work/59/8a6fb9elkfl39j # set path to the failing directory
load_singularity_container() {
  container=$(grep SINGULARITY .command.run | sed 's/.*\/dev\/shm //g' | \
  	sed 's/.img.*/.img/g')
  singularity shell -B /home/jersal --containall --cleanenv --home $PWD \
  --workdir /dev/shm $container
}
cd $WORK_DIR_WITH_BUG
load_singularity_container
cat .command.sh

