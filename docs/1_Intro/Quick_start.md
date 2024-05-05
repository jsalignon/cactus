
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Tutorial](/docs/1_Intro/tutorial.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# Quickstart scripts

Scripts are available to [install the dependencies](/docs/1_Intro/installing_dependencies.sh), [runinng the test datasets](/docs/1_Intro/running_the_test_datasets.sh) and [troubleshooting](/docs/1_Intro/troubleshooting.sh). Please observe that the script to install dependencies is provided for convenience, but the commands may not be up to date; and it is recommended to follow the installation guidelines from the respective tools (links are provided in the [Dependencies](/docs/2_Install/Dependencies.md) section of the Cactus documentation).


# Dependencies and profiles

Cactus needs two software to be installed in order to run: Nextflow and one of SingularityCE, Docker, conda or Mamba. Please read the [Dependencies](/docs/2_Install/Dependencies.md) section for details on versions and installation of these software.
Then, the *-profile* argument should be used to specicify which tools manager to use. In general, it is recommended to use SingularityCE on HPC systems since Singularity containers can be [run without sudo](https://blogs.oregonstate.edu/learningbydoing/2022/01/04/docker-and-singularity-containers-which-one-is-better/) and Singularity images are [immutable](https://singularity-docs.readthedocs.io/en/latest/) which ensures a high level of reproducibility and verification of images (see also [here](https://spiediedocs.binghamton.edu/docs/conda_singularity_modules.html)). Users should see with their administrator which of these 4 options are available and recommended.


# Install and run

The first step is to create the global configuration file *.cactus.config* located in the user's home folder. This file must indicate the path where to download the references and the singularity containers. Here is an example of a *.cactus.config* file when using Singularity with a Tower token to monitor the pipeline runs on [Nextflow Tower](https://cloud.tower.nf/) (which is highly recommended):
```
params.references_dir        = "${HOME}/workspace/cactus/references"
params.singularity_cache_dir = "${HOME}/workspace/singularity_containers"
params.tower_token           = "*"
params.enable_tower          = true
```

One can then download references and test datasets for *C. elegans* with this command:
```
NXF_VER=22.10.8 nextflow run jsalignon/cactus/scripts/download/download.nf -r main -latest --test_datasets --references -profile singularity --species worm
```

This should download the references and the containers/virtual environments in user-specified parameters (i.e., params.references_dir for reference and params.singularity_cache_dir for singularity images), and it should download the test dataset in the current folder which takes the form of 3 folders: parameters, design and data. Please read the [Inputs](/docs/3_Inputs/3_Inputs.md) section for details on what the folders contain.

>**_Note_:** Test datasets are also available for the species *d. melanogaster* (fly), *m. musculus* (mouse), and *h. sapiens* (human). They can be tested by changing the *--species* argument.  
>**_Note_:** On some platforms, Nextflow should be called with a './' before (i.e., *./nextflow run...*).

Then, one can run Cactus and download all container images/virtual environments using this command:
```
NXF_VER=22.10.8 nextflow run jsalignon/cactus -profile singularity -params-file parameters/full_test.yml -r main -latest
```

One can update the pipeline using this command:
```
nextflow pull jsalignon/cactus
```

Results are stored in the folder `results/Cactus_v${cactus_version}` (this path can be changed with the parameter *params.res_dir*).

It is recommended to use either the worm or the fly test datasets when testing Cactus on a laptop to reduce runtime. With 8 cores and 16Gb RAM the worm and fly test dataset can be run in respectively ~27 and ~56 minutes using this command:
```
NXF_VER=22.10.8 nextflow run jsalignon/cactus -r main -latest -params-file parameters/full_test.yml -profile singularity --executor_local_cpus 8 --executor_local_memory '16G' --res_dir 'results/almost_full_test'  --split__peak_assignment ['all'] --split__threshold_values [200]
```

>**_Note_:** The run parameters can be set up in a *.yml* file or in the command line (as shown just above). The latter taking priority on the former. When setting parameters on the command line, one dash indicates [Nextflow's internal parameters](https://www.nextflow.io/docs/latest/cli.html#run) (e.g. -profile) and two dashes indicate [Cactus' own parameters](/docs/3_Inputs/Parameters.md) (e.g. res_dir). 

>**_Note_:** A minimum of 6 cores is required to run Cactus. 


# Additional details

It is recommended to use [Nextflow Tower](https://tower.nf/) to easily monitor pipelines progress.

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to [limit the Nextflow Java virtual machines memory](https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html)). We recommend adding the following line to your environment (typically in ~/.bashrc or ~./bash_profile):
```
NXF_OPTS='-Xms1g -Xmx4g'
```


# Parameters

Global parameters can be changed in the nextflow.config file. This include output folder names, resources (memory and CPU usage), type of machine to run the script (local or cloud). This is up to the user to set up the optimal nextflow environment according to their analysis platform. Help can be found [here](https://www.nextflow.io/docs/latest/executor.html) for that.

Analysis parameters can be changed in the yml input file. See the [Parameters](/docs/3_Inputs/Parameters.md) section for more details on parameter files and on available parameters. 


# Reproducibility

It's a good idea to specify a [release version](https://github.com/jsalignon/cactus/releases) when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since. One can specify the version of the pipeline using the *-r* argument this way:
```
NXF_VER=22.10.8 nextflow run jsalignon/cactus -r release_tag -profile {singulariy,docker,conda,mamba} -params-file parameter_file
```
For instance:
```
NXF_VER=22.10.8 nextflow run jsalignon/cactus -r v1.0.0 -profile singularity -params-file parameters/full_test.yml
```

# Troubleshooting

In general, scrolling through [Nextflow’s documentation](https://www.nextflow.io/docs/latest/index.html) can help resolving most issues.  

The general process to resolve a crashing pipeline is to go to the folder indicated in the crash report, launch the appropriate container, and run the lines of codes indicated in the crash report. This way one can try to identify and solve the issue. For finer inspection of the code and analysis, a good idea is to run cactus in the background to get a detailled log file. Note that the [*-dump-channels* argument](https://www.nextflow.io/docs/latest/cli.html#run) can also be used to explore channel contents.

The *-bg* argument can be used to run cactus in the background like this:
```
NXF_VER=22.10.8 nextflow run jsalignon/cactus -r main -latest -profile singularity -params-file parameters/full_test.yml -bg > nf_log.txt
```

This creates a .nextflow.pid file that contains the master PID to kill to stop the run in the background. However, this does not always work. A workaround to kill all running processed from the current run folder is to use this function:
```
kill_nextflow_processes() {
  kill -9 `ps -aux | grep $(whoami) | grep "${PWD}/work" | awk '{print $2}'`
}
kill_nextflow_processes
```

Then, one can inspect/grep the nf_log.txt file to go to the folder that we want to inspect in more details. Once in the appropriate folder, the following function can be used, if one uses singularity, to open a shell with the container in the same settings as in Cactus and displaying the set of commands that were ran (in the .command.sh file): 

```
load_singularity_container() {
  container=$(grep SINGULARITY .command.run | sed 's/.*\/dev\/shm //g' | sed 's/.img.*/.img/g')
  singularity shell -B /home/jersal --containall --cleanenv --home $PWD --workdir /dev/shm $container
}
cd work/59/8a6fb9*
load_singularity_container
cat .command.sh
```

Then each command can be run to try to find the error or to inspect the code.

>**_Note_:** The first row of the *.command.sh* file indicates if the script should be run in bash or in R.

If none of these work, an issue can be created on the the cactus GitHub page to report the problem.

