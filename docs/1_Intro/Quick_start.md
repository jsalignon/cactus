

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/data.md), [Design](/docs/3_Inputs/design.md), [Options](/docs/3_Inputs/options.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



# Install and run

The first step is to create the global configuration file *.cactus.config* located in the home folder. This file must indicate the path where to download the references and the singularity containers. Here is an example of a *.cactus.config* file:
```
params.references_dir         = '/home/user/workspace/cactus/references'
params.singularity_images_dir = '/home/user/workspace/singularity_containers'
```

Downloading references and test datasets:
```
nextflow run jsalignon/cactus/scripts/download/download.nf --references --test_datasets --specie worm -r main -latest
```

>**_Note_:** Test datasets are also available for species fly, human and mouse. They can be tested by changing the *--specie* argument.  

Running Cactus (and downloading containers):
```
nextflow run jsalignon/cactus -r main -latest
```

One can update the pipeline using this command:
```
nextflow pull salignon/cactus
```

Results are stored in the folder *results/Cactus_v{version}*.


# Additional details

It is highly recommended to use the [Nextflow Tower tool](https://tower.nf/) to easily monitor pipelines progress.

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in ~/.bashrc or ~./bash_profile):
```
NXF_OPTS='-Xms1g -Xmx4g'
```


# Parameters

Global parameters can be changed in the nextflow.config file. This include output folder names, resources (memory and CPU usage), type of machine to run the script (local or cloud). This is up to the user to set up the optimal nextflow environment according to their analysis platform. Help can be found [here](https://www.nextflow.io/docs/latest/executor.html) for that.
Note: by default, the analysis is cached and will resume were it was before stopping. This can be disabled by setting the parameter *resume* to *false*. 

Analysis parameters can be changed in the run.config file. See the [Configuration](/docs/3_Inputs/Configuration.md) section for more details on available parameters and configuration files. 


# Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since [text stolen from Maxime/Sarek]. On can specify the version of the pipeline using the –version argument this way:
```
nextflow run jsalignon/cactus -r commit_id
```


# Troubleshooting

In general, scrolling through [Nextflow’s documentation](https://www.nextflow.io/docs/latest/index.html) can help resolving most issues.

The general process to resolve a crashing pipeline is to go to the folder indicated in the crash report, launch the appropriate container, and run the lines of codes indicated in the crash report. This way one can try to identify and solve the issue. Example:

```
cd work/59/8a6fb9*
container=$(grep SINGULARITY .command.run | sed 's/.*tmp //g' | sed 's/.img.*/.img/g')
singularity shell -B $cactus_dir -B "$PWD" --containall --cleanenv --home $PWD --workdir $cactus_dir/tmp $container
cat .command.sh
```

Then each command can be run to try to find the error. 

>**_Note_:** The first row of the *.command.sh* file indicates if the script should be run in bash or in R.

If none of these work, an issue can be created on the the cactus github page to report the problem.

