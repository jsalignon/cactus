
* Introduction: [Cactus](/README.md), [Quick Start](/docs/Intro__Quick_start.md), 
* Overview: [Graph](/docs/Overview__Graph.md), [Outputs structure](/docs/Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/Install__Dependencies.md), [Containers](/docs/Install__Containers.md), [Data](/docs/Install__Data.md), [Test_datasets](/docs/Install__Test_datasets.md)
* Run: [Input Data](/docs/ID__Input_data.md), [Input Files](/docs/ID__Input_files.md), [Parameters](/docs/ID__Parameters.md)
* ATAC Preprocessing: [Figures](/docs/AP__Figures.md), [MultiQC](/docs/MP__MultiQC.md)
* mRNA Preprocessing: [MultiQC](/docs/Prepro__MultiQC.md)
* Differential Abundance: [ATAC](/docs/DA__DiffBind.md), [mRNA](/docs/DA__Sleuth.md), [Figures](/docs/DA__Figures.md), [Tables](/docs/DA__Tables.md)
* Splitting peak sets: [Split](/docs/SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/Enrich__CHIP.md), [Chromatin state](/docs/Enrich__Chromatin_states.md), [Motifs](/docs/Enrich__Motifs.md), [Func. Anno.](/docs/Enrich__Functional_annotations.md), [Figures](/docs/Enrich__Figures.md), [Tables](/docs/Enrich__Tables.md)


[](END_OF_MENU)






# Install
To install Cactus, 
nextflow run salignon/cactus …

One can update the pipeline using this command:
nextflow pull salignon/cactus

And one can specify the version of the pipeline using this command:
… 
In addition to that, one should also download the R containers using this command: …

# Run

Here is the command to run the pipeline:
nextflow run salignon/cactus 
A folder “results_pipeline_version” will be created with the results.
It is highly recommended to use the tower tool (https://tower.nf/) to monitor the progress of the pipeline.

Additional details:
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.
It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in ~/.bashrc or ~./bash_profile):
NXF_OPTS='-Xms1g -Xmx4g'


# Parameters
Global parameters can be changed in the nextflow.config file. This include output folder names, resources (memory and CPU usage), type of machine to run the script (local or cloud). This is up to the user to set up the optimal nextflow environment according to their analysis platform. Help can be found here for that: https://www.nextflow.io/docs/latest/executor.html.
Note: by default, the analysis is cached and will resume were it was before stopping. This can be disabled by setting the parameter “resume” to “false”. 
Analysis parameters can be changed through the optional base.config file (see the Inputs section). 
Here are the parameters that can be changed and the values they can take: TO DO

# Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since [text stolen from Maxime/Sarek]. On can specify the version of the pipeline using the –version argument this way:

nextflow run salignon/cactus -v XXX

# Troubleshooting

In general, scrolling through Nextflow’s documentation can help resolving most issues: https://www.nextflow.io/docs/latest/index.html

The general process to resolve a crashing pipeline is to go to the folder indicated in the crash report, launch the appropriate container, and run the lines of codes indicated in the crash report. This way one can identify and perhaps solve the issue.

Here is a list of common issues and suggested solutions: …

