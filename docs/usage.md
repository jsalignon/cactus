
# Menu
* [Introduction](/README.md)
* [Pipeline overview](/docs/pipeline_overview.md)
* [Usage](/docs/usage.md)
* [Output](/docs/output.md)


# Usage Menu

* [Install](#Install)
* [Input](#Input)
* [Run](#Run)
* [Parameters](#Parameters)
* [Reproducibility](#Reproducibility)

# Install
To install Cactus, 
nextflow run salignon/cactus …

One can update the pipeline using this command:
nextflow pull salignon/cactus

And one can specify the version of the pipeline using this command:
… 
In addition to that, one should also download the R containers using this command: …

# Input

Inputs of the program are raw ATAC-Seq or mRNA-Seq fastq.gz files. If there are multiple files for the same conditions, the read do not necessarily need to be merged as this can be done within cactus. 
The configuration files are expected by Cactus are: atac_seq.config, mRNA_Seq.config, comparisons.config, regions_to_remove.config, grouped_conditions_for_plots.config, base.config. All are tab separated files. Below are instructions and examples on how to write these files [note: add an example per type of file]:
-	**atac_seq.config**: one line per ATAC-Seq sample. The first entry is the sample ID, and the second is the path were the files is stored. There can be multiple paths, for multiple files, in which cases these files will be merged by Cactus.
-	**mRNA_Seq.config**: one line per mRNA-Seq sample. Same formatting as for atac_seq.config
-	**regions_to_remove.config**: This file allows to remove all reads that map to the specified regions. This is useful in particular for experiments involving RNA interference, as this is known to induce a very strong sequencing signal for the locus that is repressed. The file is configurated this way:
-	**comparisons.config**: this file allows to determine which pairs of conditions will be compared to one another for Differential Abundance Analysis. It is formatted with one comparison per line with one entry per condition. 
-	**grouped_conditions_for_plots.config**: This file allows to define groups of comparisons to plot together on the heatmaps. The format is one line per group with the group ID as first entry and the comparisons as the remaining entries. The comparisons are named this way: comparison1_vs_comparison2.
-	**base.config** (optional): this file allows to overwrite Cactus default settings with custom settings. All parameters from the run.config file can be set changed here to determine how a given experiment is analyzed. See the parameters section below for more details.

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

Troubleshooting

In general, scrolling through Nextflow’s documentation can help resolving most issues: https://www.nextflow.io/docs/latest/index.html

The general process to resolve a crashing pipeline is to go to the folder indicated in the crash report, launch the appropriate container, and run the lines of codes indicated in the crash report. This way one can identify and perhaps solve the issue.

Here is a list of common issues and suggested solutions: …


# Usage Menu

* [Install](#Install)
* [Input](#Input)
* [Run](#Run)
* [Parameters](#Parameters)
* [Reproducibility](#Reproducibility)

