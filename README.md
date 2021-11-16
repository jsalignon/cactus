
# Introduction

CACTUS (Chromatin Accessibility and Transcriptomics Unification Software) is an open-source pipeline designed to easily analyze and visualize gene expression and/or chromatin accessibility data. It can handle any of the four Encode species (human, M. musculus, D. melanogaster and C. elegans). 
The aim of the pipeline is to make it easy for the biologists, and bioinformaticians, to generate hypothesis of which molecular factors regulate gene expression and chromatin accessibility for the studied conditions. It achieves this by providing output that are both easy to use (merged pdf, formatted Excel tables, individual or separate tables/pdf) and to interpret (multiQC, volcano/PCA plots, standardized plots for mRNA and ATAC-Seq, grouped conditions for heatmap plots, …).
Finally, a key aspect of Cactus is the emphasis on efficiency and reproducibility. This is achieved via Singularity, a containerization software, and Nextflow, a dataflow-based pipeline tool that automate parallelism. 

# Pipeline Overview

[ Add the poster of the pipeline that I made. Maybe make it more fancy by putting it in Inkscape]

The pipeline’s input are raw fastq-files and provide detailed plots and tables. ATAC-Seq data is processed following Harvard’s ATAC-seq Guidelines (the original version from 2017; link: https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html). As the pipeline was originally developed to analyze ATAC-Seq data, it contains a plethora of quality controls for checking the quality of the ATAC-Seq data. Now, though, the pipeline can be run using only mRNA-Seq data. MultiQC reports are provided for both kind of data.
For ATAC-Seq, reads are merged, trimmed and aligned to the genome via bowtie2 (ref). Aligned reads are filtered (low quality, duplicates, mitochondrial), shifted (atac shift (ref)) and extended to XX bp. Narrow peaks are then called using Macs2 (ref). These are further split, filtered (blacklist, gDNA) and annotated via CHIPseeker (ref). Finally, Differentially Accessible Regions are determined using DiffBind (ref).
For mRNA seq, transcripts are quantified using Kallisto (ref), and Differential Gene Expression Analysis is carried via the sleuth R package (ref).
Differential Abundance (DA) analysis refers to both Differentially Accessibility Analysis and Differential Gene Expression Analysis. Standardized outputs are produced for all DA results. This include Volcano plots and PCA plots, and DA results tables. DA results are subsequently filtered with user defined filters (False Discovery Rate, peak annotation, Fold Change). Split and filtered DA results are saved as tables.
Split and filtered DA results are then used as input for various enrichment analysis. For both ATAC-Seq and mRNA-Seq genes (DA peaks’ closest gene and DA genes) and genomic regions (DA peaks and DA genes’ promoter) are exported. Genes are used for ontologies and pathway enrichment analysis. Genomic regions are used for chromatin state, Transcription Factor known motifs and ENCODE CHIP enrichment analysis. Additionally, self-overlap enrichment analysis of genes and peaks is performed. Results are saved as tables, and displayed via Venn Diagrams, Barplots and Heatmaps.

For more details see the publication (manuscript under preparation).

Usage

Install
To install Cactus, 
nextflow run salignon/cactus …

One can update the pipeline using this command:
nextflow pull salignon/cactus

And one can specify the version of the pipeline using this command:
… 
In addition to that, one should also download the R containers using this command: …

Inputs

Inputs of the program are raw ATAC-Seq or mRNA-Seq fastq.gz files. If there are multiple files for the same conditions, the read do not necessarily need to be merged as this can be done within cactus. 
The configuration files are expected by Cactus are: atac_seq.config, mRNA_Seq.config, comparisons.config, regions_to_remove.config, grouped_conditions_for_plots.config, base.config. All are tab separated files. Below are instructions and examples on how to write these files [note: add an example per type of file]:
-	atac_seq.config: one line per ATAC-Seq sample. The first entry is the sample ID, and the second is the path were the files is stored. There can be multiple paths, for multiple files, in which cases these files will be merged by Cactus.
-	mRNA_Seq.config: one line per mRNA-Seq sample. Same formatting as for atac_seq.config
-	regions_to_remove.config: This file allows to remove all reads that map to the specified regions. This is useful in particular for experiments involving RNA interference, as this is known to induce a very strong sequencing signal for the locus that is repressed. The file is configurated this way:
-	comparisons.config: this file allows to determine which pairs of conditions will be compared to one another for Differential Abundance Analysis. It is formatted with one comparison per line with one entry per condition. 
-	grouped_conditions_for_plots.config: This file allows to define groups of comparisons to plot together on the heatmaps. The format is one line per group with the group ID as first entry and the comparisons as the remaining entries. The comparisons are named this way: comparison1_vs_comparison2.
-	base.config (optional): this file allows to overwrite Cactus default settings with custom settings. All parameters from the run.config file can be set changed here to determine how a given experiment is analyzed. See the parameters section below for more details.

Run

Here is the command to run the pipeline:
nextflow run salignon/cactus 
A folder “results_pipeline_version” will be created with the results.
It is highly recommended to use the tower tool (https://tower.nf/) to monitor the progress of the pipeline.

Additional details:
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through screen / tmux or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.
It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in ~/.bashrc or ~./bash_profile):
NXF_OPTS='-Xms1g -Xmx4g'


Parameters
Global parameters can be changed in the nextflow.config file. This include output folder names, resources (memory and CPU usage), type of machine to run the script (local or cloud). This is up to the user to set up the optimal nextflow environment according to their analysis platform. Help can be found here for that: https://www.nextflow.io/docs/latest/executor.html.
Note: by default, the analysis is cached and will resume were it was before stopping. This can be disabled by setting the parameter “resume” to “false”. 
Analysis parameters can be changed through the optional base.config file (see the Inputs section). 
Here are the parameters that can be changed and the values they can take: TO DO

Reproducibility
It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since [text stolen from Maxime/Sarek]. On can specify the version of the pipeline using the –version argument this way:

nextflow run salignon/cactus -v XXX

Outputs

[Note: There should be one example Figure for each output]
There are three main types of outputs: processed data, figures and tables. Both figures and tables are available either merged (for easier scrolling) and separated (for lower memory usage). Tables can be either csv files or formatted Excel files.
Here is the structure of the output file with a short description of each output file and examples when relevant: 
TO DO. [add the tree of output as clickable and with a submenu, and description and image for each item]


Troubleshooting

In general, scrolling through Nextflow’s documentation can help resolving most issues: https://www.nextflow.io/docs/latest/index.html

The general process to resolve a crashing pipeline is to go to the folder indicated in the crash report, launch the appropriate container, and run the lines of codes indicated in the crash report. This way one can identify and perhaps solve the issue.

Here is a list of common issues and suggested solutions: …

Citation

Manuscript under preparation.
