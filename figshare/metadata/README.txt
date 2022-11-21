
***********************
* GENERAL INFORMATION *
***********************

Author:         Jerome Salignon, Lluis Milan Arino, Maxime Garcia, Christian Riedel
Contact e-mail: jerome.salignon@ki.se
License:        MIT
Documentation:  https://github.com/jsalignon/cactus
DOI:            10.17044/scilifelab.20171347



***********************
* DATASET DESCRIPTION *
***********************

This datasets contains the references (*_refs.tar.gz) and the test datasets (*_test.tar.gz) files for the Cactus pipeline. 
See the documentation for details.



***********************
*     REFERENCES      *
***********************

Structure of the references folders:

.
├── bowtie2_indexes_conta
├── CHIP
├── chromatin_states
│   ├── ~one folder per state
├── genome
│   ├── annotation
│   │   ├── bed_regions
│   │   ├── filtered
│   │   └── R
│   └── sequence
│       └── bowtie2_indexes
└── homer_data
    ├── accession
    ├── genomes
    │   └── ~nickname
    ├── GO
    └── promoters
	
	
	
***********************
*    TEST DATASETS    *
***********************

Structure of the test datasets folders:

data
├── atac
    ├── *.fastq.gz
└── mrna
    ├── *.fastq.gz
design
├── atac_fastq.tsv
├── comparisons.tsv
├── genes_to_remove.tsv
├── groups.tsv
├── mrna_fastq.tsv
└── regions_to_remove.tsv
parameters
├── full_test.yml
└── no_enrich.yml
└── ...