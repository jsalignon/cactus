


//// list of processes (not in order)
// get_fasta_and_gff
// get_blacklisted_regions
// get_hihmm_chromatin_state_data_part_1
// get_hihmm_chromatin_state_data_part_1
// get_homer_data
// getting_orgdb
// get_chip_data
// indexing_genomes
// filtering_annotation_file
// get_hihmm_chromatin_state_data_part_2
// getting_chromosome_length_and_transcriptome
// generating_bowtie2_indexes
// get_hihmm_chromatin_state_data_part_2
// generating_kallisto_transcriptome_indexes
// get_bed_files_of_annotated_regions
// getting_R_annotation_files


params.number_of_cores = 8   // how to initialize if missing?
number_of_cores = params.number_of_cores

params.ensembl_release  = '102'
ncbi_orgdb_version      = '3.14'
homer_genomes_version   = '6.4'
homer_organisms_version = '6.3'
homer_promoters_version = '5.5'

params.homer_odd_score_threshold = 0.5

Assembly_Ensembl = Channel
	.from(
		[   //  0          1                     2              3                4     
		//                                     assembly      assembly          Ensembl            
		// specie    scientific_name            name           type            Release
			['worm',  'caenorhabditis_elegans',  'WBcel235', 'toplevel',         '107' ],
			['fly',   'drosophila_melanogaster', 'BDGP6.32', 'toplevel',         '107' ],
			['mouse', 'mus_musculus',            'GRCm38',   'primary_assembly', '102' ],
			['human', 'homo_sapiens',            'GRCh38',   'primary_assembly', '107' ]
		]
	)
	.multiMap { it ->
					fasta: it
					 pwms: it[0, 1]
					orgdb: it[0, 1]
	}
	.set { Assembly_Ensembl_1 }

//// checking Ensembl latest release and fly genome build name (BDGP6.XX)
// https://www.ensembl.org/index.html?redirect=no
// Ensembl Release 107 (Jul 2022)
// http://ftp.ensembl.org/pub/release-107/gff3/drosophila_melanogaster/

//// UCSC releases and details:
// https://genome.ucsc.edu/FAQ/FAQreleases.html

//// A note on the mouse genome:
// ENCODE (chip and chromatin states) and blacklisted regions are not available yet fo mm11/GRCm39. Therefore, we need to stick with release 102, which is the latest release for mm10/GRCm38. I could eventually use liftover files in the future to remedy this issue, or just wait for the correct assembly files to be released.


Assembly_nickname = Channel
	.from(
		[
			['worm',  'ce11', 'WBcel235',],
			['fly',   'dm6',  'BDGP6',],
			['mouse', 'mm10', 'GRCm38',  ],
			['human', 'hg38', 'GRCh38',  ]
		]
	)
	.multiMap { it ->
			blacklist: it
			  	homer: it[0, 1]
	}
	.set { Assembly_nickname_1 }

Assembly_ENCODE = Channel
	.from( 	
		[  
			['worm',   'ce11' ],
			['fly',    'dm6' ],
			['mouse',  'mm10' ],
			['human',  'GRCh38' ]
		]
	)
	.into{ Assembly_ENCODE_chip ; Assembly_ENCODE_chromatin_state }



HiHMM_liftover = Channel
	.from( 	
		[
			['human', 'hg19', 'hg19ToHg38'],
			[ 'worm', 'ce10', 'ce10ToCe11'],
			[  'fly',  'dm3', 'dm3ToDm6']
		]
	)

HiHMM_chromatin_states = Channel
	.from( 	
		[
			['human', 'iHMM.M1K16.human_GM'],
			['human', 'iHMM.M1K16.human_H1'],
			[  'fly', 'iHMM.M1K16.fly_EL'],
			[  'fly', 'iHMM.M1K16.fly_L3'],
			[ 'worm', 'iHMM.M1K16.worm_EE'],
			[ 'worm', 'iHMM.M1K16.worm_L3']
		]
	)
	




process get_homer_data {
	tag "${specie}"

	container = params.samtools_bedtools_perl

	publishDir path: "${specie}/homer_data", mode: 'link'
	input:
		set specie, specie_code from Assembly_nickname_1.homer

	output:
		file('*')

	shell:
	'''
		
		specie='!{specie}'
		specie_code='!{specie_code}'
		homer_genomes_version='!{homer_genomes_version}'
		homer_organisms_version='!{homer_organisms_version}'
		homer_promoters_version='!{homer_promoters_version}'
		
		URL=http://homer.ucsd.edu/homer/data
		
		wget -nc $URL/genomes/${specie_code}.v${homer_genomes_version}.zip
		wget -nc $URL/organisms/${specie}.v${homer_organisms_version}.zip
		wget -nc $URL/promoters/${specie}.v${homer_promoters_version}.zip
		
		for z in *.zip; do unzip "$z"; done
		mv data/* .
		rm -r data *.zip
		 
	'''
}



process get_pwms {

	container = params.r_basic

	publishDir path: "util/motifs_PWMS", mode: 'link'

	output:
		file('dt_cisbp_encode.rds') into cisbp_motifs_all_species
		file('*')
	
	shell:
	'''
		#!/usr/bin/env Rscript
		
		homer_odd_score_threshold = !{params.homer_odd_score_threshold}
		
		url="http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset"
		f1 = 'PWMs.zip'               ; download.file(paste0(url, '/', f1), f1)
		f1 = 'TF_Information.txt.zip' ; download.file(paste0(url, '/', f1), f1)
		
		unzip('TF_Information.txt.zip')
		unzip('PWMs.zip')
		
		system("sed -i 's/#/|/g' TF_Information.txt")
		
		library(data.table)
		library(magrittr)
		library(purrr)
	
		dt_all = fread('TF_Information.txt')
		
		encode_species = c('Caenorhabditis_elegans', 'Drosophila_melanogaster', 'Mus_musculus', 'Homo_sapiens')
		dt_cisbp = dt_all[TF_Species %in% encode_species]
		dt_cisbp = dt_cisbp[Motif_ID != '.']
		dt_cisbp = dt_cisbp[!duplicated(TF_ID)]

		# removing emtpy motifs
		nrow(dt_cisbp) # 2936
		nrows = sapply(dt_cisbp$Motif_ID, function(x) nrow(read.table(paste0('pwms/', x, '.txt'), sep = '\t')))
		dt_cisbp = dt_cisbp[nrows != 1]
		nrow(dt_cisbp) # 2779

		# loading motifs and determining odd scores, thresholds and consensus sequences
		dt_cisbp$motif = lapply(dt_cisbp$Motif_ID, function(x) read.table(paste0('pwms/', x,'.txt'), sep = '\t', header = T, stringsAsFactors = F)[, -1])
		dt_cisbp[, log2_odd_score := sapply(motif, function(motif1) sum(apply(motif1, 1, function(cur_row) log2(max(cur_row)/0.25))))]
		dt_cisbp[, homer_threshold := log2_odd_score * homer_odd_score_threshold]
		dt_cisbp[, consensus_sequence := sapply(motif, function(motif1) paste0(colnames(motif1)[apply(motif1, 1, function(x) which.max(x))], collapse = ''))]

		saveRDS(dt_cisbp, 'dt_cisbp_encode.rds')
	
		
	'''
	
}

Assembly_Ensembl_1.pwms
	.combine(cisbp_motifs_all_species)
	.set{ cisbp_motifs_all_species_1 }


process split_pwms {

	container = params.r_basic

	publishDir path: "${specie}/homer_data", mode: 'link'

	input:
		set specie, specie_long, file(dt_cisbp_encode_rds) from cisbp_motifs_all_species_1

	output:
		file('homer_motifs.txt')
	
	shell:
	'''
		#!/usr/bin/env Rscript
		
		library(data.table)
		
		dt_cisbp_encode = readRDS('!{dt_cisbp_encode_rds}')
		specie_long = '!{specie_long}'
		
		specie_long_1 = paste0(toupper(substr(specie_long, 1, 1)), substr(specie_long, 2, nchar(specie_long)))
		dt_cisbp = dt_cisbp_encode[TF_Species == specie_long_1]
		
		
		sink('homer_motifs.txt')
		
			apply(dt_cisbp, 1, function(x) {
					cat(paste0('>', x$consensus_sequence, '\t', x$TF_Name, '\t', x$homer_threshold, '\n'))
					cat(paste0(apply(x$motif, 1, paste0, collapse = '\t'), '\n'))
				})

		sink()
		
	'''
	
}




// bioconductor-universalmotif:1.8.3--r40h399db7b_0

// dt_encode = dt_all[TF_Species %in% encode_species]
// dt_encode = dt_encode[Motif_ID != '.']
// dt_encode1 = dt_encode[!duplicated(TF_ID)]
// > dt_encode$TF_Species %>% table
// .
//  Caenorhabditis_elegans Drosophila_melanogaster            Homo_sapiens
//                     402                    1301                    5422
//            Mus_musculus
//                    1503
// >
// > dt_encode1$TF_Species %>% table
// .
//  Caenorhabditis_elegans Drosophila_melanogaster            Homo_sapiens
//                     373                     425                    1200
//            Mus_musculus
//                     938

// However, there are still multiple entries per TF
// dt_encode[TF_Species == 'Caenorhabditis_elegans']$TF_Name %>% table %>% sort %>% rev %>% head(10)
  // hlh-1   skn-1  snpc-4   pha-4   efl-1   tra-1 nhr-182   mab-3  lin-39   eor-1
  //     7       4       3       3       3       2       2       2       2       2
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'skn-1']
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'hlh-1']
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'snpc-4']
// dt_encode[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'pha-4']

// I think I understand what happened: the entries that are duplicated are all from the similarity regression method. And they all have the same score. So instead of picking one arbitrarily, CISBP kept all of these entries with the same score. I could keep only the CHIP.
// Conclusions: these are probably exactly the same motifs (or very very similar). So we can just keep the first of them!


// # Identifying and removing empty motifs files
// file_is_correct = map_lgl(paste0('pwms/', dt_encode$Motif_ID, '.txt'), function(f1) {ifelse(nrow(fread(f1)) == 1, F, T)})
// file_is_correct %>% table # => all TRUE
// TRUE
// 8628

// https://github.com/TedCCLeung/PhMotif/blob/9ace03a196e02375b3024af3506b9f58e98f1a1d/data-raw/processCISBPfiles.R


// nice explanations about the cisbp files here:
// https://github.com/Danko-Lab/rtfbs_db/blob/a29f73b24e8984aad8bcd1227b06446ebef44b08/rtfbsdb/man/CisBP.getTFinformation.Rd
// Three TF information files in CisBP dataset.\cr\cr
// 1: TF_Information.txt : (direct motifs) or (no direct but inferred motifs with 90\%)\cr
// 2: TF_Information_all_motifs.txt: (direct motifs) and (inferred motifs above the threshold)\cr
// 3: F_Information_all_motifs_plus.txt: All motifs\cr

// 'TF_Information_all_motifs.txt' is a superset of 'TF_Information.txt'.  It also includes any motif that can be inferred for a given TF, given the TF family-specific threshold.  For example, if a TF has a directly determined motif, and two TFs with motifs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the direct motif, but 
// TF_Information_all_motifs.txt will include all three motifs.  Likewise, if a TF does not have a direct motif, but has two TFs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the motif from the 90% indentical TF, while TF_Information_all_motifs.txt would include both.\cr\cr



// f1 = 'TF_Information_all_motifs.txt.zip' ; download.file(paste0(url, '/', f1), f1)
// unzip('TF_Information_all_motifs.txt.zip')
// system("sed -i 's/#/|/g' TF_Information_all_motifs.txt")
// dt_all = fread('TF_Information_all_motifs.txt')







process get_blacklisted_regions {
	tag "${specie}"

	container = params.cvbio

	publishDir path: "${specie}/genome/annotation", mode: 'link'

	input:
		set specie, specie_code, ncbi_code from Assembly_nickname_1.blacklist

	output:
		file("blacklisted_regions.bed")

	shell:
	'''
			
		specie_code='!{specie_code}'
		ncbi_code='!{ncbi_code}'
		
		url_blacklist="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists"
		url_mapping="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master"

		wget -O ${specie_code}_blacklist_NCBI.bed.gz $url_blacklist/${specie_code}-blacklist.v2.bed.gz
		gunzip  ${specie_code}_blacklist_NCBI.bed.gz
		
		wget -O NCBI_to_Ensembl.txt $url_mapping/${ncbi_code}_UCSC2ensembl.txt
		
		cvbio UpdateContigNames -i ${specie_code}_blacklist_NCBI.bed -o blacklisted_regions.bed -m NCBI_to_Ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
		
	'''
}

///// blacklist files
// https://github.com/Boyle-Lab/Blacklist/tree/master/lists
// ce11, hg38, dm6, mm10
// mm11 is missing for now so we don't use the latest genome assembly version

//// chromosome mapping files
// https://github.com/dpryan79/ChromosomeMappings
// human: GRCh38_UCSC2ensembl.txt
// mouse: GRCm38_UCSC2ensembl.txt
//   fly: BDGP6_ensembl2UCSC.txt 
//  worm: WBcel235_UCSC2ensembl.txt

// this lines fails with a weird Nextflow bug of variable substitution:
// ncbi_code_1=${ncbi_code//\.*/}
// Script compilation error
// - file : /home/jersal/workspace/cactus/software/get_data/parse_genome_files.nf
// - cause: Unexpected character: '\'' @ line 96, column 4.


// AnnotationDbi::saveDb(org.Ce.eg.db, file = 'org_db_ce.sqlite')
// saveRDS(cur_seq_info, '~/workspace/cactus/data/worm/genome/sequence/cur_seqinfo.rds')
// 
// cur_seq_info = rtracklayer::SeqinfoForUCSCGenome('ce11')
// cur_seq_info@seqnames %<>% gsub('chr', '', .)






process get_chip_metadata {

	container = params.encodeexplorer

	input:
		set specie, assembly from Assembly_ENCODE_chip

	output:
		set specie, file("dt_encode_chip.rds") into Encode_chip_metadata		
		set specie, file("*__split_dt.rds") into Encode_chip_split_dt mode flatten
		

	shell:
	'''

		#!/usr/bin/env Rscript
		
		source('!{params.cactus_dir}/software/get_data/bin/encode_chip_functions.R')
		assembly = '!{assembly}'
		
		library(data.table)
		library(magrittr)
		library(purrr)
		
		url_encode = 'https://www.encodeproject.org' 
		url_search = paste0(url_encode, '/search/?')
		url_append = '&frame=object&format=json&limit=all'


		## getting all tables
		query_chip_files = paste0('type=File&file_format=bed&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&status=released&assembly=', assembly)
		df_chip_files = get_encode_df(query_chip_files) %T>% pnrow # 2716
		df_experiments = get_encode_df('type=Experiment&assay_title=TF+ChIP-seq&status=released') %T>% pnrow # 4425
		df_biosample_types = get_encode_df('type=BiosampleType') %T>% pnrow # 936
		df_genes = get_encode_df('type=Gene&targets=*') %T>% pnrow # 5989
		
		## merging the genes and target tables
		df_targets_0 = get_encode_df('type=Target') %T>% pnrow # 9941
		
		df_targets = df_targets_0
		df_targets = df_targets[!map_lgl(df_targets$genes, is.null), ] %T>% pnrow # 9885
		# few entries target multiple genes, we remove them
		map_int(df_targets$genes, length) %>% table
		df_targets$genes %<>% map_chr(~.x[1])
		dt_targets = data.table(df_targets[, c('targets', 'genes', 'organism')])
		dt_targets[, organism := organism %>% gsub('/organisms/|/', '', .)]
		dt_genes = data.table(df_genes[, c('genes', 'geneid', 'symbol')])
		colnames(dt_genes)[2:3] = c('target_id', 'target_symbol')
		dt_genes_target = dt_genes[dt_targets, , on = 'genes'][, genes := NULL] %T>% pnrow # 9885
		
		## adding the target gene annotations to the experiment table
		dt_experiment = data.table(df_experiments[, c('experiments', 'target', 'biosample_ontology', 'biosample_summary', 'life_stage_age')])
		colnames(dt_experiment)[2] = 'targets'
		dt_experiment1 = dt_genes_target[dt_experiment, , on = 'targets'][, targets := NULL]
		
		## adding the ontology annotations to the experiment table
		dt_biosample_types = data.table(
		  df_biosample_types[, c('biosample_types', 'classification', 'term_name')],
		  cell_slims  = collapse_slims(df_biosample_types$cell_slims),
		  organ_slims = collapse_slims(df_biosample_types$organ_slims),
		  developmental_slims = collapse_slims(df_biosample_types$developmental_slims),
		  system_slims = collapse_slims(df_biosample_types$system_slims)
		)
		colnames(dt_biosample_types)[1] = 'biosample_ontology'
		dt_experiment2 = dt_biosample_types[dt_experiment1, , on = 'biosample_ontology'][, biosample_ontology := NULL]
		
		
		## adding the experiment table to the file table to get the final metadata table
		dt = data.table(df_chip_files[, c('accession', 'dataset', 'assembly', 'href', 'md5sum')])
		colnames(dt)[1:2] = c('file', 'experiments')
		dt = dt_experiment2[dt, , on = 'experiments'][, experiments := NULL] %T>% pnrow
		colnames(dt)[colnames(dt) == 'term_name'] = 'ontology'
		
		# renaming the human assembly and removing entries with missing annotations
		dt[assembly == 'GRCh38', assembly := 'hg38']
		dt = dt[!is.na(target_symbol)] %T>% pnrow # 2714
		
		## adding a column indicating either the cell line or tissue (human or mice) or stage (yeast or worms)
		dt[, stage_cell := biosample_summary]
		
		# for worm
		dt[organism == 'celegans', stage_cell := biosample_summary]
		dt[organism == 'celegans', stage_cell := stage_cell %>% gsub('.*whole organism ', '', .) %>% gsub('.*hermaphrodite ', '', .) %>% capture_group('(L[1-4]).*') %>% gsub('young adult.*', 'YA', .) %>% gsub('late embryo.*', 'LE', .) %>% gsub('early embryo.*', 'EE', .) %>% gsub('mixed stage (embryo)', 'mixed stage embryo', ., fixed = T) %>% gsub('mixed stage embryo.*', 'MS', .) %>% gsub('dauer.*', 'Da', .) %>% gsub('midembryo.*', 'ME', .)]
		dt[organism == 'celegans', ]$stage_cell %>% table %>% sort %>% rev
		
		# for fly
		dt[organism == 'dmelanogaster', stage_cell := biosample_summary]
		dt[organism == 'dmelanogaster', stage_cell := stage_cell %>% gsub('.*whole organism male adult.*', 'MA', .) %>% gsub('.*whole organism female adult.*', 'FA', .) %>% gsub('.*whole organism adult.*', 'MS', .) %>% gsub('.*whole organism embryo.*', 'Em', .) %>% gsub('.*whole organism prepupa.*', 'PP', .) %>% gsub('.*whole organism wandering third.*', 'WT', .) %>% gsub('.*whole organism pupa.*', 'Pu', .) %>% gsub('.*white prepupa stage.*', 'WP', .) %>% gsub('.*strain Oregon-R S2.*', 'SO', .) %>% gsub('.*ovary female.*', 'OF', .) %>% gsub('.*whole organism larva.*', 'La', .) %>% gsub('Drosophila melanogaster ', '', .) %>% gsub('Kc167', 'Kc', .)]
		dt[organism == 'dmelanogaster' & ontology == 'ovary', stage_cell := 'OF']
		dt[organism == 'dmelanogaster', ]$stage_cell %>% table %>% sort %>% rev
		
		# for human
		dt[organism == 'human', stage_cell := ontology]
		terms_to_remove = which(table(dt[organism == 'human', ]$stage_cell) <= 5) %>% names
		dt[organism == 'human' & stage_cell %in% terms_to_remove, stage_cell := 'Other']
		dt[organism == 'human', stage_cell := stage_cell %>% gsub('neural cell', 'NC', .) %>% gsub('-| ', '', .) %>% gsub('Ishikawa', 'Ishikaw', .)]
		dt[organism == 'human', ]$stage_cell %>% table %>% sort %>% rev
		
		# for mouse
		dt[organism == 'mouse', stage_cell := ontology]
		terms_to_remove = which(table(dt[organism == 'mouse', ]$stage_cell) <= 4) %>% names
		dt[organism == 'mouse' & stage_cell %in% terms_to_remove, stage_cell := 'Other']
		dt[organism == 'mouse', stage_cell := stage_cell %>% gsub('.', '', ., fixed = T) %>% gsub('-', '', .) %>% gsub('Ishikawa', 'Ishikaw', .)]
		dt[organism == 'mouse', ]$stage_cell %>% table %>% sort %>% rev

		# creating a simplified target symbol with max 7 characters
		dt[, target_symbol_1 := target_symbol %>% gsub('-', '', .) %>% gsub('.', '', ., fixed = T) %>% gsub(')', '', ., fixed = T) %>% gsub('(', '', ., fixed = T) %>% substr(., 1, 7)]
		
		# creating final unique CHIP names
		dt[, chip_name := paste0(target_symbol_1, '_', stage_cell)]
		# most entries are unique. We add an identifyier for the duplicated entries
		dt$chip_name %>% table %>% as.integer %>% summary 
		dt[, chip_name := paste0(chip_name, '_', 1:.N), chip_name]
		dt[, chip_name := chip_name %>% gsub('_1$', '', .)]
		dt[, local_file := paste0(chip_name, '.bed.gz')]

		# reordering columns
		dt = dt[, c('chip_name', 'target_symbol', 'target_symbol_1', 'stage_cell', 'target_id', 'life_stage_age', 'ontology', 'classification', 'cell_slims', 'organ_slims', 'developmental_slims', 'system_slims', 'biosample_summary', 'organism', 'assembly', 'file', 'href', 'md5sum', 'local_file')]
		
		saveRDS(dt, 'dt_encode_chip.rds')

		ldt = split(dt, dt$chip_name)
		walk(ldt, ~saveRDS(.x, paste0(.x$chip_name, '__split_dt.rds')))


	'''
}


// Legend for the stage_cell column:

// Worm:
	// L1 YA L4 LE L3 L2 MS EE ME Da
	// 93 86 77 67 50 50 35 12  2  1

//  -  L1: L1 stage
//  -  YA: Young Adult
//  -  L4: L4 stage
//  -  LE: Late Embryo
//  -  L3: L3 stage
//  -  L2: L2 stage
//  -  MS: Mixed Stage embryo
//  -  EE: Early Embryo
//  -  ME: Mid-Embryo
//  -  Da: Dauer


// Fly:
	//  Em  PP  WT  Pu  FA  WP  MA  Kc  SO  MS  OF  La
	// 372  57  46  14  10   8   8   7   3   3   1   1
//
//  - Em: Embryo
//  - PP: Prepupa
//  - WT: Wandering third instar larval stage
//  - Pu: Puppa
//  - FA: Female Adult
//  - WP: White prepupa stage
//  - MA: Male Adult
//  - Kc: Kc167 cell line (embryonic)
//  - SO: Cell line from strain Oregon-R S2
//  - MS: Mixed Sex Adult
//  - OF: Ovary female
//  - La: larva (48 days)

// Human:
  //  K562   HepG2  HEK293 GM12878   Other    MCF7    A549      H1   liver Ishikaw
  //   432     249     193     156     137     106      55      54      35      24
  // SKNSH  HeLaS3 HEK293T   IMR90  MCF10A  HCT116    T47D GM12891      NC GM23338
  //    19      19      17      12      10      10       7       7       6       6

//  - NC: Neural Cell
//  - Other: Any ontology present less than 5 times
//  - Ishikaw: Ishikawa
//  - SKNSH: SK-N-SH

// Mouse:
   // MEL CH12LX  Other  liver   lung  heart G1EER4    G1E  ESE14
   //  49     39     36      7      5      5      5      5      5

//  - Other: Any ontology present less than 4 times
//  - CH12LX: CH12.LX
//  - G1EER4: G1E-ER4
//  -  ESE14: ES-E14




// Here are examples of json links associated to a given bed file 

// of a worm experiment:
// https://www.encodeproject.org/files/ENCFF667MVT/  => bigBed narrowPeak
// https://www.encodeproject.org/experiments/ENCSR956TMJ/?format=json
// https://www.encodeproject.org/biosamples/ENCBS172XOM/?format=json
// https://www.encodeproject.org/worm-donors/ENCDO164EGI/?format=json
// https://www.encodeproject.org/genetic-modifications/ENCGM416HQF/?format=json
// https://www.encodeproject.org/targets/npax-4-celegans/?format=json
// https://www.encodeproject.org/genes/182466/?format=json

// of a human experiment:
// https://www.encodeproject.org/files/ENCFF739GZC/?format=json
// https://www.encodeproject.org/experiments/ENCSR468DVP/?format=json
// https://www.encodeproject.org/human-donors/ENCDO000AAD/?format=json
// https://www.encodeproject.org/biosample-types/cell_line_EFO_0002067/?format=json
// https://www.encodeproject.org/biosample-types/tissue_UBERON_0000160/?format=json


// dt1 = copy(dt)
// dt1[, href1 := paste0('/files/', file, '/@@download/', file, '.bed.gz')]
// identical(dt1$href, dt1$href1) # T






process make_chip_ontology_groups {
	tag "${specie}"

	container = params.encodeexplorer

	publishDir path: "${specie}", mode: 'link', saveAs: {
    if (it.indexOf(".txt") > 0) "${it}"
		else if (it.indexOf(".tsv") > 0) "CHIP/${it}"
  }
	
	input:
		set specie, file(dt_encode_chip_rds) from Encode_chip_metadata

	output:
		file("*")

	shell:
	'''

		#!/usr/bin/env Rscript
		
		library(data.table)
		library(magrittr)
		library(purrr)
		
		specie = '!{specie}'
		dt1 = readRDS('!{dt_encode_chip_rds}')
		source('!{params.cactus_dir}/software/get_data/bin/make_chip_ontology_groups_functions.R')
		source('!{params.cactus_dir}/software/get_data/bin/get_tab.R')
		
		dt = copy(dt1)
		dt$local_file %<>% gsub('.gz', '', ., fixed = T)
		
		l_tissue_chip_names        = get_ontology_chip_names(dt, 'organ_slims', 'organ')
		l_cell_type_chip_names     = get_ontology_chip_names(dt, 'cell_slims', 'cell_type')
		l_cell_line_chip_names     = get_ontology_chip_names(dt, 'ontology', 'cell_line')
		l_system_chip_names        = get_ontology_chip_names(dt, 'system_slims', 'system')
		l_developmental_chip_names = get_ontology_chip_names(dt, 'developmental_slims', 'development')
		l_all_chip_names           = list(all = dt$local_file)
		
		l_ontology_chip_names = c(l_all_chip_names, l_tissue_chip_names, l_cell_type_chip_names, l_cell_line_chip_names, l_system_chip_names, l_developmental_chip_names)
		ontology_order = map_int(l_ontology_chip_names, length) %>% sort %>% rev %>% names
		l_ontology_chip_names %<>% .[ontology_order]
		
		dt_n_chip_by_ontology = data.table(ontology = names(l_ontology_chip_names), number_of_chip = map_int(l_ontology_chip_names, length))
		
		sink('available_chip_ontology_groups.txt')
			 print(dt_n_chip_by_ontology)
		sink()
		
		tb = tibble::enframe(l_ontology_chip_names) %>% {tidyr::unnest(., cols = c(value))}
		write.table(tb, file = 'chip_ontology_groups.tsv', col.names = F, row.names = F, quote = F, sep = get_tab())
		
		
	'''
}



process get_chip_data {
tag "${specie}"

container = params.encodeexplorer

publishDir path: "${specie}/CHIP", mode: 'link'

input:
	set specie, file(dt_rds) from Encode_chip_split_dt
	
output:
	file("*.bed")

shell:
'''

	#!/usr/bin/env Rscript

	library(data.table)
	dt = readRDS('!{dt_rds}')
	url_encode = 'https://www.encodeproject.org' 
	
	local_file = dt$local_file
	download.file(url = paste0(url_encode, dt$href), quiet = T, destfile = local_file, method = 'curl', extra = '-L' )
	md5sum_downloaded_files = tools::md5sum(local_file)
	if(dt$md5sum != md5sum_downloaded_files) stop(paste('md5 sums not equal for file', local_file))
	system(paste('gunzip', local_file))
			
'''
}


Assembly_ENCODE_chromatin_state
	.filter{ specie, assembly -> specie in ['human', 'mouse']}
	.set{ Assembly_ENCODE_chromatin_state_1 }


process get_encode_chromatin_state_metadata {
	tag "${specie}"

	container = params.encodeexplorer

	publishDir path: "${specie}", mode: 'link', pattern: "*.csv"

	input:
		set specie, assembly from Assembly_ENCODE_chromatin_state_1

	output:
		set specie, assembly, file("*.rds") into Encode_chromatin_state_split_dt mode flatten
		file('encode_chromatin_states.csv')

	shell:
	'''

		#!/usr/bin/env Rscript
		
		source('!{params.cactus_dir}/software/get_data/bin/encode_chip_functions.R')
		assembly = '!{assembly}'
		
		library(data.table)
		library(magrittr)
		library(purrr)
		
		url_encode = 'https://www.encodeproject.org' 
		url_search = paste0(url_encode, '/search/?')
		url_append = '&frame=object&format=json&limit=all'
		
    
    ## fetching tables
    
    df_ecsm_files = get_encode_df(paste0('type=File&annotation_type=chromatin+state&status=released&file_format=bed&assembly=', assembly)) %T>% pnrow
    df_biosample_types = get_encode_df('type=BiosampleType') %T>% pnrow
    df_annotations = get_encode_df('type=Annotation&annotation_type=chromatin+state&status=released') %T>% pnrow
    
    
    ## formatting and merging tables
    
		dt_biosample_types = data.table(
		  df_biosample_types[, c('biosample_types', 'classification', 'term_name')],
		  cell_slims  = collapse_slims(df_biosample_types$cell_slims),
		  organ_slims = collapse_slims(df_biosample_types$organ_slims),
		  developmental_slims = collapse_slims(df_biosample_types$developmental_slims),
		  system_slims = collapse_slims(df_biosample_types$system_slims)
		)
    colnames(dt_biosample_types)[1] = 'biosample_ontology'
  
    dt_annotation = data.table(df_annotations[, c('annotations', 'aliases', 'description', 'organism', 'relevant_timepoint', 'relevant_timepoint_units', 'relevant_life_stage')])
    dt_annotation$aliases %<>% map_chr(~ifelse(is.null(.x), NA, .x))
    dt_annotation[, alias := map_chr(aliases, ~ifelse(is.null(.x), NA, .x))]
    dt_annotation[, aliases := NULL]
    dt_annotation$organism %<>% gsub('/organisms/', '', .) %>% gsub('/', '', .)
    colnames(dt_annotation) %<>% gsub('relevant_', '', .)
    
    dt = data.table(df_ecsm_files[, c('accession', 'files', 'assembly', 'md5sum', 'href', 'dataset', 'biosample_ontology', 'file_size')])
    colnames(dt)[colnames(dt) == 'dataset'] = 'annotations'
    dt = dt_biosample_types[dt, , on = 'biosample_ontology']
    dt = dt_annotation[dt, , on = 'annotations']
		dt[, local_file := paste0(accession, '.bed.gz')]    
	
		dt = dt[, c('accession', 'assembly', 'description', 'alias', 'timepoint', 'timepoint_units', 'life_stage', 'term_name', 'classification', 'cell_slims', 'organ_slims', 'developmental_slims', 'system_slims', 'organism', 'local_file', 'href', 'md5sum', 'file_size')]
		
		if(assembly == 'mm10') dt = dt[grep('18states', dt$alias)]
		
		dt1 = dt[, 1:(ncol(dt)-4)]
		write.csv(dt1, 'encode_chromatin_states.csv')
		
		ldt = split(dt, dt$accession)
		walk(ldt, ~saveRDS(.x, paste0(.x$accession, '.rds')))
		
		
	'''
}

	



// dt[local_file == 'ENCFF686HVR.bed.gz']

// search link
// https://www.encodeproject.org/search/?type=File&annotation_type=chromatin+state&assembly=GRCh38&assembly=mm10&status=released&file_format=bed

// query results when using all species
// df_ecsm_files = get_encode_df('type=File&annotation_type=chromatin+state&status=released&file_format=bed&assembly=GRCh38&assembly=mm10') %T>% pnrow # 962
// df_biosample_types = get_encode_df('type=BiosampleType') %T>% pnrow # 936
// df_annotations = get_encode_df('type=Annotation&annotation_type=chromatin+state&status=released') %T>% pnrow # 1251

// >   df_ecsm_files = get_encode_df(paste0('type=File&annotation_type=chromatin+state&status=released&file_format=bed&assembly=', assembly)) %T>% pnrow
// Error in get_encode_df(paste0("type=File&annotation_type=chromatin+state&status=released&file_format=bed&assembly=",  :
//   url doesnt exist
// > paste0('type=File&annotation_type=chromatin+state&status=released&file_format=bed&assembly=', assembly)
// [1] "type=File&annotation_type=chromatin+state&status=released&file_format=bed&assembly=hg38"
// >


	// ## donwloading and unzipping the data
	// sapply(1:nrow(dt), function(c1) download.file(url = paste0(url_encode, dt$href[c1]), quiet = T, destfile = dt$local_file[c1], method = 'curl', extra = '-L' ))
	// dt[, md5sum_downloaded_files := tools::md5sum(local_file)]
	// if(any(dt$md5sum != dt$md5sum_downloaded_files)) stop('not all md5 sums are equal')
	// system('for z in *.gz; do gunzip "$z"; done')


// All experiments have been made in a 15 and 18 state model
// However, the 15 state model is much bigger
// > dt[alias %in% c("zhiping-weng:chromhmm-8marks-15states-mm10-heart-12.5", "zhiping-weng:chromhmm-10marks-18states-mm10-heart-12.5")][, c('description', 'file_size')]
//                                                description file_size
// 1: ChromHMM 15-state model for heart (embryonic 12.5 days) 121270263
// 2: ChromHMM 18 state model for heart (embryonic 12.5 days)   8684869

// => using the 18 state model

// dt[grep('18states', dt$alias)]$file_size %>% max #   8684869
// dt[grep('15states', dt$alias)]$file_size %>% min # 109950661







process get_encode_chromatin_state_data {
	tag "${specie}, ${dt_rds}"

	container = params.encodeexplorer

	publishDir path: "${specie}/chromatin_states", mode: 'link'

	input:
		set specie, assembly, file(dt_rds) from Encode_chromatin_state_split_dt

	output:
		file("*")
		set specie, assembly, file("*.bed") into Encode_bed_for_liftover

	shell:
	'''

		#!/usr/bin/env Rscript
		
		library(magrittr)
		library(data.table)
		dt = readRDS('!{dt_rds}')
		source('!{params.cactus_dir}/software/get_data/bin/get_tab.R')
		
		# downloading the bed file
		url_encode = 'https://www.encodeproject.org'
		local_file = dt$local_file
		download.file(url = paste0(url_encode, dt$href), quiet = T, destfile = local_file, method = 'curl', extra = '-L' )
		md5sum_downloaded_files = tools::md5sum(local_file)
		if(dt$md5sum != md5sum_downloaded_files) stop(paste('md5 sums not equal for file', local_file))
		system(paste('gunzip', local_file))
		
		
		# splitting each state in a separate file
		accession = dt$accession
		bed_file = paste0(accession, '.bed')
		df = read.table(bed_file, sep = get_tab(), stringsAsFactors = F)
		dir.create(accession)
		df$V4 %<>% gsub('/', '_', .)
		
		v_states = sort(unique(df[,4]))
		for(state in v_states){
			df1 = df[df$V4 == state, ]
			write.table(df1, paste0(accession, '/', state, '.bed'), col.names = F, row.names = F, quote = F, sep = get_tab())
		}
		
	'''
}



liftover_files_to_get = Channel
	.from( 	
		[
			['human', 'hg19', 'hg38'],
			[ 'worm', 'ce10', 'ce11'],
			[  'fly',  'dm3', 'dm6'],
			['mouse',  'mm10', 'mm39']
		]
	)



process get_liftover_files {
	tag "${specie}/${original_assembly} to ${target_assembly}"
	
	publishDir path: "util/liftover_files", mode: 'link'

	container = params.liftover

	input:
		set specie, original_assembly, target_assembly from liftover_files_to_get

	output:
		set specie, original_assembly, target_assembly, file("*.over.chain.gz") into liftover_files

	shell:
	'''
					
			original_assembly="!{original_assembly}"
			target_assembly="!{target_assembly}"
			
			target_assembly_1=$(echo $target_assembly | awk '{ printf ("%s%s", toupper (substr ($0, 1, 1)), substr ($0, 2)); }')
			liftover_name="${original_assembly}To${target_assembly_1}"
			liftover_file="${liftover_name}.over.chain.gz"

			url="https://hgdownload.cse.ucsc.edu/goldenPath/${original_assembly}/liftOver"
			wget $url/$liftover_file
			
	'''
}

//// md5sum is missing for ce10ToCe11 => so we unfortunately need to skip md5sum checking of files
// (command not needed: wget $url/$liftover_file $url/md5sum.txt)
// https://hgdownload.cse.ucsc.edu/goldenPath/ce10/liftOver/md5sum.txt
// However, I do still manually check that the file of the file we get is the same as the one indicated online. That is:
// ce10ToCe11.over.chain.gz     2015-06-23 15:39  3.2K  
// mm10ToMm39.over.chain.gz     2020-07-30 14:50   24K  
// hg19ToHg38.over.chain.gz     2013-12-31 23:08  222K
// dm3ToDm6.over.chain.gz     2014-08-28 14:29  1.8M  

// ls -sh util/liftover_files


Encode_bed_for_liftover
	.join(liftover_files, by:[0,1])
	.view()


// process convert_bed_coordinates_to_latest_assembly {
// 	tag "${specie}/${out_folder}"
// 
// 	publishDir path: "${specie}/{out_folder}", mode: 'link'
// 
// 	container = params.liftover
// 
// 	input:
// 		set specie, out_folder, original_assembly, target_assembly, file(bed_file) from channel_xx
// 
// 	output:
// 		file("*")
// 
// 	shell:
// 	'''
// 
// 			original_assembly="!{original_assembly}"
// 			liftover_name="!{liftover_name}"
// 			bed_name="!{bed_name}"
// 			bed_file="!{bed_file}"
// 
// 			liftover_file="${liftover_name}.over.chain.gz"
// 			bed_file_lifted="${bed_name}_lifted.bed"
// 
// 			liftover_path="https://hgdownload.cse.ucsc.edu/goldenPath/${original_assembly}/liftOver"
// 			wget $liftover_path/$liftover_file
// 			liftOver $bed_file  $liftover_file $bed_file_lifted unmapped.bed
// 			mkdir $bed_name
// 			awk  -v FOLDER="$bed_name" ' { print > FOLDER"/"$4".bed" }' $bed_file_lifted
// 
// 			rm $bed_name/17_Unmap.bed
// 
// 	'''
// }






HiHMM_chromatin_states
	.combine(HiHMM_liftover, by:0)
	.set{ HiHMM_chromatin_states_1 }
	


process get_hihmm_chromatin_state_data_part_1 {
	tag "${specie}"
	
	// No containers for the process as of now since I couldn't get one that works (and that is downloaded directly from within nextflow. This may be an issue with the outdated version of Singularity that is installed on my server)
	// container = "kernsuite-debian/singularity-container"
	// container = "kernsuite-debian/singularity-container"
	// container = "${params.depot_galaxy}/ucsc_tools:357--0"
	// container = params.bioconductor
	// container = 'ubuntu/library/ubuntu/22.10'
	// container = 'debian/stable-slim.sif'

	input:
		set specie, bed_name, original_assembly, liftover_name from HiHMM_chromatin_states_1

	output:
		set specie, bed_name, file('*.bed'), original_assembly, liftover_name into HiHMM_chromatin_states_2
		

	shell:
	'''
			
			bed_name="!{bed_name}"
			specie="!{specie}"
			
			bed_file="${bed_name}.bed"
			bed_file_lifted="${bed_name}_lifted.bed"
			
			hiHMM_path="http://compbio.med.harvard.edu/modencode/webpage/hihmm"
			
			wget $hiHMM_path/$bed_file
			
			if [ $specie = 'worm' ]; then
				gawk -i inplace '{print "chr"$0}' $bed_file
			fi

	'''
}


process get_hihmm_chromatin_state_data_part_2 {
	tag "${specie}"
	
	publishDir path: "${specie}/chromatin_states", mode: 'link'

	container = params.liftover

	input:
		set specie, bed_name, file(bed_file), original_assembly, liftover_name from HiHMM_chromatin_states_2

	output:
		file("${bed_name}/*")

	shell:
	'''
					
			original_assembly="!{original_assembly}"
			liftover_name="!{liftover_name}"
			bed_name="!{bed_name}"
			bed_file="!{bed_file}"

			liftover_file="${liftover_name}.over.chain.gz"
			bed_file_lifted="${bed_name}_lifted.bed"

			liftover_path="https://hgdownload.cse.ucsc.edu/goldenPath/${original_assembly}/liftOver"
			wget $liftover_path/$liftover_file
			liftOver $bed_file  $liftover_file $bed_file_lifted unmapped.bed
			mkdir $bed_name
			awk  -v FOLDER="$bed_name" ' { print > FOLDER"/"$4".bed" }' $bed_file_lifted

			rm $bed_name/17_Unmap.bed
			
	'''
}


// this command fails in the container but not outside the container:
// wget "http://compbio.med.harvard.edu/modencode/webpage/hihmm/iHMM.M1K16.human_H1.bed"
// wget: bad header line:     XeuOGalu: ptQ; path=/; Max-Age=900








// process get_hihmm_chromatin_state_data {
// 	publishDir path: "${specie}/chromatin_states", mode: 'link'
// 
// 	container = params.liftover
// 
// 	input:
// 		set specie, bed_name, original_assembly, liftover_name from HiHMM_chromatin_states_1
// 
// 	output:
// 		file(bed_name)
// 
// 	shell:
// 	'''
// 
// 			liftover_name="!{liftover_name}"
// 			bed_name="!{bed_name}"
// 			original_assembly="!{original_assembly}"
// 
// 			liftover_file="${liftover_name}.over.chain.gz"
// 			bed_file="${bed_name}.bed"
// 			bed_file_lifted="${bed_name}_lifted.bed"
// 
// 			liftover_path="https://hgdownload.cse.ucsc.edu/goldenPath/${original_assembly}/liftOver"
// 			hiHMM_path="http://compbio.med.harvard.edu/modencode/webpage/hihmm"
// 			wget $liftover_path/$liftover_file
// 			wget $hiHMM_path/$bed_file
// 			gawk -i inplace '{print "chr"$0}' $bed_file
// 			liftOver $bed_file  $liftover_file $bed_file_lifted unmapped.bed 
// 			mkdir $bed_name
// 			awk ' { print >  $bed_name"/"$4".bed" }' $bed_file_lifted
// 	'''
// }
// 
// // this command fails in the container but not outside the container:
// // wget "http://compbio.med.harvard.edu/modencode/webpage/hihmm/iHMM.M1K16.human_H1.bed"
// // wget: bad header line:     XeuOGalu: ptQ; path=/; Max-Age=900
// => all containers I tried fail with this command




// process get_bowtie2_indexes_for_checking_contaminations {
// 	tag "${specie}"
// 
// 	// container = "	kernsuite-debian/singularity-container"
// 	container = "aws-cli_latest.sif"
// 	// container = params.bioconductor
// 
// 	input:
// 		// set specie, bed_name, original_assembly, liftover_name from HiHMM_chromatin_states_1
// 
// 	output:
// 		// set specie, bed_name, file('*.bed'), original_assembly, liftover_name into HiHMM_chromatin_states_2
// 
// 
// 	shell:
// 	'''
// 
// 			aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/ ./my_refs/
// 
// 
// 	'''
// }
// 
// // https://ewels.github.io/AWS-iGenomes/
// // https://support.illumina.com/sequencing/sequencing_software/igenome.html
// note: I pulled the singularity image separately prior to run the script in this case

// => E coli OP50 strain is not present on iGenomes. Need to get the file another way.




process get_contamination_bowtie2 {
	tag "${specie}"

	container = params.bowtie2_samtools

	publishDir path: "human/bowtie2_indexes_conta", mode: 'link'
	publishDir path: "mouse/bowtie2_indexes_conta", mode: 'link'
	publishDir path:   "fly/bowtie2_indexes_conta", mode: 'link'
	publishDir path:  "worm/bowtie2_indexes_conta", mode: 'link'

	input:

	output:
		set file('README.txt'), file('*.bt2')

	shell:
	'''
			
		wget -O genome_contamination.fa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/496/595/GCF_009496595.1_ASM949659v1/GCF_009496595.1_ASM949659v1_genomic.fna.gz
		echo "fasta file: GCF_009496595.1_ASM949659v1_genomic.fna.gz" > README.txt
		
		gunzip genome_contamination.fa.gz
		
		bowtie2-build --threads !{number_of_cores} genome_contamination.fa genome_contamination
			
	'''
}

// There are 3 assemblies for OP50 (one from 2009 and 2 from 2019):
// https://www.ncbi.nlm.nih.gov/assembly/?term=Escherichia+coli+OP50+%28E.+coli%29
// The latest one is this one (ASM949659v1):
// https://www.ncbi.nlm.nih.gov/assembly/GCF_009496595.1 (file: GCF_009496595.1_ASM949659v1_genomic.fna.gz)
// (previously I was using the one from 2009: ASM17681v1, with the file GCA_000176815.1_ASM17681v1_genomic.fna)



process get_fasta_and_gff {
	tag "${specie}"

	container = params.samtools_bedtools_perl

	publishDir path: "${specie}/genome", mode: 'link', saveAs: {
    if (it.indexOf(".gff3") > 0) "annotation/${it}"
		else if (it.indexOf(".fa") > 0) "sequence/${it}"
    else if (it.indexOf(".txt") > 0) "${it}"
  }

	input:
		set specie, specie_long, genome, assembly, release from Assembly_Ensembl_1.fasta

	output:
		set specie, file('annotation.gff3'), file('genome.fa') into (Genome_and_annotation, Genome_and_annotation_1, Genome_and_annotation_2, Genome_and_annotation_3)
		set specie, file('annotation.gff3') into ( Gff3_annotation_for_filtering, Gff3_annotation_for_bed_files)

	shell:
	'''
	
	source !{params.cactus_dir}/software/get_data/bin/check_checksum_ensembl.sh
	specie='!{specie}'
	specie_long='!{specie_long}'
	genome='!{genome}'
	assembly='!{assembly}'
	release='!{release}'
	
	base_url=ftp://ftp.ensembl.org/pub/release-$release
	
	gff3_file=${specie_long^}.$genome.$release.gff3.gz
	fasta_file=${specie_long^}.$genome.dna_sm.$assembly.fa.gz
	
	url_anno=$base_url/gff3/$specie_long
	wget -O annotation.gff3.gz $url_anno/$gff3_file 
	wget -O checksums_anno $url_anno/CHECKSUMS
	check_checksum_and_gunzip annotation.gff3.gz $gff3_file checksums_anno
	
	url_fasta=$base_url/fasta/$specie_long/dna
	wget -O genome.fa.gz $url_fasta/$fasta_file 
	wget -O checksums_fasta $url_fasta/CHECKSUMS
	check_checksum_and_gunzip genome.fa.gz $fasta_file checksums_fasta
		
	'''
}

// downloaded_file=annotation.gff3.gz
// downloaded_file_name=$gff3_file
// checksums_file=checksums_anno


// https://genomespot.blogspot.com/2015/06/mapping-ngs-data-which-genome-version.html
// Which options are available? My recommendation at this point is Ensembl, for a number of reasons... it is clear to see what genome build and version just from the file names
// Repeat masking? The short answer is that for the purpose of mapping NGS reads, you don't want to use a RepeatMasked genome The files with "sm" in the name are "soft masked" which means instead of the repeat regions being given Ns, they are converted to lower-case. These are OK to use because most major aligners recognise lower-case letters as valid bases.
//  Primary assembly vs Top level? The short answer is to choose the "Primary assembly" (i.e. Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz) as it does not contain the haplotype information. The top level file contains additional sequences that are relatively common variants to the reference. 
 // Should I use the latest version? Short answer: Generally, yes. But of course this depends. in some projects we are still using the older version of the human genome (GRCh37/hg19) because the aim of those projects is to integrate with a bunch of ENCODE data which is not yet available for the new genome build.
  // => in our case we don't use the latest mouse genome version specifically for this reason

// https://bioinformatics.stackexchange.com/questions/540/what-ensembl-genome-version-should-i-use-for-alignments-e-g-toplevel-fa-vs-p

// Examples of link:
// http://ftp.ensembl.org/pub/release-102/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz
// http://ftp.ensembl.org/pub/release-102/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.102.gff3.gz
// wget -O checksums_anno.tsv http://ftp.ensembl.org/pub/release-102/gff3/drosophila_melanogaster/CHECKSUMS

// => primary_assembly is not available for worm an fly
// This means that top_level = primary_assembly in these species. See here:
// https://www.biostars.org/p/398167/
// "If the primary assembly file is not present, that indicates that there are no haplotype/patch regions, and the 'toplevel' file is equivalent."

// checksums for ensembl genomes
// https://stackoverflow.com/questions/14569281/checksum-and-md5-not-the-same-thing





process generating_bowtie2_indexes {
  tag "${specie}"

  container = params.bowtie2_samtools

  publishDir path: "${"${specie}/genome/sequence"}/bowtie2_indexes", mode: 'link'

  input:
    set specie, file(gff3), file(fasta) from Genome_and_annotation_1

  output:
		file("*.bt2")

  shell:
  '''

		bowtie2-build --threads !{number_of_cores} !{fasta} genome


  '''

}


process indexing_genomes {
  tag "${specie}"

  container = params.bowtie2_samtools

  publishDir path: "${"${specie}/genome/sequence"}", mode: 'link'

  input:
    set specie, file(gff3), file(fasta) from Genome_and_annotation_2

  output:
		set specie, file('genome.fa.fai') into Fasta_fai

  shell:
  '''

		samtools faidx !{fasta}


  '''

}

Genome_and_annotation_3
	.join(Fasta_fai)
	.set{ Fasta_fai_gff3 }





process getting_transcriptome {
  tag "${specie}"

	container = params.gffread

	publishDir path: "${specie}/genome/sequence", mode: 'link'
		
  input:
    set specie, file(gff3), file(fasta), file(fasta_indexes) from Fasta_fai_gff3

  output:
		set specie, file('transcriptome.fa') into Transcriptome_for_building_kallisto_indexes
		
  shell:
  '''
			
      gffread -C !{gff3} -g !{fasta} -w transcriptome.fa
			
  '''

}




process filtering_annotation_file {
  tag "${specie}"

	container = params.ubuntu

	publishDir path: "${specie}/genome/annotation/filtered", mode: 'link'
		
  input:
    set specie, file(gff3) from Gff3_annotation_for_filtering

  output:
		set specie, file('chromosomes_sizes.txt'), file(gff3), file('protein_coding_genes.gff3'), file('annotation_filtered_regions_and_pseudogenes.gff3') into GFF3_filtered_for_R
		file("*.txt")
		
  shell:
  '''
			
			gff3=!{gff3}
			source !{params.cactus_dir}/software/get_data/bin/get_tab.sh
		
			tab=$(get_tab)
		
			gff3_filtered=annotation_filtered_regions_and_pseudogenes.gff3

			# savings statistics on the unfiltered annotation file
			cut -f1 ${gff3} | sort | uniq -c | grep -v "#" | sort -nr > nb_of_feature_by_region__raw_gff3_file.txt
			awk -v FS="$tab" -v OFS="$tab" '{if($3 == "gene"){print $1}}' ${gff3} | sort | uniq -c | sort -nr	> nb_of_genes_by_region__raw_gff3_file.txt 
			cut -f3 ${gff3} | sort | uniq -c | sort -nr | grep -v "#" > nb_of_feature_by_feature_type__raw_gff3_file.txt
			
			# keeping only regions (contigs and chromosomes) with at least 5 genes, excluding the mitochondrial chromosome and removing all genes that are not encoding for proteins
			regions_to_keep=$(awk -v OFS="$tab" '{if($1 > 4) print $2}' nb_of_genes_by_region__raw_gff3_file.txt | sort | paste -sd ' ')
			awk -v FS="$tab" -v OFS="$tab" -v rtk="$regions_to_keep" 'BEGIN{split(rtk, rtk1, " ")} {for (c1 in rtk1) {if($1 == rtk1[c1] && $1 != "mitochondrion_genome" && $1 != "MT" && $1 != "MtDNA" && ($3 == "gene" || $3 !~ "gene")) print}}' ${gff3} > ${gff3_filtered}

			# savings statistics on the filtered annotation file
			cut -f1 ${gff3_filtered} | sort | uniq -c | grep -v "#" | sort -nr > nb_of_feature_by_region__raw_gff3_file.txt
			awk -v FS="$tab" -v OFS="$tab" '{if($3 == "gene"){print $1}}' ${gff3_filtered} | sort | uniq -c | sort -nr	> nb_of_genes_by_region__filtered_gff3_file.txt 
			cut -f3 ${gff3_filtered} | sort | uniq -c | sort -nr | grep -v "#" > nb_of_feature_by_feature_type__filtered_gff3_file.txt
			
			# exporting chromosomes sizes
			awk -v FS="$tab" -v OFS="$tab" '{if ($3 == "region" || $3 == "chromosome" || $3 == "scaffold") {print $1, $5}}' ${gff3_filtered} | sort -k3nr > chromosomes_sizes.txt
			
			# exporting protein coding genes
			grep -i protein_coding ${gff3_filtered} | awk -v FS="$tab" -v OFS="$tab" '{if($3 == "gene") print $0}' - > protein_coding_genes.gff3
			
			
			
  '''

}


// note: the important thing is to keep all contigs when mapping reads to prevent misalignement of reads to wrong location. However, discarding reads mapped to contigs after mapping is a matter of taste and depends on the downstream application. In our case, we need to avoid small contigs to avoid erroneous annotations of peaks to the closest gene. So we discard small contigs (current threshold: we keep only contigs with 5 genes or more).

// singularity pull ubuntu:22.04

// not sure how the annotation will work in chipseeker. Need to check manually. Maybe we need to remove more features. In this case this command will do it:
// awk -v FS="\t" -v OFS="\t" -v rtk="$regions_to_keep" 'BEGIN{split(rtk, rtk1, " ")} {for (c1 in rtk1) {if($1 == rtk1[c1] && $1 != "mitochondrion_genome" && $1 != "MT" && $1 != "MtDNA" && $3 !~ "ncRNA" && $3 !~ "pseudogen" && $3 !~ "transposable_element") print $1}}' ${gff3} > ${gff3_filtered}


// awk '{if($3 == "gene"){print $1}}' mouse/genome/annotation/annotation.gff3 | sort | uniq -c | sort -nr | awk '{if ($1 >= 5) {print $2}}' - | wc -l # 30
// awk '{if($3 == "gene"){print $1}}' human/genome/annotation/annotation.gff3 | sort | uniq -c | sort -nr | awk '{if ($1 >= 5) {print $2}}' - | wc -l # 25
// 
// awk '{if($3 == "gene"){print $1}}' human/genome/annotation/annotation.gff3 | sort | uniq -c | sort -nr | awk '{if ($1 >= 5) {print}}' -
// awk '{if($3 == "gene"){print $1}}' mouse/genome/annotation/annotation.gff3 | sort | uniq -c | sort -nr | awk '{if ($1 >= 5) {print}}' -
//  awk '{if($3 == "gene"){print $1}}' fly/genome/annotation/annotation.gff3 | sort | uniq -c | sort -nr 

//// outputs:
// - protein_coding_genes.gff3 is used to create the gene metadata table for all analysis
// - anno_without_pseudogenes.gff3 is used for creating a txdb database without non coding transcripts, for the annotation of ATAC seq peaks

//// mitochondrial chromosomes:
// worm  |         fly          | mouse  | human  
// MtDNA | mitochondrion_genome |   MT   |  MT

// scaffolds
// fly: 2110000*
// human: GL000*, KI270*

// source !{params.cactus_dir}/software/get_data/bin/get_tab.sh

process get_bed_files_of_annotated_regions {
	tag "${specie}"

	container = params.bedops

	publishDir path: "${specie}/genome/annotation/bed_regions", mode: 'link'

	input:
		set specie, file(gff3), file(chr_size) from Gff3_annotation_for_bed_files

	output:
		file('*.bed')

	shell:
	'''
	
			gff3=!{gff3}

			# Filter out predicted "biological regions" (= TSS, CpG islands, fosmids, ect and chromosomes), and sort the gff file
			awk -v FS="\t" -v OFS="\t" '{if (substr($1, 1, 1) != "#" && $3 != "biological_region" && $3 != "region" && $3 != "scaffold") print}' $gff3 | sort -k1,1 -k4,4n > annotation_clean.gff3
			gff2bed < annotation_clean.gff3 > annotation_clean.bed
			
			# Get already defined regions: exons, genes
			awk -v FS="\t" -v OFS="\t" '{if ($3 == "exon") print $1, $4-1, $5}' annotation_clean.gff3 > exons.bed
			awk -v FS="\t" -v OFS="\t" '{if ($3 == "gene") print $1, $4-1, $5}' annotation_clean.gff3 > genes.bed

			# Get chromosomes sizes
			awk -v FS="\t" -v OFS="\t" '{if ($3 == "region" || $3 == "chromosome" || $3 == "scaffold") {print $1, $5}}' ${gff3} | sort -k1,1 > chromosomes_sizes.txt
			awk -v FS="\t" -v OFS="\t" '{print $1, "0", $2}' chromosomes_sizes.txt | sort -k1,1 -k2,2n > chromosomes_sizes.bed
			
			# Get intergenic regions
			bedops --difference chromosomes_sizes.bed annotation_clean.bed > intergenic.bed
			
			# Get introns
			bedops --difference chromosomes_sizes.bed <(cat exons.bed intergenic.bed | sort -k1,1 -k2,2n) > introns.bed
			
			# Get promoters
			gff2bed < <(awk -v FS="\t" -v OFS="\t" '($3 == "gene") {print}' annotation_clean.gff3) > genes_detailed.bed
			awk -v FS="\t" -v OFS="\t" '($6 == "+"){ print $1, ($2 - 1), $2, $4, $5, $6; }' genes_detailed.bed | bedops --range -1500:500 --everything - > promoters_forward.bed
			awk -v FS="\t" -v OFS="\t" '($6 == "-"){ print $1, $3, ($3 + 1), $4, $5, $6; }' genes_detailed.bed | bedops --range -1500:500 --everything - > promoters_reverse.bed
			bedops --everything promoters_forward.bed promoters_reverse.bed | cut -f1-3 > promoters.bed
			
			rm genes_detailed.bed promoters_forward.bed promoters_reverse.bed annotation_clean.bed chromosomes_sizes.bed
		
	'''
}


// Not needed finally?
// # Filter out predicted "biological regions" (= TSS, CpG islands and chromosomes),  and sort the gff file
// awk -v FS="\t" -v OFS="\t" '{if (substr($1, 1, 1) != "#" && $3 != "biological_region" && $3 != "region" && $3 != "scaffold") print}' $gff3 | sort -k1,1 -k4,4n > annotation_clean.gff3
// gff2bed < annotation_clean.gff3 > annotation_clean.bed


 // awk '{if ($1 !~ "##" && $1 !~ "Scaffold"  && $1 !~ "2110000") print}' annotation.gff3  |  cut -f1 | sort | uniq

// tab=$(get_tab)
// awk -v my_var="$tab" ' BEGIN {OFS=my_var} {print $1, $2}' test_anno_2.bed

// not that to get intergenic there commands are equivalent"
// bedops --difference chrom_size.bed target_file.bed
// bedtools --complement target_file.bed chrom_size.txt

// cut -f8 annotation_clean.bed | sort | uniq | paste -sd "  " - 
// cut -f1 exons.bed | uniq	

// Here are the lines used to create the get_tab.sh script (creating within bash prevent the issue with windows line-breaks)
// cat > get_data/bin/get_tab.sh << EOL
// #!/usr/bin/env bash
// 
// get_tab () { echo "\t" ; }
// 
// EOL


// cat '!{params.cactus_dir}/software/get_data/bin/get_tab.sh'


// https://www.biostars.org/p/112251/#314840 => get intergenic, exons and introns .bed
// https://www.biostars.org/p/360842/#9531161 // => get promoter.bed

// other links
// https://github.com/TriLab-bioinf/Script_wrappers/blob/main/get_intergenic_bed_from_gtf.sh
// https://davetang.org/muse/2013/01/18/defining-genomic-regions/
// https://bedparse.readthedocs.io/en/stable/Tutorial.html // => could not make it work

// grep -v "sequence-region" annotation.bed | cut -f 8 | sort | uniq | paste -sd "  " - 
// biological_region CDS C_gene_segment chromosome D_gene_segment exon five_prime_UTR gene J_gene_segment lnc_RNA miRNA mRNA ncRNA ncRNA_gene pseudogene pseudogenic_transcript rRNA scaffold scRNA snoRNA snRNA three_prime_UTR tRNA unconfirmed_transcript V_gene_segment

// commands not needed:
// # convert gff3 to bed
// gff2bed < annotation.gff3 > annotation.bed
//
// # make a bed file of chromosome sizes
// awk 'OFS="\t" {print $1, "0", $2}' chromosome_size.txt | sort -k1,1 -k2,2n > chromosome_size.bed
// awk -v my_var="$tab" 'BEGIN {OFS=my_var} {if ($3 == "five_prime_UTR")  print $1, $4-1, $5}' annotation_clean.gff3 > five_prime_UTR.bed
// awk -v my_var="$tab" 'BEGIN {OFS=my_var} {if ($3 == "three_prime_UTR") print $1, $4-1, $5}' annotation_clean.gff3 > three_prime_UTR.bed
// # Get all regions
// bedops --everything promoters.bed exons.bed introns.bed intergenic.bed genic_regions.bed > all_regions.bed

// The gff3 file needs to be sorted:
// Error: Sorted input specified, but the file annotation.gff3 has the following out of order record
// 1       havana  pseudogenic_transcript  12010   13670   .       +       .       ID=transcript:ENST00000450305;Parent=gene:ENSG00000223972;Name=DDX11L1-201;biotype=transcribed_unprocessed_pseudogene;tag=basic;transcript_id=ENST00000450305;transcript_support_level=NA;version=2


// ## Exploring a bit the structure of the gff file
// 
// awk '{ if(substr($1, 1, 1) != "#") {print $0}}' annotation.gff3 | head
// 
// awk '{ if($3 == "gene") {print $0}}' annotation.gff3 | head -1
// // 1       ensembl_havana  gene    65419   71585   .       +       .       ID=gene:ENSG00000186092;Name=OR4F5;biotype=protein_coding;description=olfactory receptor family 4 subfamily F member 5 [Source:HGNC Symbol%3BAcc:HGNC:14825];gene_id=ENSG00000186092;logic_name=ensembl_havana_gene_homo_sapiens;version=6
// 
// awk '{ if($3 == "exon") {print $0}}' annotation.gff3 | head -1
// // 1       havana  exon    11869   12227   .       +       .       Parent=transcript:ENST00000456328;Name=ENSE00002234944;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00002234944;rank=1;version=1
// 
// awk '{ if($3 == "mRNA") {print $0}}' annotation.gff3 | head -1
// // 1       havana  mRNA    65419   71585   .       +       .       ID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-202;biotype=protein_coding;tag=basic;transcript_id=ENST00000641515;version=2
// 
// 
// awk '{ if($3 == "exon") {print $0}}' annotation.gff3 | head -10
// // 1       havana  exon    11869   12227   .       +       .       Parent=transcript:ENST00000456328;Name=ENSE00002234944;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00002234944;rank=1;version=1
// // 1       havana  exon    12010   12057   .       +       .       Parent=transcript:ENST00000450305;Name=ENSE00001948541;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001948541;rank=1;version=1
// // 1       havana  exon    14404   14501   .       -       .       Parent=transcript:ENST00000488147;Name=ENSE00001843071;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00001843071;
// 
// awk '{ if($9 ~ /ENST00000456328/) {print $0}}' annotation.gff3 | head -1
// // 1       havana  lnc_RNA 11869   14409   .       +       .       ID=transcript:ENST00000456328;Parent=gene:ENSG00000223972;Name=DDX11L1-202;biotype=processed_transcript;tag=basic;transcript_id=ENST00000456328;transcript_support_level=1;version=2
// 
// awk '{ if($9 ~ /ENST00000450305/) {print $0}}' annotation.gff3 | head -1
// // 1       havana  pseudogenic_transcript  12010   13670   .       +       .       ID=transcript:ENST00000450305;Parent=gene:ENSG00000223972;Name=DDX11L1-201;biotype=transcribed_unprocessed_pseudogene;tag=basic;transcript_id=ENST00000450305;transcript_support_level=NA;version=2
// 
// awk '{ if($9 ~ /ENST00000488147/) {print $0}}' annotation.gff3 | head -1
// // 1       havana  pseudogenic_transcript  14404   29570   .       -       .       ID=transcript:ENST00000488147;Parent=gene:ENSG00000227232;Name=WASH7P-201;biotype=unprocessed_pseudogene;tag=basic;transcript_id=ENST00000488147;transcript_support_level=NA;version=1
// 
// // => the first exons are from pseudogenes or non-coding RNAs that is why the coordinates of the first gene starts way after the first exons



// head -300 annotation.gff3

// ##sequence-region   KI270756.1 1 79590
// ##sequence-region   KI270757.1 1 71251
// ##sequence-region   MT 1 16569
// ##sequence-region   X 1 156040895
// ##sequence-region   Y 2781480 56887902
// #!genome-build  GRCh38.p13
// #!genome-version GRCh38
// #!genome-date 2013-12
// #!genome-build-accession NCBI:GCA_000001405.28
// #!genebuild-last-updated 2020-09
// 1       GRCh38  chromosome      1       248956422       .       .       .       ID=chromosome:1;Alias=CM000663.2,chr1,NC_000001.11
// ###
// 1       .       biological_region       10469   11240   1.3e+03 .       .       external_name=oe %3D 0.79;logic_name=cpg
// 1       .       biological_region       10650   10657   0.999   +       .       logic_name=eponine
// 1       .       biological_region       10655   10657   0.999   -       .       logic_name=eponine
// 1       .       biological_region       10678   10687   0.999   +       .       logic_name=eponine
// 1       .       biological_region       10681   10688   0.999   -       .       logic_name=eponine
// 1       .       biological_region       10707   10716   0.999   +       .       logic_name=eponine
// 1       .       biological_region       10708   10718   0.999   -       .       logic_name=eponine
// 1       .       biological_region       10735   10747   0.999   -       .       logic_name=eponine
// 1       .       biological_region       10737   10744   0.999   +       .       logic_name=eponine
// 1       .       biological_region       10766   10773   0.999   +       .       logic_name=eponine
// 1       .       biological_region       10770   10779   0.999   -       .       logic_name=eponine
// 1       .       biological_region       10796   10801   0.999   +       .       logic_name=eponine
// 1       .       biological_region       10810   10819   0.999   -       .       logic_name=eponine
// 1       .       biological_region       10870   10872   0.999   +       .       logic_name=eponine
// 1       .       biological_region       10889   10893   0.999   -       .       logic_name=eponine
// 1       havana  pseudogene      11869   14409   .       +       .       ID=gene:ENSG00000223972;Name=DDX11L1;biotype=transcribed_unprocessed_pseu
// dogene;description=DEAD/H-box helicase 11 like 1 (pseudogene) [Source:HGNC Symbol%3BAcc:HGNC:37102];gene_id=ENSG00000223972;logic_name=havana_hom
// o_sapiens;version=5

// => the first lines are chromosome sizes, then few lines describe the genome build, then a first line indicate the start of the chromosome 1, then there is a serie of short regions with external id "logic_name=eponine". These corresponds to predicted Transcription Start Sites. And then finally the real features appear. These predicted TSS could be removed from all analysis since these are anyway putative features (these regions would thus be considered intergenic if no other feature overlap these)
// https://www.sanger.ac.uk/tool/eponine/

// $ awk 'OFS="\t" {if (substr($1, 1, 1) != "#" && $3 == "biological_region") print}' annotation.gff3 | head -100
// 1       .       biological_region       827637  827722  1       -       .       external_name=rank %3D 2;logic_name=firstef
// 1       .       biological_region       904315  905706  1.34e+03        .       .       external_name=oe %3D 1.03;logic_name=cpg

// => firstef regions are also predicted TSS regions
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3577817/#pone.0054843-Davuluri1
// To evaluate the prediction performance of our method, we have compared our method with four other open access human promoter prediction programs. We have selected these four methods among a number of promoter prediction tools because they have performed well in the previous comparative studies. The methods are FirstEF [52], Eponine [14], Ep3 [22] and ProSOM [53]. All these methods are free and publicly available. FirstEF is based on discriminant analysis to find potential first donor sites and CpG-related and non-CpG-related promoter regions. Ep3 uses large scale structural properties of DNA to locate promoter regions in whole genome sequence and ProSOM utilizes self-organizing maps to distinguish promoters from non-promoter sequences.

// And CpG island can be intergenic as well. We remove all these


// -> the command to filter and sort the gff file:
// awk 'OFS="\t" {if (substr($1, 1, 1) != "#" && $3 != "biological_region") print}' annotation.gff3 | head
// awk 'OFS="\t" {if (substr($1, 1, 1) != "#" && $3 != "biological_region") print}' annotation.gff3 | sort -k1,1V -k4,4n -k5,5n | head


// scaffolds should be kept
// https://www.biostars.org/p/223541/
// The scaffolds are needed, they are likely parts of the reference genome which could not yet in this assembly build be confidently assigned to chromosomes. This is typical of repeat containing scaffolds. They are likely to be small and not gene rich, but excluding them would cause short reads which should map to them to be forced to map to the chromosomes and cause problemes downstream, such as false SNPs etc.




process getting_kallisto_indexes {
  tag "${specie}"

  container = params.kallisto

  publishDir path: "${"${specie}/genome/sequence"}", mode: 'link'

  input:
		set specie, file(transcriptome) from Transcriptome_for_building_kallisto_indexes

  output:
    file('transcriptome_kallisto_index')

  script:
  """
      kallisto index --make-unique -i transcriptome_kallisto_index ${transcriptome}
  """

}




process getting_orgdb {
  tag "${specie}"

  container = params.annotationhub

  publishDir path: "${specie}/genome/annotation/R", mode: 'link'

  input:
		set specie, specie_long from Assembly_Ensembl_1.orgdb

  output:
    set specie, specie_long, file('orgdb.sqlite') into Orgdb_for_annotation

  shell:
  '''

		#!/usr/bin/env Rscript

		library(AnnotationHub)

		specie = '!{specie}'
		specie_long = '!{specie_long}'
		ncbi_orgdb_version = '!{ncbi_orgdb_version}'
		annotationhub_cache = '!{params.annotationhub_cache}'

		specie_initial =  sapply(strsplit(specie_long, '_')[[1]], substr, 1, 1) %>% {paste0(toupper(.[1]), .[2])}

		ah = AnnotationHub(cache = annotationhub_cache)

		orgdb = query(ah, paste0('ncbi/standard/', ncbi_orgdb_version, '/org.', specie_initial, '.eg.db.sqlite'))[[1]]

		AnnotationDbi::saveDb(orgdb, 'orgdb.sqlite')

  '''

}


Orgdb_for_annotation
	.join(GFF3_filtered_for_R)
	.set{GFF3_filtered_for_R_1}

process getting_R_annotation_files {
  tag "${specie}"

  container = params.bioconductor

  publishDir path: "${specie}/genome/annotation/R", mode: 'link'

  input:
		set specie, specie_long, orgdb, file(chr_size), file(gff3_raw), file(gff3_genes_only), file(gff3_without_pseudogenes_and_ncRNA) from GFF3_filtered_for_R_1

  output:
    file('*')

  shell:
  '''

		#!/usr/bin/env Rscript

		library(magrittr)
		library(AnnotationDbi)

		source('!{params.cactus_dir}/software/bin/export_df_to_bed.R')
		gff3_raw = '!{gff3_raw}'
		gff3_genes_only = '!{gff3_genes_only}'
		gff3_without_pseudogenes_and_ncRNA = '!{gff3_without_pseudogenes_and_ncRNA}'
		specie = '!{specie}'
		orgdb = AnnotationDbi::loadDb('!{orgdb}')
		df_chr_size = read.table('!{chr_size}', sep = '\t')
		specie_long = '!{specie_long}'

		# export annotation dataframe
		anno_df = rtracklayer::readGFF(gff3_genes_only) 
		entrez_id = mapIds(orgdb, keys = anno_df$gene_id,  column = 'ENTREZID', keytype = 'ENSEMBL', multiVals = 'first')
		anno_df %<>% dplyr::rename(gene_name = Name)
		anno_df %<>% dplyr::mutate(width = end - start)
		anno_df$entrez_id = entrez_id
		anno_df$chr = as.character(anno_df$seqid)
		anno_df %<>% dplyr::select(chr, start, end, width, strand, gene_name, gene_id, entrez_id)
		saveRDS(anno_df, 'df_genes_metadata.rds')

		# create the seqinfo object
		nr = nrow(df_chr_size)
		seqinfo = GenomeInfoDb::Seqinfo(seqnames = df_chr_size[,1], seqlengths = df_chr_size[,2], isCircular = rep(F, nr) , genome = rep(specie_long, nr))
		saveRDS(seqinfo, 'seqinfo.rds')

		# export txdb
		txdb = GenomicFeatures::makeTxDbFromGFF(gff3_without_pseudogenes_and_ncRNA, chrominfo = seqinfo)
		AnnotationDbi::saveDb(txdb, file="txdb.sqlite")
		
		# export gene vs transcripts
		txdb1 = GenomicFeatures::makeTxDbFromGFF(gff3_raw)
		df_genes_transcripts = select(txdb1, keys = keys(txdb1), columns = c('GENEID', 'TXNAME', 'TXCHROM'), keytype = 'GENEID')
		saveRDS(df_genes_transcripts, 'df_genes_transcripts.rds')

		# extract promoters
		promoters = GenomicFeatures::promoters(GenomicFeatures::genes(txdb), upstream = 1500, downstream = 500)

		# export promoters as dataframe
		promoters_df = as.data.frame(promoters, stringsAsFactors = F)
		promoters_df$start[promoters_df$start < 0] = 0
		promoters_df$end[promoters_df$end < 0] = 0
		promoters_df %<>% dplyr::rename(gene_name = gene_id)
		colnames(promoters_df)[colnames(promoters_df) == 'gene_name'] = 'gene_id'
		promoters_df %<>% dplyr::inner_join(anno_df %>% dplyr::select(gene_id, gene_name, entrez_id), by = 'gene_id')
		promoters_df %<>% dplyr::arrange(seqnames, start)
		promoters_df %<>% dplyr::rename(chr = seqnames)
		saveRDS(promoters_df, file = 'promoters_df.rds')

		# export promoters as bed file
		promoters_df1 = promoters_df
		promoters_df1$score = 0
		promoters_df1 %<>% dplyr::select(chr, start, end, gene_name, score, strand, gene_id)
		export_df_to_bed(promoters_df1, 'promoters.bed')


  '''

}

// source('!{cactusdir}/src/R_analysis/general_functions.R')

// cur_seq_info = rtracklayer::SeqinfoForUCSCGenome('ce11')
// cur_seq_info@seqnames %<>% gsub('chr', '', .)



// rtracklayer::export(promoters, 'promoters.bed')

// unfortunately, the "promoters function" extract promoters with coordinate below zero (I found this bug for 3 mitochondrial genes in C elegans). Thus we need to put them at zero.

// promoters_df %<>% dplyr::select(chr:strand, gene_name, gene_id, entrez_id)
// promoters_df %<>% dplyr::mutate(seqnames = as.character(seqnames), strand = as.character(strand))


// Rle(promoters_df$seqnames)
// factor-Rle of length 20191 with 7 runs
//   Lengths:  2908  3520  2683  3285  4994  2789    12
//   Values :     I    II   III    IV     V     X MtDNA



// df_genes_transcripts = select(txdb, keys = keys(txdb), columns = c('GENEID', 'TXNAME'), keytype = "GENEID")
// dim(df_genes_transcripts)
// [1] 33552     2

// a = select(txdb, keys = keys(txdb)[1:5], columns = c('GENEID', 'TXID', 'TXNAME', 'TXSTART', 'TXEND', 'TXSTRAND'), keytype = "GENEID")
//   GENEID  TXID               TXNAME TXSTRAND  TXSTART    TXEND
// 1 2L52.1  5291 transcript:2L52.1a.1        +     1848     4717
// 2 2L52.1  5292 transcript:2L52.1b.1        +     3410     4663
// 3 4R79.2 22496 transcript:4R79.2b.1        - 17480396 17482329
// 4 4R79.2 22497 transcript:4R79.2a.1        - 17480396 17483332

// columns(txdb)
// [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"
// [6] "CDSSTART"   "CDSSTRAND"  "EXONCHROM"  "EXONEND"    "EXONID"
// [11] "EXONNAME"   "EXONRANK"   "EXONSTART"  "EXONSTRAND" "GENEID"
// [16] "TXCHROM"    "TXEND"      "TXID"       "TXNAME"     "TXSTART"
// [21] "TXSTRAND"   "TXTYPE"

// vec = anno_df$gene_name %>% setNames(., anno_df$gene_id)
// saveRDS(vec, 'gene_id_to_name_vector.rds', version = 2)
// 
// vec = anno_df$entrez_id %>% setNames(., anno_df$gene_id)
// saveRDS(vec, 'gene_ensembl_id_to_entrez_id_vector.rds', version = 2)


// saveDfFromGRangesAsBed(all_promoters, file = 'promoters.bed')
// this line of code was replaced by:
// rtracklayer::export(promoters, 'promoters.bed')

 
 
// a = GenomicFeatures::promoters(GenomicFeatures::genes(txdb), upstream = 1500, downstream = 500)
// rtracklayer::export(a, 'test.bed')
// rtracklayer::export(a, 'test.bed', format = 'bed')
// 
// txdb = 
// a = annotatePeak(gr, c(-200, 200), level = 'gene')



// txdb = GenomicFeatures::makeTxDbFromGFF("/home/jersal/lluis/atac/bin/igenomes/references/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gff3")
// 
// txdb = GenomicFeatures::makeTxDbFromGFF("genes.gff3")
// 
// anno_df1 = rtracklayer::readGFF(anno_file) 
// gr = GenomicRanges::makeGRangesFromDataFrame(anno_df1)



// => in fact I don't need a txdb object
// GenomicFeatures::promoters(a, upstream = 500, downstream = 500)
// GenomicFeatures::promoters(a, upstream = 1, downstream = 1)


// library(biomaRt)

// specie_initial = paste(sapply(strsplit(tolower(specie), '_')[[1]], function(x) substr(x, 1,1)), collapse = '')


// dataset = switch(specie, 
// 	Caenorhabditis_elegans = 'celegans_gene_ensembl', 
// 	Drosophila_melanogaster = 'dmelanogaster_gene_ensembl', 
// 	Mus_musculus = 'mmusculus_gene_ensembl', 
// 	Homo_sapiens = 'hsapiens_gene_ensembl'
// )
// 
// mart = useMart("ensembl", dataset = dataset)
// 
// ensembl_gene_id = gff3_df$gene_id
// df_biomart = getBM( filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id', 'entrezgene_id'), values = ensembl_gene_id, mart = mart )
// 
// gff3_df$entrez_gene_id = df_biomart$entrezgene_id[match(gff3_df$gene_id, df_biomart$ensembl_gene_id)]



// mart = useMart("ensembl")
// a = listDatasets(mart)
// a[grepl('elegans|hsapiens|musculus|melanogaster', a$dataset), ]
// 									 dataset                              description    version
// 48       hsapiens_gene_ensembl                 Human genes (GRCh38.p13) GRCh38.p13
// 53  dmelanogaster_gene_ensembl Drosophila melanogaster genes (BDGP6.28)   BDGP6.28
// 104      celegans_gene_ensembl  Caenorhabditis elegans genes (WBcel235)   WBcel235
// 130     mmusculus_gene_ensembl                  Mouse genes (GRCm38.p6)  GRCm38.p6

// ## for C elegans
// dim(gff3_df) # 31293
// gff3_df1 = gff3_df[, c('gene_id', 'entrez_gene_id')] %>% .[!duplicated(.), ]
// dim(gff3_df1) # 20447
// length(which(is.na(gff3_df1$entrez_gene_id))) # => 563 id were lost out of 20K. Not great but still reasonnable.



// on completion
 workflow.onComplete
 {
  println ""
  println "Workflow completed on: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'Succeeded' : 'Failed' }"
  println "Workflow Duration: $workflow.duration"
  println ""
 }
 
