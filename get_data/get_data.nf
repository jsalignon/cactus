

params.number_of_cores = 8   // how to initialize if missing?
number_of_cores = params.number_of_cores

params.ensembl_release  = '102'
ncbi_orgdb_version      = '3.14'
homer_genomes_version   = '6.4'
homer_organisms_version = '6.3'
homer_promoters_version = '5.5'

params.homer_odd_score_threshold = 0.5

Species_channel = Channel
	.from( 	
		[
		//                                    assembly  assembly                       assembly
		// specie    scientific_name           Ensembl  NCBI       assembly type       ENCODE
			['worm',  'caenorhabditis_elegans',  'ce11', 'WBcel235', 'toplevel',         'ce11' ],
			['fly',   'drosophila_melanogaster', 'dm6',  'BDGP6.28', 'toplevel',         'dm6' ],
			['mouse', 'mus_musculus',            'mm10', 'GRCm38',   'primary_assembly', 'mm10' ],
		  ['human', 'homo_sapiens',            'hg38', 'GRCh38',   'primary_assembly', 'GRCh38' ]
		]
	)
	.dump(tag: 'start_channel')
	.multiMap { it ->
			  	homer: it[0, 2]
					 pwms: it[0, 1]
      blacklist: it[0, 2, 3]
					 chip: it[0, 5]
chromatin_state: it[0, 5]
					fasta: it[0, 1, 3, 4]
					orgdb: it[0, 1]
	}
	.set { Start_channel }




process get_homer_data {
	tag "${specie}"

	container = params.samtools_bedtools_perl

	publishDir path: "${specie}/homer_data", mode: 'link'
	input:
		set specie, specie_code from Start_channel.homer

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

Start_channel.pwms
	.combine(cisbp_motifs_all_species)
	.set{ cisbp_motifs_all_species_1 }


process split_pwms {

	container = params.r_basic

	publishDir path: "${specie}/motifs_PWMs", mode: 'link'

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

	publishDir path: "${specie}/blacklisted_regions", mode: 'link'

	input:
		set specie, specie_code, ncbi_code from Start_channel.blacklist

	output:
		file("${specie_code}_blacklist_Ensembl.bed")

	shell:
	'''
			
		specie_code='!{specie_code}'
		ncbi_code='!{ncbi_code}'
		
		url_blacklist="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists"
		url_mapping="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master"

		wget -O ${specie_code}_blacklist_NCBI.bed.gz $url_blacklist/${specie_code}-blacklist.v2.bed.gz
		gunzip  ${specie_code}_blacklist_NCBI.bed.gz
		
		ncbi_code_1=${ncbi_code%%.*}
				
		wget -O NCBI_to_Ensembl.txt $url_mapping/${ncbi_code_1}_UCSC2ensembl.txt
		
		cvbio UpdateContigNames -i ${specie_code}_blacklist_NCBI.bed -o ${specie_code}_blacklist_Ensembl.bed -m NCBI_to_Ensembl.txt --comment-chars '#' --columns 0 --skip-missing true
		
	'''
}


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
		set specie, assembly from Start_channel.chip

	output:
		set specie, file("dt_encode_chip.rds") into (Encode_chip_metadata_channel_1, Encode_chip_metadata_channel_2)		

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

	publishDir path: "${specie}/CHIP", mode: 'link'
	
	input:
		set specie, file(dt_encode_chip_rds) from Encode_chip_metadata_channel_1

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
		
		dt = copy(dt1)
		
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

publishDir path: "${specie}/CHIP/files", mode: 'link'

input:
	set specie, file(dt_encode_chip_rds) from Encode_chip_metadata_channel_2

output:
	file("*.bed")

shell:
'''

	#!/usr/bin/env Rscript

	library(data.table)
	library(magrittr)
	library(purrr)

	specie = '!{specie}'
	dt1 = readRDS('!{dt_encode_chip_rds}')
	
	dt = copy(dt1)
	
	# donwloading the data
	url_encode = 'https://www.encodeproject.org' 
	sapply(1:nrow(dt), function(c1) download.file(url = paste0(url_encode, dt$href[c1]), quiet = T, destfile = dt$local_file[c1], method = 'curl', extra = '-L' ))

	dt[, md5sum_downloaded_files := tools::md5sum(local_file)]
	if(any(dt$md5sum != dt$md5sum_downloaded_files)) stop('not all md5 sums are equal')

	system('for z in *.gz; do gunzip "$z"; done')
	
'''
}


Start_channel.chromatin_state
	.filter{ specie, assembly -> specie in ['human', 'mouse']}
	.set{ Start_channel_chromatin_state }


process get_encode_chromatin_state_metadata {
	tag "${specie}"

	container = params.encodeexplorer

	input:
		set specie, assembly from Start_channel_chromatin_state

	output:
		set specie, file("dt_encode_chromatin_states.rds") into Encode_chromatin_state_metadata_channel

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
    
		/labs/zhiping-weng/
    
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
	
		dt = dt[, c('local_file', 'assembly', 'description', 'alias', 'timepoint', 'timepoint_units', 'life_stage', 'term_name', 'classification', 'cell_slims', 'organ_slims', 'developmental_slims', 'system_slims', 'organism', 'accession', 'href', 'md5sum', 'file_size')]
		
		if(assembly == 'mm10') dt = dt[grep('18states', dt$alias)]
		
		saveRDS(dt, 'dt_encode_chromatin_states.rds')


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

	// container = params.huge_container
	// could not find any container that works. Doing it without containers for now

	input:
		set specie, file(dt_encode_chromatin_state_rds) from Encode_chromatin_state_metadata_channel

	output:
		file("dt_encode_chip.rds") 

	shell:
	'''

		#!/usr/bin/env Rscript
		
		library(data.table)
		library(bedr)
		
		dt = readRDS('!{dt_encode_chromatin_state_rds}')
		
		url_encode = 'https://www.encodeproject.org' 
		
		## processing files one by one as they are huge to not fill up the server
		for(c1 in 1:nrow(dt)) {
			
			accession = dt$accession[c1]
			local_file = paste0(accession, '.bed.gz')
			download.file(url = paste0(url_encode, dt$href[c1]), quiet = T, destfile = local_file, method = 'curl', extra = '-L' )
			md5sum_downloaded_files = tools::md5sum(local_file)
			if(dt$md5sum[c1] != md5sum_downloaded_files) stop(paste('md5 sums not equal for file', local_file))
			cur_bed = rtracklayer::import(local_file)
			system(paste('gunzip', local_file))
		}
			grep "Quies" ENCFF786MCR.bed | bedtools merge -i - > ENCFF786MCR_merged_Quies.bed

			} function(c1) {
			
				
			dt[, md5sum_downloaded_files := tools::md5sum(local_file)]
			if(any(dt$md5sum != dt$md5sum_downloaded_files)) stop('not all md5 sums are equal')
			system('for z in *.gz; do gunzip "$z"; done')
	
			.
			      Enh    EnhLo1    EnhLo2  EnhPois1  EnhPois2   HetCons    HetFac     Quies
			    15563     44844     27291     61335    304983    155305    166918   8910113
			   QuiesG      TssA TssAFlnk1 TssAFlnk2    TssBiv       Tx1       Tx2
			  3259122     41348     11619     65583     36083     19092    508479
			>
				
				
	'''
}



process get_fasta_and_gff {
	tag "${specie}"

	container = params.samtools_bedtools_perl

	publishDir path: "${specie}/genome", mode: 'link', saveAs: {
    if (it.indexOf(".gff3") > 0) "annotation/${it}"
		else if (it.indexOf(".fa") > 0) "sequence/${it}"
    else if (it.indexOf(".txt") > 0) "${it}"
  }

	input:
		set specie, specie_long, genome, assembly from Start_channel.fasta

	output:
		set specie, file('annotation.gff3'), file('genome.fa') into (Genome_and_annotation, Genome_and_annotation_1, Genome_and_annotation_2, Genome_and_annotation_3)

	shell:
	'''
		
	specie='!{specie}'
	specie_long='!{specie_long}'
	genome='!{genome}'
	assembly='!{assembly}'
	release='!{params.ensembl_release}'
	
	URL=ftp://ftp.ensembl.org/pub/release-$release
	
	gff3_file=${specie_long^}.$genome.$release.gff3.gz
	fasta_file=${specie_long^}.$genome.dna_sm.$assembly.fa.gz
	
	wget -O annotation.gff3.gz $URL/gff3/$specie_long/$gff3_file 
	gunzip annotation.gff3.gz 
	
	wget -O genome.fa.gz $URL/fasta/$specie_long/dna/$fasta_file 
	gunzip genome.fa.gz
	
	echo "gff3 file : $gff3_file" > README.txt
	echo "fasta file: $fasta_file" >> README.txt
		
		
	'''
}



process filtering_annotation_file {
  tag "${specie}"

  container = params.samtools_bedtools_perl

  publishDir path: "${specie}/genome/annotation", mode: 'link'

  input:
    set specie, file(gff3), file(fasta) from Genome_and_annotation

  output:
		set specie, file('anno_genes_only.gff3'), file('gff3_without_pseudogenes_and_ncRNA.gff3') into Channel_gff3_filtered

  shell:
  '''

		grep -i protein_coding !{gff3} | awk  -F"\t"  '{if($3 == "gene") print $0}' - > anno_genes_only.gff3

		grep -viE "pseudogen|ncRNA" !{gff3} > gff3_without_pseudogenes_and_ncRNA.gff3


  '''

}

// Notes

// grep -i protein_coding ${gff} > protein_coding_mRNA_and_genes.gff3

// anno_without_pseudogenes.gff3 is used for creating a txdb database without non coding transcripts, for the annotation of ATAC seq peaks

// protein_coding_genes.gff3 is used to create the gene metadata table for all analysis


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



process getting_chromosome_length_and_transcriptome {
  tag "${specie}"

	container = params.gffread

  publishDir path: "${"${specie}/genome/sequence"}", mode: 'link'

  input:
    set specie, file(gff3), file(fasta), file(fasta_indexes) from Fasta_fai_gff3

  output:
    set specie, file('chromosome_size.txt') into Chromosome_size_for_seqinfo
		set specie, file('transcriptome.fa') into Transcriptome_for_building_kallisto_indexes

  script:
  """
      cut -f1-2 ${fasta_indexes} > chromosome_size.txt
      gffread -C ${gff3} -g ${fasta} -w transcriptome.fa
  """

}





process generating_kallisto_transcriptome_indexes {
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

  publishDir path: "${specie}/orgdb", mode: 'link'

  input:
		set specie, specie_long from Start_channel.orgdb

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
	.join(Chromosome_size_for_seqinfo)
	.join(Channel_gff3_filtered)
	.set{Channel_gff3_filtered_1}

process getting_R_annotation_files {
  tag "${specie}"

  container = params.bioconductor

  publishDir path: "${specie}/genome/annotation", mode: 'link'

  input:
		set specie, specie_long, orgdb, file(chr_size), file(gff3_genes_only), file(gff3_without_pseudogenes_and_ncRNA) from Channel_gff3_filtered_1

  output:
    file('*')

  shell:
  '''

		#!/usr/bin/env Rscript

		library(magrittr)
		library(AnnotationDbi)

		source('!{params.cactus_dir}/software/bin/export_df_to_bed.R')
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
		saveRDS(anno_df, 'df_genes_metadata.rds', version = 2)

		# create the seqinfo object
		nr = nrow(df_chr_size)
		seqinfo = GenomeInfoDb::Seqinfo(seqnames = df_chr_size[,1], seqlengths = df_chr_size[,2], isCircular = rep(F, nr) , genome = rep(specie_long, nr))
		saveRDS(seqinfo, 'seqinfo.rds')

		# export txdb
		txdb = GenomicFeatures::makeTxDbFromGFF(gff3_without_pseudogenes_and_ncRNA, chrominfo = seqinfo)
		AnnotationDbi::saveDb(txdb, file="txdb_without_pseudogenes_and_ncRNA.sqlite")

		# export gene vs transcripts 
		df_genes_transcripts = select(txdb, keys = keys(txdb), columns = c('GENEID', 'TXNAME'), keytype = 'GENEID')
		saveRDS(df_genes_transcripts, 'df_genes_transcripts.rds', version = 2)

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
		saveRDS(promoters_df, file = 'promoters_df.rds', version = 2)

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
 
