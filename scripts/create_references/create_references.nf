


params.number_of_cores = 8   // how to initialize if missing?
number_of_cores = params.number_of_cores

ncbi_orgdb_version      = '3.14'
homer_genomes_version   = '6.4'
homer_organisms_version = '6.3'
homer_promoters_version = '5.5'

params.homer_odd_score_threshold = 0.5

Ensembl_Assembly = Channel
	.from(
		[   //  0          1                     2              3                4     
		//                                     assembly      assembly          Ensembl            
		// species    scientific_name            name           type            Release
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
	.set { Ensembl_Assembly_1 }

//// checking Ensembl latest release and fly genome build name (BDGP6.XX)
// https://www.ensembl.org/index.html?redirect=no
// Ensembl Release 107 (Jul 2022)
// http://ftp.ensembl.org/pub/release-107/gff3/drosophila_melanogaster/

//// UCSC releases and details:
// https://genome.ucsc.edu/FAQ/FAQreleases.html

//// A note on the mouse genome:
// ENCODE (chip and chromatin states), blacklisted regions and homer genome are not available yet fo mm11/GRCm39. Therefore, we need to stick with release 102, which is the latest release for mm10/GRCm38. I tried to do the genomic shift myself using liftover of all bed files (see commits 6ba5f3a and 894c723). However, I could not succeed in generating the homer genome within a homer container. Therefore, an older genome build will be used for mouse for now.


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

Encode_Assembly = Channel
	.from( 	
		[  
			['worm',   'ce11' ],
			['fly',    'dm6' ],
			['mouse',  'mm10' ],
			['human',  'GRCh38' ]
		]
	)
	.tap{ Encode_Assembly_for_CHIP }
	.filter{ it[0] in ['human', 'mouse']}
	.set{ Encode_Assembly_for_chromatin_state }



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
	




process getting_homer_data {
	tag "${species}"

	label "samtools_bedtools_perl"

	publishDir path: "${species}/homer_data", mode: 'link'
	input:
		set species, species_code from Assembly_nickname_1.homer

	output:
		file('*')

	shell:
	'''

		species='!{species}'
		species_code='!{species_code}'
		homer_genomes_version='!{homer_genomes_version}'
		homer_organisms_version='!{homer_organisms_version}'
		homer_promoters_version='!{homer_promoters_version}'

		URL=http://homer.ucsd.edu/homer/data

		wget -nc $URL/genomes/${species_code}.v${homer_genomes_version}.zip
		wget -nc $URL/organisms/${species}.v${homer_organisms_version}.zip
		wget -nc $URL/promoters/${species}.v${homer_promoters_version}.zip

		for z in *.zip; do unzip "$z"; done
		mv data/* .
		rm -r data *.zip

	'''
}



process getting_pwms {

	label "r_basic"

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
		pnrow <- function(x) print(nrow(x))

		dt_all = fread('TF_Information.txt') %T>% pnrow # 399265

		encode_species = c('Caenorhabditis_elegans', 'Drosophila_melanogaster', 'Mus_musculus', 'Homo_sapiens')
		dt_cisbp = dt_all[TF_Species %in% encode_species] %T>% pnrow # 10329
		dt_cisbp = dt_cisbp[Motif_ID != '.'] %T>% pnrow # 8628
		dt_cisbp = dt_cisbp[!duplicated(TF_ID)] %T>% pnrow # 2936

		# removing emtpy motifs
		nrows = sapply(dt_cisbp$Motif_ID, function(x) nrow(read.table(paste0('pwms/', x, '.txt'), sep = '\t')))
		dt_cisbp = dt_cisbp[nrows != 1] %T>% pnrow # 2779

		# loading motifs and determining odd scores, thresholds and consensus sequences
		dt_cisbp$motif = lapply(dt_cisbp$Motif_ID, function(x) read.table(paste0('pwms/', x,'.txt'), sep = '\t', header = T, stringsAsFactors = F)[, -1])
		dt_cisbp[, log2_odd_score := sapply(motif, function(motif1) sum(apply(motif1, 1, function(cur_row) log2(max(cur_row)/0.25))))]
		dt_cisbp[, homer_threshold := log2_odd_score * homer_odd_score_threshold]
		dt_cisbp[, consensus_sequence := sapply(motif, function(motif1) paste0(colnames(motif1)[apply(motif1, 1, function(x) which.max(x))], collapse = ''))]
		saveRDS(dt_cisbp, 'dt_cisbp_encode.rds')


	'''

}

// # filtering the cisbp database
// nice explanations about the cisbp files here:
// https://github.com/Danko-Lab/rtfbs_db/blob/a29f73b24e8984aad8bcd1227b06446ebef44b08/rtfbsdb/man/CisBP.getTFinformation.Rd
// Three TF information files in CisBP dataset.\cr\cr
// 1: TF_Information.txt : (direct motifs) or (no direct but inferred motifs with 90\%)\cr
// 2: TF_Information_all_motifs.txt: (direct motifs) and (inferred motifs above the threshold)\cr
// 3: F_Information_all_motifs_plus.txt: All motifs\cr
// 
// 'TF_Information_all_motifs.txt' is a superset of 'TF_Information.txt'.  It also includes any motif that can be inferred for a given TF, given the TF family-specific threshold.  For example, if a TF has a directly determined motif, and two TFs with motifs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the direct motif, but 
// TF_Information_all_motifs.txt will include all three motifs.  Likewise, if a TF does not have a direct motif, but has two TFs with 90% and 80% identical DBDs to it, TF_Information.txt will only include the motif from the 90% indentical TF, while TF_Information_all_motifs.txt would include both.\cr\cr

// However, there are still multiple entries per TF
// dt_cisbp[TF_Species == 'Caenorhabditis_elegans']$TF_Name %>% table %>% sort %>% rev %>% head(10)
//   hlh-1   skn-1  snpc-4   pha-4   efl-1   tra-1 nhr-182   mab-3  lin-39   eor-1
//       7       4       3       3       3       2       2       2       2       2
// dt_cisbp[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'skn-1']
// dt_cisbp[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'hlh-1']
// dt_cisbp[TF_Species == 'Caenorhabditis_elegans' & TF_Name == 'snpc-4']
// 
// I think I understand what happened: the entries that are duplicated are all from the similarity regression method. And they all have the same score. So instead of picking one arbitrarily, CISBP kept all of these entries with the same score.
// Conclusions: these are probably exactly the same motifs (or very very similar). So we can just keep the first of them!

// removing empty motifs
// https://github.com/TedCCLeung/PhMotif/blob/9ace03a196e02375b3024af3506b9f58e98f1a1d/data-raw/processCISB


Ensembl_Assembly_1.pwms
	.combine(cisbp_motifs_all_species)
	.set{ cisbp_motifs_all_species_1 }


process splitting_pwms {
	tag "${species}"

	label "r_basic"

	publishDir path: "${species}/homer_data", mode: 'link'

	input:
		set species, species_long, file(dt_cisbp_encode_rds) from cisbp_motifs_all_species_1

	output:
		file('homer_motifs.txt')

	shell:
	'''
		#!/usr/bin/env Rscript

		library(data.table)

		dt_cisbp_encode = readRDS('!{dt_cisbp_encode_rds}')
		species_long = '!{species_long}'

		species_long_1 = paste0(toupper(substr(species_long, 1, 1)), substr(species_long, 2, nchar(species_long)))
		dt_cisbp = dt_cisbp_encode[TF_Species == species_long_1]

		sink('homer_motifs.txt')

			apply(dt_cisbp, 1, function(x) {
					cat(paste0('>', x$consensus_sequence, '\t', x$TF_Name, '\t', x$homer_threshold, '\n'))
					cat(paste0(apply(x$motif, 1, paste0, collapse = '\t'), '\n'))
			})

		sink()

	'''

}





process getting_blacklisted_regions {
	tag "${species}"

	label "cvbio"

	publishDir path: "${species}/genome/annotation", mode: 'link'

	input:
		set species, species_code, ncbi_code from Assembly_nickname_1.blacklist

	output:
		file("blacklisted_regions.bed")

	shell:
	'''

		species_code='!{species_code}'
		ncbi_code='!{ncbi_code}'

		url_blacklist="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists"
		url_mapping="https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master"
		wget -O ${species_code}_blacklist_NCBI.bed.gz $url_blacklist/${species_code}-blacklist.v2.bed.gz
		gunzip  ${species_code}_blacklist_NCBI.bed.gz

		wget -O NCBI_to_Ensembl.txt $url_mapping/${ncbi_code}_UCSC2ensembl.txt

		cvbio UpdateContigNames -i ${species_code}_blacklist_NCBI.bed -o blacklisted_regions.bed -m NCBI_to_Ensembl.txt --comment-chars '#' --columns 0 --skip-missing true

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




process getting_chip_metadata {
	tag "${species}"

	label "encodexplorer"

	publishDir path: "${species}", mode: 'link', saveAs: {
    if (it.equals("encode_chip_metadata.csv")) "${it}"
  }
	
	input:
		set species, assembly from Encode_Assembly_for_CHIP

	output:
		set species, file("dt_encode_chip.rds") into Encode_chip_metadata		
		set species, file("*__split_dt.rds") into Encode_chip_split_dt mode flatten
		file("encode_chip_metadata.csv")

	shell:
	'''
		#!/usr/bin/env Rscript
		
		source('!{params.bin_dir}/encode_chip_functions.R')
		assembly = '!{assembly}'
		
		library(data.table)
		library(magrittr)
		library(purrr)
		library(jsonlite)
		library(RCurl)
		
		url_encode = 'https://www.encodeproject.org' 
		url_search = paste0(url_encode, '/search/?')
		url_append = '&frame=object&format=json&limit=all'

		## getting all tables
		query_chip_files = paste0('type=File&file_format=bed&output_type=optimal+IDR+thresholded+peaks&assay_title=TF+ChIP-seq&status=released&assembly=', assembly)
		df_chip_files = get_encode_df(query_chip_files) %T>% pnrow # 2716
		
		query_experiments = paste0('type=Experiment&assay_title=TF+ChIP-seq&status=released&assembly=', assembly)
		df_experiments = get_encode_df(query_experiments) %T>% pnrow # 
		
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
		
		# creating final unique CHIP names
		dt[, target_symbol_1 := target_symbol %>% gsub('-', '', .) %>% gsub('.', '', ., fixed = T) %>% gsub(')', '', ., fixed = T) %>% gsub('(', '', ., fixed = T) %>% substr(., 1, 7)]
		dt[, chip_name := paste0(target_symbol_1, '_', stage_cell)]
		# most entries are unique. We add an identifyier for the duplicated entries
		dt$chip_name %>% table %>% as.integer %>% summary 
		dt[, chip_name := paste0(chip_name, '_', 1:.N), chip_name]
		dt[, chip_name := chip_name %>% gsub('_1$', '', .)]
		dt[, local_file := paste0(chip_name, '.bed.gz')]
		# reordering columns
		dt = dt[, c('chip_name', 'target_symbol', 'target_symbol_1', 'stage_cell', 'target_id', 'life_stage_age', 'ontology', 'classification', 'cell_slims', 'organ_slims', 'developmental_slims', 'system_slims', 'biosample_summary', 'organism', 'assembly', 'file', 'href', 'md5sum', 'local_file')]
		
		# saving tables
		dt1 = dt[, 1:(ncol(dt)-4)]
		write.csv(dt1, 'encode_chip_metadata.csv', row.names = F)
		saveRDS(dt, 'dt_encode_chip.rds')
		ldt = split(dt, dt$chip_name)
		walk(ldt, ~saveRDS(.x, paste0(.x$chip_name, '__split_dt.rds')))
		
	'''
}

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



process making_chip_ontology_groups {
	tag "${species}"

	label "encodexplorer"

	publishDir path: "${species}", mode: 'link', saveAs: {
    if (it.indexOf(".txt") > 0) "${it}"
		else if (it.indexOf(".tsv") > 0) "CHIP/${it}"
  }

	input:
		set species, file(dt_encode_chip_rds) from Encode_chip_metadata

	output:
		file("*")

	shell:
	'''
		#!/usr/bin/env Rscript

		library(data.table)
		library(magrittr)
		library(purrr)

		species = '!{species}'
		dt1 = readRDS('!{dt_encode_chip_rds}')
		source('!{params.bin_dir}/make_chip_ontology_groups_functions.R')
		source('!{params.bin_dir}/get_tab.R')

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

		sink('chip_ontology_groups_sizes.txt')
			 print(dt_n_chip_by_ontology)
		sink()

		tb = tibble::enframe(l_ontology_chip_names) %>% {tidyr::unnest(., cols = c(value))}
		write.table(tb, file = 'chip_ontology_groups.tsv', col.names = F, row.names = F, quote = F, sep = get_tab())


	'''
}



process getting_chip_data {
tag "${species}"

label "encodexplorer"

publishDir path: "${species}/CHIP", mode: 'link'

input:
	set species, file(dt_rds) from Encode_chip_split_dt

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


process getting_encode_chromatin_state_metadata {
	tag "${species}"

	label "encodexplorer"

	publishDir path: "${species}", mode: 'link', pattern: "*.csv"

	input:
		set species, assembly from Encode_Assembly_for_chromatin_state

	output:
		set species, assembly, file("*.rds") into Encode_chromatin_state_split_dt mode flatten
		file('encode_chromatin_states_metadata.csv')

	shell:
	'''
		#!/usr/bin/env Rscript

		source('!{params.bin_dir}/encode_chip_functions.R')
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
		write.csv(dt1, 'encode_chromatin_states_metadata.csv')

		ldt = split(dt, dt$accession)
		walk(ldt, ~saveRDS(.x, paste0(.x$accession, '.rds')))


	'''
}

// search link
// https://www.encodeproject.org/search/?type=File&annotation_type=chromatin+state&assembly=GRCh38&assembly=mm10&status=released&file_format=bed

// query results when using all species
// df_ecsm_files = get_encode_df('type=File&annotation_type=chromatin+state&status=released&file_format=bed&assembly=GRCh38&assembly=mm10') %T>% pnrow # 962
// df_biosample_types = get_encode_df('type=BiosampleType') %T>% pnrow # 936
// df_annotations = get_encode_df('type=Annotation&annotation_type=chromatin+state&status=released') %T>% pnrow # 1251

// All experiments have been made in a 15 and 18 state model
// However, the 15 state model is much bigger
// > dt[alias %in% c("zhiping-weng:chromhmm-8marks-15states-mm10-heart-12.5", "zhiping-weng:chromhmm-10marks-18states-mm10-heart-12.5")][, c('description', 'file_size')]
//                                                description file_size
// 1: ChromHMM 15-state model for heart (embryonic 12.5 days) 121270263
// 2: ChromHMM 18 state model for heart (embryonic 12.5 days)   8684869

// dt[grep('18states', dt$alias)]$file_size %>% max #   8684869
// dt[grep('15states', dt$alias)]$file_size %>% min # 109950661

// => we will use the 18 state models only

// dt[local_file == 'ENCFF686HVR.bed.gz']




process getting_encode_chromatin_state_data {
	tag "${species}, ${dt_rds}"

	label "encodexplorer"

	publishDir path: "${species}/chromatin_states", mode: 'link'

	input:
		set species, assembly, file(dt_rds) from Encode_chromatin_state_split_dt

	output:
		file("*")
		// set species, assembly, file("*.bed") into Encode_bed_for_liftover

	shell:
	'''
		#!/usr/bin/env Rscript

		library(magrittr)
		library(data.table)
		dt = readRDS('!{dt_rds}')
		source('!{params.bin_dir}/get_tab.R')

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

		file.remove(bed_file)

	'''
}




HiHMM_chromatin_states
	.combine(HiHMM_liftover, by:0)
	.set{ HiHMM_chromatin_states_1 }



process getting_hihmm_chromatin_state_data_part_1 {
	tag "${species}"

	// label "curl"
	// label "ubuntu"
	// label "hihmm"

	input:
		set species, bed_name, original_assembly, liftover_name from HiHMM_chromatin_states_1

	output:
		set species, bed_name, file('*.bed'), original_assembly, liftover_name into HiHMM_chromatin_states_2


	shell:
	'''

			bed_name="!{bed_name}"
			species="!{species}"

			bed_file="${bed_name}.bed"
			bed_file_lifted="${bed_name}_lifted.bed"

			hiHMM_path="http://compbio.med.harvard.edu/modencode/webpage/hihmm"

			wget $hiHMM_path/$bed_file

			if [ $species = 'worm' ]; then
				gawk -i inplace '{print "chr"$0}' $bed_file
			fi
			
	'''
}


// => could not get any containers to work
// with wegt:   wget: bad address 'compbio.med.harvard.edu'
// wget $hiHMM_path/$bed_file
// with curl:  FATAL:   stat /bin/bash: no such file or directory


process getting_hihmm_chromatin_state_data_part_2 {
	tag "${species}"

	publishDir path: "${species}/chromatin_states", mode: 'link'

	label "liftover"

	input:
		set species, bed_name, file(bed_file), original_assembly, liftover_name from HiHMM_chromatin_states_2

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




process getting_contamination_bowtie2 {
	tag "${species}"

	label "bowtie2_samtools"

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



process getting_fasta_and_gff {
	tag "${species}"

	label "samtools_bedtools_perl"

	publishDir path: "${species}/genome", mode: 'link', saveAs: {
    if (it.indexOf(".gff3") > 0) "annotation/${it}"
		else if (it.indexOf(".fa") > 0) "sequence/${it}"
    else if (it.indexOf(".txt") > 0) "${it}"
  }

	input:
		set species, species_long, genome, assembly, release from Ensembl_Assembly_1.fasta

	output:
		set species, file('annotation.gff3'), file('genome.fa') into (Genome_and_annotation, Genome_and_annotation_1, Genome_and_annotation_2, Genome_and_annotation_3)
		set species, file('annotation.gff3') into ( Gff3_annotation_for_filtering, Gff3_annotation_for_bed_files)

	shell:
	'''

	source !{params.bin_dir}/check_checksum_ensembl.sh
	species='!{species}'
	species_long='!{species_long}'
	genome='!{genome}'
	assembly='!{assembly}'
	release='!{release}'

	base_url=ftp://ftp.ensembl.org/pub/release-$release

	gff3_file=${species_long^}.$genome.$release.gff3.gz
	fasta_file=${species_long^}.$genome.dna_sm.$assembly.fa.gz

	url_anno=$base_url/gff3/$species_long
	wget -O annotation.gff3.gz $url_anno/$gff3_file 
	wget -O checksums_anno $url_anno/CHECKSUMS
	check_checksum_and_gunzip annotation.gff3.gz $gff3_file checksums_anno

	url_fasta=$base_url/fasta/$species_long/dna
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





process making_bowtie2_indexes {
  tag "${species}"

  label "bowtie2_samtools"

  publishDir path: "${"${species}/genome/sequence"}/bowtie2_indexes", mode: 'link'

  input:
    set species, file(gff3), file(fasta) from Genome_and_annotation_1

  output:
		file("*.bt2")

  shell:
  '''
		bowtie2-build --threads !{number_of_cores} !{fasta} genome
  '''

}


process indexing_genomes {
  tag "${species}"

  label "bowtie2_samtools"

  publishDir path: "${"${species}/genome/sequence"}", mode: 'link'

  input:
    set species, file(gff3), file(fasta) from Genome_and_annotation_2

  output:
		set species, file('genome.fa.fai') into Fasta_fai

  shell:
  '''
		samtools faidx !{fasta}
  '''

}

Genome_and_annotation_3
	.join(Fasta_fai)
	.set{ Fasta_fai_gff3 }





process making_transcriptome {
  tag "${species}"

	label "gffread"

	publishDir path: "${species}/genome/sequence", mode: 'link'

  input:
    set species, file(gff3), file(fasta), file(fasta_indexes) from Fasta_fai_gff3

  output:
		set species, file('transcriptome.fa') into Transcriptome_for_building_kallisto_indexes

  shell:
  '''

      gffread -C !{gff3} -g !{fasta} -w transcriptome.fa

  '''

}

// the -C argument tells gffread to filter the gff3 to keep only protein coding genes




process filtering_annotation_file {
  tag "${species}"

	label "ubuntu"

	publishDir path: "${species}/genome/annotation/filtered", mode: 'link'

  input:
    set species, file(gff3) from Gff3_annotation_for_filtering

  output:
		set species, file('chromosomes_sizes.txt'), file(gff3), file('protein_coding_genes.gff3'), file('annotation_filtered_regions_and_pseudogenes.gff3') into GFF3_filtered_for_R
		file("*.txt")

  shell:
  '''

			gff3=!{gff3}
			source !{params.bin_dir}/get_tab.sh

			tab=$(get_tab)

			gff3_filtered=annotation_filtered_regions_and_pseudogenes.gff3
			# savings statistics on the unfiltered annotation file
			cut -f1 ${gff3} | sort | uniq -c | grep -v "#" | sort -nr > nb_of_feature_by_region__raw_gff3_file.txt
			awk -v FS="$tab" -v OFS="$tab" '{if($3 == "gene"){print $1}}' ${gff3} | sort | uniq -c | sort -nr	> nb_of_genes_by_region__raw_gff3_file.txt 
			cut -f3 ${gff3} | sort | uniq -c | sort -nr | grep -v "#" > nb_of_feature_by_feature_type__raw_gff3_file.txt

			# keeping only regions (contigs and chromosomes) with at least 5 genes, excluding the mitochondrial chromosome and removing all genes that are not encoding for proteins
			regions_to_keep=$(awk -v OFS="$tab" '{if($1 > 4) print $2}' nb_of_genes_by_region__raw_gff3_file.txt | sort | paste -sd ' ')
			awk -v FS="$tab" -v OFS="$tab" -v rtk="$regions_to_keep" 'BEGIN{split(rtk, rtk1, " ")} {for (c1 in rtk1) {if($1 == rtk1[c1] && $1 != "mitochondrion_genome" && $1 != "MT" && $1 != "MtDNA" && $3 !~ "ncRNA" && $3 !~ "pseudogen" && $3 !~ "transposable_element") print}}' ${gff3} > ${gff3_filtered}
			
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


//// outputs:
// - protein_coding_genes.gff3 is used to create the gene metadata table for all analysis
// - anno_without_pseudogenes.gff3 is used for creating a txdb database without non coding transcripts, for the annotation of ATAC seq peaks

//// mitochondrial chromosomes:
// worm  |         fly          | mouse  | human  
// MtDNA | mitochondrion_genome |   MT   |  MT

// scaffolds
// fly: 2110000*
// human: GL000*, KI270*

// note: the important thing is to keep all contigs when mapping reads to prevent misalignement of reads to wrong location. However, discarding reads mapped to contigs after mapping is a matter of taste and depends on the downstream application. In our case, we need to avoid small contigs to avoid erroneous annotations of peaks to the closest gene. So we discard small contigs (current threshold: we keep only contigs with 5 genes or more).

// awk '{if($3 == "gene"){print $1}}' mouse/genome/annotation/annotation.gff3 | sort | uniq -c | sort -nr | awk '{if ($1 >= 5) {print $2}}' - | wc -l # 30
// awk '{if($3 == "gene"){print $1}}' human/genome/annotation/annotation.gff3 | sort | uniq -c | sort -nr | awk '{if ($1 >= 5) {print $2}}' - | wc -l # 25



process making_bed_files_of_annotated_regions {
	tag "${species}"

	label "bedops"

	publishDir path: "${species}/genome/annotation/bed_regions", mode: 'link'

	input:
		set species, file(gff3) from Gff3_annotation_for_bed_files

	output:
		file('*.bed')

	shell:
	'''

			gff3=!{gff3}

			# generating a gff3 file without chromosomes 
			awk -v FS="\t" -v OFS="\t" '{if (substr($1, 1, 1) != "#" && $3 != "chromosome" && $3 != "scaffold" && $3 != "region" && $3 != "biological_region") print}' $gff3 | sort -k1,1 -k4,4n > annotation_clean.gff3
			gff2bed < annotation_clean.gff3 > annotation_clean.bed

			# Get already defined regions: exons, genes
			awk -v FS="\t" -v OFS="\t" '{if ($3 == "exon") print $1, $4-1, $5}' annotation_clean.gff3 > exons.bed
			awk -v FS="\t" -v OFS="\t" '{if ($3 == "gene") print $1, $4-1, $5}' annotation_clean.gff3 > genes.bed
			# Get chromosomes sizes
			awk -v FS="\t" -v OFS="\t" '{if ($3 == "chromosome" || $3 == "scaffold" || $3 == "region") {print $1, $5}}' ${gff3} | sort -k1,1 > chromosomes_sizes.txt
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

			# deleting intermediate files
			rm genes_detailed.bed promoters_forward.bed promoters_reverse.bed annotation_clean.bed chromosomes_sizes.bed

	'''
}

// In this process, we make bed files for introns, exons, intergenic regions and promoters. 
// A gff3 file without chromosomes entries is used to create intergenic regions using the complement between chromosome sizes and all annotated features in the genome (excluding "biological_regions"; i.e., Predicted TSS, CpG islands, fosmids).
// Introns are defined as the complement of exonic and intergenic regions.
 
// main links
// https://www.biostars.org/p/112251/#314840 => get intergenic, exons and introns .bed
// https://www.biostars.org/p/360842/#9531161 // => get promoter.bed

// other links
// https://github.com/TriLab-bioinf/Script_wrappers/blob/main/get_intergenic_bed_from_gtf.sh
// https://davetang.org/muse/2013/01/18/defining-genomic-regions/
// https://bedparse.readthedocs.io/en/stable/Tutorial.html // => could not make it work

// scaffolds should be kept
// https://www.biostars.org/p/223541/
// The scaffolds are needed, they are likely parts of the reference genome which could not yet in this assembly build be confidently assigned to chromosomes. This is typical of repeat containing scaffolds. They are likely to be small and not gene rich, but excluding them would cause short reads which should map to them to be forced to map to the chromosomes and cause problemes downstream, such as false SNPs etc.

// Here are the lines used to create the get_tab.sh script (creating within bash prevent the issue with Windows line-breaks)
// cat > get_data/bin/get_tab.sh << EOL
// #!/usr/bin/env bash
// 
// get_tab () { echo "\t" ; }
// 
// EOL


process making_kallisto_indexes {
  tag "${species}"

  label "kallisto"

  publishDir path: "${"${species}/genome/sequence"}", mode: 'link'

  input:
		set species, file(transcriptome) from Transcriptome_for_building_kallisto_indexes

  output:
    file('transcriptome_kallisto_index')

  script:
  """
      kallisto index --make-unique -i transcriptome_kallisto_index ${transcriptome}
  """

}




process getting_orgdb {
  tag "${species}"

  label "annotationhub"

  publishDir path: "${species}/genome/annotation/R", mode: 'link'

  input:
		set species, species_long from Ensembl_Assembly_1.orgdb

  output:
    set species, species_long, file('orgdb.sqlite') into Orgdb_for_annotation

  shell:
  '''
		#!/usr/bin/env Rscript

		library(AnnotationHub)
		species = '!{species}'
		species_long = '!{species_long}'
		ncbi_orgdb_version = '!{ncbi_orgdb_version}'

		species_initial =  sapply(strsplit(species_long, '_')[[1]], substr, 1, 1) %>% {paste0(toupper(.[1]), .[2])}

		ah = AnnotationHub(localHub = F, cache = '.')
		
		pattern = paste0('ncbi/standard/', ncbi_orgdb_version, '/org.', species_initial, '.eg.db.sqlite')
		orgdb = query(ah, pattern)[[1]]
		AnnotationDbi::saveDb(orgdb, 'orgdb.sqlite')

  '''

}

// annotationhub_cache = '!{params.annotationhub_cache}'
// ah = AnnotationHub(cache = annotationhub_cache)

// if cache is corrupted, delete it so it will be rebuilt automatically
// rm $cactus_dir/data/util/annotationhub_cache/*

Orgdb_for_annotation
	.join(GFF3_filtered_for_R)
	.set{GFF3_filtered_for_R_1}

process making_R_annotation_files {
  tag "${species}"

  label "bioconductor"

  publishDir path: "${species}/genome/annotation/R", mode: 'link'

  input:
		set species, species_long, orgdb, file(chr_size), file(gff3_raw), file(gff3_genes_only), file(gff3_filtered_regions_and_pseudogenes) from GFF3_filtered_for_R_1

  output:
    file('*')

  shell:
  '''
		#!/usr/bin/env Rscript

		library(magrittr)
		library(AnnotationDbi)

		source('!{params.cactus_dir}/bin/export_df_to_bed.R')
		source('!{params.bin_dir}/get_tab.R')

		gff3_raw = '!{gff3_raw}'
		gff3_genes_only = '!{gff3_genes_only}'
		gff3_filtered_regions_and_pseudogenes = '!{gff3_filtered_regions_and_pseudogenes}'
		species = '!{species}'
		orgdb = AnnotationDbi::loadDb('!{orgdb}')
		df_chr_size = read.table('!{chr_size}', sep = get_tab())
		species_long = '!{species_long}'


		# exporting annotation dataframe
		anno_df = rtracklayer::readGFF(gff3_genes_only) 
		entrez_id = mapIds(orgdb, keys = anno_df$gene_id,  column = 'ENTREZID', keytype = 'ENSEMBL', multiVals = 'first')
		anno_df %<>% dplyr::rename(gene_name = Name)
		anno_df %<>% dplyr::mutate(width = end - start)
		anno_df$entrez_id = entrez_id
		anno_df$chr = as.character(anno_df$seqid)
		anno_df %<>% dplyr::select(chr, start, end, width, strand, gene_name, gene_id, entrez_id)
		saveRDS(anno_df, 'df_genes_metadata.rds')

		# creating the seqinfo object
		nr = nrow(df_chr_size)
		seqinfo = GenomeInfoDb::Seqinfo(seqnames = df_chr_size[,1], seqlengths = df_chr_size[,2], isCircular = rep(F, nr) , genome = rep(species_long, nr))
		saveRDS(seqinfo, 'seqinfo.rds')

		# exporting txdb
		txdb = GenomicFeatures::makeTxDbFromGFF(gff3_filtered_regions_and_pseudogenes, chrominfo = seqinfo)
		AnnotationDbi::saveDb(txdb, file="txdb.sqlite")

		# exporting gene vs transcripts
		txdb1 = GenomicFeatures::makeTxDbFromGFF(gff3_raw)
		df_genes_transcripts = select(txdb1, keys = keys(txdb1), columns = c('GENEID', 'TXNAME', 'TXCHROM'), keytype = 'GENEID')
		saveRDS(df_genes_transcripts, 'df_genes_transcripts.rds')

		# extracting promoters
		promoters = GenomicFeatures::promoters(GenomicFeatures::genes(txdb), upstream = 1500, downstream = 500)

		# exporting promoters as dataframe
		promoters_df = as.data.frame(promoters, stringsAsFactors = F)
		promoters_df$start[promoters_df$start < 0] = 0
		promoters_df$end[promoters_df$end < 0] = 0
		promoters_df %<>% dplyr::rename(gene_name = gene_id)
		colnames(promoters_df)[colnames(promoters_df) == 'gene_name'] = 'gene_id'
		promoters_df %<>% dplyr::inner_join(anno_df %>% dplyr::select(gene_id, gene_name, entrez_id), by = 'gene_id')
		promoters_df %<>% dplyr::arrange(seqnames, start)
		promoters_df %<>% dplyr::rename(chr = seqnames)
		saveRDS(promoters_df, file = 'promoters_df.rds')

		# exporting promoters as bed file
		promoters_df1 = promoters_df
		promoters_df1$score = 0
		promoters_df1 %<>% dplyr::select(chr, start, end, gene_name, score, strand, gene_id)
		export_df_to_bed(promoters_df1, 'promoters.bed')

		# exporting kegg data
		library(clusterProfiler)
		species_split = strsplit(species_long, '_')[[1]]
		kegg_code = paste0(substr(species_split[[1]], 1, 1), substr(species_split[[2]], 1, 2))
		options(clusterProfiler.download.method = 'wget')
		kegg_environment = clusterProfiler:::prepare_KEGG(kegg_code, 'KEGG', "ncbi-geneid")
		saveRDS(kegg_environment, 'kegg_environment.rds')


  '''

}

//// inputs:
// - gff3_raw (i.e., annotation.gff3) is the raw gff3 file. It is used to create the full gene vs transcripts table. All transcripts should be there as this table is used by Kallisto for alignment free transcripts quantification
// - gff3_genes_only (i.e., protein_coding_genes.gff3) only include entries with feature type being gene. It is used to create the df_genes_metadata table that contains coordinates, gene name and gene ID for Ensembl and NCBI. This file is used in multiple occasions in Cactus for correspondance between gene_id and gene name (volcano plots), for adding coordinate informations (for the res_detailed_atac tables) or for correspondance between Ensembl and NCBI ids (KEGG enrichment).
// - gff3_filtered_regions_and_pseudogenes (i.e., anno_without_pseudogenes.gff3) contains all features types excepting ncRNA_gene. It is used for creating a txdb database without non coding transcripts. The resulting txdb file omits all features that don't have a parental gene. Therefere, all transcripts and exons of ncRNA_genes (i.e., miRNA, lnc_RNA, rRNA, tRNA, ...) are excluded from the txdb file. This txdb file is used for the annotation of ATAC-Seq peaks, but also for the creation of promoter files used in the mRNA-Seq analysis.


// For enrichGO, the go terms are stored into a package and accessed liked that:  goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
// so everything is run in local within the clusterProfiler:::get_GO_data function
// for enrichKEGG, the KEGG database is fetched online. We use a preparsed one here instead
// using download.method = 'wget' is necessary as curl is not in the container. Loading the clusterProfiler at the end is necessary otherwise it interfere with some AnnotationDbi functions


// unfortunately, the "promoters function" extract promoters with coordinate below zero (I found this bug for 3 mitochondrial genes in C elegans). Thus we need to put them at zero.
// promoters_df %<>% dplyr::select(chr:strand, gene_name, gene_id, entrez_id)
// promoters_df %<>% dplyr::mutate(seqnames = as.character(seqnames), strand = as.character(strand))


///// Checking if the -C option of gffread really did filter all non-coding transcripts, and if we have the same set of genes as in my manually filtered gff3 files (here in worms)

//// in the appropriate worm work folders
// library(magrittr)
// pnrow <- function(x) print(nrow(x))
// df_genes_transcripts = readRDS('df_genes_transcripts.rds') %T>% pnrow # 44264
// df_genes_metadata = readRDS('df_genes_metadata.rds') %T>% pnrow # 19985
// length(unique(df_genes_transcripts$GENEID)) # 31571
// lw(df_genes_metadata$gene_id %in% df_genes_transcripts$GENEID)
// 
// fastaFiles  = Biostrings::readDNAStringSet("~/workspace/cactus/data/worm/genome/sequence/transcriptome.fa")
// seq_name = names(fastaFile)
// sequence = paste(fastaFile)
// df <- data.frame(seq_name, sequence)
// df1 = df[,1]
// v_genes_in_transcriptome =  df1 %>% gsub('.*gene=', '', .) %>% gsub(' .*', '', .)
// v_genes_in_transcriptome %>% length # 31768
// v_genes_in_transcriptome %>% unique %>% length # 19997
// 
// v_missing = unique(v_genes_in_transcriptome) %>% .[!. %in% df_genes_metadata$gene_name]
//  <!-- [1] "WBGene00022116" "WBGene00021837" "WBGene00019083" "WBGene00023413"
//  [5] "WBGene00010654" "WBGene00010164" "WBGene00021985" "WBGene00021822"
//  [9] "WBGene00021250" "WBGene00021916" "WBGene00021803" "nduo-6"
// [13] "WBGene00010958" "WBGene00010959" "atp-6"          "nduo-2"
// [17] "ctb-1"          "ctc-3"          "nduo-4"         "ctc-1"
// [21] "ctc-2"          "nduo-3"         "nduo-5" -->
// 
// 
// => there are minor differences. So the -C option of gff3 does a good job. However, it is good to still make our own gff3-derived gene table to make sure we work on the same genes in ATAC-Seq and mRNA-Seq (since for instance mitochondrial transcripts are removed)
// 
// 
// df_gff3 = rtracklayer::readGFF('annotation.gff3') %>% data.frame
// df_gff3[df_gff3$Name %in% v_missing, ][, c('seqid', 'Name')]
//        <!-- seqid   Name
// 407022 MtDNA nduo-6
// 407054 MtDNA  atp-6
// 407068 MtDNA nduo-2
// 407086 MtDNA  ctb-1
// 407094 MtDNA  ctc-3
// 407102 MtDNA nduo-4
// 407107 MtDNA  ctc-1
// 407125 MtDNA  ctc-2
// 407136 MtDNA nduo-3
// 407141 MtDNA nduo-5
// df_gff3[df_gff3$ID %in% v_missing, ][, c('seqid', 'Name')]
//        seqid Name
// 12852      I <NA>
// 13596      I <NA>
// 30397      I <NA>
// 49362      I <NA>
// 70404      I <NA>
// 75048      I <NA>
// 99904     II <NA>
// 110505    II <NA>
// 189304   III <NA>
// 194351   III <NA>
// 195562   III <NA>
// 407027 MtDNA <NA>
// 407050 MtDNA <NA>








process making_timestamp {
  tag "${species}"

  // label "r_basic"

  publishDir path: "${species}", mode: 'link'

	cache false

  input:
	  val species from Channel.of( 'worm', 'fly', 'mouse', 'human' )

  output:
    file('timestamp.txt')

  shell:
  '''
		#!/usr/bin/env Rscript
		
		container_engine = '!{workflow.containerEngine}'
		
		library(magrittr)
		
		wl <- function(x) writeLines(x)
		wl1 <- function(x) writeLines(paste0(x, "\n"))
		wl2 <- function(x) writeLines(paste0("\n", x, "\n"))
		
		sink("timestamp.txt")
		
		  paste0("This reference was built on ", Sys.Date(), ".") %>% wl2
		
		  "The version of the workflow and package manager tools used are:" %>% wl
		
		  system("nextflow -version", intern = T) %>% wl
		  system(paste(container_engine, "--version"), intern = T) %>% wl1
		
		sink()

	'''
}



// on completion
 workflow.onComplete
 {
  println ""
  println "Workflow completed on: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'Succeeded' : 'Failed' }"
  println "Workflow Duration: $workflow.duration"
  println ""
 }
