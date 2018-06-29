#!/usr/bin/env Rscript
##	++=======================================================++
##	|| DADA2 taxnonomy script for marker gene GLDS datasets  ||
##	++=======================================================++
##
##	Christian Erikson: Bio.Erkson@gmail.com christian.b.erikson@nasa.gov
##
library(argparse)

parser <- ArgumentParser()
	parser$add_argument("asv_table",type="character", help="Name of (or path to) the ASV table output from DADA2. In R matrix format. aka first feild of col names is missing")
	parser$add_argument("output",default="taxa.tav",type="character",help="Name of (or path to) the output BIOM file.")
	parser$add_argument("taxa_train",type="character", help="The path to the taxa trainning file")
	parser$add_argument("-s","--species_train",type="character", help="The path to the species training file")
	parser$add_argument("-d","--add_taxa_ds",help="Default: 'Sample Name'.The name or feild for the the sample name in the isa file ")
	
opt = parser$parse_args()

if (file.access(opt$taxa_train)==-1 || file.access(opt$asv_table)==-1){
	print_help(opt_parser)
	stop("asv_table, taxa_train must be supplied and exist", call.=FALSE)
}
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))


# Read in seqtab file
seqtab=as.matrix(read.table(opt[['asv_table']]))

# assignTaxonomy
taxa = assignTaxonomy(seqtab, opt[['taxa_train']], multithread=TRUE, verbose = TRUE) # main training set
if (!is.null(opt$species_train)) taxa = addSpecies(taxa, opt[['species_train']]) # Add species 
taxa = add_rownames(as.data.frame(taxa), var = 'ASV')
if (!is.null(opt[['add_taxa_ds']][[1]])){
	for (i in opt[['add_taxa_ds']][[1]]){
		add_taxa = assignTaxonomy(seqtab, i, multithread=TRUE, verbose = TRUE)
		add_taxa = add_rownames(as.data.frame(add_taxa), var = 'ASV')
		taxa = bind_rows(taxa, add_taxa) %>%
			distinct(., ASV)
	}
}

# Save
write_tsv(taxa, opt[['output']])
