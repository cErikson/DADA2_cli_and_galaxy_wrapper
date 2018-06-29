#!/usr/bin/env Rscript
##	++==============================================++
##	|| DADA2 to BIOM for marker gene GLDS datasets  ||
##	++==============================================++
##
##	Christian Erikson: Bio.Erkson@gmail.com christian.b.erikson@nasa.gov
##

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

	parser$add_argument("-i","--isa_dir",type="character",help="Name of (or path to) the sample isa directory")
	parser$add_argument("-a","--asv_table",type="character",help="Name of (or path to) the ASV table output from DADA2. In R matrix format. aka first feild of col names is missing")
	parser$add_argument("-t","--taxa_table",type="character",help="Name of (or path to) the taxa table file from assign_tax")
	parser$add_argument("-s","--isa_sample_file",default=1,help="Deafult: 1st file. The name of the ISA sample file corresponding to the data")
	parser$add_argument("-f","--isa_sample_feild",default="Sample Name",help="Default: 'Sample Name'.The name or feild for the the sample name in the isa file ")
	parser$add_argument("-o","--output",default="biom.json",type="character",help="Name of (or path to) the output BIOM file.")

opt = parser$parse_args()

if (!dir.exists(opt$isa_dir) || file.access(opt$asv_table) == -1 || file.access(opt$taxa_table) == -1){
	print_help(opt_parser)
	stop("isa_dir, asv_table, taxa_table must be supplied and acessable", call.=FALSE)
}

suppressPackageStartupMessages(library(biom))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(Risa))
suppressPackageStartupMessages(library(stringr))


# read in all the files
isa=readISAtab(path = opt$isa_dir)
counts=read.table(opt$asv_table)
taxa=read_tsv(opt$taxa_table)

# HaCk away the duplicate columns, The isa standard needs to have unique columns. Debuging this error was cryptic and took forever
dedup_isa=isa['study.files'][[opt[['isa_sample_file']]]][, !duplicated(colnames(isa['study.files'][[opt[['isa_sample_file']]]]))]  ##!!!!! This is a hack, The isa file has duplicate col names for unit, ref protocol, and terms. The isa file standard needs to have unique column names 
dedup_isa=dedup_isa[order(dedup_isa[opt[['isa_sample_feild']]]),]

# If the orders do not match, stop. 
stopifnot(colnames(counts)==taxa$ASV) # Count sequences not in same order as Taxa seq 
stopifnot(str_split_fixed(rownames(counts), '_', 2)[,1]==dedup_isa[opt[['isa_sample_feild']]]) # count sample names not in same order as isa study sample names

# Create the biom object from data
the_biom=make_biom(t(counts),observation_metadata = taxa, sample_metadata = dedup_isa )

# And write it out
write_biom(the_biom, opt[['output']])
