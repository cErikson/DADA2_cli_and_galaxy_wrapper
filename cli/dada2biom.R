#!/usr/bin/env Rscript
################################################################################
# Title: dada2biom.R
# Discription: Convert the DADA2 output files to a Biom file
# Author: Christian Erikson/Nicholas Bense
# GitHub: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
# Date: 6/26/18
################################################################################
packages = c("biom", "argparse")
if (any(!(packages %in% installed.packages()[,"Package"]))) stop(sprintf('These packages are required: %s', packages[!(packages %in% installed.packages()[,"Package"])]))
library("argparse")


parser <- ArgumentParser()


	parser$add_argument("asv_table",type="character",help="Name of (or path to) the ASV table output from DADA2. In R matrix format. aka first feild of col names is missing")
	parser$add_argument("taxa_table",type="character",help="Name of (or path to) the taxa table file from assign_tax")
	parser$add_argument("-m", "--metadata", default=FALSE, type="character",help="Name of (or path to) the metadata, samples as rows")
	parser$add_argument("-s", "--samp_col", default=NULL, help="Name of column in the metadata that uniquely matches part of the sample name used in the ASV")
	parser$add_argument("output",default="biom.json",type="character",help="Name of (or path to) the output BIOM file.")
	parser$add_argument("-R","--regex", action='store',type="character",default=FALSE,help="Regex to exctract names form ASV and exact match to metadata. Double escapes \\w")
	
args = parser$parse_args()
if(any(args$metadata != F) && !file.exists(args$metadata)) stop('The metadata file does not exist')
if(!file.exists(args$asv_table)) stop('The asv table file does not exist')
if(!file.exists(args$taxa_table)) stop('The taxa table file does not exist')

library(biom)

# read in all the files
if(args$metadata != F) metadata=read.delim(args$metadata, header = T)
asv=read.table(args$asv_table)
taxa=read.table(args$taxa_table)

# Check if the samp_col is numeric or text 
if (is.na(as.numeric(args$samp_col)) && args$metadata != F){
	if(!args$samp_col %in% colnames(metadata)){
		write(colnames(metadata), stderr())
	 	stop('The samp_col name was not found in metadata, above names were found in the file')
	}
} else {
	args$samp_col = as.numeric(args$samp_col)
	warning('The samp_col is being interpreted as a column number selection')
}

if(ncol(asv)!=ncol(asv)) stop('The length of the ASV does not match the taxa')
if(args$metadata != F && nrow(metadata)!=nrow(asv)) warning('The number of samples in the ASV does not match the metadata')

# sort the input to match taxa
asv <- asv[, order(colnames(asv))]
taxa <- taxa[order(rownames(taxa)),]

if(args$metadata != F){
	ord=numeric()
	if (args$regex != F){
		library(stringr)
		asv_names = str_match(rownames(asv), args$regex)[,-1] #extract using regex
		for (x in metadata[[args$samp_col]]){
			loc=which(x == asv_names)
			if (length(loc)>1) stop(sprintf('The name %s matched multipule asv rows in the asv table, a row that matched: %s\n',x, rownames(asv)[loc]))
			ord=c(ord,loc)
		}
	} else {
		for (x in metadata[[args$samp_col]]){
			loc=grep(x, rownames(asv))
			if (length(loc)>1) stop(sprintf('The name %s matched multipule asv rows in the asv table, a row that matched: %s\n',x, rownames(asv)[loc]))
			ord=c(ord,loc)
			}
		}
	
	if(any(duplicated(ord))) stop(sprintf('The metadata column did not uniquely match the sample names: %s\n',metadata[[args$samp_col]][ord[duplicated(ord)]]))
	asv = asv[ ord, ]
	if(nrow(metadata)!=nrow(asv)) stop(sprintf('After selecting the metadata based on %s, the number of samples in the ASV does not match the metadata',args$samp_col))
	
}
# If the orders do not match, stop. 
if(any(colnames(asv)!=rownames(taxa))) stop(sprintf('After sorting the ASV and TAXA a mismatch occured: \nASV:%s\nTAXA:%s\n', asv[colnames(asv)!=rownames(taxa)], taxa[colnames(asv)!=rownames(taxa)]))  # Count sequences not in same order as Taxa seq 
if(args$metadata != F && any(!rownames(asv)!=metadata[[args$samp_col]])) stop('After sorting ASV samplename do not match Metadata sample names')

# Create the biom object from data

if(any(args$metadata != F)){
	the_biom=make_biom(t(asv),observation_metadata = taxa, sample_metadata = metadata )
}else{
	the_biom=make_biom(t(asv),observation_metadata = taxa )}

# And write it out
write_biom(the_biom, args$output)
# 
