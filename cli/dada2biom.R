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
	parser$add_argument("-s", "--samp_col", default=NULL, help="Name or number of column in the metadata that uniquely matches part of the sample name used in the ASV")
	parser$add_argument("output",default="biom.json",type="character",help="Name of (or path to) the output BIOM file.")

args = parser$parse_args()
if(any(args$metadata != F) && !file.exists(args$metadata)) stop('The metadata file does not exist')
if(!file.exists(args$asv_table)) stop('The asv table file does not exist')
if(!file.exists(args$taxa_table)) stop('The taxa table file does not exist')

print(args$samp_col)
library(biom)

# read in all the files
if(args$metadata != F) metadata=read.delim(args$metadata)
asv=read.table(args$asv_table)
taxa=read.table(args$taxa_table)

if(!args$samp_col %in% colnames(metadata)) warning('The samp_col name was not found in metadata, ignore if selecting by number')
if(ncol(asv)!=ncol(asv)) stop('The length of the ASV does not match the taxa')
if(args$metadata != F && nrow(metadata)!=nrow(asv)) warning('The number of samples in the ASV does not match the metadata')

# sort the input
asv <- asv[, order(colnames(asv))]
asv <- asv[order(rownames(asv)),]
taxa <- taxa[order(rownames(taxa)),]


if(args$metadata != F){
	ord=numeric()
	for (x in metadata[[args$samp_col]]){
		loc=grep(x, rownames(asv))
		if (length(loc)>1) stop(sprintf('The column values are not unique to the asv table, multiple rows were matched: %s',rownames(asv)[loc]))
		ord=c(loc, ord)
	}
	if(any(duplicated(ord))) stop(sprintf('The metadata column did not uniquely match the sample names: %s\n',metadata[[args$samp_col]][ord[duplicated(ord)]]))
	metadata = metadata[ ord, ]
	if(nrow(metadata)!=nrow(asv)) stop(sprintf('After selecting the metadata based on %s, the number of samples in the ASV does not match the metadata',args$samp_col))
	
}
# If the orders do not match, stop. 
if(any(colnames(asv)!=rownames(taxa))) stop(sprintf('After sorting the ASV and TAXA a mismatch occured: \nASV:%s\nTAXA:%s\n', asv[colnames(asv)!=rownames(taxa)], taxa[colnames(asv)!=rownames(taxa)]))  # Count sequences not in same order as Taxa seq 
#if(args$metadata != F && any(rownames(asv)!=metadata[[args$samp_col]])) stop('After sorting ASV samplename do not match Metadata sample names')

# Create the biom object from data

if(any(args$metadata != F)){
	the_biom=make_biom(t(asv),observation_metadata = taxa, sample_metadata = metadata )
}else{
	the_biom=make_biom(t(asv),observation_metadata = taxa )}

# And write it out
write_biom(the_biom, args$output)
# 
