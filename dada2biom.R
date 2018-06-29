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
	parser$add_argument("-m", "--metadata", default=NULL, type="character",help="Name of (or path to) the metadata, samples as rows")
	parser$add_argument("-s", "--samp_col", default=NULL, type="interger",help="Numer of column that matches the sample name used in the ASV")
	parser$add_argument("output",default="biom.json",type="character",help="Name of (or path to) the output BIOM file.")

args = parser$parse_args(c('~/lab/genelab/dada2_cli/data/asv.tsv'))
if(!is.null(args$metadata) && !file.exists(args$metadata)) stop('The metadata file does not exist')
if(!file.exists(args$asv_table)) stop('The asv table file does not exist')
if(!file.exists(args$taxa_table)) stop('The taxa table file does not exist')

library(biom)

# read in all the files
metadata=read.table(args$metadata)
asv=read.table(args$asv_table)
taxa=read.table(args$taxa_table)

if(ncol(asv)==nrow(taxa)) stop('The length of the ASV does not match the taxa')
if(!is.null(args$metadata) && row(meta)==nrow(asv)) stop('The number of samples in the ASV does not match the metadata')

# sort the input
asv <- asv[, order(colnames(asv))]
taxa <- taxa[order(rownames(taxa)),]
# If the orders do not match, stop. 
if(any(colnames(asv)!=rownames(taxa))) stop(sprintf('After sorting the ASV and TAXA a mismatch occured: \nASV:%s\nTAXA:%s', asv[colnames(asv)!=rownames(taxa)], taxa[colnames(asv)!=rownames(taxa)]))  # Count sequences not in same order as Taxa seq 

# sort metadata
metadata = metadata[ order(metadata[[args$samp_col]]),]
if(any(rownames(asv)!=metadata[[args$samp_col]])) stop('After sorting ASV samplename do not match Metadata sample names')

# Create the biom object from data
the_biom=make_biom(t(asv),observation_metadata = taxa, sample_metadata = metadata )

# And write it out
write_biom(the_biom, args$output)
