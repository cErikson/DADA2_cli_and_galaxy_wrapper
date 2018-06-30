#!/usr/bin/env Rscript
################################################################################
# Title: dada2_assign_taxa.R
# Discription: A wrapper for DADA2 taxonomic assignment
# Author: Christian Erikson
# GitHub: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
# Date: 6/26/18
################################################################################
packages = c("dada2", "argparse")
if (any(!(packages %in% installed.packages()[,"Package"]))) stop(sprintf('These packages are required: %s', packages[!(packages %in% installed.packages()[,"Package"])]))
library(argparse)

parser <- ArgumentParser()
	parser$add_argument("asv_table",type="character", help="Name of (or path to) the ASV table output from DADA2. In R matrix format. aka first feild of col names is missing")
	parser$add_argument("output",default="taxa.tav",type="character",help="Name of (or path to) the output taxa table.")
	parser$add_argument("-t", "--taxa_train", nargs='+', type="character", help="The path to tone or more taxa trainning file(s)")
	parser$add_argument("-s","--species_train",default=FALSE,nargs='+', type="character", help="The path to  species training file")

	parser$add_argument('--minBoot', default= T, help='(Optional). Default 50. The minimum bootstrap confidence for assigning a taxonomic level.')
	parser$add_argument('--tryRC', default= F, help='(Optional). Default FALSE. If TRUE, the reverse-complement of each sequences will be used for classification if it is a better match to the reference sequences than the forward sequence.')
	parser$add_argument('--outputBootstraps', default= F, help='(Optional). Default FALSE. If TRUE, bootstrap values will be retained in an integer matrix. A named list containing the assigned taxonomies (named "taxa") and the bootstrap values (named "boot") will be returned. Minimum bootstrap confidence filtering still takes place, to see full taxonomy set minBoot=0')
	parser$add_argument('--taxLevels', default=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") , help='(Optional). Default is <Kingdom Phylum Class Order Family Genus Species>. The taxonomic levels being assigned. Truncates if deeper levels not present in training fasta.')
	parser$add_argument('--allowMultiple', default= F, help='(Optional). Default FALSE. Defines the behavior when multiple exact matches against different species are returned. By default only unambiguous identifications are return. If TRUE, a concatenated string of all exactly matched species is returned. If an integer is provided, multiple identifications up to that many are returned as a concatenated string.')
	parser$add_argument('--multithread', default= T, help='(Optional). Default is FALSE. If TRUE, multithreading is enabled and the number of available threads is automatically determined. If an integer is provided, the number of threads to use is set by passing the argument on to setThreadOptions.')
	parser$add_argument('--verbose', default= T, help='(Optional). Default FALSE. If TRUE, print status to standard output.')

args = parser$parse_args()
if (args$species_train != F && any(!file.exists(args$species_train))) stop("The species_trainning does not exist")
if (any(!file.exists(args$taxa_train))) stop("The taxa_trainning does not exist")
if (any(!file.exists(args$asv_table))) stop(" The ASV_table does not exist")

library(dada2)

# Read in seqtab file
seqtab=as.matrix(read.table(args$asv_table))

# assignTaxonomy
taxa = assignTaxonomy(seqtab, args$taxa_train, minBoot = args$minBoot, tryRC = args$tryRC, outputBootstraps = args$outputBootstraps,
					  taxLevels = args$taxLevels, multithread=args$multithread, verbose = args$verbose) # main training set
if (args$species_train!=F) taxa = addSpecies(taxa, args$species_train, allowMultiple = args$allowMultiple, verbose = args$verbose) # Add species 

# Save
write.table(taxa, args$output)
