#!/usr/bin/env Rscript
################################################################################
# Title: Dada2_ASVTable
# Author: Christian Erikson/Nicholas Bense
# GitHub: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
# Date: 6/26/18
# Description: R script to create Amplicon Sequence Variant Tables for 16s and ITS 
# microbial diversity surveys. These can be used for taxonomic alignment downstream.
################################################################################

# Import packages
packages <- c("dada2", "argparse", "stringr")
if (any(!(packages %in% installed.packages()[,"Package"]))) stop(sprintf('These packages are required: %s', packages[!(packages %in% installed.packages()[,"Package"])]))
library(dada2)
library(argparse)


##### PARSER #####
# Initialize argument parser
parser <- ArgumentParser(description='A Command line wrapper for DADA2.')

# Read arguments/parameters
parser$add_argument("-f","--fwd", nargs='+', type="character",help="Name/path to fastq or fastq.gz files for single end sequencing or forward reads if paired end sequencing")
parser$add_argument("-r","--rev", nargs='+', type="character",default=FALSE,help="Name/path to fastq or fastq.gz files for reverse reads if paired end sequencing")
parser$add_argument("-o","--asv_out",type="character", default="Dada2_ASVTable.tsv",help="Name/path to ASV table file output")
parser$add_argument("-c","--chimrm",type="character", default='consensus', help="Method to use in removeBimeraDenovo. `pooled`, `consensus`, `per-sample`, `FALSE` to skip. Default consensus")
parser$add_argument("-s","--samp_fields", nargs='+', type="integer", default=FALSE, help="The fields in the file name that should be used for sample names.")
parser$add_argument("-S","--fields_delim", type="character", default='_', help="The field delimiter used in the file name")
parser$add_argument("-R","--samp_regex", action='store',type="character",default=FALSE,help="Regex used to create samplenames. Double escapes \\w")
parser$add_argument("-l","--samp_list", nargs='+', action='store',type="character", default=FALSE, help="List of sample names that are in same order as the reads, can be used in combination with -R or -s")
##
parser$add_argument("--LEARN_READS_N", type='integer', default=1e06, help="Default 1e6. The minimum number of reads to use for error rate learning. Samples are read into memory until at least this number of reads has been reached, or all provided samples have been read in.")
parser$add_argument("--OMEGA_A", type='double', default=1e-40, help="This parameter sets the threshold for when DADA2 calls unique sequences significantly overabundant, and therefore creates a new cluster with that sequence as the center. The default value is 1e-40, which is a conservative setting to avoid making false positive inferences, but which comes at the cost of reducing the ability to identify some rare variants.")
parser$add_argument("--USE_QUALS", type='logical', default='TRUE', help="If TRUE, the dada(...) error model takes into account the consensus quality score of the dereplicated unique sequences. If FALSE, quality scores are ignored. The default is TRUE, however if applying DADA2 to pyrosequenced data it is recommended to set USE_QUALS to FALSE, as quality scores are not informative about substitution error rates in pyrosequencing.")
parser$add_argument("--USE_KMERS", type='logical', default='TRUE', help="If TRUE, a 5-mer distance screen is performed prior to performing each pairwise alignment, and if the 5mer-distance is greater than KDIST_CUTOFF, no alignment is performed. TRUE by default.")
parser$add_argument('--KDIST_CUTOFF', type='double', default=0.42, help='The default value of 0.42 was chosen to screen pairs of sequences that differ by greater 10 percent, and was calibrated on Illumina sequenced 16S amplicon data. The assumption is that sequences that differ by such a large amount cannot be linked by amplicon errors (i.e. if you sequence one, you won`t get a read of other) and so careful (and costly) alignment is unnecessary.')
parser$add_argument("--BAND_SIZE", type='integer', default=16, help = "When set, banded Needleman-Wunsch alignments are performed. Banding restricts the net cumulative number of insertion of one sequence relative to the other. The default value of BAND_SIZE is 16. If DADA is applied to marker genes with high rates of indels, such as the ITS region in fungi, the BAND_SIZE parameter should be increased. Setting BAND_SIZE to a negative number turns off banding (i.e. full Needleman-Wunsch).")
parser$add_argument("--SCORE_MATRIX", type='character', help="The path to the score matrix for the Needleman-Wunsch alignment. This is a 4x4 matrix as no ambiguous nucleotides are allowed. Default is nuc44: -4 for mismatches, +5 for matches.")
parser$add_argument("--GAP_PENALTY", type='integer', default=-8, help='The cost of gaps in the Needleman-Wunsch alignment. Default is -8.')
parser$add_argument("--HOMOPOLYMER_GAP_PENALTY", default=NULL, help='The cost of gaps in homopolymer regions (>=3 repeated bases). Default is NULL, which causes homopolymer gaps to be treated as normal gaps.')
parser$add_argument("--MIN_FOLD", type='integer', default=1, help= 'The minimum fold-overabundance for sequences to form new clusters. Default value is 1, which means this criteria is ignored.')
parser$add_argument("--MIN_HAMMING", type ='integer', default=1, help= 'The minimum hamming-separation for sequences to form new clusters. Default value is 1. which means this criteria is ignored.')
parser$add_argument("--MAX_CLUST", type='integer', default=0, help='The maximum number of clusters. Once this many clusters have been created, the algorithm terminates regardless of whether the statistical model suggests more sample sequences exist. If set to 0 this argument is ignored. Default value is 0.')
parser$add_argument("--MAX_CONSIST", type='integer', default=10, help= 'The maximum number of steps when selfConsist=TRUE. If convergence is not reached in MAX_CONSIST steps, the algorithm will terminate with a warning message. Default value is 10.')
parser$add_argument("--VERBOSE", type= 'logical', default=TRUE, help='If TRUE progress messages from the algorithm are printed. Warning: There is a lot of output. Default is FALSE.')


# Parse arguments
args <- parser$parse_args()
##### VALIDATE #####
#Ensure the chimrm selection is valid
if (!any(args$chimrm == c('pooled', 'consensus', 'per-sample', 'FALSE'))) stop(sprintf('%s is not a vaild chimrm option', args$chimrm))

 

##### FILES #####
fnFs = sort(args$fwd) # sort for pairing files 
if(!all(file.exists(fnFs))) stop(sprintf("The following file does not exist: %s\n",fnFs[!file.exists(fnFs)])) # Check to see if files exist
if(any(duplicated(fnFs))) stop(sprintf("The following files are duplicated: %s\n",fnFs[!duplicated(fnFs)])) 
# If we have paired end data
if (all(args$rev != F)){
	fnRs = sort(args$rev) # sort it
	if(!all(file.exists(fnRs))) stop(sprintf("The following file does not exist: %s\n",fnRs[!file.exists(fnRs)])) # validate
	if(any(duplicated(fnRs))) stop(sprintf("The following files are duplicated: %s\n",fnFs[!duplicated(fnRs)])) 
	if(any(fnFs %in% fnRs || fnRs %in% fnFs)) stop(sprintf("The following files are in both read streams %s, %s,\n", fnFs[fnFs %in% fnRs], fnRs[ fnRs %in% fnFs])) 
	print('File pairs to be processed in DADA2:')
	print(data.frame(fwd_files=fnFs, rev_files=fnRs)) # provide debugging
}else{
	print('File pairs to be processed in DADA2:')
	print(data.frame(fwd_files=fnFs))
}

##### SAMPLE_NAMES #####
study.name = args$prefix

if (any(args$samp_fields != F) && any(args$samp_list == F) && any(args$samp_regex == F)){ # EXTRACT_NAMES_FROM_READS
	# Extract sample names, assuming filenames have format: {ds}_{resource]_{sample}_{factor}_R1-trimmed.fastq
	sample.names = lapply(strsplit(basename(fnFs), args$fields_delim), function(x){paste(x[as.integer(args$samp_fields)],collapse = args$fields_delim)}) # grab the delimited feilds
} else if (any(args$samp_fields == F) && any(args$samp_list == F) && any(args$samp_regex == F)){ # KEEP_FULL_NAMES
	sample.names = basename(fnFs)
} else if (any(args$samp_fields == F) && any(args$samp_list == F) && any(args$samp_regex != F)){ # EXTRACT_NAMES_WITH_REGEX
	library(stringr)
	sample.names = apply(format(str_match(basename(fnFs), args$samp_regex)[,-1]), 1, paste, collapse="_") #extract using regex
} else if (any(args$samp_list != F)){ # GET_NAMES_FROM_LIST
	sample.names = args$samp_list[order(args$fwd)]
	if (any(args$samp_fields != F) && any(args$samp_regex == F)){ # SPLIT_LIST_NAMES_WITH_DELIM
		sample.names = lapply(strsplit(args$samp_list, args$fields_delim), function(x){paste(x[as.integer(args$samp_fields)],collapse = args$fields_delim)})
	} else if (any(args$samp_fields == F) && any(args$samp_regex != F)){ # LIST_WITH_REGEX
		library(stringr)
		sample.names = apply(format(str_match(basename(fnFs), args$samp_regex)[,-1]), 1, paste, collapse="_") #extract using regex
	} else if (any(args$samp_fields != F) && any(args$samp_regex != F)){
		stop('Can not use both regex and delimited at the same time for listed sample names')
	} else if(any(args$samp_fields == F) && any(args$samp_list != F) && any(args$samp_regex == F)){
		warning('Using raw names from sample list')
	} else {
		stop(sprintf('Invalid combination of --samp_feilds:%s, --samp_regex:%s, --samp_list:%s', args$samp_fields, args$samp_regex, args$samp_list))
	}
}

if (!all(!duplicated(sample.names))){
	stop(sprintf('The following Sample names are not unique: %s\n',sample.names[!duplicated(sample.names)]))
}

##### SETUP_DADA #####
setDadaOpt(OMEGA_A=args$OMEGA_A, USE_QUALS=args$USE_QUALS, USE_KMERS=args$USE_KMERS, KDIST_CUTOFF=args$KDIST_CUTOFF,
		   BAND_SIZE=as.numeric(args$BAND_SIZE, GAP_PENALTY=args$GAP_PENALTY, HOMOPOLYMER_GAP_PENALTY=args$HOMOPOLYMER_GAP_PENALTY,
		   MIN_FOLD=as.numeric(args$MIN_FOLD), MIN_HAMMING=as.numeric(args$MIN_HAMMING), MAX_CLUST=as.numeric(args$MAMAX_CLUST), MAX_CONSIST=as.numeric(args$MAX_CONSIST), VERBOSE=args$VERBOSE))
if (!is.null(args$SCORE_MATRIX)){
	setDadaOpt(SCORE_MATRIX=read.table(args$SCORE_MATRIX))
}

##### DADA #####

# Learn error rates
write('Learning error rates R1', stderr())
errF = learnErrors(fnFs, nreads = args$LEARN_READS_N, multithread=TRUE, randomize=TRUE)
#plotErrors(errF, nominalQ=TRUE)

# Dereplication
write('Dereplicating R1', stderr())
derepFs = derepFastq(fnFs, verbose=TRUE)
names(derepFs) = sample.names # Name the derep-class objects by the sample names

# sample Infrence
write('DADA2 Core R1',stderr())
dadaFs = dada(derepFs, err=errF, multithread=TRUE)

if (any(args$rev != F)){
	# Learn error rates
	write('Learning error rates R2', stderr())
	errR = learnErrors(fnRs, nreads = args$LEARN_READS_N, multithread=TRUE, randomize=TRUE)
	#plotErrors(errF, nominalQ=TRUE)
	
	# Dereplication
	write('Dereplicating', stderr())
	derepRs = derepFastq(fnRs, verbose=TRUE)
	names(derepRs) = sample.names # Name the derep-class objects by the sample names
	
	# sample Infrence
	write('DADA2 Core',stderr())
	dadaRs = dada(derepRs, err=errR, multithread=TRUE)
	
	# merge pair reads
	write('Merge pair end reads', stderr())
	dada_reads = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
} else {
	dada_reads = dadaFs
}

# construct ASV table
write('Creating ASV table', stderr())
seqtab = makeSequenceTable(dada_reads)
write('Sequence lengths of ASVs:', stdout())
write(table(nchar(getSequences(seqtab))), stdout())

if (args$chimrm != 'FALSE'){
	# Remove Chimeras
	write('Remove Chimeras', stderr())
	seqtab.chimrm = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
	write.table(seqtab.chimrm, file=args$asv_out, sep="\t")
} else {
	write.table(seqtab, file=args$asv_out, sep="\t")
}

# # Track Reads
# getN = function(x) sum(getUniques(x))
# track = data.frame(samp=unlist(sample.names), denoisedF = sapply(dadaFs, getN))
# if (args$rev != FALSE){
# 	track['denoisedR'] = sapply(dadaRs, getN)
# 	track['merged'] = sapply(mergers, getN)
# }
# if (args$chimrm != F) track['nonchim'] = rowSums(seqtab.chimrm)
# # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
# write('Reads Surviving each step', stderr())
# write.table(track, stderr())
