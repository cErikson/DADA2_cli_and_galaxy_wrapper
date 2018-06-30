#!/usr/bin/env Rscript
################################################################################
# Title: dada2_filter.R
# Discription: A wrapper for DADA2 filter
# Author: Christian Erikson
# GitHub: https://github.com/cErikson/DADA2_cli_and_galaxy_wrapper
# Date: 6/26/18
################################################################################
packages = c("dada2", "argparse")
if (any(!(packages %in% installed.packages()[,"Package"]))) stop(sprintf('These packages are required: %s', packages[!(packages %in% installed.packages()[,"Package"])]))
library(dada2)
library(argparse)

##### PARSER #####
# Initialize argument parser
parser <- ArgumentParser(description='A Command line wrapper for DADA2.', epilog='Note: When filtering paired reads... If a length 1 vector is provided, the same parameter value is used for the forward and reverse reads. If a length 2 vector is provided, the first value is used for the forward reads, and the second for the reverse reads')

parser$add_argument("-p","--prefix",type="character",default='dada2-filter_' ,help="Filtered file name prefix.")
parser$add_argument("-f","--fwd", nargs='+', type="character",help="Name/path to fastq or fastq.gz files for single end sequencing or forward reads if paired end sequencing")
parser$add_argument("-r","--rev", nargs='+', type="character",default=FALSE,help="Name/path to fastq or fastq.gz files for reverse reads if paired end sequencing")
parser$add_argument("-s","--samp_fields", nargs='+', type="integer", default=FALSE, help="The fields in the file name that should be used for sample names.")
parser$add_argument("-S","--fields_delim", nargs='+', type="character", default='_', help="The field delimiter used in the file name")
parser$add_argument("-X","--samp_regex", action='store',type="character",default=FALSE,help="Regex used to create filenames. Double escapes \\w")

parser$add_argument("-q", "--truncQ", nargs='+',default=2, help='(Optional). Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.')
parser$add_argument("-l", "--truncLen", nargs='+',default = 0, help="(Optional). Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded.")
parser$add_argument("-L", "--trimLeft", nargs='+',default = 0, help="(Optional). Default 0. The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft.")
parser$add_argument("-x", "--maxLen", nargs='+',default= 'Inf', help= "(Optional). Default Inf (no maximum). Remove reads with length greater than maxLen. maxLen is enforced before trimming and truncation.")
parser$add_argument("-n", "--minLen", nargs='+',default=20, help="(Optional). Default 20. Remove reads with length less than minLen. minLen is enforced after trimming and truncation.")
parser$add_argument("-N", "--maxN", nargs='+',default=0, help="(Optional). Default 0. After truncation, sequences with more than maxN Ns will be discarded. Note that dada does not allow Ns.")
parser$add_argument("-d", "--minQ", nargs='+',default=0, help="(Optional). Default 0. After truncation, reads contain a quality score less than minQ will be discarded.")
parser$add_argument("-E", "--maxEE", nargs='+', default='Inf',  help="(Optional). Default Inf (no EE filtering). Two are needed for pairedend. vaules After truncation, reads with higher than maxEE 'expected errors' will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))")
parser$add_argument("--rm.phix",  default=T, help="(Optional). Default TRUE. If TRUE, discard reads that match against the phiX genome, as determined by isPhiX.")
parser$add_argument("--primer.fwd", default=NULL, help="(Optional). Default NULL. Paired-read filtering only. A character string defining the forward primer. Only allows unambiguous nucleotides. The primer will be compared to the first len(primer.fwd) nucleotides at the start of the read. If there is not an exact match, the read is filtered out. For paired reads, the reverse read is also interrogated, and if the primer is detected on the reverse read, the forward/reverse reads are swapped.")
parser$add_argument("--matchIDs",default=F, help="(Optional). Default FALSE. Paired-read filtering only. Whether to enforce matching between the id-line sequence identifiers of the forward and reverse fastq files. If TRUE, only paired reads that share id fields (see below) are output. If FALSE, no read ID checking is done. Note: matchIDs=FALSE essentially assumes matching order between forward and reverse reads. If that matched order is not present future processing steps may break (in particular mergePairs).")
parser$add_argument("--id.sep", default='\\s', help="(Optional). Default '\\\\s' (white-space). Paired-read filtering only. The separator between fields in the id-line of the input fastq files. Passed to the strsplit.")
parser$add_argument("--id.field", default='NULL', help="(Optional). Default NULL (automatic detection). Paired-read filtering only. The field of the id-line containing the sequence identifier. If NULL (the default) and matchIDs is TRUE, the function attempts to automatically detect the sequence identifier field under the assumption of Illumina formatted output.")
parser$add_argument("--multithread", default=F, help="(Optional). Default is FALSE. If TRUE, input files are filtered in parallel via mclapply. If an integer is provided, it is passed to the mc.cores argument of mclapply. Note that the parallelization here is by forking, and each process is loading another fastq file into memory. If memory is an issue, execute in a clean environment and reduce the chunk size n and/or the number of threads.")
parser$add_argument("--n_reads", default= 1e5, help="(Optional). Default 1e5. The number of records (reads) to read in and filter at any one time. This controls the peak memory requirement so that very large fastq files are supported. See FastqStreamer for details.")
parser$add_argument("-v", "--verbose", default= TRUE, help="((Optional). Default TRUE. Whether to output status messages.")


args <- parser$parse_args()


#Ensure there is a sample name type
if (sum(all(args$samp_fields != F), all(args$samp_regex != F)) != 1) stop('Provide a one option to collect sample names: -S -R') 

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
	print('File pairs to be filtered in DADA2:')
	print(data.frame(fwd_files=fnFs, rev_files=fnRs)) # provide debugging
}else{
	print('File pairs to be filtered in DADA2:')
	print(data.frame(fwd_files=fnFs))
}

##### SAMPLE_NAMES #####
study.name = args$prefix

if (any(args$samp_fields != F)){
	# Extract sample names, assuming filenames have format: {ds}_{resource]_{sample}_{factor}_R1-trimmed.fastq
	sample.names = lapply(strsplit(basename(fnFs), args$fields_delim), function(x){paste(x[as.integer(args$samp_fields)],collapse = args$fields_delim)})
} else if (args$samp_regex != F){
	library(stringr)
	sample.names = apply(format(str_match(basename(fnFs), args$samp_regex)[,-1]), 1, paste, collapse="_")
}
if (!all(!duplicated(sample.names))){
	stop(sprintf('The following Sample names are not unique: %s\n',sample.names[!duplicated(sample.names)]))
}


##### FILTER #####
# Place filtered files in filtered/ subdirectory
if (any(args$rev != F)){
	write('Filtering reads',stderr())
	filtFs = file.path( paste0(study.name,sample.names, "_R1-filtered.fastq.gz"))
	filtRs = file.path( paste0(study.name,sample.names, "_R2-filtered.fastq.gz")) # create file names
	out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,
						truncQ = as.numeric(args$truncQ), truncLen = as.numeric(args$truncLen) , trimLeft = as.numeric(args$trimLeft), maxLen = as.numeric(args$maxLen), minLen = as.numeric(args$minLen),
						maxN = as.numeric(args$maxN), minQ = as.numeric(args$minQ), maxEE = as.numeric(args$maxEE) , rm.phix = args$rm.phix, primer.fwd = args$primer.fwd,
						matchIDs = args$matchIDs, id.sep = args$id.sep, id.field = args$id.field,
						multithread = args$multithread, n = args$n_reads, verbose = args$verbose)
	#
	exists = file.exists(filtFs) & file.exists(filtRs)
	filtFs = filtFs[exists]
	filtRs = filtRs[exists]
	if (any(exists==FALSE)){
		write('The fliter has removed all the reads from the following files, thus they will no longer be analyized', stderr())
		write(filtFs[exists==FALSE], stderr())
		write(filtRs[exists==FALSE], stderr())
	}
} else {
		write('Filtering reads',stderr())
		filtFs = file.path( 'data/filtered', paste0(study.name,sample.names, "_R1-filtered.fastq.gz"))
		out = filterAndTrim(fnFs, filtFs, 
							truncQ = as.numeric(args$truncQ), truncLen = as.numeric(args$truncLen) , trimLeft = as.numeric(args$trimLeft), maxLen = as.numeric(args$maxLen), minLen = as.numeric(args$minLen),
							maxN = as.numeric(args$maxN), minQ = as.numeric(args$minQ), maxEE = as.numeric(args$maxEE) , rm.phix = args$rm.phix, primer.fwd = args$primer.fwd,
							matchIDs = args$matchIDs, id.sep = args$id.sep, id.field = args$id.field,
							multithread = args$multithread, n = args$n_reads, verbose = args$verbose)
		#
		exists = file.exists(filtFs)
		filtFs = filtFs[exists]
		if (any(exists==FALSE)){
			write('The fliter has removed all the reads from the following files, thus they will no longer be analyized', stderr())
			write(filtFs[exists==FALSE], stderr())
		}
}

