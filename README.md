**+===================================+**

**| DADA2 for Command Line and Galaxy |**

**+===================================+**

# Intro
DADA2 is useful for metagenomic studies where taxonomic information of a cell can be inferred from a variable region of the genome(i.e. 16s, ITS), that is targeted by PCR amplification and sequencing. Unlike many other programs that cluster sequence reads into clumps to deal with errors in sequencing; DADA2 infers the error patterns of the sequenced reads, then uses the information to gather reads into their most likely true sequence.
In addition a Snakemake pipeline for this workflow can be found at https://github.com/cErikson/GeneLab_DADA2_snakemake_Pipeline.
DADA2: High-resolution sample inference from Illumina amplicon data: https://doi.org/10.1038/nmeth.3869 (not just for illumina)

# Requirements
R >= 3.4: DADA2, biom, argparse, stringr 
Python: Biom, argparse

# Workflow
## Raw Reads
The reads should originate from a targeted locus present in the metagenome. The reads can be from illumina single or paired end sequencing, or with some setting changes from the 454 or ion torrent sequencing platforms.

## Quality Control
### Trimming
DADA2 requires that reads have synthetic sequences(primers, adapters) removed, and be free of ambiguous base calls. In addition removal of highly erroneous reads will improve performance.
Removal of synthetic sequences can be achieved with the following cutadapt command.

`cutadapt -g {fwd_prime} -G {rev_prime}: -a {rev_revcomp_prime} -A {fwd_revcomp_prime} -o {output.fwd_trim} -p {output.rev_trim} --max-n=0 -m {min_length} --error-rate {max_error_rate_for_trimming_matched} --trim-n --pair-filter any {input.fwd} {input.rev} > {log}`

In addition,
For 454: `--maximum-length={select_based_on_chemisty}`
For Ion Torrent: `--trimLeft=15`
`--trimmed-only` can be used where primers are expected to be in the sequenced reads, which will help remove contamination. This option will not work for primers in the style of Earth Micro Biome primers, where sequencing primers target the primer region used to target the biological sequence during amplification.

### Filtering
The removal of highly erroneous reads will speed up DADA2 and improve its quality. It is recommend this is accomplished by filtering out reads, where the filter criteria is the maximum allowable probability of errors in a read, known as maximum expected error rate. This can be accomplished by the `dada2_filter.R` script, with the following command. 

`filtered.R -f {input.fwd} -r {input.rev} -p {output.directory} -s 1 -S '.' -q {Trim_bases_with_Q>} -n {min_length} -d {filter_reads with_Q>} -E {Max_error_rate} --rm_phix {True/False}`

For paired end data differing options may be set for forward and reverse read, for example `-E 0.7,1.5` will allow for the probability of 0.7 errors in the forward read, and 1.5 errors in the reverse.

## DADA2 Core
dada2.R is the will take the reads and calculate a Amplicon Sequence Variant table. Which is a count table of dereplicated sequences by sample.
Sample names used in the table my be created from the file names or provided by 1. selecting a delimited field with `-s 1,2,3 -S '_', 2. A regex `-R '(\\w*)\\.fastq` (double escapes), 3. A list of names  `--samp_list Samp_A,Samp_B`, possibly in combination with option 1 or 2. Be sure to include a unique identifier. 

`dada2.R -f {input.fwd} -r {input.rev} -o {output.asv} {samp_name}`

In addition,
for 454 and Ion Torrent: `--HOMOPOLYMER_GAP_PENALTY '-1' --BAND_SIZE 32`
if paired-end data does not overlap: `--no_overlap True`

## Assign taxonomic information
Now to assign taxonomic information to the ASVs, use the `dada2_assign_taxa.R` script. First taxonomic training files need to be downloaded. https://benjjneb.github.io/dada2/training.html. Note that species assignment may fail if `--no_overlap True` was set in dada2.R.

`dada2_assign_taxa.R {input.asv_table} {output.taxa_table} -t {taxa_train_fasta} -s {species_train_fasta]} --tryRC {try_rev_complement} --taxLevels Kingdom Phylum Class Order Family Genus Species --allowMultiple {return_multiple_species_hit} --multithread {threads}`

## Conversion to BIOM
Finally, it is possible to convert the output to the popular BIOM format. It is possible to add sample metadata to the BIOM. The metadata needs a header and a column which has an id unique to each sample`-c`. The regex `-r` is used to extract this unique id from the sample names in the ASV. 

`dada2biom.py {input.asv} {output.biom} -t {input.taxa} -s {input.metadata} -c {meta_samp_name_feild} -r {asv_samp_name_regex]} -n {Study_name}`

