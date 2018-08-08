#!/usr/bin/env python
"""
Project:
@author: Christian Erikson

TODO: 
select if taxa is qiime
taxa head space delim?    
"""
##### Debug #####
testing=False

###### Imports #####
import argparse as arg
import pdb
import sys
import numpy as np
import biom
import re


##### Args #####
if __name__ == "__main__" and testing != True:
    parser = arg.ArgumentParser()
    # Positional mandatory arguments
    parser.add_argument("asv", help="The ASV count table generated by DADA", type=str)
    parser.add_argument("out", help="Path to write the BIOM", type=str)
    # Optional arguments
    parser.add_argument("-t", "--taxa", help="The taxa table generated by DADA", type=str, default=False)
    parser.add_argument("-s", "--meta", help="Sample metadata, in which rows are samples", type=str, default=False)
    parser.add_argument("-d", "--delim", help="The delimiter used in sample metadata.", type=str, default='\t')
    parser.add_argument("-c", "--samp_col", help="The column name, which has sample ids. Ids need to have a coraspnding entry in the ASV table, which may be a substring int the ASV.", type=str, default=False)
    parser.add_argument("-r", "--regex", help="Regex used to extract sample id names from the ASV ids", type=str, default='(.*)')
    parser.add_argument("-n", "--study_id", help="Study ID to be used in the Biom ", type=str, default='DADA2 study')
    # Parse arguments
    args = parser.parse_args()
    
##### Testing #####
else:   # Else test
    sys.stderr.writelines("!!!!!___RUNNING_IN_TESTING_MODE_WITH_TEST_ARGS___!!!!!\n")
    class test_args(object):
        asv='/home/christian/lab/genelab/analyses/GLDS-191/data/GLDS-191_metagenomics_dada2-asv-matrix.tsv'
        taxa='/home/christian/lab/genelab/analyses/GLDS-191/data/GLDS-191_metagenomics_dada2-taxonomy-dataframe.tsv'
        meta='/home/christian/lab/genelab/analyses/GLDS-191/data/metadata/isa_files/s_MGP79314.txt'
        delim='\t'
        samp_col='Sample Name'
        regex= False
        study_id='GLDS-146'
    args = test_args()
    
if args.regex is not False:
    rec=re.compile(args.regex)

# Read ASV counts
with open(args.asv, 'r') as fh:
    taxa_id=[x.strip('\'\"\n') for x in fh.readline().strip().split(args.delim)]
    counts=[[y.strip('\'\"\n') for y in x.strip().split(args.delim)] for x in fh]
    samp_id=[x.pop(0).strip('\'\"\n') for x in counts]
counts=np.array(counts)
if args.regex is not False:
    ext=[]
    for ident in samp_id:
        try:
            ext.append(rec.match(ident).groups()[0])
        except AttributeError:
            sys.exit('The regex `{}` failed to match `{}`'.format(rec.pattern, ident))
    samp_id=ext

# Read Taxa
taxa={}
tax_lev={'Kingdom':'k__', 'Phylum':'p__', 'Class':'c__', 'Order':'o__', 'Family':'f__', 'Genus':'g__', 'Species':'s__'}
if args.taxa is not False:
    with open(args.taxa, 'r') as fh:
        taxa_head=[x.strip('\'\"\n') for x in fh.readline().strip().split()]  ### FIXXXXXXXXXXX
        for line in fh:
            cell=[x.strip('\n\'\"') for x in line.split(args.delim)]
            taxa[cell[0]]=[tax_lev[y]+x for x,y in zip(cell[1:],taxa_head)]
    try:
        obs_meta=[{'taxonomy':taxa.pop(x)} for x in taxa_id]
    except KeyError as e:
        pdb.set_trace()
        sys.exit('The taxonomic metadata was not found for: '+str(e)+'\n A key in the metadat is '+str(list(taxa.keys())[0])  )
    if len(taxa) is not 0:
        sys.stderr.write('After binding taxa metadata, the following entries went unused:\n'+ str(taxa.keys()) )
else:
    obs_meta=None

# Read sample metadata
samp={}
if args.meta is not False:
    with open(args.meta, 'r') as fh:
        samp_head=[x.strip('\'\"\n') for x in fh.readline().strip().split(args.delim)]
        try:
            sname=samp_head.index(args.samp_col)
        except ValueError as e:
            sys.exit('The Sample ID column was not found in the metadata.\n'+(str(e)))
        samp_head.pop(sname)
        for line in fh:
            cell=line.split(args.delim)
            name=cell.pop(sname).strip('\n\'\"')
            samp[name]={x.strip('\n\'\"'):y.strip('\n\'\"') for x,y in zip(samp_head, cell)}
    try:
        samp_meta=[samp.pop(x) for x in samp_id]
    except KeyError as e:
        sys.exit('The sample metadata was not found for key: '+str(e)+'\nThe keys are: '+str(samp.keys()) )
    if len(taxa) is not 0:
        sys.stderr.write('After binding taxa metadata, the following entries went unused:\n'+ str(samp.keys()) )
else:
    samp_meta=None
    
biom_table = biom.Table(counts.transpose(), taxa_id, samp_id, obs_meta, samp_meta, table_id='Example Table', generated_by='DADA2BIOM.py')

with open(args.out, 'w+') as oh:
    biom_table.to_json('DADA2BIOM', oh)
    
