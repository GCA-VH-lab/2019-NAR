#!/usr/bin/env python
#
#
# 
# 2018.12.27 - updated to handle gene lists
# 2019.05.21 - updated to work with py 3.6 new libraries
#

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2019"
__version__		= "0.1.6"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import sys
import argparse
import pandas as pd
import pysam
from collections import defaultdict

parser = argparse.ArgumentParser(description='Computes metagene tables around stop codons TAA, TAG and TGA.')
parser.add_argument('-i', type=str, help='input table of gene coverage in *.hd5 format')
parser.add_argument('-prefix', type=str, help='output file name first part; {prefix}_TAA.csv, {prefix}_TAG.csv, {prefix}_TGA.csv')
parser.add_argument('-annot',  type=str, help='GTF annotation file', default='0-References/genome.gtf.gz')
parser.add_argument('-genome',  type=str, help='Genome sequnece in FastA', default='0-References/Genome.fa')
parser.add_argument('-norm',  type=str, help='Normalisation: if YES genes contribute equally to final table', default='Yes')
parser.add_argument('-subsets',  type=str, help='Split to subsets based Stop codon', default='Yes')
parser.add_argument('-th',  type=float, help='Stop region total coverage - exclude if below', default=5.0)
parser.add_argument('-span',  type=int, help='Positions included to table around stop', default=120)
parser.add_argument('-col',  type=str, help='column for values - default is "sum"', default='sum')
parser.add_argument('-glist',  type=str, help='filename with gene names in interest one per line', default=None)
args = parser.parse_args()

    
print("\n\
-i      input *.h5:      {}\n\
-prefix part of oufile : {}\n\
-annot   annotation GTF: {}\n\
-genome    genome FastA: {}\n\
-col            columns: {}\n\
-th    region threshold: {}\n\
-span   included region: {}\n\
-norm  to gene coverage: {}\n\
-subsets  by stop codon: {}\n\
-glist  subsets by list: {}\n".format(args.i, args.prefix,args.annot, args.genome, args.col, args.th, args.span, args.norm, args.subsets, args.glist))

usage = "./hdf_2_metagene_tables_1_dev.py  -i WTS1_5-End_21-33_idx_assign_rpm.h5   -prefix WTS1_stop_metagene   -norm  Yes\n\n \
\t\t  assumes *.hd5 file keys structure \"Chr/For\"   - chromosome and strand \n " 


if (args.i==None)|(args.prefix==None):
     sys.exit("\n  usage:\n\t{}\n".format(usage))

infile_h5 = args.i #ok
outfile_str = args.prefix
col = args.col
Span = args.span #ok
annotation = args.annot
genomefile = args.genome
normalise = args.norm.upper()
Threshold = args.th #ok
split_by_codons = args.subsets.upper()

########## Functions        
def df_framing2(df, index, columns, strand="+"):
    # create df2
    df1 = df[columns].copy()
    df2 = pd.DataFrame(0, index=index, columns=columns)
    
    df1 = df1.add(df2, fill_value=0, axis=1)  
    if strand == "+":
        df1.reset_index(inplace=True)  # reset index
        return df1[columns]
    
    elif strand == "-":
        df1 = df1[::-1]  # reverts table
        df1.reset_index(inplace=True)  # reset index
        return df1[columns]
    else:
        # error
        print("ERROR! Expext '+'/'-' but found {} for strand".format(strand))

def get_part_from_gtf(annotation, reference=None, feature="CDS"):   
    tabixfile = pysam.TabixFile(annotation, parser=pysam.asGTF())
    return [gtf for gtf in tabixfile.fetch(reference=reference) if (gtf.feature == feature)]

def yeastChr():
    # Ordered yeast Chr list short names from ensembl
    return ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV',
            'XVI','Mito']

def complement(seq):
    basecomplement = {'A':'T','C':'G','G':'C','T':'A','Y':'R','R':'Y','M':'K','K':'M','W':'W','V':'B','B':'V','H':'D','D':'H','N':'N',
                      'a':'t','c':'g','g':'c','t':'a','y':'r','r':'y','m':'k','k':'m','w':'w','v':'b','b':'v','h':'d','d':'h','n':'n' }
    return ''.join([basecomplement[base] for base in seq])

def revcompl(seq):
    return complement(seq[::-1])

def read_FASTA(filename, splitstr='|', SplitHeader=True):
    with open(filename) as file:

        if SplitHeader:
            return [(part[0].split(splitstr),
                part[2].replace('\n', ''))
                for part in 
                [entry.partition('\n')
                    for entry in file.read().split('>')[1:]]]
        else: 
            return [(part[0],
                part[2].replace('\n', ''))
                for part in 
                [entry.partition('\n')
                    for entry in file.read().split('>')[1:]]]

def read_FASTA_dictionary(filename, splitstr='|', SplitHeader=False):
    return {info[0]: seq for info, seq in read_FASTA(filename)}
#############################

genome      = read_FASTA_dictionary(genomefile, SplitHeader=False)
stop_codons = ['TAA', 'TAG', 'TGA']

max_rpms  = 1000 # remove outliers have effect when normalisation is off
# counting & debug
c_taa = c_tag = c_tga = 0
c = c1 = c2  = 0

notStopCodon  = defaultdict() #debugging
specifiCodon  = defaultdict()
genes2highrpm = defaultdict()

# input h5
hd5 = pd.HDFStore(infile_h5, "r")

# filename
if normalise == 'YES':
    outfile_str = "{}-norm".format(outfile_str)    
outfile_str = "{}_Sp{}_th{}".format(outfile_str, Span, int(args.th))
# outfiles
outf_stop_TAA = outfile_str + "_TAA.csv"
outf_stop_TAG = outfile_str + "_TAG.csv"
outf_stop_TGA = outfile_str + "_TGA.csv"
outf_stop_sum = outfile_str + "_sum.csv"
# metagene summary df
columns          = hd5[hd5.keys()[0]].columns
columns = columns[2:]
meta_stop_df_TAA = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
meta_stop_df_TAG = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
meta_stop_df_TGA = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)
#meta_stop_df_sum = pd.DataFrame(index=range(0, 2 * Span + 1), columns=columns).fillna(0)

# genes listfile
gnames = []
if args.glist==None:
    pass
else:             # read gene names from the file    
    try:
        # open file and read the content in a list
        with open(args.glist, 'r') as filehandle:  
            gnames = [gn.rstrip() for gn in filehandle.readlines()]
    except:
        print('Could not open genes listfile "{}"'.format(args.glist))
        sys.exit(1) 

### MAIN body START ###
for ref in yeastChr():

    df_f = hd5['For_rpm/'+ ref]
    df_r = hd5['Rev_rpm/'+ ref]

    stop_gtf_l = get_part_from_gtf(annotation, reference=ref, feature="stop_codon") # gtf part for stop

    for gtf in stop_gtf_l:
        
        # using genes list to restrict genes for metagene
        if args.glist==None:              # don't do nothing
            pass
        else:
            if gtf.gene_id not in gnames: # eituse eitus viib tagasi üles kui 
                continue
        
            
        stop_codon = genome[ref][gtf.start:gtf.start+3]                      # get stop codon
        stop_codon = stop_codon if gtf.strand=='+' else revcompl(stop_codon) # revcomp rev_strand

        #test coverage
        coverage_for = df_f.loc[gtf.start - Span:gtf.start + 3, "sum"].sum()
        coverage_rev = df_r.loc[gtf.start:gtf.start + Span - 1, "sum"].sum()
        coverage_rpm = coverage_for if gtf.strand == '+' else coverage_rev

        if coverage_rpm < Threshold: #FILTER 3
            c1+=1
            continue
            
        # get regions
        df_1    = df_f.loc[gtf.start - Span:gtf.start + Span, :]     #  Forvard
        index_1 = range(gtf.start - Span, gtf.start + Span + 1)
        df_2    = df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1, :] # Reverse
        index_2 = range(gtf.end - Span - 1, gtf.end + Span)
        # chose proper index & dataframe
        df, index = (df_1,index_1) if gtf.strand == '+' else (df_2, index_2)  # select one based strand
        #
        # normalisation part can added here
        if normalise == 'YES':
            normalisation_factor = df[col].sum()
            df = df[columns]/float(normalisation_factor)

        df = df_framing2(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
        
        # sum up
        if   stop_codon == 'TAA':
            meta_stop_df_TAA = meta_stop_df_TAA + df
            c_taa+=1
        elif stop_codon == 'TAG':
            meta_stop_df_TAG = meta_stop_df_TAG + df
            specifiCodon[gtf.gene_id]=coverage_rpm
            c_tag+=1
        elif stop_codon == 'TGA':
            meta_stop_df_TGA = meta_stop_df_TGA + df
            c_tga+=1
        else:
            pass
### MAIN body END  ###


#save files
# if subsets == YES
if split_by_codons == 'YES':
    # create 'rel_Pos' index 
    index = list(range(-Span, Span + 1))
    meta_stop_df_TAA['rel_Pos'] = index 
    meta_stop_df_TAA.set_index('rel_Pos', inplace=True)

    meta_stop_df_TAG['rel_Pos'] = index 
    meta_stop_df_TAG.set_index('rel_Pos', inplace=True)

    meta_stop_df_TGA['rel_Pos'] = index 
    meta_stop_df_TGA.set_index('rel_Pos', inplace=True)
    
    
    meta_stop_df_TAG.to_csv(outf_stop_TAG, sep='\t', header=True, index=True)
    meta_stop_df_TAA.to_csv(outf_stop_TAA, sep='\t', header=True, index=True)
    meta_stop_df_TGA.to_csv(outf_stop_TGA, sep='\t', header=True, index=True)
    print("TAA: {:,}\nTAG: {:,}\nTGA: {:,}\n".format(c_taa, c_tag, c_tga))
    print("{}\n{}\n{}\n".format(outf_stop_TAA, outf_stop_TAG, outf_stop_TGA))
    

meta_stop_df_sum = meta_stop_df_TAA  + meta_stop_df_TAG + meta_stop_df_TGA 
meta_stop_df_sum.to_csv(outf_stop_sum, sep='\t', header=True, index=True)

print("Bleow threshold {}: {:,}\nIncluded stops: {}".format(Threshold, c1, c_taa+c_tag+c_tga ))
print("{}\n//".format(outf_stop_sum))

hd5.close()


