#!/usr/bin/env python
#

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2019"
__version__		= "0.0.2"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import sys
import argparse
import pandas as pd
import pysam


parser = argparse.ArgumentParser(description="Computes relative 3' UTR coverage")
parser.add_argument('-i', type=str, help='input table of gene coverage in *.hd5 format')
parser.add_argument('-o', type=str, help='output file name  output.txt', default="3utr_coverage.txt")
parser.add_argument('-annot',  type=str, help='GTF annotation file', default='0-References/genome.gtf.gz')
parser.add_argument('-th',  type=float, help='Summary coverage from 96 to 120 nt before stop', default=2)
parser.add_argument('-nt',  type=int, help="positions after stop codon - 3' UTR" , default=30)
parser.add_argument('-th3utr',  type=float, help="Summary coverage of 3' UTR - default 30nt after stop", default=0.2)
parser.add_argument('-span',  type=int, help='120 nt before stop - fixed, dont change!', default=120)
parser.add_argument('-col',  type=str, help='column name in hd5 used for coverage - default is "sum"', default="sum")
parser.add_argument('-subsets',  type=str, help='Split to subsets based Stop codon - ToDo!', default='NA')
parser.add_argument('-glist',  type=str, help='filename with gene names in interest one per line - ToDo!', default=None)
args = parser.parse_args()


print("\n\
-i           input *.hd5: {}\n\
-o          output *.txt: {}\n\
-annot    annotation GTF: {}\n\
-col       column in hd5: {}\n\
-nt        nt after stop: {}\n\
-th   gene bkg threshold: {}\n\
-th3utr  thr. after stop: {}\n".format(args.i, args.o, args.annot, args.col, args.nt, args.th, args.th3utr))


usage = "./compute_relative_3utr_coverage.py  -i S1_gene_coverage.hdf  -o relative_3utr_coverage.txt"
description = "Add some description about script. ToDo! "

if (args.i==None)|(args.o==None):
     sys.exit("\n  usage:\n\t{}\n\n\t{}\n".format(usage, description))

infile_h5  = args.i
outfile    = args.o
thres_cover= args.th 
thres_3utr = args.th3utr
after_stop = args.nt
Span       = args.span     
annotation = args.annot
col        = args.col
     
#
#
def df_framing(df1, index, columns, strand="+"):
    # create df2
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
#
#
c=c1=c2=c3=0# counting 
d = {}  # dictionary for collecting data
 
# input h5
hd5 = pd.HDFStore(infile_h5, "r")  
# metagene summary df
columns = hd5[hd5.keys()[0]].columns

for ref in yeastChr():
    #print("{}...".format(ref))
    df_f = hd5['For_rpm/'+ ref]
    df_r = hd5['Rev_rpm/'+ ref]

    stop_gtf_l = get_part_from_gtf(annotation, reference=ref, feature="stop_codon") # gtf part for stop

    for gtf in stop_gtf_l:
        #stop_codon = genome[ref][gtf.start:gtf.start+3]                      # get stop codon
        #stop_codon = stop_codon if gtf.strand=='+' else revcompl(stop_codon) # revcomp rev_strand
        # First test coverage  !! below is second test with the same threshold but more specific region !!
        coverage_for = df_f.loc[gtf.start - Span:gtf.start + 3, col].sum()
        coverage_rev = df_r.loc[gtf.start:gtf.start + Span - 1, col].sum()
        coverage_rpm = coverage_for if gtf.strand == '+' else coverage_rev

        if coverage_rpm < thres_cover: # THRESHOLD 1 REGION 1
            c1+=1
            continue

        # get regions
        df_1    = df_f.loc[gtf.start - Span:gtf.start + Span, :]     #  Forvard
        index_1 = range(gtf.start - Span, gtf.start + Span + 1)
        df_2    = df_r.loc[gtf.end - Span - 1:gtf.end + Span - 1, :] # Reverse
        index_2 = range(gtf.end - Span - 1, gtf.end + Span)

        # chose proper index & dataframe
        df, index = (df_1,index_1) if gtf.strand == '+' else (df_2, index_2)  # select one based strand
        df = df_framing(df, index=index, columns=columns, strand=gtf.strand)  # expanded & index resetted df
        df['rel_Pos'] = list(range(-Span, Span + 1))
        df.set_index('rel_Pos', inplace=True)

        #bkg_30_60 = df.loc[-57:-39,  col].mean()*3  # 9 codons coverage mean between -60 & -30 peaks
        #bkg_60_90 = df.loc[-87:-69,  col].mean()*3
        
        # sum of region
        bkg_90_120_sum = df.loc[-120:-96,col].sum()
        utr3pr_sum     = df.loc[3:after_stop, col].sum()
        # mean of region
        bkg_90_120_u   = df.loc[-120:-96,col].mean()
        utr3pr_u       = df.loc[3:after_stop, col].mean()
        ###
        if bkg_90_120_sum < thres_cover: # THRESHOLD 1 REGION 2
            c1+=1
            continue
        
        if utr3pr_sum < thres_3utr: # THRESHOLD 3
            c2+=1
            continue
               
        # mean values for readthrough
        readthrough = utr3pr_u/bkg_90_120_u
        
        ##########
        # collect data
        c3+=1
        d[gtf.gene_id] = readthrough
        

print("{:5} failed F1 120nt".format(c1))
print("{:5} failed F1 -96 to 120nt".format(c2))
print("{:5} passed!".format(c3))

#convert to series
s = pd.Series(d)
# name is not stored when saving to *.csv
#t = infile_h5.split('/')[-1].split('_')
#name = "{}_{}".format(t[0],t[1]) 
#s.name=name
s.to_csv(outfile, sep="\t")

print("{:5} genes with readthrough scores".format(s.shape[0]))
print("{:.2f}  Total readthrough scores sum".format(s.sum()))
print("=> {}\n//".format(outfile))
