#!/usr/bin/env python
#
# works fine
# Features to add:
#       a) split by stop codon and  
#       b) work with gene list 
# 
#
__author__		= "Tonu Margus"
__copyright__	= "Copyright 2019"
__version__		= "0.0.1"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import sys
import argparse
import pandas as pd
import pysam


parser = argparse.ArgumentParser(description='Computes waviness')
parser.add_argument('-i', type=str, help='input table of gene coverage in *.hd5 format')
parser.add_argument('-o', type=str, help='output file name  *.csv')
parser.add_argument('-annot',  type=str, help='GTF annotation file', default='0-References/genome.gtf.gz')
parser.add_argument('-th1',  type=float, help='Summary gene coverage 120 nt before stop', default=10)
parser.add_argument('-th2',  type=float, help='Background coverage - codon mean between expected peak areas', default=0.2)
parser.add_argument('-span',  type=int, help='Positions before stop. Can be bigger than 120 nt but not shorter', default=120)
parser.add_argument('-col',  type=str, help='column name in hd5 used for coverage - default is "sum"', default="sum")
parser.add_argument('-subsets',  type=str, help='Split to subsets based Stop codon', default='NA')
args = parser.parse_args()


print("\n\
-i           input *.h5: {}\n\
-o         output *.csv: {}\n\
-annot   annotation GTF: {}\n\
-col             column: {}\n\
-th1   region threshold: {}\n\
-th2      bkg threshold: {}\n\
-span    nt before stop: {}\n".format(args.i, args.o, args.annot, args.col, args.th1, args.th2, args.span))

usage = "./compute_waviness.py -i *.hdf  -o *.csv"

if (args.i==None)|(args.o==None):
     sys.exit("\n  usage:\n\t{}\n".format(usage))

infile_h5  = args.i
outfile    = args.o
thres_cover= args.th1 
thres_bkg  = args.th2
Span       = args.span     
annotation = args.annot
col        = args.col

replacestr = "_th{}-{}.csv".format(thres_cover, thres_bkg)
outfile = outfile.replace('.csv', replacestr)

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
        #test coverage
        coverage_for = df_f.loc[gtf.start - Span:gtf.start + 3, col].sum()
        coverage_rev = df_r.loc[gtf.start:gtf.start + Span - 1, col].sum()
        coverage_rpm = coverage_for if gtf.strand == '+' else coverage_rev

        if coverage_rpm < thres_cover: #FILTER 1
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

        bkg_30_60 = df.loc[-57:-39,  col].mean()*3  # 9 codons coverage mean between -60 & -30 peaks
        bkg_60_90 = df.loc[-87:-69,  col].mean()*3
        bkg_90_120= df.loc[-117:-109,col].mean()*3
        
        # codon instead .mean()*3 use .sum()
        peak_30z  = df.loc[-36:-34,  col].sum()
        peak_30a  = df.loc[-33:-31,  col].sum()
        peak_30b  = df.loc[-30:-28,  col].sum()
        peak_30c  = df.loc[-27:-25,  col].sum()

        peak_60z  = df.loc[-66:-64,  col].sum()
        peak_60a  = df.loc[-63:-61,  col].sum()
        peak_60b  = df.loc[-60:-58,  col].sum()
        peak_60c  = df.loc[-57:-55,  col].sum()
        
        peak_90z  = df.loc[-96:-94,  col].sum()
        peak_90a  = df.loc[-93:-91,  col].sum()
        peak_90b  = df.loc[-90:-88,  col].sum()
        peak_90c  = df.loc[-87:-85,  col].sum()
        ##########
        # collect data
        c2+=1
        d[gtf.gene_id] = [bkg_30_60,bkg_60_90,bkg_90_120, peak_30z, peak_30a, peak_30b, peak_30c, 
                          peak_60z,peak_60a,peak_60b, peak_60c, peak_90z,peak_90a, peak_90b,peak_90c ]
print(c2)
# pmake dataframe
index = ["bkg_30_60","bkg_60_90","bkg_90_120", "peak_30z","peak_30a", "peak_30b", "peak_30c","peak_60z", 
         "peak_60a", "peak_60b", "peak_60c", "peak_90z","peak_90a", "peak_90b","peak_90c" ]
df = pd.DataFrame(d, index=index).T
df['bkg'] = df.loc[:,["bkg_30_60","bkg_60_90","bkg_90_120"]].mean(axis=1)
# drop df.bkg<th2
print("{} bleow bkg threshold {}".format(df.loc[df.bkg<thres_bkg,:].shape[0], thres_bkg))
df = df.loc[df.bkg>thres_bkg,:]
df['peak30'] = df.loc[:,["peak_30z","peak_30a", "peak_30b", "peak_30c"]].max(axis=1)
df['peak60'] = df.loc[:,["peak_60z","peak_60a", "peak_60b", "peak_60c"]].max(axis=1)
df['peak90'] = df.loc[:,["peak_90z","peak_90a", "peak_90b", "peak_90c"]].max(axis=1)
df['peakSum']= df.loc[:, ['peak30','peak60','peak90']].sum(axis=1)
df['waviness']=df['peakSum']/df['bkg']

i = 10 
mask = df['waviness']>i
print("There are  {} genes in the table".format(df.shape[0]))
print("There are  {} most wavier genes {} times above bakground".format(df[mask].shape[0], i))
df.to_csv(outfile, sep='\t', header=True, index=True)
print("Output:\n   {}".format(outfile))
print("//")
