library( "DESeq2" )
library(ggplot2)

##########################################
#  Run DESeq analysis for RNA-Seq data  #  
##########################################

sampleFile = "htseq_count_RNASeq.csv" # RNA-Seq
annotation_file = "rnaseq-anno.txt"   # RNA-Seq

countData <- read.table(sampleFile,header=TRUE,row.names=1,comment.char="#",blank.lines.skip=TRUE)
cts <- as.matrix(countData,row.names="gene_id")

coldata = read.table(annotation_file, header=TRUE, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Prefiltering low coverage
dds <- dds[ rowSums(counts(dds)) > 10, ]

# define levels
dds$condition <- relevel(dds$condition, ref="WT")

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# reorder table
resOrdered <- res[order(res$padj),]
summary(res)

# write to file
write.csv(as.data.frame(resOrdered), 
          file="KO_vs_WT_RNASeq_DESeq2.csv")


##########################################
#  Run DESeq analysis for Ribo-Seq data  #  
##########################################


# comment uncomment relevant below
sampleFile = "htseq_count_RiboSeq.csv" # Ribo-Seq
annotation_file = "riboseq-anno.txt"  # Ribo-Seq

countData <- read.table(sampleFile,header=TRUE,row.names=1,comment.char="#",blank.lines.skip=TRUE)
cts <- as.matrix(countData,row.names="gene_id")

coldata = read.table(annotation_file, header=TRUE, row.names=1)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Prefiltering low coverage
dds <- dds[ rowSums(counts(dds)) > 10, ]

# define levels
dds$condition <- relevel(dds$condition, ref="WT")

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# reorder table
resOrdered <- res[order(res$padj),]
summary(res)

# write to file
write.csv(as.data.frame(resOrdered), 
          file="KO_vs_WT_RiboSeq_DESeq2.csv")


