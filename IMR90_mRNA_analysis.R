setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/")

library(dplyr)
library(readr)

samplefile <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/samplefile_IMR90_gencode.txt")

library(tximportData)

dir <- "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90"

files <- file.path(dir, samplefile$sampleID, "quant.sf")

names(files) <- paste0(samplefile$sample)

all(file.exists(files))


library(rtracklayer)
library(tximport)

gtf_file <- file.path(dir, "gencode.v34.annotation.gtf")

gtf_content <- import(gtf_file, feature.type = "gene")  

annotation <- data.frame(elementMetadata(gtf_content), stringsAsFactors = FALSE)

# set data frame row names to the gene IDs 

rownames(annotation) <- annotation$gene_id

annotation_transcript <- elementMetadata(import(gtf_file, feature.type = "transcript"))

tx2gene <- annotation_transcript[,c("transcript_id", "gene_id")]

# load the salmon dataset 

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)

head(txi.salmon$counts)

# construct a DESeqDataSet from the txi object and sample information in samples

library(DESeq2)

ddsTxi_obj <- DESeqDataSetFromTximport(txi.salmon,
                                       colData = samplefile,
                                       design = ~ batch + condition)

#################################################################################################################################

# Quality control

colData(ddsTxi_obj)

tpm <- txi.salmon$abundance

sum(assay(ddsTxi_obj)[,1])

colSums(assay(ddsTxi_obj))



# filtering non-expressed genes

is_expressed <- assay(ddsTxi_obj) >= 5

head(is_expressed)

hist(rowSums(is_expressed),main="Number of samples a gene is expressed in",xlab="Sample Count")

keep <- rowSums(assay(ddsTxi_obj) >= 5) >= 2

table(keep)

ddsTxi_obj <- ddsTxi_obj[keep,]


# visualising counts distribution

# Get log2 counts

vsd <- vst(ddsTxi_obj, blind=TRUE)

# Check distributions of samples using boxplots

png("Normalizedcounts_boxplot.png")

boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")

# Let's add a blue horizontal line that corresponds to the median logCPM

abline(h=median(assay(vsd)), col="blue")

dev.off()


#################################################################################################################################

# carry out the differential expression with apeglm log fold shrinkage

# had to relevel to make Empty the reference for comparison ... normally done alphabetically

ddsTxi_obj$condition <- relevel(ddsTxi_obj$condition, ref = "Empty")

ddsTxi <- DESeq(ddsTxi_obj)

res <- results(ddsTxi)

resultsNames(ddsTxi)


# also had to redo the binomial wald test 


ddsTxi <- nbinomWaldTest(ddsTxi)

resultsNames(ddsTxi)

OIS_IMR <- lfcShrink(ddsTxi, coef = "condition_OIS_vs_Empty", type="apeglm") 

#################################################################################################################################

# p-values and adjusted p-values

# We can order our results table by the smallest p value

resOrdered <- OIS_IMR[order(OIS_IMR$pvalue),]

# We can summarize some basic tallies using the summary function. 

summary(res)


# How many adjusted p-values were less than 0.05?

sum(OIS_IMR$padj < 0.05, na.rm=TRUE)


# MA-plot

png("MA-plot.png")

plotMA(OIS_IMR)

dev.off()

# this gives log2(n + 1)

ntd <- normTransform(ddsTxi_obj)

#################################################################################################################################

library(ggplot2)

# PCA post batch correction 

ddsTxi_obj$batch <- factor(rep(c(1, 2, 3, 4, 5, 6, 7)))

png("PCAplot_pre_batchcorrection.png")

vsd <- varianceStabilizingTransformation(ddsTxi_obj) # do a variance stablizing transformation on the model matrix 

plotPCA(vsd, "batch") +
  geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5) # this should give you the variance associated with the effect of the batch 

dev.off()

png("PCAplot_post_batchcorrection.png")

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch) # applying limma's removeBatchEffect 

plotPCA(vsd, intgroup = c("condition", "batch")) +
  geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5)+
  coord_fixed(ratio = 100)  # plotting the residuals 

dev.off()

# Replot the PCAs in ggplot 

# need to find way to replot the PCA before and after batch correction 

#################################################################################################################################

# P-VALUES and Adj p-values HISTOGRAM SANITY CHECK 

png("pvalues_histogram.png")

hist(OIS_IMR$pvalue)

dev.off()

png("padj_histogram.png")

hist(OIS_IMR$padj)

dev.off()

#################################################################################################################################
# Extract the DEseq2 Results 

library(dplyr)
library(tidyverse)

OIS_IMR_table <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("GeneID") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj)

# Exctract columns of interest from annotation

deseq_annotation <- as.data.frame(annotation) %>%
  select(-c(source, type, score, phase, gene_id, level, hgnc_id, havana_gene, tag)) %>%
  rownames_to_column("GeneID") 

OIS_IMR_table <- OIS_IMR_table %>%
  left_join(deseq_annotation, "GeneID")

write_tsv(OIS_IMR_table, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/OIS_DESeq2resultstable.csv")


#################################################################################################################################

# volcano plot (both figure and interactive)

# Enhanced Volcano 

library(ggrepel)
library(EnhancedVolcano)

png("OIS_IMR_volcano.png")

EnhancedVolcano(OIS_IMR_table,
                lab = OIS_IMR_table$gene_name,
                x = "log2FC",
                y = "FDR",
                ylab = "-Log10(FDR)",
                #xlim = c(-8, 8), 
                title = "Oncogene Induced Senescence (ER:Ras) in IMR90s", 
                border = "full",
                legendPosition = "right",
                legend = c("NS", "Log2FC", "FDR", "Log2FC & FDR"),
                legendLabSize = 16,
                legendIconSize = 5.0,
                gridlines.minor = TRUE)

dev.off()

#################################################################################################################################

# Interactive StripChart with Glimma

library(Glimma)

group <- samplefile$sample

de <- as.integer(OIS_IMR$padj <= 0.05 & !is.na(OIS_IMR$padj))

normCounts <- log2(counts(ddsTxi_obj))

statusCol <- match(samplefile$condition, c("Empty", "OIS"))

glXYPlot(
  x = OIS_IMR$log2FoldChange,
  y = -log10(OIS_IMR$pvalue),
  xlab = "log2FC",
  ylab = "-log10(FDR)",
  side.xlab = "Condition",
  side.ylab = "Log(CPM)",
  main = "Oncogene Induced Senescence (ER:RAS) in IMR90s",
  counts = normCounts,
  groups = group,
  sample.cols = statusCol,
  status = de,
  anno = OIS_IMR_table[, c("GeneID", "gene_name")],
  folder = "interactive_volcano")

#################################################################################################################################

# HEATMAPS

library(ComplexHeatmap)
library(circlize)

# get the top 200 genes
sigGenes <- as.data.frame(OIS_IMR_table) %>% 
  top_n(200, Empty=-FDR) %>% 
  pull("GeneID")

# filter the data for the top 200 by padj in the LRT test
plotDat <- varianceStabilizingTransformation(ddsTxi_obj)[sigGenes,] %>% 
  assay()
z.mat <- t(scale(t(plotDat), center=TRUE, scale=TRUE))

# colour palette
myPalette <- c("red3", "ivory", "blue3")
myRamp = colorRamp2(c(-2, 0, 2), myPalette)

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_names = FALSE,
        cluster_columns = FALSE)

# cluster the data and split the tree
hcDat <- hclust(dist(z.mat))
cutGroups <- cutree(hcDat, h=4)

ha1 = HeatmapAnnotation(df = colData(ddsTxi_obj)[,c("condition", "batch")])

png("stringent_heatmap.png")

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_name = FALSE,
        cluster_columns = FALSE,
        split=cutGroups,
        rect_gp = gpar(col = "darkgrey", lwd=0.5),
        top_annotation = ha1)

dev.off()

#################################################################################################################################

#################################################################################################################################
# Find the most abundant transcripts FOR EACH CONDITION by tpm

# import all the salmon files in

transcript_salmon <- tximport(files, type = "salmon", txOut = TRUE)

# take the abundances for each replicate 

transcript_tpm <- transcript_salmon$abundance

# change the rownames to a column called transcript_id

transcript_tpm <- as.data.frame(transcript_tpm) %>%
  rownames_to_column("transcript_id")

# add mean column ... remember, to do the SCU you need to do this for each condition

abundant_transcripts_Empty <- transcript_tpm %>%
  select(-c(Ras1, Ras2, Ras3, Ras4)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(Empty1, Empty2, Empty3))) %>%
  filter(mean_tpm != 0)

abundant_transcripts_OIS <- transcript_tpm %>%
  select(-c(Empty1, Empty2, Empty3)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(Ras1, Ras2, Ras3, Ras4))) %>%
  filter(mean_tpm != 0)

# find the mean raw counts for each of the abundant transcripts 

#meanrawcounts_Empty <- as.data.frame(transcript_salmon$counts) %>%
#  rownames_to_column("transcript_id") %>%
#  select(-c(Ras1, Ras2, Ras3, Ras4)) %>%
#  rowwise() %>%
#  mutate(mean_counts = round(mean(c(Empty1, Empty2, Empty3)))) %>%
#  select(-c(Empty1, Empty2, Empty3)) %>%
#  filter(mean_counts != 0)

#meanrawcounts_OIS <- as.data.frame(transcript_salmon$counts) %>%
#  rownames_to_column("transcript_id") %>%
#  select(-c(Empty1, Empty2, Empty3)) %>%
#  rowwise() %>%
#  mutate(mean_counts = round(mean(c(Ras1, Ras2, Ras3, Ras4)))) %>%
# select(-c(Ras1, Ras2, Ras3, Ras4)) %>%
# filter(mean_counts != 0)  


# add the gene_id from the previously constructed tx 2 gene list

tx2gene_df <- as.data.frame(tx2gene)

abundant_transcripts_Empty <- dplyr::inner_join(abundant_transcripts_Empty, tx2gene_df, by = "transcript_id")

abundant_transcripts_OIS <- dplyr::inner_join(abundant_transcripts_OIS, tx2gene_df, by = "transcript_id")

# select transcripts with the biggest mean tpm ... selecting the most abundant transcript in each gene set

abundant_transcripts_Empty_all <- abundant_transcripts_Empty %>%
  group_by(gene_id) %>%
  #top_n(n = 1, Empty = mean_tpm) %>%
  filter(mean_tpm != 0)

abundant_transcripts_OIS_all <- abundant_transcripts_OIS %>%
  group_by(gene_id) %>%
  #top_n(n = 1, Empty = mean_tpm) %>%
  filter(mean_tpm != 0)

abundant_transcripts_Empty_ab <- abundant_transcripts_Empty %>%
  group_by(gene_id) %>%
  top_n(n = 1, Empty = mean_tpm) 

abundant_transcripts_OIS_ab <- abundant_transcripts_OIS %>%
  group_by(gene_id) %>%
  top_n(n = 1, Empty = mean_tpm) 

# dont delete this until i know i dont need it !!!!!
#transcripts_byTPM <- abundant_transcripts %>%
#select(c(gene_id, transcript_id))


# sanity check for duplications ... change according to condition

duplications <- abundant_transcripts_OIS_all%>%
  add_count(transcript_id) %>%
  filter(n>1) %>%
  distinct(transcript_id) %>%
  nrow()

duplications

# filter the mean raw counts table so they only include the most abundant transcripts for each gene 

#meanrawcounts_Empty <- dplyr::inner_join(meanrawcounts_Empty, abundant_transcripts_Empty, by = "transcript_id") %>%
#  select(-c(Empty1, Empty2, Empty3, mean_tpm, gene_id))

#meanrawcounts_OIS <- dplyr::inner_join(meanrawcounts_OIS, abundant_transcripts_OIS, by = "transcript_id") %>%
#  select(-c(Ras1, Ras2, Ras3, Ras4, mean_tpm, gene_id))

# Rasate CDS fasta file with all expressed genes from salmon mapping step  

# Rasate matrix with transcript tpms

tpm_Empty_all <- transcript_tpm %>%
  select(-c(Ras1, Ras2, Ras3, Ras4)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(Empty1, Empty2, Empty3))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(Empty1, Empty2, Empty3))

#tpm_Empty_all <- as.matrix(tpm_Empty_all)

tpm_OIS_all <- transcript_tpm %>%
  select(-c(Empty1, Empty2, Empty3)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(Ras1, Ras2, Ras3, Ras4))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(Ras1, Ras2, Ras3, Ras4))

#tpm_OIS_all <- as.matrix(tpm_OIS_all)

tpm_Empty_abundant <- transcript_tpm %>%
  dplyr::inner_join(filtered_abundantTx_Empty, by = "transcript_id") %>%
  select(-c(Ras1, Ras2, Ras3, Ras4, Empty1.y, Empty2.y, Empty3.y )) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(Empty1.x, Empty2.x, Empty3.x))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(Empty1.x, Empty2.x, Empty3.x)) 

#tpm_Empty_abundant <- as.matrix(tpm_Empty_abundant)

tpm_OIS_abundant <- transcript_tpm %>%
  dplyr::inner_join(filtered_abundantTx_OIS, by = "transcript_id") %>%
  select(-c(Empty1, Empty2, Empty3, Ras1.y, Ras2.y, Ras3.y, Ras4.y )) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(Ras1.x, Ras2.x, Ras3.x, Ras4.x))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(Ras1.x, Ras2.x, Ras3.x, Ras4.x))

#tpm_OIS_abundant <- as.matrix(tpm_OIS_abundant)

library(Biostrings)
library(dplyr)

# take the abundant transcripts and remove the columns i dont need 

filtered_abundantTx_Empty <- as.data.frame(abundant_transcripts_Empty_ab) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_abundantTx_OIS <- as.data.frame(abundant_transcripts_OIS_ab) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_allTx_Empty <- as.data.frame(abundant_transcripts_Empty_all) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_allTx_OIS <- as.data.frame(abundant_transcripts_OIS_all) %>%
  select(-c(mean_tpm, gene_id)) 


# remove any rows where transcripts are not being expressed in the condition 

#keep1 <- rowSums(filtered_abundantTx_Empty) > 0
#filtered_abundantTx_Empty <- filtered_abundantTx_Empty[keep1, ]

#keep2 <- rowSums(filtered_abundantTx_OIS) > 0
#filtered_abundantTx_OIS <- filtered_abundantTx_OIS[keep2, ]

# row names to column 

#filtered_abundantTx_Empty <-as.data.frame(abundant_transcripts_Empty) %>%
#  select(-c(mean_tpm, gene_id))

#filtered_abundantTx_OIS <-as.data.frame(abundant_transcripts_OIS) %>%
#  select(-c(mean_tpm, gene_id))

# load the transcriptome (filtered CDS)
# all these CDSs have 5' and 3' UTRs, Start/Stop codon and is divisible by 3 

filtered_CDSfasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/gencode.v34.pc_transcripts_CDSs.fa")

# Rasate a data frame with only the expressed transcripts 

Empty_CDS_allTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_allTx_Empty, by = "transcript_id") %>%
  select(-c(Empty1, Empty2, Empty3))

OIS_CDS_allTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_allTx_OIS, by = "transcript_id") %>%
  select(-c(Ras1, Ras2, Ras3, Ras4))

Empty_CDS_abundantTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_Empty, by = "transcript_id") %>%
  select(-c(Empty1, Empty2, Empty3))

OIS_CDS_abundantTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_OIS, by = "transcript_id") %>%
  select(-c(Ras1, Ras2, Ras3, Ras4))

# sanity check 

#dim(Empty_CDS)
#dim(meanrawcounts_Empty)

#dim(OIS_CDS)
#dim(meanrawcounts_OIS)

# write the newly generated transcriptome dataframes into fasta and save on disk 

library(tidyverse)
library(seqinr)

# remeber to set a new wd 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/")

df1 = Empty_CDS_allTx
seqs = as.list(dplyr::pull(df1, CDS))
names = dplyr::pull(df1, `transcript_id`)
write.fasta(seqs, names, "Empty_CDS_allTx.fasta",
            open = "w", as.string = FALSE)

df2 = OIS_CDS_allTx
seqs = as.list(dplyr::pull(df2, CDS))
names = dplyr::pull(df2, `transcript_id`)
write.fasta(seqs, names, "OIS_CDS_allTx.fasta",
            open = "w", as.string = FALSE)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/")

df3 = Empty_CDS_abundantTx
seqs = as.list(dplyr::pull(df3, CDS))
names = dplyr::pull(df3, `transcript_id`)
write.fasta(seqs, names, "Empty_CDS.fasta",
            open = "w", as.string = FALSE)

df4 = OIS_CDS_abundantTx
seqs = as.list(dplyr::pull(df4, CDS))
names = dplyr::pull(df4, `transcript_id`)
write.fasta(seqs, names, "OIS_CDS.fasta",
            open = "w", as.string = FALSE)

#################################################################################################################################

# Calculate the SCU

# double check that the directory has no dodgy fasta files in there otherwise this might mess up the calculations

library("seqinr")

aa <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/codon_aa_table.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)

head(aa)

## Empty All Transcripts

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/"

out.file <- ""

file.names <- dir(path, pattern ="Empty_CDS_allTx.fasta")

table <- Empty

# calculate the codon frequency 

for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)

write.table(table, file = "Empty_allTx_SCUoutput.txt", sep="\t")

##############################################################################
## OIS All Transcripts 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/"

out.file <- ""

file.names <- dir(path, pattern ="OIS_CDS_allTx.fasta")

table <- Empty

# calculate the codon frequency 

for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)

write.table(table, file = "OIS_allTx_SCUoutput.txt", sep="\t")

###############################################################################################################################

## Empty Abundant Transcripts

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/"

out.file <- ""

file.names <- dir(path, pattern ="Empty_CDS.fasta")

table <- Empty

# calculate the codon frequency 

for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)

write.table(table, file = "Empty_abTx_SCUoutput.txt", sep="\t")

##############################################################################
## OIS Abundant Transcripts 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/"

out.file <- ""

file.names <- dir(path, pattern ="OIS_CDS.fasta")

table <- Empty

# calculate the codon frequency 

for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)

write.table(table, file = "OIS_abTx_SCUoutput.txt", sep="\t")

###############################################################################################################################
# combine the generated tables and store them as dataframes 
# remember to reset the wd directory ... all figures generated will be saved there 

### ALL TRANSCRIPTS

library(readr)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/")

letter_to_amino <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/letter_to_amino.txt")

Empty_allTx_RSCU <- read.table("Empty_allTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon")

OIS_allTx_RSCU <- read.table("OIS_allTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon") 

EvOIS_allTx_RSCU <- Empty_allTx_RSCU %>%
  inner_join(OIS_allTx_RSCU, by = "codon")

EvOIS_allTx_RSCU <- EvOIS_allTx_RSCU %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(EvOIS_allTx_RSCU, file = "EvOIS_allTx_SCUoutput.txt", sep="\t")

### ABUNDANT TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/")

Empty_abTx_RSCU <- read.table("Empty_abTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon")

OIS_abTx_RSCU <- read.table("OIS_abTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon") 

EvOIS_abTx_RSCU <- Empty_abTx_RSCU %>%
  inner_join(OIS_abTx_RSCU, by = "codon")

EvOIS_abTx_RSCU <- EvOIS_abTx_RSCU %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(EvOIS_abTx_RSCU, file = "EvOIS_allTx_SCUoutput.txt", sep="\t")

###############################################################################################################################

# Generate gingold type plots for both all and abundant transcripts 

### ALL TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/")

df_melt_all <- mutate(EvOIS_allTx_RSCU,
                      wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                         str_detect(codon, "G$|C$") ~ "GC_wobble"))

# gingold type plot

allTx_gingoldplot <- ggplot(df_melt_all, aes(x=Empty_CDS_allTx, y= OIS_CDS_allTx)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (all transcripts)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EvOIS_allTx_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)


allTx_gingoldplot

### ABUNDANT TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/")

df_melt_abundant <- mutate(EvOIS_abTx_RSCU,
                           wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                              str_detect(codon, "G$|C$") ~ "GC_wobble"))

# gingold type plot

abTx_gingoldplot <- ggplot(df_melt_abundant, aes(x=Empty_CDS, y= OIS_CDS)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (abundant transcripts only)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EvOIS_abTx_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)


abTx_gingoldplot

###############################################################################################################################
# Sperate up and down regulated genes 
# Empty (Ctrl) = downregulated
# OIS (condition) = upregulated 

# first i need to take the most abundant transcript for all genes and put them in a data frame 
# this will be used to convert the gene_id back to their transcript ID 

# find the most abundant transcript (but for all conditions instead of for each condition)

abundant_transcripts_everything <- transcript_tpm %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(Empty1, Empty2, Empty3, Ras1, Ras2, Ras3, Ras4))) %>%
  filter(mean_tpm != 0)

# get the gene id for each transcript 

abundant_transcripts_everything <- dplyr::inner_join(abundant_transcripts_everything, tx2gene_df, by = "transcript_id")

# take the most abundant transcript for each gene 

abundant_transcripts_everything <- abundant_transcripts_everything %>%
  group_by(gene_id) %>%
  top_n(n = 1, Empty = mean_tpm) %>%
  filter(mean_tpm != 0) %>%
  select(-c(Empty1, Empty2, Empty3, Ras1, Ras2, Ras3, Ras4)) %>%
  ungroup() %>%
  select(-c(gene_id))

# this is what will be used to weight the abundance for each gene 

##### set the working directory 

Empty_gene2abTx <- abundant_transcripts_Empty_ab %>%
  select(-c(Empty1, Empty2, Empty3))

OIS_gene2abTx <- abundant_transcripts_OIS_ab %>%
  select(-c(Ras1, Ras2, Ras3, Ras4))

upregulated_cutoff0 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC > 0) %>%
  dplyr::inner_join(OIS_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff0.5 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 0.5) %>%
  dplyr::inner_join(OIS_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff1 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 1) %>%
  dplyr::inner_join(OIS_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff2 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 2) %>%
  dplyr::inner_join(OIS_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff0 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC < 0) %>%
  dplyr::inner_join(Empty_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff0.5 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -0.5) %>%
  dplyr::inner_join(Empty_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff1 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -1) %>%
  dplyr::inner_join(Empty_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff2 <- as.data.frame(OIS_IMR) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -2) %>%
  dplyr::inner_join(Empty_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

###################################

Empty_differential_cutoff0 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff0, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Empty_differential_cutoff0.5 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff0.5, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Empty_differential_cutoff1 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff1, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Empty_differential_cutoff2 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff2, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

OIS_differential_cutoff0 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff0, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

OIS_differential_cutoff0.5 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff0.5, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

OIS_differential_cutoff1 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff1, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

OIS_differential_cutoff2 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff2, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

# will also use log2FC cutoffs of 0, 0.5, 1 and 2 

##########################

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_differentials")

# write the fasta files in the appropriate directory

df5 = Empty_differential_cutoff0
seqs = as.list(dplyr::pull(df5, CDS))
names = dplyr::pull(df5, `transcript_id`)
write.fasta(seqs, names, "Empty_differential_cutoff0.fasta",
            open = "w", as.string = FALSE)

df6 = Empty_differential_cutoff0.5
seqs = as.list(dplyr::pull(df6, CDS))
names = dplyr::pull(df6, `transcript_id`)
write.fasta(seqs, names, "Empty_differential_cutoff0.5.fasta",
            open = "w", as.string = FALSE)

df7 = Empty_differential_cutoff1
seqs = as.list(dplyr::pull(df7, CDS))
names = dplyr::pull(df7, `transcript_id`)
write.fasta(seqs, names, "Empty_differential_cutoff1.fasta",
            open = "w", as.string = FALSE)

df8 = Empty_differential_cutoff2
seqs = as.list(dplyr::pull(df8, CDS))
names = dplyr::pull(df8, `transcript_id`)
write.fasta(seqs, names, "Empty_differential_cutoff2.fasta",
            open = "w", as.string = FALSE)


df9 = OIS_differential_cutoff0
seqs = as.list(dplyr::pull(df9, CDS))
names = dplyr::pull(df9, `transcript_id`)
write.fasta(seqs, names, "OIS_differential_cutoff0.fasta",
            open = "w", as.string = FALSE)

df10 = OIS_differential_cutoff0.5
seqs = as.list(dplyr::pull(df10, CDS))
names = dplyr::pull(df10, `transcript_id`)
write.fasta(seqs, names, "OIS_differential_cutoff0.5.fasta",
            open = "w", as.string = FALSE)

df11 = OIS_differential_cutoff1
seqs = as.list(dplyr::pull(df11, CDS))
names = dplyr::pull(df11, `transcript_id`)
write.fasta(seqs, names, "OIS_differential_cutoff1.fasta",
            open = "w", as.string = FALSE)

df12 = OIS_differential_cutoff2
seqs = as.list(dplyr::pull(df12, CDS))
names = dplyr::pull(df12, `transcript_id`)
write.fasta(seqs, names, "OIS_differential_cutoff2.fasta",
            open = "w", as.string = FALSE)

########################################################

# calculate the RSCU

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_differentials/")
path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_differentials/"
out.file<-""
file.names <- dir(path, pattern ="OIS_differential_cutoff0.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "OIS_differential_cutoff0_SCUoutput.txt",sep="\t")

########################################################################################
out.file<-""
file.names <- dir(path, pattern ="OIS_differential_cutoff0.5.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "OIS_differential_cutoff0.5_SCUoutput.txt",sep="\t")

###############################################################################

out.file<-""
file.names <- dir(path, pattern ="OIS_differential_cutoff1.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "OIS_differential_cutoff1_SCUoutput.txt",sep="\t")

#################################################################################

out.file<-""
file.names <- dir(path, pattern ="OIS_differential_cutoff2.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "OIS_differential_cutoff2_SCUoutput.txt",sep="\t")

###############################################################################################

out.file<-""
file.names <- dir(path, pattern ="Empty_differential_cutoff0.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "Empty_differential_cutoff0_SCUoutput.txt",sep="\t")

#############################################################################################################

out.file<-""
file.names <- dir(path, pattern ="Empty_differential_cutoff0.5.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "Empty_differential_cutoff0.5_SCUoutput.txt",sep="\t")

#################################################################################################################

out.file<-""
file.names <- dir(path, pattern ="Empty_differential_cutoff1.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "Empty_differential_cutoff1_SCUoutput.txt",sep="\t")

#######################################################################################################

out.file<-""
file.names <- dir(path, pattern ="Empty_differential_cutoff2.fasta")
table<- Empty
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name <- file.names[i]
  name <- gsub(".fasta","",name)
  print(name)
  UCO <- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1 <- do.call(rbind,UCO)
  UCO1 <- as.data.frame(UCO1)
  colnames(UCO1) <- aa$AA
  head(UCO1) 
}

# multiply the codon frequency by the mean_tpm to weight for abundance 

library(tidyverse)

UCO2 <- UCO1 %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(UCO1)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm)

colnames(UCO2) <- str_replace(colnames(UCO2),"\\..*","")

UCO2 <- UCO2 %>%
  column_to_rownames(var = "transcript_id")

# find the proportion of each codon for every amino acid (codon freq / sum freq)

x<- as.matrix(UCO2, rownames.force = TRUE)

head(x)

x1 <- aggregate(t(x), by=list(rownames(t(x))), FUN=sum)

head(x1)

rownames(x1) <- x1$Group.1

x1$Group.1 <- Empty

x1 <- t(x1)

head(x1)

rownames(UCO2) = make.names(rownames(UCO2), unique=TRUE)

head(UCO2)

rownames(x1) = make.names(rownames(x1), unique=TRUE)

a <- UCO2[1,] / x1[1,] [ match( colnames(UCO2) , colnames(x1) ) ]

b <- data.frame()

for(i in 1:nrow(UCO2)){
  
  div <- UCO2[i,] / x1[i,] [ match( colnames(UCO2) , colnames(x1) ) ]
  
  b <- rbind(b,div)
}

head(b)

colnames(b)<- aa$codon

freq<- colMeans(b,na.rm = TRUE)

table<-cbind(table,freq)

colnames(table)[ncol(table)] <- paste(name)

head(table)
write.table(table, file = "Empty_differential_cutoff2_SCUoutput.txt",sep="\t")

########################################################

# combine all the generated tables 

OIS_cutoff0 <- read.table("OIS_differential_cutoff0_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
OIS_cutoff0.5 <- read.table("OIS_differential_cutoff0.5_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
OIS_cutoff1 <- read.table("OIS_differential_cutoff1_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
OIS_cutoff2 <- read.table("OIS_differential_cutoff2_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")

Empty_cutoff0 <- read.table("Empty_differential_cutoff0_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Empty_cutoff0.5 <- read.table("Empty_differential_cutoff0.5_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Empty_cutoff1 <- read.table("Empty_differential_cutoff1_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Empty_cutoff2 <- read.table("Empty_differential_cutoff2_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")

EmptyvsOIS_differential_SCUoutput <- OIS_cutoff0 %>%
  dplyr::inner_join(OIS_cutoff0.5) %>%
  dplyr::inner_join(OIS_cutoff1) %>%
  dplyr::inner_join(OIS_cutoff2) %>%
  dplyr::inner_join(Empty_cutoff0) %>%
  dplyr::inner_join(Empty_cutoff0.5) %>%
  dplyr::inner_join(Empty_cutoff1) %>%
  dplyr::inner_join(Empty_cutoff2) 

# Add AA features to the RSCU output

EmptyvsOIS_differential_SCUoutput <- EmptyvsOIS_differential_SCUoutput %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(EmptyvsOIS_differential_SCUoutput, file = "EmptyvsOIS_differential_SCUoutput.txt", sep="\t")

# plot the gingold type plots 

df_melt_differentials <- mutate(EmptyvsOIS_differential_SCUoutput,
                                wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                                   str_detect(codon, "G$|C$") ~ "GC_wobble"),
                                wobble_single = case_when(str_detect(codon, "A$") ~ "A_wobble",
                                                          str_detect(codon, "T$") ~ "T_wobble",
                                                          str_detect(codon, "G$") ~ "G_wobble",
                                                          str_detect(codon, "C$") ~ "C_wobble"))


cutoff0_gingoldplot <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff0, y= OIS_differential_cutoff0)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = 0)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff0_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0_gingoldplot

cutoff0.5_gingoldplot <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff0.5, y= OIS_differential_cutoff0.5)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±0.5)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff0.5_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0.5_gingoldplot

cutoff1_gingoldplot <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff1, y= OIS_differential_cutoff1)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±1)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff1_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff1_gingoldplot

cutoff2_gingoldplot <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff2, y= OIS_differential_cutoff2)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±2)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff2_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff2_gingoldplot

###############################################################################################################################
# now do the same but highlight by single wobble

cutoff0_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff0, y= OIS_differential_cutoff0)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = 0)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  #geom_text_repel(aes(label = aminoacid)) +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff0_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0_gingoldplot_wobble

cutoff0.5_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff0.5, y= OIS_differential_cutoff0.5)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±0.5)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff0.5_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0.5_gingoldplot_wobble

cutoff1_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff1, y= OIS_differential_cutoff1)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±1)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff1_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff1_gingoldplot_wobble

cutoff2_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Empty_differential_cutoff2, y= OIS_differential_cutoff2)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±2)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_EmptyvsOIS_differential_cutoff2_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff2_gingoldplot_wobble
###############################################################################################################################

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff0, y= OIS_differential_cutoff0)) +
    geom_point(aes(color = wobble_single)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff0.5_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyvsOIS_differential_cutoff0_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff0.5, y= OIS_differential_cutoff0.5)) +
    geom_point(aes(color = wobble_single)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log cutoff = ±0.5)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff0.5_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyvsOIS_differential_cutoff0.5_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}


for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff1_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff1, y= OIS_differential_cutoff1)) +
    geom_point(aes(color = wobble_single)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log cutoff = ±1)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff1_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyOIS_differential_cutoff1_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}



for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff2_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff2, y= OIS_differential_cutoff2)) +
    geom_point(aes(color = wobble_single)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log cutoff = ±2)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff2_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyvsOIS_differential_cutoff2_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

###############################################################################################################################
for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff0, y= OIS_differential_cutoff0)) +
    geom_point(aes(color = woble)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff0.5_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyvsOIS_differential_cutoff0_wobble_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff0.5, y= OIS_differential_cutoff0.5)) +
    geom_point(aes(color = wobble)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log cutoff = ±0.5)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff0.5_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyvsOIS_differential_cutoff0.5_wobble_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}


for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff1_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff1, y= OIS_differential_cutoff1)) +
    geom_point(aes(color = wobble)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log cutoff = ±1)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff1_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyOIS_differential_cutoff1_wobble_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}



for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff2_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Empty_differential_cutoff2, y= OIS_differential_cutoff2)) +
    geom_point(aes(color = wobble)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log cutoff = ±2)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff2_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_EmptyvsOIS_differential_cutoff2_wobble_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}





################################################################################################################################

# calculate AA freq for each gene and normalise to their abundance (TPM)

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Empty_allTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/Empty_CDS_allTx.fasta")

Ras_allTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_all/OIS_CDS_allTx.fasta")

# translate the DNA sequences into protein sequence 

Empty_allTx_translated_fasta <- Biostrings::translate(Empty_allTx_fasta) 

Ras_allTx_translated_fasta <- Biostrings::translate(Ras_allTx_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Empty_allTx_translated_freq <- Biostrings::alphabetFrequency(Empty_allTx_translated_fasta)

Ras_allTx_translated_freq <- Biostrings::alphabetFrequency(Ras_allTx_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Empty_allTx_translated_fasta <- as.data.frame(Empty_allTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Empty_allTx_translated_freq <- as.data.frame(Empty_allTx_translated_freq)

rownames(Empty_allTx_translated_freq) <- Empty_allTx_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Empty_allTx_translated_freq <- Empty_allTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Empty_allTx_translated_freq_normalized <- Empty_allTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Empty_allTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Empty_AA_freq_normalized <- Empty_allTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Empty_AA_freq_normalized <- as.data.frame(Empty_AA_freq_normalized) %>%
  rownames_to_column("AA") %>%
  rename(Empty_AA_freq = V1) 

###################################################################################################

Ras_allTx_translated_fasta <- as.data.frame(Ras_allTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Ras_allTx_translated_freq <- as.data.frame(Ras_allTx_translated_freq)

rownames(Ras_allTx_translated_freq) <- Ras_allTx_translated_fasta$transcript_id

Ras_allTx_translated_freq <- Ras_allTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Ras_allTx_translated_freq_normalized <- Ras_allTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Ras_allTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Ras_AA_freq_normalized <- Ras_allTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Ras_AA_freq_normalized <- as.data.frame(Ras_AA_freq_normalized) %>%
  rownames_to_column("AA") %>%
  rename(Ras_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

OIS_AA_freq_normalized <- Empty_AA_freq_normalized %>%
  dplyr::inner_join(Ras_AA_freq_normalized, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(OIS_AA_freq_normalized, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/codon_freq/OIS_AA_freq_normalized.txt")

###############################################################################################################################

# calculate AA freq for the most abundant transcript and normalise to their abundance (TPM)

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Empty_abTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/Empty_CDS.fasta")

Ras_abTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_abundant/OIS_CDS.fasta")

# translate the DNA sequences into protein sequence 

Empty_abTx_translated_fasta <- Biostrings::translate(Empty_abTx_fasta) 

Ras_abTx_translated_fasta <- Biostrings::translate(Ras_abTx_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Empty_abTx_translated_freq <- Biostrings::alphabetFrequency(Empty_abTx_translated_fasta)

Ras_abTx_translated_freq <- Biostrings::alphabetFrequency(Ras_abTx_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Empty_abTx_translated_fasta <- as.data.frame(Empty_abTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Empty_abTx_translated_freq <- as.data.frame(Empty_abTx_translated_freq)

rownames(Empty_abTx_translated_freq) <- Empty_abTx_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Empty_abTx_translated_freq <- Empty_abTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Empty_abTx_translated_freq_normalized <- Empty_abTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Empty_abTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Empty_AA_freq_normalized_ab <- Empty_abTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Empty_AA_freq_normalized_ab <- as.data.frame(Empty_AA_freq_normalized_ab) %>%
  rownames_to_column("AA") %>%
  rename(Empty_AA_freq = V1) 

###################################################################################################

Ras_abTx_translated_fasta <- as.data.frame(Ras_abTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Ras_abTx_translated_freq <- as.data.frame(Ras_abTx_translated_freq)

rownames(Ras_abTx_translated_freq) <- Ras_abTx_translated_fasta$transcript_id

Ras_abTx_translated_freq <- Ras_abTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Ras_abTx_translated_freq_normalized <- Ras_abTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Ras_abTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Ras_AA_freq_normalized_ab <- Ras_abTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Ras_AA_freq_normalized_ab <- as.data.frame(Ras_AA_freq_normalized_ab) %>%
  rownames_to_column("AA") %>%
  rename(Ras_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

OIS_AA_freq_normalized_ab <- Empty_AA_freq_normalized_ab %>%
  dplyr::inner_join(Ras_AA_freq_normalized_ab, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(OIS_AA_freq_normalized_ab, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/codon_freq/OIS_AA_freq_normalized_ab.txt")

###################################################################################################
# calculate AA frequency for differential genes

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Empty_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_differentials/Empty_differential_cutoff2.fasta")

Ras_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/SCU/SCU_differentials/OIS_differential_cutoff2.fasta")

# translate the DNA sequences into protein sequence 

Empty_diff_translated_fasta <- Biostrings::translate(Empty_diff_fasta) 

Ras_diff_translated_fasta <- Biostrings::translate(Ras_diff_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Empty_diff_translated_freq <- Biostrings::alphabetFrequency(Empty_diff_translated_fasta)

Ras_diff_translated_freq <- Biostrings::alphabetFrequency(Ras_diff_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Empty_diff_translated_fasta <- as.data.frame(Empty_diff_translated_fasta) %>%
  rownames_to_column("transcript_id")

Empty_diff_translated_freq <- as.data.frame(Empty_diff_translated_freq)

rownames(Empty_diff_translated_freq) <- Empty_diff_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Empty_diff_translated_freq <- Empty_diff_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Empty_diff_translated_freq_normalized <- Empty_diff_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Empty_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Empty_diff_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  na.omit()

Empty_AA_freq_normalized_diff <- Empty_diff_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Empty_AA_freq_normalized_diff <- as.data.frame(Empty_AA_freq_normalized_diff) %>%
  rownames_to_column("AA") %>%
  rename(Empty_AA_freq = V1) 

###################################################################################################

Ras_diff_translated_fasta <- as.data.frame(Ras_diff_translated_fasta) %>%
  rownames_to_column("transcript_id")

Ras_diff_translated_freq <- as.data.frame(Ras_diff_translated_freq)

rownames(Ras_diff_translated_freq) <- Ras_diff_translated_fasta$transcript_id

Ras_diff_translated_freq <- Ras_diff_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Ras_diff_translated_freq_normalized <- Ras_diff_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_OIS_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Ras_abTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  na.omit()

Ras_AA_freq_normalized_diff <- Ras_diff_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Ras_AA_freq_normalized_diff <- as.data.frame(Ras_AA_freq_normalized_diff) %>%
  rownames_to_column("AA") %>%
  rename(Ras_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

OIS_AA_freq_normalized_diff <- Empty_AA_freq_normalized_diff %>%
  dplyr::inner_join(Ras_AA_freq_normalized_diff, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(OIS_AA_freq_normalized_diff, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/codon_freq/OIS_AA_freq_normalized_differential_logFC2cutoff.txt")


###################################################################################################

# AA property plots

aa_properties <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/aa_properties.txt", skip=1)

colnames(aa_properties) <- c("aminoacid", "property")

property_plot_input <- OIS_AA_freq_normalized_diff %>%
  left_join(aa_properties, "aminoacid")
  
ggplot(data = property_plot_input, aes(x=Empty_AA_freq, y= Ras_AA_freq))+
  geom_point(aes(color = property)) +
  theme_bw(base_size = 16)+
  labs(title = "Amino Acid Frequency in IMR90 Empty vs Ras") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  xlab("Empty AA frequency")+
  ylab("Ras AA frequency") +
  geom_label_repel(aes(label = aminoacid),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()  

ggsave(filename = "IMR90_aa_properties.png", bg = "white", width = 7, height = 7, dpi = 600)











