setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/")

library(dplyr)
library(readr)

samplefile <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat/samplefile_bcat_cmyc_gencode.txt")

samplefile <- samplefile[-c(10), ] # removed bcat5 because of the original PCA 



library(tximportData)

library(tximportData)

dir <- "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat"

files <- file.path(dir, samplefile$sampleID, "quant.sf")

names(files) <- paste0(samplefile$sample)

all(file.exists(files))


library(rtracklayer)
library(tximport)

gtf_file <- file.path(dir, "gencode.vM25.annotation.gtf")

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

png("Normalizedcounts_boxplot_WTvsbcatMyc.png")

boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")

# Let's add a blue horizontal line that corresponds to the median logCPM

abline(h=median(assay(vsd)), col="blue")

dev.off()


#################################################################################################################################

# carry out the differential expression with apeglm log fold shrinkage

# had to relevel to make WT the reference for comparison ... normally done alphabetically

ddsTxi_obj$condition <- relevel(ddsTxi_obj$condition, ref = "WT")

ddsTxi <- DESeq(ddsTxi_obj)

res <- results(ddsTxi)

resultsNames(ddsTxi)


# also had to redo the binomial wald test 


ddsTxi <- nbinomWaldTest(ddsTxi)

resultsNames(ddsTxi)

WTvsbcat_cmyc <- lfcShrink(ddsTxi, coef = "condition_bcat_cmyc_vs_WT", type="apeglm") 


#################################################################################################################################

# p-values and adjusted p-values

# We can order our results table by the smallest p value

resOrdered <- WTvsbcat_cmyc[order(WTvsbcat_cmyc$pvalue),]

# We can summarize some basic tallies using the summary function. 

summary(WTvsbcat_cmyc)


# How many adjusted p-values were less than 0.05?

sum(WTvsbcat_cmyc$padj < 0.05, na.rm=TRUE)


# MA-plot

png("MA-plot_WTvsbcat_cmyc_WTvsbcatMyc.png")

plotMA(WTvsbcat_cmyc, ylim=c(-5, 5))

dev.off()


# this gives log2(n + 1)

ntd <- normTransform(ddsTxi_obj)

#################################################################################################################################

library(ggplot2)

# PCA post batch correction 

#ddsTxi_obj$batch <- factor(rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18)))
png("PCAplot_pre_batchcorrection_rmOutliers_WTvsbcatMyc.png")

vsd <- varianceStabilizingTransformation(ddsTxi_obj) # do a variance stablizing transformation on the model matrix 

plotPCA(vsd, "condition") +
  geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5)+ 
  coord_equal(clip = "off")

# this should give you the variance associated with the effect of the batch 

dev.off()

#png("PCAplot_post_batchcorrection.png")

#assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch) # applying limma's removeBatchEffect 

#plotPCA(vsd, intgroup = c("condition", "batch")) +
#geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5) +
#coord_equal(clip = "off")# plotting the residuals 

#dev.off()

# Replot the PCAs in ggplot 

# need to find way to replot the PCA before and after batch correction 

#################################################################################################################################

# P-VALUES and Adj p-values HISTOGRAM SANITY CHECK 

png("pvalues_histogram_WTvsbcatMyc.png")

hist(WTvsbcat_cmyc$pvalue)

dev.off()

png("padj_histogram_bcat_cmyc_WTvsbcatMyc.png")

hist(WTvsbcat_cmyc$padj)

dev.off()

#################################################################################################################################
# Extract the DEseq2 Results 

library(dplyr)
library(tidyverse)

WTvsbcat_cmyc_table <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("GeneID") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj)


# Exctract columns of interest from annotation

deseq_annotation <- as.data.frame(annotation) %>%
  select(-c(source, type, score, phase, gene_id, level, havana_gene, tag, havana_gene, mgi_id)) %>%
  rownames_to_column("GeneID") 

WTvsbcat_cmyc_table <- WTvsbcat_cmyc_table %>%
  left_join(deseq_annotation, "GeneID")


write_tsv(WTvsbcat_cmyc_table, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/WTvsbcat_cmyc_DESeq2resultstable.txt")


#################################################################################################################################

# volcano plot (both figure and interactive)

# Enhanced Volcano 

library(ggrepel)
library(EnhancedVolcano)

png("WTvsbcat_cmyc_volcano_WTvsbcatMyc.png")

EnhancedVolcano(WTvsbcat_cmyc_table,
                lab = WTvsbcat_cmyc_table$gene_name,
                x = "log2FC",
                y = "FDR",
                ylab = "-Log10(FDR)",
                xlim = c(-8, 8), 
                title = "WT vs b-catenin/c-Myc Volcano Plot", 
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

de <- as.integer(WTvsbcat_cmyc$padj <= 0.05 & !is.na(WTvsbcat_cmyc$padj))

normCounts <- log2(counts(ddsTxi_obj))

#statusCol <- match(samplefile$condition, c("WT", "bcat_cmyc"))

glXYPlot(
  x = WTvsbcat_cmyc$log2FoldChange,
  y = -log10(WTvsbcat_cmyc$pvalue),
  xlab = "log2FC",
  ylab = "-log10(padj)",
  side.xlab = "Condition",
  side.ylab = "Log(CPM)",
  main = "WT vs b-catenin/c-Myc in Mouse Caudate Lobe",
  counts = normCounts,
  groups = group,
  #sample.cols = statusCol,
  status = de,
  anno = WTvsbcat_cmyc_table[, c("GeneID", "gene_name")],
  folder = "interactive_volcano")

#################################################################################################################################

# HEATMAPS

library(ComplexHeatmap)
library(circlize)

# get the top 200 genes
sigGenes <- as.data.frame(WTvsbcat_cmyc_table) %>% 
  top_n(200, wt=-FDR) %>% 
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

png("stringent_heatmap_WTvsbcat_myc_WTvsbcatMyc.png")

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

abundant_transcripts_WT <- transcript_tpm %>%
  select(-c(bcat1, bcat2, bcat3, bcat4, cmyc1, cmyc2, cmyc3, cmyc4, cmyc5, bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(WT1, WT2, WT3, WT4, WT5))) %>%
  filter(mean_tpm != 0)

abundant_transcripts_bcatMyc <- transcript_tpm %>%
  select(-c(WT1, WT2, WT3, WT4, WT5, bcat1, bcat2, bcat3, bcat4, cmyc1, cmyc2, cmyc3, cmyc4, cmyc5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5))) %>%
  filter(mean_tpm != 0)

# add the gene_id from the previously constructed tx 2 gene list

tx2gene_df <- as.data.frame(tx2gene)

abundant_transcripts_WT <- dplyr::inner_join(abundant_transcripts_WT, tx2gene_df, by = "transcript_id")

abundant_transcripts_bcatMyc <- dplyr::inner_join(abundant_transcripts_bcatMyc, tx2gene_df, by = "transcript_id")

# select transcripts with the biggest mean tpm ... selecting the most abundant transcript in each gene set

abundant_transcripts_WT_all <- abundant_transcripts_WT %>%
  group_by(gene_id) %>%
  #top_n(n = 1, wt = mean_tpm) %>%
  filter(mean_tpm != 0)

abundant_transcripts_bcatMyc_all <- abundant_transcripts_bcatMyc %>%
  group_by(gene_id) %>%
  #top_n(n = 1, wt = mean_tpm) %>%
  filter(mean_tpm != 0)

abundant_transcripts_WT_ab <- abundant_transcripts_WT %>%
  group_by(gene_id) %>%
  top_n(n = 1, wt = mean_tpm) 

abundant_transcripts_bcatMyc_ab <- abundant_transcripts_bcatMyc %>%
  group_by(gene_id) %>%
  top_n(n = 1, wt = mean_tpm) 

# dont delete this until i know i dont need it !!!!!
#transcripts_byTPM <- abundant_transcripts %>%
#select(c(gene_id, transcript_id))


# sanity check for duplications ... change according to condition

duplications <- abundant_transcripts_WT_all%>%
  add_count(transcript_id) %>%
  filter(n>1) %>%
  distinct(transcript_id) %>%
  nrow()

duplications

# filter the mean raw counts table so they only include the most abundant transcripts for each gene 

#meanrawcounts_Fed <- dplyr::inner_join(meanrawcounts_Fed, abundant_transcripts_Fed, by = "transcript_id") %>%
#  select(-c(Fed1, Fed2, Fed5, mean_tpm, gene_id))

#meanrawcounts_Starved <- dplyr::inner_join(meanrawcounts_Starved, abundant_transcripts_Starved, by = "transcript_id") %>%
#  select(-c(Starved1, Starved2, Starved5, mean_tpm, gene_id))

# create CDS fasta file with all expressed genes from salmon mapping step  

# take the abundant transcripts and remove the columns i dont need 

filtered_abundantTx_WT <- as.data.frame(abundant_transcripts_WT_ab) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_abundantTx_bcatMyc <- as.data.frame(abundant_transcripts_bcatMyc_ab) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_allTx_WT <- as.data.frame(abundant_transcripts_WT_all) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_allTx_bcatMyc <- as.data.frame(abundant_transcripts_bcatMyc_all) %>%
  select(-c(mean_tpm, gene_id)) 

# create matrix with transcript tpms

tpm_WT_all <- transcript_tpm %>%
  select(-c(bcat1, bcat2, bcat3, bcat4, cmyc1, cmyc2, cmyc3, cmyc4, cmyc5, bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(WT1, WT2, WT3, WT4, WT5))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(WT1, WT2, WT3, WT4, WT5))

tpm_bcatMyc_all <- transcript_tpm %>%
  select(-c(WT1, WT2, WT3, WT4, WT5, bcat1, bcat2, bcat3, bcat4, cmyc1, cmyc2, cmyc3, cmyc4, cmyc5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5))

#tpm_Starved_all <- as.matrix(tpm_Starved_all)

tpm_WT_abundant <- transcript_tpm %>%
  dplyr::inner_join(filtered_abundantTx_WT, by = "transcript_id") %>%
  select(-c(bcat1, bcat2, bcat3, bcat4, cmyc1, cmyc2, cmyc3, cmyc4, cmyc5, bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5, WT1.y, WT2.y, WT3.y, WT4.y, WT5.y)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(WT1.x, WT2.x, WT3.x, WT4.x, WT5.x))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(WT1.x, WT2.x, WT3.x, WT4.x, WT5.x)) 

#tpm_Fed_abundant <- as.matrix(tpm_Fed_abundant)

tpm_bcatMyc_abundant <- transcript_tpm %>%
  dplyr::inner_join(filtered_abundantTx_bcatMyc, by = "transcript_id") %>%
  select(-c(WT1, WT2, WT3, WT4, WT5, bcat1, bcat2, bcat3, bcat4, cmyc1, cmyc2, cmyc3, cmyc4, cmyc5, bcat_cmyc1.y, bcat_cmyc2.y, bcat_cmyc3.y, bcat_cmyc4.y, bcat_cmyc5.y)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(bcat_cmyc1.x, bcat_cmyc2.x, bcat_cmyc3.x, bcat_cmyc4.x, bcat_cmyc5.x))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(bcat_cmyc1.x, bcat_cmyc2.x, bcat_cmyc3.x, bcat_cmyc4.x, bcat_cmyc5.x))

#tpm_Starved_abundant <- as.matrix(tpm_Starved_abundant)

library(Biostrings)
library(dplyr)

# remove any rows where transcripts are not being expressed in the condition 

#keep1 <- rowSums(filtered_abundantTx_Fed) > 0
#filtered_abundantTx_Fed <- filtered_abundantTx_Fed[keep1, ]

#keep2 <- rowSums(filtered_abundantTx_Starved) > 0
#filtered_abundantTx_Starved <- filtered_abundantTx_Starved[keep2, ]

# row names to column 

#filtered_abundantTx_Fed <-as.data.frame(abundant_transcripts_Fed) %>%
#  select(-c(mean_tpm, gene_id))

#filtered_abundantTx_Starved <-as.data.frame(abundant_transcripts_Starved) %>%
#  select(-c(mean_tpm, gene_id))

# load the transcriptome (filtered CDS)
# all these CDSs have 5' and 3' UTRs, Start/Stop codon and is divisible by 3 

filtered_CDSfasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat/gencode.vM25.pc_transcripts_CDSs.fa")

# create a data frame with only the expressed transcripts 
LL
WT_CDS_allTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_allTx_WT, by = "transcript_id") %>%
  select(-c(WT1, WT2, WT3, WT4, WT5))

bcatMyc_CDS_allTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_allTx_bcatMyc, by = "transcript_id") %>%
  select(-c(bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5))

WT_CDS_abundantTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_WT, by = "transcript_id") %>%
  select(-c(WT1, WT2, WT3, WT4, WT5))

bcatMyc_CDS_abundantTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_bcatMyc, by = "transcript_id") %>%
  select(-c(bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5))

# sanity check 

#dim(Fed_CDS)
#dim(meanrawcounts_Fed)

#dim(Starved_CDS)
#dim(meanrawcounts_Starved)

# write the newly generated transcriptome dataframes into fasta and save on disk 

library(tidyverse)
library(seqinr)

# remeber to set a new wd 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_all/")

df1 = WT_CDS_allTx
seqs = as.list(dplyr::pull(df1, CDS))
names = dplyr::pull(df1, `transcript_id`)
write.fasta(seqs, names, "WT_CDS_allTx.fasta",
            open = "w", as.string = FALSE)

df2 = bcatMyc_CDS_allTx
seqs = as.list(dplyr::pull(df2, CDS))
names = dplyr::pull(df2, `transcript_id`)
write.fasta(seqs, names, "bcatMyc_CDS_allTx.fasta",
            open = "w", as.string = FALSE)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_abundant/")

df3 = WT_CDS_abundantTx
seqs = as.list(dplyr::pull(df3, CDS))
names = dplyr::pull(df3, `transcript_id`)
write.fasta(seqs, names, "WT_CDS.fasta",
            open = "w", as.string = FALSE)

df4 = bcatMyc_CDS_abundantTx
seqs = as.list(dplyr::pull(df4, CDS))
names = dplyr::pull(df4, `transcript_id`)
write.fasta(seqs, names, "bcatMyc_CDS.fasta",
            open = "w", as.string = FALSE)

#################################################################################################################################

# Calculate the SCU

# double check that the directory has no dodgy fasta files in there otherwise this might mess up the calculations

library("seqinr")

aa <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/codon_aa_table.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)

head(aa)

## WT All Transcripts

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_all/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_all/"

out.file <- ""

file.names <- dir(path, pattern ="WT_CDS_allTx.fasta")

table <- NULL

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
  left_join(tpm_WT_all, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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

write.table(table, file = "WT_allTx_SCUoutput.txt", sep="\t")

##############################################################################
## B-cat/c-myc All Transcripts 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_all/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_all/"

out.file <- ""

file.names <- dir(path, pattern ="bcatMyc_CDS_allTx.fasta")

table <- NULL

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
  left_join(tpm_bcatMyc_all, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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

write.table(table, file = "bcatMyc_allTx_SCUoutput.txt", sep="\t")

###############################################################################################################################

## WT Abundant Transcripts

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_abundant/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_abundant/"

out.file <- ""

file.names <- dir(path, pattern ="WT_CDS.fasta")

table <- NULL

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
  left_join(tpm_WT_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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

write.table(table, file = "WT_abTx_SCUoutput.txt", sep="\t")

##############################################################################
## B-cat/c-myc Abundant Transcripts 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_abundant/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_abundant/"

out.file <- ""

file.names <- dir(path, pattern ="bcatMyc_CDS.fasta")

table <- NULL

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
  left_join(tpm_bcatMyc_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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

write.table(table, file = "bcatMyc_abTx_SCUoutput.txt", sep="\t")

###############################################################################################################################
# combine the generated tables and store them as dataframes 
# remember to reset the wd directory ... all figures generated will be saved there 

### ALL TRANSCRIPTS

library(readr)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_all/")

letter_to_amino <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/letter_to_amino.txt")

WT_allTx_RSCU <- read.table("WT_allTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon")

bcatMyc_allTx_RSCU <- read.table("bcatMyc_allTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon") 

WTvsbcatMyc_allTx_RSCU <- WT_allTx_RSCU %>%
  inner_join(bcatMyc_allTx_RSCU, by = "codon")

WTvsbcatMyc_allTx_RSCU <- WTvsbcatMyc_allTx_RSCU %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(WTvsbcatMyc_allTx_RSCU, file = "WTvsbcatMyc_allTx_SCUoutput.txt", sep="\t")

### ABUNDANT TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_abundant/")

WT_abTx_RSCU <- read.table("WT_abTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon")

bcatMyc_abTx_RSCU <- read.table("bcatMyc_abTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon") 

WTvsbcatMyc_abTx_RSCU <- WT_abTx_RSCU %>%
  inner_join(bcatMyc_abTx_RSCU, by = "codon")

WTvsbcatMyc_abTx_RSCU %>%
  dplyr::rename()

WTvsbcatMyc_abTx_RSCU <- WTvsbcatMyc_abTx_RSCU %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(WTvsbcatMyc_abTx_RSCU, file = "WTvsbcatMyc_allTx_SCUoutput.txt", sep="\t")

###############################################################################################################################

# Generate gingold type plots for both all and abundant transcripts 

### ALL TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_all/")

df_melt_all <- mutate(WTvsbcatMyc_allTx_RSCU,
                      wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                         str_detect(codon, "G$|C$") ~ "GC_wobble"),
                      wobble_single = case_when(str_detect(codon, "A$") ~ "A_wobble",
                                                str_detect(codon, "T$") ~ "T_wobble",
                                                str_detect(codon, "G$") ~ "G_wobble",
                                                str_detect(codon, "C$") ~ "C_wobble"))

# gingold type plot

allTx_gingoldplot <- ggplot(df_melt_all, aes(x=WT_CDS_allTx, y= bcatMyc_CDS_allTx)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (all transcripts)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_allTx_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)


allTx_gingoldplot

allTx_gingoldplot_wobble <- ggplot(df_melt_all, aes(x=WT_CDS_allTx, y= bcatMyc_CDS_allTx)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (all transcripts)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") 

  ggsave(filename = "RSCU_WTvsbcatMyc_allTx_gingold_wobblesingle.png", bg = "white", width = 7, height = 7, dpi = 600)


allTx_gingoldplot_wobble

### ABUNDANT TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_abundant/")

df_melt_abundant <- mutate(WTvsbcatMyc_abTx_RSCU,
                           wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                              str_detect(codon, "G$|C$") ~ "GC_wobble"),
                           wobble_single = case_when(str_detect(codon, "A$") ~ "A_wobble",
                                                     str_detect(codon, "T$") ~ "T_wobble",
                                                     str_detect(codon, "G$") ~ "G_wobble",
                                                     str_detect(codon, "C$") ~ "C_wobble"))

# gingold type plot

abTx_gingoldplot <- ggplot(df_melt_abundant, aes(x=WT_CDS, y= bcatMyc_CDS)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (abundant transcripts only)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_abTx_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)


abTx_gingoldplot

abTx_gingoldplot_wobble <- ggplot(df_melt_abundant, aes(x=WT_CDS, y= bcatMyc_CDS)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (abundant transcripts only)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") 

  ggsave(filename = "RSCU_WTvsbcatMyc_abTx_gingold_wobblesingle.png", bg = "white", width = 7, height = 7, dpi = 600)


abTx_gingoldplot_wobble

###############################################################################################################################
# Sperate up and down regulated genes 
# WT (Ctrl) = downregulated
# b-catenine/c-myc (condition) = upregulated 

# first i need to take the most abundant transcript for all genes and put them in a data frame 
# this will be used to convert the gene_id back to their transcript ID 

##### set the working directory 

WT_gene2abTx <- abundant_transcripts_WT_ab %>%
  select(-c(WT1, WT2, WT3, WT4, WT5))

bcatMyc_gene2abTx <- abundant_transcripts_bcatMyc_ab %>%
  select(-c(bcat_cmyc1, bcat_cmyc2, bcat_cmyc3, bcat_cmyc4, bcat_cmyc5))

upregulated_cutoff0 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC > 0) %>%
  dplyr::inner_join(bcatMyc_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff0.5 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 0.5) %>%
  dplyr::inner_join(bcatMyc_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff1 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 1) %>%
  dplyr::inner_join(bcatMyc_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff2 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 2) %>%
  dplyr::inner_join(bcatMyc_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff0 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC < 0) %>%
  dplyr::inner_join(WT_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff0.5 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -0.5) %>%
  dplyr::inner_join(WT_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff1 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -1) %>%
  dplyr::inner_join(WT_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff2 <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -2) %>%
  dplyr::inner_join(WT_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

###################################

WT_differential_cutoff0 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff0, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

WT_differential_cutoff0.5 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff0.5, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

WT_differential_cutoff1 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff1, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

WT_differential_cutoff2 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff2, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

bcatMyc_differential_cutoff0 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff0, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

bcatMyc_differential_cutoff0.5 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff0.5, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

bcatMyc_differential_cutoff1 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff1, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

bcatMyc_differential_cutoff2 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff2, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

# will also use log2FC cutoffs of 0, 0.5, 1 and 2 

##########################
###############################################################################################################################
# Sperate up and down regulated genes 
# WT (Ctrl) = downregulated
# bcat/c-myc (condition) = upregulated 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_differentials/")

# write the fasta files in the appropriate directory

df5 = WT_differential_cutoff0
seqs = as.list(dplyr::pull(df5, CDS))
names = dplyr::pull(df5, `transcript_id`)
write.fasta(seqs, names, "WT_differential_cutoff0.fasta",
            open = "w", as.string = FALSE)

df6 = WT_differential_cutoff0.5
seqs = as.list(dplyr::pull(df6, CDS))
names = dplyr::pull(df6, `transcript_id`)
write.fasta(seqs, names, "WT_differential_cutoff0.5.fasta",
            open = "w", as.string = FALSE)

df7 = WT_differential_cutoff1
seqs = as.list(dplyr::pull(df7, CDS))
names = dplyr::pull(df7, `transcript_id`)
write.fasta(seqs, names, "WT_differential_cutoff1.fasta",
            open = "w", as.string = FALSE)

df8 = WT_differential_cutoff2
seqs = as.list(dplyr::pull(df8, CDS))
names = dplyr::pull(df8, `transcript_id`)
write.fasta(seqs, names, "WT_differential_cutoff2.fasta",
            open = "w", as.string = FALSE)


df9 = bcatMyc_differential_cutoff0
seqs = as.list(dplyr::pull(df9, CDS))
names = dplyr::pull(df9, `transcript_id`)
write.fasta(seqs, names, "bcatMyc_differential_cutoff0.fasta",
            open = "w", as.string = FALSE)

df10 = bcatMyc_differential_cutoff0.5
seqs = as.list(dplyr::pull(df10, CDS))
names = dplyr::pull(df10, `transcript_id`)
write.fasta(seqs, names, "bcatMyc_differential_cutoff0.5.fasta",
            open = "w", as.string = FALSE)

df11 = bcatMyc_differential_cutoff1
seqs = as.list(dplyr::pull(df11, CDS))
names = dplyr::pull(df11, `transcript_id`)
write.fasta(seqs, names, "bcatMyc_differential_cutoff1.fasta",
            open = "w", as.string = FALSE)

df12 = bcatMyc_differential_cutoff2
seqs = as.list(dplyr::pull(df12, CDS))
names = dplyr::pull(df12, `transcript_id`)
write.fasta(seqs, names, "bcatMyc_differential_cutoff2.fasta",
            open = "w", as.string = FALSE)

########################################################

# calculate the RSCU

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_differentials/")
path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_differentials/"
out.file<-""
file.names <- dir(path, pattern ="bcatMyc_differential_cutoff0.fasta")
table<- NULL
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
  left_join(tpm_bcatMyc_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "bcatMyc_differential_cutoff0_SCUoutput.txt",sep="\t")

########################################################################################
out.file<-""
file.names <- dir(path, pattern ="bcatMyc_differential_cutoff0.5.fasta")
table<- NULL
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
  left_join(tpm_bcatMyc_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "bcatMyc_differential_cutoff0.5_SCUoutput.txt",sep="\t")

###############################################################################

out.file<-""
file.names <- dir(path, pattern ="bcatMyc_differential_cutoff1.fasta")
table<- NULL
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
  left_join(tpm_bcatMyc_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "bcatMyc_differential_cutoff1_SCUoutput.txt",sep="\t")

#################################################################################

out.file<-""
file.names <- dir(path, pattern ="bcatMyc_differential_cutoff2.fasta")
table<- NULL
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
  left_join(tpm_bcatMyc_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "bcatMyc_differential_cutoff2_SCUoutput.txt",sep="\t")

###############################################################################################

out.file<-""
file.names <- dir(path, pattern ="WT_differential_cutoff0.fasta")
table<- NULL
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
  left_join(tpm_WT_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "WT_differential_cutoff0_SCUoutput.txt",sep="\t")

#############################################################################################################

out.file<-""
file.names <- dir(path, pattern ="WT_differential_cutoff0.5.fasta")
table<- NULL
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
  left_join(tpm_WT_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "WT_differential_cutoff0.5_SCUoutput.txt",sep="\t")

#################################################################################################################

out.file<-""
file.names <- dir(path, pattern ="WT_differential_cutoff1.fasta")
table<- NULL
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
  left_join(tpm_WT_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "WT_differential_cutoff1_SCUoutput.txt",sep="\t")

#######################################################################################################

out.file<-""
file.names <- dir(path, pattern ="WT_differential_cutoff2.fasta")
table<- NULL
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
  left_join(tpm_WT_abundant, by = "transcript_id") %>% 
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

x1$Group.1 <- NULL

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
write.table(table, file = "WT_differential_cutoff2_SCUoutput.txt",sep="\t")

########################################################################################################################

# combine all the generated tables 

bcatMyc_cutoff0 <- read.table("bcatMyc_differential_cutoff0_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
bcatMyc_cutoff0.5 <- read.table("bcatMyc_differential_cutoff0.5_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
bcatMyc_cutoff1 <- read.table("bcatMyc_differential_cutoff1_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
bcatMyc_cutoff2 <- read.table("bcatMyc_differential_cutoff2_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")

WT_cutoff0 <- read.table("WT_differential_cutoff0_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
WT_cutoff0.5 <- read.table("WT_differential_cutoff0.5_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
WT_cutoff1 <- read.table("WT_differential_cutoff1_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
WT_cutoff2 <- read.table("WT_differential_cutoff2_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")

WTvsbcatMyc_differential_SCUoutput <- bcatMyc_cutoff0 %>%
  dplyr::inner_join(bcatMyc_cutoff0.5) %>%
  dplyr::inner_join(bcatMyc_cutoff1) %>%
  dplyr::inner_join(bcatMyc_cutoff2) %>%
  dplyr::inner_join(WT_cutoff0) %>%
  dplyr::inner_join(WT_cutoff0.5) %>%
  dplyr::inner_join(WT_cutoff1) %>%
  dplyr::inner_join(WT_cutoff2) 

# Add AA features to the RSCU output

WTvsbcatMyc_differential_SCUoutput <- WTvsbcatMyc_differential_SCUoutput %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(WTvsbcatMyc_differential_SCUoutput, file = "WTvsbcatMyc_differential_SCUoutput.txt", sep="\t")

# plot the gingold type plots 

df_melt_differentials <- mutate(WTvsbcatMyc_differential_SCUoutput,
                                wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                                   str_detect(codon, "G$|C$") ~ "GC_wobble"),
                                wobble_single = case_when(str_detect(codon, "A$") ~ "A_wobble",
                                                          str_detect(codon, "T$") ~ "T_wobble",
                                                          str_detect(codon, "G$") ~ "G_wobble",
                                                          str_detect(codon, "C$") ~ "C_wobble"))


cutoff0_gingoldplot <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff0, y= bcatMyc_differential_cutoff0)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = 0)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff0_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0_gingoldplot

cutoff0.5_gingoldplot <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff0.5, y= bcatMyc_differential_cutoff0.5)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±0.5)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff0.5_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0.5_gingoldplot

cutoff1_gingoldplot <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff1, y= bcatMyc_differential_cutoff1)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±1)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff1_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff1_gingoldplot

cutoff2_gingoldplot <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff2, y= bcatMyc_differential_cutoff2)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±2)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff2_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff2_gingoldplot

###############################################################################################################################
# now do the same but highlight by single wobble

cutoff0_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff0, y= bcatMyc_differential_cutoff0)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = 0)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  #geom_text_repel(aes(label = aminoacid)) +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff0_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0_gingoldplot_wobble

cutoff0.5_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff0.5, y= bcatMyc_differential_cutoff0.5)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±0.5)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff0.5_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0.5_gingoldplot_wobble

cutoff1_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff1, y= bcatMyc_differential_cutoff1)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±1)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff1_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff1_gingoldplot_wobble

cutoff2_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=WT_differential_cutoff2, y= bcatMyc_differential_cutoff2)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log cutoff = ±2)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_WTvsbcatMyc_differential_cutoff2_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff2_gingoldplot_wobble
###############################################################################################################################

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=WT_differential_cutoff0, y= bcatMyc_differential_cutoff0)) +
    geom_point(aes(color = wobble)) +
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
  
  ggsave(filename = paste(this_aminoacid, "RSCU_WTvsbcatmyc_differential_cutoff0_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=WT_differential_cutoff0.5, y= bcatMyc_differential_cutoff0.5)) +
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
  
  ggsave(filename = paste(this_aminoacid, "RSCU_WTvsbcatmyc_differential_cutoff0.5_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}


for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff1_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=WT_differential_cutoff1, y= bcatMyc_differential_cutoff1)) +
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
  
  ggsave(filename = paste(this_aminoacid, "RSCU_WTbcatMyc_differential_cutoff1_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}



for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff2_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=WT_differential_cutoff2, y= bcatMyc_differential_cutoff2)) +
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
  
  ggsave(filename = paste(this_aminoacid, "RSCU_WTvsbcatMyc_differential_cutoff2_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

######################################################################################################################################################
setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/")

# going to reorganise the data at the gene level to get enterez ID / medianTxLength to do GOseq and KEGG analysis 
# using bioMart to get the annotation details and will create a new gene level results table

# view the available databases on biomarts

library(biomaRt)
library(DESeq2)
library(tidyverse)

WTvsbcat_cmyc_GENElevel <- as.data.frame(WTvsbcat_cmyc) %>%
  rownames_to_column("GeneID")

# remove the decimal to convert to just normal geneID rather that having the version no. too 

WTvsbcat_cmyc_GENElevel$GeneID <- sub('\\.[0-9]*$', "", WTvsbcat_cmyc_GENElevel$GeneID)

WTvsbcat_cmyc_GENElevel <- WTvsbcat_cmyc_GENElevel %>%
  column_to_rownames(var = "GeneID")
  
WTvsbcat_cmyc_genelevel_GO <- WTvsbcat_cmyc_GENElevel %>%
  rownames_to_column("GeneID") %>%
  dplyr::select(-c(baseMean, lfcSE)) %>%
  drop_na()

write_tsv(WTvsbcat_cmyc_genelevel_GO, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/WTvsbcat_cmyc_genelevel_GO.txt")

# view the available databases

listMarts(host="uswest.ensembl.org")

## set up connection to ensembl database

ensembl=useMart("ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")

# list the available datasets (species)

listDatasets(ensembl) %>% 
  filter(str_detect(description, "Mouse"))

# specify a data set to use

ensembl = useDataset("mmusculus_gene_ensembl", mart=ensembl)

# check the available "filters" - things you can filter for

listFilters(ensembl) %>% 
  filter(str_detect(name, "ensembl"))

# Set the filter type and values

ourFilterType <- "ensembl_gene_id"
filterValues <- rownames(WTvsbcat_cmyc_GENElevel)

# check the available "attributes" - things you can retreive

listAttributes(ensembl) %>% 
  head(20)

# Set the list of attributes

# you can use this to find a particular attribute of interest e.g. length etc

listAttributes(ensembl) %>%
  filter(str_detect(name, "length"))

#list attributes you want

attributeNames <- c('ensembl_gene_id', 'entrezgene_id')

# run the query

annot <- getBM(attributes=attributeNames, 
               filters = ourFilterType, 
               values = filterValues, 
               mart = ensembl)

head(annot)

# get the median transcript length from the salmon output 

lengthData <- as.data.frame(txi.salmon$length) %>%
  mutate(medianTxlength = rowMeans(txi.salmon$length))

rownames(lengthData) <- rownames(txi.salmon$length)

# remove the unwanted columns

lengthData <- lengthData %>%
  dplyr::select(-c(WT1:bcat_cmyc5)) %>%
  rownames_to_column("GeneID")

# remove the deimal place to get just the Gene ID without version number 

lengthData$GeneID <- sub('\\.[0-9]*$', "", lengthData$GeneID)

 
# check the no. of duplicate GeneIDs you have, that would most probably be because of multiple entrezIDs for each gene 

dpulicate_entrez <- annot %>%  
  add_count(ensembl_gene_id) %>%  
  filter(n>1)

# concatenate entrezIDs

annot <- annot %>%
  group_by(ensembl_gene_id, entrezgene_id) %>%
  summarise(entrezgene_id = paste0(entrezgene_id, collapse = ";")) %>%
  ungroup()


# change column names of the BioMarts annotations

annot <- annot %>%
  rename(GeneID = ensembl_gene_id, Entrez = entrezgene_id)

# sanity checks for any duplicates ... should be 0

annot %>%  
  add_count(GeneID) %>%  
  filter(n>1) %>%
  distinct(GeneID) %>%
  nrow()

# join the DESeq2 data / annotation / median TX length 

Genelevel_annot_WTvsbcatMyc <- WTvsbcat_cmyc_GENElevel %>%
  rownames_to_column("GeneID") %>%
  left_join(annot, "GeneID") %>%
  left_join(lengthData, "GeneID") %>%
  rename(logFC=log2FoldChange, FDR=padj)


  
  
  
library(goseq)

# The input for goseq is a vector that indicates, for each gene, whether or not it is significantly differentially expressed. This should be a named vector, where the names are the gene ids and the values are 1 if the gene is significant and a 0 if it is not.

# I'm using the gene symbol for input since these have been mapped to the transcriptome instead of the Ensembl gene IDs.

supportedOrganisms() %>% filter(str_detect(Genome, "mm"))

signifdata <- as.integer( Genelevel_annot_WTvsbcatMyc$FDR < 0.05 & !is.na(Genelevel_annot_WTvsbcatMyc$FDR) )

names(signifdata) <- Genelevel_annot_WTvsbcatMyc$GeneID

pwf <- nullp(signifdata, "mm10", "ensGene", bias.data = Genelevel_annot_WTvsbcatMyc$medianTxlength)

###################################################################################################

# calculate AA frequency for differential genes

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Null_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_differentials/WT_differential_cutoff2.fasta")

Cre_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/SCU_differentials/bcatMyc_differential_cutoff2.fasta")

# translate the DNA sequences into protein sequence 

Null_diff_translated_fasta <- Biostrings::translate(Null_diff_fasta) 

Cre_diff_translated_fasta <- Biostrings::translate(Cre_diff_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Null_diff_translated_freq <- Biostrings::alphabetFrequency(Null_diff_translated_fasta)

Cre_diff_translated_freq <- Biostrings::alphabetFrequency(Cre_diff_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Null_diff_translated_fasta <- as.data.frame(Null_diff_translated_fasta) %>%
  rownames_to_column("transcript_id")

Null_diff_translated_freq <- as.data.frame(Null_diff_translated_freq)

rownames(Null_diff_translated_freq) <- Null_diff_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Null_diff_translated_freq <- Null_diff_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Null_diff_translated_freq_normalized <- Null_diff_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_WT_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Null_diff_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  na.omit()

Null_AA_freq_normalized_diff <- Null_diff_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Null_AA_freq_normalized_diff <- as.data.frame(Null_AA_freq_normalized_diff) %>%
  rownames_to_column("AA") %>%
  rename(Null_AA_freq = V1) 

###################################################################################################

Cre_diff_translated_fasta <- as.data.frame(Cre_diff_translated_fasta) %>%
  rownames_to_column("transcript_id")

Cre_diff_translated_freq <- as.data.frame(Cre_diff_translated_freq)

rownames(Cre_diff_translated_freq) <- Cre_diff_translated_fasta$transcript_id

Cre_diff_translated_freq <- Cre_diff_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Cre_diff_translated_freq_normalized <- Cre_diff_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_bcatMyc_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Cre_diff_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  na.omit()

Cre_AA_freq_normalized_diff <- Cre_diff_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Cre_AA_freq_normalized_diff <- as.data.frame(Cre_AA_freq_normalized_diff) %>%
  rownames_to_column("AA") %>%
  rename(Cre_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

bcatMyc_AA_freq_normalized_diff <- Null_AA_freq_normalized_diff %>%
  dplyr::inner_join(Cre_AA_freq_normalized_diff, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(bcatMyc_AA_freq_normalized_diff, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/bcat_reanalysis/SCU_WTvsbcatMyc/bcatMyc_AA_freq_normalized_differential_logFC2cutoff.txt")


###################################################################################################

# AA property plots

aa_properties <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/aa_properties.txt", skip=1)

colnames(aa_properties) <- c("aminoacid", "property")

property_plot_input <- bcatMyc_AA_freq_normalized_diff %>%
  left_join(aa_properties, "aminoacid")

ggplot(data = property_plot_input, aes(x=Null_AA_freq, y= Cre_AA_freq))+
  geom_point(aes(color = property)) +
  theme_bw(base_size = 16)+
  labs(title = "Amino Acid Frequency in bcat/Myc Null vs Cre") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  xlab("Null AA frequency")+
  ylab("Cre AA frequency") +
  geom_label_repel(aes(label = aminoacid),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()  

ggsave(filename = "bcatMyc_aa_properties.png", bg = "white", width = 7, height = 7, dpi = 600)
