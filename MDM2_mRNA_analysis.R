setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/")

library(dplyr)
library(readr)

samplefile <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/samplefile_MDM2_maleNullvsCre_gencode.txt")

#samplefile$batch <- factor(samplefile$batch)

library(tximportData)

library(tximportData)

dir <- "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/"

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

png("Normalizedcounts_boxplot_MALE.png")

boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")

# Let's add a blue horizontal line that corresponds to the median logCPM

abline(h=median(assay(vsd)), col="blue")

dev.off()


#################################################################################################################################

# carry out the differential expression with apeglm log fold shrinkage

# had to relevel to make WT the reference for comparison ... normally done alphabetically

ddsTxi_obj$condition <- relevel(ddsTxi_obj$condition, ref = "Null")

ddsTxi <- DESeq(ddsTxi_obj)

res <- results(ddsTxi)

resultsNames(ddsTxi)


# also had to redo the binomial wald test 


ddsTxi <- nbinomWaldTest(ddsTxi)

resultsNames(ddsTxi)

NullvsCre <- lfcShrink(ddsTxi, coef = "condition_Cre_vs_Null", type="apeglm") 


# carry out the differential expression with apeglm log fold shrinkage

#ddsTxi <- DESeq(ddsTxi_obj)

#res <- results(ddsTxi)

#resultsNames(ddsTxi)

#apeglm_ddsTxi <- lfcShrink(ddsTxi, coef="condition_Starved_vs_Fed", type="apeglm")

#NullvsCre <- results(ddsTxi, contCret = c("condition", "Cre", "Null"), alpha = 0.05)

#################################################################################################################################

# p-values and adjusted p-values

# We can order our results table by the smallest p value

resOrdered <- NullvsCre[order(NullvsCre$pvalue),]

# We can summarize some basic tallies using the summary function. 

summary(res)


# How many adjusted p-values were less than 0.05?

sum(NullvsCre$padj < 0.05, na.rm=TRUE)


# MA-plot

png("MA-plot_MALE.png")

plotMA(NullvsCre, ylim=c(-5, 5))

dev.off()

# this gives log2(n + 1)

ntd <- normTransform(ddsTxi_obj)

#################################################################################################################################

library(ggplot2)

# PCA post batch correction 

ddsTxi_obj$batch <- factor(rep(c("1","2","3","4","5","6","7","8","9","10")))

png("PCAplot_pre_batchcorrection_MALE.png")

vsd <- varianceStabilizingTransformation(ddsTxi_obj) # do a variance stablizing transformation on the model matrix 

plotPCA(vsd, "batch") +
  geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5)+
  ylim(-30,15)

# this should give you the variance associated with the effect of the batch 

dev.off()

png("PCAplot_post_batchcorrection_MALE.png")

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch) # applying limma's removeBatchEffect 

plotPCA(vsd, intgroup = c("condition", "batch")) +
  geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5)+
  ylim(-15,15)+
  xlim(-35,35)

# plotting the residuals 

dev.off()

# Replot the PCAs in ggplot 

# need to find way to replot the PCA before and after batch correction 

#################################################################################################################################

# P-VALUES and Adj p-values HISTOGRAM SANITY CHECK 

png("pvalues_histogram_MALE.png")

hist(NullvsCre$pvalue)

dev.off()

png("padj_histogram_MALE.png")

hist(NullvsCre$padj)

dev.off()

#################################################################################################################################
# Extract the DEseq2 Results 

library(dplyr)
library(tidyverse)

NullvsCre_table <- as.data.frame(NullvsCre) %>%
  rownames_to_column("GeneID") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj)

# Exctract columns of interest from annotation

deseq_annotation <- as.data.frame(annotation) %>%
  select(-c(source, type, score, phase, gene_id, level, havana_gene, tag, havana_gene, mgi_id)) %>%
  rownames_to_column("GeneID") 

NullvsCre_table <- NullvsCre_table %>%
  left_join(deseq_annotation, "GeneID")

write_tsv(NullvsCre_table, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/NullvsCre_MALE_DESeq2resultstable.csv")


#################################################################################################################################

# volcano plot (both figure and interactive)

# Enhanced Volcano 

library(ggrepel)
library(EnhancedVolcano)

png("NullvsCre_volcano_MALE.png")

EnhancedVolcano(NullvsCre_table,
                lab = NullvsCre_table$gene_name,
                x = "log2FC",
                y = "FDR",
                ylab = "-Log10(FDR)",
                xlim = c(-8, 8), 
                title = "MDM2 Male Null vs Cre Volcano Plot", 
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

de <- as.integer(NullvsCre$padj <= 0.05 & !is.na(NullvsCre$padj))

normCounts <- log2(counts(ddsTxi_obj))

statusCol <- match(samplefile$condition, c("Null", "Cre"))

glXYPlot(
  x = NullvsCre$log2FoldChange,
  y = -log10(NullvsCre$pvalue),
  xlab = "log2FC",
  ylab = "-log10(padj)",
  side.xlab = "Condition",
  side.ylab = "Log(CPM)",
  main = "Male MDM2 Null vs Cre",
  counts = normCounts,
  groups = group,
  sample.cols = statusCol,
  status = de,
  anno = NullvsCre_table[, c("GeneID", "gene_name")],
  folder = "interactive_volcano")

#################################################################################################################################

# HEATMAPS

library(ComplexHeatmap)
library(circlize)

# get the top 200 genes
sigGenes <- as.data.frame(NullvsCre_table) %>% 
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

png("stringent_heatmap_MALE.png")

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

abundant_transcripts_Null <- transcript_tpm %>%
  select(-c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5))) %>%
  filter(mean_tpm != 0)

abundant_transcripts_Cre <- transcript_tpm %>%
  select(-c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5))) %>%
  filter(mean_tpm != 0)

# add the gene_id from the previously constructed tx 2 gene list

tx2gene_df <- as.data.frame(tx2gene)

abundant_transcripts_Null <- dplyr::inner_join(abundant_transcripts_Null, tx2gene_df, by = "transcript_id")

abundant_transcripts_Cre <- dplyr::inner_join(abundant_transcripts_Cre, tx2gene_df, by = "transcript_id")

# select transcripts with the biggest mean tpm ... selecting the most abundant transcript in each gene set

abundant_transcripts_Null_all <- abundant_transcripts_Null %>%
  group_by(gene_id) %>%
  #top_n(n = 1, wt = mean_tpm) %>%
  filter(mean_tpm != 0)

abundant_transcripts_Cre_all <- abundant_transcripts_Cre %>%
  group_by(gene_id) %>%
  #top_n(n = 1, wt = mean_tpm) %>%
  filter(mean_tpm != 0)

abundant_transcripts_Null_ab <- abundant_transcripts_Null %>%
  group_by(gene_id) %>%
  top_n(n = 1, wt = mean_tpm) 

abundant_transcripts_Cre_ab <- abundant_transcripts_Cre %>%
  group_by(gene_id) %>%
  top_n(n = 1, wt = mean_tpm) 

# dont delete this until i know i dont need it !!!!!
#transcripts_byTPM <- abundant_transcripts %>%
#select(c(gene_id, transcript_id))


# sanity check for duplications ... change according to condition

duplications <- abundant_transcripts_Null_all%>%
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

filtered_abundantTx_Null <- as.data.frame(abundant_transcripts_Null_ab) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_abundantTx_Cre <- as.data.frame(abundant_transcripts_Cre_ab) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_allTx_Null <- as.data.frame(abundant_transcripts_Null_all) %>%
  select(-c(mean_tpm, gene_id)) 

filtered_allTx_Cre <- as.data.frame(abundant_transcripts_Cre_all) %>%
  select(-c(mean_tpm, gene_id)) 

# create matrix with transcript tpms

tpm_Null_all <- transcript_tpm %>%
  select(-c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5))

#tpm_Fed_all <- as.matrix(tpm_Fed_all)

tpm_Cre_all <- transcript_tpm %>%
  select(-c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5))

#tpm_Starved_all <- as.matrix(tpm_Starved_all)

tpm_Null_abundant <- transcript_tpm %>%
  dplyr::inner_join(filtered_abundantTx_Null, by = "transcript_id") %>%
  select(-c(AAV_Null1.y, AAV_Null2.y, AAV_Null3.y, AAV_Null4.y, AAV_Null5.y, AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(AAV_Null1.x, AAV_Null2.x, AAV_Null3.x, AAV_Null4.x, AAV_Null5.x))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(AAV_Null1.x, AAV_Null2.x, AAV_Null3.x, AAV_Null4.x, AAV_Null5.x)) 

#tpm_Fed_abundant <- as.matrix(tpm_Fed_abundant)

tpm_Cre_abundant <- transcript_tpm %>%
  dplyr::inner_join(filtered_abundantTx_Cre, by = "transcript_id") %>%
  select(-c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5, AAV_Cre1.y, AAV_Cre2.y, AAV_Cre3.y, AAV_Cre4.y, AAV_Cre5.y)) %>%
  rowwise() %>%
  mutate(mean_tpm = mean(c(AAV_Cre1.x, AAV_Cre2.x, AAV_Cre3.x, AAV_Cre4.x, AAV_Cre5.x))) %>%
  filter(mean_tpm != 0) %>%
  select(-c(AAV_Cre1.x, AAV_Cre2.x, AAV_Cre3.x, AAV_Cre4.x, AAV_Cre5.x))

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

filtered_CDSfasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/gencode.vM25.pc_transcripts_CDSs.fa")

# create a data frame with only the expressed transcripts 

Null_CDS_allTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_allTx_Null, by = "transcript_id") %>%
  select(-c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5))

Cre_CDS_allTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_allTx_Cre, by = "transcript_id") %>%
  select(-c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5))

Null_CDS_abundantTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_Null, by = "transcript_id") %>%
  select(-c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5))

Cre_CDS_abundantTx <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_Cre, by = "transcript_id") %>%
  select(-c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5))

# sanity check 

#dim(Fed_CDS)
#dim(meanrawcounts_Fed)

#dim(Starved_CDS)
#dim(meanrawcounts_Starved)

# write the newly generated transcriptome dataframes into fasta and save on disk 

library(tidyverse)
library(seqinr)

# remeber to set a new wd 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_all/")

df1 = Null_CDS_allTx
seqs = as.list(dplyr::pull(df1, CDS))
names = dplyr::pull(df1, `transcript_id`)
write.fasta(seqs, names, "Null_CDS_allTx.fasta",
            open = "w", as.string = FALSE)

df2 = Cre_CDS_allTx
seqs = as.list(dplyr::pull(df2, CDS))
names = dplyr::pull(df2, `transcript_id`)
write.fasta(seqs, names, "Cre_CDS_allTx.fasta",
            open = "w", as.string = FALSE)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_abundant/")

df3 = Null_CDS_abundantTx
seqs = as.list(dplyr::pull(df3, CDS))
names = dplyr::pull(df3, `transcript_id`)
write.fasta(seqs, names, "Null_CDS.fasta",
            open = "w", as.string = FALSE)

df4 = Cre_CDS_abundantTx
seqs = as.list(dplyr::pull(df4, CDS))
names = dplyr::pull(df4, `transcript_id`)
write.fasta(seqs, names, "Cre_CDS.fasta",
            open = "w", as.string = FALSE)

#################################################################################################################################

# Calculate the SCU

# double check that the directory has no dodgy fasta files in there otherwise this might mess up the calculations

library("seqinr")

aa <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/codon_aa_table.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)

head(aa)

## Null All Transcripts

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_all/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_all/"

out.file <- ""

file.names <- dir(path, pattern ="Null_CDS_allTx.fasta")

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
  left_join(tpm_Null_all, by = "transcript_id") %>% 
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

write.table(table, file = "Null_allTx_SCUoutput.txt", sep="\t")

##############################################################################
## Cre All Transcripts 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_all/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_all/"

out.file <- ""

file.names <- dir(path, pattern ="Cre_CDS_allTx.fasta")

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
  left_join(tpm_Cre_all, by = "transcript_id") %>% 
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

write.table(table, file = "Cre_allTx_SCUoutput.txt", sep="\t")

###############################################################################################################################

## Null Abundant Transcripts

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_abundant/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_abundant/"

out.file <- ""

file.names <- dir(path, pattern ="Null_CDS.fasta")

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
  left_join(tpm_Null_abundant, by = "transcript_id") %>% 
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

write.table(table, file = "Null_abTx_SCUoutput.txt", sep="\t")

##############################################################################
## Cre Abundant Transcripts 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_abundant/")

path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_abundant/"

out.file <- ""

file.names <- dir(path, pattern ="Cre_CDS.fasta")

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
  left_join(tpm_Cre_abundant, by = "transcript_id") %>% 
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

write.table(table, file = "Cre_abTx_SCUoutput.txt", sep="\t")

###############################################################################################################################
# combine the generated tables and store them as dataframes 
# remember to reset the wd directory ... all figures generated will be saved there 

### ALL TRANSCRIPTS

library(readr)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_all/")

letter_to_amino <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/letter_to_amino.txt")

Null_allTx_RSCU <- read.table("Null_allTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon")

Cre_allTx_RSCU <- read.table("Cre_allTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon") 

NullvsCre_allTx_RSCU <- Null_allTx_RSCU %>%
  inner_join(Cre_allTx_RSCU, by = "codon")

NullvsCre_allTx_RSCU <- NullvsCre_allTx_RSCU %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(NullvsCre_allTx_RSCU, file = "NullvsCre_allTx_SCUoutput.txt", sep="\t")

### ABUNDANT TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_abundant/")

Null_abTx_RSCU <- read.table("Null_abTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon")

Cre_abTx_RSCU <- read.table("Cre_abTx_SCUoutput.txt") %>%
  as.data.frame() %>%
  rownames_to_column(var = "codon") 

NullvsCre_abTx_RSCU <- Null_abTx_RSCU %>%
  inner_join(Cre_abTx_RSCU, by = "codon")

NullvsCre_abTx_RSCU <- NullvsCre_abTx_RSCU %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(NullvsCre_abTx_RSCU, file = "NullvsCre_allTx_SCUoutput.txt", sep="\t")

###############################################################################################################################

# Generate gingold type plots for both all and abundant transcripts 

### ALL TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_all/")

df_melt_all <- mutate(NullvsCre_allTx_RSCU,
                      wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                         str_detect(codon, "G$|C$") ~ "GC_wobble"),
                      wobble_single = case_when(str_detect(codon, "A$") ~ "A_wobble",
                                                str_detect(codon, "T$") ~ "T_wobble",
                                                str_detect(codon, "G$") ~ "G_wobble",
                                                str_detect(codon, "C$") ~ "C_wobble"))

# gingold type plot

allTx_gingoldplot <- ggplot(df_melt_all, aes(x=Null_CDS_allTx, y= Cre_CDS_allTx)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (all transcripts)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_allTx_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)


allTx_gingoldplot


allTx_gingoldplot_wobble <- ggplot(df_melt_all, aes(x=Null_CDS_allTx, y= Cre_CDS_allTx)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (all transcripts)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") 

  ggsave(filename = "RSCU_NullvsCre_allTx_gingold_wobblesingle.png", bg = "white", width = 7, height = 7, dpi = 600)


allTx_gingoldplot_wobble

### ABUNDANT TRANSCRIPTS

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_abundant/")

df_melt_abundant <- mutate(NullvsCre_abTx_RSCU,
                           wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                              str_detect(codon, "G$|C$") ~ "GC_wobble"),
                           wobble_single = case_when(str_detect(codon, "A$") ~ "A_wobble",
                                                     str_detect(codon, "T$") ~ "T_wobble",
                                                     str_detect(codon, "G$") ~ "G_wobble",
                                                     str_detect(codon, "C$") ~ "C_wobble"))

# gingold type plot

abTx_gingoldplot <- ggplot(df_melt_abundant, aes(x=Null_CDS, y= Cre_CDS)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (abundant transcripts only)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_abTx_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)


abTx_gingoldplot

abTx_gingoldplot_single <- ggplot(df_melt_abundant, aes(x=Null_CDS, y= Cre_CDS)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "RSCU weighted by TPM (abundant transcripts only)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") 

  ggsave(filename = "RSCU_NullvsCre_abTx_gingold_single.png", bg = "white", width = 7, height = 7, dpi = 600)


abTx_gingoldplot_single

###############################################################################################################################
# Sperate up and down regulated genes 
# Null (Ctrl) = downregulated
# Cre (condition) = upregulated 

# first i need to take the most abundant transcript for all genes and put them in a data frame 
# this will be used to convert the gene_id back to their transcript ID 

##### set the working directory 

Null_gene2abTx <- abundant_transcripts_Null_ab %>%
  select(-c(AAV_Null1, AAV_Null2, AAV_Null3, AAV_Null4, AAV_Null5))

Cre_gene2abTx <- abundant_transcripts_Cre_ab %>%
  select(-c(AAV_Cre1, AAV_Cre2, AAV_Cre3, AAV_Cre4, AAV_Cre5))

upregulated_cutoff0 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC > 0) %>%
  dplyr::inner_join(Cre_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff0.5 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 0.5) %>%
  dplyr::inner_join(Cre_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff1 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 1) %>%
  dplyr::inner_join(Cre_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

upregulated_cutoff2 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC >= 2) %>%
  dplyr::inner_join(Cre_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff0 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC < 0) %>%
  dplyr::inner_join(Null_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff0.5 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -0.5) %>%
  dplyr::inner_join(Null_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff1 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -1) %>%
  dplyr::inner_join(Null_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

downregulated_cutoff2 <- as.data.frame(NullvsCre) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC <= -2) %>%
  dplyr::inner_join(Null_gene2abTx, by = "gene_id") %>%
  select(-c(gene_id))

###################################

Null_differential_cutoff0 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff0, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Null_differential_cutoff0.5 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff0.5, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Null_differential_cutoff1 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff1, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Null_differential_cutoff2 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated_cutoff2, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Cre_differential_cutoff0 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff0, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Cre_differential_cutoff0.5 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff0.5, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Cre_differential_cutoff1 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff1, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

Cre_differential_cutoff2 <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated_cutoff2, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, pvalue, FDR))

# will also use log2FC cutoffs of 0, 0.5, 1 and 2 

##########################
###############################################################################################################################
# Sperate up and down regulated genes 
# Null (Ctrl) = downregulated
# Cre (condition) = upregulated 

# first i need to take the most abundant transcript for all genes and put them in a data frame 
# this will be used to convert the gene_id back to their transcript ID 

# find the most abundant transcript (but for all conditions instead of for each condition)

#abundant_transcripts_everything <- transcript_tpm %>%
  #rowwise() %>%
  #mutate(mean_tpm = mean(c(Resec1, Resec2, Resec3, Resec4, Resec5, Regen1, Regen2, Regen3, Regen4, Regen5))) %>%
  #filter(mean_tpm != 0)

# get the gene id for each transcript 

#abundant_transcripts_everything <- dplyr::inner_join(abundant_transcripts_everything, tx2gene_df, by = "transcript_id")

# take the most abundant transcript for each gene 

#abundant_transcripts_everything <- abundant_transcripts_everything %>%
  #group_by(gene_id) %>%
  #top_n(n = 1, wt = mean_tpm) %>%
  #filter(mean_tpm != 0) %>%
  #select(-c(Resec1, Resec2, Resec3, Resec4, Resec5, Regen1, Regen2, Regen3, Regen4, Regen5)) %>%
  #ungroup() %>%
  #select(-c(gene_id))

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_differentials/")

# write the fasta files in the appropriate directory

df5 = Null_differential_cutoff0
seqs = as.list(dplyr::pull(df5, CDS))
names = dplyr::pull(df5, `transcript_id`)
write.fasta(seqs, names, "Null_differential_cutoff0.fasta",
            open = "w", as.string = FALSE)

df6 = Null_differential_cutoff0.5
seqs = as.list(dplyr::pull(df6, CDS))
names = dplyr::pull(df6, `transcript_id`)
write.fasta(seqs, names, "Null_differential_cutoff0.5.fasta",
            open = "w", as.string = FALSE)

df7 = Null_differential_cutoff1
seqs = as.list(dplyr::pull(df7, CDS))
names = dplyr::pull(df7, `transcript_id`)
write.fasta(seqs, names, "Null_differential_cutoff1.fasta",
            open = "w", as.string = FALSE)

df8 = Null_differential_cutoff2
seqs = as.list(dplyr::pull(df8, CDS))
names = dplyr::pull(df8, `transcript_id`)
write.fasta(seqs, names, "Null_differential_cutoff2.fasta",
            open = "w", as.string = FALSE)


df9 = Cre_differential_cutoff0
seqs = as.list(dplyr::pull(df9, CDS))
names = dplyr::pull(df9, `transcript_id`)
write.fasta(seqs, names, "Cre_differential_cutoff0.fasta",
            open = "w", as.string = FALSE)

df10 = Cre_differential_cutoff0.5
seqs = as.list(dplyr::pull(df10, CDS))
names = dplyr::pull(df10, `transcript_id`)
write.fasta(seqs, names, "Cre_differential_cutoff0.5.fasta",
            open = "w", as.string = FALSE)

df11 = Cre_differential_cutoff1
seqs = as.list(dplyr::pull(df11, CDS))
names = dplyr::pull(df11, `transcript_id`)
write.fasta(seqs, names, "Cre_differential_cutoff1.fasta",
            open = "w", as.string = FALSE)

df12 = Cre_differential_cutoff2
seqs = as.list(dplyr::pull(df12, CDS))
names = dplyr::pull(df12, `transcript_id`)
write.fasta(seqs, names, "Cre_differential_cutoff2.fasta",
            open = "w", as.string = FALSE)

########################################################

# calculate the RSCU

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_differentials/")
path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_differentials/"
out.file<-""
file.names <- dir(path, pattern ="Cre_differential_cutoff0.fasta")
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
  left_join(tpm_Cre_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Cre_differential_cutoff0_SCUoutput.txt",sep="\t")

########################################################################################
out.file<-""
file.names <- dir(path, pattern ="Cre_differential_cutoff0.5.fasta")
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
  left_join(tpm_Cre_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Cre_differential_cutoff0.5_SCUoutput.txt",sep="\t")

###############################################################################

out.file<-""
file.names <- dir(path, pattern ="Cre_differential_cutoff1.fasta")
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
  left_join(tpm_Cre_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Cre_differential_cutoff1_SCUoutput.txt",sep="\t")

#################################################################################

out.file<-""
file.names <- dir(path, pattern ="Cre_differential_cutoff2.fasta")
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
  left_join(tpm_Cre_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Cre_differential_cutoff2_SCUoutput.txt",sep="\t")

###############################################################################################

out.file<-""
file.names <- dir(path, pattern ="Null_differential_cutoff0.fasta")
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
  left_join(tpm_Null_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Null_differential_cutoff0_SCUoutput.txt",sep="\t")

#############################################################################################################

out.file<-""
file.names <- dir(path, pattern ="Null_differential_cutoff0.5.fasta")
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
  left_join(tpm_Null_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Null_differential_cutoff0.5_SCUoutput.txt",sep="\t")

#################################################################################################################

out.file<-""
file.names <- dir(path, pattern ="Null_differential_cutoff1.fasta")
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
  left_join(tpm_Null_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Null_differential_cutoff1_SCUoutput.txt",sep="\t")

#######################################################################################################

out.file<-""
file.names <- dir(path, pattern ="Null_differential_cutoff2.fasta")
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
  left_join(tpm_Null_abundant, by = "transcript_id") %>% 
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
write.table(table, file = "Null_differential_cutoff2_SCUoutput.txt",sep="\t")

########################################################################################################################

# combine all the generated tables 

Cre_cutoff0 <- read.table("Cre_differential_cutoff0_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Cre_cutoff0.5 <- read.table("Cre_differential_cutoff0.5_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Cre_cutoff1 <- read.table("Cre_differential_cutoff1_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Cre_cutoff2 <- read.table("Cre_differential_cutoff2_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")

Null_cutoff0 <- read.table("Null_differential_cutoff0_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Null_cutoff0.5 <- read.table("Null_differential_cutoff0.5_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Null_cutoff1 <- read.table("Null_differential_cutoff1_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")
Null_cutoff2 <- read.table("Null_differential_cutoff2_SCUoutput.txt") %>%
  rownames_to_column(var = "codon")

NullvsCre_differential_SCUoutput <- Cre_cutoff0 %>%
  dplyr::inner_join(Cre_cutoff0.5) %>%
  dplyr::inner_join(Cre_cutoff1) %>%
  dplyr::inner_join(Cre_cutoff2) %>%
  dplyr::inner_join(Null_cutoff0) %>%
  dplyr::inner_join(Null_cutoff0.5) %>%
  dplyr::inner_join(Null_cutoff1) %>%
  dplyr::inner_join(Null_cutoff2) 

# Add AA features to the RSCU output

NullvsCre_differential_SCUoutput <- NullvsCre_differential_SCUoutput %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(NullvsCre_differential_SCUoutput, file = "NullvsCre_differential_SCUoutput.txt", sep="\t")

# plot the gingold type plots 

df_melt_differentials <- mutate(NullvsCre_differential_SCUoutput,
                                wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                                   str_detect(codon, "G$|C$") ~ "GC_wobble"),
                                wobble_single = case_when(str_detect(codon, "A$") ~ "A_wobble",
                                                          str_detect(codon, "T$") ~ "T_wobble",
                                                          str_detect(codon, "G$") ~ "G_wobble",
                                                          str_detect(codon, "C$") ~ "C_wobble"))


cutoff0_gingoldplot <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff0, y= Cre_differential_cutoff0)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = 0)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff0_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0_gingoldplot

cutoff0.5_gingoldplot <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff0.5, y= Cre_differential_cutoff0.5)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = ±0.5)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff0.5_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0.5_gingoldplot

cutoff1_gingoldplot <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff1, y= Cre_differential_cutoff1)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = ±1)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff1_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff1_gingoldplot

cutoff2_gingoldplot <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff2, y= Cre_differential_cutoff2)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = ±2)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff2_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff2_gingoldplot

###############################################################################################################################
# now do the same but highlight by single wobble

cutoff0_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff0, y= Cre_differential_cutoff0)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = 0)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  #geom_text_repel(aes(label = aminoacid)) +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff0_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0_gingoldplot_wobble

cutoff0.5_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff0.5, y= Cre_differential_cutoff0.5)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = ±0.5)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff0.5_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff0.5_gingoldplot_wobble

cutoff1_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff1, y= Cre_differential_cutoff1)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = ±1)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff1_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff1_gingoldplot_wobble

cutoff2_gingoldplot_wobble <- ggplot(df_melt_differentials, aes(x=Null_differential_cutoff2, y= Cre_differential_cutoff2)) +
  geom_point(aes(color = wobble_single)) +
  theme_bw(base_size = 16)+
  labs(title = "mRNA RSCU of Differential Genes (Log2FC cutoff = ±2)") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  ggsave(filename = "RSCU_NullvsCre_differential_cutoff2_wobble_gingold.png", bg = "white", width = 7, height = 7, dpi = 600)

cutoff2_gingoldplot_wobble
###############################################################################################################################

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff0.5, y=Cre_differential_cutoff0.5)) +
    geom_point(aes(color = wobble_single)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log = ±0.5)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff0.5_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff0.5_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff0, y=Cre_differential_cutoff0)) +
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
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff0_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}


for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff1_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff1, y=Cre_differential_cutoff1)) +
    geom_point(aes(color = wobble_single)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log = ±1)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff1_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff1_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}



for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff2_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff2, y=Cre_differential_cutoff2)) +
    geom_point(aes(color = wobble_single)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log = ±2)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff2_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff2_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

###############################################################################################################################

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0.5_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff0.5, y=Cre_differential_cutoff0.5)) +
    geom_point(aes(color = wobble)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log = ±0.5)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff0.5_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff0.5_AT_GC_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff0_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff0, y=Cre_differential_cutoff0)) +
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
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff0_AT_GC_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}


for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff1_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff1, y=Cre_differential_cutoff1)) +
    geom_point(aes(color = wobble)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log = ±1)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff1_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff1_AT_GC_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}



for (this_aminoacid in unique(df_melt_differentials$aminoacid)) {
  
  cutoff2_gingold_loop <- ggplot(
    data = subset(df_melt_differentials, aminoacid == this_aminoacid),
    aes(x=Null_differential_cutoff2, y=Cre_differential_cutoff2)) +
    geom_point(aes(color = wobble)) +
    theme_bw(base_size = 16)+
    labs(title = paste(this_aminoacid , "RSCU of Differential Genes (Log = ±2)")) +
    geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
    ylim(0,1) +
    xlim(0,1)
  
  cutoff2_gingold_loop +
    geom_label_repel(aes(label = codon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +
    theme_classic()
  
  ggsave(filename = paste(this_aminoacid, "RSCU_maleNullvsCre_differential_cutoff2_AT_GC_gingold.png", sep = "_"), bg = "white", width = 7, height = 7, dpi = 600)
}

###############################################################################################################################

# calculate AA freq for each gene and normalise to their abundance (TPM)

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Null_allTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/codon_freq/Null_CDS_allTx.fasta")

Cre_allTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/codon_freq/Cre_CDS_allTx.fasta")

# translate the DNA sequences into protein sequence 

Null_allTx_translated_fasta <- Biostrings::translate(Null_allTx_fasta) 

Cre_allTx_translated_fasta <- Biostrings::translate(Cre_allTx_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Null_allTx_translated_freq <- Biostrings::alphabetFrequency(Null_allTx_translated_fasta)

Cre_allTx_translated_freq <- Biostrings::alphabetFrequency(Cre_allTx_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Null_allTx_translated_fasta <- as.data.frame(Null_allTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Null_allTx_translated_freq <- as.data.frame(Null_allTx_translated_freq)

rownames(Null_allTx_translated_freq) <- Null_allTx_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Null_allTx_translated_freq <- Null_allTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Null_allTx_translated_freq_normalized <- Null_allTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Null_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Null_allTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Null_AA_freq_normalized <- Null_allTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Null_AA_freq_normalized <- as.data.frame(Null_AA_freq_normalized) %>%
  rownames_to_column("AA") %>%
  rename(Null_AA_freq = V1) 

###################################################################################################

Cre_allTx_translated_fasta <- as.data.frame(Cre_allTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Cre_allTx_translated_freq <- as.data.frame(Cre_allTx_translated_freq)

rownames(Cre_allTx_translated_freq) <- Cre_allTx_translated_fasta$transcript_id

Cre_allTx_translated_freq <- Cre_allTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Cre_allTx_translated_freq_normalized <- Cre_allTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Cre_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Cre_allTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Cre_AA_freq_normalized <- Cre_allTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Cre_AA_freq_normalized <- as.data.frame(Cre_AA_freq_normalized) %>%
  rownames_to_column("AA") %>%
  rename(Cre_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

MDM2_AA_freq_normalized <- Null_AA_freq_normalized %>%
  dplyr::inner_join(Cre_AA_freq_normalized, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(MDM2_AA_freq_normalized, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/codon_freq/MDM2_AA_freq_normalized.txt")

###############################################################################################################################

# calculate AA freq for the most abundant transcript and normalise to their abundance (TPM)

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Null_abTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/codon_freq/Null_CDS.fasta")

Cre_abTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/codon_freq/Cre_CDS.fasta")

# translate the DNA sequences into protein sequence 

Null_abTx_translated_fasta <- Biostrings::translate(Null_abTx_fasta) 

Cre_abTx_translated_fasta <- Biostrings::translate(Cre_abTx_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Null_abTx_translated_freq <- Biostrings::alphabetFrequency(Null_abTx_translated_fasta)

Cre_abTx_translated_freq <- Biostrings::alphabetFrequency(Cre_abTx_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Null_abTx_translated_fasta <- as.data.frame(Null_abTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Null_abTx_translated_freq <- as.data.frame(Null_abTx_translated_freq)

rownames(Null_abTx_translated_freq) <- Null_abTx_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Null_abTx_translated_freq <- Null_abTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Null_abTx_translated_freq_normalized <- Null_abTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Null_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Null_abTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Null_AA_freq_normalized_ab <- Null_abTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Null_AA_freq_normalized_ab <- as.data.frame(Null_AA_freq_normalized_ab) %>%
  rownames_to_column("AA") %>%
  rename(Null_AA_freq = V1) 

###################################################################################################

Cre_abTx_translated_fasta <- as.data.frame(Cre_abTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Cre_abTx_translated_freq <- as.data.frame(Cre_abTx_translated_freq)

rownames(Cre_abTx_translated_freq) <- Cre_abTx_translated_fasta$transcript_id

Cre_abTx_translated_freq <- Cre_abTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Cre_abTx_translated_freq_normalized <- Cre_abTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Cre_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Cre_abTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Cre_AA_freq_normalized_ab <- Cre_abTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Cre_AA_freq_normalized_ab <- as.data.frame(Cre_AA_freq_normalized_ab) %>%
  rownames_to_column("AA") %>%
  rename(Cre_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

MDM2_AA_freq_normalized_ab <- Null_AA_freq_normalized_ab %>%
  dplyr::inner_join(Cre_AA_freq_normalized_ab, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(MDM2_AA_freq_normalized_ab, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/codon_freq/MDM2_AA_freq_normalized_ab.txt")

###################################################################################################
###################################################################################################

# Radar Plots to visualise AA frequency 

library(fmsb)

Radar_allTx <- MDM2_AA_freq_normalized %>%
  column_to_rownames(var = "aminoacid") %>%
  t(.)

# convert to percentages

Radar_allTx_input <- as.data.frame(Radar_allTx) %>%
  mutate_if(is.numeric, ~ . * 100)

rownames(Radar_allTx_input) <- row.names(Radar_allTx)

# add min and max rows ... here only using 15% since abundance doesnt go over 15% for any AA and makes plot easier to visualise

Radar_allTx_input <- rbind(rep(15, 21), rep(0, 21), Radar_allTx_input)

# Customize and plot the Radar plots

# Color vector

colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

png(filename ="Radar_allTx_plot.png")

Radar_allTx_plot <- radarchart(as.data.frame(Radar_allTx_input),
           axistype = 2,
           seg = 3,
           #custom polygon
           pcol=colors_border, 
           pfcol=colors_in, 
           plwd=4, 
           plty=1,
           #custom the grid
           cglcol="grey", 
           cglty=1, 
           axislabcol="black", 
           caxislabels=seq(0,20,5), 
           cglwd=0.8,
           #custom labels
           vlcex=1.4,
           title = "Amino acid abundance of all transcripts in MDM2 Null vs Cre",
           paxislabels =  "15%",
           palcex = 0.9)

legend(x=1.3, y=1, legend = rownames(Radar_allTx_input[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.2, pt.cex=3)

dev.off()

###################################################################################################

Radar_abTx <- MDM2_AA_freq_normalized_ab %>%
  column_to_rownames(var = "aminoacid") %>%
  t(.)

# convert to percentages

Radar_abTx_input <- as.data.frame(Radar_abTx) %>%
  mutate_if(is.numeric, ~ . * 100)

rownames(Radar_abTx_input) <- row.names(Radar_abTx)

# add min and max rows ... here only using 15% since abundance doesnt go over 15% for any AA and makes plot easier to visualise

Radar_abTx_input <- rbind(rep(15, 21), rep(0, 21), Radar_abTx_input)

png(filename ="Radar_abTx_plot.png")

Radar_abTx_plot <- radarchart(as.data.frame(Radar_abTx_input),
                               axistype = 2,
                               seg = 3,
                               #custom polygon
                               pcol=colors_border, 
                               pfcol=colors_in, 
                               plwd=4, 
                               plty=1,
                               #custom the grid
                               cglcol="grey", 
                               cglty=1, 
                               axislabcol="black", 
                               caxislabels=seq(0,20,5), 
                               cglwd=0.8,
                               #custom labels
                               vlcex=1.4,
                               title = "Amino acid abundance of most abundant transcript per gene in MDM2 Null vs Cre",
                               paxislabels =  "15%",
                               palcex = 0.9)

legend(x=1.3, y=1, legend = rownames(Radar_allTx_input[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.2, pt.cex=3)

dev.off()

###################################################################################################

# calculate AA frequency for differential genes

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Null_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_differentials/Null_differential_cutoff2.fasta")

Cre_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/SCU_differentials/Cre_differential_cutoff2.fasta")

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
  left_join(tpm_Null_abundant, by = "transcript_id") %>% 
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
  left_join(tpm_Cre_abundant, by = "transcript_id") %>% 
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

MDM2_AA_freq_normalized_diff <- Null_AA_freq_normalized_diff %>%
  dplyr::inner_join(Cre_AA_freq_normalized_diff, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(MDM2_AA_freq_normalized_diff, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/MDM2/SCU_male/codon_freq/MDM2_AA_freq_normalized_differential_logFC2cutoff.txt")


###################################################################################################

# AA property plots

aa_properties <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/aa_properties.txt", skip=1)

colnames(aa_properties) <- c("aminoacid", "property")

property_plot_input <- MDM2_AA_freq_normalized_diff %>%
  left_join(aa_properties, "aminoacid")

ggplot(data = property_plot_input, aes(x=Null_AA_freq, y= Cre_AA_freq))+
  geom_point(aes(color = property)) +
  theme_bw(base_size = 16)+
  labs(title = "Amino Acid Frequency in MDM2 Null vs Cre") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  xlab("Null AA frequency")+
  ylab("Cre AA frequency") +
  geom_label_repel(aes(label = aminoacid),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()  

ggsave(filename = "MDM2_aa_properties.png", bg = "white", width = 7, height = 7, dpi = 600)

