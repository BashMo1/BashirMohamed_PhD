setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/")

library(dplyr)
library(readr)

samplefile <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/samplefile_FvS_gencode.txt")

library(tximportData)

dir <- "/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding"

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

# carry out the differential expression 

ddsTxi <- DESeq(ddsTxi_obj)

res <- results(ddsTxi)

resultsNames(ddsTxi)

#apeglm_ddsTxi <- lfcShrink(ddsTxi, coef="condition_Starved_vs_Fed", type="apeglm")

FedvsStarved <- results(ddsTxi, contStarvedt = c("condition", "Starved", "Fed"), alpha = 0.05)

#################################################################################################################################

# p-values and adjusted p-values

# We can order our results table by the smallest p value

resOrdered <- FedvsStarved[order(FedvsStarved$pvalue),]

# We can summarize some basic tallies using the summary function. 

summary(res)


# How many adjusted p-values were less than 0.05?

sum(FedvsStarved$padj < 0.05, na.rm=TRUE)


# MA-plot

png("MA-plot.png")

plotMA(FedvsStarved, ylim=c(-5, 5))

dev.off()

# this gives log2(n + 1)

ntd <- normTransform(ddsTxi_obj)

#################################################################################################################################

library(ggplot2)

# PCA post batch correction 

ddsTxi_obj$batch <- factor(rep(c("Batch1", "Batch2", "Batch3")))

png("PCAplot_pre_batchcorrection.png")

vsd <- varianceStabilizingTransformation(ddsTxi_obj) # do a variance stablizing transformation on the model matrix 

plotPCA(vsd, "batch") +
  geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5) # this should give you the variance associated with the effect of the batch 

dev.off()

png("PCAplot_post_batchcorrection.png")

assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch) # applying limma's removeBatchEffect 

plotPCA(vsd, intgroup = c("condition", "batch")) +
  geom_text(aes(label = ddsTxi_obj$sample), vjust = 1.5)        # plotting the residuals 

dev.off()

# Replot the PCAs in ggplot 

# need to find way to replot the PCA before and after batch correction 

#################################################################################################################################

# P-VALUES and Adj p-values HISTOGRAM SANITY CHECK 

png("pvalues_histogram.png")

hist(FedvsStarved$pvalue)

dev.off()

png("padj_histogram.png")

hist(FedvsStarved$padj)

dev.off()

#################################################################################################################################
# Extract the DEseq2 Results 

library(dplyr)
library(tidyverse)

FedvsStarved_table <- as.data.frame(FedvsStarved) %>%
  rownames_to_column("GeneID") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj)

# Exctract columns of interest from annotation

deseq_annotation <- as.data.frame(annotation) %>%
  select(-c(source, type, score, phase, gene_id, level, hgnc_id, havana_gene, tag)) %>%
  rownames_to_column("GeneID") 

FedvsStarved_table <- FedvsStarved_table %>%
  left_join(deseq_annotation, "GeneID")

write_tsv(FedvsStarved_table, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/FvS_DESeq2resultstable.csv")


#################################################################################################################################

# volcano plot (both figure and interactive)

# Enhanced Volcano 

library(ggrepel)
library(EnhancedVolcano)

png("FedvsStarved_volcano.png")

EnhancedVolcano(FedvsStarved_table,
                lab = FedvsStarved_table$gene_name,
                x = "log2FC",
                y = "FDR",
                ylab = "-Log10(FDR)",
                xlim = c(-8, 8), 
                title = "Fed vs Starved Volcano Plot", 
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

de <- as.integer(FedvsStarved$padj <= 0.05 & !is.na(FedvsStarved$padj))

normCounts <- log2(counts(ddsTxi_obj))

statusCol <- match(samplefile$condition, c("Fed", "Starved"))

glXYPlot(
  x = FedvsStarved$log2FoldChange,
  y = -log10(FedvsStarved$pvalue),
  xlab = "log2FC",
  ylab = "-log10(FDR)",
  side.xlab = "Condition",
  side.ylab = "Log(CPM)",
  main = "Fed vs Starved",
  counts = normCounts,
  groups = group,
  sample.cols = statusCol,
  status = de,
  anno = FedvsStarved_table[, c("GeneID", "gene_name")],
  folder = "interactive_volcano")

#################################################################################################################################

# HEATMAPS

library(ComplexHeatmap)
library(circlize)

# get the top 200 genes
sigGenes <- as.data.frame(FedvsStarved_table) %>% 
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

# Find the most abundant transcripts FOR EACH CONDITION

# import all the salmon files in

transcript_salmon <- tximport(files, type = "salmon", txOut = TRUE)

# take the abundances for each replicate 

transcript_tpm <- transcript_salmon$counts

# change the rownames to a column called transcript_id

transcript_tpm <- as.data.frame(transcript_tpm) %>%
  rownames_to_column("transcript_id")

# add mean column ... remember, to do the SCU you need to do this for each condition

abundant_transcripts_Fed <- transcript_tpm %>%
  select(-c(Starved1, Starved2, Starved5)) %>%
  rowwise() %>%
  mutate(mean_tpm = round(mean(c(Fed1, Fed2, Fed5)))) %>%
  filter(mean_tpm != 0)

abundant_transcripts_Starved <- transcript_tpm %>%
  select(-c(Fed1, Fed2, Fed5)) %>%
  rowwise() %>%
  mutate(mean_tpm = round(mean(c(Starved1, Starved2, Starved5)))) %>%
  filter(mean_tpm != 0)

# find the mean raw counts for each of the abundant transcripts 

meanrawcounts_Fed <- as.data.frame(transcript_salmon$counts) %>%
  rownames_to_column("transcript_id") %>%
  select(-c(Starved1, Starved2, Starved5)) %>%
  rowwise() %>%
  mutate(mean_counts = round(mean(c(Fed1, Fed2, Fed5)))) %>%
  select(-c(Fed1, Fed2, Fed5)) %>%
  filter(mean_counts != 0)
  
meanrawcounts_Starved <- as.data.frame(transcript_salmon$counts) %>%
  rownames_to_column("transcript_id") %>%
  select(-c(Fed1, Fed2, Fed5)) %>%
  rowwise() %>%
  mutate(mean_counts = round(mean(c(Starved1, Starved2, Starved5)))) %>%
  select(-c(Starved1, Starved2, Starved5)) %>%
  filter(mean_counts != 0)  


# add the gene_id from the previously constructed tx 2 gene list

tx2gene_df <- as.data.frame(tx2gene)

abundant_transcripts_Fed <- dplyr::inner_join(abundant_transcripts_Fed, tx2gene_df, by = "transcript_id")

abundant_transcripts_Starved <- dplyr::inner_join(abundant_transcripts_Starved, tx2gene_df, by = "transcript_id")

# select transcripts with the biggest mean tpm ... selecting the most abundant transcript in each gene set

abundant_transcripts_Fed <- abundant_transcripts_Fed %>%
  group_by(gene_id) %>%
  top_n(n = 1, wt = mean_tpm) %>%
  filter(mean_tpm != 0)

abundant_transcripts_Starved <- abundant_transcripts_Starved %>%
  group_by(gene_id) %>%
  top_n(n = 1, wt = mean_tpm) %>%
  filter(mean_tpm != 0)

# dont delete this until i know i dont need it !!!!!
#transcripts_byTPM <- abundant_transcripts %>%
  #select(c(gene_id, transcript_id))


# sanity check for duplications ... change according to condition

duplications <- abundant_transcripts_Starved %>%
  add_count(transcript_id) %>%
  filter(n>1) %>%
  distinct(transcript_id) %>%
  nrow()

duplications

# filter the mean raw counts table so they only include the most abundant transcripts for each gene 

meanrawcounts_Fed <- dplyr::inner_join(meanrawcounts_Fed, abundant_transcripts_Fed, by = "transcript_id") %>%
  select(-c(Fed1, Fed2, Fed5, mean_tpm, gene_id))

meanrawcounts_Starved <- dplyr::inner_join(meanrawcounts_Starved, abundant_transcripts_Starved, by = "transcript_id") %>%
  select(-c(Starved1, Starved2, Starved5, mean_tpm, gene_id))

# create CDS fasta file with all expressed genes from salmon mapping step  

library(Biostrings)
library(dplyr)

# take the abundant transcripts and remove the columns i dont need 

#filtered_abundantTx_Fed <- as.data.frame(abundant_transcripts_Fed) %>%
  #select(-c(mean_tpm, gene_id)) %>%
  #column_to_rownames(var = "transcript_id")

#filtered_abundantTx_Starved <- as.data.frame(abundant_transcripts_Starved) %>%
  #select(-c(mean_tpm, gene_id)) %>%
  #column_to_rownames(var = "transcript_id")


# remove any rows where transcripts are not being expressed in the condition 

#keep1 <- rowSums(filtered_abundantTx_Fed) > 0
#filtered_abundantTx_Fed <- filtered_abundantTx_Fed[keep1, ]

#keep2 <- rowSums(filtered_abundantTx_Starved) > 0
#filtered_abundantTx_Starved <- filtered_abundantTx_Starved[keep2, ]

# row names to column 

filtered_abundantTx_Fed <-as.data.frame(abundant_transcripts_Fed) %>%
  select(-c(mean_tpm, gene_id))

filtered_abundantTx_Starved <-as.data.frame(abundant_transcripts_Starved) %>%
  select(-c(mean_tpm, gene_id))

# load the transcriptome (filtered CDS)
# all these CDSs have 5' and 3' UTRs, Start/Stop codon and is divisible by 3 

filtered_CDSfasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/gencode.v34.pc_transcripts_CDSs.fa")

# create a data frame with only the expressed transcripts 

Fed_CDS <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_Fed, by = "transcript_id") %>%
  select(-c(Fed1, Fed2, Fed5))

Starved_CDS <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(filtered_abundantTx_Starved, by = "transcript_id") %>%
  select(-c(Starved1, Starved2, Starved5))

# sanity check 

dim(Fed_CDS)
dim(meanrawcounts_Fed)

dim(Starved_CDS)
dim(meanrawcounts_Starved)

# multiple the CDS string by the mean counts ... this 

#Fed_CDS <- as.data.frame(Fed_CDS) %>%
  #dplyr::inner_join(meanrawcounts_Fed, by = "transcript_id") %>%
  #mutate_at(vars(CDS), ~ strrep(Fed_CDS$CDS, Fed_CDS$mean_counts))
  


# write the newly generated transcriptome dataframes into fasta and save on disk 

library(tidyverse)
library(seqinr)

# remeber to set a new wd 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/")

df1 = Fed_CDS
seqs = as.list(dplyr::pull(df1, CDS))
names = dplyr::pull(df1, `transcript_id`)
write.fasta(seqs, names, "Fed_CDS.fasta",
            open = "w", as.string = FALSE)

df2 = Starved_CDS
seqs = as.list(dplyr::pull(df2, CDS))
names = dplyr::pull(df2, `transcript_id`)
write.fasta(seqs, names, "Starved_CDS.fasta",
            open = "w", as.string = FALSE)

# write table wih the mean raw counts for each condition ... make sure there are no headers and tab delimited 

write.table(meanrawcounts_Fed, file="meanrawcounts_Fed.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(meanrawcounts_Starved, file="meanrawcounts_Starved.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



#################################################################################################################################

# Calculate the SCU

# double check that the directory has no dodgy fasta files in there otherwise this might mess up the calculations

library("seqinr")

aa <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/codon_aa_table.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)

head(aa)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/")
path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/"
out.file<-""
file.names <- dir(path, pattern =".fasta")
table<- NULL
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name<-file.names[i]
  name<-gsub(".fasta","",name)
  print(name)
  UCO<- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1<-do.call(rbind,UCO)
  UCO1<- as.data.frame(UCO1)
  colnames(UCO1)<- aa$AA
  head(UCO1)
  
  x<- as.matrix(UCO1, rownames.force = TRUE)
  head(x)
  x1<-aggregate(t(x),by=list(rownames(t(x))),FUN=sum)
  head(x1)
  rownames(x1)<- x1$Group.1
  x1$Group.1<-NULL
  x1<- t(x1)
  head(x1)
  rownames(UCO1) = make.names(rownames(UCO1), unique=TRUE)
  head(UCO1)
  rownames(x1) = make.names(rownames(x1), unique=TRUE)
  a <- UCO1[1,] / x1[1,] [ match( colnames(UCO1) , colnames(x1) ) ]
  
  b<- data.frame()
  for(i in 1:nrow(UCO1)){
    div<-UCO1[i,] / x1[i,] [ match( colnames(UCO1) , colnames(x1) ) ]
    
    b<- rbind(b,div)
  }
  head(b)
  colnames(b)<- aa$codon
  freq<- colMeans(b,na.rm = TRUE)
  table<-cbind(table,freq)
  colnames(table)[ncol(table)] <- paste(name)
}
head(table)
write.table(table, file = "FvS_CDS_SCUoutput.txt",sep="\t")


# upload the output file and plot the SCU

letter_to_amino <- read_tsv("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/letter_to_amino.txt")

FvS_CDS_SCUoutput <- as.data.frame(table) %>%
  rownames_to_column(var = "codon") %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(FvS_CDS_SCUoutput, file = "full_FvS_CDS_SCUoutput.txt", sep="\t")

#ggplot(FvS_CDS_SCUoutput, aes(x=Fed_CDS, y= Starved_CDS)) + geom_point()

df_melt <- mutate(FvS_CDS_SCUoutput, 
                  #Condition = case_when(variable == "Fed_differential" ~ "Fed",
                  #variable == "Starved_differential" ~ "Starved"),
                  wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                     str_detect(codon, "G$|C$") ~ "GC_wobble"))

# gingold type plot

everything_gingoldplot <- ggplot(df_melt, aes(x=Fed_CDS, y= Starved_CDS)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "Relative Synonymous Codon Usage (all transcripts)") + 
  #ggrepel::geom_text_repel(aes(label = codon))+
  ggsave(filename = "RSCU_FvS_gingold.png", bg = "white", width = 15, height = 7, dpi = 700)


everything_gingoldplot

melted<- reshape2::melt(FvS_CDS_SCUoutput)

ggplot(melted, aes(variable, y= value)) + geom_point() + facet_wrap(~aminoacid)

RSCU_FvS_Tx <- ggplot(melted, aes(variable, y= value)) + 
  geom_line(aes(group=codon))+ 
  geom_point() + 
  facet_wrap(~aminoacid) +
  ggsave(filename = "RSCU_FvS_Tx.png", bg = "white", width = 15, height = 7, dpi = 700)+
  #geom_text(data = melted, aes(label = codon),
   #         size = 3, hjust = -.5) +
  #geom_text(data = melted, aes(label = codon),
   #         size = 3, hjust = 1.5)+ 
  theme_bw()
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))

RSCU_FvS_Tx

#library(cowplot)

#save_plot(RSCU_FvS_Tx, filename = "RSCU_FvS_Tx.png")

#################################################################################################################################

# filter the up and donregulated genes to get list of differntially expressed RSCU

# i will treat all upregulated genes as Fed genes and all downregulated genes as Starved gene

upregulated <- as.data.frame(FedvsStarved) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC > 0) %>%
  dplyr::inner_join(tx2gene_df, by = "gene_id") %>%
  select(-c(gene_id))
  
downregulated <- as.data.frame(FedvsStarved) %>%
  rownames_to_column("gene_id") %>%
  dplyr::rename(log2FC = log2FoldChange, FDR = padj) %>%
  filter(FDR <= 0.05 & log2FC < 0) %>%
  dplyr::inner_join(tx2gene_df, by = "gene_id") %>%
  select(-c(gene_id))

################################################################################

Fed_differential <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(downregulated, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, stat, pvalue, FDR))

Starved_differential <- as.data.frame(filtered_CDSfasta) %>%
  rownames_to_column(var = "transcript_id") %>%
  dplyr::rename(CDS = x) %>%
  dplyr::inner_join(upregulated, by = "transcript_id") %>%
  select(-c(baseMean, log2FC, lfcSE, stat, pvalue, FDR)) 

# calculate the RSCU

# write the newly generated transcriptome dataframes into fasta and save on disk 

library(tidyverse)
library(seqinr)

# remeber to set a new wd 

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/differentials/")

df1 = Fed_differential
seqs = as.list(dplyr::pull(df1, CDS))
names = dplyr::pull(df1, `transcript_id`)
write.fasta(seqs, names, "Fed_differential.fasta",
            open = "w", as.string = FALSE)

df2 = Starved_differential
seqs = as.list(dplyr::pull(df2, CDS))
names = dplyr::pull(df2, `transcript_id`)
write.fasta(seqs, names, "Starved_differential.fasta",
            open = "w", as.string = FALSE)

#################################################################################################################################

# Calculate the SCU

# double check that the directory has no dodgy fasta files in there otherwise this might mess up the calculations

library("seqinr")

aa <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/codon_aa_table.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)

head(aa)

setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/differentials/")
path = "/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/differentials/"
out.file<-""
file.names <- dir(path, pattern =".fasta")
table<- NULL
for(i in 1:length(file.names)){
  print(i)	
  file <- read.fasta(file.names[i],as.string=FALSE)
  name<-file.names[i]
  name<-gsub(".fasta","",name)
  print(name)
  UCO<- lapply(file, function(x) uco(x, frame = 0, index = "freq", as.data.frame = FALSE, NA.rscu = NA))
  UCO1<-do.call(rbind,UCO)
  UCO1<- as.data.frame(UCO1)
  colnames(UCO1)<- aa$AA
  head(UCO1)
  
  x<- as.matrix(UCO1, rownames.force = TRUE)
  head(x)
  x1<-aggregate(t(x),by=list(rownames(t(x))),FUN=sum)
  head(x1)
  rownames(x1)<- x1$Group.1
  x1$Group.1<-NULL
  x1<- t(x1)
  head(x1)
  rownames(UCO1) = make.names(rownames(UCO1), unique=TRUE)
  head(UCO1)
  rownames(x1) = make.names(rownames(x1), unique=TRUE)
  a <- UCO1[1,] / x1[1,] [ match( colnames(UCO1) , colnames(x1) ) ]
  
  b<- data.frame()
  for(i in 1:nrow(UCO1)){
    div<-UCO1[i,] / x1[i,] [ match( colnames(UCO1) , colnames(x1) ) ]
    
    b<- rbind(b,div)
  }
  head(b)
  colnames(b)<- aa$codon
  freq<- colMeans(b,na.rm = TRUE)
  table<-cbind(table,freq)
  colnames(table)[ncol(table)] <- paste(name)
}
head(table)
write.table(table, file = "FvS_differential_SCUoutput.txt",sep="\t")

# plot the RSCU

FvS_differential_SCUoutput <- as.data.frame(table) %>%
  rownames_to_column(var = "codon") %>%
  dplyr::inner_join(aa, by = "codon") %>%
  dplyr::inner_join(letter_to_amino, by = "AA")

write.table(FvS_differential_SCUoutput, file = "full_FvS_differential_SCUoutput.txt", sep="\t")



df_melt <- mutate(FvS_differential_SCUoutput, 
                  #Condition = case_when(variable == "Fed_differential" ~ "Fed",
                                                 #variable == "Starved_differential" ~ "Starved"),
                  wobble = case_when(str_detect(codon, "A$|T$") ~ "AT_wobble",
                                     str_detect(codon, "G$|C$") ~ "GC_wobble"))


gingoldplot <- ggplot(df_melt, aes(x=Fed_differential, y= Starved_differential)) +
  geom_point(aes(color = wobble)) +
  theme_bw(base_size = 16)+
  labs(title = "Relative Synonymous Codon Usage in differential genes") + 
  #ggrepel::geom_text_repel(aes(label = codon))+
  ggsave(filename = "RSCU_diff_FvS_gingold.png", bg = "white", width = 15, height = 7, dpi = 700)


gingoldplot

for (aminoacid in c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "iMet", "STOP")) {
  
  df_loop <- subset(df_melt, aminoacid == aminoacid)
  
  gingold_loop <- ggplot(df_loop, aes(x=Fed_differential, y= Starved_differential, group = aminoacid)) +
    geom_point(aes(color = wobble)) +
    theme_bw(base_size = 16)+
    labs(title = aminoacid, "Relative Synonymous Codon Usage in differential genes") + 
    ggrepel::geom_text_repel(aes(label = codon))+
    ggsave(paste(aminoacid, "RSCU_diff_FvS_gingold.png", sep = "_"), bg = "white", width = 15, height = 7, dpi = 700)
  
}

#df_melt <- tidyr::gather(FvS_differential_SCUoutput, variable, value, Fed_differential:Starved_differential)


melted<- reshape2::melt(FvS_differential_SCUoutput)

ggplot(melted, aes(variable, y= value)) + geom_point() + facet_wrap(~aminoacid)

right_label <- melted %>%
  filter(variable == "Starved_differential")

left_label <- melted %>%
  filter(variable == "Fed_differential")

RSCU_diff_FvS_Tx <- ggplot(melted, aes(variable, y= value)) + 
  geom_line(aes(group=codon))+ 
  geom_point() +
  #geom_text(hjust = 0, nudge_x = -1.5, aes(label = codon))+
  facet_wrap(~aminoacid) +
  ggsave(filename = "RSCU_diff_FvS_Tx.png", bg = "white", width = 15, height = 7, dpi = 700)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = "none")+
  geom_text(data = right_label, aes(label = codon),
            size = 3, hjust = -.5) +
  geom_text(data = left_label, aes(label = codon),
            size = 3, hjust = 1.5) +
  geom_text_repel()+
  labs(title = "geom_text_repel()")
#theme_bw()+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

RSCU_diff_FvS_Tx




#################################################################################################################################

df_melt <- FvS_CDS_SCUoutput

df_melt <- tidyr::gather(df_melt, variable, value, Fed_CDS:Starved_CDS)

df_melt <- mutate(df_melt, Condition = case_when(variable == "Fed_CDS" ~ "Fed",
                                                  variable == "Starved_CDS" ~ "Starved"))


ggplot(df_melt, aes(Condition)) +
  geom_line(aes(group = variable)) #+
  #geom_point(aes(color = AA))


#for amino in c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "iMet", "STOP")) {
  
  df_melt <- subset(df_melt, aminoacid == amino)
  
  ggplot(df_melt, aes(codon)) + 
    geom_bar(aes(fill = Condition), position = "fill") + 
    coord_flip() + 
    theme_bw()+
    ggtitle(aminoacid, "Fed vs Starved - Relative Synonymous Codon Usage")
#}  

# stat in identity as opposed to fill 

ggplot(df_melt, aes(codon, y = variable)) + 
  geom_bar(aes(fill = variable), stat = "identity") + 
  coord_flip()




ggplot(df_melt, aes(Condition)) +
  geom_point() #+
  #scale_color_manual(values = c("red", "orange", "yellow", "green", "blue", "purple", "violet"))
  #geom_line(aes(group = codon))

#################################################################################################################################
setwd("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/")
#################################################################################################################################

###############################################################################################################################

# calculate AA freq for each gene and normalise to their abundance (TPM)

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Fed_allTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/SCU_all/Fed_CDS_allTx.fasta")

Starved_allTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/SCU_all/Starved_CDS_allTx.fasta")

# translate the DNA sequences into protein sequence 

Fed_allTx_translated_fasta <- Biostrings::translate(Fed_allTx_fasta) 

Starved_allTx_translated_fasta <- Biostrings::translate(Starved_allTx_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Fed_allTx_translated_freq <- Biostrings::alphabetFrequency(Fed_allTx_translated_fasta)

Starved_allTx_translated_freq <- Biostrings::alphabetFrequency(Starved_allTx_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Fed_allTx_translated_fasta <- as.data.frame(Fed_allTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Fed_allTx_translated_freq <- as.data.frame(Fed_allTx_translated_freq)

rownames(Fed_allTx_translated_freq) <- Fed_allTx_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Fed_allTx_translated_freq <- Fed_allTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Fed_allTx_translated_freq_normalized <- Fed_allTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Resec_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Fed_allTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Fed_AA_freq_normalized <- Fed_allTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Fed_AA_freq_normalized <- as.data.frame(Fed_AA_freq_normalized) %>%
  rownames_to_column("AA") %>%
  rename(Fed_AA_freq = V1) 

###################################################################################################

Starved_allTx_translated_fasta <- as.data.frame(Starved_allTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Starved_allTx_translated_freq <- as.data.frame(Starved_allTx_translated_freq)

rownames(Starved_allTx_translated_freq) <- Starved_allTx_translated_fasta$transcript_id

Starved_allTx_translated_freq <- Starved_allTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Starved_allTx_translated_freq_normalized <- Starved_allTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Regen_all, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Starved_allTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id")

Starved_AA_freq_normalized <- Starved_allTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Starved_AA_freq_normalized <- as.data.frame(Starved_AA_freq_normalized) %>%
  rownames_to_column("AA") %>%
  rename(Starved_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

FedvsStarved_AA_freq_normalized <- Fed_AA_freq_normalized %>%
  dplyr::inner_join(Starved_AA_freq_normalized, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(FedvsStarved_AA_freq_normalized, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/FedvsStarved_AA_freq_normalized.txt")

###############################################################################################################################

# calculate AA freq for the most abundant transcript and normalise to their abundance (TPM)

library(Biostrings)

tpm_Fed_abundant <- abundant_transcripts_Fed %>%
  select(-c("Fed1", "Fed2", "Fed5", "gene_id"))


tpm_Starved_abundant <- abundant_transcripts_Starved %>%
  select(-c(Starved1, Starved2, Starved5, gene_id))

# read in the fasta file of all expressed transcripts from earlier 

Fed_abTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/SCU_abundant/Fed_CDS.fasta")

Starved_abTx_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/SCU_abundant/Starved_CDS.fasta")

# translate the DNA sequences into protein sequence 

Fed_abTx_translated_fasta <- Biostrings::translate(Fed_abTx_fasta) 

Starved_abTx_translated_fasta <- Biostrings::translate(Starved_abTx_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Fed_abTx_translated_freq <- Biostrings::alphabetFrequency(Fed_abTx_translated_fasta)

Starved_abTx_translated_freq <- Biostrings::alphabetFrequency(Starved_abTx_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Fed_abTx_translated_fasta <- as.data.frame(Fed_abTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Fed_abTx_translated_freq <- as.data.frame(Fed_abTx_translated_freq)

rownames(Fed_abTx_translated_freq) <- Fed_abTx_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Fed_abTx_translated_freq <- Fed_abTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Fed_abTx_translated_freq_normalized <- Fed_abTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Fed_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Fed_abTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  select(-c(gene_id)) %>%
  na.omit()

Fed_AA_freq_normalized_ab <- Fed_abTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # divide sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Fed_AA_freq_normalized_ab <- as.data.frame(Fed_AA_freq_normalized_ab) %>%
  rownames_to_column("AA") %>%
  rename(Fed_AA_freq = V1) 

###################################################################################################

Starved_abTx_translated_fasta <- as.data.frame(Starved_abTx_translated_fasta) %>%
  rownames_to_column("transcript_id")

Starved_abTx_translated_freq <- as.data.frame(Starved_abTx_translated_freq)

rownames(Starved_abTx_translated_freq) <- Starved_abTx_translated_fasta$transcript_id

Starved_abTx_translated_freq <- Starved_abTx_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Starved_abTx_translated_freq_normalized <- Starved_abTx_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Starved_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Starved_abTx_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  select(-c(gene_id)) %>%
  na.omit()

Starved_AA_freq_normalized_ab <- Starved_abTx_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Starved_AA_freq_normalized_ab <- as.data.frame(Starved_AA_freq_normalized_ab) %>%
  rownames_to_column("AA") %>%
  rename(Starved_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

FedvsStarved_AA_freq_normalized_ab <- Fed_AA_freq_normalized_ab %>%
  dplyr::inner_join(Starved_AA_freq_normalized_ab, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(FedvsStarved_AA_freq_normalized_ab, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/FedvsStarved_AA_freq_normalized_ab.txt")

###################################################################################################

library(Biostrings)

# read in the fasta file of all expressed transcripts from earlier 

Fed_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/SCU_differentials/Fed_differential_cutoff2.fasta")

Starved_diff_fasta <- readDNAStringSet("/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/SCU/SCU_differentials/Starved_differential_cutoff2.fasta")

# translate the DNA sequences into protein sequence 

Fed_diff_translated_fasta <- Biostrings::translate(Fed_diff_fasta) 

Starved_diff_translated_fasta <- Biostrings::translate(Starved_diff_fasta)

# calculate the frequency of each amino acid in the sequence for each transcript 

Fed_diff_translated_freq <- Biostrings::alphabetFrequency(Fed_diff_translated_fasta)

Starved_diff_translated_freq <- Biostrings::alphabetFrequency(Starved_diff_translated_fasta)

###################################################################################################

# normalize the AA frequency to abundance by using the TPM values for each transcript 

Fed_diff_translated_fasta <- as.data.frame(Fed_diff_translated_fasta) %>%
  rownames_to_column("transcript_id")

Fed_diff_translated_freq <- as.data.frame(Fed_diff_translated_freq)

rownames(Fed_diff_translated_freq) <- Fed_diff_translated_fasta$transcript_id

# remove the excess columns (other AA's or unidentifiable AA's)

Fed_diff_translated_freq <- Fed_diff_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

# multiply frequencies by the TPM to normalize to abundance 

Fed_diff_translated_freq_normalized <- Fed_diff_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Fed_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Fed_diff_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  na.omit()

Fed_AA_freq_normalized_diff <- Fed_diff_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>% # multiply sum of each column by sum of matrix to get composition %
  t(.) # transpose matrix 

# clean up column and rownames

Fed_AA_freq_normalized_diff <- as.data.frame(Fed_AA_freq_normalized_diff) %>%
  rownames_to_column("AA") %>%
  rename(Fed_AA_freq = V1) 

###################################################################################################

Starved_diff_translated_fasta <- as.data.frame(Starved_diff_translated_fasta) %>%
  rownames_to_column("transcript_id")

Starved_diff_translated_freq <- as.data.frame(Starved_diff_translated_freq)

rownames(Starved_diff_translated_freq) <- Starved_diff_translated_fasta$transcript_id

Starved_diff_translated_freq <- Starved_diff_translated_freq %>%
  select(-c("*", "-", "+", ".", other, U, O, B, X, J, Z))

Starved_diff_translated_freq_normalized <- Starved_diff_translated_freq %>% 
  rownames_to_column("transcript_id") %>% 
  left_join(tpm_Starved_abundant, by = "transcript_id") %>% 
  mutate_at(vars(colnames(Starved_diff_translated_freq)), list(~ . * mean_tpm)) %>% 
  select(-mean_tpm) %>%
  column_to_rownames(var = "transcript_id") %>%
  na.omit()

Starved_AA_freq_normalized_diff <- Starved_diff_translated_freq_normalized %>%
  select(where(is.numeric)) %>%
  summarize(across(.fn = sum) / sum(across())) %>%
  t(.)

Starved_AA_freq_normalized_diff <- as.data.frame(Starved_AA_freq_normalized_diff) %>%
  rownames_to_column("AA") %>%
  rename(Starved_AA_freq = V1) 

###################################################################################################

# combine newly generated tables and add further details for plotting 

FvS_AA_freq_normalized_diff <- Fed_AA_freq_normalized_diff %>%
  dplyr::inner_join(Starved_AA_freq_normalized_diff, by = "AA") %>%
  dplyr::inner_join(letter_to_amino, by = "AA") %>%
  select(- full_amino, - AA)

write_tsv(FvS_AA_freq_normalized_diff, "/mnt/data/BMOHAMED/Total_RNAseq/salmon/protein_coding/FvS_AA_freq_normalized_differential_logFC2cutoff.txt")


###################################################################################################

# AA property plots

aa_properties <- read.table("/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90/aa_properties.txt", skip=1)

colnames(aa_properties) <- c("aminoacid", "property")

property_plot_input <- FvS_AA_freq_normalized_diff %>%
  left_join(aa_properties, "aminoacid")

ggplot(data = property_plot_input, aes(x=Fed_AA_freq, y= Starved_AA_freq))+
  geom_point(aes(color = property)) +
  theme_bw(base_size = 16)+
  labs(title = "Amino Acid Frequency in BJ5TA Fed vs Starved") +
  geom_abline(slope = 1, intercept = 0, linetype= "dashed") +
  xlab("Fed AA frequency")+
  ylab("Starved AA frequency") +
  geom_label_repel(aes(label = aminoacid),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()  

ggsave(filename = "FvS_aa_properties.png", bg = "white", width = 7, height = 7, dpi = 600)




















































