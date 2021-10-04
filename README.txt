The code in this repository is a compilation of different programs used in the processing, mapping, analysis and visualization of ALICE-tRNA-seq and mRNA-seq data sets used in my PhD. ALICE-tRNA-seq is a high-throuput RNA sequencing methodology I developed during my PhD that allows us to decipher the relative changes in tRNA expression at the gene level, in varying cellular and disease conditions. The mRNA-seq analysis is done to calculate the relative synonymous codon usage of translating genes.

tRNAs are small adapter RNA molecules used by the cell to decode the mRNA, bringing amino acids to the ribosomes, aiding in protein synthesis. Molecular and cancer biologists have long tried to decipher the role of the tRNAs in protein synthesis, and how the tRNA pool may change according to differing cellular states, and how they may affect protein synthesis in disease state. 

This methodology utilizes Illumina High-throughput sequencing, which allows the generation of pricise and large datasets, utilized in the study of relative tRNA expression. tRNA-omics, a subfield of genomics, fits in beautifully with the already established transcriptomics, proteomics and metabolomics  

Although raw data will not be found on this page (due to intellectual property), the programs found here were either written by myself, or adapted from publications (and this will be cited in text). 

# Adapter_trimming.sh 

representative script that trims the Illumina sequencing adapters from the raw fastq input files, generating the processed files with prefix trimmed_

# processsamples.py

This can be used for mapping but can also use mapreads.py (or .bash or .pyc), using the following options: 

time python processsamples.py --experimentname=XXX --databasename=hg38_nov2019_new- --samplefile=samplefile_fedvsstarved.txt --ensemblgtf=hg38-genes_ensembl2019.gtf --exppairs=pairfile_fedvsstarved.txt   

xxx = can be named as anything

adapted from: Holmes AD, Howard JM, Chan PP, and Lowe TM. tRNA Analysis of eXpression (tRAX): A tool for integrating analysis of tRNAs, tRNA-derived small RNAs, and tRNA modifications. 2020.

# countreads.py

can be used to count the reads using the following options:

time python countreads.py --samplefile=samplefile_fedvsstarved.txt --ensemblgtf=hg38-genes_ensembl2019.gtf --nomultimap --trnatable=hg38_nov2019_new-trnatable.txt --maturetrnas=hg38_nov2019_new-maturetRNAs.bed

adapted from: Holmes AD, Howard JM, Chan PP, and Lowe TM. tRNA Analysis of eXpression (tRAX): A tool for integrating analysis of tRNAs, tRNA-derived small RNAs, and tRNA modifications. 2020.


# files with prefix hg38_nov2019_new-

required by the processsamples.py, countreads.py and mapping packages. These are generated from the tRNA-Scan-SE programme (http://lowelab.ucsc.edu/tRNAscan-SE/) with the genome of interest being used as input.

# generation of sample and pair files

samplefile = tab delimited .txt file with unique sample names, condition and path to fastq input files in three respective columns 

e.g. 
Replicate1  Condition1 /path/to/Replicate1_fastq
Replicate2  Condition1 /path/to/Replicate2_fastq
Replicate3  Condition2 /path/to/Replicate3_fastq
Replicate4  Condition2 /path/to/Replicate4_fastq

pairfile = tab delimited .txt file with the two conditions you want to analyze. The first column should be the control sample

e.g. 
Condition1 Condition2

# .GTF file

Attached is an Ensemble GTF file containing non-coding RNA annotations for Humans from the hg38 annotation, with the chromosome names being converted to the UCSC Genome Browser style. If you have a an up to date GTF file or require one for a different organism, then download your GTF file of interest from the Ensemble database and convert to the correct formatt:

wget -q -O - ftp://ftp.ensembl.org/pub/release/of/interest.gtf.gz | gzip -cd | grep -v '^#' | awk '{print "chr" $0;}' | grep -e Mt_rRNA -e miRNA -e misc_RNA -e rRNA -e snRNA -e snoRNA -e ribozyme -e sRNA -e scaRNA > XXX-genes.gtf

where XXX can be names at your convenience 

# System requirements

 Due to the large size of sequencing datasets, the programs listed here are to be run on a Linux/Unix system with at least 8 cores and 16 GB memory.
