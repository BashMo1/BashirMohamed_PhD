# Adapter_trimming.sh 

trims the Illumina sequencing adapters from the raw fastq input files, generating the processed files with prefix trimmed_

bash Adapter_trimming.sh

# processsamples.py

This can be used for mapping but can also use mapreads.py (or .bash or .pyc), using the following options: 

time python processsamples.py --experimentname=XXX --databasename=hg38_nov2019_new- --samplefile=samplefile_fedvsstarved.txt --ensemblgtf=hg38-genes_ensembl2019.gtf --exppairs=pairfile_fedvsstarved.txt   

xxx = can be named as anything

Holmes AD, Howard JM, Chan PP, and Lowe TM. tRNA Analysis of eXpression (tRAX): A tool for integrating analysis of tRNAs, tRNA-derived small RNAs, and tRNA modifications. 2020.

# countreads.py

can be used to count the reads using the following options:

time python countreads.py --samplefile=samplefile_fedvsstarved.txt --ensemblgtf=hg38-genes_ensembl2019.gtf --nomultimap --trnatable=hg38_nov2019_new-trnatable.txt --maturetrnas=hg38_nov2019_new-maturetRNAs.bed

Holmes AD, Howard JM, Chan PP, and Lowe TM. tRNA Analysis of eXpression (tRAX): A tool for integrating analysis of tRNAs, tRNA-derived small RNAs, and tRNA modifications. 2020.

# FedvsStarved_countsanalysis.R 

counts processing and differential expression analysis. Also has codes for visualization. sampleinformation.txt file and counts text file needed as input.

# files with prefix hg38_nov2019_new-

required by the processsamples.py, countreads.py and mapping packages 


