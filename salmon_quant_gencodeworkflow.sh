#!/bin/bash

maindirectory="/mnt/data/BMOHAMED/Total_RNAseq/salmon/IMR90" 

Samples="Empty1_S1 Empty2_S2 Empty3_S3 Ras1_S4 Ras2_S5 Ras3_S6 Ras4_S7"

for sample in $Samples

# salmon quantification at against the transcriptome 

do 

    salmon quant -i gencode.v34.pc_transcripts_f_index/ -l A -r ${maindirectory}/${sample}_R1_001.fastq -p 12 -o ${sample}_filteredproteincoding

done 
