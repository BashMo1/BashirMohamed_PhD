#!/bin/bash

maindirectory="/mnt/data/BMOHAMED/ALICE-tRNA-seq/MDM2/201006_NS500360_0181_AH3CYVBGXG/fastq_files"

outdir=${maindirectory}/trimmed

mkdir $outdir 

Samples="M_MDM2_Null1_S1 M_MDM2_Null2_S2 M_MDM2_Null3_S3 M_MDM2_Null4_S4 M_MDM2_Null5_S5 M_MDM2_Cre1_S6 M_MDM2_Cre2_S7 M_MDM2_Cre3_S8 M_MDM2_Cre4_S9 M_MDM2_Cre5_S10 F_MDM2_Null1_S11 F_MDM2_Null2_S12 F_MDM2_Null3_S13 F_MDM2_Null4_S14 F_MDM2_Cre2_S15 F_MDM2_Cre3_S16 F_MDM2_Cre4_S17"

for sample in $Samples

#3' Adapter trimming 

do 

    cutadapt -a 'CAGATCGGAAGAGCACACGTCT' -o trimmed_${sample}_R1_001.fastq ${maindirectory}/${sample}_R1_001.fastq > ${sample}_trimming_summary.txt & 

done 

