#!/bin/bash

maindirectory="/mnt/data/BMOHAMED/ALICE-tRNA-seq/fedvsstarved_rvcomp"

outdir=${maindirectory}/trimmed

mkdir $outdir 

Samples="Fed1_rvcomp_S1 Fed2_rvcomp_S2 Fed5_rvcomp_S3 Starved1_rvcomp_S4 Starved2_rvcomp_S5 Starved5_rvcomp_S6"

for sample in $Samples

#3' Adapter trimming 

do 

    cutadapt -a 'CAGATCGGAAGAGCACACGTCT' -o trimmed_${sample}_R1_001.fastq ${maindirectory}/${sample}_R1_001.fastq > ${sample}_trimming_summary.txt & 

done 

