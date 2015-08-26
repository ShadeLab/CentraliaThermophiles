#! /bin/bash
for file in $(<Cen13_R1_Headers.txt)
do
	grep ${file} Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.qc.fastq -A 3 >> Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.qc.matches.fastq
done
