#! /bin/bash
for file in $(<Cen13_Matched_Headers.txt)
do
	grep ${file} Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.qc.fastq -A 3 >> Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.qc.matches.fastq
done
