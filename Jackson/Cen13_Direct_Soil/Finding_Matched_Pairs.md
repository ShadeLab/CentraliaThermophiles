I think megahit requires paired end reads to have an exact match in their matching files, if they do not it freezes up and does not work. This is problematic because when using fastq_quality_filter, we do so separately on the forward and reverse read, and thus need to adjust our two input files so there are not reads without a mate. Start with the following command. 
```
grep "@HWI" Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.qc.fastq | awk '{print $1}' > Cen13_R1_Headers.txt
```
This command will take the headers from the quality filtered forward reads, and writes the 1st portion, which we will use to match the reverse reads, to a file called `Cen13_R1_Headers.txt`. Then we can find the matches that are present in the reverse read by doing the following. 
```
#! /bin/bash
for file in $(<Cen13_R1_Headers.txt)
do
	grep ${file} Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.qc.fastq -A 3 >> Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.qc.matches.fastq
done
```


Next we need to get the headers from the matched reads so we can grab the correct reads from the 1st read.
```
grep "@HWI" Cen13_mgDNA_Pooled_CTTGTA_L002_R2_001.qc.matches.fastq | awk '{print $1}' > Cen13_Matched_Headers.txt
```

Make the matched forward read fastq file
```
#! /bin/bash
for file in $(<Cen13_Matched_Headers.txt)
do
	grep ${file} Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.qc.fastq -A 3 >> Cen13_mgDNA_Pooled_CTTGTA_L002_R1_001.qc.matches.fastq
done
```