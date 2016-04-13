Genome_Phyla<-read.table("/Users/JSorensen/Desktop/Genome_Phyla.txt", header=TRUE, sep="\t")
tail(Genome_Phyla)
Therm_Genome_Phyla<-read.table("/Users/JSorensen/Desktop/Genome_Phyla_Thermophiles.txt", header=TRUE, sep="\t")

Therm_Phyla_Counts<-as.data.frame(table(Therm_Genome_Phyla$Phylum))


Phyla_Counts<-as.data.frame(table(Genome_Phyla$Phylum))
sum(Phyla_Counts$Freq)
slices<- Phyla_Counts$Freq
lbls<- Phyla_Counts$Var1
pie(slices, labels=lbls, main="Pie Chart of Phyla")

ggplot(Phyla_Counts, aes(Phyla_Counts$Freq, fill=cut)) + geom_bar()