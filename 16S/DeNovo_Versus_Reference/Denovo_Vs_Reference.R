setwd("~/Ashley_Analysis/")
data<- read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_collapse_even321798.txt", header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")

rdp<- data[,19]
data<-data[,-19]
data=data[,order(colnames(data))]

data.dn <- data[grepl("OTU", row.names(data)),]
Percentdn<-colSums(data.dn)/colSums(data)

#read in mapping file with soil data
map=read.table("Centralia_Collapsed_Map_forR.txt", header=TRUE, sep="\t")
temp<- map$SoilTemperature_to10cm

plot(temp, colSums(data.dn)/colSums(data), xlab = "Sample Temperature", ylab="% De novo OTUs", main = "Fraction De Novo versus Samnple Temperature")

# Correlation Tests for % denovo and Sample Temperature
cor.test(temp, Percentdn, method="pearson")
cor.test(temp, Percentdn, method="spearman")
cor.test(temp, Percentdn, method="kendall")

# Format data for ggplot2 stacked barcharts
Plot_data<-cbind(colnames(data.dn),rep("Denovo",18), colSums(data.dn))
Plot_data<-rbind(Plot_data, cbind(colnames(data.dn),rep("Reference", 18), (colSums(data)- colSums(data.dn))))
colnames(Plot_data)<- c("Sample", "OTU_Type", "Count")

library(ggplot2)
Plot_data<- as.data.frame(Plot_data)
Counts <- c(colSums(data.dn), (colSums(data) -colSums(data.dn)))
Plot_data$Count<- Counts


p=ggplot(data=Plot_data, aes(x=Sample, y=Count, fill=OTU_Type)) + geom_bar(stat="identity", aes(x=Sample, y=Count)) 
p

# Test % denovo versus fireaffected and cool sites
t.test(Percentdn~(map$SoilTemperature_to10cm>=15)
)


### t test for % denovo in Recovered versus Fireaffected
Percent_dn_norefs<- Percentdn[-17]
Percent_dn_norefs<- Percent_dn_norefs[-8]
Classification<- map$Classification
Classification<- Classification[-17]
Classification<- Classification[-8]

t.test(Percent_dn_norefs~Classification)



### Trying to figure out how to order the samples on the bottom by sample temperature