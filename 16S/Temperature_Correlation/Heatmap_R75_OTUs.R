load("~/Desktop/VennAnalysis/Temp_Correlation2.RData")

map=read.table("Centralia_Full_Map_Fixed.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE, stringsAsFactors=FALSE)#read in the mapping file
sorted_map<- map[order(map$SoilTemperature_to10cm),]
Temp<-sorted_map$SoilTemperature_to10cm 
Temp<-unique(Temp)
Temp
Full<-c(Temp[1:3],13.3,Temp[4:17])
Column_order<-unique(sorted_map$Sample)

OTUs_R75=coretest.out[coretest.out[,2]>0.75,]
Soils_75=collapsed_soils_rel[OTUs_R75[,1],]
hmap<-as.matrix(Soils_75)

ordered_hmap=NULL
for (i in 1:18){ordered_hmap<- cbind(ordered_hmap, hmap[,Column_order[i]])}
colnames(ordered_hmap)<- Column_order

heatmap(ordered_hmap, Rowv = NA, Colv = NA  , scale = "row", main = "R=0.75 OTUs HeatMap")
library(gplots)
heatmap.2(ordered_hmap, Colv=FALSE, Rowv = TRUE, col=pal(20), dendrogram="row" , scale = "row", main = "R>0.75 OTUs HeatMap", key=TRUE, trace="none", labRow=c(rep("",217)), labCol=Full, xlab= "Soil Temperature (C)", ylab="OTUs (97%)" )

library(grDevices)
library(RColorBrewer)
display.brewer.all()


help(brewer.pal)
colors<- brewer.pal(3, "YlOrRd")
pal <- colorRampPalette(colors)

soils=read.table("subsamplingtable_even73419_R.txt", header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE)
soils=as.matrix(soils)
subset_otus<- soils[rownames(hmap),]


taxonomy.df=NULL
for (i in 1:nrow(subset_otus)){
taxonomy.df<-rbind(taxonomy.df, as.vector(unlist(strsplit(subset_otus[i,55],";")),mode="list"))
}

Phyla=unique(taxonomy.df[,2])
str(unlist(Phyla[1]))
Phyla_OTUs=NULL
for (i in 1:length(Phyla)){
  Phyla_OTUs=rbind(Phyla_OTUs, c(Phyla[i], sum(1*(taxonomy.df[,2]==unlist(Phyla[i])))))
}

barplot(as.numeric(Phyla_OTUs$V2), names.arg=Phyla_Names)


Phyla_Names=c("Acidobacteria", "Proteobacteria", "WS3", "Verrucomicrobia", "GOUTA4", "Bacteroidetes", "Planctomycetes", "Chlorobi", "Chloroflexi", "Nitrospirae", "Gemmatimonadetes", "Spirochaetes", "Firmicutes", "Actinobacteria", "GN04", "Cyanobacteria")
x<- barplot(as.numeric(Phyla_OTUs$V2), main="Strongly Heat Responsive OTUs", ylab="Number of OTUs", xlab="Phyla")
labs<- paste(Phyla_Names)
text(cex=.75, x=x+.25, y=-1.25, labs, xpd=TRUE, srt=60, pos=2)


sorted_Phyla<- Phyla_OTUs[order(-as.numeric(Phyla_OTUs$V2)),]

x_ordered<- barplot(as.numeric(sorted_Phyla$V2), main="Strongly Temperature Responsive OTUs", ylab="Number of OTUs", xlab="Phyla")
labs_ordered<-paste(sorted_Phyla$V1)
text(cex=.75, x=x+.25, y=-1.25, labs, xpd=TRUE, srt=60, pos=2)

