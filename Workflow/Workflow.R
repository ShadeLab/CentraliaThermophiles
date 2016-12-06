## Workflow Script for Sorensen_INPREP
library(ggplot2)
library(vegan)
library(outliers)
library(gplots)
library(colorRamps)
### Reading in Data Files and Manipulating Datafile

# Mapping File
map <- read.table("Centralia_Collapsed_Map_forR.txt",sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)

map_MG <- map[c(1,3,4,5,6,7,10,12,14,15,16,17),]

# 16S OTU Table
comm <- read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)
rdp<- comm[,19]
rdp.sigs <- rdp
comm<-comm[,-19]
comm=comm[,order(colnames(comm))]

#designate a full dataset
comm.sigs=comm

#remove OTUs with an abundance = 1, across the entire dataset (singleton OTUs)
comm=comm[rowSums(comm)>1,]
sum(colSums(comm))
comm_rel <- decostand(x=comm, method="total", MARGIN=2)
colSums(comm_rel)

rdp <- rdp[rowSums(comm.sigs)>1]
# KEGG Ortholog Table
KO <- read.table("KO_minus2col_09-27-2016.tab.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
colnames(KO) <- map_MG$Sample

rpoB <- KO["KO:K03043",] + KO["KO:K13798",] ### add the bacterial and archaeal copies together.
#Relativze to rpoB
KO.rpob <- NULL
for (i in 1:nrow(KO)){
  KO.rpob <- rbind(KO.rpob, KO[i,]/rpoB)
}

### Venn Analysis 16S 
source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
library("limma")

fireclass=map[,"Classification"]

active.pa=1*(comm[,fireclass=="FireAffected"]>0)
recov.pa=1*(comm[,fireclass=="Recovered"]>0)
ref.pa=1*(comm[,fireclass=="Reference"]>0)

#summarize and combine
venndata=cbind(1*rowSums(active.pa>0),1*rowSums(recov.pa>0),1*rowSums(ref.pa>0))
colnames(venndata)=c("Fire-affected", "Recovered", "Reference")

#apply venn analysis
v=vennCounts(venndata)
v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)
v[,"Counts"]<-v2
vennDiagram(v)


#Write out the results of venncounts

write.table(v, "VennCounts.txt", quote=FALSE, sep="\t")
par(mfrow=c(1,1))
setEPS()

postscript("../Figures/Fig1.eps", width = 5, height=5, pointsize=8,paper="special")
par(mar=c(5,3,2,2)+0.1)
fig1=vennDiagram(v)
dev.off()

### DeNovo Vs Reference OTU Analysis FINISHED
library(reshape2)
library(ggplot2)
Percent_DN <- colSums(comm[grepl("dn",row.names(comm)),])/colSums(comm)
map$PercentDeNovo <- Percent_DN

map.long=melt(map, id.vars=c("Sample", "SoilTemperature_to10cm", "Classification"), measure.vars=c("NO3N_ppm","NH4N_ppm","pH","SulfateSulfur_ppm","K_ppm","Ca_ppm","Mg_ppm","OrganicMatter_500","Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","PercentDeNovo"))

p <- ggplot(map, aes(SoilTemperature_to10cm, PercentDeNovo))
p+ geom_point() + labs(y="Community Novelty", x="Soil Temperature")
ggsave("../Figures/Fig2.eps", width=178, units="mm")

#Correlation test between Novelty and Soil Temperature
cor.test(map$SoilTemperature_to10cm, map$PercentDeNovo)

fit <- lm(map$PercentDeNovo ~ map$SoilMoisture_Per)
summary(fit)
anova(fit)


### Broken
#fig2=ggplot(map.long, aes(y=as.numeric(SoilTemperature_to10cm), x=value))
#add points layer
#geom_point(aes(y=as.numeric(SoilTemperature_to10cm), x=value, shape=Classification, color=as.numeric(SoilTemperature_to10cm)))
### Broken



### OTU abundance and Soil Temperature Correlations
comm_rel.mat <- as.matrix(comm_rel)

Temp_Full <- as.numeric(map$SoilTemperature_to10cm)
#Perform Correlation tests 
coretest.out=NULL
for(i in 1:nrow(comm_rel.mat)){
  results=cor.test(comm_rel.mat[i,],Temp_Full)
  coretest.out=rbind(coretest.out,c(row.names(comm_rel.mat)[i],results$estimate,results$p.value))
}

hist(as.numeric(coretest.out[,2]))
# Only Significant Correlations with Temperature
sigcor <- coretest.out[coretest.out[,3]<.05,]
comm_sigcor <- cbind(comm_rel.mat[coretest.out[,3]<.05,],rdp[coretest.out[,3]<.05])
#Only Postive Significant Correlations with Temp
sigcor_pos <- sigcor[sigcor[,2]>0,]
comm_sigcor_pos <- comm_sigcor[sigcor[,2]>0,]
comm_sigcor_pos <- as.data.frame(comm_sigcor_pos)
hist(as.numeric(sigcor_pos[,2]))

# Designate a dataset with read counts (as opposed to relative abundance) 
Acomm_sigcor <- cbind(comm[coretest.out[,3]<.05,],rdp[coretest.out[,3]<.05])
Acomm_sigcor_pos <- Acomm_sigcor[sigcor[,2]>0,]
Acomm_sigcor_pos <- as.data.frame(Acomm_sigcor_pos)

# Getting Taxonomy for temperature correlated OTUs
Taxonomy <- NULL
for (i in 1:nrow(comm_sigcor_pos)){
  Taxonomy <- rbind(Taxonomy,unlist(strsplit(as.character(comm_sigcor_pos[i,19]), ";")))
}
#Phyla with Temperature Correlated (AKA Thermophilic) OTUs
Phyla_Therm <- unique(Taxonomy[,2])
Phyla_Therm <- sub("k__","",Phyla_Therm)
Phyla_Therm <- sub(" p__", "", Phyla_Therm)
Phyla_Therm <- gsub("\\[|\\]","", Phyla_Therm)

Counts<- NULL
Num_OTUs <- NULL
for (i in 1:length(Phyla_Therm)){
  Object <- Acomm_sigcor_pos[grep(Phyla_Therm[i],Taxonomy[,2]),] 
  Counts <- c(Counts, sum(Object[,1:18]))
  Num_OTUs <- c(Num_OTUs, nrow(Object))
}

Phylum_Counts_Therm <- data.frame(matrix(NA, nrow = length(Phyla_Therm), ncol = 2))
Phylum_Counts_Therm[,1] <- Phyla_Therm
Phylum_Counts_Therm[,2] <- Counts
colnames(Phylum_Counts_Therm) <- c("Phylum", "Reads")
Phylum_Counts_Therm$Percentage_Therm <- Phylum_Counts_Therm[,2]/sum(Phylum_Counts_Therm[,2])
sum(Phylum_Counts_Therm[,3])

Phylum_OTUs_Therm <- data.frame(matrix(NA, nrow = length(Phyla_Therm), ncol = 2))
Phylum_OTUs_Therm[,1] <- Phyla_Therm
Phylum_OTUs_Therm[,2] <- Num_OTUs
colnames(Phylum_OTUs_Therm) <- c("Phylum", "Number_of_Temperature_Responsive_OTUs")
Phylum_OTUs_Therm$PercentageOTUs <- Phylum_OTUs_Therm[,2]/sum(Phylum_OTUs_Therm[,2])
sum(Phylum_OTUs_Therm$PercentageOTUs)

Phylum_Therm <- cbind(Phylum_Counts_Therm, Phylum_OTUs_Therm[,2:3])
Phylum_Therm[,1] <- sub("p__", "", Phylum_Therm[,1])
Phylum_Therm[,1] <- sub("k__", "", Phylum_Therm[,1])
Phylum_Therm[,1] <- sub(" ", "", Phylum_Therm[,1])
Phylum_Therm[5,1] <- "Unclassified_Bacteria"
Phylum_Therm[42,1] <- "Unclassified_Archaea"

# Read in Phylum Counts of genomes from (hyper)thermophiles in IMG 
JGI_Therm <- read.table("JGI_Thermophile_Phylum_Counts_November9th.txt", sep="\t", stringsAsFactors = FALSE, header = FALSE)
JGI_Therm$V3 <- JGI_Therm$V2/sum(JGI_Therm$V2)

#Make combined Dataset
Combined <- data.frame(matrix(NA, nrow = length(Phyla_Therm)+length(Phyla_Therm)+nrow(JGI_Therm), ncol = 3))
Combined[,1] <- c(Phylum_Therm$Phylum, Phylum_Therm$Phylum, JGI_Therm[,1])
Combined[,2] <- c(Phylum_Therm$Percentage_Therm, Phylum_Therm$PercentageOTUs, JGI_Therm[,3])
Combined[,3]<-c(rep("RelativeAbundanceTherm",42),rep("PercentageThermophileOTUs",42),rep("JGI",28))
Combined[,4]<-c(rep("A",42),rep("B",42),rep("C",28))

colnames(Combined) <- c("Phylum", "Proportion", "Category", "Text")
Combined$Category_f <- factor(Combined$Text, levels=c("A","B","C"))

p2 <- ggplot(data=Combined, aes(x=Phylum, y=Proportion, fill=Phylum)) + geom_bar(stat="identity",aes(x=Phylum,y=Proportion)) + facet_wrap(~Category_f, nrow=3,) + theme(strip.background = element_blank(), strip.text.x = element_blank(),axis.text.x = element_text(angle=90, hjust=1, vjust=.5,  face="bold", lineheight = rel(2), size=rel(1.2)), axis.title.x=element_text(vjust=-1), strip.text.x=element_text(size=rel(2))) + guides(fill=FALSE)
p2
ggsave("../Figures/Fig3.eps", width=200, units="mm")

### Venn Meta_Analysis
map.f <- read.table("Centralia_Full_Map_Fixed.txt", sep="\t", header= TRUE, row.names=1, stringsAsFactors = FALSE)
meta <- read.table("Carini_RDP_rmCM.txt", sep="\t", header = TRUE, row.names=1, stringsAsFactors = FALSE)

Cen <- meta[,grepl("C",colnames(meta))]
Car <- meta[,grepl("SRR", colnames(meta))]

Samples<- unique(map.f$Sample)

collapse <- NULL

for (i in 1:length(Samples)){
  x <- Cen[,grepl(Samples[i], colnames(Cen))]
  collapse <- cbind(collapse,rowSums(x))
}
colnames(collapse) <- Samples

Reference <- collapse[,map$Classification=="Reference"]
FireAffected <- collapse[,map$Classification=="FireAffected"]
Recovered <- collapse[,map$Classification=="Recovered"]

sum(Reference)
sum(FireAffected)
sum(Recovered)
sum(Car)
data.f <-cbind(rowSums(FireAffected),rowSums(Recovered), rowSums(Reference), rowSums(Car))
colnames(data.f) <- c("FireAffected", "Recovered", "Reference", "Carini")
library(vegan)

library(limma)
?rrarefy
data.rare<-rrarefy(t(data.f), 900542)
data.rare <- t(data.rare)
colSums(data.rare)


data.rare.pa <- 1*(data.rare>0)
data.rare.pa.nz <- data.rare.pa[rowSums(data.rare.pa)>0,]
v=vennCounts(data.rare.pa.nz)
v2=round(v[,"Counts"]/sum(v[,"Counts"]),2)
v[,"Counts"]<-v2
vennDiagram(v)

dev.off()
setEPS()
postscript("../Figures/Fig4.eps", width = 4.000, height=4.000, pointsize=8,paper="special")
vennDiagram(v)
dev.off()
### End 16S Analysis


### PCoA Relativized to RPOB WIP
library(calibrate)
class <- rep("black", nrow(map_MG))
class[map_MG$Classification=="Recovered"]='yellow'
class[map_MG$Classification=="FireAffected"]='red'
class[map_MG$Classification=="Reference"]='green'

# KEGG Orthology distanc matrix and PCoA
K.dist <- vegdist(t(KO.rpob), method="bray" )
K.dist
k.pcoa <-cmdscale(K.dist, eig=TRUE)

ax1.v=k.pcoa$eig[1]/sum(k.pcoa$eig)
ax2.v=k.pcoa$eig[2]/sum(k.pcoa$eig)
# Read in Weighted UniFrac Distance MAtrix and Calculate PCoA
uf=read.table("weighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)

#sort by rows, columns (so they are in the consecutive order)
uf=uf[order(row.names(uf)),order(colnames(uf))]
uf.d=as.dist(uf)

uf.pcoa=cmdscale(uf.d, eig=TRUE)
#calculate percent variance explained, then add to plot
ax1.v.uf=uf.pcoa$eig[1]/sum(uf.pcoa$eig)
ax2.v.uf=uf.pcoa$eig[2]/sum(uf.pcoa$eig)


class.uf <- rep("black", nrow(map))
class.uf[map$Classification=="Recovered"]='yellow'
class.uf[map$Classification=="FireAffected"]='red'
class.uf[map$Classification=="Reference"]='green'
?plot
dev.off()
setEPS()
postscript("../Figures/Figure5.eps", width = 8, height=4, pointsize=8,paper="special")
par(mfrow=c(1,2))
plot(k.pcoa$points[,1], k.pcoa$points[,2],cex=1.5, bg=class, pch=21, main= "rpoB Relativized Bray Curtis KEGG PCoA", xlab= paste("PCoA1: ",100*round(ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(ax2.v,3),"% var. explained",sep=""))
textxy(X=k.pcoa$points[,1], Y=k.pcoa$points[,2],labs=map_MG$Sample, cex=1)

env=map_MG[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]
envEF=envfit(k.pcoa, env)
plot(envEF, p.max=0.05, col="black")

plot(uf.pcoa$points[,1],uf.pcoa$points[,2] ,cex=1.5,pch=21,bg=class.uf,main="Weighted UniFrac PCoA",xlab= paste("PCoA1: ",100*round(ax1.v.uf,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(ax2.v.uf,3),"% var. explained",sep=""))
#textxy is from the calibrate library
textxy(X=uf.pcoa$points[,1], Y=uf.pcoa$points[,2],labs=map$Sample, cex=1)
legend('bottomleft',c('Fire Affected','Recovered','Reference'),pch=21,pt.bg=c("red", "yellow", "green"),lty=0)

#Add env vectors to plot that are significant p < 0.10
env.uf=map[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]
envEF.uf=envfit(uf.pcoa, env.uf)
plot(envEF.uf, p.max=0.05, col="black")
dev.off()

Class2=sub("green", "yellow", class)
### rpoB relativized hypothesis testing
a=adonis(K.dist~Class2, distance=TRUE, permutations=1000)
a
### rpoB relativized
b=betadisper(K.dist, group=Class2)
TukeyHSD(b, which = "group", ordered = FALSE,conf.level = 0.95)


### Using Coefficient of Variation for KEGG Orthologs to identify "responding" genes.
# Reference for single copy genes http://www.nature.com/articles/srep14840
KO <- as.matrix(KO)
KO.rpob <- as.matrix(KO.rpob)

KO_NZ <- KO[rowSums(KO)>0,]
KO.rpob_NZ <- KO.rpob[rowSums(KO.rpob)>0,]

CV<- NULL
for(i in 1:nrow(KO_NZ)){
  CV <- c(CV, sd(KO_NZ[i,])/mean(KO_NZ[i,]))
}

plot(CV)

CV.rpob <- NULL
for(i in 1:nrow(KO.rpob_NZ)){
  CV.rpob <- c(CV.rpob, sd(KO.rpob_NZ[i,])/mean(KO.rpob_NZ[i,]))
}

hist(CV)
hist(CV.rpob)

#Single Copy Genes
SCG <- read.table("Single_Copy_Genes.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

SCG_CV <- NULL 
for(i in 1:nrow(SCG)){
  SCG_CV <- c(SCG_CV,CV[grep(SCG[i,3], rownames(KO_NZ))])
}
SCG <- cbind(SCG, SCG_CV)

SCG_CV.rpob <- NULL 
for(i in 1:nrow(SCG)){
  SCG_CV.rpob <- c(SCG_CV.rpob,CV.rpob[grep(SCG[i,3], rownames(KO.rpob_NZ))])
}
SCG <- cbind(SCG, SCG_CV.rpob)

SCG_Order <- NULL
for(i in 1:nrow(SCG)){
  SCG_Order <- c(SCG_Order,grep(SCG$KO[i],rownames(KO_NZ)))
}
SCG_BaKO<- SCG[SCG$Class=="bacteria",]
SCG_BoKO <- SCG[SCG$Class=="both",]
SCG_ArKO<- SCG[SCG$Class=="archaea",]

SCG_Bacteria <- NULL
for(i in 1:nrow(SCG_BaKO)){
  SCG_Bacteria <- c(SCG_Bacteria, grep(SCG_BaKO$KO[i],rownames(KO_NZ)))
}

SCG_Both <- NULL
for(i in 1:nrow(SCG_BoKO)){
  SCG_Both <- c(SCG_Both, grep(SCG_BoKO$KO[i],rownames(KO_NZ)))
}

SCG_Archaea <- NULL
for(i in 1:nrow(SCG_ArKO)){
  SCG_Archaea <- c(SCG_Archaea, grep(SCG_ArKO$KO[i],rownames(KO_NZ)))
}

Weights <- rep(.5, length(CV))
Weights[SCG_Order]=5
colors <- rep('black',length(CV))
colors[SCG_Bacteria]="blue"
colors[SCG_Both]="green"
colors[SCG_Archaea]="red"
dev.off()
setEPS()
postscript("../Figures/Figure6.eps", width = 8, height=3, pointsize=8,paper="special")
par(mfrow=c(1,2))
plot(CV, col=colors,lwd=Weights)
hist(CV, prob=TRUE)
lines(density(CV))
abline(v = mean(SCG_CV[SCG$Class=="bacteria"]), col = "blue", lwd = 2)
abline(v= mean(SCG_CV[SCG$Class=="archaea"]), col="red", lwd=2)
abline(v= mean(SCG_CV[SCG$Class=="both"]), col="green", lwd=2)
dev.off()

setEPS()
postscript("../Figures/Figure6.eps", width = 8, height=3, pointsize=8,paper="special")
hist(CV.rpob, prob=TRUE)
lines(density(CV.rpob))
abline(v = mean(SCG_CV.rpob[SCG$Class=="bacteria"]), col = "blue", lwd = 2)
abline(v= mean(SCG_CV.rpob[SCG$Class=="archaea"]), col="red", lwd=2)
abline(v= mean(SCG_CV.rpob[SCG$Class=="both"]), col="green", lwd=2)
dev.off()
# Find KOs with CV greater than archaeal Single Copy Housekeeping Genes
High_CVA <- KO_NZ[CV.rpob>mean(SCG_CV.rpob[SCG$Class=="archaea"]),]
High_CVB <- KO_NZ[CV.rpob>mean(SCG_CV.rpob[SCG$Class=="bacteria"]),]
row.names(KO_NZ) <- sub("KO:", "", row.names(KO_NZ))

row.names(High_CVA) <- sub("KO:", "", row.names(High_CVA))
write.table(x=row.names(High_CVA), file="High_CVA_KO.txt", sep="\t")

row.names(High_CVB) <- sub("KO:", "", row.names(High_CVB))
write.table(x=row.names(High_CVB), file="High_CVB_KO.txt", sep="\t")

#KO Responses (using z-score across samples to look for similar patterns) CV greater than Archaeal Single Copy Genes

High_CVA.zs <- decostand(High_CVA, method="standardize", MARGIN=1)
High_CVA.zs <- as.matrix(High_CVA.zs)

#Create Dendrogram
HCVA_Den <- hclust(dist(High_CVA.zs))

#create color pallette; see: http://colorbrewer2.org/ 
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

Figure7A<-heatmap.2(High_CVA.zs,col=hc(100),scale="row",key=TRUE,symkey=FALSE,main="Archael Cutoff", trace="none", density.info="none",dendrogram="both", labRow=row.names(High_CVA.zs), colCol=c(class), margins=c(5,13), srtCol=90)

# Ordering Columns by temperature
High_CVA.zs.temp <- High_CVA.zs[,order(map_MG$SoilTemperature_to10cm)]


Figure7B<-heatmap.2(High_CVA.zs.temp,col=hc(100),scale="row",key=TRUE,symkey=FALSE, Colv= FALSE, trace="none", density.info="none",dendrogram="row", main="Archael Cutoff Ordered by Temp", labRow=row.names(High_CVA.zs), colCol=c("green","yellow","yellow","yellow","yellow","yellow","red","red","red","red","red","red"), margins=c(5,13), srtCol=90)


#KO Responses (using z-score across samples to look for similar patterns) CV greater than Bacterial Single Copy Genes

High_CVB.zs <- decostand(High_CVB, method="standardize", MARGIN=1)
High_CVB.zs <- as.matrix(High_CVB.zs)

#Creat Dendrogram
HCVB_Den <- hclust(dist(High_CVB.zs))

#create color pallette; see: http://colorbrewer2.org/ 
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

Figure8A<-heatmap.2(High_CVB.zs,col=hc(100),scale="row",key=TRUE,symkey=FALSE, trace="none", main="Bacterial Cutoff", density.info="none",dendrogram="column", labRow=row.names(High_CVB.zs), colCol=c(class), margins=c(5,13), srtCol=90,)
Figure8A<-heatmap.2(High_CVB.zs,col=hc(100),scale="row",key=TRUE,symkey=FALSE, trace="none", main="Bacterial Cutoff", density.info="none",dendrogram="both", labRow=row.names(High_CVB.zs), colCol=c(class), margins=c(5,13), srtCol=90,)


# Ordering Columns by temperature
High_CVB.zs.temp <- High_CVB.zs[,order(map_MG$SoilTemperature_to10cm)]


Figure8B<-heatmap.2(High_CV.zs.temp,col=hc(100),scale="row",key=TRUE,symkey=FALSE, Colv= FALSE, trace="none", main="Bacterial Cutoff Ordered by Temp", density.info="none",dendrogram="row", labRow=row.names(High_CVB.zs), colCol=c("green","yellow","yellow","yellow","yellow","yellow","red","red","red","red","red","red"), margins=c(5,13), srtCol=90)



### How to define meaningful clusters of KEGG Orthologs based on abundance patterns across the gradient
HCVA_Den <- hclust(dist(High_CVA.zs))
plot(HCVA_Den)
rect.hclust(HCVA_Den, h=5.5)
cutree(HCVA_Den, h=5.5)
Output <- NULL
for(i in 1:13){
  cluster_Cor <- NULL
  cluster <- High_CVA[cutree(HCVA_Den, h=5.5)==i,]
  for (n in 1:nrow(cluster)){
    results<-cor.test(cluster[n,],apply(cluster, 2, median))
    cluster_Cor<- rbind(cluster_Cor, c(results$estimate,results$p.value))
    row.names(cluster_Cor[n,]) <- row.names(cluster[n,])
  }
  Output <- rbind(Output, c(apply(cluster_Cor, 2, mean), nrow(cluster_Cor)))
}

High_CVA_Clusters <- NULL
High_CVA_Clusters <- as.list(High_CVA_Clusters)
for (i in 1:length(unique(cutree(HCVA_Den, h=5.5)))){
  High_CVA_Clusters[[length(High_CVA_Clusters)+1]] <- as.data.frame(High_CVA[cutree(HCVA_Den, h=5.5)==i,])
}

write.table(row.names(High_CVA_Clusters[[1]]),"HCVA_H55_C1.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[2]]),"HCVA_H55_C2.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[3]]),"HCVA_H55_C3.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[4]]),"HCVA_H55_C4.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[5]]),"HCVA_H55_C5.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[6]]),"HCVA_H55_C6.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[7]]),"HCVA_H55_C7.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[8]]),"HCVA_H55_C8.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[9]]),"HCVA_H55_C9.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[10]]),"HCVA_H55_C10.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[11]]),"HCVA_H55_C11.txt", sep="\t",quote=FALSE, col.names=FALSE)

write.table(row.names(High_CVA_Clusters[[12]]),"HCVA_H55_C12.txt", sep="\t",quote=FALSE, col.names=FALSE)


### KO T-test (difference in functional richness)

KO_NZ_PA <- 1*(KO_NZ>0)
colSums(KO_NZ_PA)


class3 <- class
class3
class3[12] <- "yellow"
t.test(colSums(KO_NZ_PA[,class3=="yellow"]),colSums(KO_NZ_PA[,class3=="red"]))
# No Significant difference in richness of KEGG Orthologs found between Recovered/Reference and FireAffected





plot(hclust(dist(t(KO.rpob_NZ))))

### t test of KO abudnace in recovered vs reference sites
output.out<- NULL
for(i in 1:nrow(KO.rpob_NZ)){
  active <- KO.rpob_NZ[i,class3=="red"]
  recref <- KO.rpob_NZ[i,class3=="yellow"]
  output <- t.test(active,recref)
  output.out <- rbind(output.out, c(rownames(KO.rpob_NZ)[i], output$statistic, output$parameter, output$p.value))
}

hist(output.out[,3])
sig_KO <- output.out[output.out[,4]<0.05,]

KO_FA_SIG.rpob <- KO.rpob_NZ[output.out[,4]<0.05,]

FA_KOs <- sig_KO[sig_KO[,2]>0,]
KO_FA_SIG.rpob <- KO_FA_SIG.rpob[sig_KO[,2]>0,]

RR_KOs <- sig_KO[sig_KO[,2]<0,]

row.names(KO_FA_SIG.rpob) <- sub("KO:", "", row.names(KO_FA_SIG.rpob))

FA_KOs[,1] <- sub("KO:", "", FA_KOs[,1])
RR_KOs[,1] <- sub("KO:", "", RR_KOs[,1])

FAKO.zs <- decostand(KO_FA_SIG.rpob, method="standardize", MARGIN=1)
FAKO.zs.temp <- FAKO.zs[,order(map_MG$SoilTemperature_to10cm)]

Figure9<-heatmap.2(FAKO.zs.temp,col=hc(100),scale="row",key=TRUE,symkey=FALSE, Colv= FALSE, trace="none", main="Fire Affected KOs", density.info="none",dendrogram="row", labRow=row.names(FAKO.zs.temp), colCol=c("green","yellow","yellow","yellow","yellow","yellow","red","red","red","red","red","red"), margins=c(5,13), srtCol=90)


# Outputs for MinPath (http://omics.informatics.indiana.edu/MinPath/run.php)
write.table(FA_KOs[,1],"Significant_FireAffected_KOs.txt", sep="\t",quote=FALSE, col.names=FALSE)
write.table(RR_KOs[,1], "Significant_RecRef_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)


### What about all of our KOs
KO.rpob_NZ
row.names(KO.rpob_NZ)<-sub("KO:", "", row.names(KO.rpob_NZ))

write.table(row.names(KO.rpob_NZ),"All_KOs.txt", sep="\t",quote=FALSE, col.names=FALSE)


MinPath_AllKOs <- read.table("MinPath_AllKOs_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,1]>0,]), "C01_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,2]>0,]), "C03_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,3]>0,]), "C04_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,4]>0,]), "C05_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,5]>0,]), "C06_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)


C01 <- 
write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,6]>0,]), "C07_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,7]>0,]), "C10_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,8]>0,]), "C12_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,9]>0,]), "C14_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,10]>0,]), "C15_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,11]>0,]), "C16_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)

write.table(row.names(KO.rpob_NZ[KO.rpob_NZ[,12]>0,]), "C17_KOs.txt", sep="\t", quote=FALSE, col.names=FALSE)


C01 <- read.table("C01_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C03 <- read.table("C03_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C04 <- read.table("C04_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C05 <- read.table("C05_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C06 <- read.table("C06_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C07 <- read.table("C07_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C10 <- read.table("C10_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C12 <- read.table("C12_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C14 <- read.table("C14_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C15 <- read.table("C15_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C16 <- read.table("C16_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
C17 <- read.table("C17_MinPath_Formatted.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

#Shared_Paths <- unique(c(C01$Path, C03$Path, C04$Path, C05$Path, C06$Path, C07$Path, C10$Path, C12$Path, C14$Path, C15$Path, C16$Path, C17$Path))
Shared_Paths <- Reduce(intersect, list(C01$Path, C03$Path, C04$Path, C05$Path, C06$Path, C07$Path, C10$Path, C12$Path, C14$Path, C15$Path, C16$Path, C17$Path))

grep(Shared_Paths[1],C01$Path)
MinPath_Data <- NULL
for (i in 1:length(Shared_Paths)){
  x <- c(C01[C01$Path==Shared_Paths[i],4],C03[C03$Path==Shared_Paths[i],4],C04[C04$Path==Shared_Paths[i],4],C05[C05$Path==Shared_Paths[i],4],C06[C06$Path==Shared_Paths[i],4],C07[C07$Path==Shared_Paths[i],4],C10[C10$Path==Shared_Paths[i],4],C12[C12$Path==Shared_Paths[i],4],C14[C14$Path==Shared_Paths[i],4],C15[C15$Path==Shared_Paths[i],4],C16[C16$Path==Shared_Paths[i],4],C17[C17$Path==Shared_Paths[i],4])
  MinPath_Data <- rbind(MinPath_Data, x)
}



MinPath_RecRef <- MinPath_Data[,map_MG$Classification!="FireAffected"]
MinPath_FA <- MinPath_Data[,map_MG$Classification=="FireAffected"]

colnames(MinPath_RecRef) <- map_MG[map_MG$Classification!="FireAffected",][,1]
colnames(MinPath_FA) <- map_MG[map_MG$Classification=="FireAffected",][,1]
MinPath_tTest <- NULL
for (i in 1:nrow(MinPath_Data)){
  result <- t.test(MinPath_RecRef[i,],MinPath_FA[i,])
  MinPath_tTestt <- rbind(MinPath_tTest, c(Shared_Paths[i],result$statistic, result$parameter, result$p.value))
}

plot(apply(MinPath_RecRef, 1, mean), apply(MinPath_FA,1,mean))
setdiff(Shared_Paths, Pathways)
