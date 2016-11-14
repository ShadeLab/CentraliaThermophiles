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

hist(CV.rpob, prob=TRUE)
lines(density(CV.rpob))
abline(v = mean(SCG_CV.rpob[SCG$Class=="bacteria"]), col = "blue", lwd = 2)
abline(v= mean(SCG_CV.rpob[SCG$Class=="archaea"]), col="red", lwd=2)
abline(v= mean(SCG_CV.rpob[SCG$Class=="both"]), col="green", lwd=2)

# Find KOs with CV greater than archaeal Single Copy Housekeeping Genes
High_CV <- KO_NZ[CV.rpob>mean(SCG_CV.rpob[SCG$Class=="archaea"]),]
row.names(KO_NZ) <- sub("KO:", "", row.names(KO_NZ))

row.names(High_CV) <- sub("KO:", "", row.names(High_CV))
write.table(x=row.names(High_CV), file="High_CV_KO.txt", sep="\t")

#KO Responses (using z-score across samples to look for similar patterns)
comm.phylum.oc=decostand(comm.phylum, method="standardize", margin=1)
comm.phylum.oc=as.matrix(comm.phylum.oc)

High_CV.zs <- decostand(High_CV, method="standardize", MARGIN=1)
High_CV.zs <- as.matrix(High_CV.zs)

#create color pallette; see: http://colorbrewer2.org/ 
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

Figure7<-heatmap.2(High_CV.zs,col=hc(100),scale="row",key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", labRow=row.names(High_CV.zs), colCol=c(class), margins=c(5,13), srtCol=90)


# Ordering Columns by temperature
High_CV.zs.temp <- High_CV.zs[,order(map_MG$SoilTemperature_to10cm)]
comm=comm[,order(colnames(comm))]

Figure7B<-heatmap.2(High_CV.zs.temp,col=hc(100),scale="row",key=TRUE,symkey=FALSE, Colv= FALSE, trace="none", density.info="none",dendrogram="row", labRow=row.names(High_CV.zs), colCol=c("green","yellow","yellow","yellow","yellow","yellow","red","red","red","red","red","red"), margins=c(5,13), srtCol=90)


### How to define meaningful clusters of KEGG Orthologs based on abundance patterns across the gradient
HCV_Den <- hclust(dist(High_CV.zs))
rect.hclust(HCV_Den, k=13)
cutree(HCV_Den, k=13)
Output <- NULL
for(i in 1:13){
  cluster_Cor <- NULL
  cluster <- High_CV.zs[cutree(HCV_Den, k=13)==i,]
  for (n in 1:nrow(cluster)){
    results<-cor.test(cluster[n,],apply(cluster, 2, median))
    cluster_Cor<- rbind(cluster_Cor, c(results$estimate,results$p.value))
    row.names(cluster_Cor[n,]) <- row.names(cluster[n,])
  }
  Output <- rbind(Output, c(apply(cluster_Cor, 2, mean), nrow(cluster_Cor)))
}

### KO T-test (difference in functional richness)

KO_NZ_PA <- 1*(KO_NZ>0)
colSums(KO_NZ_PA)


class3 <- class
class3
class3[12] <- "yellow"
t.test(colSums(KO_NZ_PA[,class3=="yellow"]),colSums(KO_NZ_PA[,class3=="red"]))
# No Significant difference in richness of KEGG Orthologs found between Recovered/Reference and FireAffected
