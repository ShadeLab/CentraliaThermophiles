## Workflow Script for Sorensen_INPREP
library(ggplot2)
library(vegan)
library(outliers)
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

# Calculating Coefficient of Variation:
http://www.nature.com/articles/srep14840
KO <- as.matrix(KO)

KO_NZ <- KO[rowSums(KO)>0,]
CV<- NULL

for(i in 1:nrow(KO_NZ)){
  CV <- c(CV, sd(KO_NZ[i,])/mean(KO_NZ[i,]))
}
warnings()
plot(CV)


CV[grep("K03043", rownames(KO_NZ),)]

#Single Copy Genes
SCG <- read.table("Single_Copy_Genes.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

SCG_CV <- NULL 
for(i in 1:nrow(SCG)){
  SCG_CV <- c(SCG_CV,CV[grep(SCG[i,3], rownames(KO_NZ),)])
}
SCG <- cbind(SCG, SCG_CV)

SCG_Order <- NULL
for(i in 1:nrow(SCG)){
  SCG_Order <- c(SCG_Order,grep(SCG$KO[i],rownames(KO_NZ)))
}

colors <- rep('black',length(CV))
colors[SCG_Order]='red'
plot(CV, col=colors)

#KO Responses (using z-score across samples to look for similar patterns)
comm.phylum.oc=decostand(comm.phylum, method="standardize", margin=1)
comm.phylum.oc=as.matrix(comm.phylum.oc)

KO.rpob.zs <- decostand(KO.rpob, method="standardize", margin=1)
KO.rpob.zs <- as.matrix(KO.rpob.zs)
