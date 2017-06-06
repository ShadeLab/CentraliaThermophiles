## Workflow Script for Sorensen_INPREP
library(ggplot2)
library(vegan)
library(outliers)
library(gplots)
library(colorRamps)
library(dplyr)
### Reading in Data Files and Manipulating Datafile
setwd("~/GitHub_Repos/ShadeLab/CentraliaThermophiles/Workflow/")
# Mapping File
map <- read.table("Centralia_Collapsed_Map_forR.txt",sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)

map_MG <- map[c(1,3,4,5,6,7,10,12,14,15,16,17),]

sum_stats <- read.table("../Centralia_Metagenome_Summary_Stats.txt", sep="\t", header=TRUE, row.names=1)
sum_stats <- sum_stats[,-10]
sum_stats$PercentMapped <- sum_stats$Aligned.Reads/sum_stats$Quality.Reads

sum_stats$Sample <- map_MG$Sample

CO2_2015 <- c(525,455,500,512,484,779,581,503,6094,17112,15760,13969,7390,723,1624,597,480,477)
CO2_MG_2015 <- CO2_2015[c(1,3,4,5,6,7,10,12,14,15,16,17)]
t.test(CO2_MG_2015[map_MG$Classification=="FireAffected"], CO2_MG_2015[map_MG$Classification!="FireAffected"])

# Alpha Diversity
alpha <- read.table("MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_alphadiv.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names = 1)

alpha <- alpha[order(row.names(alpha)),]
alpha_MG <- alpha[c(1,3,4,5,6,7,10,12,14,15,16,17),]

setEPS()
postscript("../Figures/MGStatsVs16Sstats.eps", width = 5, height=5, pointsize=8,paper="special")
par(mfrow=c(2,2))
fig = plot(alpha_MG$PD_whole_tree, sum_stats$PercentMapped, ylab="Percent Reads Mapped", xlab="PD") 
plot(alpha_MG$observed_otus, sum_stats$PercentMapped, ylab="Percent Reads Mapped", xlab="OTUs")
plot(alpha_MG$PD_whole_tree, sum_stats$Assembled.Length, ylab="Assembled Length", xlab="PD")
plot(alpha_MG$observed_otus, sum_stats$Assembled.Length, ylab="Assembled Length", xlab="OTUs")
dev.off()
cor.test(alpha_MG$PD_whole_tree, sum_stats$PercentMapped)
cor.test(alpha_MG$observed_otus, sum_stats$PercentMapped)
cor.test(alpha_MG$PD_whole_tree, sum_stats$Assembled.Length)
cor.test(alpha_MG$observed_otus, sum_stats$Assembled.Length)

cor.test(sum_stats$Quality.Reads, sum_stats$Assembled.Length)
plot(sum_stats$Quality.Reads, sum_stats$Assembled.Length)
cor.test(alpha_MG$observed_otus, sum_stats$Quality.Reads)
library(outliers)

grubbs.test(sum_stats$Quality.Reads, type=10)
grubbs.test(sum_stats$Quality.Reads, type=11)
grubbs.test(sum_stats$Quality.Reads, type=20)

grubbs.test(sum_stats$Assembled.Length, type=10)
grubbs.test(sum_stats$Assembled.Length, type=11)
grubbs.test(sum_stats$Assembled.Length, type=20)


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
library(vegan)
## Assembled and Unassembled
KO <- read.table("4-13-17abundance_ko_126107.tab.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)
# Assembled Only
##KO <- read.table("KO_minus2col_09-27-2016.tab.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
colnames(KO) <- map_MG$Sample
row.names(KO) <- sub("KO:", "", row.names(KO))
KO<- KO[rowSums(KO)>0,]

# So whenever in R, make sure its not accidentally trying to add strings...
#rpoB <- KO["KO:K03043",] + KO["KO:K13798",] ### add the bacterial and archaeal copies together.
rpoB <- as.numeric(KO["K03043",]) + as.numeric(KO["K13798",])
# informative Plots
 plot(rpoB, map_MG$SoilTemperature_to10cm)
 plot(as.numeric( KO["K03043",]), map_MG$SoilTemperature_to10cm)
 plot(as.numeric(KO["K13798",]), map_MG$SoilTemperature_to10cm)
 plot(as.numeric(KO["K03043",]), as.numeric(KO["K13798",]))
#
# Remove KOs with zeroes across the entire dataset
#Rarefy for Unassembled and Assembled
KO.rare <- t(rrarefy(t(KO),11855473))
#Rarefy for Assembled Only
#KO <- t(rrarefy(t(KO),10064577))
rpoB <- as.numeric(KO["K03043",]) + as.numeric(KO["K13798",])
#Relativze to rpoB
#KO.rpob <- NULL
#for (i in 1:nrow(KO)){
#  KO.rpob <- rbind(KO.rpob, KO[i,]/rpoB)
#}

#row.names(KO.rpob) <- sub("KO:", "", row.names(KO))

# Relativize to rpsM
#KO.rpsM <- NULL
#for(i in 1:nrow(KO)){
#  KO.rpsM<- rbind(KO.rpsM, KO[i,]/KO["K02952",])
#}
#row.names(KO.rpsM) <- row.names(KO.rpob)
#Relativized to KO count
KO.rel <- decostand(x=KO, method="total", MARGIN=2)
row.names(KO.rel) <- sub("KO:","",row.names(KO.rel))
#Single Copy KO's
library(readr)
COG_Key <- read_delim("~/GitHub_Repos/ShadeLab/CentraliaThermophiles/Workflow/Supplemental/He_et_al_COG_to_KEGG.txt","\t", escape_double = FALSE, trim_ws = TRUE)
COG_Key <- as.data.frame(COG_Key)
sc <- COG_Key[,3]
SCG_Rel <- KO.rel[sc,] 
SCG_Rel <- SCG_Rel[complete.cases(SCG_Rel),]
SCG_Absolute <- KO[sc,]
SCG_Absolute <- SCG_Absolute[complete.cases(SCG_Absolute),]

# For use in Odds Ratios
Average_MG_SCG <- apply(SCG_Rel, 1, mean)
#For use in Normalizing KEGG Data
Average_SCG <- apply(SCG_Absolute, 2, median)

# KOs relatvized to average Single Copy Gene abundance

KO.sr<- NULL
for(i in 1:nrow(KO)){
  KO.sr <- rbind(KO.sr, KO[i,]/Average_SCG)
}

KO <- as.data.frame(KO)
Odds_Ratio <- NULL
for (i in 1:nrow(SCG_Rel)){
  Odds_Ratio <- rbind(Odds_Ratio, KO.rel[row.names(SCG_Rel[i,]),]/Average_MG_SCG[i])
}

Odds_Ratio <- Odds_Ratio[complete.cases(Odds_Ratio),]

SCG_Correlations <- NULL
for(i in 1:nrow(Odds_Ratio)){
  result <- cor.test(map_MG$SoilTemperature_to10cm, as.numeric(Odds_Ratio[i,]))
  SCG_Correlations <- rbind(SCG_Correlations,unlist(result[1:4]))
}
SCG_Correlations <- as.data.frame(SCG_Correlations)
SCG_Correlations$Adjusted.p.value <- p.adjust(SCG_Correlations$p.value, method="fdr")
SCG_SigCorrelations <- SCG_Correlations[SCG_Correlations$Adjusted.p.value<0.05,]

row.names(SCG_Correlations) <- row.names(Odds_Ratio)
par(mfrow=c(1,1))
SCG_Correlations <- as.data.frame(SCG_Correlations)
SCG_Correlations <- SCG_Correlations[order(row.names(SCG_Correlations)),]
SCG_Correlations$KEGG <- row.names(SCG_Correlations)
New_Table <- inner_join(COG_Key, SCG_Correlations, "KEGG")

### ggplots Odds_Ratio
library(reshape2)
Odds_Ratio


melt(Odds_Ratio, id.vars=colnames(Odds_Ratio), )

MO <- cbind(Odds_Ratio, row.names(Odds_Ratio))
MO <- melt(MO, id.vars="row.names(Odds_Ratio)",variable.name="Sample", value.name="Measurement")
colnames(MO) <- c("KO", "Sample", "Measurement")
library(dplyr)
?inner_join

fuller_map <- inner_join(map_MG, sum_stats, by="Sample")
Joined_Data <- inner_join(MO, fuller_map, by="Sample")
library(ggplot2)
ggplot(Joined_Data, aes(x=SoilTemperature_to10cm, y=PercentMapped, color=KO)) + geom_point() +geom_smooth(method="lm", alpha=0)
odr <- ggplot(Joined_Data, aes(x=SoilTemperature_to10cm, y=Measurement, color=KO)) + geom_point() + geom_smooth(method="lm", alpha=0) + guides(color=FALSE)

#odr <- ggplot(Joined_Data, aes(x=SoilTemperature_to10cm, y=Measurement, color=KO)) + geom_point() + geom_smooth(method="lm", alpha=0) + facet_wrap(~KO, nrow=6)
ggsave("Supplemental/SCG_OddsRation.png", odr)

### Venn Analysis 16S 
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("limma")
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

### Who belongs to each category
# f= FireAffected, r=Reference, c=Recovered
comm.rdp <- cbind(comm, rdp)

#FireAffected Only
f<- venndata[venndata[,1]>0,]
f<- f[f[,2]==0,]
f <- f[f[,3]==0,]
FireOnly <- comm.rdp[row.names(f),]
f_tax <- c(1:7)
for (i in 1:nrow(FireOnly)){
  f_tax <- rbind(f_tax,unlist(strsplit(as.character(FireOnly[i,19]), ";")))
}
f_tax <- f_tax[-1,]
# FireAffected + Recovered
fc <- venndata[venndata[,1]>0,]
fc <- fc[fc[,2]>0,]
fc <- fc[fc[,3]==0,]
Fire_Recovered <- comm.rdp[row.names(fc),]
fc_tax <- c(1:7)
for (i in 1:nrow(Fire_Recovered)){
  fc_tax <- rbind(fc_tax,unlist(strsplit(as.character(Fire_Recovered[i,19]), ";")))
}
fc_tax <- fc_tax[-1,]
# FireAffected + Reference
fr <- venndata[venndata[,1]>0,]
fr <- fr[fr[,3]>0,]
fr <- fr[fr[,2]==0,]
Fire_Reference <- comm.rdp[row.names(fr),]
fr_tax <- c(1:7)
for (i in 1:nrow(Fire_Reference)){
  fr_tax <- rbind(fr_tax,unlist(strsplit(as.character(Fire_Reference[i,19]), ";")))
}
fr_tax <- fr_tax[-1,]
#Recovered Only
c <- venndata[venndata[,2]>0,]
c<- c[c[,1]==0,]
c <- c[c[,3]==0,]
RecoveredOnly <- comm.rdp[row.names(c),]
c_tax <- c(1:7)
for (i in 1:nrow(RecoveredOnly)){
  c_tax <- rbind(c_tax,unlist(strsplit(as.character(RecoveredOnly[i,19]), ";")))
}
c_tax <- c_tax[-1,]
#Recovered + Reference
cf <- venndata[venndata[,2]>0,]
cf<- cf[cf[,3]>0,]
cf <- cf[cf[,1]==0,]
Recovered_Reference <- comm.rdp[row.names(cf),]
cf_tax <- NULL
for (i in 1:nrow(Recovered_Reference)){
  cf_tax <- rbind(cf_tax,unlist(strsplit(as.character(Recovered_Reference[i,19]), ";")))
}
# Reference Only
r <- venndata[venndata[,3]>0,]
r <- r[r[,2]==0,]
r <- r[r[,1]==0,]
ReferenceOnly <- comm.rdp[row.names(r),]
r_tax <- NULL
for (i in 1:nrow(ReferenceOnly)){
  r_tax <- rbind(r_tax,unlist(strsplit(as.character(ReferenceOnly[i,19]), ";")))
}


# All Classes
a <- venndata[venndata[,1]>0,]
a <- a[a[,2]>0,]
a <- a[a[,3]>0,]
All_Classes <- comm.rdp[row.names(a),]
ac_tax <- NULL
for (i in 1:nrow(All_Classes)){
  ac_tax <- rbind(ac_tax,unlist(strsplit(as.character(All_Classes[i,19]), ";")))
}

Taxa <- NULL
for (i in 1:nrow(comm.rdp)){
  Taxa <- rbind(Taxa,unlist(strsplit(as.character(comm.rdp[i,19]), ";")))
}
Phyla <- unique(Taxa[,2])
Phyla <- gsub("\\[|\\]","", Phyla)
Phyla <- gsub("k__","", Phyla)
Phyla <- gsub(" p__","", Phyla)
Phyla <- Phyla[-57]
Phyla
Venn_Phyla_Counts <- NULL
Venn_Phyla_OTUs<- NULL
for (i in 1:length(Phyla)){
  tryCatch({Venn_Phyla_Counts <- rbind(Venn_Phyla_Counts,c(sum(Fire_Recovered[grepl(Phyla[i],fc_tax[,2]),1:18]),sum(FireOnly[grepl(Phyla[i],f_tax[,2]),1:18]),sum(All_Classes[grepl(Phyla[i],ac_tax[,2]),1:18]),sum(Fire_Reference[grepl(Phyla[i],fr_tax[,2]),1:18]),sum(Recovered_Reference[grepl(Phyla[i],cf_tax[,2]),1:18]),sum(RecoveredOnly[grepl(Phyla[i],c_tax[,2]),1:18]),sum(ReferenceOnly[grepl(Phyla[i],r_tax[,2]),1:18])))},error=function(e){})
  Venn_Phyla_OTUs <- rbind(Venn_Phyla_OTUs, c(nrow(Fire_Recovered[grepl(Phyla[i],fc_tax[,2]),]),nrow(FireOnly[grepl(Phyla[i],f_tax[,2]),]),nrow(All_Classes[grepl(Phyla[i],ac_tax[,2]),]),nrow(Fire_Reference[grepl(Phyla[i],fr_tax[,2]),]),nrow(Recovered_Reference[grepl(Phyla[i],cf_tax[,2]),]),nrow(RecoveredOnly[grepl(Phyla[i],c_tax[,2]),]),nrow(ReferenceOnly[grepl(Phyla[i],r_tax[,2]),])))
}
row.names(Venn_Phyla_OTUs) <- Phyla
colnames(Venn_Phyla_OTUs) <- c("FireRecovered", "FireOnly","AllClasses","FireReference","RecoveredReference","RecoveredOnly","ReferenceOnly")
#Write out the results of venncounts

write.table(v, "VennCounts.txt", quote=FALSE, sep="\t")
par(mfrow=c(1,1))
setEPS()

postscript("../Figures/Fig1.eps", width = 5, height=5, pointsize=8,paper="special")
par(mar=c(5,3,2,2)+0.1)
fig1=vennDiagram(v)
dev.off()


### KEGG Ortholog Eveness
s=specnumber(t(KO))
h=diversity(t(KO), index="shannon")
pielou = h/log(s)

### DeNovo Vs Reference OTU Analysis FINISHED
library(reshape2)
library(ggplot2)
Percent_DN <- colSums(comm[grepl("dn",row.names(comm)),])/colSums(comm)
map$PercentDeNovo <- Percent_DN

map.long=melt(map, id.vars=c("Sample", "SoilTemperature_to10cm", "Classification"), measure.vars=c("NO3N_ppm","NH4N_ppm","pH","SulfateSulfur_ppm","K_ppm","Ca_ppm","Mg_ppm","OrganicMatter_500","Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","PercentDeNovo"))

p <- ggplot(map, aes(SoilTemperature_to10cm, PercentDeNovo))
p+ geom_point(size=3) + labs(y="Community Novelty", x="Soil Temperature") + geom_smooth(method="lm", alpha=0,colour="black") + theme(axis.text= element_text(size=15), axis.ticks = element_line(size=1)) + ylim(0,0.7)
ggsave("../Figures/Fig2.eps", width=178, units="mm")


# De Novo Otus
dn_OTUs <- cbind(comm[grep("dn", row.names(comm)),], rdp[grep("dn", row.names(comm))])

dn_tax <- NULL
for (i in 1:nrow(dn_OTUs)){
  dn_tax <- rbind(dn_tax,unlist(strsplit(as.character(dn_OTUs[i,19]), ";")))
}

gg_OTUs <- cbind(comm[grep("dn", row.names(comm), invert=TRUE),], rdp[grep("dn", row.names(comm), invert=TRUE)])
gg_tax <- NULL
for (i in 1:nrow(gg_OTUs)){
  gg_tax <- rbind(gg_tax,unlist(strsplit(as.character(gg_OTUs[i,19]), ";")))
}

All_Phyla <- unique(c(gg_tax[,2], dn_tax[,2]))

All_Phyla <- sub("k__","",All_Phyla)
All_Phyla <- sub(" p__", "", All_Phyla)
All_Phyla <- gsub("\\[|\\]","", All_Phyla)
All_Phyla <- All_Phyla[-49]

dn_Phyla_Stats <- NULL
gg_Phyla_Stats <- NULL
for (i in 1:length(All_Phyla)){
  x <- tryCatch({sum(dn_OTUs[grep(All_Phyla[i], dn_tax[,2]),1:18])}, error=function(e){return(0)})
  y <- tryCatch({nrow(dn_OTUs[grep(All_Phyla[i], dn_tax[,2]),1:18])}, error=function(e){return(0)})
  a <- tryCatch({sum(gg_OTUs[grep(All_Phyla[i], gg_tax[,2]),1:18])}, error=function(e){return(0)})
  b <- tryCatch({nrow(gg_OTUs[grep(All_Phyla[i], gg_tax[,2]),1:18])}, error=function(e){return(0)})
  dn_Phyla_Stats <- rbind(dn_Phyla_Stats, c(x, y))
  gg_Phyla_Stats <- rbind(gg_Phyla_Stats, c(a, b))
}


Results <- NULL
Percent_Denovo_Phyla <- NULL
for(i in 1:length(All_Phyla)){
  x <- tryCatch({colSums(dn_OTUs[grep(All_Phyla[i], dn_tax[,2]),1:18])}, error=function(e){return(rep(0,18))})
  y <- tryCatch({colSums(gg_OTUs[grep(All_Phyla[i], gg_tax[,2]),1:18])}, error=function(e){return(rep(0,18))})
  z <- x/(x+y)
  Percent_Denovo_Phyla <- rbind(Percent_Denovo_Phyla, z)
  a <- tryCatch({cor.test(z, map$SoilTemperature_to10cm)}, error=function(e){cor.test(c(1,2,3,4), c(1,2,3,4))})
  Results <- rbind(Results, c(a$statistic, a$parameter, a$p.value, a$estimate))
}
Results <- as.data.frame(Results)
Results$Phylum <- All_Phyla
row.names(Percent_Denovo_Phyla) <- All_Phyla
Results<- Results[Results$cor!=1,]
Results$fdr <- p.adjust(Results$V3, method="fdr")
subset <- Percent_Denovo_Phyla[c("Acidobacteria","Proteobacteria","Verrucomicrobia"),]
subset <- rbind(subset, map$SoilTemperature_to10cm)
row.names(subset) <- c("Acidobacteria","Proteobacteria","Verrucomicrobia", "Soil_temperature")
melted_subset <- melt(t(subset))
colnames(melted_subset) <- c("Sample", "Phylum", "Percent_Denovo")
melted_subset$SoilTemperature <- rep(map$SoilTemperature_to10cm,3)


odr <- ggplot(Joined_Data, aes(x=SoilTemperature_to10cm, y=Measurement, color=KO)) + geom_point() + geom_smooth(method="lm", alpha=0) + guides(color=FALSE)
PPDN <- ggplot(melted_subset, aes(x=SoilTemperature, y=Percent_Denovo, shape=Phylum)) + geom_point(size=3) + geom_smooth(colour="black", method="lm", alpha=0, aes(linetype=factor(Phylum))) + scale_linetype_manual("", values=c(1,2,3)) + scale_shape_manual("", values=c(1,2,3)) +guides(linetype=FALSE, shape=FALSE) +ylim(0,0.7) +theme(axis.text= element_text(size=15), axis.ticks = element_line(size=1))
ggsave("../Figures/PPDN.eps", plot=PPDN, width=178, units="mm")
setEPS()
postscript("../Figures/Phyla.eps", width = 5.000, height=10.000, pointsize=10,paper="special")
par(mfrow=c(3,1))
plot( map$SoilTemperature_to10cm,Percent_Denovo_Phyla["Acidobacteria",], main="Acidobacteria", xaxt=NULL, xlab="", yaxt=NULL, ylab="")
plot(map$SoilTemperature_to10cm, Percent_Denovo_Phyla["Proteobacteria",], main="Proteobacteria", xaxt=NULL, xlab="",yaxt=NULL, ylab="")
plot(map$SoilTemperature_to10cm, Percent_Denovo_Phyla["Verrucomicrobia",], main="Verrucomicrobia",xaxt=NULL, xlab="",yaxt=NULL, ylab="")
dev.off()

dn_Phyla_Stats<- as.data.frame(dn_Phyla_Stats)
gg_Phyla_Stats <- as.data.frame(gg_Phyla_Stats)
Total_Phyla_Stats <- cbind(All_Phyla, dn_Phyla_Stats, gg_Phyla_Stats)
Total_Phyla_Stats[,1] <- as.character(Total_Phyla_Stats[,1])
colnames(Total_Phyla_Stats) <- c("Phylum", "DeNovoReads", "DeNovoOTUs", "GGReads", "GGOTUs")

TPS_Abundant <- Total_Phyla_Stats[(Total_Phyla_Stats[,3]+Total_Phyla_Stats[,5])>1000,]
TPS_Rare <- Total_Phyla_Stats[(Total_Phyla_Stats[,3]+Total_Phyla_Stats[,5])<1000,]

Other <- TPS_Rare[(TPS_Rare[,3]+TPS_Rare[,5])<50,]
Other <- c(0, colSums(Other[,2:5]))
TPS_Rare <- TPS_Rare[(TPS_Rare[,3]+TPS_Rare[,5])>50,]
TPS_Rare <- rbind(TPS_Rare, Other )
TPS_Rare[21,1]<- "Other" 

LD_A <- melt(TPS_Abundant, id.vars="Phylum")
LD_R <- melt(TPS_Rare, id.vars = "Phylum")
LD_A[,4] <- c(rep("Reads",8), rep("OTUs",8), rep("Reads",8), rep("OTUs",8))
LD_A[,2] <- c(rep("DeNovo",16), rep("GG", 16))

LD_R[,4] <- c(rep("Reads",21), rep("OTUs",21), rep("Reads",21), rep("OTUs",21))
LD_R[,2] <- c(rep("DeNovo",42), rep("GG", 42))

colnames(LD_A) <- c("Phylum", "Classification", "Count", "Type")
LD_A$Type <- factor(LD_A$Type, levels=c("OTUs", "Reads"))

colnames(LD_R) <- c("Phylum", "Classification", "Count", "Type")
LD_R$Type <- factor(LD_R$Type, levels=c("OTUs", "Reads"))


#Long_Data <- melt(Total_Phyla_Stats, id.vars = "Phylum")
#Long_Data[,4] <- c(rep("Reads", 63), rep("OTUs", 63), rep("Reads", 63), rep("OTUs", 63))
#Long_Data[,2] <- c(rep("DeNovo", 126), rep("GG", 126))

#colnames(Long_Data) <- c("Phylum", "Classification", "Count", "Type")
#Long_Data$Type <- factor(Long_Data$Type, levels=c("OTUs","Reads"))

cbPalette <- c("#bdbdbd", "#636363")

ggplot(LD_A, aes(x=Phylum, y=Count, fill=Classification)) + geom_bar(stat="identity", aes(x=Phylum, y=Count)) + facet_wrap(~Type, nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_text(size=15)) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)
ggsave("../Figures/Phylum_DNvsGG_Abundant.jpg", width=75, units="mm")
ggplot(LD_R, aes(x=Phylum, y=Count, fill=Classification)) + geom_bar(stat="identity", aes(x=Phylum, y=Count)) + facet_wrap(~Type, nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_text(size=15)) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)
ggsave("../Figures/Phylum_DNvsGG_Rare.eps", width=131, units="mm")

ggplot(Long_Data, aes(x=Phylum, y=Count, fill=Classification)) + geom_bar(stat="identity", aes(x=Phylum, y=Count)) + facet_wrap(~Type, nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_text(size=15)) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)


cbPalette <- c("#f0f0f0", "#636363")
pa <- ggplot(data=PlotData_A, aes(x=Phylum_Column, y=Proportion, fill=Fill_Categories)) + geom_bar(stat="identity", aes(x=Phylum_Column, y=Proportion), colour="black") + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_text(size=15)) + guides(fill=FALSE) + scale_fill_manual(values=cbPalette)

pa <- ggplot(data=CA_P, aes(x=Phylum, y=Proportion)) + geom_bar(stat="identity", aes(x=Phylum,y=Proportion)) + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(),strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_blank()) + guides(fill=FALSE)


DNR_Rel <- (Total_Phyla_Stats$DeNovoReads/(Total_Phyla_Stats$DeNovoReads+Total_Phyla_Stats$GGReads))
GGR_Rel <- (Total_Phyla_Stats$GGReads/(Total_Phyla_Stats$DeNovoReads+Total_Phyla_Stats$GGReads))
DNR_Rel + GGR_Rel



full <- full_join(dn_Phyla_Stats, gg_Phyla_Stats, by= V1)
?melt

#Correlation test between Novelty and Soil Temperature
cor.test(map$SoilTemperature_to10cm, map$PercentDeNovo)

fit <- lm(map$PercentDeNovo ~ map$SoilMoisture_Per)
summary(fit)
anova(fit)

grepl()
### Broken
#fig2=ggplot(map.long, aes(y=as.numeric(SoilTemperature_to10cm), x=value))
#add points layer
#geom_point(aes(y=as.numeric(SoilTemperature_to10cm), x=value, shape=Classification, color=as.numeric(SoilTemperature_to10cm)))
### Broken



### OTU abundance and Soil Temperature Correlations
comm_rel.mat <- as.matrix(comm_rel)

Temp_Full <- as.numeric(map$SoilTemperature_to10cm)
#Perform Correlation tests 
fireclass_bin <- gsub(pattern = "Reference", replacement = "Recovered", x = fireclass)
ThermSites <- rep(1,18)
ThermSites[Temp_Full<30]<-0 

sum(1*rowSums(comm_rel.mat[,Temp_Full>30])>0)
comm_therm.mat <- comm_rel.mat[rowSums(comm_rel.mat[,Temp_Full>30])>0,] 

coretest.out=NULL
ttest.out = NULL
kendall.out <- NULL
spearman.out <- NULL
kw.out <- NULL
for(i in 1:nrow(comm_therm.mat)){
  results=cor.test(comm_therm.mat[i,],Temp_Full)
  results.t=t.test(comm_therm.mat[i,]~ThermSites)
  results.k <- kruskal.test(comm_therm.mat[i,]~ThermSites)
  spearman <- cor.test(comm_therm.mat[i,], Temp_Full , method="spearman")
  kendall <- cor.test(comm_therm.mat[i,], Temp_Full, method="kendall")
  coretest.out=rbind(coretest.out,c(row.names(comm_therm.mat)[i],results$estimate,results$p.value))
  ttest.out = rbind(ttest.out, c(row.names(comm_therm.mat)[i], results.t[1], results.t[2], results.t[3]))
  kendall.out <- rbind(kendall.out, c(row.names(comm_therm.mat)[i], kendall[1], kendall[2], kendall[3], kendall[4]))
  spearman.out <- rbind(spearman.out, c(row.names(comm_therm.mat)[i], spearman[1], spearman[2], spearman[3], spearman[4]))
 kw.out <- rbind(kw.out, c(row.names(comm_therm.mat)[i], results.k[1], results.k[2], results.k[3])) 
}

sum(1*(p.adjust(ttest.out[,4], method="fdr")<0.05))
sum(1*(p.adjust(coretest.out[,3], method="fdr")<0.05))
spearman.out <- cbind(spearman.out, p.adjust(spearman.out[,4], method="fdr"))
range(spearman.out[,6])
sum(1*(spearman.out[,6]<0.05))
spearsig <- spearman.out[spearman.out[,6]<0.05,]
spearsig[spearsig[,5]>0,]

range(p.adjust(ttest.out[,4], method="fdr"))
hist(p.adjust(ttest.out[,4], method="fdr"))


fo.out=NULL
for(i in 1:nrow(FireOnly)){
  results <- cor.test(as.numeric(FireOnly[i,1:18]), Temp_Full)
  fo.out <- rbind(fo.out, c(row.names(FireOnly)[i], results[1], results[2], results[3]))
}

# OTUs present in C10 and/or C13
therm.otus <- comm_rel[,map$SoilTemperature_to10cm>50]
therm.otus <- therm.otus[rowSums(therm.otus)>0,]
therm.comm <- comm_rel[which(row.names(therm.otus)%in%row.names(comm_rel)),]
therm.out <- NULL
thermkendall.out <- NULL
thermspearman.out <- NULL
for(i in 1:nrow(therm.comm)){
  results <- cor.test(as.numeric(therm.comm[i,1:18]), Temp_Full)
  results.spearman <- cor.test(as.numeric(therm.comm[i,1:18]), Temp_Full, method="spearman")
  results.kendall <- cor.test(as.numeric(therm.comm[i,1:18]), Temp_Full, method="kendall")
  therm.out <- rbind(therm.out, c(row.names(therm.comm)[i], results[1], results[2], results[3]))
  thermspearman.out <- rbind(thermspearman.out, c(row.names(therm.comm)[i], results.spearman[1], results.spearman[2], results.spearman[3], results.spearman[4]))
  thermkendall.out <- rbind(thermkendall.out, c(row.names(therm.comm)[i], results.kendall[1], results.kendall[2], results.kendall[3], results.kendall[4]))
}
therm.out <- cbind(therm.out, p.adjust(therm.out[,4], method="fdr"))
thermspearman.out <- cbind(thermspearman.out, p.adjust(thermspearman.out[,4], method="fdr"))
thermkendall.out <- cbind(thermkendall.out, p.adjust(thermkendall.out[,4], method="fdr"))
sig_spearman <- thermspearman.out[thermspearman.out[,6]<0.05,]
hist(as.numeric(sig_spearman[,5]))
sig_spearman[sig_spearman[,5]>0,]

# All OTUs present in Fire Affected Sites
fire_otus <- comm_rel[rowSums(comm_rel[,fireclass_bin=="FireAffected"])>0,]
fo.out=NULL
for(i in 1:nrow(fire_otus)){
  results <- cor.test(as.numeric(fire_otus[i,1:18]), Temp_Full)
  fo.out <- rbind(fo.out, c(row.names(fire_otus)[i], results[1], results[2], results[3]))
}
hist(as.numeric(fo.out[,4]))
range(p.adjust(as.numeric(fo.out[,4]), method="fdr"))



cor.test(comm_rel.mat[1,], Temp_Full)
cor.test(as.numeric(comm[1,]), Temp_Full)

hist(as.numeric(coretest.out[,3]))

#Only look at OTUs present in 9 or more samples
#Half_or_More <- coretest.out[rowSums(1*(comm_rel>0))>8,]
#range(p.adjust(Half_or_More[,3], method="fdr"))
hist(p.adjust(as.numeric(coretest.out[,3]), method="fdr"))
range(p.adjust(as.numeric(coretest.out[,3]), method="fdr"))

# Only Significant Correlations with Temperature
sigcor <- coretest.out[as.numeric(coretest.out[,3])<.05,]
comm_sigcor <- cbind(comm_rel.mat[as.numeric(coretest.out[,3])<.05,],rdp[as.numeric(coretest.out[,3])<.05])
#Only Postive Significant Correlations with Temp
sigcor_pos <- sigcor[sigcor[,2]>0,]
comm_sigcor_pos <- comm_sigcor[sigcor[,2]>0,]
comm_sigcor_pos <- as.data.frame(comm_sigcor_pos)
hist(as.numeric(sigcor_pos[,2]))

# Designate a dataset with read counts (as opposed to relative abundance) 
Acomm_sigcor <- cbind(comm[as.numeric(coretest.out[,3])<.05,],rdp[as.numeric(coretest.out[,3])<.05])
Acomm_sigcor_pos <- Acomm_sigcor[as.numeric(sigcor[,2])>0,]
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

#Counts<- NULL
#Num_OTUs <- NULL
#for (i in 1:length(Phyla_Therm)){
 # Object <- Acomm_sigcor_pos[grep(Phyla_Therm[i],Taxonomy[,2]),] 
  #Counts <- c(Counts, sum(Object[,1:18]))
  #Num_OTUs <- c(Num_OTUs, nrow(Object))
#}

Counts_dn <- NULL
Num_OTUs_dn <- NULL
Counts_gg<- NULL
Num_OTUs_gg <- NULL
for (i in 1:length(Phyla_Therm)){
  Object <- Acomm_sigcor_pos[grep(Phyla_Therm[i],Taxonomy[,2]),] 
  Object_dn <- NULL
  tryCatch({Object_dn <- Object[grep("dn",row.names(Object)),]},error=function(e){})
  Object_gg <- NULL
  tryCatch({Object_gg <- Object[grep("dn",row.names(Object),invert=TRUE),]},error=function(e){})
  Counts_dn <- c(Counts_dn, tryCatch({sum(Object_dn[,1:18])},error=function(e){0}))
  Counts_gg <- c(Counts_gg, tryCatch({sum(Object_gg[,1:18])},error=function(e){0}))
  Num_OTUs_dn <- c(Num_OTUs_dn, nrow(Object_dn))
  Num_OTUs_gg <- c(Num_OTUs_gg, nrow(Object_gg))
}

Proportion_OTUs <- c(Num_OTUs_gg,Num_OTUs_dn)/sum(c(Num_OTUs_gg,Num_OTUs_dn))
Proportion_RA <- c(Counts_gg,Counts_dn)/sum(c(Counts_gg,Counts_dn))

JGI_Therm <- read.table("JGI_Thermophile_Phylum_Counts_November9th.txt", sep="\t", stringsAsFactors = FALSE, header = FALSE)
JGI_Therm$V3 <- JGI_Therm$V2/sum(JGI_Therm$V2)

Proportion <- c(Proportion_OTUs, Proportion_RA, JGI_Therm[,3])

Facet_Categories <- c(rep("A",84), rep("B",84), rep("C", nrow(JGI_Therm)))
Fill_Categories <- c(rep("GG", 42), rep("DN",42), rep("GG", 42), rep("DN",42), rep("Genomes", nrow(JGI_Therm)))

Phylum_Column <- c(rep(Phyla_Therm,4), JGI_Therm[,1])
Phylum_Column <- sub("Thermi", "Deinococcus-Thermus",Phylum_Column)

Plot_Data <- cbind(Phylum_Column,Fill_Categories)
Plot_Data <- as.data.frame(Plot_Data)

Plot_Data$Proportion <- Proportion
Plot_Data$Facet_Categories <- Facet_Categories
data.class(Plot_Data$Proportion)

Plot_Data$Category_f <- factor(Plot_Data$Facet_Categories, levels=c("A","B","C"))

Lumped_OTUs <- (Num_OTUs_gg+Num_OTUs_dn)/sum(c(Num_OTUs_gg,Num_OTUs_dn))
Lumped_RA <- (Counts_gg+Counts_dn)/sum(c(Counts_gg,Counts_dn))

Col1 <- c(rep(Phyla_Therm,2),JGI_Therm$V1)
Col2 <- c(Lumped_OTUs,Lumped_RA, JGI_Therm$V3)
Dummy_Matrix <- cbind(Col1,Col2)

Phyla_A <- NULL
Phyla_R <- NULL
for(i in 1:length(unique(Col1))){
  x <- Col2[grep(unique(Col1)[i], Col1)]
  if(sum(1*(x<=0.05))==length(x)){
    Phyla_R <- c(Phyla_R,unique(Col1)[i])
  }else{
    Phyla_A <- c(Phyla_A, unique(Col1)[i])
  }
}



PlotData_A <- NULL
for(i in 1:length(Phyla_A)){
  x <- Plot_Data[grep(Phyla_A[i],Plot_Data[,1]),]
  PlotData_A <- rbind(PlotData_A,x)
}

PlotData_R <- NULL
for(i in 1:length(Phyla_R)){
  x <- Plot_Data[grep(Phyla_R[i],Plot_Data[,1]),]
  PlotData_R <- rbind(PlotData_R,x)
}

data.class(PlotData_A$Proportion)

cbPalette <- c("#f0f0f0", "#bdbdbd", "#636363")
pa <- ggplot(data=PlotData_A, aes(x=Phylum_Column, y=Proportion, fill=Fill_Categories)) + geom_bar(stat="identity", aes(x=Phylum_Column, y=Proportion), colour="black") + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_text(size=15)) + guides(fill=FALSE) + scale_fill_manual(values=cbPalette)
pa
ggsave("../Figures/Fig3A.eps", width=50, units="mm")

pr <- ggplot(data=PlotData_R, aes(x=Phylum_Column, y=Proportion, fill=Fill_Categories)) + geom_bar(stat="identity", aes(x=Phylum_Column, y=Proportion), colour="black") + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_text(size=15)) + guides(fill=FALSE) + scale_fill_manual(values=cbPalette)
pr
ggsave("../Figures/Fig3B.eps", width=200, units="mm")

### Trying to Even out Bars
pa <- ggplot(data=CA_P, aes(x=Phylum, y=Proportion)) + geom_bar(stat="identity", aes(x=Phylum,y=Proportion)) + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(),strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_blank()) + guides(fill=FALSE)
pa
ggsave("../Figures/Fig3A.png", width=50, units="mm")

pr <- ggplot(data=CR_P, aes(x=Phylum, y=Proportion)) + geom_bar(stat="identity", aes(x=Phylum,y=Proportion)) + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_blank()) + guides(fill=FALSE)
pr
ggsave("../Figures/Fig3B.png", width=200, units="mm")

### Venn Meta_Analysis
map.f <- read.table("Supplemental/Centralia_Full_Map_Fixed.txt", sep="\t", header= TRUE, row.names=1, stringsAsFactors = FALSE)
meta <- read.table("supplemental/Carini_RDP_rmCM.txt", sep="\t", header = TRUE, row.names=1, stringsAsFactors = FALSE)

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


### KEGG Module Analysis
library(stringr)
Modules <- read.table("kmodlist47982_23-nov-2016.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

MO_Per_M <- rep(0,nrow(Modules))
for (i in 1:nrow(Modules)){
  MO_Per_M[i] <- str_count(Modules$Definition[i],"M")
}

Modules_w_Modules <- Modules[MO_Per_M>0,]

Dictionary=NULL
Dictionary <- vector(mode="list", length=nrow(Modules))
names(Dictionary) <- row.names(Modules)

### For Some Module definitions, they have other modules in them, this is formating those definitions so that there are only K0 and no M0's in the Definition. As a side note these are all "Signature" module types
#Acetogen
Modules["M00618",3] <- paste(c(Modules["M00377",3],Modules["M00579",3]), collapse=" ")
#Anoxygenic photosynthesis in green nonsulfur bacteria
Modules["M00613",3] <- paste(c(Modules["M00597",3],Modules["M00376",3]), collapse=" ")
#Anoxygenic photosynthesis in green sulfur bacteria
Modules["M00614",3] <- paste(c(Modules["M00598",3],Modules["M00173",3]), collapse=" ")
#Anoxygenic photosynthesis in purple bacteria
Modules["M00612",3] <- paste(c(Modules["M00597",3],Modules["M00165",3]), collapse=" ")
#Methanogen
Modules["M00617",3] <- paste(c(Modules["M00567",3],Modules["M00357",3],Modules["M00356",3],Modules["M00563",3]), collapse=" ")
#Nitrate assimilations
Modules["M00615",3] <- paste(c(Modules["M00438",3],Modules["M00531",3]), collapse=" ")
#Oxygenic photosynthesis in plants and cyanobacteria 
Modules["M00611",3] <- paste(c(Modules["M00161",3],Modules["M00163",3],Modules["M00165",3]), collapse=" ")
#Sulfate-sulfur assimilation
Modules["M00616",3] <- paste(c(Modules["M00185",3],Modules["M00176",3]), collapse=" ")

#Produces a list, each item in this list is a dataframe with the KO's in our dataset 
#Relativized to rpoB
#for (y in 1:nrow(KO.rpob)){
 # KO_M <- grep(row.names(KO.rpob)[y], Modules$Definition)
  #if (length(KO_M)>0){
   # for (x in 1:length(KO_M)){
    #  Dictionary[[KO_M[x]]] <- rbind(Dictionary[[KO_M[x]]], KO.rpob[y,])
    #}
  #}
#}
# Relativized to rpsM
#for (y in 1:nrow(KO.rpsM)){
 # KO_M <- grep(row.names(KO.rpsM)[y], Modules$Definition)
  #if (length(KO_M)>0){
   # for (x in 1:length(KO_M)){
    #  Dictionary[[KO_M[x]]] <- rbind(Dictionary[[KO_M[x]]], KO.rpsM[y,])
  #  }
#  }
#}

# Relativized to total KO Count
#for (y in 1:nrow(KO.rel)){
 # KO_M <- grep(row.names(KO.rel)[y], Modules$Definition)
  #if (length(KO_M)>0){
   # for (x in 1:length(KO_M)){
    #  Dictionary[[KO_M[x]]] <- rbind(Dictionary[[KO_M[x]]], KO.rel[y,])
    #}
  #}
#}

# Relativized to Average Single Copy Gener count
for (y in 1:nrow(KO.sr)){
  KO_M <- grep(row.names(KO.sr)[y], Modules$Definition)
  if (length(KO_M)>0){
    for (x in 1:length(KO_M)){
      Dictionary[[KO_M[x]]] <- rbind(Dictionary[[KO_M[x]]], KO.sr[y,])
    }
  }
}

### Number of modules not represented by dataset
count <- 0
for (y in 1:length(Dictionary)){
  if (data.class(Dictionary[[y]])=="NULL"){
    count <- count + 1
  }
}
count

#Counts the number of KOs in each module
library(stringr)
KO_Per_M <- rep(0, nrow(Modules))
for (i in 1:nrow(Modules)){
  KO_Per_M[i] <- str_count(Modules$Definition[i],"K")
}
# Counts number of KOs from each module that are present in our dataset
KO_Per_M_Data <- rep(0, nrow(Modules))
for (i in 1:nrow(Modules)){
  if(data.class(Dictionary[[i]])!="NULL"){
    KO_Per_M_Data[[i]] <- nrow(Dictionary[[i]])
  }
}

Missing_KOs <- KO_Per_M - KO_Per_M_Data

KO_Per_M_Present <- KO_Per_M[KO_Per_M_Data>0]
KO_Per_M_Data_Present <- KO_Per_M_Data[KO_Per_M_Data>0]
### Number of Modules that are less than 50% complete 
sum(1*((KO_Per_M_Data_Present/KO_Per_M_Present)<0.5))

### Adding zeroes for K0 missing in our data set from some modules (doesnt make a big difference, commenting it out because of this)
#Dictionary_AddedZeroes <- Dictionary
#for( i in 1:length(Missing_KOs)){
 # if (Missing_KOs[i]>0){
  #  Addition <- matrix(0,Missing_KOs[i],12)
   # colnames(Addition) <- colnames(Dictionary_AddedZeroes[[i]])
    #Dictionary_AddedZeroes[[i]] <- rbind(Dictionary_AddedZeroes[[i]],Addition)
  #}
#}
### Number of incomplete modules 367 modules that are not 100% complete
length(Missing_KOs[Missing_KOs>0])

KO_Per_M_Data/KO_Per_M


### Getting Rid of Modules that are not present at all in our data. 
Dictionary_Subset <- Dictionary[lapply(Dictionary, data.class) == "data.frame"]
#Dictionary_AZ_Subset <- Dictionary_AddedZeroes[lapply(Dictionary_AddedZeroes, data.class) == "data.frame"]
Modules_Subset <- Modules[lapply(Dictionary,data.class) == "data.frame",]
KO_Per_M_Subset <- KO_Per_M[lapply(Dictionary, data.class) == "data.frame"]
KO_Per_M_Data_Subset <- KO_Per_M_Data[lapply(Dictionary, data.class) == "data.frame"]

#Fraction of KOs from each module present in our dataset
Module_Completeness <- KO_Per_M_Data_Subset/KO_Per_M_Subset

### Figure out correlation coefficients? Just do t.tests? 
# the real question is whether to use the worst, average, or best correlation/t.test result
#Dummy Code
#  for every module 
#  for every KO in that module
#    Calculate the correlation coefficient
#    calculate t.test result
#    write correlation coefficient and t. test results to an item in a list
#  Calculate average and standard deviation of correlation and cor pvalue
#  calculate average and standard deviation of t.test statistic and pvalue

#Calculating correlation coefficient and T.test results for all KOs in all modules  
Summary_Output <- NULL
avg_mod <- NULL
mid_mod <- NULL
max_mod <- NULL
min_mod <- NULL
stdev_mod <- NULL
for (m in 1:length(Dictionary_Subset)){
  temp_data <- Dictionary_Subset[[m]]
  K_out <- NULL
  agg <- NULL
  ### Calculates the AVG, Median, Min, Max, and STDEV for the KOs of a module in each given site.
  for (s in 1:ncol(temp_data)){
    avg <- mean(temp_data[,s])
    mid <- median(temp_data[,s])
    minimum <- min(temp_data[,s])
    maximum <- max(temp_data[,s])
    stdev <- sd(temp_data[,s])
    agg <- cbind(agg, c(avg, mid, minimum, maximum, stdev))
  }
  avg_mod <- rbind(avg_mod, agg[1,])
  mid_mod <- rbind(mid_mod, agg[2,])
  min_mod <- rbind(min_mod, agg[3,])
  max_mod <- rbind(max_mod, agg[4,])
  stdev_mod <- rbind(stdev_mod, agg[5,])
}

###colnames(Summary_Output) <- c("Number KOs", "Average Pearsons", "Average Pearsons Pvalue", "Average T Statistic", "Average T P value", "SD Pearsons", "SD Pearsons Pvalue", "SD T Statistic", "SD T Pvalue") 
colnames(avg_mod) <- map_MG$Sample
colnames(mid_mod) <- map_MG$Sample
colnames(stdev_mod) <- map_MG$Sample
colnames(mid_mod)<- map_MG$Sample
colnames(max_mod) <- map_MG$Sample
 
row.names(avg_mod) <- names(Dictionary_Subset)
row.names(mid_mod) <- names(Dictionary_Subset)
row.names(max_mod) <- names(Dictionary_Subset)
row.names(min_mod) <- names(Dictionary_Subset)
row.names(stdev_mod) <- names(Dictionary_Subset)


KO_Per_M_Data_Present <- KO_Per_M_Data_Present[rowSums(mid_mod)!=0]
KO_Per_M_Present <- KO_Per_M_Present[rowSums(mid_mod)!=0]
mid_mod <- mid_mod[rowSums(mid_mod)!=0,]


###   ***Multivariate Analysis based on Modules, KO, and 16S*** 
# Distance Matrix based on Weighted unifrac
uf=read.table("weighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)
uf=uf[order(row.names(uf)),order(colnames(uf))]
# Subsetting to only the samples we have metagenomes for...
uf_mg <- uf[c(1,3,4,5,6,7,10,12,14,15,16,17),c(1,3,4,5,6,7,10,12,14,15,16,17)]
uf_mg.d <- as.dist(uf_mg)

# Distance matrix based on Modules
Mod.d <- vegdist(t(mid_mod),method="bray")
Class <- rep("Red", 12)
Class[map_MG$Classification!="FireAffected"] <- "Green"
a=adonis(Mod.d~Class, distance=TRUE, permutations=1000)
a

b=betadisper(Mod.d, group=Class)
TukeyHSD(b, which = "group", ordered = FALSE,conf.level = 0.95)

# Distance Matrix based on KO rel
KO.d <- vegdist(t(KO.sr),method="bray")

# Setting up Classifications
Class <- rep("Red", 12)
Class[map_MG$Classification!="FireAffected"] <- "Green"

#PERMANOVA For KO, Module, and UniFrac
Module.a=adonis(Mod.d~Class, distance=TRUE, permutations=1000)
Module.a # Significant 
KO.a <- adonis(KO.d~Class, distance=TRUE, permutations = 1000)
KO.a # Significant
uf.a <- adonis(uf_mg.d~Class, distance=TRUE, permutations= 1000)
uf.a # Significant


# Beta Dispersion for KO, Module, and UniFrac
Module.b=betadisper(Mod.d, group=Class)
TukeyHSD(Module.b, which = "group", ordered = FALSE,conf.level = 0.95)

KO.b=betadisper(KO.d, group=Class)
TukeyHSD(KO.b, which = "group", ordered = FALSE,conf.level = 0.95)

uf.b=betadisper(uf_mg.d, group=Class)
TukeyHSD(uf.b, which = "group", ordered = FALSE,conf.level = 0.95)

Module.pcoa <- cmdscale(Mod.d, eig=TRUE)
KO.pcoa <- cmdscale(KO.d, eig=TRUE)
uf.pcoa <- cmdscale(uf_mg.d, eig=TRUE)
# Based on the help file for PROTEST it looks like it requires a PCoA to run the test, although maybe using a distance matrix can work as well (it just has # of dimensions== # of samples instead of # of PCoA axes)

# Module VS KO
M_v_K.mantel <- mantel(Mod.d, KO.d)
M_v_K.protest <- protest(Module.pcoa, KO.pcoa)

M_v_U.mantel <- mantel(Mod.d, uf_mg.d)
M_v_U.protest <- protest(Module.pcoa,uf.pcoa)

K_v_U.mantel <- mantel(KO.d, uf_mg.d)
K_v_U.protest <- protest(KO.pcoa, uf.pcoa)

space.d <- read.table("spatialdistancematrix.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
space.d <- as.matrix(space.d)

uf.m <- as.matrix(uf)
mantel(space.d, uf.m)
as.numeric(space.d)
plot(as.numeric(space.d),as.numeric(as.matrix(uf)))

space.d.mg <- space.d[c(1,3,4,5,6,7,10,12,14,15,16,17),c(1,3,4,5,6,7,10,12,14,15,16,17)]

plot(as.numeric(space.d.mg), as.matrix(KO.d))
plot(as.numeric(space.d.mg), as.matrix(uf_mg.d))

M_v_Space.mantel <- mantel(Mod.d, space.d.mg)
K_v_Space.mantel <- mantel(KO.d, space.d.mg)
U_v_Space.mantel <- mantel(uf_mg.d, space.d.mg)

Distance1 <-c("Module", "Module", "Module", "KO", "KO", "UniFrac")
Distance2<- c("KO", "UniFrac", "Space","UniFrac", "Space", "Space")
Mantel_Names<- c("Distance1", "Distance2","Mantel_R", "p-value")
Mantel_R <- c(M_v_K.mantel[3], M_v_U.mantel[3], M_v_Space.mantel[3], K_v_U.mantel[3], K_v_Space.mantel[3], U_v_Space.mantel[3])
Mantel_Pvalues <- c(M_v_K.mantel[4], M_v_U.mantel[4], M_v_Space.mantel[4], K_v_U.mantel[4], K_v_Space.mantel[4], U_v_Space.mantel[4])
Mantel_Summary <- cbind(Distance1, Distance2, Mantel_R, Mantel_Pvalues)
write.table(as.matrix(Mantel_Summary),"Supplemental/Supplemental_Mantel.txt", sep="\t", row.names=FALSE, quote=FALSE)
# Protest and Mantel test reveal all three datasets resemble each other


#pcoa_pro <- protest(Module.pcoa,uf.pcoa)
#dist_pro <- protest(Mod.d, uf_mg.d)
#summary(dist_pro)
#summary(pcoa_pro)
#plot(pcoa_pro)
#plot(pcoa_pro, kind=1)
#plot(dist_pro)



### *** Module Temperature Correlations, T-Tests, and Heatmaps *** 
median_Module_TempCorrelations <- apply(mid_mod, 1, function(x) cor.test(x, map_MG$SoilTemperature_to10cm))


# T-Test for each Module
median_Module_t.test <- apply(mid_mod, 1, function(x) t.test(x[map_MG$Classification=="FireAffected"], x[map_MG$Classification!="FireAffected"]))


### 229 Modules significantly correlated with temperature
Med_Mod_TempCor <- NULL
for (i in 1:length(median_Module_TempCorrelations)){
  Med_Mod_TempCor <- rbind (Med_Mod_TempCor, c(median_Module_TempCorrelations[[i]][4],median_Module_TempCorrelations[[i]][2],median_Module_TempCorrelations[[i]][3]))
}
Med_Mod_TempCor <- as.data.frame(Med_Mod_TempCor)
Med_Mod_TempCor$KEGG <- row.names(mid_mod)
Med_Mod_TempCor[,1] <- unlist(Med_Mod_TempCor[,1])
Med_Mod_TempCor[,2] <- unlist(Med_Mod_TempCor[,2])
Med_Mod_TempCor[,3] <- unlist(Med_Mod_TempCor[,3])
Med_Mod_TempCor$Completeness <- KO_Per_M_Data_Present/KO_Per_M_Present
Med_Mod_TempCor$Adjusted.p.value <- p.adjust(Med_Mod_TempCor$p.value, "fdr")
row.names(Med_Mod_TempCor)<- Med_Mod_TempCor$KEGG

Sig_Temp_Cor <- Med_Mod_TempCor[Med_Mod_TempCor$Adjusted.p.value<0.05,]

Med_Mod_Ttest <- NULL
for (i in 1:length(median_Module_t.test)){
  Med_Mod_Ttest <- rbind (Med_Mod_Ttest, c(median_Module_t.test[[i]][1],median_Module_t.test[[i]][2],median_Module_t.test[[i]][3]))
}
Med_Mod_Ttest <- as.data.frame(Med_Mod_Ttest)
Med_Mod_Ttest$KEGG <- row.names(mid_mod)
Med_Mod_Ttest[,1] <- unlist(Med_Mod_Ttest[,1])
Med_Mod_Ttest[,2] <- unlist(Med_Mod_Ttest[,2])
Med_Mod_Ttest[,3] <- unlist(Med_Mod_Ttest[,3])
Med_Mod_Ttest$Adjusted.p.value <- p.adjust(Med_Mod_Ttest$p.value, "fdr")
row.names(Med_Mod_Ttest) <- Med_Mod_Ttest$KEGG
Sig_Ttest_Modules <- Med_Mod_Ttest[Med_Mod_Ttest$Adjusted.p.value<0.05,]
colnames(Sig_Ttest_Modules) <- c("T-statistic", "T-parameter", "T-p.value","KEGG", "T-Adjusted.p.value")
library(plyr)
Combined_Sig_Module_Results <- join(Sig_Temp_Cor, Sig_Ttest_Modules, by ="KEGG", type="full")
Temp_Cor_Vector <- rep(0, nrow(Combined_Sig_Module_Results))
Temp_Cor_Vector[Combined_Sig_Module_Results$Adjusted.p.value<0.05]=1
Temp_Cor_Vector[is.na(Temp_Cor_Vector)]<-0

Ttest_Sig_Vector <- rep(0,nrow(Combined_Sig_Module_Results))
Ttest_Sig_Vector[Combined_Sig_Module_Results$`T-Adjusted.p.value`<0.05]=2

Combined_Sig_Modules <- as.data.frame(mid_mod)[Combined_Sig_Module_Results$KEGG,]


Sig_Code <- Temp_Cor_Vector + Ttest_Sig_Vector
### Blue == T Test   Red= Temp Cor   Green == Both
Classification_Test <- rep("Blue", nrow(Combined_Sig_Modules)) 
Classification_Test[Sig_Code==3] <- "Green"
Classification_Test[Sig_Code==2] <- "Blue"
Classification_Test[Sig_Code==1] <- "Red"
### Heat map of Modules that are either significant based in T-Test or Temperature Correlation
library(colorRamps)
library(gplots)
library(vegan)
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")


Combined_Sig_Modules.zs <- decostand(as.matrix(Combined_Sig_Modules), method="standardize", MARGIN=1)
par(mfrow=c(1,1))
setEPS()

png("../Figures/Heatmap_Combined_averageSCG.png", width = 2000, height=4000, pointsize=8)
fig10 <- heatmap.2(Combined_Sig_Modules.zs, col=hc(100), key=FALSE, symkey=FALSE, trace="none", density.info="none", colsep=c(1:12),rowsep=c(1:nrow(Combined_Sig_Modules.zs)), sepcolor="black", sepwidth=c(0.01,0.00001), dendrogram="row",cexRow = 2, labRow=FALSE, margins=c(5,13), srtCol=90, colRow = Classification_Test, RowSideColors=Classification_Test, lhei = c(1,100))
dev.off()
### For the Legend
png("../Figures/Heatmap_Legend_rarefied_rpoB.png", width = 2000, height=4000, pointsize=8)
fig10 <- heatmap.2(Combined_Sig_Modules.zs, col=hc(100), key=TRUE, symkey=FALSE, trace="none", density.info="none", colsep=c(1:12),rowsep=c(1:nrow(Combined_Sig_Modules.zs)), sepcolor="black", sepwidth=c(0.01,0.00001), dendrogram="row",cexRow = 2, labRow=FALSE, margins=c(5,13), srtCol=90, colRow = Classification_Test, RowSideColors=Classification_Test, lhei = c(1,10))
dev.off()

### Create Supplemental Table of Summaries
SignificantModules_Summary <-join(Med_Mod_TempCor[row.names(Combined_Sig_Modules),], Med_Mod_Ttest[row.names(Combined_Sig_Modules),], by="KEGG")
SignificantModules_Summary$ModuleDescription <- Modules[row.names(Combined_Sig_Modules),1]


SignificantModules_Summary<-SignificantModules_Summary[,c(4,11,5,1,2,3,6,7,8,9,10)]
colnames(SignificantModules_Summary) <- c("Module","Module Description","Completeness","Pearson's Rho", "Pearson's Degrees Freedom", "Peason p value", "Pearson Adjusted p value", "T Statistic", "T Degrees Freedom", "T-test p value", "T Adjusted p value")
row.names(SignificantModules_Summary) <- row.names(Combined_Sig_Modules)
write.table(x = SignificantModules_Summary, file="Supplemental/SupplementalTable2_Averaged_SCG.txt", sep="\t", quote=FALSE)









### Indicator Module Analysis
library(indicspecies)
IndicClusters <- rep(1,12)
IndicClusters[map_MG$Classification=="Recovered"]=2
IndicClusters[map_MG$Classification=="Reference"]=3
tmidmod<- t(mid_mod)
tmidmod <- as.data.frame(tmidmod)
B=strassoc(tmidmod, cluster=IndicClusters,func="B")
sel=which(B[,1]>0.2)

indicators(tmidmod[,sel], cluster=IndicClusters, group="1",At=0.5,Bt=0.2, verbose=TRUE)
wetpt = multipatt(tmidmod, IndicClusters, )

### Finding Maximum Site for each Module
hist(apply(Combined_Sig_Modules.zs, 1, which.max))
hist(apply(mid_mod,1,which.max))

### PCoA's for Module, Ortholog, and UF distacnes
library(calibrate)
class <- rep("black", nrow(map_MG))
class[map_MG$Classification=="Recovered"]='yellow'
class[map_MG$Classification=="FireAffected"]='red'
class[map_MG$Classification=="Reference"]='green'

# % variance explained on 1st two axes of Module, KO, and UF PCoA
M_ax1.v <- Module.pcoa$eig[1]/sum(Module.pcoa$eig)
M_ax2.v <- Module.pcoa$eig[2]/sum(Module.pcoa$eig)

KO_ax1.v=KO.pcoa$eig[1]/sum(KO.pcoa$eig)
KO_ax2.v=KO.pcoa$eig[2]/sum(KO.pcoa$eig)

uf_ax1.v=uf.pcoa$eig[1]/sum(uf.pcoa$eig)
uf_ax2.v=uf.pcoa$eig[2]/sum(uf.pcoa$eig)


#calculate percent variance explained, then add to plot
ax1.v.uf=uf.pcoa$eig[1]/sum(uf.pcoa$eig)
ax2.v.uf=uf.pcoa$eig[2]/sum(uf.pcoa$eig)

env=map_MG[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]

dev.off()
setEPS()
postscript("../Figures/Figure5.eps", width = 15, height=15, pointsize=12,paper="special")
#par(mfrow=c(1,3))
plot(.pcoa$points[,1], Module.pcoa$points[,2],cex=4, bg=class, pch=21, main= "Single Copy Gene Median Relativized Bray Curtis KEGG Module PCoA", xlab= paste("PCoA1: ",100*round(M_ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(M_ax2.v,3),"% var. explained",sep=""))
textxy(X=Module.pcoa$points[,1],Y=Module.pcoa$points[,2], lab=map_MG$Sample,cex=2)
M_env<- envfit(Module.pcoa,env)
M_moduleenv <- envfit(Module.pcoa, t(mid_mod))
plot(M_env, p.max=0.05, col="black")

dev.off()
setEPS()
postscript("../Figures/Figure5.eps", width = 15, height=15, pointsize=12,paper="special")
plot(KO.pcoa$points[,1], KO.pcoa$points[,2],cex=3, bg=class, pch=21, main= "SCG Median Relativized Bray Curtis KEGG Ortholog PCoA", xlab= paste("PCoA1: ",100*round(KO_ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(KO_ax2.v,3),"% var. explained",sep=""))
textxy(X=KO.pcoa$points[,1], Y=KO.pcoa$points[,2],labs=map_MG$Sample, cex=2)

KO_env=envfit(KO.pcoa, env)
plot(KO_env, p.max=0.05, col="black")
dev.off()
#plot(uf.pcoa$points[,1],uf.pcoa$points[,2] ,cex=3,pch=21,bg=class,main="Weighted UniFrac PCoA",xlab= paste("PCoA1: ",100*round(uf_ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(uf_ax2.v,3),"% var. explained",sep=""))
#textxy is from the calibrate library
#textxy(X=uf.pcoa$points[,1], Y=uf.pcoa$points[,2],labs=map_MG$Sample, cex=2)
#uf_env <- envfit(uf.pcoa, env)
#plot(uf_env, p.max=0.05, col="black")
dev.off()




#List of Unclassified OTUs

sum(1*grepl("Unclassified", rdp))
grep("Bacteria", rdp)


FULL_Taxonomy <- NULL
for (i in 1:nrow(comm.rdp)){
  FULL_Taxonomy <- rbind(FULL_Taxonomy,unlist(strsplit(as.character(comm.rdp[i,19]), ";")))
}
row.names(comm.rdp)[grep("Bacteria", FULL_Taxonomy[,2])]
write.table(row.names(comm.rdp)[grep("Bacteria", FULL_Taxonomy[,2])],"Unclassified_Bacteria_OTUs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
?write.table
write.table(row.names(comm.rdp)[grep("Archaea", FULL_Taxonomy[,2])], "Unclassified_Archaea_OTUs.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
