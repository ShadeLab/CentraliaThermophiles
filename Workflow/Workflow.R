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


### DeNovo Vs Reference OTU Analysis -WIP-
library(reshape2)
library(ggplot2)
Percent_DN <- colSums(comm[grepl("dn",row.names(comm)),])/colSums(comm)
map$PercentDeNovo <- Percent_DN

map.long=melt(map, id.vars=c("Sample", "SoilTemperature_to10cm", "Classification"), measure.vars=c("NO3N_ppm","NH4N_ppm","pH","SulfateSulfur_ppm","K_ppm","Ca_ppm","Mg_ppm","OrganicMatter_500","Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","PercentDeNovo"))

p <- ggplot(map, aes(SoilTemperature_to10cm, PercentDeNovo))
p+ geom_point()
ggsave("../Figures/Fig1.eps", width=178, units="mm")

fit <- lm(map$PercentDeNovo ~ map$SoilMoisture_Per)
summary(fit)
anova(fit)

p2 <- ggplot(map, aes(SoilTemperature_to10cm, log(PercentDeNovo)))
p2+ geom_point()

fit2 <- lm(log(map$PercentDeNovo)~map$SoilMoisture_Per)
summary(fit2)

p3 <- ggplot(map, aes(log(SoilTemperature_to10cm), PercentDeNovo))
p3+ geom_point()

fit3 <- lm(map$PercentDeNovo~log(map$SoilMoisture_Per))
summary(fit3)
### Broken
#fig2=ggplot(map.long, aes(y=as.numeric(SoilTemperature_to10cm), x=value))
#add points layer
#geom_point(aes(y=as.numeric(SoilTemperature_to10cm), x=value, shape=Classification, color=as.numeric(SoilTemperature_to10cm)))
### Broken



#Positively Temperature Correlated OTUs
comm_rel.mat <- as.matrix(comm_rel)

Temp_Full <- as.numeric(map$SoilTemperature_to10cm)
coretest.out=NULL
TC_OTUs_R=NULL
TC_OTUs_A=NULL
for(i in 1:nrow(comm_rel)){
  results=cor.test(comm_rel.mat[i,],Temp_Full)
  if(results[4]>0){
    if(results[3]<0.05){
      TC_OTUs_R=rbind(TC_OTUs_R,c(comm_rel[i,],rdp[i]))
      TC_OTUs_A=rbind(TC_OTUs_A,c(comm[i,],rdp[i]))
      coretest.out=rbind(coretest.out,c(row.names(comm_rel)[i],results$estimate,results$p.value))
    }
  }
  
}

rownames(TC_OTUs_R)=coretest.out[,1]
rownames(TC_OTUs_A)=coretest.out[,1]



#All Significant(p<0.05) Correlations
full.out=NULL
TC_OTUs_RA=NULL
TC_OTUs_AA=NULL
for(i in 1:nrow(comm_rel)){
  results=cor.test(comm_rel.mat[i,],Temp_Full)
    if(results[3]<0.05){
      TC_OTUs_RA=rbind(TC_OTUs_RA,c(comm_rel[i,],rdp[i]))
      TC_OTUs_AA=rbind(TC_OTUs_AA,c(comm[i,],rdp[i]))
      full.out=rbind(full.out,c(row.names(comm_rel)[i],results$estimate,results$p.value))
  }
  
}


