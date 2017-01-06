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
# Remove KOs with zeroes across the entire dataset
KO <- KO[rowSums(KO)>0,]

#Relativze to rpoB
KO.rpob <- NULL
for (i in 1:nrow(KO)){
  KO.rpob <- rbind(KO.rpob, KO[i,]/rpoB)
}

row.names(KO.rpob) <- sub("KO:", "", row.names(KO.rpob))

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

Combined$Phylum <- sub("Thermi", "Deinococcus-Thermus", Combined$Phylum)

Names <- unique(Combined$Phylum)
CA_P <- NULL
CR_P <- NULL
for(i in 1:length(Names)){
  x <- Combined[grep(Names[i],Combined$Phylum, fixed=TRUE),]
  if(sum(1*(x[,2]<=0.05))==nrow(x)){
    CR_P <- rbind(CR_P, x)
  }else{
    CA_P <- rbind(CA_P,x)
  }
    
  }

row.names(CA_P) <- c(1:nrow(CA_P))
# Add the unclassifieds to the Rare table
CR_P <- rbind(CR_P, CA_P[36,])
CR_P <- rbind(CR_P, CA_P[33,])

# Remove unclassified archaea, unclassifieds, and duplicated unclassified_bacteria from the abundant table. 
CA_P <- CA_P[-37,]
CA_P <- CA_P[-36,]
CA_P <- CA_P[-35,]
CA_P <- CA_P[-34,]
CA_P <- CA_P[-33,]
CA_P <- CA_P[-32,]

nrow(CA_P) + nrow(CR_P) == nrow(Combined)

pa <- ggplot(data=CA_P, aes(x=Phylum, y=Proportion)) + geom_bar(stat="identity", aes(x=Phylum,y=Proportion)) + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(),strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=rel(1)), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2))) + guides(fill=FALSE)
pa
ggsave("../Figures/Fig3A.eps", width=200, units="mm")

pr <- ggplot(data=CR_P, aes(x=Phylum, y=Proportion)) + geom_bar(stat="identity", aes(x=Phylum,y=Proportion)) + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=rel(1)), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2))) + guides(fill=FALSE)
pr
ggsave("../Figures/Fig3B.eps", width=200, units="mm")

##p2 <- ggplot(data=Combined, aes(x=Phylum, y=Proportion, fill=Phylum)) + geom_bar(stat="identity",aes(x=Phylum,y=Proportion)) + facet_grid(~Category_f, fac) + theme(strip.background = element_blank(),axis.text.x = element_text(angle=90, hjust=1, vjust=.5,  face="bold", lineheight = rel(2), size=rel(1.2)), axis.title.x=element_text(vjust=-1), strip.text.x=element_text(size=rel(2))) + guides(fill=FALSE)
##p2
#ggsave("../Figures/Fig3.eps", width=200, units="mm")

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


### KEGG Module Analysis
Modules <- read.table("kmodlist47982_23-nov-2016.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

MO_Per_M <- rep(0,nrow(Modules))
for (i in 1:nrow(Modules)){
  MO_Per_M[i] <- str_count(Modules$Definition[i],"M")
}

Modules_w_Modules <- Modules[MO_Per_M>0,]

Dictionary <- vector(mode="list", length=nrow(Modules))
names(Dictionary) <- row.names(Modules)

Dictionary[[1]]
grep(row.names(KO.rpob)[10], Modules$Definition)
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
for (KO in 1:nrow(KO.rpob)){
  KO_M <- grep(row.names(KO.rpob)[KO], Modules$Definition)
  if (length(KO_M)>0){
    for (x in 1:length(KO_M)){
      Dictionary[[KO_M[x]]] <- rbind(Dictionary[[KO_M[x]]], KO.rpob[KO,])
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
#mid_mod_az <- NULL
#avg_mod_az <- NULL
max_mod <- NULL
min_mod <- NULL
stdev_mod <- NULL
for (m in 1:length(Dictionary_Subset)){
  temp_data <- Dictionary_Subset[[m]]
 # temp_data_az <- Dictionary_AZ_Subset[[m]]
  K_out <- NULL
  agg <- NULL
  #agg_az <- NULL
  ### This loop calculates the correlation coefficients and t test for each KO in a module. We've decided to go ahead and calculate these for the median values 
  #for(k in 1:nrow(temp_data)){
  # cor <- cor.test(as.numeric(temp_data[k,]), map_MG$SoilTemperature_to10cm) 
  #test <- t.test(as.numeric(temp_data[k,map_MG$Class_Color=="red"], as.numeric(temp_data[k, map_MG$Class_Color!="red"])))
  #K_out <- rbind(K_out, c(as.numeric(cor[1]), as.numeric(cor[3]), as.numeric(test[1]), as.numeric(test[3])))
  #}
  ### Calculates the AVG, Median, Min, Max, and STDEV for the KOs of a module in each given site.
  for (s in 1:ncol(temp_data)){
    avg <- mean(temp_data[,s])
    mid <- median(temp_data[,s])
    minimum <- min(temp_data[,s])
    maximum <- max(temp_data[,s])
    stdev <- sd(temp_data[,s])
    agg <- cbind(agg, c(avg, mid, minimum, maximum, stdev))
  }
  #for (t in 1:ncol(temp_data_az)){
   # avg_az <- mean(temp_data_az[,t])
    #mid_az <- median(temp_data_az[,t])
    #agg_az <- cbind(agg_az, c(avg_az, mid_az))
#  }
  #Summary_Output <- rbind(Summary_Output, c(nrow(K_out),apply(K_out, 2, mean), apply(K_out, 2, sd)))
  avg_mod <- rbind(avg_mod, agg[1,])
  mid_mod <- rbind(mid_mod, agg[2,])
  min_mod <- rbind(min_mod, agg[3,])
  max_mod <- rbind(max_mod, agg[4,])
  stdev_mod <- rbind(stdev_mod, agg[5,])
  #avg_mod_az <- rbind(avg_mod_az, agg_az[1,])
 # mid_mod_az <- rbind(mid_mod_az, agg_az[2,])
}

###colnames(Summary_Output) <- c("Number KOs", "Average Pearsons", "Average Pearsons Pvalue", "Average T Statistic", "Average T P value", "SD Pearsons", "SD Pearsons Pvalue", "SD T Statistic", "SD T Pvalue") 
colnames(avg_mod) <- map_MG$Sample
colnames(mid_mod) <- map_MG$Sample
colnames(stdev_mod) <- map_MG$Sample
colnames(mid_mod)<- map_MG$Sample
colnames(max_mod) <- map_MG$Sample
#colnames(avg_mod_az) <- map_MG$Sample
#colnames(mid_mod_az) <- map_MG$Sample

#row.names(Summary_Output) <- names(Dictionary_Subset)  
row.names(avg_mod) <- names(Dictionary_Subset)
row.names(mid_mod) <- names(Dictionary_Subset)
row.names(max_mod) <- names(Dictionary_Subset)
row.names(min_mod) <- names(Dictionary_Subset)
row.names(stdev_mod) <- names(Dictionary_Subset)
#row.names(avg_mod_az) <- names(Dictionary_AZ_Subset)
#row.names(mid_mod_az) <- names(Dictionary_AZ_Subset)

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

# Distance Matrix based on KO
KO.d <- vegdist(t(KO.rpob),method="bray")

# Settin gup Classifications
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

#med_mod_az_TC <- apply(mid_mod_az, 1, function(x) cor.test(x, map_MG$SoilTemperature_to10cm))

# T-Test for each Module
median_Module_t.test <- apply(mid_mod, 1, function(x) t.test(x[map_MG$Classification=="FireAffected"], x[map_MG$Classification!="FireAffected"]))

#med_mod_az_t.test <- apply(mid_mod_az, 1, function(x) t.test(x[map_MG$Classification=="FireAffected"], x[map_MG$Classification!="FireAffected"]))

### 229 Modules significantly correlated with temperature
Med_Mod_TempCor <- NULL
for (i in 1:607){
  Med_Mod_TempCor <- rbind (Med_Mod_TempCor, c(median_Module_TempCorrelations[[i]][4],median_Module_TempCorrelations[[i]][2],median_Module_TempCorrelations[[i]][3]))
}
Med_Mod_TempCor <- as.data.frame(Med_Mod_TempCor)
row.names(Med_Mod_TempCor) <- names(Dictionary_Subset)
Med_Mod_TempCor[,1] <- unlist(Med_Mod_TempCor[,1])
Med_Mod_TempCor[,2] <- unlist(Med_Mod_TempCor[,2])
Med_Mod_TempCor[,3] <- unlist(Med_Mod_TempCor[,3])
Med_Mod_TempCor <- cbind(Med_Mod_TempCor, KO_Per_M_Data_Present/KO_Per_M_Present)

### 
#Med_Mod_AZ_TempCor <- NULL
#for (i in 1:607){
#  Med_Mod_AZ_TempCor <- rbind(Med_Mod_AZ_TempCor, c(med_mod_az_TC[[i]][4],med_mod_az_TC[[i]][2],med_mod_az_TC[[i]][3] ))
#}
#Med_Mod_AZ_TempCor <- as.data.frame(Med_Mod_AZ_TempCor)
#row.names(Med_Mod_AZ_TempCor) <- names(Dictionary_Subset)
#Med_Mod_AZ_TempCor[,1] <- unlist(Med_Mod_AZ_TempCor[,1])
#Med_Mod_AZ_TempCor[,2] <- unlist(Med_Mod_AZ_TempCor[,2])
#Med_Mod_AZ_TempCor[,3] <- unlist(Med_Mod_AZ_TempCor[,3])
#Med_Mod_AZ_TempCor <- cbind(Med_Mod_AZ_TempCor, KO_Per_M_Data_Present/KO_Per_M_Present)
### Counting number of significant modules (188 significant modules out of 607 Modules)
#sum(1*grepl("TRUE",Med_Mod_AZ_TempCor$p.value<0.05))

### Removing Modules less than 50% completeness across dataset (leaves 538 Modules out of 607 Modules, removing 69 modules in total)
#Med_Mod_AZ_TempCor <- Med_Mod_AZ_TempCor[Med_Mod_AZ_TempCor$`KO_Per_M_Data_Present/KO_Per_M_Present`>=0.5,]

# Removing Modules for which Temperature correlation could not be calculated
#Med_Mod_AZ_TempCor <- Med_Mod_AZ_TempCor[complete.cases(Med_Mod_AZ_TempCor),]

### Counting significant Modules (188 Modules with significant temperature correlations) 
#sum(1*(Med_Mod_AZ_TempCor$p.value<0.05))

#Sig_Temp_Cor_Zeroes <- Med_Mod_AZ_TempCor[Med_Mod_AZ_TempCor$p.value <0.05,]

### Finding Modules that are significant if not adding zeroes that become insignifcant when adding zeroes. 
#sum(1*(row.names(Med_Mod_TempCor[Med_Mod_TempCor$p.value<0.05,]) %nin% row.names(Med_Mod_AZ_TempCor[Med_Mod_AZ_TempCor$p.value<0.05,])))

Sig_Temp_Cor <- Med_Mod_TempCor[Med_Mod_TempCor$p.value<0.05,]
#Non_Overlapping <- Sig_Temp_Cor[row.names(Med_Mod_TempCor[Med_Mod_TempCor$p.value<0.05,]) %nin% row.names(Med_Mod_AZ_TempCor[Med_Mod_AZ_TempCor$p.value<0.05,]),]
#hist(Non_Overlapping$`KO_Per_M_Data_Present/KO_Per_M_Present`)
# 17 Modules are signifcant w/o zeroes because they are less than 50% complete
#sum(1*(Non_Overlapping$`KO_Per_M_Data_Present/KO_Per_M_Present`<0.5))
# 33 Modules are signifcant w/o zeroes but insignifant with zeroes despite being >50% complete
#sum(1*(Non_Overlapping$`KO_Per_M_Data_Present/KO_Per_M_Present`>=0.5))
#Weirdos <- row.names(Non_Overlapping[Non_Overlapping$`KO_Per_M_Data_Present/KO_Per_M_Present`>=0.5,])
#mid_mod[Weirdos,]
#Modules[Weirdos,1:2]

# Finding Modules that become significant w/zeroes that are not significant w/o zeroes. 
#Non_Overlapping_Zeroes <- Sig_Temp_Cor_Zeroes[row.names(Sig_Temp_Cor_Zeroes)%nin%row.names(Sig_Temp_Cor),]
#Modules[row.names(Non_Overlapping_Zeroes),1:2]

Med_Mod_Ttest <- NULL
for (i in 1:607){
  Med_Mod_Ttest <- rbind (Med_Mod_Ttest, c(median_Module_t.test[[i]][1],median_Module_t.test[[i]][2],median_Module_t.test[[i]][3]))
}
Med_Mod_Ttest <- as.data.frame(Med_Mod_Ttest)
row.names(Med_Mod_Ttest) <- names(Dictionary_Subset)
Med_Mod_Ttest[,1] <- unlist(Med_Mod_Ttest[,1])
Med_Mod_Ttest[,2] <- unlist(Med_Mod_Ttest[,2])
Med_Mod_Ttest[,3] <- unlist(Med_Mod_Ttest[,3])

#Med_Mod_AZ_Ttest <- NULL
#for(i in 1:607){
 # Med_Mod_AZ_Ttest <- rbind(Med_Mod_AZ_Ttest, c(med_mod_az_t.test[[i]][1],med_mod_az_t.test[[i]][2],med_mod_az_t.test[[i]][3]))
#}
#Med_Mod_AZ_Ttest <- as.data.frame(Med_Mod_AZ_Ttest)
#row.names(Med_Mod_AZ_Ttest) <- names(Dictionary_Subset)
#Med_Mod_AZ_Ttest[,1] <- unlist(Med_Mod_AZ_Ttest[,1])
#Med_Mod_AZ_Ttest[,2] <- unlist(Med_Mod_AZ_Ttest[,2])
#Med_Mod_AZ_Ttest[,3] <- unlist(Med_Mod_AZ_Ttest[,3])

Sig_Ttest_Modules <- Modules_Subset[Med_Mod_Ttest[,3]<0.05,]
Sig_TempCor_Modules <- Modules_Subset[Med_Mod_TempCor[,3]<0.05,]

Combined_Sig_Module_Names <- unique(c(row.names(Sig_Ttest_Modules), row.names(Sig_TempCor_Modules)))

Combined_Sig_Modules <- mid_mod[Combined_Sig_Module_Names,]


TempCor_Sig_vector <-1*(Combined_Sig_Module_Names%nin%row.names(Sig_TempCor_Modules)=="FALSE")
Ttest_Sig_vector <- 2*(Combined_Sig_Module_Names%nin%row.names(Sig_Ttest_Modules)=="FALSE")

Sig_Code <- TempCor_Sig_vector + Ttest_Sig_vector
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

Combined_Sig_Modules.zs <- decostand(Combined_Sig_Modules, method="standardize", MARGIN=1)
par(mfrow=c(1,1))
setEPS()

png("../Figures/Heatmap_Combined.png", width = 2000, height=4000, pointsize=8)
fig10 <- heatmap.2(Combined_Sig_Modules.zs, col=hc(100), key=FALSE, symkey=FALSE, trace="none", density.info="none", colsep=c(1:12),rowsep=c(1:nrow(Combined_Sig_Modules.zs)), sepcolor="black", sepwidth=c(0.01,0.00001), dendrogram="row",cexRow = 2, labRow=FALSE, margins=c(5,13), srtCol=90, colRow = Classification_Test, RowSideColors=Classification_Test, lhei = c(1,100))
dev.off()

# Investigating the Clusters
Temp_Sensitive <- c("M00591","M00644", "M00339","M00480","M00474","M00298","M00251","M00090","M00235","M00204","M00722","M00216","M00230","M00607","M00662","M00446","M00040","M00208","M00042","M00641","M00672","M00627","M00509","M00639","M00136","M00025","M00328","M00593","M00631")
Modules[Temp_Sensitive,]
write.table(Modules[Temp_Sensitive,1:2], file="TempSensitiveModules.txt", quote=FALSE)

C1 <- c("M00597", "M00435","M00011","M00320","M00638","M00249","M00153","M00163","M00609","M00161","M00145","M00346","M00061","M00580","M00345","M00172","M00171","M00027")
Modules[C1,]
write.table(Modules[C1,1:2], file="Cluster1_Modules.txt", quote=FALSE)

C2 <- c("M00479","M00116","M00196","M00343","M00342","M00309","M00482","M00113","M00256","M00243","M00283","M00550","M00151","M00540","M00416","M00548","M00604","M00623","M00203","M00219","M00001","M00167","M00033","M00081","M00162","M00545","M00569","M00436","M00510","M00010","M00176","M00165","M00612")
Modules[C2,]
write.table(Modules[C2,1:2], file="Cluster2_Modules.txt", quote=FALSE)

C3 <- c("M00668","M00080","M00543","M00091","M00331","M00613","M00099","M00186","M00250","M00190","M00602","M00544","M00135")
Modules[C3,]
write.table(Modules[C3,1:2], file="Cluster3_Modules.txt", quote=FALSE)

C4 <- c("M00400","M00403","M00423","M00413","M00179","M00184","M00391","M00390","M00288","M00177","M00425","M00181","M00180","M00182","M00529","M00120","M00095","M00344","M00530","M00367","M00365","M00166","M00374", "M00032","M00401","M00072","M00264","M00261","M00633","M00596","M00763","M00031","M00159","M00582","M00338")
Modules[C4,]
write.table(Modules[C4,1:2], file="Cluster4_Modules.txt", quote=FALSE)

C5 <- c("M00013", "M00717","M00442","M00290","M00473")
Modules[C5,]
write.table(Modules[C5,1:2], file="Cluster5_Modules.txt", quote=FALSE)

C6<- c("M00168","M00740", "M00368","M00018","M00012","M00052","M00178","M00023","M00019","M00570","M00140","M00222","M00360","M00359","M00378","M00237","M00053","M00016","M00114","M00254","M00141","M00364","M00003","M00005","M00045","M00239","M00207","M00170","M00489","M00490","M00611", "M00392","M00267","M00275","M00507","M00552","M00614","M00051","M00366","M00007","M00121","M00002", "M00149", "M00021", "M00029")
Modules[C6,]
write.table(Modules[C6,1:2], file="Cluster6_Modules.txt", quote=FALSE)

C7 <- c("M00133","M00620","M00616","M00484","M00087","M00088","M00375")
Modules[C7,]
write.table(Modules[C7,1:2], file="Cluster7_Modules.txt", quote=FALSE)

C8 <- c("M00189","M00471","M00155","M00376","M00483","M00434","M00236","M00193","M00028","M00036","M00549","M00506","M00260","M00454","M00505","M00246","M00245","M00093","M00060","M00572")
Modules[C8,]
write.table(Modules[C8,1:2], file="Cluster8_Modules.txt", quote=FALSE)

C9 <- c("M00008","M00066","M00333","M00308","M00452","M00554","M00632","M00173","M00535","M00432","M00082","M00006","M00361", "M00004", "M00335")
Modules[C9,]
write.table(Modules[C9,1:2], file="Cluster9_Modules.txt", quote=FALSE)

C10 <- c("M00125","M00022","M00049","M00048","M00026","M00157","M00336","M00188","M00096","M00127","M00131","M00209","M00526","M00525","M00527","M00531","M00117","M00064","M00009","M00259","M00144","M00729", "M00115", "M00183")
Modules[C10,]
write.table(Modules[C10,1:2], file="Cluster10_Modules.txt", quote=FALSE)

C11 <- c("M00129","M00394","M00122","M00119","M00362","M00670","M00669")
Modules[C11,]
write.table(Modules[C11,1:2], file="Cluster11_Modules.txt", quote=FALSE)

C12 <- c("M00567","M00465","M00647","M00725","M00242","M00358","M00458","M00481")
Modules[C12,]
write.table(Modules[C12,1:2], file="Cluster12_Modules.txt", quote=FALSE)

C13<- c("M00357","M00077","M00450","M00646","M00035","M00460")
Modules[C13,]
write.table(Modules[C13,1:2], file="Cluster13_Modules.txt", quote=FALSE)

C14 <- c("M00124","M00050","M00565","M00210","M00240","M00519","M00741","M00123","M00577","M00573")
Modules[C14,]
write.table(Modules[C14,1:2], file="Cluster14_Modules.txt", quote=FALSE)


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
postscript("../Figures/Figure5.eps", width = 12, height=4, pointsize=8,paper="special")
par(mfrow=c(1,3))
plot(Module.pcoa$points[,1], Module.pcoa$points[,2],cex=1.5, bg=class, pch=21, main= "rpoB Relativized Bray Curtis KEGG Module PCoA", xlab= paste("PCoA1: ",100*round(M_ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(M_ax2.v,3),"% var. explained",sep=""))
textxy(X=Module.pcoa$points[,1],Y=Module.pcoa$points[,2], lab=map_MG$Sample,cex=1)
M_env<- envfit(Module.pcoa,env)
plot(M_env, p.max=0.05, col="black")

plot(KO.pcoa$points[,1], KO.pcoa$points[,2],cex=1.5, bg=class, pch=21, main= "rpoB Relativized Bray Curtis KEGG Ortholog PCoA", xlab= paste("PCoA1: ",100*round(KO_ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(KO_ax2.v,3),"% var. explained",sep=""))
textxy(X=KO.pcoa$points[,1], Y=KO.pcoa$points[,2],labs=map_MG$Sample, cex=1)

KO_env=envfit(KO.pcoa, env)
plot(KO_env, p.max=0.05, col="black")

plot(uf.pcoa$points[,1],uf.pcoa$points[,2] ,cex=1.5,pch=21,bg=class,main="Weighted UniFrac PCoA",xlab= paste("PCoA1: ",100*round(uf_ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(uf_ax2.v,3),"% var. explained",sep=""))
#textxy is from the calibrate library
textxy(X=uf.pcoa$points[,1], Y=uf.pcoa$points[,2],labs=map_MG$Sample, cex=1)
uf_env <- envfit(uf.pcoa, env)
plot(uf_env, p.max=0.05, col="black")
dev.off()