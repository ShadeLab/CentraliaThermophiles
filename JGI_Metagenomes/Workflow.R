## Workflow Script for Sorensen_INPREP
library(ggplot2)
library(vegan)
library(outliers)
library(gplots)
library(colorRamps)
library(dplyr)
### Reading in Data Files and Manipulating Datafile
setwd("~/GitHub_Repos/ShadeLab/CentraliaThermophiles/JGI_Metagenomes/")
# Mapping File
map <- read.table("Input_Files/Centralia_Collapsed_Map_forR.txt",sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)

map_MG <- map[c(1,3,4,5,6,7,10,12,14,15,16,17),]

sum_stats <- read.table("Input_Files/Centralia_Metagenome_Summary_Stats.txt", sep="\t", header=TRUE, row.names=1)
sum_stats <- sum_stats[,-10]
sum_stats$PercentMapped <- sum_stats$Aligned.Reads/sum_stats$Quality.Reads

sum_stats$Sample <- map_MG$Sample

CO2_2015 <- c(525,455,500,512,484,779,581,503,6094,17112,15760,13969,7390,723,1624,597,480,477)
CO2_MG_2015 <- CO2_2015[c(1,3,4,5,6,7,10,12,14,15,16,17)]
t.test(CO2_MG_2015[map_MG$Classification=="FireAffected"], CO2_MG_2015[map_MG$Classification!="FireAffected"])

# Alpha Diversity
alpha <- read.table("Input_Files/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_alphadiv.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names = 1)

alpha <- alpha[order(row.names(alpha)),]
alpha_MG <- alpha[c(1,3,4,5,6,7,10,12,14,15,16,17),]

setEPS()
postscript("Figures/MGStatsVs16Sstats.eps", width = 5, height=5, pointsize=8,paper="special")
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
comm <- read.table("Input_Files/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)
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
KO <- read.table("Input_Files/4-13-17abundance_ko_126107.tab.txt", sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)
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
outl<-NULL
for (i in 1:12){
  outl <- c(outl, row.names(MO[MO[,i]==outlier(MO[,i]),])) 
}
#MO <- melt(MO, id.vars="row.names(Odds_Ratio)",variable.name="Sample", value.name="Measurement")
# Remove K01519 as it is an outlier in 10/12 samples
y<-MO[MO$C12!=max(MO$C12),]
y <- melt(y, id.vars="row.names(Odds_Ratio)",variable.name="Sample", value.name="Measurement")
colnames(y) <- c("KO", "Sample", "Measurement")
library(dplyr)
?inner_join

fuller_map <- inner_join(map_MG, sum_stats, by="Sample")
Joined_Data <- inner_join(y, fuller_map, by="Sample")
colnames(Joined_Data)[10] <- "Temperature"
ggplot(Joined_Data, aes(x=Temperature, y=PercentMapped, color=KO)) + geom_point() +geom_smooth(method="lm", alpha=0) + theme_bw(base_size=12)

z <- cor.test(Joined_Data$Temperature, Joined_Data$Measurement)
odr <- ggplot(Joined_Data, aes(x=Temperature, y=Measurement)) + geom_point(size=1.5)  + guides(color=FALSE) + theme_bw(base_size=12) + theme(text=element_text(size=8), axis.text = element_text(size=8)) +labs(x=expression("Temperature " ( degree~C)), y="Odds Ratio") + annotate("text", x=44, y=1.4, label=paste("Pearson's r =",round(as.numeric(z[4]),3)), size=2.9) + scale_y_continuous(limits=c(.5,1.5))
odr
#Code for adding a trendline 
#geom_smooth(aes(x=Temperature, y=Measurement), method="lm", alpha=0, group=1, colour="black")
cor.test(Joined_Data$Temperature, Joined_Data$Measurement)
#odr <- ggplot(Joined_Data, aes(x=SoilTemperature_to10cm, y=Measurement, color=KO)) + geom_point() + geom_smooth(method="lm", alpha=0) + facet_wrap(~KO, nrow=6)


JGI_Archaea <- read.table("Input_Files/JGI_Archaea_06192017.txt", sep="\t", header=TRUE, row.names=NULL, stringsAsFactors = FALSE, fill=TRUE, quote="")

Permanent_Archaea <- JGI_Archaea[grepl("Permanent", JGI_Archaea$Status),]
Finished_Archaea <- JGI_Archaea[grepl("Finished",JGI_Archaea$Status),]


JGI_Bacteria <- read.table("Input_Files/JGI_Bacteria_06192017.txt", sep="\t", header=TRUE, row.names=NULL, stringsAsFactors = FALSE, fill=TRUE, quote="")

Permanent_Bacteria <- JGI_Bacteria[grepl("Permanent", JGI_Bacteria$Status),]
Finished_Bacteria <- JGI_Bacteria[grepl("Finished", JGI_Bacteria$Status),]


PF_Archaea <- rbind(Permanent_Archaea, Finished_Archaea)
PF_Bacteria<- rbind(Permanent_Bacteria,Finished_Bacteria)

PF_All <- rbind(PF_Archaea, PF_Bacteria)

rdp <- rdp[rowSums(comm.sigs)>1]
rdp <- gsub("\\[|\\]","", rdp)

Taxonomy <- NULL
for(i in 1:length(rdp)){
  Taxonomy <- rbind(Taxonomy, unlist(strsplit(as.character(rdp[i]), ";")))
}

Phylogeny <- Taxonomy[,2]
Phyla <- Taxonomy[,2]
Phyla<- gsub(" p__", "", Phyla)
Phyla <- gsub("k__", "", Phyla)
Phyla <- unique(Phyla)
Phyla
Phyla[which(!Phyla%in%PF_All$Phylum)]

Phylogeny <- Taxonomy[,2]
# Bacteria
# TM6, WS3, OD1, WPS-2, FCPU426, Parvarchaeota, AD3, OP11, GOUTA4, OP3, GN04, GAL15, TM7, Unclassified, SBR1093, NKB19, FBP, GN02, WWE1, Archaea, NC10, BHI80-139, WS4, WS5, Caldithrix, Thermi, OctSpA1-106, MVP-21, SC4, WS2, "Blank", SAR406, PAUC34f, LD1, AC1, "MVS-104, SR1 

Not_Present <- c("AC1", "AD3", "FCPU426", "GN02", "GN04", "GOUTA4", "LD1", "MVP-21", "MVS-104", "PAUC34f", "SAR406", "SBR1093", "SC4", "WS2","WS3", "WS4", "WS5", "WWE1")
for(i in 1:length(Not_Present)){
  Phylogeny<-gsub(Not_Present[i], "Bacteria", Phylogeny )
}

GG_to_JGI <- read.table("Input_Files/GG_to_JGI.txt", header=FALSE, row.names=NULL, sep="\t", stringsAsFactors = FALSE)

for(i in 1:nrow(GG_to_JGI)){
  Phylogeny <- gsub(GG_to_JGI[i,1], GG_to_JGI[i,2], Phylogeny)
}

Phylogeny<- gsub(" p__", "", Phylogeny)
Phylogeny <- gsub("k__", "", Phylogeny)
v <- unique(Phylogeny)
v<- v[-42]
Fixed_Phylum <- unique(Phylogeny)
Fixed_Phylum[22] <- "Candidatus Parvarchaeota"
Fixed_Phylum <- Fixed_Phylum[complete.cases(Fixed_Phylum)]
Fixed_Phylum <- Fixed_Phylum[-42]  
library(outliers)

# For each phylum calculate
# Total,Average, IQR(Q3-Q1), Q1, Q3, Lower inner fence(Q1-1.5*IQR), Upper inner fence(Q3+1.5*IQR), Outliers, New Mean, SD  

Necessary_JGI <- NULL
for (i in 1:length(Fixed_Phylum)){
  x <- PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14]
  Total <- as.numeric(length(x))
  avg_bad <- mean(x) 
  Q1<- as.numeric(quantile(x)[2])
  Q3 <- as.numeric(quantile(x)[4])
  IQR <- Q3-Q1
  LIF <- Q1-(1.5*IQR)
  UIF <- Q3 + (1.5*IQR)
  O <- Total-sum(1*(LIF<x & x<UIF))
  x_good <- x[LIF<x & x<UIF]
  avg_good <- mean(x_good)
  sd_good <- sd(x_good)
  Necessary_JGI <- rbind(Necessary_JGI, c(Total, avg_bad, Q1, Q3, IQR, LIF, UIF, O, avg_good, sd_good, sd_good/avg_good))
}
row.names(Necessary_JGI)<- Fixed_Phylum


plot_data<- NULL
for (i in 1:length(Fixed_Phylum)){
  x <- cbind(rep(Fixed_Phylum[i], length(PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14])), PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14])
  plot_data <- rbind(plot_data, x)             
}
plot_data <- as.data.frame(plot_data)

#ggplot(plot_data, aes(V1, V2)) + geom_boxplot() + facet_wrap(~factor(V1), nrow=6)

x <- PF_All[grepl(Fixed_Phylum[i], PF_All$Phylum),14]

#creating Phylum level table for Centralia Data
comm_12<-comm[,c(1,3,4,5,6,7,10,12,14,15,16,17)]
Phyla_Comm<- NULL
for(i in 1:length(Fixed_Phylum)){
  x <- comm_12[grepl(Fixed_Phylum[i],Phylogeny),]
  y<- as.numeric(colSums(x))
  Phyla_Comm <- rbind(Phyla_Comm, y)
}
row.names(Phyla_Comm) <- Fixed_Phylum
Phyla_Comm_Rel <- decostand(Phyla_Comm, method= "total", MARGIN=2 )
colSums(Phyla_Comm_Rel)



Phyla_Comm_Whole<- NULL
for(i in 1:length(Fixed_Phylum)){
  x <- comm[grepl(Fixed_Phylum[i],Phylogeny),]
  y<- as.numeric(colSums(x))
  Phyla_Comm_Whole <- rbind(Phyla_Comm_Whole, y)
}
row.names(Phyla_Comm_Whole) <- Fixed_Phylum
Phyla_Comm_Whole_Rel <- decostand(Phyla_Comm_Whole, method= "total", MARGIN=2 )
colSums(Phyla_Comm_Whole_Rel)

sum(1*(PF_Bacteria$Genome.Size.....assembled<100000))
nrow(PF_Bacteria)
Necessary_JGI["Bacteria", 9] <- mean(PF_Bacteria$Genome.Size.....assembled)
Necessary_JGI["Unclassified", 9] <- mean(PF_All$Genome.Size.....assembled)
Necessary_JGI["OP11", 9] <- mean(PF_Bacteria$Genome.Size.....assembled)
Necessary_JGI["FBP", 9] <- mean(PF_Bacteria$Genome.Size.....assembled)
Necessary_JGI["Archaea", 9] <- mean(PF_Archaea$Genome.Size.....assembled)
Necessary_JGI["candidate division TM6", 9] <- Necessary_JGI["candidate division TM6",2]
Necessary_JGI["WS1", 9] <- Necessary_JGI["WS1",2]
Necessary_JGI["candidate division SR1", 9] <- Necessary_JGI["candidate division SR1",2]


Output2 <- NULL
for(i in 1:ncol(Phyla_Comm_Rel)){
  x <- Phyla_Comm_Rel[,i]*Necessary_JGI[,9]
  z <- sum(x)
  Output2 <- c(Output2, z)
}

Output2_Whole <- NULL
for(i in 1:ncol(Phyla_Comm_Whole_Rel)){
  x <- Phyla_Comm_Whole_Rel[,i]*Necessary_JGI[,9]
  z <- sum(x)
  Output2_Whole <- c(Output2_Whole, z)
}


plot(map_MG$SoilTemperature_to10cm, Output2)
cor.test(map_MG$SoilTemperature_to10cm, Output2)
plot_data_2 <- cbind(Output2, map_MG$SoilTemperature_to10cm)
plot_data_2 <- as.data.frame(plot_data_2)
colnames(plot_data_2) <- c("Estimate", "Temperature")
z <- cor.test(plot_data_2$Temperature, plot_data_2$Estimate)
rRNA_Size_Estimate <- ggplot(plot_data_2, aes(x=Temperature, y =Estimate )) + geom_point(size=1.5) + theme_bw(base_size=8) + theme(text=element_text(size=8), axis.text = element_text(size=8)) +labs(x=expression("Temperature " ( degree~C)), y="Average Genome Size (bp)") + annotate("text", x=45, y=3600000, label=paste("Pearson's r =",round(as.numeric(z[4]),3)), size=2.9) + scale_y_continuous(labels = function(y) format(y, scientific = TRUE), limits=c(3000000, 3800000))
rRNA_Size_Estimate

M_Census <- c(6705016, 6046631, 6174721, 6213575, 6088619, 5810701, 4136147, 5295248, 4506566, 4881807, 5421668, 6144962)
plot_data_3 <- cbind(M_Census, map_MG$SoilTemperature_to10cm)
plot_data_3 <- as.data.frame(plot_data_3)
colnames(plot_data_3) <- c("Estimate", "Temperature")

z <- cor.test(plot_data_3$Temperature, plot_data_3$Estimate)
MCensus_Size_Estimate <- ggplot(plot_data_3, aes(x=Temperature, y =Estimate )) + geom_point(size=1.5) + theme_bw(base_size=8) + theme(text=element_text(size=8), axis.text = element_text(size=8)) +labs(x=expression("Temperature " ( degree~C)), y="Average Genome Size (bp)")+ annotate("text", x=45, y=6000000, label=paste("Pearson's r =",round(as.numeric(z[4]),3)), size=2.9) + scale_y_continuous(labels = function(y) format(y, scientific = TRUE),limits=c(4000000, 7000000))
MCensus_Size_Estimate
Dummy_plot <- plot.new()

### Multiplot code taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

mp <- multiplot(odr, rRNA_Size_Estimate, Dummy_plot, MCensus_Size_Estimate, cols=2)
setEPS()
postscript("Figures/Figure1_AverageGenomeSize.eps",  width= 3, height=6)
par(ps = 8, cex = 1, cex.main = 1)
mp <- multiplot(odr, MCensus_Size_Estimate,rRNA_Size_Estimate, cols=1)
dev.off()


### KEGG Ortholog Eveness
s=specnumber(t(KO))
h=diversity(t(KO), index="shannon")
pielou = h/log(s)

### KEGG Module Analysis
library(stringr)
Modules <- read.table("Input_Files/kmodlist47982_23-nov-2016.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

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

# Relativized to Average Single Copy Gene count
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
uf=read.table("Input_Files/weighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", header=TRUE, row.names=1)
uf=uf[order(row.names(uf)),order(colnames(uf))]
# Subsetting to only the samples we have metagenomes for...
uf_mg <- uf[c(1,3,4,5,6,7,10,12,14,15,16,17),c(1,3,4,5,6,7,10,12,14,15,16,17)]
uf_mg.d <- as.dist(uf_mg)

# JGI Enzyme, pfam, COG datatables
enzyme <- read.table("Input_Files/abundance_enzyme_79244.txt", sep="\t", header=TRUE, row.names = 1, stringsAsFactors = FALSE)
colnames(enzyme) <- map_MG$Sample
enzyme.rel <- decostand(enzyme, MARGIN = 2, method="total")
enzyme.d <- vegdist(t(enzyme.rel), method = "bray")

cog <- read.table("Input_Files/abundance_cog_79133.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
colnames(cog) <- map_MG$Sample
cog.rel <- decostand(cog, MARGIN=2, method="total")
cog.d<- vegdist(t(cog.rel), method="bray")

pfam <- read.table("Input_Files/abundance_pfam_72927.txt", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
colnames(pfam) <- map_MG$Sample
pfam.rel <- decostand(pfam, MARGIN=2, method="total")
pfam.d <- vegdist(t(pfam.rel), method="bray")


# Distance matrix based on Modules
Mod.d <- vegdist(t(mid_mod),method="bray")
Mod_sub.d <- vegdist(t(mid_mod[,-12]), method="bray")
Class <- rep("Red", 12)
Class[map_MG$Classification!="FireAffected"] <- "Green"
a=adonis(Mod.d~Class, distance=TRUE, permutations=1000)
a

b=betadisper(Mod.d, group=Class)
TukeyHSD(b, which = "group", ordered = FALSE,conf.level = 0.95)

b.sub= betadisper(Mod_sub.d, group=Class[-12])
TukeyHSD(b.sub, which="group", ordered=FALSE, conf.level = 0.95)
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

M_v_cog.mantel <- mantel(Mod.d, cog.d)
M_v_enzyme.mantel <- mantel(Mod.d, enzyme.d)
M_v_pfam.mantel <- mantel(Mod.d, pfam.d)

# Database Comparisons 
dbs <- read.table("Input_Files/DB_Comparison.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
colnames(dbs) <- c("COG", "pfam", "Enzyme", "KO")
dbs$Site <- row.names(dbs)
dbs_melted <- melt(dbs, id.vars = "Site")
dbc <- ggplot(dbs_melted, aes(x=variable, y=value)) + geom_boxplot() +labs(x="Database", y="Percentage of Genes")
ggsave("Supplemental/DB_Comparison.png", device = "png", dbc)


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
row.names(Combined_Sig_Modules.zs)

Complete <- Combined_Sig_Module_Results[Combined_Sig_Module_Results$Completeness>=0.5,]
Complete_Combined_Sig_Modules <- Combined_Sig_Modules[Combined_Sig_Module_Results$Completeness>=0.5,]


Pos_Temp_Cor_Modules <- Complete_Combined_Sig_Modules[Complete$estimate>0,]
Pos_Temp_Cor_Modules <- Pos_Temp_Cor_Modules[complete.cases(Pos_Temp_Cor_Modules),]
#Positively Correlated Modules z-scored
PCM.zs <- decostand(as.matrix(Pos_Temp_Cor_Modules), method="standardize", MARGIN=1)

PCM.zs <- PCM.zs[,order(map_MG$SoilTemperature_to10cm)]

# Negative Temp Cor Modules
Neg_Temp_Cor_Modules <- Complete_Combined_Sig_Modules[Complete$estimate<0,]
Neg_Temp_Cor_Modules <- Neg_Temp_Cor_Modules[complete.cases(Neg_Temp_Cor_Modules),]
# Z-score Negative Correlated Modules
NCM.zs <- decostand(as.matrix(Neg_Temp_Cor_Modules), method="standardize", MARGIN = 1)
NCM.zs <- NCM.zs[,order(map_MG$SoilTemperature_to10cm)]


# Plot and Save PCM Heatmap
setEPS()
png("../Figures/PCM.png", width = 500, height=1000, pointsize=8)
heatmap.2(PCM.zs, col=hc(100), key=TRUE, symkey=TRUE, trace="none", colsep=c(1:12),rowsep=c(1:nrow(PCM.zs)), sepcolor="black", Colv=FALSE, sepwidth=c(0.01,0.00001),labRow = row.names(PCM.zs), dendrogram="row",cexRow = 1, margins=c(5,13), srtCol=90, lhei=c(1,100))
dev.off()

# Plot and Save PCM heatmap Key
setEPS()
postscript("../Figures/PCM_Key.eps", width = 4, height=8)
par(ps = 12, cex = 1, cex.main = 1)
heatmap.2(PCM.zs, col=hc(100), key=TRUE, symkey=TRUE, trace="none",density.info = "none", colsep=c(1:12),rowsep=c(1:nrow(PCM.zs)), sepcolor="black", Colv=FALSE, sepwidth=c(0.01,0.00001),labRow = row.names(PCM.zs), dendrogram="row",cexRow = 1, srtCol=90, lmat= rbind(c(3,4),c(2,1)), lhei=c(1,4))
dev.off()


# Plot and Save NCM heatmap
setEPS()
postscript("../Figures/NCM.eps", width = 4, height=8)
par(ps = 12, cex = 1, cex.main = 1)
heatmap.2(NCM.zs, col=hc(100),key=TRUE, symkey=TRUE, trace="none", colsep=c(1:12),density.info="none", sepcolor="black", Colv=FALSE, sepwidth=c(0.01,0.0000001), dendrogram="row",cexRow = 1, labRow=FALSE, srtCol=90,lmat=rbind(c(3,4),c(2,1)), lhei=c(1,4))
dev.off()
par(mfrow=c(1,1))
# Plot and Sace NCM heatmap Key
setEPS()
png("../Figures/NCM_Key.png", width = 500, height=500, pointsize=8)
heatmap.2(NCM.zs, col=hc(100), key=TRUE, symkey=TRUE, trace="none",density.info = "none", colsep=c(1:12),rowsep=c(1:nrow(NCM.zs)), sepcolor="black", Colv=FALSE, sepwidth=c(0.01,0.00001),labRow = FALSE, dendrogram="row",cexRow = 1, margins=c(5,13), srtCol=90)
dev.off()

Complete <- Complete[complete.cases(Complete$KEGG),] 
row.names(Complete) <- Complete$KEGG

nrow(Complete[Complete$estimate>0,])
nrow(Complete[Complete$estimate<0,])


### Create Supplemental Table of Summaries
SignificantModules_Summary <-join(Med_Mod_TempCor[row.names(Combined_Sig_Modules),], Med_Mod_Ttest[row.names(Combined_Sig_Modules),], by="KEGG")

SignificantModules_Summary$ModuleDescription <- Modules[row.names(Combined_Sig_Modules),1]

row.names(SignificantModules_Summary) <- SignificantModules_Summary$KEGG 

SignificantModules_Summary <- SignificantModules_Summary[row.names(Complete),]

SignificantModules_Summary<-SignificantModules_Summary[,c(4,11,5,1,2,3,6,7,8,9,10)]
colnames(SignificantModules_Summary) <- c("Module","Module Description","Completeness","Pearson's Rho", "Pearson's Degrees Freedom", "Peason p value", "Pearson Adjusted p value", "T Statistic", "T Degrees Freedom", "T-test p value", "T Adjusted p value")
write.table(x = SignificantModules_Summary, file="Supplemental/SupplementalTable2_Median_SCG.txt", sep="\t", quote=FALSE)

PosCorModules <- mid_mod[row.names(SignificantModules_Summary[SignificantModules_Summary$`Pearson's Rho`>0,]),]
NegCorModules <- mid_mod[row.names(SignificantModules_Summary[SignificantModules_Summary$`Pearson's Rho`<0,]),]

P_sd <- apply(PCM.zs[,map_MG$Classification=="FireAffected"],1,sd)
N_sd <- apply(NCM.zs[,map_MG$Classification=="FireAffected"],1,sd)
t.test(P_sd^2, N_sd^2)
t.test(c(P_sd, N_sd)^2~c(rep("Positive", length(P_sd)), rep("Negative", length(N_sd))))




# Percent KO analysis
PKO <- read.table("Centralia_MGKOPercent.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
PKO <- PKO[order(PKO$Genome.Name...Sample.Name),]
PKO$Genome.Name...Sample.Name <- map_MG$Sample

plot(map_MG$SoilTemperature_to10cm,PKO$Not.KO.......both)
cor.test(map_MG$SoilTemperature_to10cm,PKO$Not.KO.......both)
# % of genes not in KO is negatively and significantly correlated with temperature...

# Plotting DeNitrification and Sulfate Reduction Modules against Temperature and Sulfur and Nitrogen Measurements

# Denitrification
DNitrate <- Combined_Sig_Modules["M00529",]

#Dissimilatory Nitrate Reduction
DNR <- Combined_Sig_Modules["M00530",]

# Dissimilalry Sulfate Reduction
DSR <- Combined_Sig_Modules["M00596",]

plot(map_MG$SulfateSulfur_ppm, as.numeric(DSR))

ModuleAnalysis <- rbind(DNitrate, DNR, DSR, map_MG$SoilTemperature_to10cm, map_MG$NO3N_ppm, map_MG$NH4N_ppm, map_MG$SulfateSulfur_ppm)
row.names(ModuleAnalysis) <- c("M00529","M00530", "M00596", "Temperature", "Nitrate", "Ammonium", "Sulfate_Sulfur")
ModuleAnalysis <- t(ModuleAnalysis)
ModuleAnalysis<- as.data.frame(ModuleAnalysis)

ggplot(ModuleAnalysis, aes(x=Temperature, y=M00529)) + geom_point()
ggplot(ModuleAnalysis, aes(x=Temperature, y=M00530)) + geom_point()
ggplot(ModuleAnalysis, aes(x=Temperature, y=M00596)) + geom_point()
ggplot(ModuleAnalysis, aes(x=Nitrate, y=M00529)) + geom_point() 

ggsave("../Figures/DNitrate_Nitrate.eps", width=100, height=75, units="mm")
ggplot(ModuleAnalysis, aes(x=Ammonium, y=M00529)) + geom_point() 

cor.test(x = ModuleAnalysis$Temperature, y=ModuleAnalysis$M00529)
cor.test(x=ModuleAnalysis$Nitrate, y=ModuleAnalysis$M00529)

ggplot(ModuleAnalysis, aes(x=Nitrate, y=M00530)) + geom_point() 
ggsave("../Figures/DNR_Nitrate.eps", width=100, height=75, units="mm")
ggplot(ModuleAnalysis, aes(x=Ammonium, y=M00530)) + geom_point()
cor.test(x=ModuleAnalysis$Nitrate, y=ModuleAnalysis$M00530)
cor.test(x=ModuleAnalysis$Ammonium, y=ModuleAnalysis$M00530)

ggplot(ModuleAnalysis, aes(x=Sulfate_Sulfur, y=M00596)) + geom_point() 
cor.test(x=ModuleAnalysis$Sulfate_Sulfur, y=ModuleAnalysis$M00596)
ggsave("../Figures/DSR.eps", width=100, height=75, units="mm")


# Two-component Regulatory System Results
TCRS <- SignificantModules_Summary[grep("two-component", SignificantModules_Summary$`Module Description`),]
TCRS_M <- Combined_Sig_Modules[TCRS[,1],]
med_TCRS <- apply(TCRS_M, 2, median)

TCRS_M$KM <- row.names(TCRS_M)
TCRS_med_correlation<- NULL
for (i in 1:nrow(TCRS_M)){
  tc <- cor.test(med_TCRS, as.numeric(TCRS_M[i,1:12]))
  TCRS_med_correlation <- rbind(TCRS_med_correlation, c(tc[3],tc[4]))
}

#TCRS_M <- TCRS_M[TCRS_M[,1]<1,]
plot_data_TC <- melt(TCRS_M)
SoilTemp <- data.frame(Temp = map_MG$SoilTemperature_to10cm)  
SoilTemp$Site <- map_MG$Sample
colnames(plot_data_TC) <- c("KM", "Site", "Measurement")
TCRS_M <- join(plot_data_TC, SoilTemp, by = "Site")



TCRS_M$Color <- rep("black",nrow(TCRS_M))

TCRS_plot <- ggplot(TCRS_M, aes(x=Temp, y=Measurement, color=KM)) + geom_smooth(method="lm", inherit.aes=TRUE, se=0) + geom_point(inherit.aes = TRUE) + guides(color=FALSE) + theme_bw(base_size=8) + theme(text=element_text(size=8), axis.text = element_text(size=8)) + labs(x="Temperature (Celsius)", y="Abundance (Copies per Genome)") 
TCRS_plot

ggsave("../Figures/Figure5_TwoComponenetRegulatorySystesms.eps", TCRS_plot, width=4, height=4, units="in")

# Drug resistance Results
DR <- read.table("../Drug_Resistance.txt", stringsAsFactors = FALSE)
DR_Stats <- SignificantModules_Summary[DR[,1],]
DR_M <- Combined_Sig_Modules[DR[,1],]
DR_M$KM <- row.names(DR_M)
DR_M <- DR_M[DR_Stats$Completeness>.5,]
#DR_M <- DR_M[DR_M[,1]<1,]
plot_data_DR <- melt(DR_M)
colnames(plot_data_DR) <- c("KM", "Site", "Measurement")
plot_data_DR <- join(plot_data_DR, SoilTemp, by= "Site")

DR_plot <- ggplot(plot_data_DR, aes(x=Temp, y=Measurement, color=KM)) +geom_point() + geom_smooth(method="lm", se=0) +guides(color=FALSE) + theme_bw(base_size=8) + theme(text=element_text(size=8), axis.text = element_text(size=8)) + labs(x="Temperature (Celsius)", y="Abundance (Copies per Genome)")
DR_plot
ggsave("../Figures/Figure6_DrugResistance_Biosynthesis.eps", DR_plot, width=4, height=4, units="in")


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
postscript("../Figures/Figure5_KEGGModulePCoA.eps", width = 4, height=4, paper="special")
par(ps = 8, cex = 1, cex.main = 1)
plot(Module.pcoa$points[,1], Module.pcoa$points[,2],cex=1, bg=class, pch=21, xlab= paste("PCoA1: ",100*round(M_ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(M_ax2.v,3),"% var. explained",sep=""))
#textxy(X=Module.pcoa$points[,1],Y=Module.pcoa$points[,2], lab=map_MG$Sample,cex=.5)
M_env<- envfit(Module.pcoa,env)
plot(M_env, p.max=0.05, col="black", lwd=3)
dev.off()


# Trying to plot in ggplot2 BROKEN
module_plot_data <- cbind.data.frame(Module.pcoa$points[,1], Module.pcoa$points[,2])
colnames(module_plot_data) <- c("PCoA1", "PCoA2")
module_plot_data$Sample <- map_MG$Sample

module_plot_data <- inner_join(module_plot_data, map_MG, by="Sample")
module_plot_data$Class_Color <- c("Yellow","Yellow","Yellow","Yellow","Red", "Yellow","Red","Red","Red","Red","Red","Green"  )
ggplot(module_plot_data, aes(x=PCoA1, y=PCoA2, color=Class_Color)) + geom_point(colour=module_plot_data$Class_Color, size=3) +geom_text(aes(label=module_plot_data$Sample), hjust=1, vjust=-1)


#textxy(X=Module.pcoa$points[,1],Y=Module.pcoa$points[,2], lab=map_MG$Sample,cex=0)

# https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2
M_env.df <- as.data.frame(M_env$vectors$arrows*sqrt(M_env$vectors$r))
M_env.df$Metadata <- row.names(M_env.df) 

M_env.df <- M_env.df[M_env$vectors$pvals<0.05,]

ggplot(module_plot_data, aes(x=PCoA1, y=PCoA2, color=Class_Color)) + geom_point(colour=module_plot_data$Class_Color, size=3) +geom_text(aes(label=module_plot_data$Sample), hjust=1, vjust=-1) + geom_segment(data=M_env.df,aes(x=0,xend=Dim1,y=0,yend=Dim2), inherit.aes = FALSE,
             arrow = arrow(length = unit(0.25, "cm")),colour="grey",inherit_aes=FALSE) + 
  geom_text(data=M_env.df,aes(x=Dim1,y=Dim2,label=Metadata),size=5, inherit.aes = FALSE)

vec.sp<-envfit(sol$points, NMDS.log, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)



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





# Presence Absence Module Table
mod_PA <- 1*(mid_mod>0)
colSums(mod_PA)
mod_FA <- rowSums(mid_mod[,map_MG$Classification=="FireAffected"])
mod_RecR <- rowSums(mid_mod[,map_MG$Classification!="FireAffected"])
mod_FA_PA <- 1*(mod_FA>0)
mod_RecR_PA <- 1*(mod_RecR>0)

### Indicator Module Analysis
library(indicspecies)
IndicClusters <- rep(1,12)
IndicClusters[map_MG$Classification!="FireAffected"]=2
tmidmod<- t(mid_mod)
tmidmod <- as.data.frame(tmidmod)
B=strassoc(tmidmod, cluster=IndicClusters,func="B")
sel=which(B[,1]>0.2)

indicators(tmidmod[,sel], cluster=IndicClusters, group="1",At=0.5,Bt=0.2, verbose=TRUE)
wetpt = multipatt(tmidmod, IndicClusters, control=how(nperm=999))
summary(wetpt)
