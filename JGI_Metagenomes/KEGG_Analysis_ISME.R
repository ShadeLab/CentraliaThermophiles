setwd("~/GitHub_Repos/ShadeLab/CentraliaThermophiles/JGI_Metagenomes/")
KEGG<- read.table("Gradient_KEGG.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)### Read in KEGG Table from JGI
colnames(KEGG)<-c("C01","C03","C04","C05","C06","C07", "C10", "C12", "C14", "C15", "C16", "C17")### Replace the headers from JGI with Sample names.

map=read.table("Centralia_Collapsed_Map_forR.txt", header=TRUE, sep="\t")
row.names(map) <- map$Sample
map_sub <- map[colnames(KEGG),]


rpoB <- KEGG["KO:K03043",] + KEGG["KO:K13798",] ### add the bacterial and archaeal copies together.
KEGG <- KEGG[-grep("KO:K03043", row.names(KEGG)),]
KEGG <- KEGG[-grep("KO:K13798",row.names(KEGG)),]

KEGG.rpob <- NULL

for (i in 1:nrow(KEGG)){
  KEGG.rpob <- rbind(KEGG.rpob, KEGG[i,]/rpoB)
}


class <- rep("Black", nrow(map_sub))
class[map_sub$Classification=="Recovered"]='yellow'
class[map_sub$Classification=="FireAffected"]='red'
class[map_sub$Classification=="Reference"]='green'

library(calibrate)
library(ggplot2)
library(vegan)
### PCoA Presence Absence 
KEGG.pa <- 1*(KEGG.rpob>0)
K.dist.pa <- vegdist(t(KEGG.pa), method="bray")
K.pa.pcoa<- cmdscale(K.dist.pa, eig=TRUE)
ax1.pa.v=K.pa.pcoa$eig[1]/sum(K.pa.pcoa$eig)
ax2.pa.v=K.pa.pcoa$eig[2]/sum(K.pa.pcoa$eig)
?plot
dev.off()
setEPS()
postscript("KO_PcoA_PA.eps", width = 6.770, height=3.385, pointsize=8,paper="special")
plot(K.pa.pcoa$points[,1], K.pa.pcoa$points[,2],cex=1.5, bg=class, pch=21,main= "Presence/Absence Bray Curtis KEGG PCoA",xlab= paste("PCoA1: ",100*round(ax1.pa.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(ax2.pa.v,3),"% var. explained",sep=""))
###textxy(X=K.pa.pcoa$points[,1], Y=K.pa.pcoa$points[,2],labs=map_sub$Sample, cex=1)

env=map_sub[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]
envEF=envfit(K.pa.pcoa, env)
plot(envEF, p.max=0.05, col="black")
dev.off()

### PCoA Relativized to RPOB
K.dist <- vegdist(t(KEGG.rpob), method="bray" )
K.dist
k.pcoa <-cmdscale(K.dist, eig=TRUE)

ax1.v=k.pcoa$eig[1]/sum(k.pcoa$eig)
ax2.v=k.pcoa$eig[2]/sum(k.pcoa$eig)
?plot
dev.off()
setEPS()
postscript("KO_PcoA_Relativized.eps", width = 6.770, height=3.385, pointsize=8,paper="special")
plot(k.pcoa$points[,1], k.pcoa$points[,2],cex=1.5, bg=class, pch=21, main= "rpoB Relativized Bray Curtis KEGG PCoA", xlab= paste("PCoA1: ",100*round(ax1.v,3),"% var. explained",sep=""), ylab= paste("PCoA2: ",100* round(ax2.v,3),"% var. explained",sep=""))
###textxy(X=k.pcoa$points[,1], Y=k.pcoa$points[,2],labs=map_sub$Sample, cex=1)

env=map_sub[,c("SoilTemperature_to10cm", "NO3N_ppm", "pH", "K_ppm", "Mg_ppm", "OrganicMatter_500", "NH4N_ppm", "SulfateSulfur_ppm", "Ca_ppm", "Fe_ppm", "As_ppm", "P_ppm", "SoilMoisture_Per","Fire_history")]
envEF=envfit(k.pcoa, env)
plot(envEF, p.max=0.05, col="black")
dev.off()

Class2=sub("green", "yellow", class)
### rpoB relativized hypothesis testing
a=adonis(K.dist~Class2, distance=TRUE, permutations=1000)
a

### PA hypothesis testing
a.pa=adonis(K.dist.pa~Class2, distance=TRUE, permutations=1000)
a.pa

#multivariate dispersion with Tukey HSD
### rpoB relativized
b=betadisper(K.dist, group=Class2)
TukeyHSD(b, which = "group", ordered = FALSE,conf.level = 0.95)
### PA 
b.pa=betadisper(K.dist.pa, group=Class2)
TukeyHSD(b.pa, which = "group", ordered = FALSE,conf.level = 0.95)

### Alpha Diversity
# Boxplot for KO content
library(reshape2)
richness <- as.numeric(colSums(KEGG.pa))
alpha <- cbind(colnames(KEGG.pa),Class2)
alpha<- as.data.frame(alpha)
alpha$Richness <- richness
str(alpha$Richness)
colnames(alpha) <- c("SampleID", "Classification", "Richness")
#reshape the data
alpha.long=melt(alpha, id.vars=c("SampleID", "Classification"))
#plot a facet
#comment toggle for color v. bw
colors=c("red", "yellow")

fig3 <- ggplot(data=alpha.long, aes(x=Classification, y=value)) +
  geom_boxplot() + 
  geom_jitter()+
  #geom_jitter(aes(shape=Classification))+
  #geom_jitter(aes(color=Classification, cex=2))+
  #scale_shape(guide=FALSE)+
  scale_size(guide=FALSE)+
  #scale_color_manual(values=colors)+
  scale_x_discrete(name="Fire classification")+
  theme_bw(base_size=10)
fig3
###ggsave("KO_Richness.eps", width=86, units="mm")
ggsave("KO_Richness.png", width=86, units="mm")
active <- alpha[alpha$Classification=="red",colnames(alpha)=="Richness"]
recref <- alpha[alpha$Classification=="yellow",colnames(alpha)=="Richness"]

test1=t.test(as.numeric(active), as.numeric(recref), paired=FALSE, var.equal = FALSE)
### Fireaffected versus nonfireaffected sites don't vary in the # of KOs present. 

### Indicator "Species" Analysis
library(indicspecies)

KEGG.rpob.t=as.data.frame(t(KEGG.rpob))
class.ind=as.vector(map_sub$Classification)
ind=multipatt(KEGG.rpob.t, class.ind, control=how(nperm=999), duleg = TRUE, func="IndVal.g")
summary(ind, alpha=0.001, indvalcomp=TRUE, At=0.95, Bt=0.95)



### Determining Shared KEGGS between Fire Affected and NonFireAffected sites.

sum(1*rowSums(KEGG.pa[,Class2=="yellow"])>0)
  ### 7216 KOs in Non Fire Affected Sites

sum(1*rowSums(KEGG.pa[,Class2=="red"])>0)
  ### 7236 KOs in Fire Affected Sites 

dim(KEGG.pa[(1*rowSums(KEGG.pa[,Class2=="red"])>0)!=(1*rowSums(KEGG.pa[,Class2=="yellow"])>0),])
  ### 424 KEGGS that are either only in FA sites or only in Non FA sites

Exclusive <- KEGG.pa[(1*rowSums(KEGG.pa[,Class2=="red"])>0)!=(1*rowSums(KEGG.pa[,Class2=="yellow"])>0),]
FA_Exclusive <- Exclusive[rowSums(Exclusive[,Class2=="red"])>0,]
RecRef_Exclusive <- Exclusive[rowSums(Exclusive[,Class2=="yellow"])>0,]

write.table(row.names(FA_Exclusive), "FA_Exclusive_KOs.txt", sep="\t")                 
