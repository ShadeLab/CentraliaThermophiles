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
row.names(subset) <- c("Acidobacteria","Proteobacteria","Verrucomicrobia")
melted_subset <- melt(t(subset))
colnames(melted_subset) <- c("Sample", "Phylum", "Percent_Denovo")
melted_subset$SoilTemperature <- rep(map$SoilTemperature_to10cm,3)


odr <- ggplot(Joined_Data, aes(x=SoilTemperature_to10cm, y=Measurement, color=KO)) + geom_point() + geom_smooth(method="lm", alpha=0) + guides(color=FALSE)
PPDN <- ggplot(melted_subset, aes(x=SoilTemperature, y=Percent_Denovo, shape=Phylum)) + geom_point(size=3) + geom_smooth(colour="black", method="lm", alpha=0, aes(linetype=factor(Phylum))) + scale_linetype_manual("", values=c(1,2,3)) + scale_shape_manual("", values=c(1,2,3)) +ylim(0,0.7) +theme(axis.text= element_text(size=15), axis.ticks = element_line(size=1))
ggsave("../Figures/PPDN.eps", plot=PPDN, width=178, units="mm")

Total <- data.frame(1:4,1:4)

Total$Count <- c(nrow(dn_OTUs),nrow(gg_OTUs),sum(dn_OTUs[,1:18]),sum(gg_OTUs[,1:18]))
Total$DataType <- c("OTUs", "OTUs", "Reads", "Reads")
Total$OTUType <- c("DN","GG", "DN", "GG")
Total$Total <- rep("total",4)

cbPalette <- c("#bdbdbd", "#636363")
ggplot(Total, aes(x=total, y=Count, fill=OTUType)) + geom_bar(stat="identity", aes(x=Total, y=Count)) + facet_wrap(~factor(DataType, levels=c("OTUs","Reads")), nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)),axis.title.y=element_blank()) +  scale_fill_manual(values=cbPalette) + guides(fill=FALSE)
ggplot(Total, aes(x=total, y=Count, fill=OTUType,width=.5)) + geom_bar(stat="identity", aes(x=Total, y=Count)) + facet_wrap(~factor(DataType, levels=c("OTUs","Reads")), nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_blank(), axis.title.x= element_blank(),strip.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank()) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)

ggsave("../Figures/TotalDNvsGG_Abundant.jpg", width=20, units="mm")


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

ggplot(LD_A, aes(x=Phylum, y=Count, fill=Classification)) + geom_bar(stat="identity", aes(x=Phylum, y=Count)) + facet_wrap(~Type, nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=75, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)),axis.text.y=element_blank(),axis.title.y=element_blank()) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)
ggplot(LD_A, aes(x=Phylum, y=Count, fill=Classification)) + geom_bar(stat="identity", aes(x=Phylum, y=Count)) + facet_wrap(~Type, nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_blank(), axis.title.x= element_blank(),strip.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank()) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)
ggsave("../Figures/Phylum_DNvsGG_Abundant.jpg", width=80, units="mm")
ggplot(LD_R, aes(x=Phylum, y=Count, fill=Classification)) + geom_bar(stat="identity", aes(x=Phylum, y=Count)) + facet_wrap(~Type, nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_blank(), axis.title.x= element_blank(),strip.text.x=element_text(size=rel(2)), axis.text.y=element_blank(),axis.title.y=element_blank()) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)
ggsave("../Figures/Phylum_DNvsGG_Rare.jpg", width=210, units="mm")

ggplot(LD_R, aes(x=Phylum, y=Count, fill=Classification)) + geom_bar(stat="identity", aes(x=Phylum, y=Count)) + facet_wrap(~Type, nrow=2, strip.position="right", scales ="free_y") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=75, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_blank(),axis.title.y=element_blank()) + scale_fill_manual(values=cbPalette) + guides(fill=FALSE)




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
  #results.k <- kruskal.test(comm_therm.mat[i,]~ThermSites)
  #spearman <- cor.test(comm_therm.mat[i,], Temp_Full , method="spearman")
  #kendall <- cor.test(comm_therm.mat[i,], Temp_Full, method="kendall")
  coretest.out=rbind(coretest.out,c(row.names(comm_therm.mat)[i],results$estimate,results$p.value))
  ttest.out = rbind(ttest.out, c(row.names(comm_therm.mat)[i], results.t[1], results.t[2], results.t[3]))
  #kendall.out <- rbind(kendall.out, c(row.names(comm_therm.mat)[i], kendall[1], kendall[2], kendall[3], kendall[4]))
  #spearman.out <- rbind(spearman.out, c(row.names(comm_therm.mat)[i], spearman[1], spearman[2], spearman[3], spearman[4]))
  #kw.out <- rbind(kw.out, c(row.names(comm_therm.mat)[i], results.k[1], results.k[2], results.k[3])) 
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
Acomm_sigcor_pos <- Acomm_sigcor[as.numeric(sigcor[,3])>0,]
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
Phyla_Therm <- Phyla_Therm[-46]
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

Facet_Categories <- c(rep("A",length(Proportion_OTUs)), rep("B",length(Proportion_RA)), rep("C", nrow(JGI_Therm)))
Fill_Categories <- c(rep("GG", length(Proportion_OTUs)/2), rep("DN",length(Proportion_OTUs)/2), rep("GG", length(Proportion_OTUs)/2), rep("DN",length(Proportion_OTUs)/2), rep("Genomes", nrow(JGI_Therm)))

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

pr <- ggplot(data=PlotData_R, aes(x=Phylum_Column, y=Proportion, fill=Fill_Categories)) + geom_bar(stat="identity", aes(x=Phylum_Column, y=Proportion), colour="black") + facet_wrap(~Category_f, nrow=3, strip.position="right") + theme(strip.background = element_blank(), strip.text.y=element_blank(), axis.text.x = element_text(angle=60, hjust=1, vjust=1, lineheight = rel(2), size=10), axis.title.x= element_text(vjust=1),strip.text.x=element_text(size=rel(2)), axis.text.y=element_text(size=15))  + scale_fill_manual(values=cbPalette)
pr
ggsave("../Figures/Fig3B.eps", width=200, units="mm")

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