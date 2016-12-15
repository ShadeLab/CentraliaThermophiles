
map <- read.table("Centralia_Collapsed_Map_forR.txt",sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)

map_MG <- map[c(1,3,4,5,6,7,10,12,14,15,16,17),]

KO <- read.table("KO_minus2col_09-27-2016.tab.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactors = FALSE)
colnames(KO) <- map_MG$Sample

rpoB <- KO["KO:K03043",] + KO["KO:K13798",] ### add the bacterial and archaeal copies together.
#Relativze to rpoB

KO.rpob <- NULL
for (i in 1:nrow(KO)){
  KO.rpob <- rbind(KO.rpob, KO[i,]/rpoB)
}

KO_NZ <- KO[rowSums(KO)>0,]
KO.rpob_NZ <- KO.rpob[rowSums(KO.rpob)>0,]
row.names(KO_NZ) <- sub("KO:", "", row.names(KO_NZ))
row.names(KO.rpob_NZ) <- sub("KO:", "", row.names(KO.rpob_NZ))

Modules <- read.table("kmodlist47982_23-nov-2016.txt", sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)

Dictionary <- vector(mode="list", length=nrow(Modules))
names(Dictionary) <- row.names(Modules)

Dictionary[[1]]
grep(row.names(KO.rpob_NZ)[10], Modules$Definition)


grep("K18304", Modules$Definition)


#Produces a list, each item in this list is a dataframe with the KO's in our dataset 
for (KO in 1:nrow(KO.rpob_NZ)){
  KO_M <- grep(row.names(KO.rpob_NZ)[KO], Modules$Definition)
  if (length(KO_M)>0){
    for (x in 1:length(KO_M)){
      Dictionary[[KO_M[x]]] <- rbind(Dictionary[[KO_M[x]]], KO.rpob_NZ[KO,])
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

KO_Per_M_Data/KO_Per_M

Dictionary_Subset <- Dictionary[lapply(Dictionary, data.class) != "NULL"]
Modules_Subset <- Modules[lapply(Dictionary,data.class) != "NULL",]
KO_Per_M_Subset <- KO_Per_M[lapply(Dictionary, data.class) != "NULL"]
KO_Per_M_Data_Subset <- KO_Per_M_Data[lapply(Dictionary, data.class) != "NULL"]

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
  #Summary_Output <- rbind(Summary_Output, c(nrow(K_out),apply(K_out, 2, mean), apply(K_out, 2, sd)))
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

#row.names(Summary_Output) <- names(Dictionary_Subset)  
row.names(avg_mod) <- names(Dictionary_Subset)
row.names(mid_mod) <- names(Dictionary_Subset)
row.names(max_mod) <- names(Dictionary_Subset)
row.names(min_mod) <- names(Dictionary_Subset)
row.names(stdev_mod) <- names(Dictionary_Subset)

### PERMANOVA based on the Module level aggregation Median values
library(vegan)
Mod.d <- vegdist(t(mid_mod),method="bray")
Class <- rep("Red", 12)
Class[map_MG$Classification!="FireAffected"] <- "Green"
a=adonis(Mod.d~Class, distance=TRUE, permutations=1000)
a

# Temperature Correlations for Each Module
median_Module_TempCorrelations <- apply(mid_mod, 1, function(x) cor.test(x, map_MG$SoilTemperature_to10cm))
# T-Test for each Module
median_Module_t.test <- apply(mid_mod, 1, function(x) t.test(x[map_MG$Classification=="FireAffected"], x[map_MG$Classification!="FireAffected"]))
Med_Mod_TempCor <- NULL
for (i in 1:600){
  Med_Mod_TempCor <- rbind (Med_Mod_TempCor, c(median_Module_TempCorrelations[[i]][4],median_Module_TempCorrelations[[i]][2],median_Module_TempCorrelations[[i]][3]))
}
Med_Mod_Ttest <- NULL
for (i in 1:600){
  Med_Mod_Ttest <- rbind (Med_Mod_Ttest, c(median_Module_t.test[[i]][1],median_Module_t.test[[i]][2],median_Module_t.test[[i]][3]))
}


Med_Mod_TempCor <- as.data.frame(Med_Mod_TempCor)
row.names(Med_Mod_TempCor) <- names(Dictionary_Subset)
Med_Mod_TempCor[,1] <- unlist(Med_Mod_TempCor[,1])
Med_Mod_TempCor[,2] <- unlist(Med_Mod_TempCor[,2])
Med_Mod_TempCor[,3] <- unlist(Med_Mod_TempCor[,3])

Med_Mod_Ttest <- as.data.frame(Med_Mod_Ttest)
row.names(Med_Mod_Ttest) <- names(Dictionary_Subset)
Med_Mod_Ttest[,1] <- unlist(Med_Mod_Ttest[,1])
Med_Mod_Ttest[,2] <- unlist(Med_Mod_Ttest[,2])
Med_Mod_Ttest[,3] <- unlist(Med_Mod_Ttest[,3])

Sig_Ttest_Modules <- Modules_Subset[Med_Mod_Ttest[,3]<0.05,]
Sig_TempCor_Modules <- Modules_Subset[Med_Mod_TempCor[,3]<0.05,]

Module_Completeness_Ttest <- Module_Completeness[Med_Mod_Ttest[,3]<0.05]
Module_Completeness_TempCor <- Module_Completeness[Med_Mod_TempCor[,3]<0.05]

Complete_Sig_Ttest_Modules <- Sig_Ttest_Modules[Module_Completeness_Ttest==1,]
Complete_TempCor_Modules <- Sig_TempCor_Modules[Module_Completeness_TempCor==1,]

x <- Complete_Sig_Ttest_Modules[order(Complete_Sig_Ttest_Modules$Module.Type),]

write.table(Complete_Sig_Ttest_Modules[order(Complete_Sig_Ttest_Modules$Module.Type),], "T_Test_Sig_Modules.txt", sep="\t", quote=FALSE)
write.table(Complete_TempCor_Modules[order(Complete_TempCor_Modules$Module.Type),], "TempCor_Sig_Modules.txt", sep="\t", quote=FALSE)

#Summary_Output <- cbind(Summary_Output,KO_Per_M_Subset)
#Sig_Cor_Modules <- Summary_Output[Summary_Output[,3]<0.05,]
#Pos_Sig_Cor_Modules <- Sig_Cor_Modules[Sig_Cor_Modules[,2]>0,]

Complete_Pos_Sig_Cor_Modules <- Pos_Sig_Cor_Modules[(Pos_Sig_Cor_Modules[,1]/Pos_Sig_Cor_Modules[,10])==1,]


sum(1*((Pos_Sig_Cor_Modules[,1]/Pos_Sig_Cor_Modules[,10])==1))
Modules[row.names(Complete_Pos_Sig_Cor_Modules),]

library(colorRamps)
library(gplots)
hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

mid_mod_Ttest <- mid_mod[Med_Mod_Ttest[,3]<0.05,]
mid_mod_Ttest.zs <- decostand(mid_mod_Ttest, method="standardize", MARGIN=1)
mid_mod_Ttest.zs.temp <- mid_mod_Ttest.zs[,order(map_MG$SoilTemperature_to10cm)]

mid_mod_TempCor <- mid_mod[Med_Mod_TempCor[,3]<0.05,]
mid_mod_TempCor.zs <- decostand(mid_mod_TempCor, method="standardize", MARGIN=1)
mid_mod_TempCor.zs.temp <- mid_mod_TempCor.zs[,order(map_MG$SoilTemperature_to10cm)]  

setEPS()

postscript("../Figures/Heatmap_Ttest.eps", width = 5, height=5, pointsize=8,paper="special")
fig10 <- heatmap.2(mid_mod_Ttest.zs, col=hc(100), key=TRUE, symkey=FALSE, trace="none", density.info="none", dendrogram="both", labRow=FALSE, margins=c(5,13), srtCol=90)

dev.off()

setEPS()
postscript("../Figures/Heatmap_Ttest_TempOrdered.eps", width = 5, height=5, pointsize=8,paper="special")
figure11 <- heatmap.2(mid_mod_Ttest.zs.temp, col=hc(100), key=TRUE, symkey=FALSE, trace="none", density.info="none", dendrogram="row", Colv= FALSE, labRow=FALSE, margins=c(5,13), srtCol=90, main ="T-test Significant Modules")
dev.off()

setEPS()
postscript("../Figures/Heatmap_TempCor.eps", width = 5, height=5, pointsize=8,paper="special")
Figure12 <- heatmap.2(mid_mod_TempCor.zs, col=hc(100), key=TRUE, symkey=FALSE, trace="none", density.info="none", dendrogram="both", labRow=FALSE, margins=c(5,13), srtCol=90)
dev.off()

setEPS()
postscript("../Figures/Heatmap_TempCor_TempOrdered.eps", width = 5, height=5, pointsize=8,paper="special")
Figure13 <-heatmap.2(mid_mod_TempCor.zs.temp, col=hc(100), key=TRUE, Colv=FALSE, symkey=FALSE, trace="none", density.info="none", dendrogram="row", labRow=FALSE, margins=c(5,13), srtCol=90, main ="Temperature Correlated Modules")
dev.off()


### Area for Trial Code etc. 
install.packages("pcvlust")
library(pvclust)

result <- pvclust(t(mid_mod.zs[Med_Mod_Cor[,3]<0.05,]), method.dist="cor", method.hclust="average", nboot=1000)

plot(result)
pvrect(result, alpha=0.6)



setEPS()

postscript("../Figures/Fig1.eps", width = 5, height=5, pointsize=8,paper="special")
par(mar=c(5,3,2,2)+0.1)
fig1=vennDiagram(v)
dev.off()




