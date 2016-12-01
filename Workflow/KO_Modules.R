
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
KO_Per_M_Subset <- KO_Per_M[lapply(Dictionary, data.class) != "NULL"]

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
stdev_mod <- NULL
for (m in 1:length(Dictionary_Subset)){
  temp_data <- Dictionary_Subset[[m]]
  K_out <- NULL
  agg <- NULL
  for(k in 1:nrow(temp_data)){
    cor <- cor.test(as.numeric(temp_data[k,]), map_MG$SoilTemperature_to10cm) 
    test <- t.test(as.numeric(temp_data[k,map_MG$Class_Color=="red"], as.numeric(temp_data[k, map_MG$Class_Color!="red"])))
    K_out <- rbind(K_out, c(as.numeric(cor[1]), as.numeric(cor[3]), as.numeric(test[1]), as.numeric(test[3])))
  }
  for (s in 1:ncol(temp_data)){
    avg <- mean(temp_data[,s])
    mid <- median(temp_data[,s])
    stdev <- sd(temp_data[,s])
    agg <- cbind(agg, c(avg, mid, stdev))
  }
  Summary_Output <- rbind(Summary_Output, c(nrow(K_out),apply(K_out, 2, mean), apply(K_out, 2, sd)))
  avg_mod <- rbind(avg_mod, agg[1,])
  mid_mod <- rbind(mid_mod, agg[2,])
  stdev_mod <- rbind(stdev_mod, agg[3,])
  }

colnames(Summary_Output) <- c("Number KOs", "Average Pearsons", "Average Pearsons Pvalue", "Average T Statistic", "Average T P value", "SD Pearsons", "SD Pearsons Pvalue", "SD T Statistic", "SD T Pvalue") 
colnames(avg_mod) <- map_MG$Sample
colnames(mid_mod) <- map_MG$Sample
colnames(stdev_mod) <- map_MG$Sample

row.names(Summary_Output) <- names(Dictionary_Subset)  
row.names(avg_mod) <- names(Dictionary_Subset)
row.names(mid_mod) <- names(Dictionary_Subset)
row.names(stdev_mod) <- names(Dictionary_Subset)

median_Module_Correlations <- apply(mid_mod, 1, function(x) cor.test(x, map_MG$SoilTemperature_to10cm))
Med_Mod_Cor <- NULL
for (i in 1:600){
  Med_Mod_Cor <- rbind (Med_Mod_Cor, c(median_Module_Correlations[[i]][1],median_Module_Correlations[[i]][2],median_Module_Correlations[[i]][3]))
}
Med_Mod_Cor <- as.data.frame(Med_Mod_Cor)
row.names(Med_Mod_Cor) <- names(Dictionary_Subset)
Pos_Med_Mod_Cor <- Med_Mod_Cor[Med_Mod_Cor[,1]>0,]
Sig_Pos_Med_Mod_Cor <- Pos_Med_Mod_Cor[Pos_Med_Mod_Cor[,3]<0.05,]

Summary_Output <- cbind(Summary_Output,KO_Per_M_Subset)
Sig_Cor_Modules <- Summary_Output[Summary_Output[,3]<0.05,]
Pos_Sig_Cor_Modules <- Sig_Cor_Modules[Sig_Cor_Modules[,2]>0,]

Complete_Pos_Sig_Cor_Modules <- Pos_Sig_Cor_Modules[(Pos_Sig_Cor_Modules[,1]/Pos_Sig_Cor_Modules[,10])==1,]


sum(1*((Pos_Sig_Cor_Modules[,1]/Pos_Sig_Cor_Modules[,10])==1))
Modules[row.names(Complete_Pos_Sig_Cor_Modules),]

library(colorRamps)
library(gplots)
mid_mod.zs <- decostand(mid_mod, method="standardize", MARGIN=1)
mid_mod.zs.temp <- mid_mod.zs[,order(map_MG$SoilTemperature_to10cm)]


hc=colorRampPalette(c("#91bfdb","white","#fc8d59"), interpolate="linear")

heatmap.2(mid_mod.zs,col=hc(100), key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", labRow=row.names(mid_mod.zs), margins=c(5,13), srtCol=90)

heatmap.2(mid_mod.zs[Med_Mod_Cor[,3]<0.05,],col=hc(100), key=TRUE,symkey=FALSE, trace="none", density.info="none",dendrogram="both", labRow=row.names(mid_mod.zs[Med_Mod_Cor[,3]<0.05,]), margins=c(5,13), srtCol=90)

heatmap.2(mid_mod.zs.temp[Med_Mod_Cor[,3]<0.05,],col=hc(100),scale="row",key=TRUE,symkey=FALSE, Colv = FALSE, main="Temp Correlated Modules", trace="none", density.info="none",dendrogram="row", labRow=row.names(mid_mod.zs.temp[Med_Mod_Cor[,3]<0.05,]), margins=c(5,13), srtCol=90)
