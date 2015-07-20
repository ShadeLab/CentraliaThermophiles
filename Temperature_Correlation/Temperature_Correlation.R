### Reading in the Data
soils=read.table("subsamplingtable_even73419_R.txt", header=TRUE, sep="\t", row.names=1, stringsAsFactors=FALSE)#read in the otu table
rdp=soils[,ncol(soils)]# save the taxonomy as a vector
map=read.table("Centralia_Full_Map_Fixed.txt", header=TRUE, sep="\t", row.names=1, check.names=FALSE, stringsAsFactors=FALSE)#read in the mapping file
soils=soils[,-ncol(soils)] #remove consensus lineage column
library(reshape2)
tmap=t(map) #Transpose the map so that sample names are column names and metadata labels are row names 
temp=tmap["SoilTemperature_to10cm",]
Sample_Temp=as.numeric(temp[c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52)])#Creat Vector of Sample Temperatures


### Collapsing the OTU table by summing across the replicates
u=unique(map[,"Sample"]) #Grab the unique sample names from the column named "Sample" in the map file
collapsed_soils=NULL
for(i in 1:length(u)){
  new=soils[,map[,"Sample"]==u[i]]
  sNew=rowSums(new)
  collapsed_soils=cbind(collapsed_soils,sNew)
}
colnames(collapsed_soils)=u
### Relativizing the collapsed data set
library(vegan)
collapsed_soils_rel=decostand(collapsed_soils, method="total", MARGIN=2)

head(collapsed_soils_rel)
length(Sample_Temp)
###Calculate Pearson OTU vs Temperature (+) Correlation p-value<0.05
coretest.out=NULL
TempCorOTUs_rel=NULL
TempCorOTUs_ABS=NULL
for(i in 1:nrow(collapsed_soils_rel)){
  results=cor.test(collapsed_soils_rel[i,],Sample_Temp)
  if(results[4]>0){
    if(results[3]<0.05){
      TempCorOTUs_rel=rbind(TempCorOTUs_rel,c(collapsed_soils_rel[i,],rdp[i]))
      TempCorOTUs_ABS=rbind(TempCorOTUs_ABS,c(collapsed_soils[i,],rdp[i]))
      coretest.out=rbind(coretest.out,c(row.names(collapsed_soils_rel)[i],results$estimate,results$p.value))
                        }
  }
  
}
colnames(TempCorOTUs_rel)=c(colnames(collapsed_soils_rel),"Taxonomy")
colnames(coretest.out)=c("OTU_ID","Correlation_Coefficient","p-value")
colnames(TempCorOTUs_ABS)=c(colnames(collapsed_soils_rel),"Taxonomy")
rownames(TempCorOTUs_rel)=coretest.out[,1]
rownames(TempCorOTUs_ABS)=coretest.out[,1]
range(coretest.out[,2])
range(coretest.out[,3])
summary(coretest.out[,2])
sum(1*(coretest.out[,2]>0.75))

OTUs_R75=coretest.out[coretest.out[,2]>0.75,]

### Make and write out a summed table of otus for temperature correlated otus
summed=NULL
summed_row=NULL
for (i in 1:nrow(TempCorOTUs_ABS)){
  summed_row=sum(as.double(TempCorOTUs_ABS[i,1:18]))
  summed=rbind(summed,c(summed_row,TempCorOTUs_ABS[i,19]))
}
write.table(summed, file="/Users/JSorensen/Desktop/VennAnalysis/Summed_TempCor_OTUs.txt", sep="\t")

### Calculate % of community per sample that is from otus that respond positively to temperature
Percent_Thermophile=c(rep(0,18))
for(i in 1:18){Percent_Thermophile[i]=sum(as.numeric(TempCorOTUs[,i]))}
Percent_Thermophile
plot(Sample_Temp,Percent_Thermophile)

Positivie_sig_PA=1*(TempCorOTUs_ABS>0)

Positive_sig_otus=colSums(Positivie_sig_PA)
Positive_sig_otus_Percent=NULL
Positive_sig_otus_Percent=Positive_sig_otus/dim(collapsed_soils_rel)[1]
Positive_sig_otus_Percent
plot(Sample_Temp,Positive_sig_otus_Percent[1:18])

###Calculate Pearson OTU vs Temperature Correlation p-value>0.05
coretest.out_nonsig=NULL
TempCorOTUs_nonsig=NULL
for(i in 1:nrow(collapsed_soils_rel)){
  results=cor.test(collapsed_soils_rel[i,],Sample_Temp)
  if(results[4]>0){
    if(results[3]>0.05){
      TempCorOTUs_nonsig=rbind(TempCorOTUs_nonsig,c(collapsed_soils_rel[i,],rdp[i]))
      coretest.out_nonsig=rbind(coretest.out_nonsig,c(row.names(collapsed_soils_rel)[i],results$estimate,results$p.value))
    }
  }
  
}
colnames(TempCorOTUs_nonsig)=c(colnames(collapsed_soils_rel),"Taxonomy")
colnames(coretest.out_nonsig)=c("OTU_ID","Correlation_Coefficient","p-value")
rownames(TempCorOTUs_nonsig)=coretest.out_nonsig[,1]

### Calculate % of community per sample that is from otus that respond positively to temperature pvalue>0.05
Percent_Thermophile_nonsig=c(rep(0,18))
for(i in 1:18){Percent_Thermophile_nonsig[i]=sum(as.numeric(TempCorOTUs_nonsig[,i]))}
Percent_Thermophile_nonsig
plot(Sample_Temp,Percent_Thermophile_nonsig)

### Correlation for Negative Temperature Response p-value<0.05
coretest.out_neg=NULL
TempCorOTUs_neg=NULL
for(i in 1:nrow(collapsed_soils_rel)){
  results=cor.test(collapsed_soils_rel[i,],Sample_Temp)
  if(results[4]<0){
    if(results[3]<0.05){
      TempCorOTUs_neg=rbind(TempCorOTUs_neg,c(collapsed_soils_rel[i,],rdp[i]))
      coretest.out_neg=rbind(coretest.out_neg,c(row.names(collapsed_soils_rel)[i],results$estimate,results$p.value))
    }
  }
  
}
colnames(TempCorOTUs_neg)=c(colnames(collapsed_soils_rel),"Taxonomy")
colnames(coretest.out_neg)=c("OTU_ID","Correlation_Coefficient","p-value")
rownames(TempCorOTUs_neg)=coretest.out_neg[,1]

### Calculate % of community per sample that is from otus that respond negatively to temperature pvalue <0.05
Percent_NonThermophile=c(rep(0,18))
for(i in 1:18){Percent_NonThermophile[i]=sum(as.numeric(TempCorOTUs_neg[,i]))}
Percent_NonThermophile
plot(Sample_Temp,Percent_NonThermophile)

### Correlation for Negative Temperature Response p-value>0.05
coretest.out_neg_nonsig=NULL
TempCorOTUs_neg_nonsig=NULL
for(i in 1:nrow(collapsed_soils_rel)){
  results=cor.test(collapsed_soils_rel[i,],Sample_Temp)
  if(results[4]<0){
    if(results[3]>0.05){
      TempCorOTUs_neg_nonsig=rbind(TempCorOTUs_neg_nonsig,c(collapsed_soils_rel[i,],rdp[i]))
      coretest.out_neg_nonsig=rbind(coretest.out_neg_nonsig,c(row.names(collapsed_soils_rel)[i],results$estimate,results$p.value))
    }
  }
  
}
colnames(TempCorOTUs_neg_nonsig)=c(colnames(collapsed_soils_rel),"Taxonomy")
colnames(coretest.out_neg_nonsig)=c("OTU_ID","Correlation_Coefficient","p-value")
rownames(TempCorOTUs_neg_nonsig)=coretest.out_neg_nonsig[,1]

### Calculate % of community per sample that is from otus that respond negatively to temperature pvalue > 0.05
Percent_NonThermophile_nonsig=c(rep(0,18))
for(i in 1:18){Percent_NonThermophile_nonsig[i]=sum(as.numeric(TempCorOTUs_neg_nonsig[,i]))}
Percent_NonThermophile_nonsig
plot(Sample_Temp,Percent_NonThermophile_nonsig)


Percent_NonThermophile_nonsig + Percent_NonThermophile + Percent_Thermophile_nonsig + Percent_Thermophile

### Correlation for Negative Temperature Response p-value>0.05
No_coretest.out=NULL
NoCorOTUs=NULL
for(i in 1:nrow(collapsed_soils_rel)){
  results=cor.test(collapsed_soils_rel[i,],Sample_Temp)
  results[4]
  if(results[4]==0){
      NoCorOTUs=rbind(NoCorOTUs,c(collapsed_soils_rel[i,],rdp[i]))
      No_coretest.out=rbind(No_coretest.out,c(row.names(collapsed_soils_rel)[i],results$estimate,results$p.value))
    }
  }
  
}
colnames(NoCorOTUs)=c(colnames(collapsed_soils_rel),"Taxonomy")
colnames(No_coretest.out)=c("OTU_ID","Correlation_Coefficient","p-value")
rownames(NoCorOTUs)=No_coretest.out[,1]






row.names(combined)[1]

u=unique(map[,"Sample"])

u
out=NULL
for(i in 1:length(u)){
  new=soils[,map[,"Sample"]==u[i]]
  sNew=rowSums(new)
  out=cbind(out,sNew)
}
colnames(out)=u








library(vegan)
out.REL=decostand(out, method="total", MARGIN=2)

head(out.REL)

colSums(out.REL)
new=data(,data(u[i,]))


### Begin Experimentation
plot(combined[1,],combined["SoilTemperature_to10cm",])
combined=rbind(soils,tmap) # add metadata to the end of the otu table
combined[1,]
combined["SoilTemperature_to10cm",]
result=cor.test(as.numeric(combined[1,]),as.numeric(combined["SoilTemperature_to10cm",]))
result[3]
result[4]
### End Experimentation

### Try using relative-ized data
output=c(rep(0,nrow(soils)))### Maybe look at a Kendall since it is rank based, especially if pearson's gives too many