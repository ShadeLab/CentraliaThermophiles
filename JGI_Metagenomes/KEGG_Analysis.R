setwd("/Users/JSorensen/Gradient_Metagenome/")### Set the working directory
KEGG<- read.table("Gradient_KEGG.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)### Read in KEGG Table from JGI
colnames(KEGG)<-c("C01","C03","C04","C05","C06","C07", "C10", "C12", "C14", "C15", "C16", "C17")### Replace the headers from JGI with Sample names.


Context<-read.table("/Users/JSorensen/Gradient_Metagenome/Centralia_Full_Map_Fixed.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)### Read in Mapping file

###Subset the Mapping file
Samples<-colnames(KEGG)
Context_Sub=NULL
for(i in 1:length(Samples)){
  new=Context[Context[,"Sample"]==Samples[i],]
  Context_Sub=rbind(Context_Sub,new)
}

###Collapse the mapping file
Numbers<-c(1,4,7,10,13,16,19,22,25,28,31,34)
Context_Collapsed<-Context_Sub[Numbers,]


FA_Sites<- Context_Collapsed[Context_Collapsed$Classification=="FireAffected",]

OnlyFireAffected<-KEGG[,FA_Sites$Sample]
OnlyFireAffted_1<- OnlyFireAffected[rowSums(OnlyFireAffected)>0,]
sum(1*(rowSums(OnlyFireAffected)>0))
write.table(row.names(OnlyFireAffted_1), "OnlyFireAffectedKOs.txt", sep="\t")

KEGG<-as.matrix(KEGG)
colSums(KEGG)

rpoB <- KEGG["KO:K03043",] + KEGG["KO:K13798",] ### add the bacterial and archaeal copies together.

### Make Data set relativized to the number of rpoB genes.
KEGG_rel_Rpob<- NULL
for(i in 1:nrow(KEGG)){
  row <- KEGG[i,]/rpoB
  KEGG_rel_Rpob<- rbind(KEGG_rel_Rpob, row)
}
row.names(KEGG_rel_Rpob)<- row.names(KEGG)
library(vegan)
KEGG_rel=decostand(KEGG, method="total", MARGIN=2)


### Updated 1/22/2016 to use the rpoB relativized data instead of the the % relativized
coretest.out=NULL
noncoretest.out=NULL
TempCorKEGG_rel=NULL
NonCorKEGG=NULL
NoSigcoretest.out=NULL
for(i in 1:nrow(KEGG_rel_Rpob)){
  results=cor.test(as.numeric(KEGG_rel_Rpob[i,]),Context_Collapsed$SoilTemperature_to10cm)
  if(results[4]>0){
    if(results[3]<0.05){
      TempCorKEGG_rel=rbind(TempCorKEGG_rel,c(KEGG_rel_Rpob[i,],row.names(KEGG_rel_Rpob)[i]))
      coretest.out=rbind(coretest.out,c(row.names(KEGG_rel_Rpob)[i],results$estimate,results$p.value))
    }
    if(results[3]>0.05){
      NoSigcoretest.out=rbind(NoSigcoretest.out,c(row.names(KEGG_rel_Rpob)[i],results$estimate,results$p.value))
    }    }
  if(results[4]<0){
    NonCorKEGG=rbind(NonCorKEGG,KEGG_rel_Rpob[i,])
    noncoretest.out=rbind(noncoretest.out,c(row.names(KEGG_rel_Rpob)[i],results$estimate,results$p.value))
  }

}


### Write the list of temperature correlated kegg orthologs
write.table(coretest.out[,1], "Temperature_Correlated_KOs_rpobrel.txt")
