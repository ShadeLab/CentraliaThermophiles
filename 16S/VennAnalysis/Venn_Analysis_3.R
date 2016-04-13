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

### Identifying Cold and Hot Sites
Cold_Soils<-unique(map$Sample[map$SoilTemperature_to10cm<20])
Hot_Soils<-unique(map$Sample[map$SoilTemperature_to10cm>50])

### Identifying Hot Site Exclusive OTUs
Hot_Exclusive_OTUs=NULL
index=NULL
for (i in 1:nrow(collapsed_soils_rel)){
  if (sum(collapsed_soils_rel[i,Hot_Soils])>0 & sum(collapsed_soils_rel[i,Cold_Soils])==0){
    index<-c(index,i)
    Hot_Exclusive_OTUs<-rbind(Hot_Exclusive_OTUs, collapsed_soils_rel[i,])
  }
}

rownames(Hot_Exclusive_OTUs)<- rownames(collapsed_soils_rel[index,])
colnames(Hot_Exclusive_OTUs)<-colnames(collapsed_soils_rel)

soils["806726",]
rdp[index[1]]

###Parsing Taxonomy
taxonomy.df=NULL
for (i in 1:nrow(Hot_Exclusive_OTUs)){
  taxonomy.df<-rbind(taxonomy.df, as.vector(unlist(strsplit(Hot_Exclusive_OTUs[i,19],";")),mode="list"))
}

Family<-grep("f__", taxonomy.df[,5], value=TRUE)
Order<-grep("o__", taxonomy.df[,4], value=TRUE)
Class<-grep("c__", taxonomy.df[,3], value=TRUE)
Phylum<-grep("p__", taxonomy.df[,2], value=TRUE)



taxonomy_FULL.df=NULL
for (i in 1:nrow(soils)){
  taxonomy_FULL.df<-rbind(taxonomy_FULL.df, as.vector(unlist(strsplit(soils[i,55],";")),mode="list"))
}


Phylum_Full<-grep("p__", taxonomy_FULL.df[,2], value=TRUE)
unique(Phylum_Full)



