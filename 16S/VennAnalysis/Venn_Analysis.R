#Reading in the OTU Table
OTUs<-read.table("subsamplingtable_even73419_R.txt", header=TRUE, row.names=1, sep="\t")

#Merging/Summing the Replicates of the Samples
Names<-c("C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11","C12","C13","C14","C15","C16","C17","C18")
Merged<-data.frame(C01=rowSums(OTUs[,1:3]))
for (x in 2:18) {
  Merged[,Names[x]]<-rowSums(OTUs[,(3*x-2):(3*x)])
}
#Adding on the Taxonomy
Merged$Taxonomy<-OTUs[,55]
Temperature<-c(14.1,13.6,14.7,13.3,14.0,24.1,13.5,12.6,34.2,54.2,21.1,32.0,57.4,34.1,38.9,21.7,12.1,13.3,0)
FireFront<-c(2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,0)
RowNames<-row.names(Merged)
RowNames<-c(RowNames,"Temperature")
row.names(Merged_Metadata)<-RowNames
Guess<-Merged_Metadata["Temperature",]>50
Merged_Metadata["Temperature",]

#Getting Samples with Temperature >50C
Hot_Only<-Merged_Metadata[,which(as.numeric(Merged_Metadata["Temperature",])>50)]
Merged_Metadata$HotSoils<-rowSums(Merged_Metadata[,which(as.numeric(Merged_Metadata["Temperature",])>50)])
Merged_Metadata$ColdSoils<-rowSums(Merged_Metadata[,which(as.numeric(Merged_Metadata["Temperature",])<15)])
Hot_Only$Taxonomy<-Merged_Metadata[,19]

#Getting Samples with Temperature <15C
Cold_Only<-Merged_Metadata[,which(as.numeric(Merged_Metadata["Temperature",])<15)]
Cold_Only$Taxonomy<-Merged_Metadata[,19]


#Hot Soil Exclusive Taxa:
Thermophilic_Taxa_Exclusive<-Merged_Metadata[Merged_Metadata$HotSoils>0 & Merged_Metadata$ColdSoils<1,] 
#Write out themophilic exclusive Taxa
write.table(Thermophilic_Taxa_Exclusive, file="Thermophilic_Taxa_Exclusive.txt", col.names=TRUE, row.names=TRUE, sep="\t")


#Taxa that are 10X more prevalent in soils >50C than in soils <15C
No_Zeros<-Merged_Metadata[Merged_Metadata$HotSoils>0 & Merged_Metadata$ColdSoils>0,]
Thermophilic_Taxa_10X<-No_Zeros[No_Zeros$HotSoils>10*No_Zeros$ColdSoils,]
#Write out differentially abundant taxa
write.table(Thermophilic_Taxa_10X, file="Thermophilic_Taxa_10X.txt", sep="\t", col.names=TRUE, row.names=TRUE)