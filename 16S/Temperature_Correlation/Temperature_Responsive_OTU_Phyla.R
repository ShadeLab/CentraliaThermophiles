### Phyla Counts of all Heat Responsive OTUs
OTUs_R75=coretest.out[coretest.out[,2]>0.75,]
Temp_Responsive_OTUs=soils[coretest.out[,1],]
hmap<-as.matrix(Soils_75)

Temp_Responsive_taxonomy.df=NULL
for (i in 1:nrow(Temp_Responsive_OTUs)){
  Temp_Responsive_taxonomy.df<-rbind(Temp_Responsive_taxonomy.df, as.vector(unlist(strsplit(Temp_Responsive_OTUs[i,55],";")),mode="list"))
}


Temp_Responsive_Phyla=unique(Temp_Responsive_taxonomy.df[,2])

### Calculating Phyla counts of Heat Responsive OTUs
Temp_Responsive_OTUs_Phyla=NULL
for (i in 1:length(Temp_Responsive_Phyla)){
  Temp_Responsive_OTUs_Phyla=rbind(Temp_Responsive_OTUs_Phyla, c(Temp_Responsive_Phyla[i], sum(1*(Temp_Responsive_taxonomy.df[,2]==unlist(Temp_Responsive_Phyla[i])))))
}

Temp_Responsive_OTUs_Phyla=as.data.frame(Temp_Responsive_OTUs_Phyla)
### Bar chart Makings
sorted_Temp_Resonsive_OTUs_Phyla<- Temp_Responsive_OTUs_Phyla[order(-as.numeric(Temp_Responsive_OTUs_Phyla$V2)),]
adjustment=c(0,1)

x_ordered<- barplot(as.numeric(sorted_Temp_Resonsive_OTUs_Phyla$V2), main="Temperature Responsive OTUs", ylab="Number of OTUs", xlab="Phyla")
labs_ordered<-paste(Trimmed_Phyla_Names)
text(cex=.75, x=x_ordered+.25, y=-1.25, adj=c(0,.5), labels=labs_ordered, xpd=TRUE, srt=60, pos=2)

### Trials
sum(1*(Temp_Responsive_taxonomy.df[,2]==unlist(Temp_Responsive_Phyla[1])))

Trimmed_Phyla_Names=NULL
for (i in 1:length(sorted_Temp_Resonsive_OTUs_Phyla$V1)){
  Trimmed_Phyla_Names<- c(Trimmed_Phyla_Names, unlist(strsplit(unlist(sorted_Temp_Resonsive_OTUs_Phyla$V1[i]), split='__', fixed=TRUE))[2])
}

labs_ordered<-paste(sorted_Temp_Resonsive_OTUs_Phyla$V1)