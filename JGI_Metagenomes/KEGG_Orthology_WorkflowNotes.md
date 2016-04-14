# KEGG Orthology workflow
## Rationale/Approach
We are trying to discover novel taxonomic and functional diversity of thermophiles in Centralia, PA. By using the metagenomes across a gradient of temperatures the goal is to identify specific funtions that are more prevalent in the metagenomes from hot sites.

### Accessing KEGG data from IMG
Login to [IMG-JGI](img.jgi.doe.gov). Click on `Find Genomes` ---> `Genome Search`. Search for "Centralia" under the filer "Genome Name" and press `Go`. Press `Select all` and then `Add Selected to Genome Cart`. Then click `Compare Genomes` ---> `Abundance Profiles` ---> `Overview (All Functions)`. Click `Matrix`, `Estimated gene copies`, and `KO`. Select `Assembled` from MER-FS Metagenome, `All Finished, Permanent Draft and Draft` from Sequencing Status and `Genome Cart` from Domain. Click `show` and select all 12 metagenomes and press `Go`. Wait for IMG to do its thing and once the page is loaded with results click "Download tab-delimited file for Excel". Silly IMG this is clearly for R.

### Temperature Correlation of KO's in R
Look at the R script file `KEGG_Analysis.R`.
```R
setwd("/Users/JSorensen/Gradient_Metagenome/")
KEGG<- read.table("Gradient_KEGG.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
colnames(KEGG)<-c("C01","C03","C04","C05","C06","C07", "C10", "C12", "C14", "C15", "C16", "C17")


Context<-read.table("/Users/JSorensen/Gradient_Metagenome/Centralia_Full_Map_Fixed.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

Samples<-colnames(KEGG)
Context_Sub=NULL
for(i in 1:length(Samples)){
  new=Context[Context[,"Sample"]==Samples[i],]
  Context_Sub=rbind(Context_Sub,new)
}

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

write(coretest.out[,1], "Temperature_Correlated_KOs_rpobrel.txt")
```

Need to alter the file so that it works in the KEGG pathway tool. Open the output file from the Rscript in a plain text editior and find and replace `KO:` with nothing. Also find and replace all `"` with nothing. Finally, replace the first line with `#IdentifierforyourKEGGS`. The example file [Temperature_Correlated_KOs_rpobrel_ForPathway.txt](Temperature_Correlated_KOs_rpobrel_ForPathway.txt) is included in the directory.

### Using the KEGG Pathway Reconstruct tool
Go to the [KEGG Mapper Tool](http://www.genome.jp/kegg/tool/map_pathway.html). Upload the formatted output file from the R script(containing one column of numbers and a second column of KO numbers), check `Include global/overview maps` and press `exec`. If you want to look at multiple sets of KOs at a time, like positively correlated temperature KOs and negatively correlated temperature KOs, you need them both in the same file but the lists separated by a header like as is the first line(IE header for group 1, all KOs for group 1, header for group 2, all KOs for group 2).   
