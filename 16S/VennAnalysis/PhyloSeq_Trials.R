library("phyloseq")
library("ggplot2")
library("scales")
library("grid")

mapfile="/Users/JSorensen/Desktop/VennAnalysis/Centralia_full_map_corrected.txt"
map=import_qiime_sample_data(mapfile)
class(map)

biomfile="/Users/JSorensen/Desktop/VennAnalysis/Subsampling_otu_table_even73419.biom"
Cen_biom=import_biom(biomfile, parseFunction=parse_taxonomy_greengenes)
class(Cen_biom)