library(vegan)

#import data
setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 3 - FC interactions/FC data/")

#THIS SCRIPT IS TOO DIFFICULT TO REWRITE, SO JUST CHANGE THE TABLE FOR 16S INPUT
otu_table_ITS  <- read.delim("./ITS/rarefaction/FCITS_rare6767_notax_aus.txt", row.names = 1)
mapping_file <- read.delim("./ITS/metadata_nonam.txt")
#the infamous 16s table
otu_table_ITS <- read.delim("./16S/rarefaction/summarize_taxa1405/FC16S_tarare1405_L2_top20.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata_all.txt")

#Transpose the table so you can do adonis and do sqrt transformation
new_table <- as.data.frame(t(otu_table_ITS))

new_table <- vegdist(new_table)
#marge the table and mapping file
#combi_table <- cbind(mapping_file,new_table)
attach(mapping_file)

#test for depth and location and source
depth <- adonis(t(otu_table_ITS)~Depth)
source <- adonis(t(otu_table_ITS)~Source)
location <- adonis(t(otu_table_ITS)~Location)

#or as shown below, make a much more robust permanova and include all the vars
adonis(new_table~Depth*Source)


# you may have to subset this as your have multiple confounding factors
# look at how low the r^2 values are
# this is an indication that you are not explaining much of the variation
# you need to hcange your dataset to only look at the features you want to look at

#first you order by location
new_combi_frame <- combi_table[order(combi_table$Location),]

#then select features you want
australia_combi <- new_combi_frame[1:54,]
australia_biol <- australia_combi[,length(mapping_file):length(australia_combi)]
nam_combi <- new_combi_frame[55:94, ]
nam_biol <- nam_combi[,length(mapping_file):length(nam_combi)]
dune_combi <- nam_combi[1:20,]
dune_biol <- dune_combi[,length(mapping_file):length(dune_combi)]
gp_combi <- nam_combi[21:40,]
gp_biol <- gp_combi[,length(mapping_file):length(gp_combi)]

#now for some permanovas
#these are all single parameter permanovas and capture a small portion of the variance
adonis(australia_biol~australia_combi$Source)
adonis(nam_biol~nam_combi$Depth)
adonis(dune_biol~dune_combi$Depth)

#this is a much more realistic permanova, showing only source to be significant
adonis(gp_biol~gp_combi$Depth*gp_combi$Source)


betadispp <- betadisper(new_table, Source)
vecs <- betadispp$vectors
vect12 <- vecs[,1]
ANOVA <- aov(vect12~Source)
