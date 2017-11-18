#process the biom and otu tables
libs <- c("phyloseq", "RDPutils", "plyr", "ggplot2", "wesanderson")
lapply(libs, require, character.only = TRUE)

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/CHapter 2 - FC interactions/AMplicon analysis/otu_tables")

biom <- import_biom("nam.16S.tax.biom")
tax <- read.delim("nam.16S_perfect_tax.txt")

otu_table <- merge_phyloseq(biom, tax)

otu_table(biom, taxa_are_rows)

