#Script to get otu table into phyloseq

libs <- c("ggplot2", "ggthemes", "SpiecEasi", "plyr",  "reshape", "Matrix", "phyloseq", "ape")
lapply(libs, require, character.only = TRUE)

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 2 - FC interactions/FC data")


#read the otu table, mapping file and reference taxonomy sets
tax <- read.delim("./16S/rarefaction/FC16S_tarare13208_mrdna_reorder_aus_t_avelarger1_AFD_nosmaller20.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata.txt", row.names = 1)
taxonomy_mapping <- read.delim("./16S/rarefaction/networks/FC16S_rep_numbered_tax_assignments_phyloseq.txt", row.names = 1)

#transpose the otu table and make the reference taxa a matrix
phyloseq_otu <- t(tax)
taxonomy_mapping <- as.matrix(taxonomy_mapping)

#create otu table and tax for phyloseq objects and make the combined object
OTU <- otu_table(phyloseq_otu, taxa_are_rows = T)
TAX <- tax_table(taxonomy_mapping)
SAMPL <- sample_data(mapping_file)
physeq = phyloseq(OTU, TAX, SAMPL)
physeq


#some phyloseq stuff
#barplot
plot_bar(physeq, fill = "Phylum")

#diversity measures
plot_richness(physeq)

#tree
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq1 = merge_phyloseq(physeq, random_tree)
plot_tree(physeq1, color = "Loc_dep", label.tips = "taxa_names", size = "abundance")

#network
spiec.out=spiec.easi(physeq, method="mb",icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(physeq)))
plot_network(spiec.graph, physeq, type='taxa', color="Phylum", label=NULL)



