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

#This is the actual modelling step for sparcc
sparcc.est <- sparcc(tax)
sparcc.graph <- abs(sparcc.est$Cor) >= 0.6
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)

ig.sparcc <- adj2igraph(sparcc.graph, vertex.attr=list(name=taxa_names(physeq)))


#Visualize the plots
library(igraph)
## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(tax, 1))+6
am.coord <- layout.fruchterman.reingold(ig.sparcc)


#get clusters
clusters <- cluster_fast_greedy(ig.sparcc)

#using phyloseq
plot_network(ig.sparcc, physeq, type='taxa', color="Phylum", label=NULL)
#Export to cytoscape
write.graph(ig.sparcc,file="sparcc.ncol.txt",format="ncol") 
write.table(TAX,file="taxonomy_network.txt",sep="\t", quote=FALSE)
