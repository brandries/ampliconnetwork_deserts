#Script to get otu table into phyloseq

libs <- c("ggplot2", "seqtime", "SpiecEasi", "plyr",  "reshape", "Matrix", "phyloseq", "ape", "igraph", "BIOMASS")
lapply(libs, require, character.only = TRUE)

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 2 - FC interactions/FC data")


#read the otu table, mapping file and reference taxonomy sets
tax <- read.delim("./16S/rarefaction/FC16S_tarare13208_mrdna_reorder_aus_sur_t_avelarger1_top500.txt", row.names = 1)
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


#network
spiec.out=spiec.easi(physeq, method="mb", lambda.min.ratio = 1e-2, nlambda = 20,icov.select.params=list(rep.num=20))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(physeq)))
plot_network(spiec.graph, physeq, type='taxa', color="Phylum", label=NULL)

#get positive and negative interaction
betaMat=as.matrix(symBeta(getOptBeta(spiec.out)))


otu.ids=colnames(spiec.out$data)
edges=E(spiec.graph)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(spiec.graph,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,"forestgreen")
  }else if(beta<0){
    edge.colors=append(edge.colors,"red")
  }
}
E(spiec.graph)$color=edge.colors

#now plot using red and green for positive and negative
#this is not very effective and should not be used
spiec.graph.b=spiec.graph
nodenames=V(spiec.graph.b)$name
V(spiec.graph.b)$name=getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
E(spiec.graph.b)$arrow.size=5
V(spiec.graph.b)$color="white"
V(spiec.graph.b)$frame.color="black"
tkplot(spiec.graph.b)

#cluster the network
clusters=cluster_fast_greedy(spiec.graph)
clusterOneIndices=which(clusters$membership==1)
clusterOneOtus=clusters$names[clusterOneIndices]
clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]
clusterThreeIndices=which(clusters$membership==3)
clusterThreeOtus=clusters$names[clusterThreeIndices]
clusterFourIndices=which(clusters$membership==4)
clusterFourOtus=clusters$names[clusterFourIndices]
clusterFiveIndices=which(clusters$membership==5)
clusterFiveOtus=clusters$names[clusterFiveIndices]

sort(table(getTaxonomy(clusterOneOtus,taxa.f,useRownames = TRUE)),decreasing = TRUE)

#Export to cytoscape
write.graph(spiec.graph,file="spieceasi.ncol.txt",format="ncol") 
write.table(TAX,file="taxonomy_network.txt",sep="\t", quote=FALSE)
