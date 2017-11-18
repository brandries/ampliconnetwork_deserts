libs <- c("ggplot2", "ggthemes", "SpiecEasi", "plyr",  "reshape", "Matrix", "igraph")
lapply(libs, require, character.only = TRUE)

#Import data
#Your data has to be a table with first column being the one of axis, and the second being the other axis.
#the third column contains the values used to draw the bar chart. 
#alternatively, you can input a OTU table style table, and change the format using the reshape::melt function.

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 2 - FC interactions/FC data")



tax <- read.delim("./16S/rarefaction/FC16S_tarare13208_mrdna_reorder_aus_t_avelarger1_AFD_nosmaller20.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata.txt")

#THIS IS ACTUALLY PERFORMED IN THE ALGORITH, AND IS ONLY DOCUMENTED HERE IF YOU WANT TO USE IT FOR SOMETHING ELSE
#This is a replacement for rarefaction, getting the sums for each sample, 
# normalizing by them, and multiplying by the smallest sample

depths <- rowSums(tax)
tax.filt.n <- t(apply(tax, 1, norm_to_total))
tax.filt.cs <- round(tax.filt.n * min(depths))
d <- ncol(tax)
n <- nrow(tax)
e <- d

set.seed(10010)
graph <- SpiecEasi::make_graph("cluster", d, e)
Prec <- graph2prec(graph)
Cor <- cov2cor(prec2cov(Prec))

synth <- synth_comm_from_counts(tax.filt.cs, mar = 2, distr = 'zinegbin', Sigma = Cor, n = n)
se.synth.est <- spiec.easi(synth, method='mb', lambda.min.ratio=0.5, nlambda=15)
#check your model
huge::huge.roc(se.synth.est$path, graph, verbose = T)


#This is the actual modelling step, do it for both spieceasi methods  
#You need to check what each parameter means
#IT SEEMS LIKE THE LAMBDA IS NB TO CONSIDER
se.mb.est <- spiec.easi(tax, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20, icov.select.params=list(rep.num=20))
se.gl.est <- spiec.easi(tax, method = 'glasso', lambda.min.ratio = 1e-2, nlambda = 20, icov.select.params=list(rep.num=20))

#This is the actual modelling step for sparcc
sparcc.est <- sparcc(tax)
sparcc.graph <- abs(sparcc.est$Cor) >= 0.8
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)

## Create igraph objects
ig.mb <- adj2igraph(se.mb.est$refit)
ig.gl <- adj2igraph(se.gl.est$refit)
ig.sparcc <- adj2igraph(sparcc.graph, vertex.attr=list(name=taxa_names(physeq)))
ig.synth <- adj2igraph(se.synth.est$refit)

#Visualize the plots
library(igraph)
## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(tax, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

#other layouts
am.coord <- layout_nicely(ig.mb)
am.coord <- layout.auto(ig.sparcc)

#get clusters
clusters <- cluster_fast_greedy(ig.sparcc)


par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")
plot(ig.synth, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="Synth.mb")

#using phyloseq
plot_network(ig.sparcc, physeq, type='taxa', color="Phylum", label=NULL)
#Export to cytoscape
write.graph(ig.sparcc,file="sparcc.ncol.txt",format="ncol") 
write.table(TAX,file="taxonomy_network.txt",sep="\t", quote=FALSE)


#GET THE GRAPHS OUT
write.graph(ig.sparcc, file = "edgelist_AFD_mb.gml", format = "gml")

#Evaluate the edge weights
library(Matrix)
elist.gl <- summary(triu(cov2cor(se.gl.est$opt.cov)*se.gl.est$refit, k=1))
elist.mb <- summary(symBeta(getOptBeta(se.mb.est), mode='maxabs'))
elist.sparcc <- summary(sparcc.graph*sparcc.est$Cor)

hist(elist.sparcc[,3], main="", xlab="edge weights")
hist(elist.mb[,3], add=TRUE, col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')

#the degree distributions

dd.gl <- degree.distribution(ig.gl)
dd.mb <- degree.distribution(ig.mb)
dd.sparcc <- degree.distribution(ig.sparcc)

plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"), 
       col=c("forestgreen", "red", "black"), pch=1, lty=1)