#THIS IS AN ATTEMPT AT AN TSNE SCRIPT

libs <- c("ggplot2", "vegan", "tsne", "RColorBrewer")
lapply(libs, require, character.only = TRUE)

#import data
setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 2 - FC interactions/FC data/")
compare_genes_tax  <- read.delim("./16S/rarefaction/summarize_taxa13208/FC16S_tarare13208_L6_tsne.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata.txt")

otu_transf <- decostand(otu_table_ITS, "hellinger")
compare_genes_tax <- cbind(mapping_file, t(otu_transf))
#TRANSPOSE YOUR MARIX

biol.data <- compare_genes_tax[,10:297]

#compare_genes_tax <- compare_genes_tax[1:94,]
dist_mat <- vegdist(biol.data)



#this is a function you set, one for number of colors
#the other for the names
#then you plot the tsne and you get it to giev you a plot at each iteration
#colors = rainbow(length(unique(compare_genes_tax$Loc_sour)))
colors = c("#1b4f72", "#aed6f1", "#1d8348","#58d68d", "#943126",  "#ec7063")
colors = brewer.pal(length(unique(compare_genes_tax$Loc_dep)), "Paired")
names(colors) = unique(compare_genes_tax$Loc_dep)
ecb = function(x,y){ plot(x,t='p', col=colors[compare_genes_tax$Loc_dep], lwd=10) }
tsne_iris = tsne(dist_mat, epoch_callback = ecb, perplexity=20, max_iter = 3500, epoch = 50)


#PLOT
tsne_plot <- as.data.frame(tsne_iris)

attach(mapping_file)
plot_1 <- ggplot() +
  #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Site, group = Site), alpha = 0.30) +
  #scale_fill_brewer(type = div, palette = 'Set1')+
    geom_point(data = tsne_plot, aes(x = V1, y = V2, fill = Loc_dep, shape = Source), size = 6) +  #add your points
  scale_fill_brewer(type = "div", palette = "Paired")+
  scale_shape_manual(values = shps) +
  coord_fixed(ratio = 1.05) + 
  theme_bw()  +
  #  guides(fill = FALSE) +
  #  guides(colour = FALSE) +
  #  guides(shape = FALSE) +
  #  guides(fill = guide_legend(title = NULL))+  #remove elements of legend
  #  guides(colour = guide_legend(title = NULL))+
  #  guides(shape = guide_legend(title = NULL))+
  theme(axis.text.x = element_blank(), #most of arguments to follow is to remove features of the graph
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 18), #increase the size of x label text
        axis.title.y = element_text(size = 18), #increase size of y label text
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18))


#save your plot
svg(filename = "tsne_.svg", width = 12, height = 6, pointsize = 12)
plot_1
dev.off()
