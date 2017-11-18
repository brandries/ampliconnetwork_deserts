library(vegan)
library(ggplot2)

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 2 - FC interactions/FC data/")

otu_table_ITS <- read.delim("./ITS/rarefaction/summarize_taxa/FCITS_rare6767_L6_nomargin.txt", row.names = 1)
mapping_file <- read.delim("./ITS/metadata_nomargin.txt")
  
  
otu_transf <- decostand(otu_table_ITS, "hellinger")

#bray curtis in vegan
vegdist(t(otu_transf), "bray") -> d

#metaMDS
fit <- metaMDS(d, "bray", k = 2, trymax = 999)

#take a look at the results
fit

#extract the scores
data.scores <- as.data.frame(scores(fit))
data.scores



#add hulls to your plot
#grp.dco <- data.scores[data.scores$Site == "Dune Control", ][chull(data.scores[data.scores$Site == "Dune Control", c("NMDS1", "NMDS2")]),]
#grp.dce <- data.scores[data.scores$Site == "Dune Centre", ][chull(data.scores[data.scores$Site == "Dune Centre", c("NMDS1", "NMDS2")]),]
#grp.gpco <- data.scores[data.scores$Site == "Gravel Plain Control", ][chull(data.scores[data.scores$Site == "Gravel Plain Control", c("NMDS1", "NMDS2")]),]
#grp.gpce <- data.scores[data.scores$Site == "Gravel Plain Centre", ][chull(data.scores[data.scores$Site == "Gravel Plain Centre", c("NMDS1", "NMDS2")]),]

#hull.data <- rbind(grp.dco, grp.dce, grp.gpco, grp.gpce)
#hull.data


shps = c(22, 21, 24)
attach(mapping_file)
#plot using ggplots
plot_1 <- ggplot() +
  #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Site, group = Site), alpha = 0.30) +
  #scale_fill_brewer(type = div, palette = 'Set1')+
  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, fill = Loc_dep, shape = Source), size = 6) +  #add your points
  scale_fill_brewer(type = div, palette = 'Paired')+
  scale_shape_manual(values = shps) +
  coord_fixed(ratio = 1.05) + 
  theme_bw()  +
#  guides(fill = FALSE) +
#  guides(colour = FALSE) +
#  guides(shape = FALSE) +
#  guides(fill = guide_legend(title = NULL))+  #remove elements of legend
#  guides(colour = guide_legend(title = NULL))+
#  guides(shape = guide_legend(title = NULL))+
  annotate("text", x = 0.23, y = -0.41, label = "stress = 0.16", size = 8) +
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

plot_1
#save your plot
svg(filename = "ITS_6767_nmds_nomargin.svg", width = 12, height = 6, pointsize = 12)
plot_1
dev.off()

#use multiplot
multiplot(plot_1, plot_2, plot_3, cols = 2)
