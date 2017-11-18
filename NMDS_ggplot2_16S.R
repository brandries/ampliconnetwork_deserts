library(vegan)
library(ggplot2)

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 3 - FC interactions/FC data/")


#RAREFIED TO 1405 

otu_table_16S <- read.delim("./16S/rarefaction/summarize_taxa1405/FC16S_tarare1405_L6_nonam.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata_all_nonam.txt")
  
  
otu_transf <- decostand(otu_table_16S, "hellinger")

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
grp.dco <- data.scores[data.scores$Site == "Dune Control", ][chull(data.scores[data.scores$Site == "Dune Control", c("NMDS1", "NMDS2")]),]
grp.dce <- data.scores[data.scores$Site == "Dune Centre", ][chull(data.scores[data.scores$Site == "Dune Centre", c("NMDS1", "NMDS2")]),]
grp.gpco <- data.scores[data.scores$Site == "Gravel Plain Control", ][chull(data.scores[data.scores$Site == "Gravel Plain Control", c("NMDS1", "NMDS2")]),]
grp.gpce <- data.scores[data.scores$Site == "Gravel Plain Centre", ][chull(data.scores[data.scores$Site == "Gravel Plain Centre", c("NMDS1", "NMDS2")]),]

hull.data <- rbind(grp.dco, grp.dce, grp.gpco, grp.gpce)
hull.data

attach(mapping_file)
#plot using ggplot

#set the colors and shapes
cols = c("#1b4f72", "#aed6f1", "#1d8348","#58d68d", "#943126",  "#ec7063")
#change to the ones which has the black edges
shps = c(22, 21, 24)

plot_1 <- ggplot() +
  #geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Site, group = Site), alpha = 0.30) +
  #scale_fill_brewer(type = div, palette = 'Set1')+
  geom_point(data = data.scores, aes(x = NMDS1, y = NMDS2, fill = Loc_dep, shape = Source), size = 6) +  #add your points
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
  annotate("text", x = 0.15, y = -0.55, label = "stress = 0.14", size = 7) +
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
svg(filename = "16S_1405_nonam.svg", width = 12, height = 6, pointsize = 12)
plot_1
dev.off()

#THIS IS FOR LESS DEEP
#THIS IS THE SCRIPT USED FOR ONLY LOOKING AT THE DIFFS BETWEEN AUS AND NAM

#RAREFIED TO 13208 

otu_table_16S <- read.delim("./16S/rarefaction/summarize_taxa13208/FC16S_tarare13208_L6.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata.txt")


otu_transf <- decostand(otu_table_16S, "hellinger")

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
grp.dco <- data.scores[data.scores$Site == "Dune Control", ][chull(data.scores[data.scores$Site == "Dune Control", c("NMDS1", "NMDS2")]),]
grp.dce <- data.scores[data.scores$Site == "Dune Centre", ][chull(data.scores[data.scores$Site == "Dune Centre", c("NMDS1", "NMDS2")]),]
grp.gpco <- data.scores[data.scores$Site == "Gravel Plain Control", ][chull(data.scores[data.scores$Site == "Gravel Plain Control", c("NMDS1", "NMDS2")]),]
grp.gpce <- data.scores[data.scores$Site == "Gravel Plain Centre", ][chull(data.scores[data.scores$Site == "Gravel Plain Centre", c("NMDS1", "NMDS2")]),]

hull.data <- rbind(grp.dco, grp.dce, grp.gpco, grp.gpce)
hull.data
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
  annotate("text", x = 0.35, y = -0.55, label = "stress = 0.13", size = 8) +
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
svg(filename = "16S_13208_nmds.svg", width = 12, height = 6, pointsize = 12)
plot_1
dev.off()


#use multiplot
multiplot(plot_1, plot_2, plot_3, cols = 2)
