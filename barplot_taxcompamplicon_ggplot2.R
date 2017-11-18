#Barplots for taxonomic composition
libs <- c("ggplot2", "ggthemes", "extrafont", "plyr", "scales", "wesanderson")
lapply(libs, require, character.only = TRUE)

#Import data
#Your data has to be a table with first column being the one of axis, and the second being the other axis.
#the third column contains the values used to draw the bar chart. 
#alternatively, you can input a OTU table style table, and change the format using the reshape::melt function.

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/CHapter 2 - FC interactions/AMplicon analysis/otu_tables/16S/summarized_rarefied")

#for old
setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/CHapter 2 - FC interactions/AMplicon analysis/otu_tables/16S/Raw otus/old_qiime")


tax <- read.delim("nam_16S_rare_L2_newnames.txt")
tax_old <- read.delim("filtered_otu_with_tax_L2.txt")

melt(tax) -> m_tax

#set up locations 
location <- c(rep("Dune", 390), rep(c(rep("GP Deep", 78), rep("GP Surface", 78)), 5))

location <- as.data.frame(location)

#set up source
source <- c(rep(c(rep("FC", 39), rep("Control", 39)),15))
source <- as.data.frame(source)

#set up depth
depth <- c(rep("Surface", 390), rep(c(rep("Deep", 78), rep("Surface", 78)), 5))

depth <- as.data.frame(depth)


#add to dataframe
m_tax <- cbind(m_tax,location, source, depth)

m_tax$location <- factor(m_tax$location, levels = c("Dune", "GP Deep", "GP Surface"))
m_tax$source <- factor(m_tax$source, levels = c("FC", "Control"))
m_tax$depth <- factor(m_tax$depth, levels = c("Surface", "Deep"))
##########################################################
#FOR THE ODL ONE
##########################################################

melt(tax_old) -> m_tax_old

#set up locations 
location <- c(rep("Dune", 200), rep(c(rep("GP Deep", 40), rep("GP Surface", 40)), 5))

location <- as.data.frame(location)

#set up source
source <- c(rep(c(rep("FC", 20), rep("Control", 20)),15))
source <- as.data.frame(source)

#set up depth
depth <- c(rep("Surface", 200), rep(c(rep("Deep", 40), rep("Surface", 40)), 5))

depth <- as.data.frame(depth)


#add to dataframe
m_tax_old <- cbind(m_tax_old,location, source, depth)

m_tax$location <- factor(m_tax$location)
#order the figure - pur in the phyla in your own dataset
#corr_order <- c("Other", "Ascomycota", "Acidobacteria", "Actinobacteria ", "Bacteroidetes", "Chloroflexi", "Cyanobacteria",
#                "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota")
#corr_nam <- c("Other", "Ascomycota", "Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Cyanobacteria",
#                "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota")
#compare_genes_tax$Phylum <- factor(compare_genes_tax$Phylum, levels = corr_order, 
#                            labels = corr_nam)

#create palette for plot
filcols <- wes_palette("Darjeeling", 39, type = c("continuous"))
filcols <- scale_color_manual()
#create plot with barplot, x and y labels, remove background and change theme to black lines
plot1 <- ggplot() + 
  geom_bar(aes(y = value, x = variable, fill = X.OTU.ID), colour = "black", data = m_tax, stat = "identity") +
  labs(x = "Sample Site", y = "Relative abundance (%)") +
  scale_fill_manual(values = filcols) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
facet_grid(~location)

#plot in jpeg
jpeg(filename = "barplot_genes_tax_final.jpeg", width = 7000, height = 5000, res = 800)
plot1
dev.off()
