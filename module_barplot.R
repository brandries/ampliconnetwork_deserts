#Barplots for taxonomic composition
libs <- c("ggplot2", "ggthemes", "extrafont", "plyr", "scales", "wesanderson", "reshape")
lapply(libs, require, character.only = TRUE)

#Import data
#Your data has to be a table with first column being the one of axis, and the second being the other axis.
#the third column contains the values used to draw the bar chart. 
#alternatively, you can input a OTU table style table, and change the format using the reshape::melt function.
setwd("~/../Dropbox/Masters 2017/Thesis/Chapter 3 - FC interactions/FC data/Combined ITS-16S MENA/Nam network files/")
module_abundance <- read.delim("./Control_nam_dune_module_abundance.txt")
module_abundance <- melt(module_abundance)
names(module_abundance) <- c("Phylum", "ModuleNum", "Count")


#create palette for plot
filcols <- rainbow(25)
attach(module_abundance)
#create plot with barplot, x and y labels, remove background and change theme to black lines
plot1 <- ggplot() + 
  geom_bar(aes(y = Count, x = ModuleNum, fill = Phylum), colour = "black", data = module_abundance, stat = "identity") +
  labs(x = "Module number", y = "Relative abundance (%)") +
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
  coord_flip()

plot1


#plot in svg
svg(filename = "barplot_module_abundance_dune_control.svg", width = 10, height = 5, pointsize = 12)
plot1
dev.off()
