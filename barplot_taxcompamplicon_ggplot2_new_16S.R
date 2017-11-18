#Barplots for taxonomic composition
libs <- c("ggplot2", "ggthemes", "extrafont", "plyr", "scales", "wesanderson", "reshape")
lapply(libs, require, character.only = TRUE)

#Import data
#Your data has to be a table with first column being the one of axis, and the second being the other axis.
#the third column contains the values used to draw the bar chart. 
#alternatively, you can input a OTU table style table, and change the format using the reshape::melt function.

setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 3 - FC interactions/FC data")


#ITS
tax <- read.delim("./ITS/rarefaction/summarize_taxa/FCITS_rare6767_L2_ave.txt", row.names = 1)
mapping_file <- read.delim("./ITS/metadata_tax.txt")
tax <- cbind(mapping_file, tax)

#16S
tax <- read.delim("./16S/rarefaction/summarize_taxa13208/FC16S_tarare13208_L2_top20_ave.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata_tax.txt")
tax <- cbind(mapping_file, tax)


m_tax <- melt(tax, id.vars = c("SampleID", "Name", "Location", "Source", "Loc_sour", "Loc_dep", "Depth", "sequencing", "pseudorep", "Repnumber"))


#order the figure - pur in the phyla in your own dataset
#corr_order <- c("Other", "Ascomycota", "Acidobacteria", "Actinobacteria ", "Bacteroidetes", "Chloroflexi", "Cyanobacteria",
#                "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota")
#corr_nam <- c("Other", "Ascomycota", "Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chloroflexi", "Cyanobacteria",
#                "Firmicutes", "Gemmatimonadetes", "Planctomycetes", "Proteobacteria", "Thaumarchaeota")
#compare_genes_tax$Phylum <- factor(compare_genes_tax$Phylum, levels = corr_order, 
#                            labels = corr_nam)

#create palette for plot
filcols <- wes_palette("Darjeeling", 21, type = c("continuous"))
filcols <- rainbow(19)
#create plot with barplot, x and y labels, remove background and change theme to black lines
plot1 <- ggplot() + 
  geom_bar(aes(y = value, x = SampleID, fill = variable), colour = "black", data = m_tax, stat = "identity") +
  labs(x = "Sample Site", y = "Relative abundance (%)") +
  scale_fill_manual(values = filcols) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 10, angle = 90, hjust = 1),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())

plot1
#plot in jpeg
svg(filename = "16S_barplot.svg", width = 12, height = 8, pointsize = 12)
plot1
dev.off()
tr