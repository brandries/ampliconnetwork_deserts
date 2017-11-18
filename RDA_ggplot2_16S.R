library(vegan)
library(ggplot2)


setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 2 - FC interactions/FC data/")


#RAREFIED TO 1405 

otu_table_16S <- read.delim("./16S/rarefaction/FC16S_tarare13208_mrdna_reorder_nomargin.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata_nomargin_mrdna_ordered.txt")
  
#Import soilchem
setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/Soil Chemistry")
soil_chem <- read.delim("./soil_chem_final_all_reorder.txt", dec = ".", row.names = 1)

otu_transf <- decostand(otu_table_16S, "hellinger")
soil_chem <- decostand(soil_chem, "log")
###########################################################
#DO RDA
###########################################################

#construct rda points
rda.biol <- rda(t(otu_transf), soil_chem)

smry <- summary(rda.biol)

df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2
df1 <- cbind(df1, mapping_file)


#check if your vectors are significant
envfit(t(otu_transf), soil_chem)#output gives you levels of significance

#set new sets of variables for vectors
rda(biol.trans, working_FCSOIL_CHEM_log[,c(1, 3, 4, 6,  9, 12, 13, 14, 15)]) -> rda.new.vars

#set vectors
df2  <- data.frame(smry$biplot[,1:2]) 
vect_soil <- df2[c(1,3,4,6,9,10,12,13,14), ]  #select only those that are significant


#set the colors and shapes
cols = c("#1b4f72", "#aed6f1", "#1d8348","#58d68d", "#943126",  "#ec7063")
#change to the ones which has the black edges
shps = c(22, 21, 24)


attach(mapping_file)
#draw rda plot points
plot_2 <- ggplot(df1, aes(x=RDA1, y=RDA2)) + 
  geom_point(aes(label=rownames(df1), fill = Loc_dep, shape = Source),size=6) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()+
  guides(colour = FALSE) +
  guides(shape = FALSE) +
  #  guides(colour = guide_legend(title = NULL))+
  #  guides(shape = guide_legend(title = NULL))+
  scale_fill_brewer(type = "div", palette = "Paired")+
  scale_shape_manual(values = shps) +
  theme(axis.text.x = element_blank(), #most of arguments to follow is to remove features of the graph
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 18), #increase the size of x label text
        axis.title.y = element_text(size = 18), #increase size of y label text
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18))+
  #Add the vectors
  geom_segment(data=df2, aes(x=0, xend=RDA1, y=0, yend=RDA2), size = 0.5,
               color="navy", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df2, 
            aes(x=RDA1,y=RDA2,label=rownames(df2),
                hjust=0.2*(1-sign(RDA1)),vjust=0.5*(1-sign(RDA2))), 
            color="navy", size=5)


svg(filename = "rda_16S.svg", width = 12, height = 6, pointsize = 12)
plot_2
dev.off()

#save your plot
jpeg(filename = "archaea_rda.jpeg", width = 6000, height = 5000, res = 600)
plot_rda
dev.off()


# to get the amount of variation explained
head(summary(rda.new.vars))  #the proportion constrained variance is the variance explained
RsquareAdj(rda.biol)       #get R^2 and R^2 adj     


#############################################################################################




