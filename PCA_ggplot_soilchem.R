library(vegan)
library(ggplot2)
library(RColorBrewer)

setwd("C:/Users/Andries van der Walt/Google Drive/Masters/Data 2017/Soil Chemistry")

soil_chem <- read.table("./soil_chem_transformed.txt", dec = ".", row.names = 1)
mapping_file <- read.delim("./metadata.txt")
#trnasform using the correc t srcipt
soil_dist <- vegdist(soil_chem, method = "euclidean")
pca1 = prcomp(soil_chem, scale. = TRUE)    # Set principal components

attach(mapping_file)
scores = as.data.frame(pca1$x)
summary(pca1)                 #gives cumulative variance explained

# plot of observations
shps = c(22, 21, 24)  
    
plot_pca <- ggplot() +
  geom_point(data = scores, aes(x = PC1, y = PC2, label=rownames(scores), fill = Loc_dep, shape = Source), size = 6) +  #add your points
    scale_fill_brewer(type = "div", palette = "Paired")+
    scale_shape_manual(values = shps) +
    coord_fixed(ratio = 1.85) + 
    theme_bw()  +
    xlab("PC1 (58.67%)") +
    ylab("PC2 (19.92%)") +
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
  
#Test significance
adonis(soil_chem~ Location*Source*Depth)
  
  
#save your plot
svg(filename = "soil_chem_pca.svg", width = 12, height = 6, pointsize = 12)
plot_pca
dev.off()



###########################################
#DO ANOVAS
##########################################

soil_model <- lm(soil_chem$pH~Loc_sour)

aov_soil <- anova(soil_model)


aov_soil2 <- aov(soil_chem$Carbon~Depth)
TukeyHSD(aov_soil2)












#################################################
# to do this in a tsne
#################################################
library(tsne)
compare_soil_chem <- cbind(mapping_file, soil_chem)
#TRANSPOSE YOUR MARIX

biol.data <- compare_soil_chem[,7:21]


#this is a function you set, one for number of colors
#the other for the names
#then you plot the tsne and you get it to giev you a plot at each iteration

colors = brewer.pal(length(unique(compare_soil_chem$Loc_dep)), "Paired")
names(colors) = unique(compare_soil_chem$Loc_dep)
ecb = function(x,y){ plot(x,t='p', col=colors[compare_soil_chem$Loc_dep], lwd=10) }

#TSNE intro. set perplexity based on your sample size, and how big you expect the groups to be. This is the perplexity value.
#rerun at high iteration for final figures
Sys.sleep(10); tsne_iris = tsne(biol.data, epoch_callback = ecb, perplexity=7, max_iter = 7000, epoch = 10)
#with this you can make a recording of the tsne

plot(1, type="n", axes=F, xlab="", ylab="")
