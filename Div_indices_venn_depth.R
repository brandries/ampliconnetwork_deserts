library(vegan)

#import data
setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 3 - FC interactions/FC data/")

#THIS SCRIPT IS TOO DIFFICULT TO REWRITE, SO JUST CHANGE THE TABLE FOR 16S INPUT
otu_table_ITS  <- read.delim("./ITS/rarefaction/FCITS_rare6767_reorder_subsample3_nomar.txt", row.names = 1)
mapping_file <- read.delim("./ITS/metadata_reorder_subsample3_nomar.txt")
#the infamous 16s table
otu_table_ITS <- read.delim("./16S/rarefaction/FC16S_tarare13208_mrdna_reorder_subsample3_nomargin.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata_mrdna_subsample3_nomargin.txt")

attach(mapping_file)
#minus one to remove the taxonomy or which ever coloms you need to remove
otu_table <-as.data.frame(t(otu_table_ITS[,1:(length(otu_table_ITS))]))


#THIS IS TO REORDER IF YOU HAVE TO 
#otu_table_comb <- cbind(mapping_file,otu_table)
#otu_table <- otu_table_comb[order(otu_table_comb$Location),(length(mapping_file)+1):length(otu_table_comb)]
#mapping_file <- mapping_file[order(mapping_file$Location), ]

#Diversity metrics
#**********************
alpha <-data.frame(specnumber((otu_table)))

alpha.mean <-tapply(specnumber(otu_table), Depth, mean)
alpha.sd <-tapply(specnumber(otu_table), Depth, sd)
names(alpha) <-c("alpha1")
attach(alpha)
ANOVA1 <-aov(alpha1~Depth)
summary(ANOVA1)
TukeyHSD(ANOVA1)

gamma <-specnumber(otu_table, Depth)
gamma


#I think it is the easiest to do this manually
beta.whittaker <-gamma/alpha.mean

#make new dataframe with deep and surface and mapping files
deep <- cbind(c(alpha[1:3,], alpha[7:9,], alpha[13:15,], alpha[19:21,],alpha[25:27,],alpha[31:33,]))
surface <- cbind(c(alpha[4:6,], alpha[10:12,], alpha[16:18,], alpha[22:24,],alpha[28:30,],alpha[34:36,]))
deep_m <- rbind(mapping_file[1:3,], mapping_file[7:9,], mapping_file[13:15,], mapping_file[19:21,],mapping_file[25:27,],mapping_file[31:33,])
Surface_m <- rbind(mapping_file[4:6,], mapping_file[10:12,], mapping_file[16:18,], mapping_file[22:24,],mapping_file[28:30,],mapping_file[34:36,])

beta_deep <- gamma[1]/deep
beta_surface <- gamma[2]/surface

beta_mapping <- cbind(c(deep_m$Depth,Surface_m$Depth))
  
beta_all <-data.frame(c(beta_deep, beta_surface), beta_mapping)
names(beta_all)<-c("beta","location")
anova.beta <-aov(beta ~Depth, data=beta_all)
summary(anova.beta)
TukeyHSD(anova.beta)

new_beta <- cbind(beta_all, rbind(deep_m,Surface_m))
#combine all these results

write.csv(alpha, "alpha_16S_depth.txt")
write.csv(new_beta, "beta_16S_depth.txt")
write.csv(gamma, "gamma_16S_depth.txt")


#Venn diagram
#*****************
library(gplots)


otu_table <- t(otu_table)
deep <- rbind(otu_table[1:3,], otu_table[7:9,], otu_table[13:15,], otu_table[19:21,],otu_table[25:27,],otu_table[31:33,])
surface <- rbind(otu_table[4:6,], otu_table[10:12,], otu_table[16:18,], otu_table[22:24,],otu_table[28:30,],otu_table[34:36,])
otu_table <- rbind(deep, surface)
otu_table_t <- t(otu_table)

comb <-data.frame(cbind(rowSums(otu_table_t[,1:18]), rowSums(otu_table_t[,19:36]))>0)
names(comb)<-c("Deep", "Surface")
venn(comb)

#Or draw using venneuler
#*************************
library(venneuler)
plot<- venneuler(comb > 0)
plot(plot)


#Indicator species
#***********************

library(labdsv)

indval <-indval(rarefied_34000_ITS_ws, area, too.many = 400)
summary(indval, p=0.05)		#p-value
