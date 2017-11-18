library(vegan)

#import data
setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 3 - FC interactions/FC data/")

#THIS SCRIPT IS TOO DIFFICULT TO REWRITE, SO JUST CHANGE THE TABLE FOR 16S INPUT
otu_table_ITS  <- read.delim("./ITS/rarefaction/FCITS_rare6767.txt", row.names = 1)
mapping_file <- read.delim("./ITS/metadata.txt")
#the infamous 16s table
otu_table_ITS <- read.delim("./16S/rarefaction/FC16S_tarare1405_nonam.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata_all_nonam_aut.txt")

attach(mapping_file)
#minus one to remove the taxonomy or which ever coloms you need to remove
otu_table <-as.data.frame(t(otu_table_ITS[,1:(length(otu_table_ITS)-6)]))


#THIS IS TO REORDER IF YOU HAVE TO 
#otu_table_comb <- cbind(mapping_file,otu_table)
#otu_table <- otu_table_comb[order(otu_table_comb$Location),(length(mapping_file)+1):length(otu_table_comb)]
#mapping_file <- mapping_file[order(mapping_file$Location), ]

#Diversity metrics
#**********************
alpha <-data.frame(specnumber((otu_table)))

alpha.mean <-tapply(specnumber(otu_table), Source, mean)
alpha.sd <-tapply(specnumber(otu_table), Source, sd)
names(alpha) <-c("alpha1")
attach(alpha)
ANOVA1 <-aov(alpha1~Source)
summary(ANOVA1)
TukeyHSD(ANOVA1)

gamma <-specnumber(otu_table, Source)
gamma


#I think it is the easiest to do this manually
beta.whittaker <-gamma/alpha.mean

#make new dataframe with deep and surface and mapping files
centre <- cbind(c(alpha[1:18,],  alpha[55:64,], alpha[75:84,]))
control <- cbind(c(alpha[19:36,], alpha[65:74,], alpha[85:94,]))
margin <- cbind(c(alpha[37:54,]))
centre_m <- rbind(mapping_file[1:18,],  mapping_file[55:64,], mapping_file[75:84,])
control_m <- rbind(mapping_file[19:36,], mapping_file[65:74,], mapping_file[85:94,])
margin_m <- rbind(mapping_file[37:54,])

beta_centre <- gamma[1]/centre
beta_control <- gamma[2]/control
beta_margin <- gamma[3]/margin

beta_mapping <- cbind(c(centre_m$Source,control_m$Source, margin_m$Source))
  
beta_all <-data.frame(c(beta_centre, beta_control, beta_margin), beta_mapping)
names(beta_all)<-c("beta","location")
anova.beta <-aov(beta ~Source, data=beta_all)
summary(anova.beta)
TukeyHSD(anova.beta)

new_beta <- cbind(beta_all, rbind(centre_m,control_m, margin_m))
#combine all these results

write.csv(alpha, "alpha_16S_source.txt")
write.csv(new_beta, "beta_16S_source.txt")
write.csv(gamma, "gamma_16S_source.txt")


#Venn diagram
#*****************
library(gplots)

centre <- rbind(otu_table[1:18,],  otu_table[55:64,], otu_table[75:84,])
control <- rbind(otu_table[19:36,], otu_table[65:74,], otu_table[85:94,])
margin <- rbind(otu_table[37:54,])

otu_table <- rbind(centre, control, margin)
otu_table_t <- t(otu_table)
comb <-data.frame(cbind(rowSums(otu_table_t[,1:38]), rowSums(otu_table_t[,39:76]), rowSums(otu_table_t[,77:94]))>0)
names(comb)<-c("Centre", "Control", "Margin")


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
