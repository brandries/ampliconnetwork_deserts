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
centre <- cbind(c(alpha[1:6,],  alpha[13:18,], alpha[25:30,]))
control <- cbind(c(alpha[7:12,], alpha[19:24,], alpha[31:36,]))

centre_m <- rbind(mapping_file[1:6,],  mapping_file[13:18,], mapping_file[25:30,])
control_m <- rbind(mapping_file[7:12,], mapping_file[19:24,], mapping_file[31:36,])


beta_centre <- gamma[1]/centre
beta_control <- gamma[2]/control


beta_mapping <- cbind(c(centre_m$Source,control_m$Source))
  
beta_all <-data.frame(c(beta_centre, beta_control), beta_mapping)
names(beta_all)<-c("beta","location")
anova.beta <-aov(beta ~Source, data=beta_all)
summary(anova.beta)
TukeyHSD(anova.beta)

new_beta <- cbind(beta_all, rbind(centre_m,control_m))
#combine all these results

write.csv(alpha, "alpha_16S_source.txt")
write.csv(new_beta, "beta_16S_source.txt")
write.csv(gamma, "gamma_16S_source.txt")


#Venn diagram
#*****************
library(gplots)

centre <- rbind(otu_table[1:6,],  otu_table[13:18,], otu_table[25:30,])
control <- rbind(otu_table[7:12,], otu_table[19:24,], otu_table[31:36,])
margin <- rbind(otu_table[37:54,])

otu_table <- rbind(centre, control)
otu_table_t <- t(otu_table)
comb <-data.frame(cbind(rowSums(otu_table_t[,1:18]), rowSums(otu_table_t[,19:36]))>0)
names(comb)<-c("Centre", "Control")


venn(comb)

#################
#VENN FOR CENTRE ONLY
##################
auscentre <- otu_table[1:6,]
dunecentre <- otu_table[13:18,]
gpcentre <- otu_table[25:30,]
auscontrol <- otu_table[7:12,]
dunecontrol <- otu_table[19:24,]
gpcontrol <- otu_table[31:36,]
otu_table <- rbind(auscentre, dunecentre, gpcentre, auscontrol, dunecontrol, gpcontrol)
otu_table_t <- t(otu_table)
comb <-data.frame(cbind(rowSums(otu_table_t[,1:6]), rowSums(otu_table_t[,7:12]), rowSums(otu_table_t[,13:18]), rowSums(otu_table_t[,19:24]), 
                        rowSums(otu_table_t[,25:30]), rowSums(otu_table_t[,31:36]))>0)
names(comb)<-c("Aus", "Dune", "GP")


write.csv(comb, "venn_6 lobes_16s.txt")
venn(comb)

#Or draw using venneuler
#*************************
library(venneuler)
plot<- venneuler(comb > 0)
plot(plot)

#VENN DIAGRAM
#######################
library(VennDiagram)
venn.diagram(comb, "venn_6lobes.tiff")

#Indicator species
#***********************

library(labdsv)

indval <-indval(rarefied_34000_ITS_ws, area, too.many = 400)
summary(indval, p=0.05)		#p-value
