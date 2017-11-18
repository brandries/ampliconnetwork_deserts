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

Shannon <- data.frame(diversity(otu_table))
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

beta_aus <- gamma[1]/alpha[1:12,]
beta_dune <- gamma[2]/alpha[13:24,]
beta_gp <- gamma[3]/alpha[25:36,]



beta_all <-data.frame(c(beta_aus, beta_dune, beta_gp), Location)
names(beta_all)<-c("beta","location")
anova.beta <-aov(beta ~location, data=beta_all)
summary(anova.beta)
TukeyHSD(anova.beta)

#combine all these results
alpha_beta <- cbind(alpha,Shannon, beta_all)
write.csv(alpha_beta, "alpha_16S_final.txt")
write.csv(gamma, "gamma_16S_final.txt")


#Venn diagram
#*****************
library(gplots)

otu_table_t <- t(otu_table)
comb <-data.frame(cbind(rowSums(otu_table_t[,1:12]), rowSums(otu_table_t[,13:24]), rowSums(otu_table_t[,25:36]))>0)
names(comb)<-c("Australia", "Dune", "GP")
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
