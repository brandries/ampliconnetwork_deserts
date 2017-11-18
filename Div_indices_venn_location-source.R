library(vegan)

#import data
setwd("C:/Users/Andries van der Walt/Dropbox/Masters 2017/Thesis/Chapter 2 - FC interactions/FC data/")

#THIS SCRIPT IS TOO DIFFICULT TO REWRITE, SO JUST CHANGE THE TABLE FOR 16S INPUT
otu_table_ITS  <- read.delim("./ITS/rarefaction/FCITS_rare6767.txt", row.names = 1)
mapping_file <- read.delim("./ITS/metadata.txt")
#the infamous 16s table
otu_table_ITS <- read.delim("./16S/rarefaction/FC16S_tarare13208_mrdna_reorder.txt", row.names = 1)
mapping_file <- read.delim("./16S/metadata_mrdna.txt")

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

alpha.mean <-tapply(specnumber(otu_table), Loc_sour, mean)
alpha.sd <-tapply(specnumber(otu_table), Loc_sour, sd)
names(alpha) <-c("alpha1")
attach(alpha)
ANOVA1 <-aov(alpha1~Loc_sour)
summary(ANOVA1)
TukeyHSD(ANOVA1)

gamma <-specnumber(otu_table, Loc_sour)
gamma


#I think it is the easiest to do this manually
beta.whittaker <-gamma/alpha.mean

beta_ausce <- gamma[1]/alpha[1:18,]
beta_ausco <- gamma[2]/alpha[19:36,]
beta_ausma <- gamma[3]/alpha[37:54,]
beta_dunce <- gamma[4]/alpha[55:64,]
beta_dunco <- gamma[5]/alpha[65:74,]
beta_gpce <- gamma[6]/alpha[75:84,]
beta_gpco <- gamma[7]/alpha[85:94,]

beta_all <-data.frame(c(beta_ausce, beta_ausco, beta_ausma, beta_dunce, beta_dunco, beta_gpce, beta_gpco), Loc_sour)
names(beta_all)<-c("beta","location")
anova.beta <-aov(beta ~Loc_sour, data=beta_all)
summary(anova.beta)
TukeyHSD(anova.beta)

#combine all these results
alpha_beta <- cbind(alpha, beta_all)
write.csv(alpha_beta, "alpha_16S_loc_sour.txt")
write.csv(gamma, "gamma_16S_loc_sour.txt")


#Venn diagram
#*****************
library(gplots)

otu_table_t <- t(otu_table)
comb <-data.frame(cbind(rowSums(otu_table_t[,1:18]), rowSums(otu_table_t[,19:36]), rowSums(otu_table_t[,37:54]),
         rowSums(otu_table_t[,55:64]),rowSums(otu_table_t[,65:74]), rowSums(otu_table_t[,75:84]) , 
         rowSums(otu_table_t[,85:94]) )>0)
names(comb)<-c("Australia Centre", "Australia Control", "Australia Margin", "Dune Centre", "Dune Control", 
               "GP Centre", "GP Control")
venn(comb)

#ALTERNATIVE VENNS WITH OTHER GOURPINGS
#only centre
comb <-data.frame(cbind(rowSums(otu_table_t[,1:18]), rowSums(otu_table_t[,55:64]), rowSums(otu_table_t[,75:84]) )>0)
names(comb)<-c("Australia Centre",  "Dune Centre", "GP Centre")
venn(comb)

#only classic sites
comb <-data.frame(cbind(rowSums(otu_table_t[,1:18]), rowSums(otu_table_t[,19:36]), 
                        rowSums(otu_table_t[,55:64]),rowSums(otu_table_t[,65:74]), rowSums(otu_table_t[,75:84]) , 
                        rowSums(otu_table_t[,85:94]) )>0)
names(comb)<-c("Australia Centre", "Australia Control", "Dune Centre", "Dune Control", 
               "GP Centre", "GP Control")
venn(comb)

#nam all
comb <-data.frame(cbind( rowSums(otu_table_t[,55:64]),rowSums(otu_table_t[,65:74]), rowSums(otu_table_t[,75:84]) , 
                        rowSums(otu_table_t[,85:94]) )>0)
names(comb)<-c("Dune Centre", "Dune Control", 
               "GP Centre", "GP Control")
venn(comb)

#nam centres
comb <-data.frame(cbind( rowSums(otu_table_t[,55:64]), rowSums(otu_table_t[,75:84]) )>0)
names(comb)<-c("Dune Centre",  "GP Centre")
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
