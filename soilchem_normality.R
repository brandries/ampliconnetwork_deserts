libs <- c("ggplot2", "ggthemes", "BiodiversityR", "vegan", "reshape")
lapply(libs, require, character.only = TRUE)

setwd("C:/Users/Andries van der Walt/Documents/Masters/Data 2017/Soil Chemistry")

soil_chem <- read.delim("./soil_chem_final.txt", dec = ".", row.names = 1)
mapping_file <- read.delim("./metadata.txt")

#do this for each variable
with(soil_chem, qqPlot(pH, dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
shapiro.test(soil_chem$pH)
#not transformation

with(soil_chem, qqPlot(EC..mS.m., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$EC..mS.m. <- decostand(soil_chem$EC..mS.m., method = "log")
shapiro.test(soil_chem$EC..mS.m.)
#improved normality

with(soil_chem, qqPlot(soil_chem$P..mg.l., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$P..mg.l. <- decostand(soil_chem$P..mg.l., method = "log")
shapiro.test(soil_chem$P..mg.l.)
#improved normality

with(soil_chem, qqPlot(soil_chem$Na..mg.l., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$Na..mg.l. <- decostand(soil_chem$Na..mg.l., method = "log")
shapiro.test(soil_chem$Na..mg.l.)
#improved normality

with(soil_chem, qqPlot(soil_chem$K..mg.l., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$K..mg.l. <- decostand(soil_chem$K..mg.l., method = "log")
shapiro.test(soil_chem$K..mg.l.)
#improved normality

with(soil_chem, qqPlot(soil_chem$Ca..mg.l., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$Ca..mg.l. <- decostand(soil_chem$Ca..mg.l., method = "log")
shapiro.test(soil_chem$Ca..mg.l.)
#improved normality

with(soil_chem, qqPlot(soil_chem$Mg..mg.l., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$Mg..mg.l. <- decostand(soil_chem$Mg..mg.l., method = "log")
shapiro.test(soil_chem$Mg..mg.l.)
#improved normality

with(soil_chem, qqPlot(soil_chem$Cl..mg.l., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$Cl..mg.l. <- decostand(soil_chem$Cl..mg.l., method = "log")
shapiro.test(soil_chem$Cl..mg.l.)
#improved normality

with(soil_chem, qqPlot(soil_chem$SO4..mg.., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$SO4..mg.. <- decostand(soil_chem$SO4..mg.., method = "log")
shapiro.test(soil_chem$SO4..mg..)
#improved normality

with(soil_chem, qqPlot(soil_chem$NH4.N..mgl., dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$NH4.N..mgl. <- decostand(soil_chem$NH4.N..mgl., method = "log")
shapiro.test(soil_chem$NH4.N..mgl.)
#improved normality

with(soil_chem, qqPlot(soil_chem$NO3.N, dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$NO3.N <- decostand(soil_chem$NO3.N, method = "log")
shapiro.test(soil_chem$NO3.N)
#improved normality

with(soil_chem, qqPlot(soil_chem$Sand, dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
shapiro.test(soil_chem$Sand)
#improved normality

with(soil_chem, qqPlot(soil_chem$Clay, dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
shapiro.test(soil_chem$Clay)
#improved normality

with(soil_chem, qqPlot(soil_chem$Silt, dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
shapiro.test(soil_chem$Silt)
#improved normality

with(soil_chem, qqPlot(soil_chem$Carbon, dist="norm", id.method="y", id.n=2, labels=rownames(soil_chem)))
soil_chem$Carbon <- decostand(soil_chem$Carbon, method = "log")
shapiro.test(soil_chem$Carbon)

write.table(soil_chem, "soil_chem_transformed.txt")


