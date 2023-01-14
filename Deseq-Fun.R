#Reset R's Brain
rm(list=ls())

#Set working directory............................................................
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence/abundance_barplot")

#Load Packages-------------------------------------------------------------------------------
library(qiime2R)
library(DESeq2)#log2fold analysis
library(ggplot2)#Plotting
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(ochRe)# color paletteibrary(multcomp)#multiple comparisons--for glm 
library(tidyverse)
library(tidyr)

#LOAD QIIME DATA-RAREFIED ---------------------------------------------------------------------
Metadata<-read_tsv("metadata_fungi-merged_NC-2.tsv")#Rare metadata
Table<-read_qza("rarefied_table.qza", tmp = "C:/tmp")#Rarefied table
Tree<-read_qza("Fungal-rooted-tree-NC.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Fungal-Taxonomy.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
            c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))


  
#CREATE PHYLOSEQ ARTIFACT-----------------------------------------------------------------------------
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("sample-id"))
)



#----
#----
##----
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("k__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("p__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("c__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("o__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("f__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("g__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("s__", "", tax_table(physeq)[, "Species"])


sample_data(physeq)$severity<-factor(sample_data(physeq)$severity) # make as factor

#SUBSAMPLE DATA TO LOOK ONLY AT SOIL SAMPLES ----------------------------------------------------------------------
physeqSoil<-subset_samples(physeq,severity!="low");physeqSoil#11595x302
sample_names(physeqSoil); sample_data(physeqSoil)$severity

#physeqRain<-subset_samples(physeq,SampleType=="Rain");physeqRain#11595x33
#sample_names(physeqRain); sample_data(physeqRain)$TSFdays


#Look at the data to make sure there are no "none" or NA = remove
head(sample_data(physeqSoil)$severity, 50)
head(sample_data(physeqSoil)$fire, 50)
head(sample_data(physeqSoil)$depth, 50)

#Check the structure of the data------------------------------------------------------
str(sample_data(physeqSoil)$severity)#chr
str(sample_data(physeqSoil)$Plot)#num
str(sample_data(physeqSoil)$Timepoint)#chr

#Convert all variables to objects as required per Deseq-------------------------------
sample_data(physeqSoil)$Plot<-factor(sample_data(physeqSoil)$Plot)

#Ran in case there are any n/a in data------------------------------------------------
physeqSoil<-subset_samples(physeqSoil, Treatment!= "None")
physeqSoil<-subset_samples(physeqSoil, Plot!= "None")
physeqSoil<-subset_samples(physeqSoil, Timepoint!= "None")



#
#
#*********************************************************************************************-----
#--- SUBSET DATA BY Site FOR ANALYSIS --------------------------------------------------------
#*********************************************************************************************-----

#Subset samples per Fire to see what the data looks like at T1--------------
physeqSoilCP<-subset_samples(physeqSoil,fire=="church park")
sample_names(physeqSoilCP)

physeqSoilBV<-subset_samples(physeqSoil,fire=="beaver creek")
sample_names(physeqSoilBV)

physeqSoilBC<-subset_samples(physeqSoil,fire=="badger creek")
sample_names(physeqSoilBC)

physeqSoilR<-subset_samples(physeqSoil,fire=="ryan")
sample_names(physeqSoilR)

physeqSoilM<-subset_samples(physeqSoil,fire=="mullen")
sample_names(physeqSoilM)



#
#
#****************************************************************************************************-----
#--- RUN DESEQ ANALYSIS TREATMENT PER Site----------------------------------------------------
#****************************************************************************************************-----
#

#Church Park **********************************************************************----
head(sample_data(physeqSoilCP)$severity, 35)

#Ran in case there are any n/a in data
physeqSoilCP<-subset_samples(physeqSoilCP, severity!= "None")

#Import phyloseq data to create a deseq object.....................
dsCP<-phyloseq_to_deseq2(physeqSoilCP, ~ severity)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeansCP = apply(counts(dsCP), 1, gm_mean1)
dsCP<-estimateSizeFactors(dsCP, geoMeans = geoMeansCP)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#fdsCP<-rowSums(counts(dsCP, normalize = TRUE) >= 200)>=2
#dsCP <- dsCP[fdsCP,]

#Set reference level...............................................
dsCP$severity<-relevel(dsCP$severity,ref="control")
levels(dsCP$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
dsCP<-DESeq(dsCP, fitType="local")



#Beaver Creek**********************************************************************----
head(sample_data(physeqSoilBV)$severity, 35)

#Ran in case there are any n/a in data
physeqSoilBV<-subset_samples(physeqSoilBV, severity!= "None")

#Import phyloseq data to create a deseq object...................
dsBV<-phyloseq_to_deseq2(physeqSoilBV, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeansBV = apply(counts(dsBV), 1, gm_mean2)
dsBV<-estimateSizeFactors(dsBV, geoMeans = geoMeansBV)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#fdsBV<-rowSums(counts(dsBV, normalize = TRUE) >= 200)>=2
#dsBV <- dsBV[fdsBV,]

#Set reference level.............................................
dsBV$severity<-relevel(dsBV$severity,ref="control")
levels(dsBV$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
dsBV<-DESeq(dsBV, fitType="local")


#Badger Creek **********************************************************************----

head(sample_data(physeqSoilBC)$severity, 35)

#Ran in case there are any n/a in data
physeqSoilBC<-subset_samples(physeqSoilBC, severity!= "None")

#Import phyloseq data to create a deseq object...................
dsBC<-phyloseq_to_deseq2(physeqSoilBC, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean3 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeansBC = apply(counts(dsBC), 1, gm_mean3)
dsBC<-estimateSizeFactors(dsBC, geoMeans = geoMeansBC)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#fdsBC<-rowSums(counts(dsBC, normalize = TRUE) >= 200)>=2
#dsBC <- dsBC[fdsBC,]

#Set reference level.............................................
dsBC$severity<-relevel(dsBC$severity,ref="control")
levels(dsBC$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
dsBC<-DESeq(dsBC, fitType="local")


#Ryan **********************************************************************----
head(sample_data(physeqSoilR)$severity, 35)

#Ran in case there are any n/a in data
physeqSoilR<-subset_samples(physeqSoilR, severity!= "None")

#Import phyloseq data to create a deseq object...................
dsR<-phyloseq_to_deseq2(physeqSoilR, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean4 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeansR = apply(counts(dsR), 1, gm_mean4)
dsR<-estimateSizeFactors(dsR, geoMeans = geoMeansR)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#fdsR<-rowSums(counts(dsR, normalize = TRUE) >= 200)>=2
#dsR <- dsR[fdsR,]

#Set reference level.............................................
dsR$severity<-relevel(dsR$severity,ref="control")
levels(dsR$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
dsR<-DESeq(dsR, fitType="local")


#Mullen**********************************************************************----
head(sample_data(physeqSoilM)$severity, 35)

#Ran in case there are any n/a in data
physeqSoilM<-subset_samples(physeqSoilM, severity!= "None")

#Import phyloseq data to create a deseq object...................
dsM<-phyloseq_to_deseq2(physeqSoilM, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean5 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeansM = apply(counts(dsM), 1, gm_mean5)
dsM<-estimateSizeFactors(dsM, geoMeans = geoMeansM)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#fdsM<-rowSums(counts(dsM, normalize = TRUE) >= 200)>=2
#dsM <- dsM[fdsM,]

#Set reference level.............................................
dsM$severity<-relevel(dsM$severity,ref="control")
levels(dsM$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
dsM<-DESeq(dsM, fitType="local")


#
#****************************************************************************************************-----
#--- RUN RESULT ------------------------------------------------------------------------------------------
#****************************************************************************************************-----


#Church Park **********************************************************************----
resCP<-results(dsCP)
resCP<-resCP[order(resCP$padj, na.last=NA), ]
alpha<-0.001
sigtabCP<-resCP[(resCP$padj < alpha), ];sigtabCP
sigtabCP<-cbind(as(sigtabCP, "data.frame"), 
          as(tax_table(physeqSoilCP)[rownames(sigtabCP), ], "matrix"));sigtabCP

#Create ta table of most significant based on alpha and log2fold change.....
posigtabCP<-sigtabCP[sigtabCP[, "log2FoldChange"] > 0, ]
posigtabCP<-posigtabCP[, c("baseMean", "log2FoldChange", "lfcSE", "padj", 
            "Phylum", "Class", "Family", "Genus","Species")];posigtabCP


#Beaver Creek *********************************************************************----
resBV<-results(dsBV)
resBV<-resBV[order(resBV$padj, na.last=NA), ]
alpha<-0.001
sigtabBV<-resBV[(resBV$padj < alpha), ];sigtabBV
sigtabBV<-cbind(as(sigtabBV, "data.frame"), 
                as(tax_table(physeqSoilBV)[rownames(sigtabBV), ], "matrix"));sigtabBV

#Create ta table of most significant based on alpha and log2fold change.....
posigtabBV<-sigtabBV[sigtabBV[, "log2FoldChange"] > 0, ]
posigtabBV<-posigtabBV[, c("baseMean", "log2FoldChange", "lfcSE", "padj", 
                           "Phylum", "Class", "Family", "Genus","Species")];posigtabBV


#Badger Creek3 **********************************************************************----
resBC<-results(dsBC)
resBC<-resBC[order(resBC$padj, na.last=NA), ]
alpha<-0.001
sigtabBC<-resBC[(resBC$padj < alpha), ];sigtabBC
sigtabBC<-cbind(as(sigtabBC, "data.frame"), 
                as(tax_table(physeqSoilBC)[rownames(sigtabBC), ], "matrix"));sigtabBC

#Create ta table of most significant based on alpha and log2fold change.....
posigtabBC<-sigtabBC[sigtabBC[, "log2FoldChange"] > 0, ]
posigtabBC<-posigtabBC[, c("baseMean", "log2FoldChange", "lfcSE", "padj", 
                           "Phylum", "Class", "Family", "Genus","Species")];posigtabBC


#Ryan**********************************************************************----
resR<-results(dsR)
resR<-resR[order(resR$padj, na.last=NA), ]
alpha<-0.001
sigtabR<-resR[(resR$padj < alpha), ];sigtabR
sigtabR<-cbind(as(sigtabR, "data.frame"), 
                as(tax_table(physeqSoilR)[rownames(sigtabR), ], "matrix"));sigtabR

#Create ta table of most significant based on alpha and log2fold change.....
posigtabR<-sigtabR[sigtabR[, "log2FoldChange"] > 0, ]
posigtabR<-posigtabR[, c("baseMean", "log2FoldChange", "lfcSE", "padj", 
                           "Phylum", "Class", "Family", "Genus","Species")];posigtabR


#Mullen**********************************************************************----
resM<-results(dsM)
resM<-resM[order(resM$padj, na.last=NA), ]
alpha<-0.001
sigtabM<-resM[(resM$padj < alpha), ];sigtabM
sigtabM<-cbind(as(sigtabM, "data.frame"), 
               as(tax_table(physeqSoilM)[rownames(sigtabM), ], "matrix"));sigtabM

#Create ta table of most significant based on alpha and log2fold change.....
posigtabM<-sigtabM[sigtabM[, "log2FoldChange"] > 0, ]
posigtabM<-posigtabM[, c("baseMean", "log2FoldChange", "lfcSE", "padj", 
                         "Phylum", "Class", "Family", "Genus","Species")];posigtabM

#
#
#****************************************************************************************************-----
#--- EXPORT FILES ---------------------------------------------------------------------------------------
#****************************************************************************************************-----

write.csv(posigtabCP,"posigtabCP.csv")
write.csv(posigtabBV,"posigtabBV.csv")
write.csv(posigtabBC,"posigtabBC.csv")
write.csv(posigtabR,"posigtabR.csv")
write.csv(posigtabM,"posigtabM.csv")

write.csv(sigtabCP,"sigtabCP.csv")
write.csv(sigtabBV,"sigtabBV.csv")
write.csv(sigtabBC,"sigtabBC.csv")
write.csv(sigtabR,"sigtabR.csv")
write.csv(sigtabM,"sigtabM.csv")


dir.create(file.path("1-Analysis/Deseq/Fungi/Tables"), recursive = TRUE)
#write.csv(sigtab1,"1-Analysis/Deseq/Fungi/Tables/sigtab1.csv")
#write.csv(sigtab2,"1-Analysis/Deseq/Fungi/Tables/sigtab2.csv")
#write.csv(sigtab3,"1-Analysis/Deseq/Fungi/Tables/sigtab3.csv")
#write.csv(sigtab4,"1-Analysis/Deseq/Fungi/Tables/sigtab4.csv")
#write.csv(sigtab5,"1-Analysis/Deseq/Fungi/Tables/sigtab5.csv")
#write.csv(sigtab6,"1-Analysis/Deseq/Fungi/Tables/sigtab6.csv")
#write.csv(sigtab7,"1-Analysis/Deseq/Fungi/Tables/sigtab7.csv")
#write.csv(sigtab8,"1-Analysis/Deseq/Fungi/Tables/sigtab8.csv")
#write.csv(sigtab9,"1-Analysis/Deseq/Fungi/Tables/sigtab9.csv")


#write.csv(posigtab1,"1-Analysis/Deseq/Fungi/Tables/posigtab1.csv")
#write.csv(posigtab2,"1-Analysis/Deseq/Fungi/Tables/posigtab2.csv")
#write.csv(posigtab3,"1-Analysis/Deseq/Fungi/Tables/posigtab3.csv")
#write.csv(posigtab4,"1-Analysis/Deseq/Fungi/Tables/posigtab4.csv")
#write.csv(posigtab5,"1-Analysis/Deseq/Fungi/Tables/posigtab5.csv")
#write.csv(posigtab6,"1-Analysis/Deseq/Fungi/Tables/posigtab6.csv")
#write.csv(posigtab7,"1-Analysis/Deseq/Fungi/Tables/posigtab7.csv")
#write.csv(posigtab8,"1-Analysis/Deseq/Fungi/Tables/posigtab8.csv")
#write.csv(posigtab9,"1-Analysis/Deseq/Fungi/Tables/posigtab9.csv")


#
#
#****************************************************************************************************-----
#--- PHYLUM TO GENUS GRAPHS---------------------------------------------------------------------------------------
#****************************************************************************************************-----

#ChurchPark........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgenCP = subset(sigtabCP, !is.na(Genus))
# Phylum order
x = tapply(sigtabgenCP$log2FoldChange, sigtabgenCP$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgenCP$Phylum = factor(as.character(sigtabgenCP$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgenCP$log2FoldChange, sigtabgenCP$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgenCP$Genus = factor(as.character(sigtabgenCP$Genus), levels = names(x))

sigtabgenCP<-subset_samples(sigtabgenCP, Genus!= "unidentified")
  
CP<-ggplot(sigtabgenCP, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
      geom_point(size=4)+ 
     #scale_color_manual(values=c("#b0022a","#7aa0c2")) +
      theme(axis.text.x = element_text(angle = -90, 
      hjust = 0, vjust = 0.5))





#Timepoint 2 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen2 = subset(sigtab2, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Phylum = factor(as.character(sigtabgen2$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Genus = factor(as.character(sigtabgen2$Genus), levels = names(x))


T2<-ggplot(sigtabgen2, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T2


#Timepoint 3 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen3= subset(sigtab3, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Phylum = factor(as.character(sigtabgen3$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Genus = factor(as.character(sigtabgen3$Genus), levels = names(x))


T3<-ggplot(sigtabgen3, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T3




#Timepoint 4 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen4= subset(sigtab4, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Phylum = factor(as.character(sigtabgen4$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Genus = factor(as.character(sigtabgen4$Genus), levels = names(x))


T4<-ggplot(sigtabgen4, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T4


#Timepoint 5 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen5= subset(sigtab5, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen5$log2FoldChange, sigtabgen5$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen5$Phylum = factor(as.character(sigtabgen5$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen5$log2FoldChange, sigtabgen5$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen5$Genus = factor(as.character(sigtabgen5$Genus), levels = names(x))


T5<-ggplot(sigtabgen5, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T5




#Timepoint 6 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen6= subset(sigtab6, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen6$log2FoldChange, sigtabgen6$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen6$Phylum = factor(as.character(sigtabgen6$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen6$log2FoldChange, sigtabgen6$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen6$Genus = factor(as.character(sigtabgen6$Genus), levels = names(x))


T6<-ggplot(sigtabgen6, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T6




#Timepoint 7 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen7= subset(sigtab7, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen7$log2FoldChange, sigtabgen7$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen7$Phylum = factor(as.character(sigtabgen7$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen7$log2FoldChange, sigtabgen7$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen7$Genus = factor(as.character(sigtabgen7$Genus), levels = names(x))


T7<-ggplot(sigtabgen7, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T7




#Timepoint 8 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen8= subset(sigtab8, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen8$log2FoldChange, sigtabgen8$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen8$Phylum = factor(as.character(sigtabgen8$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen8$log2FoldChange, sigtabgen8$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen8$Genus = factor(as.character(sigtabgen8$Genus), levels = names(x))


T8<-ggplot(sigtabgen8, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T8






#
#
#***************************************************************************************-----
#--------------------- EXPORT GRAPHS --------------------------------------------------------
#***************************************************************************************-----
dir.create(file.path("1-Analysis/Deseq/Fungi/Graphs/Timepoints"), recursive = TRUE)

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T1-Treatment.pdf", height=8, width=10)
T1
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T2-Treatment.pdf", height=8, width=10)
T2
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T3-Treatment.pdf", height=8, width=10)
T3
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T4-Treatment.pdf", height=8, width=10)
T4
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T5-Treatment.pdf", height=8, width=10)
T5
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T6-Treatment.pdf", height=8, width=10)
T6
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T7-Treatment.pdf", height=8, width=10)
T7
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T8-Treatment.pdf", height=8, width=10)
T8
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/T9-Treatment.pdf", height=8, width=10)
T9
dev.off()




#panel graphs....................................................................
library(ggpubr)
pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/All-Timepoints.pdf", height=12, width=14)
ggarrange(T1,T2,T3,T4,T5,T6,T7,T8,T9, legend = TRUE, labels = c("T1","T2","T3",
"T4","T5","T6","T7","T8","T9"))
dev.off()


#
#
#*****************************************************************************************----
# CREATE GRAPHS------ FAMILY TO GENUS --------------------------------------------------------
#*****************************************************************************************----


#Timepoint 1 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen1 = subset(sigtab1, !is.na(Genus))
# Family order
x = tapply(sigtabgen1$log2FoldChange, sigtabgen1$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen1$Family = factor(as.character(sigtabgen1$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen1$log2FoldChange, sigtabgen1$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen1$Genus = factor(as.character(sigtabgen1$Genus), levels = names(x))


T1Gen<-ggplot(sigtabgen1, aes(x = Genus, y = log2FoldChange, color = Family)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#701C00","pink","#62a6b3",
        "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC","#00aba2")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T1Gen





#Timepoint 2 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen2 = subset(sigtab2, !is.na(Genus))
# Family order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Family = factor(as.character(sigtabgen2$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Genus = factor(as.character(sigtabgen2$Genus), levels = names(x))


T2Gen<-ggplot(sigtabgen2, aes(x = Genus, y = log2FoldChange, color = Family)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#701C00","pink","#62a6b3",
          "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC",
          "#00aba2","#18330d","#1102e6","#876641","#c7eaeb","#ffd500",
          "#4b6080","#c710c1","#24702d")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T2Gen


#Timepoint 3 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen3= subset(sigtab3, !is.na(Genus))
# Family order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Family = factor(as.character(sigtabgen3$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Genus = factor(as.character(sigtabgen3$Genus), levels = names(x))


T3Gen<-ggplot(sigtabgen3, aes(x = Genus, y = log2FoldChange, color = Family)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#701C00","pink","#62a6b3",
          "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC",
           "#00aba2","#18330d","#1102e6","#876641","#c7eaeb","#ffd500","#4b6080")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T3Gen




#Timepoint 4 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen4= subset(sigtab4, !is.na(Genus))
# Family order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Family = factor(as.character(sigtabgen4$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Genus = factor(as.character(sigtabgen4$Genus), levels = names(x))


T4Gen<-ggplot(sigtabgen4, aes(x = Genus, y = log2FoldChange, color = Family)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=c("#701C00","pink","#62a6b3",
    "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC",
     "#00aba2","#18330d","#1102e6","#876641","#c7eaeb",
    "#ffd500","#4b6080")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T4Gen


#Timepoint 5 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen5= subset(sigtab5, !is.na(Genus))
# Family order
x = tapply(sigtabgen5$log2FoldChange, sigtabgen5$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen5$Family = factor(as.character(sigtabgen5$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen5$log2FoldChange, sigtabgen5$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen5$Genus = factor(as.character(sigtabgen5$Genus), levels = names(x))


T5Gen<-ggplot(sigtabgen5, aes(x = Genus, y = log2FoldChange, color = Family)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#701C00","pink","#62a6b3",
       "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC",
       "#00aba2","#18330d","#1102e6","#876641","#c7eaeb","#ffd500","#4b6080",
       "#7355a3","#a8a476","#6e500a","#b81286","#ededca","red")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T5Gen




#Timepoint 6 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen6= subset(sigtab6, !is.na(Genus))
# Family order
x = tapply(sigtabgen6$log2FoldChange, sigtabgen6$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen6$Family = factor(as.character(sigtabgen6$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen6$log2FoldChange, sigtabgen6$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen6$Genus = factor(as.character(sigtabgen6$Genus), levels = names(x))


T6Gen<-ggplot(sigtabgen6, aes(x = Genus, y = log2FoldChange, color = Family)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#701C00","pink","#62a6b3",
   "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC",
   "#00aba2","#18330d","#1102e6","#876641","#c7eaeb","#ffd500",
   "#4b6080","#03fce3")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T6Gen



#Timepoint 7 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen7= subset(sigtab7, !is.na(Genus))
# Family order
x = tapply(sigtabgen7$log2FoldChange, sigtabgen7$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen7$Family = factor(as.character(sigtabgen7$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen7$log2FoldChange, sigtabgen7$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen7$Genus = factor(as.character(sigtabgen7$Genus), levels = names(x))


T7Gen<-ggplot(sigtabgen7, aes(x = Genus, y = log2FoldChange, color = Family)) +  
    geom_point(size=4)+ scale_color_manual(values=c("#701C00","pink","#62a6b3",
      "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC",
      "#00aba2","#18330d","#1102e6","#876641","#c7eaeb","#ffd500","#4b6080",
      "#7355a3","#a8a476","#6e500a","#b81286","#ededca","#03fce3","red","green")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T7Gen




#Timepoint 8 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen8= subset(sigtab8, !is.na(Genus))
# Family order
x = tapply(sigtabgen8$log2FoldChange, sigtabgen8$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen8$Family = factor(as.character(sigtabgen8$Family), levels = names(x))
# Genus order
x = tapply(sigtabgen8$log2FoldChange, sigtabgen8$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen8$Genus = factor(as.character(sigtabgen8$Genus), levels = names(x))


T8Gen<-ggplot(sigtabgen8, aes(x = Genus, y = log2FoldChange, color = Family)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#701C00","pink","#62a6b3",
     "#0c6e05","#e6da02","#D59D3A","#E0D2B9","#c90e24","#84A9BC",
     "#00aba2","#18330d","#1102e6","#876641","#c7eaeb","#ffd500","#4b6080",
     "#7355a3","#a8a476","#6e500a","#b81286","#ededca","#03fce3",
     "#4107e0","red","#0f0f06")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5));T8Gen








#
#
#****************************************************************************************************-----
#--- CREATE GRAPHS Phylum TO SPECIES ---------------------------------------------------------------------------------------
#****************************************************************************************************-----

#Load the files with the manually blasted species data
sigtab1<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab1.csv")
sigtab2<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab2.csv")
sigtab3<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab3.csv")
sigtab4<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab4.csv")
sigtab5<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab5.csv")
sigtab6<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab6.csv")
sigtab7<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab7.csv")
sigtab8<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab8.csv")
sigtab9<-read.csv("1-Analysis/Deseq/Fungi/Tables/sigtab9.csv")



#Timepoint 1 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen1 = subset(sigtab1, !is.na(Species))
# Phylum order
x = tapply(sigtabgen1$log2FoldChange, sigtabgen1$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen1$Phylum = factor(as.character(sigtabgen1$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen1$log2FoldChange, sigtabgen1$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen1$Species = factor(as.character(sigtabgen1$Species), levels = names(x))


T1sp<-ggplot(sigtabgen1, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1, face ="italic" ));T1sp





#Timepoint 2 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen2 = subset(sigtab2, !is.na(Species))
# Phylum order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Phylum = factor(as.character(sigtabgen2$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Species = factor(as.character(sigtabgen2$Species), levels = names(x))


T2sp<-ggplot(sigtabgen2, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,face ="italic"));T2sp


#Timepoint 3 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen3= subset(sigtab3, !is.na(Species))
# Phylum order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Phylum = factor(as.character(sigtabgen3$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Species = factor(as.character(sigtabgen3$Species), levels = names(x))


T3sp<-ggplot(sigtabgen3, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,face ="italic"));T3sp




#Timepoint 4 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen4= subset(sigtab4, !is.na(Species))
# Phylum order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Phylum = factor(as.character(sigtabgen4$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Species = factor(as.character(sigtabgen4$Species), levels = names(x))


T4sp<-ggplot(sigtabgen4, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,face ="italic"));T4sp


#Timepoint 5 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen5= subset(sigtab5, !is.na(Species))
# Phylum order
x = tapply(sigtabgen5$log2FoldChange, sigtabgen5$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen5$Phylum = factor(as.character(sigtabgen5$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen5$log2FoldChange, sigtabgen5$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen5$Species = factor(as.character(sigtabgen5$Species), levels = names(x))


T5sp<-ggplot(sigtabgen5, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,face ="italic"));T5sp




#Timepoint 6 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen6= subset(sigtab6, !is.na(Species))
# Phylum order
x = tapply(sigtabgen6$log2FoldChange, sigtabgen6$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen6$Phylum = factor(as.character(sigtabgen6$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen6$log2FoldChange, sigtabgen6$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen6$Species = factor(as.character(sigtabgen6$Species), levels = names(x))


T6sp<-ggplot(sigtabgen6, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,face ="italic"));T6sp



#Timepoint 7 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen7= subset(sigtab7, !is.na(Species))
# Phylum order
x = tapply(sigtabgen7$log2FoldChange, sigtabgen7$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen7$Phylum = factor(as.character(sigtabgen7$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen7$log2FoldChange, sigtabgen7$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen7$Species = factor(as.character(sigtabgen7$Species), levels = names(x))


T7sp<-ggplot(sigtabgen7, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1));T7sp




#Timepoint 8 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen8= subset(sigtab8, !is.na(Species))
# Phylum order
x = tapply(sigtabgen8$log2FoldChange, sigtabgen8$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen8$Phylum = factor(as.character(sigtabgen8$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen8$log2FoldChange, sigtabgen8$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen8$Species = factor(as.character(sigtabgen8$Species), levels = names(x))


T8sp<-ggplot(sigtabgen8, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
 theme(axis.text.x = element_text(angle = 45,hjust = 1,face ="italic"));T8sp



#Timepoint 9 ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
sigtabgen9= subset(sigtab9, !is.na(Species))
# Phylum order
x = tapply(sigtabgen9$log2FoldChange, sigtabgen9$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen9$Phylum = factor(as.character(sigtabgen9$Phylum), levels = names(x))
# Species order
x = tapply(sigtabgen9$log2FoldChange, sigtabgen9$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen9$Species = factor(as.character(sigtabgen9$Species), levels = names(x))


T9sp<-ggplot(sigtabgen9, aes(x = Species, y = log2FoldChange, color = Phylum)) +  
      geom_point(size=4)+ scale_color_manual(values=c("#b0022a","#7aa0c2","#cca002")) +
 theme(axis.text.x = element_text(angle = 45,hjust = 1,face ="italic"));T9sp







#
#
#***************************************************************************************-----
#--------------------- EXPORT GRAPHS --------------------------------------------------------
#***************************************************************************************-----
dir.create(file.path("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp"), recursive = TRUE)

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T1-Treatment-PhylumSp.pdf", height=8, width=10)
T1sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T2-Treatment-PhylumSp.pdf", height=8, width=10)
T2sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T3-Treatment-PhylumSp.pdf", height=8, width=10)
T3sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T4-Treatment-PhylumSp.pdf", height=8, width=10)
T4sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T5-Treatment-PhylumSp.pdf", height=8, width=10)
T5sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T6-Treatment-PhylumSp.pdf", height=8, width=10)
T6sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T7-Treatment-PhylumSp.pdf", height=8, width=10)
T7sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T8-Treatment-PhylumSp.pdf", height=8, width=10)
T8sp
dev.off()

pdf("1-Analysis/Deseq/Fungi/Graphs/Timepoints/PhylumSpp/T9-Treatment-PhylumSp.pdf", height=8, width=10)
T9sp
dev.off()













