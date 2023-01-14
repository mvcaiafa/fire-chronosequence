#Reset R's Brain
rm(list=ls())

#Set working directory............................................................
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence/abundance_barplot")

#Load Packages-------------------------------------------------------------------------------
library(qiime2R)
library(DESeq2)#log2fold analysis
library(ggplot2)#Plotting
library(phyloseq)# to load in qiime files and work w them based on .tza extension
library(tidyverse)
library(tidyr)


#LOAD QIIME rarified table -------------------------------------------------------------------------------------------
Metadata_B<-read_tsv("metadata_bacteria-merged-NC.tsv")#Rare metadata
Bac_Table<-read_qza("B_rarefied_table.qza", tmp = "C:/tmp")#Rarefied table
Bac_Tree<-read_qza("Bacteria-rooted-tree-NC.qza",tmp = "C:/tmp" )
Bac_Taxonomy<-read_qza("Bacteria-Taxonomy132.qza",tmp = "C:/tmp")
Bac_Taxtable<-Bac_Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                                                             c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

#CREATE PHYLOSEQ ARTIFACT
Bac_physeq<-phyloseq(
  otu_table(Bac_Table$data, taxa_are_rows = TRUE), 
  phy_tree(Bac_Tree$data), 
  tax_table(as.data.frame(Bac_Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata_B%>% as.data.frame() %>% column_to_rownames("sample-id"))
)

#----
#----
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(Bac_physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number..............................................................
colnames(tax_table(Bac_physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(Bac_physeq)[, "Kingdom"] <- gsub("D_0__", "", tax_table(Bac_physeq)[, "Kingdom"])
tax_table(Bac_physeq)[, "Phylum"] <- gsub("D_1__", "", tax_table(Bac_physeq)[, "Phylum"])
tax_table(Bac_physeq)[, "Class"] <- gsub("D_2__", "", tax_table(Bac_physeq)[, "Class"])
tax_table(Bac_physeq)[, "Order"] <- gsub("D_3__", "", tax_table(Bac_physeq)[, "Order"])
tax_table(Bac_physeq)[, "Family"] <- gsub("D_4__", "", tax_table(Bac_physeq)[, "Family"])
tax_table(Bac_physeq)[, "Genus"] <- gsub("D_5__", "", tax_table(Bac_physeq)[, "Genus"])
tax_table(Bac_physeq)[, "Species"] <- gsub("D_6__", "", tax_table(Bac_physeq)[, "Species"])

sample_data(Bac_physeq)$severity<-factor(sample_data(physeq)$severity) # make as factor

#SUBSAMPLE DATA TO LOOK ONLY AT SOIL SAMPLES ----------------------------------------------------------------------
Bac_physeqSoil<-subset_samples(Bac_physeq,severity!="low");Bac_physeqSoil#29621x187
sample_names(Bac_physeqSoil); sample_data(Bac_physeqSoil)$severity

#Look at the data to make sure there are no "none" or NA = remove
head(sample_data(Bac_physeqSoil)$severity, 50)
head(sample_data(Bac_physeqSoil)$fire, 50)
head(sample_data(Bac_physeqSoil)$depth, 50)

#Check the structure of the data------------------------------------------------------
str(sample_data(Bac_physeqSoil)$severity)#chr
str(sample_data(Bac_physeqSoil)$fire)#num
str(sample_data(Bac_physeqSoil)$depth)#chr

#
#
#*********************************************************************************************-----
#--- SUBSET DATA BY Site FOR ANALYSIS --------------------------------------------------------
#*********************************************************************************************-----

#Subset samples per Fire--------------
Bac_physeqSoilCP<-subset_samples(Bac_physeqSoil,fire=="church park")
sample_names(Bac_physeqSoilCP)

Bac_physeqSoilBV<-subset_samples(Bac_physeqSoil,fire=="beaver creek")
sample_names(Bac_physeqSoilBV)

Bac_physeqSoilBC<-subset_samples(Bac_physeqSoil,fire=="badger creek")
sample_names(Bac_physeqSoilBC)

Bac_physeqSoilR<-subset_samples(Bac_physeqSoil,fire=="ryan")
sample_names(Bac_physeqSoilR)

Bac_physeqSoilM<-subset_samples(Bac_physeqSoil,fire=="mullen")
sample_names(Bac_physeqSoilM)

#
#
#****************************************************************************************************-----
#--- RUN DESEQ ANALYSIS TREATMENT PER Site----------------------------------------------------
#****************************************************************************************************-----
#

#Church Park **********************************************************************----
head(sample_data(Bac_physeqSoilCP)$severity, 35)

#Ran in case there are any n/a in data
Bac_physeqSoilCP<-subset_samples(Bac_physeqSoilCP, severity!= "None")

#Import phyloseq data to create a deseq object.....................
Bac_dsCP<-phyloseq_to_deseq2(Bac_physeqSoilCP, ~ severity)

# calculate geometric means prior to estimate size factors.........
gm_mean1 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
Bac_geoMeansCP = apply(counts(Bac_dsCP), 1, gm_mean1)
Bac_dsCP<-estimateSizeFactors(Bac_dsCP, geoMeans = Bac_geoMeansCP)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#Bac_fdsCP<-rowSums(counts(Bac_dsCP, normalize = TRUE) >= 200)>=2
#Bac_dsCP <- Bac_dsCP[Bac_fdsCP,]

#Set reference level...............................................
Bac_dsCP$severity<-relevel(Bac_dsCP$severity,ref="control")
levels(Bac_dsCP$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ...............................
# (a.k.a. Gamma-Poisson) distribution
Bac_dsCP<-DESeq(Bac_dsCP, fitType="local")

#Beaver Creek**********************************************************************----
head(sample_data(Bac_physeqSoilBV)$severity, 35)

#Ran in case there are any n/a in data
Bac_physeqSoilBV<-subset_samples(Bac_physeqSoilBV, severity!= "None")

#Import phyloseq data to create a deseq object...................
Bac_dsBV<-phyloseq_to_deseq2(Bac_physeqSoilBV, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean2 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
Bac_geoMeansBV = apply(counts(Bac_dsBV), 1, gm_mean2)
Bac_dsBV<-estimateSizeFactors(Bac_dsBV, geoMeans = Bac_geoMeansBV)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#Bac_fdsBV<-rowSums(counts(Bac_dsBV, normalize = TRUE) >= 200)>=2
#Bac_dsBV <- Bac_dsBV[Bac_fdsBV,]

#Set reference level.............................................
Bac_dsBV$severity<-relevel(Bac_dsBV$severity,ref="control")
levels(Bac_dsBV$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
Bac_dsBV<-DESeq(Bac_dsBV, fitType="local")


#Badger Creek **********************************************************************----

head(sample_data(Bac_physeqSoilBC)$severity, 35)

#Ran in case there are any n/a in data
Bac_physeqSoilBC<-subset_samples(Bac_physeqSoilBC, severity!= "None")

#Import phyloseq data to create a deseq object...................
Bac_dsBC<-phyloseq_to_deseq2(Bac_physeqSoilBC, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean3 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
Bac_geoMeansBC = apply(counts(Bac_dsBC), 1, gm_mean3)
Bac_dsBC<-estimateSizeFactors(Bac_dsBC, geoMeans = Bac_geoMeansBC)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#Bac_fdsBC<-rowSums(counts(Bac_dsBC, normalize = TRUE) >= 200)>=2
#Bac_dsBC <- Bac_dsBC[Bac_fdsBC,]

#Set reference level.............................................
Bac_dsBC$severity<-relevel(Bac_dsBC$severity,ref="control")
levels(Bac_dsBC$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
Bac_dsBC<-DESeq(Bac_dsBC, fitType="local")


#Ryan **********************************************************************----
head(sample_data(Bac_physeqSoilR)$severity, 35)

#Ran in case there are any n/a in data
Bac_physeqSoilR<-subset_samples(Bac_physeqSoilR, severity!= "None")

#Import phyloseq data to create a deseq object...................
Bac_dsR<-phyloseq_to_deseq2(Bac_physeqSoilR, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean4 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
Bac_geoMeansR = apply(counts(Bac_dsR), 1, gm_mean4)
Bac_dsR<-estimateSizeFactors(Bac_dsR, geoMeans = Bac_geoMeansR)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#Bac_fdsR<-rowSums(counts(Bac_dsR, normalize = TRUE) >= 200)>=2
#Bac_dsR <- Bac_dsR[Bac_fdsR,]

#Set reference level.............................................
Bac_dsR$severity<-relevel(Bac_dsR$severity,ref="control")
levels(Bac_dsR$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
Bac_dsR<-DESeq(Bac_dsR, fitType="local")


#Mullen**********************************************************************----
head(sample_data(Bac_physeqSoilM)$severity, 35)

#Ran in case there are any n/a in data
Bac_physeqSoilM<-subset_samples(Bac_physeqSoilM, severity!= "None")

#Import phyloseq data to create a deseq object...................
Bac_dsM<-phyloseq_to_deseq2(Bac_physeqSoilM, ~ severity)

# calculate geometric means prior to estimate size factors.......
gm_mean5 = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
Bac_geoMeansM = apply(counts(Bac_dsM), 1, gm_mean5)
Bac_dsM<-estimateSizeFactors(Bac_dsM, geoMeans = Bac_geoMeansM)

#Filter out ASVs that have less than 200 reads and are no present in at least 2 samples
#Bac_fdsM<-rowSums(counts(Bac_dsM, normalize = TRUE) >= 200)>=2
#Bac_dsM <- Bac_dsM[Bac_fdsM,]

#Set reference level.............................................
Bac_dsM$severity<-relevel(Bac_dsM$severity,ref="control")
levels(Bac_dsM$severity)

# RUN DESEQ-BASED ON NEG BINOMIAL  ..............................
# (a.k.a. Gamma-Poisson) distribution
Bac_dsM<-DESeq(Bac_dsM, fitType="local")


#
#****************************************************************************************************-----
#--- RUN RESULT ------------------------------------------------------------------------------------------
#****************************************************************************************************-----


#Church Park **********************************************************************----
Bac_resCP<-results(Bac_dsCP)
Bac_resCP<-Bac_resCP[order(Bac_resCP$padj, na.last=NA), ]
alpha<-0.01
Bac_sigtabCP<-Bac_resCP[(Bac_resCP$padj < alpha), ];Bac_sigtabCP
Bac_sigtabCP<-cbind(as(Bac_sigtabCP, "data.frame"), 
                as(tax_table(Bac_physeqSoilCP)[rownames(Bac_sigtabCP), ], "matrix"));Bac_sigtabCP


#Beaver Creek *********************************************************************----
Bac_resBV<-results(Bac_dsBV)
Bac_resBV<-Bac_resBV[order(Bac_resBV$padj, na.last=NA), ]
alpha<-0.01
Bac_sigtabBV<-Bac_resBV[(Bac_resBV$padj < alpha), ];Bac_sigtabBV
Bac_sigtabBV<-cbind(as(Bac_sigtabBV, "data.frame"), 
                as(tax_table(Bac_physeqSoilBV)[rownames(Bac_sigtabBV), ], "matrix"));Bac_sigtabBV


#Badger Creek**********************************************************************----
Bac_resBC<-results(Bac_dsBC)
Bac_resBC<-Bac_resBC[order(Bac_resBC$padj, na.last=NA), ]
alpha<-0.01
Bac_sigtabBC<-Bac_resBC[(Bac_resBC$padj < alpha), ];Bac_sigtabBC
Bac_sigtabBC<-cbind(as(Bac_sigtabBC, "data.frame"), 
                as(tax_table(Bac_physeqSoilBC)[rownames(Bac_sigtabBC), ], "matrix"));Bac_sigtabBC


#Ryan**********************************************************************----
Bac_resR<-results(Bac_dsR)
Bac_resR<-Bac_resR[order(Bac_resR$padj, na.last=NA), ]
alpha<-0.01
Bac_sigtabR<-Bac_resR[(Bac_resR$padj < alpha), ];Bac_sigtabR
Bac_sigtabR<-cbind(as(Bac_sigtabR, "data.frame"), 
               as(tax_table(Bac_physeqSoilR)[rownames(Bac_sigtabR), ], "matrix"));Bac_sigtabR



#Mullen**********************************************************************----
Bac_resM<-results(Bac_dsM)
Bac_resM<-Bac_resM[order(Bac_resM$padj, na.last=NA), ]
alpha<-0.01
Bac_sigtabM<-Bac_resM[(Bac_resM$padj < alpha), ];Bac_sigtabM
Bac_sigtabM<-cbind(as(Bac_sigtabM, "data.frame"), 
               as(tax_table(Bac_physeqSoilM)[rownames(Bac_sigtabM), ], "matrix"));Bac_sigtabM

#****************************************************************************************************-----
#--- EXPORT FILES ---------------------------------------------------------------------------------------
#****************************************************************************************************-----

write.csv(Bac_sigtabCP,"Bac_sigtabCP.csv")
write.csv(Bac_sigtabBV,"Bac_sigtabBV.csv")
write.csv(Bac_sigtabBC,"Bac_sigtabBC.csv")
write.csv(Bac_sigtabR,"Bac_sigtabR.csv")
write.csv(Bac_sigtabM,"Bac_sigtabM.csv")


###Make Figures

#Load dataset .....................................................................
Bac_sigtabCP<-read.csv("Bac_sigtabCP.csv", na.strings = "NA")
Bac_sigtabBV<-read.csv("Bac_sigtabBV.csv", na.strings = "NA")
Bac_sigtabBC<-read.csv("Bac_sigtabBC.csv", na.strings = "NA")
Bac_sigtabR<-read.csv("Bac_sigtabR.csv", na.strings = "NA")
Bac_sigtabM<-read.csv("Bac_sigtabM.csv", na.strings = "NA")


Bac_sigtabCP<-Bac_sigtabCP[(Bac_sigtabCP$Genus!="unidentified"&Bac_sigtabCP$Genus!="uncultured"&Bac_sigtabCP$Genus!="uncultured bacterium"),]
Bac_sigtabBV<-Bac_sigtabBV[(Bac_sigtabBV$Genus!="unidentified"&Bac_sigtabBV$Genus!="uncultured"&Bac_sigtabBV$Genus!="uncultured bacterium"),]
Bac_sigtabBC<-Bac_sigtabBC[(Bac_sigtabBC$Genus!="unidentified"&Bac_sigtabBC$Genus!="uncultured"&Bac_sigtabBC$Genus!="uncultured bacterium"),]
Bac_sigtabR<-Bac_sigtabR[(Bac_sigtabR$Genus!="unidentified"&Bac_sigtabR$Genus!="uncultured"&Bac_sigtabR$Genus!="uncultured bacterium"),]
Bac_sigtabM<-Bac_sigtabM[(Bac_sigtabM$Genus!="unidentified"&Bac_sigtabM$Genus!="uncultured"&Bac_sigtabM$Genus!="uncultured bacterium"),]

#****************************************************************************************************-----
#--- CREATE GRAPHS PHYLUM TO Genus ---------------------------------------------------------------------------------------
#****************************************************************************************************-----
Bac_colores<-c("Acidobacteria"="#01665E","Actinobacteria"="#7F3B08","Bacteroidetes"="#A50026",
               "Chloroflexi"="#35978F","Firmicutes"="#8073AC","Gemmatimonadetes"="#E7D4E8","Planctomycetes"="#5AAE61",
               "Proteobacteria"="#FDB863","Verrucomicrobia"="#FEE0B6")

#Church Park ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
Bac_sigtabgenCP = subset(Bac_sigtabCP, !is.na(Genus))
# Phylum order
x = tapply(Bac_sigtabgenCP$log2FoldChange, Bac_sigtabgenCP$Phylum, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenCP$Phylum = factor(as.character(Bac_sigtabgenCP$Phylum), levels = names(x))
# Genus order
x = tapply(Bac_sigtabgenCP$log2FoldChange, Bac_sigtabgenCP$Genus, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenCP$Genus = factor(as.character(Bac_sigtabgenCP$Genus), levels = names(x))


Bac_CP<-ggplot(Bac_sigtabgenCP, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=Bac_colores) +
  theme(axis.text.x = element_text(angle = 45, hjust = .98, vjust=1,size=8, colour = "black"),
        axis.text.y = element_text(colour = "black", size=8, face = "italic"),
        axis.title.x = element_blank())

Bac_CP

pdf("DESEQ_CPBac.pdf", height=4, width=8)
Bac_CP
dev.off()


#Beaver Creek ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
Bac_sigtabgenBV = subset(Bac_sigtabBV, !is.na(Genus))
# Phylum order
x = tapply(Bac_sigtabgenBV$log2FoldChange, Bac_sigtabgenBV$Phylum, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenBV$Phylum = factor(as.character(Bac_sigtabgenBV$Phylum), levels = names(x))
# Genus order
x = tapply(Bac_sigtabgenBV$log2FoldChange, Bac_sigtabgenBV$Genus, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenBV$Genus = factor(as.character(Bac_sigtabgenBV$Genus), levels = names(x))

Bac_BV<-ggplot(Bac_sigtabgenBV, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=Bac_colores) +
  theme(axis.text.x = element_text(angle = 45, hjust = .98, vjust=1,size=8, colour = "black"),
        axis.text.y = element_text(colour = "black", size=8, face = "italic"),
        axis.title.x = element_blank())

Bac_BV

pdf("DESEQ_BVbac.pdf", height=4, width=8)
Bac_BV
dev.off()

#Badger Creek ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
Bac_sigtabgenBC = subset(Bac_sigtabBC, !is.na(Genus))
# Phylum order
x = tapply(Bac_sigtabgenBC$log2FoldChange, Bac_sigtabgenBC$Phylum, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenBC$Phylum = factor(as.character(Bac_sigtabgenBC$Phylum), levels = names(x))
# Genus order
x = tapply(Bac_sigtabgenBC$log2FoldChange, Bac_sigtabgenBC$Genus, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenBC$Genus = factor(as.character(Bac_sigtabgenBC$Genus), levels = names(x))

Bac_BC<-ggplot(Bac_sigtabgenBC, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=Bac_colores) +
  theme(axis.text.x = element_text(angle = 45, hjust = .98, vjust=1,size=8, colour = "black"),
        axis.text.y = element_text(colour = "black", size=8, face = "italic"),
        axis.title.x = element_blank())

Bac_BC


pdf("DESEQ_BCbac.pdf", height=4, width=8)
Bac_BC
dev.off()

#Ryan ........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
Bac_sigtabgenR = subset(Bac_sigtabR, !is.na(Genus))
# Phylum order
x = tapply(Bac_sigtabgenR$log2FoldChange, Bac_sigtabgenR$Phylum, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenR$Phylum = factor(as.character(Bac_sigtabgenR$Phylum), levels = names(x))
# Genus order
x = tapply(Bac_sigtabgenR$log2FoldChange, Bac_sigtabgenR$Genus, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenR$Genus = factor(as.character(Bac_sigtabgenR$Genus), levels = names(x))

Bac_R<-ggplot(Bac_sigtabgenR, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=Bac_colores) +
  theme(axis.text.x = element_text(angle = 45, hjust = .98, vjust=1,size=8, colour = "black"),
        axis.text.y = element_text(colour = "black", size=8, face = "italic"),
        axis.title.x = element_blank())

Bac_R


pdf("DESEQ_Rbac.pdf", height=4, width=8)
Bac_R
dev.off()


#Mullen........................................................................
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palette, ...)
}
Bac_sigtabgenM = subset(Bac_sigtabM, !is.na(Genus))
# Phylum order
x = tapply(Bac_sigtabgenM$log2FoldChange, Bac_sigtabgenM$Phylum, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenM$Phylum = factor(as.character(Bac_sigtabgenM$Phylum), levels = names(x))
# Genus order
x = tapply(Bac_sigtabgenM$log2FoldChange, Bac_sigtabgenM$Genus, function(x) max(x))
x = sort(x, TRUE)
Bac_sigtabgenM$Genus = factor(as.character(Bac_sigtabgenM$Genus), levels = names(x))

Bac_M<-ggplot(Bac_sigtabgenM, aes(x = Genus, y = log2FoldChange, color = Phylum)) +  
  geom_point(size=4)+ 
  scale_color_manual(values=Bac_colores) +
  theme(axis.text.x = element_text(angle = 45, hjust = .98, vjust=1,size=8, colour = "black"),
        axis.text.y = element_text(colour = "black", size=8, face = "italic"),
        axis.title.x = element_blank())

Bac_M

pdf("DESEQ_Mbac.pdf", height=4, width=8)
Bac_M
dev.off()

Bac_deseqfigure <- ggarrange(Bac_CP,Bac_BV,Bac_BC,Bac_R,Bac_M, labels = c("A", "B","C","D","E"), ncol = 2, nrow = 3, common.legend=TRUE)
Bac_deseqfigure

pdf("DESEQ_BACTERIA.pdf", height=10, width=10)
Bac_deseqfigure
dev.off()

