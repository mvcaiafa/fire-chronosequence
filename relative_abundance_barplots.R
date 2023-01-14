#Set working directory............................................................
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence/abundance_barplot")

##INstall packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")


#Load librarires..................................................................
library(phyloseq)#load qiime files and work w them based on .tza extension
library(qiime2R)
library(ape)#to build tree but I imported from QIIME
library(tidyverse) #required to load tax table
library(ggplot2)
library(ggpubr) 
library(ochRe)
library(MicEco)#tp create venndiagram


###########Fungi#####################
#####################################

#LOAD QIIME DATA-NON-NORMALIZED -------------------------------------------------------------------------------------------
Metadata<-read_tsv("metadata_fungi-merged_NC-2.tsv")#Rare metadata
Metadata<-read.csv("metadata_fungi-merged_NC-2.csv")

Table<-read_qza("rarefied_table.qza", tmp = "C:/tmp")#Rarefied table
Tree<-read_qza("Fungal-rooted-tree-NC.qza",tmp = "C:/tmp" )
Taxonomy<-read_qza("Fungal-Taxonomy.qza",tmp = "C:/tmp")
Taxtable<-Taxonomy$data %>% as_tibble() %>% separate(Taxon, sep = ";",
                                                     c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))


#CREATE PHYLOSEQ ARTIFACT
physeq<-phyloseq(
  otu_table(Table$data, taxa_are_rows = TRUE), 
  phy_tree(Tree$data), 
  tax_table(as.data.frame(Taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), 
  sample_data(Metadata%>% as.data.frame() %>% column_to_rownames("sample.id"))
)

#----
#----
#********************************************************************************************************************----
#------------------------------------     QUALITY CONTROL        ---------------------------------------------
#********************************************************************************************************************----
rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number..............................................................
colnames(tax_table(physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(physeq)[, "Kingdom"] <- gsub("k__", "", tax_table(physeq)[, "Kingdom"])
tax_table(physeq)[, "Phylum"] <- gsub("p__", "", tax_table(physeq)[, "Phylum"])
tax_table(physeq)[, "Class"] <- gsub("c__", "", tax_table(physeq)[, "Class"])
tax_table(physeq)[, "Order"] <- gsub("o__", "", tax_table(physeq)[, "Order"])
tax_table(physeq)[, "Family"] <- gsub("f__", "", tax_table(physeq)[, "Family"])
tax_table(physeq)[, "Genus"] <- gsub("g__", "", tax_table(physeq)[, "Genus"])
tax_table(physeq)[, "Species"] <- gsub("s__", "", tax_table(physeq)[, "Species"])

#--SUBSET DATA IN DIFFERENT SEVERITY SAMPLES ...........................................................

df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)
str(sample_data(physeq))

# * * Subset All to maintain only th.................................
physeqAllControl<-subset_samples(physeq,severity=="control");physeqAllControl
physeqAllHigh<-subset_samples(physeq,severity=="high");physeqAllHigh
physeqAllLow<-subset_samples(physeq,severity=="low");physeqAllLow


#----
#----
#**************************************************************************************************************-----
# RELATIVE ABUNDANCE SOIL SAMPLES ONLY--------------------------------------------------------------------------
#**************************************************************************************************************-----
#Phylum 
Phylum <- tax_glom(physeq, taxrank = 'Phylum') # agglomerate taxa
(PhySev = merge_samples(Phylum, "treatment")) # merge samples on sample variable of interest
RelPhySev <- transform_sample_counts(PhySev, function(x) x/sum(x)) #get abundance in %
RelPhySev1 <- psmelt(RelPhySev) # create dataframe from phyloseq object
RelPhySev1$Phylum <- as.character(RelPhySev1$Phylum) #convert to character
RelPhySev1$Phylum[RelPhySev1$Abundance < 0.05] 

AbunSev<-ggplot(data=RelPhySev1, aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(aes(), stat="identity", position="stack") + 
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 38, colour = "black"), 
        axis.title.x = element_text(size = 36, colour = "black"), 
        axis.text.y = element_text(size=30, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=30, colour = "black"),
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=28),
        legend.title=element_text(size=30, face="bold"))+ 
  xlab("Severity")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 4))

AbunSev


# * * * LOOK AT SOIL SAMPLES ONLY  -----------------------------------------------------------------------
#** Treatment analysis.................................................................................
Genus <- tax_glom(physeq, taxrank = 'Genus') # agglomerate taxa
(GenSev = merge_samples(Genus, "treatment")) # merge samples on sample variable of interest
RelGenSev <- transform_sample_counts(GenSev, function(x) x/sum(x)) #get abundance in %
RelGenSev1 <- psmelt(RelGenSev) # create dataframe from phyloseq object
RelGenSev1$Genus <- as.character(RelGenSev1$Genus) #convert to character
RelGenSev1$Genus[RelGenSev1$Abundance < 0.05] <- "< 5% abund." #rename genera with < 3% abundance


library(RColorBrewer)
n <- 24
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div'& brewer.pal.info$colorblind ,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
f=sample(col_vector, n)
pie(rep(1,n), col=f)

f

FunTrt<-c("< 5% abund."="#dee0df", "Calyptrozyma"="#ABD9E9","Collophora"="#8073AC","Coniochaeta" =  "#5AAE61" ,
          "Cortinarius"="#F6E8C3","Gymnostellatospora"="#003C30", "Hyphodiscus"="#80CDC1","Inocybe"="#A6DBA0",
          "Leohumicola"="#74ADD1","Leucosporidium"="#7FBC41","Mallocybe"="#B35806","Mycena"="#542788",
          "Naganishia"="#2D004B", "Oidiodendron"="#762A83","Paratritirachium"= "#C7EAE5", "Penicillium"="#8C510A",
          "Phialocephala"= "#8E0152","Pholiota"="#D9F0D3","Piloderma"="#053061", "Serpula"="#F46D43",
          "Tetracladium"="#F4A582","Tricholoma"="#D6604D","unidentified"="#4575B4","Venturia"=  "#00441B")
          
           
       "#1B7837"
 #TREATMENT...............................................................................
AbunSev<-ggplot(data=RelGenSev1, aes(x=Sample, y=Abundance, fill=Genus))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = FunTrt)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 38, colour = "black"), 
        axis.title.x = element_text(size = 36, colour = "black"), 
        axis.text.y = element_text(size=30, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=30, colour = "black"),
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=28),
        legend.title=element_text(size=30, face="bold"))+ 
  xlab("Severity")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 4))

AbunSev


#####################################
###########Bacteria##################
#####################################

#LOAD QIIME DATA-NON-NORMALIZED -------------------------------------------------------------------------------------------
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
rank_names(physeq)# Look at rank names

#Quality control: Remove the g__ from each rank number..............................................................
colnames(tax_table(Bac_physeq))= c("Kingdom","Phylum","Class","Order","Family","Genus","Species", "Confidence")
tax_table(Bac_physeq)[, "Kingdom"] <- gsub("D_0__", "", tax_table(Bac_physeq)[, "Kingdom"])
tax_table(Bac_physeq)[, "Phylum"] <- gsub("D_1__", "", tax_table(Bac_physeq)[, "Phylum"])
tax_table(Bac_physeq)[, "Class"] <- gsub("D_2__", "", tax_table(Bac_physeq)[, "Class"])
tax_table(Bac_physeq)[, "Order"] <- gsub("D_3__", "", tax_table(Bac_physeq)[, "Order"])
tax_table(Bac_physeq)[, "Family"] <- gsub("D_4__", "", tax_table(Bac_physeq)[, "Family"])
tax_table(Bac_physeq)[, "Genus"] <- gsub("D_5__", "", tax_table(Bac_physeq)[, "Genus"])
tax_table(Bac_physeq)[, "Species"] <- gsub("D_6__", "", tax_table(Bac_physeq)[, "Species"])


#--SUBSET DATA IN DIFFERENT SEVERITY SAMPLES ...........................................................

df <- as.data.frame(lapply(sample_data(Bac_physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(Bac_physeq)
sample_data(Bac_physeq) <- sample_data(df)
str(sample_data(Bac_physeq))

#----
#----
#**************************************************************************************************************-----
# RELATIVE ABUNDANCE SOIL SAMPLES -------------------------------------------------------------------------
#**************************************************************************************************************-----

# * * * LOOK AT SOIL SAMPLES -----------------------------------------------------------------------
#** Treatment analysis.................................................................................
Phylum <- tax_glom(Bac_physeq, taxrank = 'Phylum') # agglomerate taxa
(PhySev = merge_samples(Phylum, "treatment")) # merge samples on sample variable of interest
RelPhySev <- transform_sample_counts(PhySev, function(x) x/sum(x)) #get abundance in %
RelPhySev1 <- psmelt(RelPhySev) # create dataframe from phyloseq object
RelPhySev1$Phylum <- as.character(RelPhySev1$Phylum) #convert to character
RelPhySev1$Phylum[RelPhySev1$Abundance < 0.02] <- "< 2% abund." #rename genera with < 2% abundance

library(RColorBrewer)
n <- 12
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div'& brewer.pal.info$colorblind ,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
f=sample(col_vector, n)
pie(rep(1,n), col=f)


BacTrt<-c("< 2% abund."="#dee0df","Acidobacteria"="#01665E","Actinobacteria"="#7F3B08","Bacteroidetes"="#A50026",
          "Chloroflexi"="#35978F","Firmicutes"="#8073AC","Gemmatimonadetes"="#E7D4E8","Planctomycetes"="#5AAE61",
          "Proteobacteria"="#FDB863","Verrucomicrobia"="#FEE0B6")

#TREATMENT...............................................................................
RelPhySev1$Sample<-with(RelPhySev1,  reorder(Sample, year))

BacAbunSev<-ggplot(data=RelPhySev1, aes(x=Sample, y=Abundance, fill=Phylum))+
  geom_bar(aes(), stat="identity", position="stack") + 
  scale_fill_manual(values = BacTrt)+
  theme_bw()+ 
  theme(plot.margin = margin(1.5, 1.5, 0.8, 0.8, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 38, colour = "black"), 
        axis.title.x = element_text(size = 36, colour = "black"), 
        axis.text.y = element_text(size=30, hjust=0.5, colour = "black"),
        axis.text.x = element_text(size=30, colour = "black"),
        legend.position="bottom",
        legend.text = element_text(face = "italic", size=28),
        legend.title=element_text(size=30, face="bold"))+ 
  xlab("Severity")+
  ylab("Relative Abundance")+ 
  guides(fill=guide_legend(ncol = 3))

BacAbunSev

RelPhySev1$Sample<-with(RelPhySev1,  reorder(Sample, year))

