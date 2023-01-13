setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence/Diversity_analysis_B")

#Packages to load 
library(devtools)
library(tidyverse)
library(EcolUtils)
library(SPECIES)
library(vegan)
library(tcltk)
library(BiodiversityR)
library(scales)
library(ggpubr)
#load library for nonlinear mixed effect models
library(nlme)
library(multcomp)
library("stringr")
library("plyr")
library("dplyr")
library("ggplot2")

#read in dataframes
asv_Bacteria <- read.csv("Bacteria-table-with-taxonomy.csv", row.names=1)
head(asv_Bacteria)

####Filter Taxonomy
#make last column (taxonomy) into a data frame
lastcolumn <- ncol(asv_Bacteria)
taxon <- asv_Bacteria[ ,lastcolumn]
id <- row.names(asv_Bacteria)
dd <- data.frame(id,taxon)

#divide last row into subsets
asv_taxa1 <- ldply(str_split(string = dd$taxon, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(asv_taxa1) <- c("Kingdom", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6")
asv_taxa2 <- as.data.frame(lapply(asv_taxa1, gsub, pattern=" ", replacement=""))
dd <- cbind(dd[,1:2 ],asv_taxa2)

rownames(dd) <- id

dim(dd)
dd

dim(asv_Bacteria)

#remove taxonomy columns
asv_Bacteria_notax  <- asv_Bacteria[ ,1:308]
dim(asv_Bacteria_notax)

#look at no DNA controls
noDNA <- which(colnames(asv_Bacteria_notax) %in%c("Ext_Blank1_B","Ext_Blank10_B","Ext_Blank2_B","Ext_Blank3_B","Ext_Blank4_B",
                                              "Ext_Blank5_B","Ext_Blank7_B","Ext_Blank8_B", "IndexBlank2_B",
                                              "PCR_Blank1_B", "PCR_Blank2_B", "PCR_Blank3_B", "PCR_Blank4_B", "PCR_Blank5_B","Mock1_B","Mock2_B"))

#make otu table of no DNA controls
ASV_noDNA <-asv_Bacteria_notax[ ,noDNA]
colSums(ASV_noDNA)

##Extraction and PCR blanks have very low read counts. No contaminated!! :)

ASV_Bacteria_notax <- asv_Bacteria_notax[ ,-noDNA]

########alpha diversity analysis
#transform table first
dim(ASV_Bacteria_notax)
ASV_Bacteria_trans <- t(ASV_Bacteria_notax[ ,1:292])

#employ getrowsums
rowSums(ASV_Bacteria_trans)

#remove samples with not enough sequences
ASV_Bacteria_trans_2 <- ASV_Bacteria_trans[rowSums(ASV_Bacteria_trans)>=5321,]
rowSums(ASV_Bacteria_trans_2)
Z<-ASV_Bacteria_trans[rowSums(ASV_Bacteria_trans)<5321,]

#rarefy ASV table to remove low abundance sequences but keep rest, normalize to 5321 seq/sample, normalize 1x
ASV_Bacteria.e5321 <- rrarefy(ASV_Bacteria_trans_2 , 5321)

#use EcoUtils function to normalize 100 times and get mean
ASV_Bacteria.e5321_mean <- rrarefy.perm(ASV_Bacteria_trans_2 , sample =5321, n = 100, round.out = T)

#test if this matters - 97% correlated, it doesnt matter
mantel(ASV_Bacteria.e5321_mean, ASV_Bacteria.e5321)
plot(ASV_Bacteria.e5321_mean ~ ASV_Bacteria.e5321)

#Mantel statistic r: 0.9805 
#Significance: 0.001 
#Upper quantiles of permutations (null model):
#  90%     95%   97.5%     99% 
#  0.00531 0.00952 0.01526 0.03118 
#Permutation: free
#Number of permutations: 999
#Seems like we are fine with rarefying 1x

#remove ASVs that got zero reads after normalization
zeroes <- which(colSums(ASV_Bacteria.e5321)==0)

#remove these columns 
ASV_Bacteria.e5321_nozeroes <- ASV_Bacteria.e5321[ ,-zeroes]
dim(ASV_Bacteria.e5321_nozeroes) #281 Samples 29601 ASVs
dim(ASV_Bacteria_trans) #292 Samples 35698 ASVs

write.csv(ASV_Bacteria.e5321_nozeroes, "ASV_Bacteria.e5321.csv")

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(ASV_Bacteria.e5321_nozeroes, method="random")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", ylim=c(0,10000))

#figure out how to plot number of ASVs against number of sequences

#now calculate richness with BioDiversityR package
#run function once to see what it does 
ASV_Bacteria_richness <- estimateR(ASV_Bacteria.e5321_nozeroes)
#run function on transformed table to switch rows and columns
ASV_Bacteria_richness<- t(estimateR(ASV_Bacteria.e5321_nozeroes))


#invoke pvaluelegend to put pvalue on figure (Sydney's GitHub)
pvaluelegend <- function(stat,pval){
  rp2 = vector('expression',2)
  
  rp2[1] = substitute(expression("Pearson's r"== MYVALUE),
                      list(MYVALUE = format(stat, dig=2)))[2]
  rp2[2] = substitute(expression(italic(P)== MYOTHERVALUE),
                      list(MYOTHERVALUE = (ifelse(pval>0.0001,format(pval, digits = 3), "0.0001"))))[2]
  
  legend("topleft", legend = rp2, bty = 'n', text.col="black", cex=1.3) 
}

#make a function to get estimates and make plots
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_fungal_richness_correlations_MiSeq",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
  #perform correlation test
  cortest1 <- cor.test(estimates2[,2],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest1$estimate, cortest1$p.value)
  mtext("A",side=3,adj=0)
  #plot S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",ylab="ACE",col=alpha("black", 0.5),pch=16)
  #perform correlation test
  cortest2 <- cor.test(estimates2[,4],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest2$estimate, cortest2$p.value)
  mtext("B",side=3,adj=0)
  dev.off()
  
}

#run function
estimates_plot_function(ASV_Bacteria.e5321_nozeroes,"ASV_e5321")

##############################
#other ways to calculate species richness
##############################

# get species richness fo EMF not rarefied
asv.H <- diversity(ASV_Bacteria.e5321_nozeroes) # Shannon entropy
asv.N1 <- exp(asv.H ) ## Shannon number of diversity
asv.Simp <- diversity(ASV_Bacteria.e5321_nozeroes, "inv") ## Simpson Diversity

#make data frame of shannon entropy, shannon diversity, simpson diversity
asv.richness <- data.frame(asv.H,asv.N1,asv.Simp)

#add these to S obs, chao1, ACE
ASV_Bacteria_e5321_richness <- cbind(asv.richness,ASV_Bacteria_richness)

write.csv(ASV_Bacteria_e5321_richness, "ASV_Bacteria_e5321_richness.csv")

#test if data are normal distributed, if significant not normal
shapiro.test(ASV_Bacteria_e5321_richness$S.obs)
histogram(ASV_Bacteria_e5321_richness$S.obs)
#Data not normally distributed!! Use GLMM 

########################################################################################
########Make some figures
########################################################################################
#read in dataframe with METADATA
metadata_B <- read.csv("metadata_diversity_Bacteria.csv", row.names=1)

#give columns an id column for the sample ID to be able to join them by
ASV_Bacteria_e5321_richness$sample_id <- row.names(ASV_Bacteria_e5321_richness) 

#give columns the same name
metadata_B$sample_id <- row.names(metadata_B)

#reorder sample ID alphabetically
ASV_Bacteria_e5321_richness_match<-ASV_Bacteria_e5321_richness[order(ASV_Bacteria_e5321_richness$sample_id),]

#check that names match, if they don't match, use a matching function
row.names(ASV_Bacteria_e5321_richness_match) == row.names(metadata_B)

#do a left join to join by id to maintain only rows in rarefied asv table and make sure rows match
bacterialdata <- left_join(ASV_Bacteria_e5321_richness_match, metadata_B, by="sample_id")

dim(bacterialdata) #281, 15 (281 samples, 15 is joined meta and richness data)
dim(ASV_Bacteria_e5321_richness_match) #281, 9 (rows match, 9 is just richness data)

#write as csv 
write.csv(bacterialdata, "ColoradoBacteria_richnessandmetadata.csv")
bacterialdata<-read.csv("ColoradoBacteria_richnessandmetadata.csv", row.names=1)

names(bacterialdata)
head(bacterialdata)

str(bacterialdata)

#reorder by year
bacterialdata<-bacterialdata[order(bacterialdata$year),]

bacterialdata$fire<-with(bacterialdata,  reorder(fire, year))

Fig_Sobs <- ggplot(bacterialdata, aes(x=fire, y=S.obs, group=severity, col=severity, label="Sample"))+
  theme_bw() + #make black and white
  ylab("S.obs") +  #change yaxis label
  scale_colour_manual(values = c("#4575B4","#D73027","#FDAE61"))+
  stat_summary(fun=mean,geom="point", size=4)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=.8, alpha=0.7,position = position_dodge(0.01))+
  theme(legend.position = "bottom", #put legend to the right of the graph
        legend.title = element_blank(), #remove legend titles
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=12, angle=20,hjust=1), #make y axis text larger and angled
        axis.text.x = element_text(size=12,angle=45, hjust=1)) 

Fig_Sobs

Fig_Sobs_Facet <- Fig_Sobs+ facet_wrap(~depth)

Fig_Sobs_Facet 

pdf("SobsBac.pdf", height=6, width=10)
Fig_Sobs_Facet
dev.off()

Fig_Sobs_stats<-Fig_Sobs_Facet + stat_compare_means(method = "kruskal.test", label="p.signif")

Fig_Sobs_stats


ggexport(Fig_Sobs_stats,"Fig_Sobs_stats.pdf", width=8, height = 5)


#######################################
###########NMDS#######################
#######################################

ASV_Bacteria.e5321<-read.csv("ASV_Bacteria.e5321.csv", row.names=1)
metadata_B <- read.csv("metadata_diversity_Bacteria.csv", row.names=1)

#verify both tables are in teh same order
row.names(ASV_Bacteria.e5321)==row.names(metadata_B)

#reorder  alphabetically
ASV_Bacteria.e5321_order<-ASV_Bacteria.e5321[order(row.names(ASV_Bacteria.e5321)),]

row.names(ASV_Bacteria.e5321_order)==row.names(metadata_B)


##############################################################################
##############################################################################
####################NMDS Bray-curtis dissimilarity############################

nmds.bray.bac<- metaMDS(ASV_Bacteria.e5321_order,distance = "bray", k=4, trymax=300)

###Other ways to do NMDS

#make dataframe of NMDS scores
scores_bac <- data.frame(scores(nmds.bray.bac$points))

#add scores to metadata
bray_all_bac <- cbind(metadata_B, scores_bac)

nmds.bray.bac$stress

NMDS_all_bac<-ggplot(bray_all_bac, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "All fires")+
  #xlim(-0.2,0.7)+
  ylim(-0.3,0.3)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=0.25, x=-0.3, label= "Stress=0.103")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_all_bac

##Try dividing matrix by fires####
##################################
##Church Park

ASV_B_CP<- ASV_Bacteria.e5321_order[grep("CP",row.names(ASV_Bacteria.e5321_order)),]
meta_B_CP<- metadata_B[which(metadata_B$fire=="church park"),]

row.names(ASV_B_CP)==row.names(meta_B_CP)

nmds.bac.CP<- metaMDS(ASV_B_CP,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scores_bac_CP <- data.frame(scores(nmds.bac.CP$points))

#add scores to metadata
bray_bac_CP <- cbind(meta_B_CP, scores_bac_CP)

nmds.bac.CP$stress

NMDS_Bac_ChurchPark<-ggplot(bray_bac_CP, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Church Park")+
  #xlim(-0.2,0.7)+
  ylim(-0.8,1)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=1, x=1, label= "Stress= 0.113")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_Bac_ChurchPark

##Beaver Creek

ASV_B_BV<- ASV_Bacteria.e5321_order[grep("BV",row.names(ASV_Bacteria.e5321_order)),]
meta_B_BV<- metadata_B[which(metadata_B$fire=="beaver creek"),]

row.names(ASV_B_BV)==row.names(meta_B_BV)

nmds.bac.BV<- metaMDS(ASV_B_BV,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scores_bac_BV <- data.frame(scores(nmds.bac.BV$points))

#add scores to metadata
bray_bac_BV <- cbind(meta_B_BV, scores_bac_BV)

nmds.bac.BV$stress

NMDS_Bac_BeaverCreek<-ggplot(bray_bac_BV, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Beaver Creek")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=1.5, x=2, label= "Stress= 0.106")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_Bac_BeaverCreek


##Badged Creek

ASV_B_BC<- ASV_Bacteria.e5321_order[grep("BC",row.names(ASV_Bacteria.e5321_order)),]
meta_B_BC<- metadata_B[which(metadata_B$fire=="badger creek"),]

row.names(ASV_B_BC)==row.names(meta_B_BC)

nmds.bac.BC<- metaMDS(ASV_B_BC,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scores_bac_BC <- data.frame(scores(nmds.bac.BC$points))

#add scores to metadata
bray_bac_BC <- cbind(meta_B_BC, scores_bac_BC)

nmds.bac.BC$stress

NMDS_Bac_BadgerCreek<-ggplot(bray_bac_BC, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Badger Creek")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=0.8, x=-0.5, label= "Stress= 0.120")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_Bac_BadgerCreek

##Ryan

ASV_B_R<- ASV_Bacteria.e5321_order[grep("R",row.names(ASV_Bacteria.e5321_order)),]
meta_B_R<- metadata_B[which(metadata_B$fire=="ryan"),]

row.names(ASV_B_R)==row.names(meta_B_R)

nmds.bac.R<- metaMDS(ASV_B_R,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scores_bac_R <- data.frame(scores(nmds.bac.R$points))

#add scores to metadata
bray_bac_R <- cbind(meta_B_R, scores_bac_R)

nmds.bac.R$stress

NMDS_Bac_Ryan<-ggplot(bray_bac_R, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Ryan")+
  #xlim(-0.2,0.7)+
  ylim(-0.5,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=0.5, x=-0.5, label= "Stress= 0.100")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_Bac_Ryan

##Mullen

ASV_B_M<- ASV_Bacteria.e5321_order[grep("M",row.names(ASV_Bacteria.e5321_order)),]
meta_B_M<- metadata_B[which(metadata_B$fire=="mullen"),]

row.names(ASV_B_M)==row.names(meta_B_M)

nmds.bac.M<- metaMDS(ASV_B_M,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scores_bac_M <- data.frame(scores(nmds.bac.M$points))

#add scores to metadata
bray_bac_M <- cbind(meta_B_M, scores_bac_M)

nmds.bac.M$stress

NMDS_Bac_Mullen<-ggplot(bray_bac_M, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Mullen")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=0.8, x=-0.5, label= "Stress= 0.082")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_Bac_Mullen

#Arrange figures

library(ggpubr)

NMDS_figures<- ggarrange(NMDS_all_bac,NMDS_Bac_ChurchPark, NMDS_Bac_BeaverCreek,NMDS_Bac_BadgerCreek,NMDS_Bac_Ryan,NMDS_Bac_Mullen, labels = c("A", "B","C","D", "E","F"), ncol = 3, nrow = 2, common.legend=TRUE)
NMDS_figures



##Try dividing matrix by severity####
##################################
##Control

ASV_Bac_Control<- ASV_Bacteria.e5321_order[grep("_C",row.names(ASV_Bacteria.e5321_order)),]
meta_B_Control<- metadata_B[which(metadata_B$severity=="control"),]

row.names(ASV_Bac_Control)==row.names(meta_B_Control)

nmds.Bac.Control<- metaMDS(ASV_Bac_Control,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresBacControl <- data.frame(scores(nmds.Bac.Control$points))

#add scores to metadata
bray_Bac_Control <- cbind(meta_B_Control, scoresBacControl)

nmds.Bac.Control$stress

NMDS_Bac_Control<-ggplot(bray_Bac_Control, aes(x=MDS1, y=MDS2, col=fire, shape=depth, group=fire))+ 
  geom_point(size=3) +
  theme_bw() +
  labs(title = "Control")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=0.6, x=0.6, label= "Stress= 0.117")+
  scale_color_manual(labels=c("Church Park","Beaver Creek","Badger Creek", "Ryan","Mullen"), #manual labels for legend
                     values=c("#74ADD1","#542788","#D73027","#FDAE61","#5AAE61" )) 

NMDS_Bac_Control

##Low

ASV_Bac_Low<- ASV_Bacteria.e5321_order[grep("_L",row.names(ASV_Bacteria.e5321_order)),]
meta_B_Low<- metadata_B[which(metadata_B$severity=="low"),]

row.names(ASV_Bac_Low)==row.names(meta_B_Low)

nmds.Bac.Low<- metaMDS(ASV_Bac_Low,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresBacLow <- data.frame(scores(nmds.Bac.Low$points))

#add scores to metadata
bray_Bac_Low <- cbind(meta_B_Low, scoresBacLow)

nmds.Bac.Low$stress

NMDS_Bac_Low<-ggplot(bray_Bac_Low, aes(x=MDS1, y=MDS2, col=fire, shape=depth, group=fire))+ 
  geom_point(size=3) +
  theme_bw() +
  labs(title = "Low")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=0.07, x=0.1, label= "Stress= 0.139")+
  scale_color_manual(labels=c("Church Park","Beaver Creek","Badger Creek", "Ryan","Mullen"), #manual labels for legend
                     values=c("#74ADD1","#542788","#D73027","#FDAE61","#5AAE61" )) 

NMDS_Bac_Low


##High

ASV_Bac_High<- ASV_Bacteria.e5321_order[grep("_H",row.names(ASV_Bacteria.e5321_order)),]
meta_B_High<- metadata_B[which(metadata_B$severity=="high"),]

row.names(ASV_Bac_High)==row.names(meta_B_High)

nmds.Bac.High<- metaMDS(ASV_Bac_High,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresBacHigh <- data.frame(scores(nmds.Bac.High$points))

#add scores to metadata
bray_Bac_High <- cbind(meta_B_High, scoresBacHigh)

nmds.Bac.High$stress

NMDS_Bac_High<-ggplot(bray_Bac_High, aes(x=MDS1, y=MDS2, col=fire, shape=depth, group=fire))+ 
  geom_point(size=3) +
  theme_bw() +
  labs(title = "High")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  annotate(geom = "text", y=0.5, x=0.5, label= "Stress= 0.132")+
  scale_color_manual(labels=c("Church Park","Beaver Creek","Badger Creek", "Ryan","Mullen"), #manual labels for legend
                     values=c("#74ADD1","#542788","#D73027","#FDAE61","#5AAE61" )) 

NMDS_Bac_High

NMDS_Bac_figures<- ggarrange(NMDS_Bac_Control,NMDS_Bac_Low, NMDS_Bac_High, labels = c("A", "B","C"), ncol = 3, nrow = 1, common.legend=TRUE)
NMDS_Bac_figures



##Difference in Depth by severity

NMDS_depth_by_Severity<-ggplot(bray_all_B, aes(x=NMDS1, y=NMDS2, col=depth, shape=depth, group=depth))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Bacterial Bray-Curtis Dissimilarity")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 16, colour = "black"), #make 
        axis.text=element_text(size=16,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  scale_color_manual(labels=c("Shallow","Deep"), #manual labels for legend
                     values=c("#4575B4", "#D73027" )) 

NMDS_depth_by_Severity

NMDS_severity <-NMDS_depth_by_Severity + facet_wrap(~severity, scales="free")

NMDS_severity 

citation()


