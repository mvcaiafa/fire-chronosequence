setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence/Diversity_analysis")

#Packages to install
install.packages("SPECIES")
install.packages("vegan")
install.packages("spaa") #needed to get EcolUtils
install.packages("tcltk")
install.packages("BiodiversityR")
install.packages("scales")
#install EcoUtils #to get the mean of the rarefied OTU table many times
install.packages("devtools")
devtools::install_github("GuillemSalazar/EcolUtils")
install.packages ("ggpubr")
install.packages("nlme")
install.packages("multcomp")
#install.packages("EcolUtils")
install.packages('tidyverse')


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
asv_Fungi <- read.csv("Fungal-table-with-taxonomy.csv", row.names=1)

head(asv_Fungi)

####Filter Taxonomy
#make last column (taxonomy) into a data frame
lastcolumn <- ncol(asv_Fungi)
taxon <- asv_Fungi[ ,lastcolumn]
id <- row.names(asv_Fungi)
dd <- data.frame(id,taxon)

#divide last row into subsets
asv_taxa1 <- ldply(str_split(string = dd$taxon, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
names(asv_taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
asv_taxa2 <- as.data.frame(lapply(asv_taxa1, gsub, pattern=" ", replacement=""))
dd <- cbind(dd[,1:2 ],asv_taxa2)

rownames(dd) <- id

dim(dd)
dd

dim(asv_Fungi)

#attach to the otu table
asv_Fungi_2 <- cbind(asv_Fungi, dd)

names(asv_Fungi_2)

#remember to remove anything not assigned to Kingdom Fungi
asv_Fungi_2$Kingdom
asv_Fungi_only<- asv_Fungi_2[which(asv_Fungi_2$Kingdom=="k__Fungi"),]

dim(asv_Fungi_only)

#remove taxonomy columns
asv_Fungi_notax  <- asv_Fungi_only[ ,1:308]

#write as csv 
write.csv(asv_Fungi_notax, "asv_Fungi_notax.csv")

asv_Fungi_notax<- read.csv("asv_Fungi_notax.csv", row.names = 1)


#look at no DNA controls
noDNA <- which(colnames(asv_Fungi_only) %in%c("Ext_Blank1_F","Ext_Blank10_F","Ext_Blank2_F","Ext_Blank3_F","Ext_Blank4_F",
                                          "Ext_Blank5_F","Ext_Blank7_F","Ext_Blank8_F","Ext_Blank9_F", "IndexBlank1_F",
                                          "PCR_Blank1_F", "PCR_Blank2_F", "PCR_Blank3_F", "PCR_Blank4_F", "Mock1_F","Mock2_F"))
                                           
###Redoing for beta diversity############################################################################################33
noDNA <- which(colnames(asv_Fungi_notax) %in%c("Ext_Blank1_F","Ext_Blank10_F","Ext_Blank2_F","Ext_Blank3_F","Ext_Blank4_F",
                                               "Ext_Blank5_F","Ext_Blank7_F","Ext_Blank8_F","Ext_Blank9_F", "IndexBlank1_F",
                                               "PCR_Blank1_F", "PCR_Blank2_F", "PCR_Blank3_F", "PCR_Blank4_F", "Mock1_F","Mock2_F"))                                           

#make otu table of no DNA controls
ASV_noDNA <-asv_Fungi_notax[ ,noDNA]
colSums(ASV_noDNA)

##Extraction and PCR blanks have very low read counts. No contaminated!! :)

ASV_Fungi_notax <- asv_Fungi_notax[ ,-noDNA]


########alpha diversity analysis
#transform table first
dim(ASV_Fungi_notax)
ASV_Fungi_trans <- t(ASV_Fungi_notax[ ,1:292])

#employ getrowsums
rowSums(ASV_Fungi_trans)

#remove samples with not enough sequences
ASV_Fungi_trans_2 <- ASV_Fungi_trans[rowSums(ASV_Fungi_trans)>=6124,]


#rarefy ASV table to remove low abundance sequences but keep rest, normalize to 6124 seq/sample, normalize 1x
ASV_Fungi.e6124 <- rrarefy(ASV_Fungi_trans_2 , 6124)

#use EcoUtils function to normalize 100 times and get mean
ASV_Fungi.e6124_mean <- rrarefy.perm(ASV_Fungi_trans_2 , sample =6124, n = 100, round.out = T)

#test if this matters - 97% correlated, it doesnt matter
mantel(ASV_Fungi.e6124_mean, ASV_Fungi.e6124)
plot(ASV_Fungi.e6124_mean ~ ASV_Fungi.e6124)
#Mantel statistic r: 0.9987 
#Significance: 0.001 
#Upper quantiles of permutations (null model):
#  90%      95%    97.5%      99% 
# 0.000384 0.001528 0.003224 0.007781 
#Permutation: free
#Number of permutations: 999
#Seems like we are fine with rarefying 1x

write.csv(ASV_Fungi.e6124, "ASV_Fungi_rarified.csv")

#remove ASVs that got zero reads after normalization
zeroes <- which(colSums(ASV_Fungi.e6124)==0)

#renove these columns 
ASV_Fungi.e6124_nozeroes <- ASV_Fungi.e6124[ ,-zeroes]
dim(ASV_Fungi.e6124_nozeroes) #281 Samples 7754 ASVs
dim(ASV_Fungi_trans) #292 Samples 11053 ASVs

#write csv

write.csv(ASV_Fungi.e6124_nozeroes, "ASV_Fungi.e6124.csv")

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(ASV_Fungi.e6124_nozeroes, method="exact")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", ylim=c(0,10000))

#figure out how to plot number of ASVs against number of sequences

#now calculate richness with BioDiversityR package
#run function once to see what it does 
ASV_Fungi_richness <- estimateR(ASV_Fungi.e6124_nozeroes)
#run function on transformed table to switch rows and columns
ASV_Fungi_richness<- t(estimateR(ASV_Fungi.e6124_nozeroes))


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
estimates_plot_function(ASV_Fungi.e6124_nozeroes,"ASV_e6124")

##############################
#other ways to calculate species richness
##############################

# get species richness 
asv.H <- diversity(ASV_Fungi.e6124_nozeroes) # Shannon entropy
asv.N1 <- exp(otu.H ) ## Shannon number of diversity
asv.Simp <- diversity(ASV_Fungi.e6124_nozeroes, "inv") ## Simpson Diversity

#make data frame of shannon entropy, shannon diversity, simpson diversity
asv.richness <- data.frame(asv.H,asv.N1,asv.Simp)

#add these to S obs, chao1, ACE
ASV_Fungi_e6124_richness <- cbind(asv.richness,ASV_Fungi_richness)

write.csv(ASV_Fungi_e6124_richness, "ASV_Fungi_e6124_richness.csv")

#test if data are normal distributed, if significant not normal
shapiro.test(ASV_Fungi_e6124_richness$S.obs)
histogram(ASV_Fungi_e6124_richness$S.obs)
#Data not normally distributed!! Use GLMM 


########################################################################################
########Make some figures
########################################################################################
#read in dataframe with METADATA
metadata <- read.csv("metadata_diversity.csv", row.names=1)

#give columns an id column for the sample ID to be able to join them by
ASV_Fungi_e6124_richness$sample_id <- row.names(ASV_Fungi_e6124_richness) 

#reorder sample ID alphabetically
ASV_Fungi_e6124_richness_match<-ASV_Fungi_e6124_richness[order(ASV_Fungi_e6124_richness$sample_id),]

#check that names match, if they don't match, use a matching function
row.names(ASV_Fungi_e6124_richness_match) == row.names(metadata)

#check classes
class(ASV_Fungi_e6124_richness$sample_id)

#give columns the same name
metadata$sample_id <- row.names(metadata)

#do a left join to join by id to maintain only rows in rarefied asv table and make sure rows match
fungidata <- left_join(ASV_Fungi_e6124_richness_match, metadata, by="sample_id")

dim(fungidata) #281, 15 (281 samples, 15 is joined meta and richness data)
dim(ASV_Fungi_e6124_richness_match) #281, 9 (rows match, 9 is just richness data)

#write as csv 
write.csv(fungidata, "ColoradoFungi_richnessandmetadata.csv")
fungidata<-read.csv("ColoradoFungi_richnessandmetadata.csv", row.names=1)

names(fungidata)
head(fungidata)

str(fungidata)
names(fungidata)

#reorder by year
fungidata<-fungidata[order(fungidata$year),]

fungidata$fire<-with(fungidata,  reorder(fire, year))

Fig_Sobs <- ggplot(fungidata, aes(x=fire, y=S.obs, group=severity, col=severity, label="Sample"))+
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

pdf("Sobsfun.pdf", height=6, width=10)
Fig_Sobs_Facet
dev.off()

Fig_Sobs_Facet 

Fig_Sobs_stats<-Fig_Sobs_Facet + stat_compare_means(method = "kruskal.test", label="p.signif")

Fig_Sobs_stats

ggexport(Fig_Sobs_stats,"Fig_Sobs_stats.pdf", width=8, height = 5)

#######################################
###########NMDS#######################
#######################################

ASV_Fungi.e6124_nozeroes<-read.csv("ASV_Fungi.e6124.csv", row.names=1)

#verify both tables are in teh same order
row.names(ASV_Fungi.e6124_nozeroes)==row.names(metadata)

#reorder  alphabetically
ASV_Fungi.e6124_order<-ASV_Fungi.e6124_nozeroes[order(row.names(ASV_Fungi.e6124_nozeroes)),]

row.names(ASV_Fungi.e6124_order)==row.names(metadata)

#######################################
####NMDS Raup-Crick dissimilarity#####

ASVp<-1*(ASV_Fungi.e6124_order!=0)

nmds.m<- metaMDS(ASVp,distance = "raup", k=3, trymax=200)

plot(nmds.m$points[,1:2], type = "n" )
text(nmds.m$points, labels=metadata$severity)

par(mar=c(5,4,2,2),xpd=T)

plot(nmds.m$points[,1:2], type = "n" , xlab = "NMDS1", ylab = "NMDS2")

#By Fire
points(nmds.m$points[metadata$fire == "badger_creek",],  pch = 19, col = "#D73027") 
points(nmds.m$points[metadata$fire == "ryan",],  pch = 19,col = "#FDAE61") 
points(nmds.m$points[metadata$fire == "beaver_creek",],  pch = 19,col = "#4575B4") 
points(nmds.m$points[metadata$fire == "mullen",],  pch = 19,col = "green")
points(nmds.m$points[metadata$fire == "church_park",],  pch = 19,col = "#313695")

#By severity

plot(nmds.m$points[,1:2], type = "n" , xlab = "NMDS1", ylab = "NMDS2")

points(nmds.m$points[metadata$severity == "high",],  pch = 19, col = "#D73027") 
points(nmds.m$points[metadata$severity == "low",],  pch = 19,col = "#FDAE61") 
points(nmds.m$points[metadata$severity == "control",],  pch = 19,col = "#4575B4") 


##############################################################################
##############################################################################
###functions for error bars, put right at the top of your code################
##############################################################################
##############################################################################

error.bar.x <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x+upper,y, x-lower, y, angle=90, code=3, length=length, ...)
}

error.bar.y <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

##############################################################################
##############################################################################


####NMDS Bray-curtis dissimilarity#####

nmds.bray<- metaMDS(ASV_Fungi.e6124_order,distance = "bray", k=4, trymax=300)

text(-0.55,0.6, "Stress = 0.1542", cex=1)

###Other ways to do NMDS

#make dataframe of NMDS scores
scores <- data.frame(scores(nmds.bray$points))

#add scores to metadata
bray_all <- cbind(metadata, scores)

nmds.bray$stress

NMDS_all<-ggplot(bray_all, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "All fires")+
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
  annotate(geom = "text", y=0.8, x=-0.4, label= "Stress=0.125")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_all

##Try dividing matrix by fires####
##################################
##Church Park

ASV_CP<- ASV_Fungi.e6124_order[grep("CP",row.names(ASV_Fungi.e6124_order)),]
meta_CP<- metadata[which(metadata$fire=="church_park"),]

row.names(ASV_CP)==row.names(meta_CP)

nmds.CP<- metaMDS(ASV_CP,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresCP <- data.frame(scores(nmds.CP$points))

#add scores to metadata
bray_CP <- cbind(meta_CP, scoresCP)

nmds.CP$stress

NMDS_ChurchPark<-ggplot(bray_CP, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Church Park")+
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
  annotate(geom = "text", y=1.2, x=0.5, label= "Stress= 0.112")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_ChurchPark

##Beaver Creek

ASV_BV<- ASV_Fungi.e6124_order[grep("BV",row.names(ASV_Fungi.e6124_order)),]
meta_BV<- metadata[which(metadata$fire=="beaver_creek"),]

row.names(ASV_BV)==row.names(meta_BV)

nmds.BV<- metaMDS(ASV_BV,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresBV <- data.frame(scores(nmds.BV$points))

#add scores to metadata
bray_BV <- cbind(meta_BV, scoresBV)

nmds.BV$stress

NMDS_BeaverCreek<-ggplot(bray_BV, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
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
  annotate(geom = "text", y=1.5, x=2, label= "Stress= 0.130")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_BeaverCreek


##Badged Creek

ASV_BC<- ASV_Fungi.e6124_order[grep("BC",row.names(ASV_Fungi.e6124_order)),]
meta_BC<- metadata[which(metadata$fire=="badger_creek"),]

row.names(ASV_BC)==row.names(meta_BC)

nmds.BC<- metaMDS(ASV_BC,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresBC <- data.frame(scores(nmds.BC$points))

#add scores to metadata
bray_BC <- cbind(meta_BC, scoresBC)

nmds.BC$stress

NMDS_BadgerCreek<-ggplot(bray_BC, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
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
  annotate(geom = "text", y=0.8, x=-0.5, label= "Stress= 0.119")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_BadgerCreek

##Ryan

ASV_R<- ASV_Fungi.e6124_order[grep("R",row.names(ASV_Fungi.e6124_order)),]
meta_R<- metadata[which(metadata$fire=="ryan"),]

row.names(ASV_R)==row.names(meta_R)

nmds.R<- metaMDS(ASV_R,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresR <- data.frame(scores(nmds.R$points))

#add scores to metadata
bray_R <- cbind(meta_R, scoresR)

nmds.R$stress

NMDS_Ryan<-ggplot(bray_R, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Ryan")+
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
  annotate(geom = "text", y=0.8, x=-0.5, label= "Stress= 0.106")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_Ryan

##Mullen

ASV_M<- ASV_Fungi.e6124_order[grep("M",row.names(ASV_Fungi.e6124_order)),]
meta_M<- metadata[which(metadata$fire=="mullen"),]

row.names(ASV_M)==row.names(meta_M)

nmds.M<- metaMDS(ASV_M,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresM <- data.frame(scores(nmds.M$points))

#add scores to metadata
bray_M <- cbind(meta_M, scoresM)

nmds.M$stress

NMDS_Mullen<-ggplot(bray_M, aes(x=MDS1, y=MDS2, col=severity, shape=depth, group=severity))+ 
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
  annotate(geom = "text", y=0.8, x=-0.5, label= "Stress= 0.101")+
  scale_color_manual(labels=c("Control","High","Low"), #manual labels for legend
                     values=c("#4575B4", "#D73027","#FDAE61" )) 

NMDS_Mullen

#Arrange figures

NMDS_figures<- ggarrange(NMDS_all,NMDS_ChurchPark, NMDS_BeaverCreek,NMDS_BadgerCreek,NMDS_Ryan,NMDS_Mullen, labels = c("A", "B","C","D", "E","F"), ncol = 3, nrow = 2, common.legend=TRUE)
NMDS_figures



##Try dividing matrix by severity####
##################################
##Control

ASV_Control<- ASV_Fungi.e6124_order[grep("_C",row.names(ASV_Fungi.e6124_order)),]
meta_Control<- metadata[which(metadata$severity=="control"),]

row.names(ASV_Control)==row.names(meta_Control)

nmds.Control<- metaMDS(ASV_Control,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresControl <- data.frame(scores(nmds.Control$points))

#add scores to metadata
bray_Control <- cbind(meta_Control, scoresControl)

#reorder  alphabetically
bray_Control<-bray_Control[order(bray_Control$year),]

nmds.Control$stress

NMDS_Control<-ggplot(bray_Control, aes(x=MDS1, y=MDS2, col=fire, shape=depth, group=fire))+ 
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
  annotate(geom = "text", y=0.6, x=0.6, label= "Stress= 0.175")+
  scale_color_manual(values=c("#D73027","#542788","#74ADD1","#5AAE61","#FDAE61"),
                     labels=c("Badger Creek","Beaver Creek","Church Park","Mullen","Ryan")) #manual labels for legend


NMDS_Control

##Low

ASV_Low<- ASV_Fungi.e6124_order[grep("_L",row.names(ASV_Fungi.e6124_order)),]
meta_Low<- metadata[which(metadata$severity=="low"),]

row.names(ASV_Low)==row.names(meta_Low)

nmds.Low<- metaMDS(ASV_Low,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresLow <- data.frame(scores(nmds.Low$points))

#add scores to metadata
bray_Low <- cbind(meta_Low, scoresLow)

nmds.Low$stress

NMDS_Low<-ggplot(bray_Low, aes(x=MDS1, y=MDS2, col=fire, shape=depth, group=fire))+ 
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
  stat_ellipse()+  #include ellipse around the group
  annotate(geom = "text", y=0.5, x=1, label= "Stress= 0.171")+
  scale_color_manual(labels=c("Church Park","Beaver Creek","Badger Creek", "Ryan","Mullen"), #manual labels for legend
                     values=c("#74ADD1","#542788","#D73027","#FDAE61","#5AAE61" )) 

NMDS_Low


##High

ASV_High<- ASV_Fungi.e6124_order[grep("_H",row.names(ASV_Fungi.e6124_order)),]
meta_High<- metadata[which(metadata$severity=="high"),]

row.names(ASV_High)==row.names(meta_High)

nmds.High<- metaMDS(ASV_High,distance = "bray", k=3, trymax=200)

#make dataframe of NMDS scores
scoresHigh <- data.frame(scores(nmds.High$points))

#add scores to metadata
bray_High <- cbind(meta_High, scoresHigh)
#reorder  alphabetically
bray_High<-bray_High[order(bray_High$year),]

nmds.High$stress

NMDS_High<-ggplot(bray_High, aes(x=MDS1, y=MDS2, col=fire, shape=depth, group=fire))+ 
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
  stat_ellipse()+  #include ellipse around the group
  annotate(geom = "text", y=0.5, x=-0.5, label= "Stress= 0.165")+
  scale_color_manual(labels=c("Church Park","Beaver Creek","Badger Creek", "Ryan","Mullen"), #manual labels for legend
                   values=c("#74ADD1","#542788","#D73027","#FDAE61","#5AAE61" )) 

NMDS_High


NMDS_figures<- ggarrange(NMDS_Control,NMDS_Low, NMDS_High, labels = c("A", "B","C"), ncol =3, nrow = 1, common.legend=TRUE)
NMDS_figures


##Difference in Depth by severity

NMDS_depth_by_Severity<-ggplot(bray_all, aes(x=NMDS1, y=NMDS2, col=depth, shape=depth, group=depth))+ 
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Fungal Bray-Curtis Dissimilarity")+
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



