setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence/Diversity_analysis")
############################################################################################
#...................................Richness GLMM..........................................#
############################################################################################
library(lme4)
library(MASS)
library(nlme)

fungidata<-read.csv("ColoradoFungi_richnessandmetadata.csv", row.names=1)

# Faster way to identify the number of levels in your data set
str(fungidata)
names(fungidata)

fungidata$year <- as.numeric(fungidata$year)
fungidata$S.obs <- as.numeric(fungidata$S.obs)

fungidata$year<-scale(fungidata$year)

#test if data are normal distributed, if significant not normal
shapiro.test(fungidata$S.obs)
hist(fungidata$S.obs)
#Data not normally distributed!! Use GLMM 

#* * * TREATMENT EFFECTS .......................................................
#Checking mean and variance of the data to see if its NegBinom..................
df <- fungidata %>%
  group_by(severity) %>%
  summarise(means = mean(S.obs),
            variance = var(S.obs));df #suggest negative binomial

##Start with a poisson model ....................................................
Pois1 <- glm(S.obs ~ severity, data = fungidata, family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion  ods=23.72845
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

#Zero inflaited model ..........................................................
sum(fungidata$S.obs == 0) #0

#Test different identity link to see if disperssion decreases ..................

Pois2 <- glm(S.obs ~ severity , data = fungidata,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #ods=23.72845

Pois3 <- glm(S.obs ~ severity , data = fungidata,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#ods=23.72845

Gam <- glm(S.obs ~ severity ,data = fungidata,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam#0.174121

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
nb1 <- glm.nb(S.obs ~ severity,data = fungidata)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#1.035907
plot(nb1)

AIC(Pois1, Pois2, Pois3, Gam, nb1)#gamma better 

#Test negative model
qpcr_neg1 <- glmer(S.obs ~ (1|plot), family= "Gamma",data =fungidata);summary(qpcr_neg1)


#BACKWARD MODEL SELECTION........................................................................
# * * FULL MODEL ..............................................................................
Full<-glmer( S.obs~fire*year + fire*severity + fire*depth + year*severity + year*depth + 
               severity*depth +fire+year+severity+depth+(1|plot), family= "Gamma", data = fungidata);summary(Full)



Null<-glmer(S.obs ~ (1|plot), family= "Gamma",data =fungidata);summary(Null)
anova(Full,Null)#Full model better, AIC = 2927.1

#Start backwards model selection..............................................................

#Remove fire*year

Mod2<-glmer( S.obs~ fire*severity + fire*depth + year*severity + year*depth + 
                severity*depth +fire+year+severity+depth+(1|plot), family= "Gamma", data = fungidata);summary(Mod2)

anova(Full, Mod2)##Same AIC, not significant, remove

#Remove year*severity

Mod3<-glmer(S.obs~ fire*severity + fire*depth +  year*depth + 
               severity*depth +fire+year+severity+depth+(1|plot), family= "Gamma", data = fungidata); summary(Mod3)

anova(Mod3,Mod2)##Same AIC, not significant, remove

#Remove year*depth
Mod4<-glmer(S.obs~ fire*severity + fire*depth + severity*depth +fire+year+severity+depth+
              (1|plot), family= "Gamma", data = fungidata); summary(Mod4)

anova(Mod3,Mod4)##Same AIC, not significant, remove

#REmove severity*depth +

Mod5<-glmer(S.obs~ fire*severity + fire*depth + fire+year+severity+depth+
              (1|plot), family= "Gamma", data = fungidata); summary(Mod5)

anova(Mod4,Mod5)##Mod5 better, significant, DO NOT remove

#Remove fire*depth +

Mod6<-glmer(S.obs~ fire*severity + severity*depth +fire+year+severity+depth+
              (1|plot), family= "Gamma", data = fungidata); summary(Mod6)

anova(Mod4, Mod6)##Mod6 better, not significant, remove

#Remove fire*severity

Mod7<-glmer(S.obs~severity*depth +fire+year+severity+depth+
              (1|plot), family= "Gamma", data = fungidata); summary(Mod7)

anova(Mod6, Mod7)##Mod6 better, significant, DO NOT remove

##Remove year+

Mod8<-glmer(S.obs~ fire*severity + severity*depth +fire+severity+depth+
              (1|plot), family= "Gamma", data = fungidata); summary(Mod8)

anova(Mod6, Mod8)###Same AIC, not significant, remove

##Remove fire

Mod9<-glmer(S.obs~ fire*severity + severity*depth +severity+depth+
              (1|plot), family= "Gamma", data = fungidata); summary(Mod9)

anova(Mod8, Mod9)###Same AIC, not significant, remove

#Remove severity

Mod10<-glmer(S.obs~ fire*severity + severity*depth +depth+(1|plot), family= "Gamma", data = fungidata); summary(Mod10)

anova(Mod9, Mod10)###Same AIC, not significant, remove


#Remove depth

Mod11<-glmer(S.obs~ fire*severity + severity*depth +(1|plot), family= "Gamma", data = fungidata); summary(Mod11)

anova(Mod11, Mod10)###Same AIC, not significant, remove


#Compare all models....................................................................
anova(Full, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10,Mod11)
AIC(Full, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10,Mod11)




FunRichModel<-glmer(S.obs~ fire*severity + severity*depth +(1|plot), family= "Gamma", data = fungidata);summary(FunRichModel)
anova(FunRichModel)

FunRichModel_stat<-drop1(FunRichModel, test = "Chi")
write_csv(FunRichModel_stat, "FunRichStat.csv")


library(r2glmm)

R2FunRichModel<- r2beta(model=FunRichModel,partial=TRUE,method='sgv')
plot(FunRichModel)
Anova(FunRichModel)
summary(FunRichModel)

#######
##Percent Change
##
library(dplyr)

Fun_ave <- fungidata %>%
  group_by(fire, severity,depth) %>%
  summarise(means = mean(S.obs),
            variance = var(S.obs))
  
write.csv(Fun_ave, "FunRich_average.csv")

##Change<-fungidata %>%
 # group_by(fire, severity,depth) %>%
#  summarise(means = mean(S.obs),variance = var(S.obs)) %>%
 # mutate(change= (means/sum(means))*100)


############################################################################################
#..................................Community Permanova..... ...............................#
############################################################################################

library(vegan)
library(ape)

ASV_Fungi<-read.csv("ASV_Fungi.e6124.csv", row.names=1)
meta_fungi<- read.csv("metadata_diversity.csv", row.names=1)

row.names(ASV_Fungi)==row.names(meta_fungi)

#reorder  alphabetically
ASV_Fungi<-ASV_Fungi[order(row.names(ASV_Fungi)),]

##ADONIS (Bray-curtis)

str(meta_fungi)
meta_fungi$year<-as.factor(meta_fungi$year)

AdonisFungi<-adonis2(ASV_Fungi~meta_fungi$severity + meta_fungi$fire+meta_fungi$depth,permutations=9999, method="bray")

AdonisFungi2<-adonis(ASV_Fungi~meta_fungi$severity*meta_fungi$fire+meta_fungi$severity*meta_fungi$depth,permutations=9999, method="bray")

AdonisFungi2

##Betadisper

dist_mat<-vegdist(ASV_Fungi,distance="bray")
beta_severity<-betadisper(dist_mat, meta_fungi$severity)
anova(beta_severity)

beta_site<-betadisper(dist_mat, meta_fungi$fire)
anova(beta_site)

beta_depth<-betadisper(dist_mat, meta_fungi$depth)
anova(beta_depth)

##########################
#####ADONIS by site#######
##########################

###Church Park

meta_fungi_CP<-subset(meta_fungi, fire=="church_park")

meta_fungi_CP$sample_id<-row.names(meta_fungi_CP)

ASV_Fungi$sample_id<-row.names(ASV_Fungi)


ASV_Fungi_CP<-semi_join(ASV_Fungi,meta_fungi_CP, by="sample_id", copy=F )
dim(ASV_Fungi_CP)

ASV_Fungi_CP<-ASV_Fungi_CP[,-7755]

row.names(ASV_Fungi_CP)==row.names(meta_fungi_CP)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_CP)==0)

#renove these columns 
ASV_Fungi_ChurchPark <- ASV_Fungi_CP[ ,-zeroes]
dim(ASV_Fungi_CP) #52 Samples 7754 ASVs
dim(ASV_Fungi_ChurchPark) #52 Samples 3101 ASVs

AdonisFungiCP<-adonis2(ASV_Fungi_ChurchPark~meta_fungi_CP$severity*meta_fungi_CP$depth+meta_fungi_CP$severity+meta_fungi_CP$depth,permutations=9999, method="bray")
AdonisFungiCP
capture.output(AdonisFungiCP, file = "AdonisFungi_CP.csv")

###Beaver Creek

meta_fungi_BV<-subset(meta_fungi, fire=="beaver creek")

meta_fungi_BV$sample_id<-row.names(meta_fungi_BV)

ASV_Fungi$sample_id<-row.names(ASV_Fungi)

ASV_Fungi_BV<-semi_join(ASV_Fungi,meta_fungi_BV, by="sample_id", copy=F )
dim(ASV_Fungi_BV)

ASV_Fungi_BV<-ASV_Fungi_BV[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_BV)==0)

#renove these columns 
ASV_Fungi_BeaverCreek <- ASV_Fungi_BV[ ,-zeroes]
dim(ASV_Fungi_BV) #49 Samples 7754 ASVs
dim(ASV_Fungi_BeaverCreek) #49 Samples 8500 ASVs

AdonisFungiBV<-adonis(ASV_Fungi_BeaverCreek~meta_fungi_BV$severity*meta_fungi_BV$depth+meta_fungi_BV$severity+meta_fungi_BV$depth,permutations=9999, method="bray")
AdonisFungiBV
capture.output(AdonisFungiBV, file = "AdonisFungi_BV.csv")

###Badger Creek
meta_fungi_BC<-subset(meta_fungi, fire=="badger creek")

meta_fungi_BC$sample_id<-row.names(meta_fungi_BC)


ASV_Fungi_BC<-semi_join(ASV_Fungi,meta_fungi_BC, by="sample_id", copy=F )
dim(ASV_Fungi_BC)

ASV_Fungi_BC<-ASV_Fungi_BC[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_BC)==0)

#renove these columns 
ASV_Fungi_BadgerCreek <- ASV_Fungi_BC[ ,-zeroes]
dim(ASV_Fungi_BC) #72 Samples 7754 ASVs
dim(ASV_Fungi_BadgerCreek) #72 Samples 2403 ASVs

AdonisFungiBC<-adonis(ASV_Fungi_BadgerCreek~meta_fungi_BC$severity*meta_fungi_BC$depth+meta_fungi_BC$severity+meta_fungi_BC$depth,permutations=9999, method="bray")
AdonisFungiBC
capture.output(AdonisFungiBC, file = "AdonisFungi_BC.csv")

###Ryan

meta_fungi_R<-subset(meta_fungi, fire=="ryan")

meta_fungi_R$sample_id<-row.names(meta_fungi_R)

ASV_Fungi_R<-semi_join(ASV_Fungi,meta_fungi_R, by="sample_id", copy=F )
dim(ASV_Fungi_R)

ASV_Fungi_R<-ASV_Fungi_R[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_R)==0)

#renove these columns 
ASV_Fungi_Ryan <- ASV_Fungi_R[ ,-zeroes]
dim(ASV_Fungi_R) #54 Samples 7754 ASVs
dim(ASV_Fungi_Ryan) #54 Samples 1818 ASVs

AdonisFungiR<-adonis(ASV_Fungi_Ryan~meta_fungi_R$severity*meta_fungi_R$depth+meta_fungi_R$severity+meta_fungi_R$depth,permutations=9999, method="bray")
AdonisFungiR

capture.output(AdonisFungiR, file = "AdonisFungi_R.csv")


###Mullen

meta_fungi_M<-subset(metadata_B, fire=="mullen")

meta_fungi_M$sample_id<-row.names(meta_fungi_M)

ASV_Fungi_M<-semi_join(ASV_Fungi,meta_fungi_M, by="sample_id", copy=F )
dim(ASV_Fungi_M)

ASV_Fungi_M<-ASV_Fungi_M[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_M)==0)

#renove these columns 
ASV_Fungi_Mullen <- ASV_Fungi_M[ ,-zeroes]
dim(ASV_Fungi_M) #53 Samples 7754 ASVs
dim(ASV_Fungi_Mullen) #54 Samples 2430 ASVs

AdonisFungiM<-adonis(ASV_Fungi_Mullen~meta_fungi_M$severity*meta_fungi_M$depth+meta_fungi_M$severity+meta_fungi_M$depth,permutations=9999, method="bray")
AdonisFungiM

capture.output(AdonisFungiM, file = "AdonisFungi_M.csv")

##############################
#####ADONIS by severity#######
##############################

###Control

meta_fungi_Control<-subset(meta_fungi, severity=="control")

meta_fungi_Control$sample_id<-row.names(meta_fungi_Control)

ASV_Fungi$sample_id<-row.names(ASV_Fungi)


ASV_Fungi_Control<-semi_join(ASV_Fungi,meta_fungi_Control, by="sample_id", copy=F )
dim(ASV_Fungi_Control)

ASV_Fungi_Control<-ASV_Fungi_Control[,-7755]

row.names(ASV_Fungi_Control)==row.names(meta_fungi_Control)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_Control)==0)

#renove these columns 
ASV_Fungi_C <- ASV_Fungi_Control[ ,-zeroes]
dim(ASV_Fungi_Control) #95 Samples 7754 ASVs
dim(ASV_Fungi_C) #95 Samples 4109 ASVs

AdonisFungiControl<-adonis2(ASV_Fungi_C~meta_fungi_Control$fire*meta_fungi_Control$depth+meta_fungi_Control$fire+meta_fungi_Control$depth,permutations=9999, method="bray")
AdonisFungiControl
capture.output(AdonisFungiControl, file = "AdonisFungi_Control.csv")

##Betadisper

dist_matC<-vegdist(ASV_Fungi_C,distance="bray")

beta_site<-betadisper(dist_matC, meta_fungi_Control$fire)
anova(beta_site)

beta_depth<-betadisper(dist_matC, meta_fungi_Control$depth)
anova(beta_depth)

#Low

meta_fungi_Low<-subset(meta_fungi, severity=="low")

meta_fungi_Low$sample_id<-row.names(meta_fungi_Low)

ASV_Fungi$sample_id<-row.names(ASV_Fungi)


ASV_Fungi_Low<-semi_join(ASV_Fungi,meta_fungi_Low, by="sample_id", copy=F )
dim(ASV_Fungi_Low)

ASV_Fungi_Low<-ASV_Fungi_Low[,-7755]

row.names(ASV_Fungi_Low)==row.names(meta_fungi_Low)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_Low)==0)

#renove these columns 
ASV_Fungi_L <- ASV_Fungi_Low[ ,-zeroes]
dim(ASV_Fungi_Low) #94 Samples 7754 ASVs
dim(ASV_Fungi_L) #94 Samples 3717 ASVs

AdonisFungiLow<-adonis2(ASV_Fungi_L~meta_fungi_Low$fire*meta_fungi_Low$depth+meta_fungi_Low$fire+meta_fungi_Low$depth,permutations=9999, method="bray")
AdonisFungiLow
capture.output(AdonisFungiLow, file = "AdonisFungi_Low.csv")

#HIgh

meta_fungi_High<-subset(meta_fungi, severity=="high")

meta_fungi_High$sample_id<-row.names(meta_fungi_High)

ASV_Fungi$sample_id<-row.names(ASV_Fungi)


ASV_Fungi_High<-semi_join(ASV_Fungi,meta_fungi_High, by="sample_id", copy=F )
dim(ASV_Fungi_High)

ASV_Fungi_High<-ASV_Fungi_High[,-7755]

row.names(ASV_Fungi_High)==row.names(meta_fungi_High)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Fungi_High)==0)

#renove these columns 
ASV_Fungi_H <- ASV_Fungi_High[ ,-zeroes]
dim(ASV_Fungi_High) #92 Samples 7754 ASVs
dim(ASV_Fungi_H) #94 Samples 2589 ASVs

AdonisFungiHigh<-adonis2(ASV_Fungi_H~meta_fungi_High$fire*meta_fungi_High$depth+meta_fungi_High$fire+meta_fungi_High$depth,permutations=9999, method="bray")
AdonisFungiHigh
capture.output(AdonisFungiHigh, file = "AdonisFungi_High.csv")

#############R squared values by site##############

library(ggplot2)
library(ggpubr)


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence")

Rsquare<-read.csv("Rsq_plot.csv")
Rsquare

Rsquare<-Rsquare[order(Rsquare$TSF),]

Fig_Rfungi <- ggplot(Rsquare, aes(x=Site, y=R.fungi,  group=R.fungi))+
  theme_bw() + #make black and white
  ylab("R-squared") +  #change yaxis label
  stat_summary(fun=mean,geom="point", size=4)+
  theme(legend.position = "bottom", #put legend to the right of the graph
        legend.title = element_blank(), #remove legend titles
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=12, angle=20,hjust=1), #make y axis text larger and angled
        axis.text.x = element_text(size=12,angle=45, hjust=1)) 

Fig_Rfungi

Fig_Rbacteria <- ggplot(Rsquare, aes(x=Site, y=R.bacteria,  group=R.bacteria))+
  theme_bw() + #make black and white
  ylab("R-squared") +  #change yaxis label
  stat_summary(fun=mean,geom="point", size=4)+
  theme(legend.position = "bottom", #put legend to the right of the graph
        legend.title = element_blank(), #remove legend titles
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=12, angle=20,hjust=1), #make y axis text larger and angled
        axis.text.x = element_text(size=12,angle=45, hjust=1)) 

Fig_Rbacteria

NMDS_R_figures<- ggarrange(NMDS_all,Fig_Rfungi,NMDS_all_bac,Fig_Rbacteria, labels = c("A", "B","C","D"), ncol = 2, nrow = 2, common.legend=TRUE)
NMDS_R_figures

