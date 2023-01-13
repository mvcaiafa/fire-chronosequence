setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Postdoc_stuff/Fire_Chronosequence/Diversity_analysis_B")

############################################################################################
#...................................Richness GLMM..........................................#
############################################################################################
library(lme4) 
library(MASS)
library(nlme)
library(MuMIn)
citation("lme4")

bacterialdata<-read.csv("ColoradoBacteria_richnessandmetadata.csv", row.names=1)

# Faster way to identify the number of levels in your data set
str(bacterialdata)
names(bacterialdata)

bacterialdata$year <- as.numeric(bacterialdata$year)
bacterialdata$S.obs <- as.numeric(bacterialdata$S.obs)

bacterialdata$year<-scale(bacterialdata$year)

#test if data are normal distributed, if significant not normal
shapiro.test(bacterialdata$S.obs)
hist(bacterialdata$S.obs)
#Data not normally distributed!! Use GLMM 

#* * * TREATMENT EFFECTS .......................................................
#Checking mean and variance of the data to see if its NegBinom..................
df <- bacterialdata %>%
  group_by(severity) %>%
  summarise(means = mean(S.obs),
            variance = var(S.obs));df #

##Start with a poisson model ....................................................
Pois1 <- glm(S.obs ~ severity, data = bacterialdata, family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion  ods=62.65752
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

#Zero inflaited model ..........................................................
sum(bacterialdata$S.obs == 0) #0

#Test different identity link to see if disperssion decreases ..................

Pois2 <- glm(S.obs ~ severity , data = bacterialdata,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #ods=62.65752

Pois3 <- glm(S.obs ~ severity , data = bacterialdata,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#ods=62.65752

Gam <- glm(S.obs ~ severity ,data = bacterialdata,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam#0.1482164

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
nb1 <- glm.nb(S.obs ~ severity,data = bacterialdata)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb# 1.03526 ##acceptable
plot(nb1)

AIC(Pois1, Pois2, Pois3, Gam, nb1)#Negative binomial better 

#Test negative model
Sobs_bac_neg1 <- glmer.nb(S.obs ~ (1|plot),data =bacterialdata);summary(Sobs_bac_neg1)

#BACKWARD MODEL SELECTION........................................................................
# * * FULL MODEL ..............................................................................
Full_Sobs_Bac<-glmer.nb( S.obs~fire*year + fire*severity + fire*depth + year*severity + year*depth + 
               severity*depth +fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Full_Sobs_Bac)

Null_Sobs_Bac<-glmer.nb(S.obs ~ (1|plot),data =bacterialdata);summary(Null_Sobs_Bac)

anova(Full_Sobs_Bac,Null_Sobs_Bac)#Full model better, AIC = 3618.2

#Start backwards model selection..............................................................

#Remove fire*year

Mod2_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + year*severity + year*depth + 
                           severity*depth +fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Mod2_Bac)

anova(Full_Sobs_Bac,Mod2_Bac)##Same AIC, not significant, remove

#Remove fire*severity +

Mod3_Bac<-glmer.nb( S.obs~  fire*depth + year*severity + year*depth + 
                      severity*depth +fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Mod3_Bac)

anova(Mod2_Bac,Mod3_Bac)##Mod2_Bac better, significant, DO NOT remove

#Remove fire*depth +
Mod4_Bac<-glmer.nb( S.obs~ fire*severity + year*severity + year*depth + 
                      severity*depth +fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Mod4_Bac)
anova(Mod2_Bac,Mod4_Bac)##Mod2_Bac better, significant, DO NOT remove

#Remove year*severity + 
Mod5_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + year*depth + 
                      severity*depth +fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Mod5_Bac)
anova(Mod2_Bac,Mod5_Bac)##Same AIC, not significant, remove

#Remove year*depth + 
Mod6_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + 
                      severity*depth +fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Mod6_Bac)
anova(Mod5_Bac,Mod6_Bac)##Same AIC, not significant, remove

#Remove severity*depth
Mod7_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + 
                      fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Mod7_Bac)
anova(Mod6_Bac,Mod7_Bac)##Mod6_Bac better, significant, DO NOT remove

#Remove fire
Mod8_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + severity*depth +year+severity+depth+(1|plot), data = bacterialdata);summary(Mod8_Bac)
anova(Mod6_Bac,Mod8_Bac)##Same AIC, not significant, remove

#Remove year
Mod9_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + severity*depth +severity+depth+(1|plot), data = bacterialdata);summary(Mod9_Bac)
anova(Mod8_Bac,Mod9_Bac)##Same AIC, not significant, remove

#Remove severity
Mod10_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + severity*depth +depth+(1|plot), data = bacterialdata);summary(Mod10_Bac)
anova(Mod9_Bac,Mod10_Bac)##Same AIC, not significant, remove

#Remove depth
Mod11_Bac<-glmer.nb( S.obs~ fire*severity + fire*depth + severity*depth + (1|plot), data = bacterialdata);summary(Mod11_Bac)
anova(Mod10_Bac,Mod11_Bac)##Same AIC, not significant, remove

#Compare all models....................................................................
anova(Full_Sobs_Bac, Mod2_Bac, Mod3_Bac, Mod4_Bac, Mod5_Bac, Mod6_Bac, Mod7_Bac, Mod8_Bac, Mod9_Bac,Mod10_Bac,Mod11_Bac)
AICc(Full_Sobs_Bac, Mod2_Bac, Mod3_Bac, Mod4_Bac, Mod5_Bac, Mod6_Bac, Mod7_Bac, Mod8_Bac, Mod9_Bac,Mod10_Bac,Mod11_Bac)

BacRichModel<-glmer.nb( S.obs~ fire*severity + fire*depth + severity*depth + (1|plot), data = bacterialdata)
summary(BacRichModel)
anova(BacRichModel)
drop1(BacRichModel, test= "Chi" )


Anova(BacRichModel)

R2Bac<-r.squaredGLMM(BacRichModel)

#Model Diagnostic & Validation plots ...........................................
par(mfrow=c(1,1));plot(BacRichModel)
Fitted <- fitted(BacRichModel) #command fitted gives already e^model
Res <- resid(BacRichModel, type = "pearson")

plot(x = Fitted,  y = Res,xlab = "Fitted values",ylab = "Pearson residuals")
abline(h = 0, lty = 2, col="red")

drop1(Full_Sobs_Bac, test = "Chi")

drop1(BacRichModel, test = "Chi")

Modx_Bac<-glmer.nb( S.obs~  fire+year+severity+depth+(1|plot), data = bacterialdata);summary(Modx_Bac)




#######
##Percent Change
##
library(dplyr)
Bac_ave <- bacterialdata %>%
  group_by(fire, severity,depth) %>%
  summarise(means = mean(S.obs),
            variance = var(S.obs))

write.csv(Bac_ave, "BacRich_average.csv")


############################################################################################
#..................................Community Permanova..... ...............................#
############################################################################################

library(vegan)
library(ape)

ASV_Bacteria.e5321<-read.csv("ASV_Bacteria.e5321.csv", row.names=1)
metadata_B <- read.csv("metadata_diversity_Bacteria.csv", row.names=1)

#verify both tables are in teh same order
row.names(ASV_Bacteria.e5321)==row.names(metadata_B)

#reorder  alphabetically
ASV_Bacteria.e5321_order<-ASV_Bacteria.e5321[order(row.names(ASV_Bacteria.e5321)),]

row.names(ASV_Bacteria.e5321_order)==row.names(metadata_B)

##ADONIS (Bray-curtis)

str(metadata_B)
metadata_B$year<-as.numeric(meta_fungi$year)

AdonisBacteria<-adonis(ASV_Bacteria.e5321~metadata_B$fire*metadata_B$year + metadata_B$fire*metadata_B$severity + metadata_B$fire*metadata_B$depth + metadata_B$year*metadata_B$severity + metadata_B$year*metadata_B$depth + 
  metadata_B$severity*metadata_B$depth,permutations=9999, method="bray") 
AdonisBacteria

AdonisBacteria2<-adonis(ASV_Bacteria.e5321~metadata_B$severity*metadata_B$fire+metadata_B$severity*metadata_B$depth+metadata_B$year,permutations=9999, method="bray")


##Betadisper

dist_mat_B<-vegdist(ASV_Bacteria.e5321,distance="bray")
beta_severity<-betadisper(dist_mat_B, metadata_B$severity)
anova(beta_severity)

beta_site<-betadisper(dist_mat_B, metadata_B$fire)
anova(beta_site)

beta_depth<-betadisper(dist_mat_B, metadata_B$depth)
anova(beta_depth)

##########################
#####ADONIS by site#######
##########################

###Church Park

meta_bacteria_CP<-subset(metadata_B, fire=="church park")

meta_bacteria_CP$sample_id<-row.names(meta_bacteria_CP)

ASV_Bacteria.e5321_order$sample_id<-row.names(ASV_Bacteria.e5321_order)

ASV_Bacteria_CP<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_CP, by="sample_id", copy=F )
dim(ASV_Bacteria_CP)

ASV_Bacteria_CP<-ASV_Bacteria_CP[,-29602]

row.names(ASV_Bacteria_CP)==row.names(meta_bacteria_CP)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_CP)==0)

#renove these columns 
ASV_Bacteria_ChurchPark <- ASV_Bacteria_CP[ ,-zeroes]
dim(ASV_Bacteria_CP) #53 Samples 29601 ASVs
dim(ASV_Bacteria_ChurchPark) #53 Samples 10840 ASVs

AdonisBacteriaCP<-adonis(ASV_Bacteria_ChurchPark~meta_bacteria_CP$severity*meta_bacteria_CP$depth+meta_bacteria_CP$severity+meta_bacteria_CP$depth,permutations=9999, method="bray")
AdonisBacteriaCP
capture.output(AdonisBacteriaCP, file = "AdonisBacteria_CP.csv")

###Beaver Creek

meta_bacteria_BV<-subset(metadata_B, fire=="beaver creek")

meta_bacteria_BV$sample_id<-row.names(meta_bacteria_BV)

ASV_Bacteria.e5321_order$sample_id<-row.names(ASV_Bacteria.e5321_order)

ASV_Bacteria_BV<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_BV, by="sample_id", copy=F )
dim(ASV_Bacteria_BV)

ASV_Bacteria_BV<-ASV_Bacteria_BV[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_BV)==0)

#renove these columns 
ASV_Bacteria_BeaverCreek <- ASV_Bacteria_BV[ ,-zeroes]
dim(ASV_Bacteria_BV) #49 Samples 7754 ASVs
dim(ASV_Bacteria_BeaverCreek) #49 Samples 8500 ASVs

AdonisBacteriaBV<-adonis(ASV_Bacteria_BeaverCreek~meta_bacteria_BV$severity*meta_bacteria_BV$depth+meta_bacteria_BV$severity+meta_bacteria_BV$depth,permutations=9999, method="bray")
AdonisBacteriaBV
capture.output(AdonisBacteriaBV, file = "AdonisBacteria_BV.csv")

###Badger Creek
meta_bacteria_BC<-subset(metadata_B, fire=="badger creek")

meta_bacteria_BC$sample_id<-row.names(meta_bacteria_BC)


ASV_Bacteria_BC<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_BC, by="sample_id", copy=F )
dim(ASV_Bacteria_BC)

ASV_Bacteria_BC<-ASV_Bacteria_BC[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_BC)==0)

#renove these columns 
ASV_Bacteria_BadgerCreek <- ASV_Bacteria_BC[ ,-zeroes]
dim(ASV_Bacteria_BC) #72 Samples 7754 ASVs
dim(ASV_Bacteria_BadgerCreek) #72 Samples 2403 ASVs

AdonisBacteriaBC<-adonis(ASV_Bacteria_BadgerCreek~meta_bacteria_BC$severity*meta_bacteria_BC$depth+meta_bacteria_BC$severity+meta_bacteria_BC$depth,permutations=9999, method="bray")
AdonisBacteriaBC
capture.output(AdonisBacteriaBC, file = "AdonisBacteria_BC.csv")

###Ryan

meta_bacteria_R<-subset(metadata_B, fire=="ryan")

meta_bacteria_R$sample_id<-row.names(meta_bacteria_R)

ASV_Bacteria_R<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_R, by="sample_id", copy=F )
dim(ASV_Bacteria_R)

ASV_Bacteria_R<-ASV_Bacteria_R[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_R)==0)

#renove these columns 
ASV_Bacteria_Ryan <- ASV_Bacteria_R[ ,-zeroes]
dim(ASV_Bacteria_R) #54 Samples 7754 ASVs
dim(ASV_Bacteria_Ryan) #54 Samples 1818 ASVs

AdonisBacteriaR<-adonis(ASV_Bacteria_Ryan~meta_bacteria_R$severity*meta_bacteria_R$depth+meta_bacteria_R$severity+meta_bacteria_R$depth,permutations=9999, method="bray")
AdonisBacteriaR

capture.output(AdonisBacteriaR, file = "AdonisBacteria_R.csv")


###Mullen

meta_bacteria_M<-subset(metadata_B, fire=="mullen")

meta_bacteria_M$sample_id<-row.names(meta_bacteria_M)

ASV_Bacteria_M<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_M, by="sample_id", copy=F )
dim(ASV_Bacteria_M)

ASV_Bacteria_M<-ASV_Bacteria_M[,-29602]

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_M)==0)

#renove these columns 
ASV_Bacteria_Mullen <- ASV_Bacteria_M[ ,-zeroes]
dim(ASV_Bacteria_M) #53 Samples 7754 ASVs
dim(ASV_Bacteria_Mullen) #54 Samples 2430 ASVs

AdonisBacteriaM<-adonis(ASV_Bacteria_Mullen~meta_bacteria_M$severity*meta_bacteria_M$depth+meta_bacteria_M$severity+meta_bacteria_M$depth,permutations=9999, method="bray")
AdonisBacteriaM

capture.output(AdonisBacteriaM, file = "AdonisBacteria_M.csv")

##############################
#####ADONIS by severity#######
##############################

###Control

meta_bacteria_Control<-subset(metadata_B, severity=="control")

meta_bacteria_Control$sample_id<-row.names(meta_bacteria_Control)

ASV_Bacteria.e5321_order$sample_id<-row.names(ASV_Bacteria.e5321_order)


ASV_Bacteria_Control<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_Control, by="sample_id", copy=F )
dim(ASV_Bacteria_Control)

ASV_Bacteria_Control<-ASV_Bacteria_Control[,-29602]

row.names(ASV_Bacteria_Control)==row.names(meta_bacteria_Control)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_Control)==0)

#renove these columns 
ASV_Bacteria_C <- ASV_Bacteria_Control[ ,-zeroes]
dim(ASV_Bacteria_Control) #95 Samples 29601 ASVs
dim(ASV_Bacteria_C) #95 Samples 13237 ASVs

AdonisBacteriaControl<-adonis2(ASV_Bacteria_C~meta_bacteria_Control$fire*meta_bacteria_Control$depth+meta_bacteria_Control$fire+meta_bacteria_Control$depth,permutations=9999, method="bray")
AdonisBacteriaControl
capture.output(AdonisBacteriaControl, file = "AdonisBacteria_Control.csv")

#Low

meta_bacteria_Low<-subset(metadata_B, severity=="low")

meta_bacteria_Low$sample_id<-row.names(meta_bacteria_Low)

ASV_Bacteria.e5321_order$sample_id<-row.names(ASV_Bacteria.e5321_order)


ASV_Bacteria_Low<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_Low, by="sample_id", copy=F )
dim(ASV_Bacteria_Low)

ASV_Bacteria_Low<-ASV_Bacteria_Low[,-29602]

row.names(ASV_Bacteria_Low)==row.names(meta_bacteria_Low)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_Low)==0)

#renove these columns 
ASV_Bacteria_L <- ASV_Bacteria_Low[ ,-zeroes]
dim(ASV_Bacteria_Low) #93 Samples 29601 ASVs
dim(ASV_Bacteria_L) #95 Samples 15527 ASVs

AdonisBacteriaLow<-adonis2(ASV_Bacteria_L~meta_bacteria_Low$fire*meta_bacteria_Low$depth+meta_bacteria_Low$fire+meta_bacteria_Low$depth,permutations=9999, method="bray")
AdonisBacteriaLow
capture.output(AdonisBacteriaLow, file = "AdonisBacteria_Low.csv")

#HIgh

meta_bacteria_High<-subset(metadata_B, severity=="high")

meta_bacteria_High$sample_id<-row.names(meta_bacteria_High)

ASV_Bacteria.e5321_order$sample_id<-row.names(ASV_Bacteria.e5321_order)


ASV_Bacteria_High<-semi_join(ASV_Bacteria.e5321_order,meta_bacteria_High, by="sample_id", copy=F )
dim(ASV_Bacteria_High)

ASV_Bacteria_High<-ASV_Bacteria_High[,-29602]

row.names(ASV_Bacteria_High)==row.names(meta_bacteria_High)

#remove ASVs that got zero reads after deleting low severity
zeroes <- which(colSums(ASV_Bacteria_High)==0)

#renove these columns 
ASV_Bacteria_H <- ASV_Bacteria_High[ ,-zeroes]
dim(ASV_Bacteria_High) #92 Samples 29601 ASVs
dim(ASV_Bacteria_H) #92 Samples 12695 ASVs

AdonisBacteriaHigh<-adonis2(ASV_Bacteria_H~meta_bacteria_High$fire*meta_bacteria_High$depth+meta_bacteria_High$fire+meta_bacteria_High$depth,permutations=9999, method="bray")
AdonisBacteriaHigh
capture.output(AdonisBacteriaHigh, file = "AdonisBacteria_High.csv")
