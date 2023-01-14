
#set working directory
setwd("~/Desktop/qPCR")


#load libraries
library("tidyverse")
library("scales")
library("gridExtra")
library("ggpubr")

#load in data
qpcr <- read.csv ("Colorado_chronosequences.csv")

# Faster way to identify the number of levels in your data set
str(Bacteria)


names(qpcr)

qpcr$Kingdom <- as.factor(qpcr$Kingdom)
qpcr$Fire <- as.factor(qpcr$Fire)

qpcr$Severity <- as.factor(qpcr$Severity)
qpcr$Severity <- try(relevel(qpcr$Severity,"Control")) #change base level


qpcr$Depth <- as.factor(qpcr$Depth)
qpcr$Plot <- as.factor(qpcr$Plot)
qpcr$Year<- as.numeric(qpcr$Year)

qpcr$Year<-scale(qpcr$Year)
qpcr$Replicate<- as.factor(qpcr$Replicate)
# Testing if data set is normnalized or not (Large p-values, Normal - anova, Small p-values, Non-normal = Kruskal)
shapiro.test(qpcr$CopyNumber)




#########***************************############
#########.........STATS.............############
################################################

library(lme4)
library(MASS)
library(nlme)


Bacteria<- qpcr[which(qpcr$Kingdom=="Bacteria"),]
hist(Bacteria$CopyNumber,breaks = 100)
shapiro.test(Bacteria$CopyNumber)#Check normality

###Log transform Copy number. Reduces overdispersion
Bacteria$LogCopy<-log(Bacteria$CopyNumber)
shapiro.test(Bacteria$LogCopy)

attach(Bacteria)

Bac_Biomass <- Bacteria %>%
  group_by(Fire,Severity,Depth) %>%
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber))

write.csv(Bac_Biomass, "Bacbio_average.csv")

#* * * TREATMENT EFFECTS .......................................................
#Checking mean and variance of the data to see if its NegBinom..................
df <- Bacteria %>%
  group_by(Fire) %>%
  summarise(means = mean(LogCopy),
            variance = var(LogCopy));df #

##Start with a poisson model ....................................................
Pois1_B <- glm(LogCopy ~ Fire, data = Bacteria, family = poisson(link = log))
ods <- Pois1_B$deviance / Pois1_B$df.residual;ods#0verdispersion  ods= 0.2607024
par(mfrow=c(2,2))
plot(Pois1_B)#look at plots, cooks, levenes etc

#Zero inflaited model ..........................................................
sum(Bacteria$LogCopy == 0) #0

#Test different identity link to see if disperssion decreases ..................

Pois2_B <- glm(LogCopy ~ Fire , data = Bacteria,family = poisson(link = sqrt))
ods2<- Pois2_B$deviance / Pois2_B$df.residual;ods2 #0.2607024

Pois3_B <- glm(LogCopy ~ Fire , data = Bacteria,family = poisson(link = identity))
ods3<- Pois3_B$deviance / Pois3_B$df.residual;ods3#0.2607024

Gam_B <- glm(LogCopy ~ Fire,data = Bacteria,family = "Gamma")
odsGam_B<-Gam_B$deviance/Gam_B$df.residual;odsGam_B#0.02961916

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
library(MASS)
nb1_B <- glm.nb(LogCopy ~ Fire,data = Bacteria)
ods_nb1_B <- (nb1_B$deviance / nb1_B$df.residual);ods_nb1_B#0.2606998
plot(nb1_B)

#Compare models.................................................................
AIC(Pois1_B, Pois2_B, Pois3_B, Gam_B, nb1_B)#Gamma better


#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
#Test the null models to see what level of nestedness we should use
qpcr_neg1_B <- glmer(LogCopy ~ (1|Plot), family = "Gamma",data =Bacteria);summary(qpcr_neg1_B)
qpcr_neg2_B <- glmer(LogCopy ~ (1|Plot)+(1|Replicate), family = "Gamma", data = Bacteria);summary(qpcr_neg2_B)

#MuMIN for model selection..................................................
#lower AIC the better......................................
AIC(qpcr_neg1_B,qpcr_neg2_B)#qpcr_neg1_B better, AIC=6739.998

#BACKWARD MODEL SELECTION........................................................................
# * * FULL MODEL ..............................................................................
Full_B<-glmer(LogCopy~Fire*Year + Fire*Severity + Fire*Depth + Year*Severity + Year*Depth + 
                  Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Full_B)


Null_B<-glmer(LogCopy ~ (1|Plot), family = "Gamma",data =Bacteria);summary(Null_B)
anova(Full_B,Null_B)#Null model better, AIC = 1209

#Start backwards model selection..............................................................

#Remove Fire*Year

Mod2_B<- glmer(LogCopy~ Fire*Severity + Fire*Depth + Year*Severity + Year*Depth + 
               Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod2_B)

anova(Full_B, Mod2_B)##Same AIC, not significant, remove

#Remove Fire*Severity

Mod3_B<- glmer(LogCopy~  Fire*Depth + Year*Severity + Year*Depth + 
               Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod3_B)

anova(Mod2_B,Mod3_B)##Mod3_B better, not significant, remove

# Remove Fire*Depth

Mod4_B<- glmer(LogCopy~ Year*Severity + Year*Depth + 
               Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod4_B)


anova(Mod3_B,Mod4_B)## Mod4_B better, ot significant, remove

#  Remove Year*Severity

Mod5_B<- glmer(LogCopy~  Year*Depth + Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod5_B)

anova(Mod4_B,Mod5_B)## Mod5_B better, not significant, remove

# Remove Year*Depth

Mod6_B<- glmer(LogCopy~ Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod6_B)

anova(Mod5_B,Mod6_B)##Mod6_B better, not significant, remove 

# Remove Severity*Depth 

Mod7_B<- glmer(LogCopy~ Fire+Year+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod7_B)

anova(Mod6_B, Mod7_B)##Mod7_B better, not significant, remove 


# Remove Year

Mod8_B<- glmer(LogCopy~ Fire+Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod8_B)

anova(Mod7_B, Mod8_B)##Same AIC, not significant, remove

#Remove Fire

Mod9_B<- glmer(LogCopy~Severity+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod9_B)

anova(Mod8_B, Mod9_B)##Mod8_B better, signifficant, DO NOT REMOVE 

#Remove Severity

Mod10_B<- glmer(LogCopy~ Fire+Depth+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod10_B)

anova(Mod8_B,Mod10_B)#Mod8_B better, not significant, remove 

#Remove Depth

Mod11_B<- glmer(LogCopy~ Fire+(1|Plot), family = "Gamma",data = Bacteria);summary(Mod11_B)

anova(Mod11_B,Mod8_B)#Mod11_B better, not significant, remove 

Mod12_B<- glmer(LogCopy~ Fire+Severity+(1|Plot), family = "Gamma",data = Bacteria)

#Compare all models....................................................................
anova(Full_B, Null_B, Mod2_B, Mod6_B, Mod7_B, Mod8_B, Mod9_B,Mod10_B,Mod11_B, Mod12_B)
AICc(Full_B, Null_B, Mod2_B, Mod3_B, Mod4_B, Mod5_B, Mod6_B, Mod7_B, Mod8_B, Mod9_B,Mod10_B,Mod11_B)



BacBiomass_model<- glmer(LogCopy~ Fire+Severity+(1|Plot), family = "Gamma",data = Bacteria)

drop1(BacBiomass_model, test="Chi")
