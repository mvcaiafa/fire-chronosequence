###Log transform Copy number. Reduces overdispersion
Fungi$LogCopy<-log(Fungi$CopyNumber)
shapiro.test(Fungi$LogCopy)

#* * * TREATMENT EFFECTS .......................................................
#Checking mean and variance of the data to see if its NegBinom..................
df <- Fungi %>%
  group_by(Fire) %>%
  summarise(means = mean(CopyNumber),
            variance = var(CopyNumber));df #suggest negative binomial

##Start with a poisson model ....................................................
Pois1 <- glm(LogCopy ~ Fire, data = Fungi, family = poisson(link = log))
ods <- Pois1$deviance / Pois1$df.residual;ods#0verdispersion  ods=116938.3
par(mfrow=c(2,2))
plot(Pois1)#look at plots, cooks, levenes etc

#Zero inflaited model ..........................................................
sum(Fungi$LogCopy == 0) #0

#Test different identity link to see if disperssion decreases ..................

Pois2 <- glm(LogCopy ~ Fire , data = Fungi,family = poisson(link = sqrt))
ods2<- Pois2$deviance / Pois2$df.residual;ods2 #

Pois3 <- glm(LogCopy ~ Fire , data = Fungi,family = poisson(link = identity))
ods3<- Pois3$deviance / Pois3$df.residual;ods3#

Gam <- glm(LogCopy ~ Fire,data = Fungi,family = "Gamma")
odsGam<-Gam$deviance/Gam$df.residual;odsGam# 0.04150967

#Negative binomial... ..........................................................
#* mean and variance suggest this is the model to use
library(MASS)
nb1 <- glm.nb(LogCopy ~ Fire,data = Fungi)
ods_nb <- (nb1$deviance / nb1$df.residual);ods_nb#0.3761
plot(nb1)

#Compare models.................................................................
AIC(Pois1, Pois2, Pois3, Gam, nb1)#gamma better 


##Re-scale variable

Fungi$years <- scale(Fungi$Year)

#Explore Temporal correlations to ensure that we are capturing as much of the variance............... 
#Test the null models to see what level of nestedness we should use
qpcr_neg1 <- glmer(LogCopy ~ (1|Plot), family= "Gamma",data =Fungi);summary(qpcr_neg1)
qpcr_neg2 <- glmer(LogCopy ~ (1|Plot)+(1|Replicate),family= "Gamma", data = Fungi);summary(qpcr_neg2)

#MuMIN for model selection..................................................
#lower AIC the better......................................
AIC(qpcr_neg1,qpcr_neg2)#qpcr_neg1 better, AIC=6739.998

#BACKWARD MODEL SELECTION........................................................................
# * * FULL MODEL ..............................................................................
Full<-glmer( LogCopy~Fire*Year + Fire*Severity + Fire*Depth + Year*Severity + Year*Depth + 
                  Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Full)



Null<-glmer( CopyNumber~ (1|Plot),family= "Gamma", data =Fungi);summary(Null)
anova(Full,Null)#Full model better, AIC = 1115.7

#Start backwards model selection..............................................................

#Remove Fire*Year

Mod2<- glmer(LogCopy~ Fire*Severity + Fire*Depth + Year*Severity + Year*Depth + 
                Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod2)

anova(Full, Mod2)##Same AIC, not significant, remove

#Remove Fire*Severity

Mod3<- glmer(LogCopy~  Fire*Depth + Year*Severity + Year*Depth + 
               Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod3)

anova(Mod2,Mod3)##Mod3 better, not significant, remove

# Remove Fire*Depth

Mod4<- glmer(LogCopy~ Year*Severity + Year*Depth + Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod4)


anova(Mod3,Mod4)## Same AIC, not significant, remove

#  Remove Year*Severity

Mod5<- glmer(LogCopy~ Year*Depth + Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod5)

anova(Mod4,Mod5)## Mod5 better, not significant, remove

# Remove Year*Depth

Mod6<- glmer(LogCopy~ Severity*Depth +Fire+Year+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod6)

anova(Mod5,Mod6)##Mod6, not significant, remove 

# Remove Severity*Depth 

Mod7<- glmer(LogCopy~ Fire+Year+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod7)

anova(Mod6, Mod7)##Mod6 better, significant, do not remove


# Remove Year

Mod8<- glmer(LogCopy~ Fire+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod8)
anova(Mod8, Mod7)##Same AIC, not significant, remove


# Remove Fire
Mod9<-  glmer(LogCopy~ Severity+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod9)
anova(Mod8,Mod9)##Mod8 better, significant, DO NOT remove

### Remove Severity
Mod10<-  glmer(LogCopy~ Fire+Depth+(1|Plot), family= "Gamma", data = Fungi);summary(Mod10)
anova(Mod8,Mod10)##Mod8 better, significant, DO NOT remove

### Remove Depth
Mod11<- glmer(LogCopy~ Fire+Severity+(1|Plot), family= "Gamma", data = Fungi);summary(Mod11)
anova(Mod8,Mod11)##Mod8 better, significant, DO NOT remove

#Compare all models....................................................................
anova(Full, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10,Mod11)
AICc(Full, Mod2, Mod3, Mod4, Mod5, Mod6, Mod7, Mod8, Mod9,Mod10,Mod11)


FunModel<-glmer(LogCopy~ Fire+Severity+Depth+(1|Plot), family= "Gamma", data = Fungi)
summary(FunModel)
anova(FunModel)

drop1(FunModel, test = "Chi")

library(r2glmm)

R2FunModel<- r2beta(model=FunModel,partial=TRUE,method='sgv')
plot(R2FunModel)
Anova(FunModel)


