#clear workspace
rm(list = ls())

#load packages
library(data.table)
library(ggplot2)
library(DescTools) #for CCC
library(lubridate)
library(lme4)
library(DHARMa)
library(effects)
library(hms)
library(parameters)
library(igraph)
library(emmeans)
library(MuMIn) #R^2 calculations
library(glmmTMB)
library(broom.mixed)

#load functions
source("functions.R")
source("prepareTracking.R")

##### Loading and preparing Data ###########
#load tracking data
trackingData <- fread("filtered_30s_full_combined.csv")
#load social data
observData = fread('socialData.csv')
#load comb sizes
combData = fread('focalSelection.csv')
combData[,ID := paste0(Pen,ID)]
combData[, Date := as_date(Date, format = "%d.%m.%Y")]
combData = combData[Date < as_date("2020-01-01")]
combData = combData[, .SD, .SDcols = c("Pen", "ID", "Date", "HenID", "TagID", "Mean5_7")]
colnames(combData)[6] <- "Comb" 

#prepare WoA table 
tableWoA = data.table(Date = seq(ymd("2019-10-11"), ymd("2020-07-02"), by = "day"),
                      WoA = rep(18:55, each = 7))
#load health data
healthData = fread('HA_all.csv')
healthData[,ID := paste0(pen,backpack)]
#healthData = healthData[(date == "28.10.2019"| date == "09.12.2019" | date == "29.06.2020" | date == "30.06.2020")
#                        & ID %in% unique(observData$ID),]
healthData = healthData[ID %in% unique(observData$ID),]
healthData[, feathers := neck + wings +tail + cloaca+ breast]
healthData[, footproblems:= r_podo + r_bumble +r_injure + l_podo + l_bumble + l_injure]
healthData[, date := as_date(date, format = "%d.%m.%Y")]
healthData = healthData[tableWoA, on = c(date = "Date"), nomatch = NULL]
healthDataWide = dcast(healthData, formula = ID ~ WoA, value.var = list("weight", "feathers", "wounds","comb", "bare_head", "footproblems"))

KBF = fread("KBF_scores.csv")
#KBF = KBF[(Date == "28.10.2019" |Date == "09.12.2019" | Date == "29.06.2020" | Date == "30.06.2020")
#                        & HenID %in% unique(observData$ID),]
KBF = KBF[HenID %in% unique(observData$ID),]
KBF[, Date := as_date(Date, format = "%d.%m.%Y")]
KBF = KBF[tableWoA, on = "Date", nomatch = NULL]
KBFWide = dcast(KBF, formula = HenID ~ WoA, value.var = "Severity")
colnames(KBFWide)[2:length(colnames(KBFWide))] = paste("KBF", colnames(KBFWide)[2:length(colnames(KBFWide))], sep = "_")

healthData = KBF[, .(WoA, HenID, Severity)][healthData[, .(WoA, ID, weight, feathers, wounds, comb, bare_head, footproblems)], on = c(HenID = "ID",WoA = "WoA")]
healthDataWide = KBFWide[healthDataWide, on = c(HenID = "ID")]

#delete invalid rows
observData = observData[!(Exclusion == 1),]

##### Social Data #####
###### Observer reliability ##### 

reliability = observData[Reliability == 1]

Observer1 = c(reliability$Sum_Sub[reliability$Observer == 'Masha'], 
              reliability$Sum_Dom[reliability$Observer == 'Masha'])
Observer2 = c(reliability$Sum_Sub[reliability$Observer == 'Tatiana'], 
              reliability$Sum_Dom[reliability$Observer == 'Tatiana'])

res = CCC(Observer1, Observer2, ci = "z-transform", conf.level = 0.95) 
#result: CCC = 0.89 [0.81, 0.94]

# lab <- paste("CCC: ", round(res$rho.c[,1], digits = 2), " (95% CI ", 
#              round(res$rho.c[,2], digits = 2), " - ",
#              round(res$rho.c[,3], digits = 2), ")", sep = "")
# z <- lm(Observer1 ~ Observer2)
# 
# par(pty = "s")
# plot(Observer1, Observer2, xlab = "Observer 1", 
#      ylab = "Observer 2", pch = 16)
# abline(a = 0, b = 1, lty = 2)
# abline(z, lty = 1)
# legend(x = "topleft", legend = c("Line of perfect concordance", 
#                                  "Reduced major axis"), lty = c(2,1), lwd = c(1,1), bty = "n")
# text(x = 1.55, y = 3.8, labels = lab)
# 
# tmp.mean <- mean(res$blalt$delta)
# tmp.sd <- sqrt(var(res$blalt$delta))
# 
# plot(res$blalt$mean, res$blalt$delta, pch = 16, 
#      xlab = "Average PEFR ", 
#      ylab = "Difference in PEFR") 
# abline(h = tmp.mean, lty = 1, col = "gray")
# abline(h = tmp.mean - (2 * tmp.sd), lty = 2, col = "gray")
# abline(h = tmp.mean + (2 * tmp.sd), lty = 2, col = "gray")
# legend(x = "topleft", legend = c("Mean difference", 
#                                  "Mean difference +/ 2SD"), lty = c(1,2), bty = "n")
# legend(x = 0, y = 125, legend = c("Difference"), pch = 16, 
#        bty = "n")


###### Social Data prep & Index #####
socialData = observData[!(Observer == 'Masha' & Reliability == 1),]

#Descriptives
socialData[, ratioInteractions := Sum_Actions/Hen_in_sight*60]
#mean number of interactions per bird per minute
socialData[, .(MeanInter = mean(ratioInteractions), sdInter = sd(ratioInteractions),
               MedianInter = median(ratioInteractions), IQRInter = IQR(ratioInteractions))] 

#aggregate by individual
dataIndex = socialData[, .(Dom=sum(Sum_Dom), Sub=sum(Sum_Sub), Sum_Actions = sum(Sum_Actions),
                           Affil_given = sum(Affiliative_given),
                           Affil_rec = sum(Affiliative_received)),  by = .(ID, Pen)]

#median sum of interactions per bird per observation
dataIndex[, .(MedianSum = mean(Sum_Actions), sdInter = sd(Sum_Actions), min = min(Sum_Actions), max = max(Sum_Actions))] 

dataIndex[, Ratio := Dom/(Dom+Sub)]

dataIndex = dataIndex[order(Ratio),]

socialData = combData[dataIndex, on = "ID"]
socialData[, Hen := paste0("Hen_", HenID)]

socialData[,ID := factor(ID, levels = ID)]

ggplot(data = socialData, aes(x = ID, y = Ratio, colour = as.factor(Pen)))+ 
  geom_point(size = 6)+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  theme_classic(base_size = 18)


henData = socialData[healthDataWide, on = c(ID = "HenID")]
henDataLong = socialData[healthData, on = c(ID = "HenID")]

henData[, Pen := as.factor(Pen)]

rm(dataIndex)
###### Comb Size #######

#comb size by dominance
ggplot(data = henData, aes(x = Comb, y = Ratio))+ 
  # geom_violin()+
  geom_point(aes(colour = as.factor(Pen)),size = 6)+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  facet_grid(~Pen)
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'Comb Size', y = "Dominance ratio")

#weight at 26 by dominance
ggplot(data = henData, aes(x = weight_20, y = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  geom_point(aes(x = weight_55), size = 3, colour = "grey")+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  geom_smooth(aes(x = weight_55), method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'Weight', y = "Dominance ratio")

#weight at 55 by comb size
ggplot(data = henData, aes(x = Comb, y = weight_55))+ 
  # geom_violin()+
  geom_point(size = 6, aes(colour = Ratio))+
  #geom_hline(yintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'Comb Size', y = "Weight")

cor.test(henData$Comb, henData$weight_20)
# correlated -> R = 0.4, p = 0.015
cor.test(henData$Comb, henData$weight_26)
# correlated -> R = 0.55, p <0.001
cor.test(henData$Comb, henData$weight_55)
# correlated -> R = 0.54, p <0.001
cor.test(henData$weight_55, henData$weight_20)
# correlated -> R = 0.5, p = 0.001
cor.test(henData$weight_55, henData$weight_26)
# correlated -> R = 0.82, p < 0.001

model.Comb = lmer(Ratio ~ Comb +(1|Pen), data = henData)
null.Comb =  lmer(Ratio ~ 1 +(1|Pen), data = henData)
resid.Comb = simulateResiduals(model.Comb, 1000)
plot(resid.Comb)
plotResiduals(resid.Comb, form = henData$Comb)
summary(model.Comb)
anova(model.Comb, null.Comb)
parameters(model.Comb)
r.squaredGLMM(model.Comb, null.Comb)

henData[, PredictRatio:= predict(model.Comb)]
ggplot(data = henData, aes(y = Ratio, x = Comb, ))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  geom_line(aes(y = PredictRatio, colour = as.factor(Pen)), size = 1.5)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

###### Weight by dominance ####
#descriptives
henDataLong[, mean(weight), by = WoA]
henDataLong[Ratio < 0.17 , max(weight), by = ID][,mean(V1)]
henDataLong[Ratio > 0.7 , max(weight), by = ID][,mean(V1)]

#weight gain
ggplot(data = henDataLong, aes(x = WoA, y = weight))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio), size = 3)+
  geom_line(aes(group = ID, colour = Ratio))+
  geom_line(data = henDataLong[, .(meanW = mean(weight)), by = WoA], aes(y = meanW), size = 1.5, "-")+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'WoA', y = "Weight")+
  scale_color_gradient(low = "blue", high = "gold")

henData[, gain := weight_55-weight_20]

#gain by dominance
ggplot(data = henData, aes(y = gain, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'weight gain', y = 'Weight gain', x = "Dominance ratio")


hist(henDataLong$weight)#near-normal with small negative skew
model.Weight = lmer(weight ~ WoA*Ratio+(1|Pen), data = henDataLong)
null.Weight =  lmer(weight ~ 1 +(1|Pen), data = henDataLong)
resid.Weight = simulateResiduals(model.Weight, 1000)
plot(resid.Weight) #looks like polynomial effect
plotResiduals(resid.Weight, form = henDataLong$Ratio)
plotResiduals(resid.Weight, form = henDataLong$WoA) #polynomial effect

model.Weight = lmer(weight ~ poly(WoA, 2)+ WoA:Ratio + Ratio +(1|Pen), data = henDataLong)
resid.Weight = simulateResiduals(model.Weight, 1000)
plot(resid.Weight) #good
plotResiduals(resid.Weight, form = henDataLong$Ratio)
plotResiduals(resid.Weight, form = henDataLong$WoA)
summary(model.Weight)
anova(model.Weight, null.Weight)
parameters(model.Weight)
plot(allEffects(model.Weight))
#variance explained by fixed factors and entire model
r.squaredGLMM(model.Weight, null.Weight)

#TODO: starting at which WoA is does a weight difference appear?
trend.Weight = emtrends(model.Weight, "WoA", va = "Ratio", at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))
summary(trend.Weight)
test.Weight = lmer(weight ~ poly(WoA, 2)+ WoA:Ratio + Ratio +(1|Pen), data = henDataLong[WoA <51,])


plotData = as.data.table(emmeans(model.Weight, ~ pairwise ~ WoA*Ratio, at =  list(Ratio = round(quantile(henDataLong$Ratio), digits = 2),
                                                         WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))$emmeans)

henDataLong[, PredictWeight:= predict(model.Weight)]
henDataLong[, RatioSplit := ifelse(Ratio >0.5,"Dom" , "Sub")]

ggplot()+ 
  geom_point(data = henDataLong, aes(x = WoA, y = weight), size = 3)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = lower.CL, ymax = upper.CL, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)


#Gain descriptive

henData[, gain := weight_55 - weight_20, by = HenID]
henData[Ratio < 0.17 , mean(gain)]
henData[Ratio > 0.7 , mean(gain)]



###### Health & social data ####

henDataLong[, mean(Severity), by = WoA]

#KBF and dominance
ggplot(data = henData, aes(y = KBF_55, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Health & Dominance', y = 'KBF Severity', x = "Dominance ratio")

#KBF development
ggplot(data = henDataLong, aes(x = WoA, y = Severity, colour = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs(x = 'WoA', y = "KBF Severity")+
  scale_color_gradient(low = "blue", high = "gold")

#two lowest ranking hens with really low severity -> not laying???

hist(henDataLong$Severity) #normal but with a lot of zeros
model.KBF = lmer(Severity ~ Ratio*WoA +(1|Pen), data = henDataLong)
null.KBF =  lmer(Severity ~ 1 +(1|Pen), data = henDataLong)
red.KBF =  lmer(Severity ~ WoA +(1|Pen), data = henDataLong)
resid.KBF = simulateResiduals(model.KBF, 1000)
plot(resid.KBF) #problem in left corner
plotResiduals(resid.KBF, form = henDataLong[, Ratio]) #okay
plotResiduals(resid.KBF, form = henDataLong[, WoA]) #zero-inflation?
testZeroInflation(resid.KBF) #zero-inflation problem
hist(henDataLong$Severity[henDataLong$WoA >22]) #still many
hist(henDataLong$Severity[henDataLong$WoA >26]) #maybe possible

#to compare estimates direction
summary(model.KBF) #estimate direction: Ratio positive, WoA positive, interaction: negative
#tested by running model without 20 WoA 
test.KBF = lmer(Severity ~ Ratio*WoA +(1|Pen), data = henDataLong[WoA > 26,])
resid.KBF = simulateResiduals(test.KBF, 1000)
plot(resid.KBF) #better
plotResiduals(resid.KBF, form = henDataLong[WoA >26, Ratio]) #okay
plotResiduals(resid.KBF, form = henDataLong[WoA >26, WoA]) #good
summary(test.KBF) # estimate direction the same -> stick to full model

anova(model.KBF, null.KBF)
AIC(model.KBF, red.KBF)
parameters(red.KBF)
r.squaredGLMM(red.KBF)

henDataLong[, mean(feathers), by = WoA]
henDataLong[, diffFeath := feathers - shift(feathers), by = ID][, mean(diffFeath), by = WoA]

#feather cover
ggplot(data = henData, aes(y = feathers_55, x = Ratio, colour = Pen))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, se = F)+
  theme_classic(base_size = 18)+
  labs(y = 'Feather cover', x = "Dominance ratio")

#feather development
ggplot(data = henDataLong, aes(x = WoA, y = feathers,colour = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs( x = 'WoA', y = "Feather loss")+
  scale_color_gradient(low = "blue", high = "gold")

hist(henDataLong$feathers) 
model.Feathers = glmer(feathers ~ Ratio*WoA+(1|Pen), data = henDataLong, family = poisson)
resid.Feathers = simulateResiduals(model.Feathers, 1000)
plot(resid.Feathers) 
testOverdispersion(resid.Feathers)#overdispersed -> negative binomial
model.Feathers = glmmTMB(feathers ~ Ratio + Ratio:WoA+ WoA +(1|Pen), data = henDataLong, family = "nbinom2") 
resid.Feathers = simulateResiduals(model.Feathers, 1000) #looks like polynomial effect
plot(resid.Feathers)
plotResiduals(resid.Feathers, form = henDataLong$Ratio) #good
plotResiduals(resid.Feathers, form = henDataLong$WoA) #polynomial

model.Feathers = glmmTMB(feathers ~ Ratio + Ratio:WoA+ poly(WoA,2) +(1|Pen), data = henDataLong, family = "nbinom2") 
resid.Feathers = simulateResiduals(model.Feathers, 1000) #looks like polynomial effect
plot(resid.Feathers) #okay
plotResiduals(resid.Feathers, form = henDataLong$Ratio) #good
plotResiduals(resid.Feathers, form = henDataLong$WoA)#good
null.Feathers =  glmmTMB(feathers ~ 1 +(1|Pen), data = henDataLong, family = "nbinom2")
red.Feathers = glmmTMB(feathers ~ poly(WoA,2) +(1|Pen), data = henDataLong, family = "nbinom2") 
summary(model.Feathers)
AIC(model.Feathers, null.Feathers)
AIC(model.Feathers, red.Feathers)
plot(allEffects(model.Feathers))
tidy(model.Feathers, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
r.squaredGLMM(red.Feathers, null.Feathers)


#foot problems
median(henDataLong$footproblems)
ggplot(data = henDataLong, aes(x = WoA, y = footproblems,colour = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs( x = 'WoA', y = "Foot problems")+
  scale_color_gradient(low = "blue", high = "gold")

#injuries
median(henDataLong$wounds)
ggplot(data = henDataLong, aes(x = WoA, y = wounds,colour = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs( x = 'WoA', y = "Feather loss")+
  scale_color_gradient(low = "blue", high = "gold")

#comb pecks
median(henDataLong$comb)
ggplot(data = henDataLong, aes(x = WoA, y = comb,colour = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs( x = 'WoA', y = "Feather loss")+
  scale_color_gradient(low = "blue", high = "gold")


#health pca
health.pca <- prcomp(henData[,.(feathers_55, footproblems_55, 
                                KBF_55, weight_55, comb_55, wounds_55, Comb)], scale = TRUE, center = T)
summary(health.pca)
library(ggbiplot)
henData[, Tend := ifelse(Ratio >0.3, ifelse(Ratio >0.54, "Dom", "Mid"), "Sub")]
ggbiplot(health.pca, labels=henData[,HenID], groups = henData[, Tend], ellipse=TRUE)

health.pca$rotation
###### Gentle beak pecking behaviour##########

ggplot(data = henData, aes(y = Affil_given, x = Ratio))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 2)+
  geom_smooth(method = lm)+
  theme_classic(base_size = 18)

ggplot(data = henData, aes(y = Affil_rec, x = Ratio))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 2)+
  geom_smooth(method = lm)+
  theme_classic(base_size = 18)

##### Tracking data ####

#relevant hens
hens = sort(unique(socialData$HenID))

trackingData = prepareTrackingData(trackingData, hens)
#include Week of age
trackingData = trackingData[tableWoA, on = "Date", nomatch = NULL]

#fwrite(trackingData, "trackingData.csv", sep = ";")
trackingData <- fread("trackingData.csv")
rm(trackingData)

##### Sequence similarity ####

#sequence similarity calculated by how many seconds two timelines agree and how many they disagree
#calculated within birds -> comparing consecutive (or further away?) days
#only look at day hours
#Problem: need to deal with time shift for within? I think not will resolve itself  

###### Between individuals ####
betweenIndividuals = similarityBetween(trackingData[Date < as.Date("2020-03-15"),])
betweenIndividuals2 = similarityBetween(trackingData[Date > as.Date("2020-03-14"),])

meanData = data.table(Date = unique(trackingData$Date),
                      Mean = c(unlist(lapply(betweenIndividuals, function(x){ mean(x, na.rm = TRUE)})),
                               unlist(lapply(betweenIndividuals2, function(x){ mean(x, na.rm = TRUE)}))),
                      Median = c(unlist(lapply(betweenIndividuals, function(x){ median(x, na.rm = TRUE)})),
                               unlist(lapply(betweenIndividuals2, function(x){ median(x, na.rm = TRUE)}))),
                      Max = c(unlist(lapply(betweenIndividuals, function(x){ max(x, na.rm = TRUE)})),
                              unlist(lapply(betweenIndividuals, function(x){ max(x, na.rm = TRUE)}))),
                      SD = c(unlist(lapply(betweenIndividuals, function(x){ sd(x, na.rm = TRUE)})),
                               unlist(lapply(betweenIndividuals2, function(x){ sd(x, na.rm = TRUE)}))),
                      Above = c(unlist(lapply(betweenIndividuals, function(x){ sum(x>0.6, na.rm = TRUE)})),
                              unlist(lapply(betweenIndividuals2, function(x){ sum(x>0.6, na.rm = TRUE)}))),
                      Below = c(unlist(lapply(betweenIndividuals, function(x){ sum(x<0.2, na.rm = TRUE)})),
                                unlist(lapply(betweenIndividuals2, function(x){ sum(x<0.2, na.rm = TRUE)}))),
                      Q1 = c(unlist(lapply(betweenIndividuals, function(x){ quantile(x, 0.25,na.rm = TRUE)})),
                                unlist(lapply(betweenIndividuals2, function(x){ quantile(x, 0.25, na.rm = TRUE)}))),
                      Q3 = c(unlist(lapply(betweenIndividuals, function(x){ quantile(x, 0.75,na.rm = TRUE)})),
                             unlist(lapply(betweenIndividuals2, function(x){ quantile(x, 0.75, na.rm = TRUE)})))
)

ggplot(meanData, aes(x = Date))+
  geom_pointrange(aes(y = Median, ymin = Q1, ymax = Q3))+
  geom_point(aes(y = Mean), color = "red")+
  #geom_point(aes(y = Max))+
  theme_classic(base_size = 18)
ggplot(meanData, aes(x = Date))+
  #geom_point(aes(y = Below))+
  #geom_smooth(aes(y = Below))+
  geom_point(aes(y = Above))+
  geom_smooth(aes(y = Above))+
  theme_classic(base_size = 18)

# similSimply1 = lapply(betweenIndividuals, function(x){ x[x < 0.7|is.na(x)] = 0; x[(x > 0.7)|(x == 0.7)] = 1; x})
# similSimply1 = Reduce('+', similSimply1)
# similSimply2 = lapply(betweenIndividuals2, function(x){ x[x < 0.7|is.na(x)] = 0; x[(x > 0.7)|(x == 0.7)] = 1; x})
# similSimply2 = Reduce('+', similSimply2)
# similSimply = similSimply1 + similSimply2
# 
# similSimply1 = lapply(betweenIndividuals, function(x){ x[is.na(x)] = 0; x})
# similSimply1 = Reduce('+', similSimply1)
# similSimply2 = lapply(betweenIndividuals2, function(x){ x[is.na(x)] = 0; x})
# similSimply2 = Reduce('+', similSimply2)
# similSimply = similSimply1 + similSimply2
# noData1 = lapply(betweenIndividuals, function(x){ x[!is.na(x)] = 0; x[is.na(x)] = 1; x})
# noData1 = Reduce('+', noData1)
# noData2 = lapply(betweenIndividuals2, function(x){ x[!is.na(x)] = 0; x[is.na(x)] = 1; x})
# noData2 = Reduce('+', noData2)
# noData = (noData1 + noData2 - 179)*-1
# 
# similSimply = similSimply/noData
# 
# mean(similSimply, na.rm = TRUE)
# max(similSimply, na.rm = TRUE)
# 
# similSimply[is.na(similSimply)] = 0
# similSimply = similSimply + t(similSimply)
# 
# 
# #[similSimply <0.45]= 0
# 
# g  <- graph.adjacency(similSimply,mode = "undirected", weighted = TRUE, diag = F)
# bc <- edge.betweenness.community(g)
# #par(mfrow=c(1,2))
# plot(as.dendrogram(bc))
# #network vertex names
# V(g)$name
# #inspect network edge and node attributes
# edge_attr(g)
# vertex_attr(g)
# V(g)$domIndex <- socialData$Ratio[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
# V(g)$Pen <- socialData$Pen[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
# plot(g, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
#      vertex.color=c( "pink", "skyblue")[1+(V(g)$domIndex>0.5)] ) 
# assortativity(g, V(g)$domIndex, types1 = graph.strength(g))
# E(g)$width <- E(g)$weight*3
# plot(g, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
#      vertex.color=V(g)$Pen) 
# assortativity(g, V(g)$Pen)
# library(assortnet)
# assortment.continuous(similSimply, V(g)$domIndex, weighted = TRUE, SE = FALSE, M = 1)


#reformat lists into datatables containing cols of HenID1, HenID2, Date, Similarity
sample = lapply(betweenIndividuals, function(x){ x[upper.tri(x) & is.na(x)] = -1; 
                                                  x[lower.tri(x)] = NA; diag(x) = NA; return(x)})
sample = na.omit(rbindlist(lapply(sample, function(x){as.data.table(as.table(x))})))
colnames(sample) = c("Hen1", "Hen2", "Similarity")
sample[Similarity == -1, Similarity := NA] 
nPairs = length(unique(paste(sample$Hen1, sample$Hen2)))
sample[, Date := as.IDate(rep(names(betweenIndividuals), each = nPairs))]

sample2 = lapply(betweenIndividuals2, function(x){ x[upper.tri(x) & is.na(x)] = -1; 
                                                  x[lower.tri(x)] = NA; diag(x) = NA; return(x)})
sample2 = na.omit(rbindlist(lapply(sample2, function(x){as.data.table(as.table(x))})))
colnames(sample2) = c("Hen1", "Hen2", "Similarity")
sample2[Similarity == -1, Similarity := NA] 
nPairs = length(unique(paste(sample2$Hen1, sample2$Hen2)))
sample2[, Date := as.IDate(rep(names(betweenIndividuals2), each = nPairs))]

dataBetween = na.omit(rbind(sample, sample2))

dataBetween[, Hen1_Pen := socialData$Pen[match(Hen1, socialData$Hen)]]
dataBetween[, Hen2_Pen := socialData$Pen[match(Hen2, socialData$Hen)]]
dataBetween[, Hen1_Ratio := socialData$Ratio[match(Hen1, socialData$Hen)]]
dataBetween[, Hen2_Ratio := socialData$Ratio[match(Hen2, socialData$Hen)]]
dataBetween[, PenBool := Hen1_Pen == Hen2_Pen]
dataBetween[, DiffRatio := abs(Hen1_Ratio - Hen2_Ratio)]

dataBetween= dataBetween[tableWoA, on = "Date", nomatch = NULL]

#fwrite(dataBetween, "dataBetween.csv", sep = ";")
dataBetween <- fread("dataBetween.csv")

#build Model
hist(dataBetween$Similarity)

model.Between = lmer(Similarity ~ PenBool+ DiffRatio*WoA + (1|Hen1) + (1|Hen2), data = dataBetween)
null.Between = lmer(Similarity ~ 1+ (1|Hen1) + (1|Hen2), data = dataBetween)
anova(model.Between, null.Between)
resid.Between = simulateResiduals(model.Between, 1000)
plot(resid.Between)
plotResiduals(resid.Between, form = dataBetween$DiffRatio)
plotResiduals(resid.Between, form = dataBetween$PenBool)
plotResiduals(resid.Between, form = dataBetween$WoA)
summary(model.Between)
parameters(model.Between)
plot(allEffects(model.Between))
#variance explained by fixed factors and entire model
r.squaredGLMM(model.Between, null.Between)

dataBetween[, PredictSimilarity := predict(model.Between)]

ggplot(dataBetween[WoA %in% quantile(WoA),], aes(x = DiffRatio, y = Similarity)) +
  geom_point(size=2) + 
  geom_smooth(data = dataBetween[WoA %in% quantile(WoA), mean(PredictSimilarity), by = .(WoA, DiffRatio)], 
            aes(y =V1), formula = y~x, size=1.5, colour = "red") + 
  theme_classic(base_size = 18)+ 
  facet_grid(~WoA)+
  ylab("daily between-individual similarity")
#+ 
#test random Ratio allocation
#TODO: randomise similarity within day to take out day and ratio effect
#TODO: model on difference in ratio or actual ratios? and if actual then how interaction
#TODO: how to deal with pen?

betaTable = data.table(Run = "Orig", t(fixef(model.Between)))

set.seed(42)
for (i in 1:1000){
  cat("Run:",i, "\n")
  dataBetween[, Similarity_rand := sample(Similarity), by = Date]
  model.Between.rand = lmer(Similarity_rand ~ DiffRatio*PenBool+ DiffRatio*WoA + (1|Hen1) + (1|Hen2), data = dataBetween)
  entry = data.table(Run = as.character(i), t(fixef(model.Between.rand)))
  betaTable= rbind(betaTable, entry)
}

ggplot(data = betaTable, aes(x= `DiffRatio:WoA`))+
  geom_histogram()+
  geom_vline(xintercept = betaTable[Run == "Orig", `DiffRatio:WoA`], colour = "red")

p = (sum(betaTable[Run != "Orig", `DiffRatio:WoA`] <= betaTable[Run == "Orig", `DiffRatio:WoA`]) + 1)/dim(betaTable)[1]

#TODO: run again, wrong randomisation...

### Isolate by Nestbox behaviour

btwIndivNest = similarityBetween(trackingData, zone = "Ramp_Nestbox", interval = "morning")

assortNest = data.table(Date = as_date(colnames(btwIndivNest)[-c(1,2, .N)]))
assortNest[, Assort := 0]
# high degree assortativity is a measure of preferential attachment in organizations, 
# where highly connected vertices are connected with each other and a large number 
# of vertices with low degree make up the remainder of the network.
assortNest[, Degree := 0]
dailySimilNest = vector(mode='list', length= length(colnames(btwIndivNest)[-c(1,2, .N)]))
i = 1
for(day in colnames(btwIndivNest)[-c(1,2, .N)]){
  cat(day)
  dailySimilNest[[i]] = as.matrix(rbind(cbind(data.table(Hen_108 = rep(NA, 34)), 
                                          dcast(btwIndivNest[,.SD, .SDcols = c("V1", "V2", (day))], V1 ~ V2)[,-1]),
                                    data.table(Hen_108 = NA), fill = T))
  rownames(dailySimilNest[[i]]) = colnames(dailySimilNest[[i]])
  simMatrix = dailySimilNest[[i]]
  simMatrix = simMatrix *10
  #maybe set threshold smarter by checking 
  simMatrix[simMatrix < 9] = 0
  g  <- graph.adjacency(simMatrix,mode = "undirected", weighted = T, diag = F)
  #add network node attribute
  V(g)$domIndex <- socialData$Ratio[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
  assortNest[Date == (day), Assort := assortativity(g, V(g)$domIndex)]
  assortNest[Date == (day), Degree := assortativity_degree(g)]
  i = i+1
}

similSimplyNest = lapply(dailySimilNest, function(x){ x[x < 0.9] = 0; x[(x > 0.9)|(x == 0.9)] = 1; x})

similSimplyNest = Reduce('+', similSimplyNest)
similSimplyNest = similSimplyNest/10
g  <- graph.adjacency(similSimplyNest,mode = "undirected", diag = F)
bc <- edge.betweenness.community(g)
#par(mfrow=c(1,2))
plot(as.dendrogram(bc))
#network vertex names
V(g)$name
#inspect network edge and node attributes
edge_attr(g)
vertex_attr(g)
V(g)$domIndex <- socialData$Ratio[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
plot(g, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     vertex.color=c( "pink", "skyblue")[1+(V(g)$domIndex>0.5)] ) 
assortativity(g, V(g)$domIndex)


###### Within individuals ####

withinIndividuals = similarityWithin(trackingData)
# Reshape data from wide to long format
withinIndividualsL = melt(withinIndividuals,  
                           id.vars     = c("Date"),
                           variable.name = "Hen", 
                           value.name = "Similarity")
withinIndividualsL[, HenID := as.numeric(unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen))))]
withinIndividualsL = withinIndividualsL[socialData[,c(1,4,13)], on = "HenID"]
#include Week of age
withinIndividualsL= withinIndividualsL[tableWoA, on = "Date", nomatch = NULL]
withinIndividualsL[, Highlight := "Any"]
withinIndividualsL[HenID == 82 | HenID == 97 | HenID == 33, Highlight := "Dom"]
withinIndividualsL[HenID == 77 | HenID == 5 | HenID == 108, Highlight := "Sub"]


ggplot(data = na.omit(withinIndividualsL)[Highlight != "Any",], aes(x = Date, y = Similarity, colour = Hen))+
  geom_line(aes(group = Hen))
ggplot(data = na.omit(withinIndividualsL), aes(x = Date, y = Similarity, colour = Ratio))+
  geom_line(aes(group = Hen))

#test model
hist(withinIndividualsL$Similarity)
dataWithin = na.omit(withinIndividualsL)
dataWithin[, NonSimilar := 1-Similarity]
dataWithin[, RatioSplit := "Dom"]
dataWithin[Ratio < 0.5, RatioSplit := "Sub"]

#fwrite(dataWithin, "dataWithin.csv", sep = ";")
dataWithin <- fread("dataWithin.csv")
#rm(dataWithin)

model.Within = lmer(Similarity ~ Ratio*WoA + (1|Pen/HenID), data = dataWithin)
#model.Within = lmer(Similarity ~ Ratio + (WoA|HenID), data = dataWithin)
null.Within = lmer(Similarity ~ 1 + (1|Pen/HenID), data = dataWithin)
resid.Within = simulateResiduals(model.Within, 1000)
plot(resid.Within) 
plotResiduals(resid.Within, form = dataWithin$Ratio)
plotResiduals(resid.Within, form = dataWithin$WoA)
summary(model.Within)
parameters(model.Within)
plot(allEffects(model.Within))
r.squaredGLMM(model.Within, null.Within)

betaTableWithin = data.table(Run = "Orig", t(fixef(model.Within)))

set.seed(42)
for (i in 1:1000){
  cat("Run:",i, "\n")
  dataWithin[, Similarity_rand := sample(Similarity), by = Date]
  model.Within.rand = lmer(Similarity_rand ~ Ratio*WoA + (1|Pen/HenID), data = dataWithin)
  entry = data.table(Run = as.character(i), t(fixef(model.Within.rand)))
  betaTableWithin= rbind(betaTableWithin, entry)
}

ggplot(data = betaTableWithin, aes(x= `Ratio:WoA`))+
  geom_histogram()+
  geom_vline(xintercept = betaTableWithin[Run == "Orig", `Ratio:WoA`], colour = "red")

p = (sum(betaTableWithin[Run != "Orig", `Ratio:WoA`] <= betaTableWithin[Run == "Orig", `Ratio:WoA`]) + 1)/dim(betaTableWithin)[1]


dataWithin[, PredictWithin := predict(model.Within)]

ggplot(dataWithin, aes(x = Date, y = Similarity, color = Ratio)) +
  geom_jitter(size=2) + 
  geom_smooth(aes(group = Hen), se = F)+
  #geom_line(data = dataWithin[, mean(PredictWithin), by = .(Date, RatioSplit)], 
  #          aes(x = Date, y =V1, colour = RatioSplit),size=1.5) + 
  theme_classic(base_size = 18)+ 
  ylab("daily within-individual similarity")#+ 

ggplot(dataWithin[WoA %in% quantile(WoA),], aes(x = Ratio, y = Similarity)) +
  geom_point(size=2) + 
  geom_smooth(data = dataWithin[WoA %in% quantile(WoA), mean(PredictWithin), by = .(WoA, Ratio)], 
              aes(y =V1), formula = y~x, method = lm, size=1.5, colour = "red") + 
  theme_classic(base_size = 18)+ 
  facet_grid(~WoA)+
  ylab("daily within-individual similarity")


#plot individual variation
ggplot(data = dataWithin, aes(x = WoA, y = PredictWithin)) + 
  geom_line(aes(group = HenID, color = RatioSplit), size=1)+
  #geom_line(data = dataWithin[Highlight != "Any",], aes(group = HenID, colour = as.factor(HenID)), size=2)+
  labs(y = "daily within-individual similarity")+ 
  theme_classic(base_size = 18)

#plot against actual data of some examples
ggplot(data = dataWithin[Highlight != "Any",], 
       aes(x = WoA, colour = RatioSplit)) + 
  geom_jitter(aes(y = Similarity))+
  geom_line(aes(y = PredictWithin, group = HenID), size=1)+
  labs(y = "daily within-individual similarity",color = "Hen ID")+ 
  theme_classic(base_size = 18)+ 
  guides(color = guide_legend(nrow = 4))
ggplot(data = dataWithin[Highlight != "Any",], 
       aes(x = WoA, colour = Hen)) + 
  geom_jitter(aes(y = Similarity))+
  geom_line(aes(y = PredictWithin, group = HenID), size=1)+
  labs(y = "daily within-individual similarity",color = "Hen ID")+ 
  theme_classic(base_size = 18)+ 
  guides(color = guide_legend(nrow = 4))


#divide the variance explained by animal_id by the total phenotypic variance 
#(animal_id + month:year + year + residual variance)
print(VarCorr(model.Within), comp = c("Variance", "Std.Dev."))
VarCorr(model.Within)$"HenID:Pen"[1] / (VarCorr(model.Within)$"HenID:Pen"[1] + 
                                          VarCorr(model.Within)$"Pen"[1] + 
                                          attr(VarCorr(model.Within), "sc")^2)
#-> within -individual repeatability of 0.3 (within a pen)

#to get uncertainity estimate: simulate model 1000 times
set.seed(1) 
library(arm)
simulated <- sim(model.Within, n.sim = 1000)
posterior_HenID <- apply(simulated@ranef$"HenID:Pen"[ , , 1],1,var) 
posterior_Pen <- apply(simulated@ranef$"Pen"[ , , 1],1,var) 
posterior_residual <- simulated@sigma^2
quantile(posterior_HenID / (posterior_HenID + posterior_Pen + posterior_residual), prob=c(0.025, 0.5, 0.975))

#coefficient of variation for between individual variance
CVi <- sqrt(posterior_HenID) / summary(model.Within)$coefficients[1] 
quantile(CVi,prob=c(0.025, 0.5, 0.975))

#behavioural type
#the behavioral type is the best linear unbiased prediction (BLUP) of the 
#random effect, i.e. the prediction for mean behavioral expression for each individual.
library(merTools)
randomSims <- as.data.table(REsim(model.Within, n.sims = 1000)) 
head(randomSims[randomSims$groupFctr=="HenID:Pen",])
# add the dominance ratio of the individual 
randomSims[, HenID := as.numeric(unlist(strsplit(groupID,':'))[2*(1:dim(randomSims)[1])-1])] 
randomSims <- merge(randomSims[randomSims$groupFctr=="HenID:Pen",], 
                    dataWithin[!duplicated(HenID),c("HenID","Ratio", "Pen", "RatioSplit")])
# add identifier to color individuals uniquely 
randomSims[, ID := ifelse(HenID %in% c(77, 5, 108, 82, 97, 33), as.character(HenID), "Other individuals")]

#add population intercept for easier interpretation on transition scale
randomSims[, meanSimil := mean + fixef(model.Within)["(Intercept)"]]

# order by transitions
randomSims = randomSims[order(meanSimil),]
randomSims[,HenID := factor(HenID, levels = as.character(HenID))]

#plot
ggplot(data = randomSims, aes(x = as.factor(HenID), y = meanSimil))+
  geom_errorbar(aes(ymin = meanSimil-sd, ymax = meanSimil+sd, color = RatioSplit), size = 2)+
  geom_point(aes(fill = ID), shape = 21, size = 5)+
  theme_classic(base_size = 18)+
  scale_fill_manual(values = c("red","blue","orange","yellow", "lightblue","white", "gray"))+
  scale_color_manual(values = c("grey", "black"))

#### Time series plots #########################
#relevant times
times = list(ymd_hms(c("2019-11-10 04:00:00", "2019-11-10 17:30:00")),
             ymd_hms(c("2020-01-10 02:00:00", "2020-01-10 17:30:00")),
             ymd_hms(c("2020-03-13 02:00:00", "2020-03-13 17:30:00")),
             ymd_hms(c("2020-05-08 02:00:00", "2020-05-08 17:30:00")),
             ymd_hms(c("2020-06-26 02:00:00", "2020-06-26 17:30:00")))

hen_list <- vector(mode='list', length=length(hens))
splitHen = splitHenData(trackingData)
for (j in 1:length(hens)) {
  hen_list[[j]] <- vector(mode='list', length=5)
  for (i in 1:5){
    hen_list[[j]][[i]] = extractInterval(splitHen[[j]], times[[i]][1], times[[i]][2])
  }
} 


#make overview plots per hen
for (i in 1:length(hens)){
  plotData = mergeHenData(hen_list[[i]])
  plotData[, Date := as.factor(as_date(Time))]
  plotData[, Time_x := as_hms(Time)]
  plotData[, Zone := factor(Zone, levels= c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
  ggplot(plotData, aes(x = Time_x, y = Zone)) + 
    geom_step(group = 1) + 
    geom_point() + 
    facet_grid(Date~.)+
    labs(x = "time", y = "Zones") +
    ggtitle(hens[i])
  #scale_x_discrete(labels = x_times) + 
  #theme(axis.ticks = element_blank(), axis.text.y = element_text(size = 6))+
  #scale_y_discrete(labels = c("WG", "1", "2", "3", "4"))
  
  ggsave(filename = paste0(hens[i],".png"), plot = last_plot(), 
         width = 18, height = 20, dpi = 300, units = "cm")
}

plotData = trackingData
plotData[, Date := as.factor(Date)]
plotData[, Time_x := as_hms(Time)]
plotData[, Zone := factor(Zone, levels= c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]

# #ideas for plots: for each individual an average step function plot
# #requires: full sequence for every day
# 
# splitHen = splitHenData(trackingData)
# hens = unique(trackingData$HenID)
# average_hen <- vector(mode='list', length=length(hens))
# i = 1
# for (i in 1:length(hens)) {
#   cat("Hen", i)
#   list = fillSeqHenData(splitHen[[i]])
#   #cat("done")
#   average_hen[[i]] = averageDailyZone(list)
# } 
# 
# 
# #make overview plots per hen
# for (i in 1:length(hens)){
#   plotData = average_hen[[i]]
#   plotData[, Zone := factor(Zone, levels= c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
#   ggplot(plotData, aes(x = Time, y = Zone)) + 
#     geom_step(group = 1) + 
#     geom_point() + 
#     labs(x = "time", y = "Zones") +
#     ggtitle(hens[i])+
#     theme_bw(base_size = 18)
#   
#   
#   ggsave(filename = paste0(hens[i],".png"), plot = last_plot(), 
#          width = 40, height = 18, dpi = 300, units = "cm")
# }

### Test asnipe Gaussian mixture model approach ####


#take out Day and Light indicators from data!

library(asnipe)
# Generate GMM data
gmm_data <- gmmevents(time=trackingData[TrueTransition == T, as.numeric(as.ITime(Time))],
                      identity=trackingData[TrueTransition == T, Hen],
                      location=trackingData[TrueTransition == T,paste(Date, Zone, sep = "_")],
                      global_ids=trackingData[, unique(Hen)])

save(gmm_data, file = "gmm_data.Rda")
# Extract output
gbi <- as.data.table(gmm_data$gbi)
events <- as.data.table(gmm_data$metadata)
observations_per_event <- as.data.table(gmm_data$B)

# Split up location and date data
events[, Date := tstrsplit(Location, "_", fixed=TRUE, keep=c(1))]
events[, c("Location1", "Location2") := tstrsplit(Location, "_", fixed=TRUE, keep=c(2,3))]
events[, Location := ifelse(is.na(Location2),Location1,paste0(Location1,"_", Location2))]
events[, c("Location1", "Location2") := NULL]
events[, Date := as.IDate(Date)]

events = tableWoA[events, on = "Date"]

matricesWoA = vector(mode='list', length= length(unique(events$WoA)))

i = 1
for (woa in unique(events$WoA)){
  matricesWoA[[i]] = get_network(gbi[which(events$WoA == woa),], data_format="GBI",
                                 association_index="SRI")
  i = i+1
}


par(mfrow = c(1,2))
g1  <- graph.adjacency(matricesWoA[[1]],mode = "undirected", weighted = TRUE, diag = F)
V(g1)$domIndex <- socialData$Ratio[match(V(g1)$name, paste0("Hen_",as.character(socialData$HenID)))]
V(g1)$Pen <- socialData$Pen[match(V(g1)$name, paste0("Hen_",as.character(socialData$HenID)))]
plot(g1, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
      vertex.color=c( "pink", "skyblue")[1+(V(g1)$domIndex>0.5)] ) 

g2  <- graph.adjacency(matricesWoA[[32]],mode = "undirected", weighted = TRUE, diag = F)
V(g2)$domIndex <- socialData$Ratio[match(V(g2)$name, paste0("Hen_",as.character(socialData$HenID)))]
V(g2)$Pen <- socialData$Pen[match(V(g2)$name, paste0("Hen_",as.character(socialData$HenID)))]
plot(g2, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     vertex.color=c( "pink", "skyblue")[1+(V(g2)$domIndex>0.5)] ) 


meanData = data.table(WoA = unique(events$WoA),
                      Mean = unlist(lapply(matricesWoA, function(x){ mean(x, na.rm = TRUE)})),
                      Median = unlist(lapply(matricesWoA, function(x){ median(x, na.rm = TRUE)})),
                      Max = unlist(lapply(matricesWoA, function(x){ max(x, na.rm = TRUE)})),
                      SD = unlist(lapply(matricesWoA, function(x){ sd(x, na.rm = TRUE)})),
                      Q1 =unlist(lapply(matricesWoA, function(x){ quantile(x, 0.25, na.rm = TRUE)})),
                      Q3 = unlist(lapply(matricesWoA, function(x){ quantile(x, 0.75, na.rm = TRUE)})))

ggplot(meanData, aes(x = WoA))+
  geom_pointrange(aes(y = Median, ymin = Q1, ymax = Q3))+
  geom_point(aes(y = Mean), color = "red")+
  #geom_point(aes(y = Max))+
  theme_classic(base_size = 18)


dataBetween2 = na.omit(rbindlist(lapply(matricesWoA, function(x){as.data.table(as.table(x))})))
colnames(dataBetween2) = c("Hen1", "Hen2", "Similarity")
nPairs = length(unique(paste(dataBetween2$Hen1, dataBetween2$Hen2)))
dataBetween2[, WoA := rep(unique(events$WoA), each = nPairs)]
dataBetween2 = dataBetween2[Hen1 != Hen2,]

dataBetween2[, Hen1_Pen := socialData$Pen[match(Hen1, socialData$Hen)]]
dataBetween2[, Hen2_Pen := socialData$Pen[match(Hen2, socialData$Hen)]]
dataBetween2[, Hen1_Ratio := socialData$Ratio[match(Hen1, socialData$Hen)]]
dataBetween2[, Hen2_Ratio := socialData$Ratio[match(Hen2, socialData$Hen)]]
dataBetween2[, PenBool := Hen1_Pen == Hen2_Pen]
dataBetween2[, DiffRatio := abs(Hen1_Ratio - Hen2_Ratio)]


#build Model
hist(dataBetween2$Similarity)

model.Between = lmer(Similarity ~ DiffRatio*PenBool+ DiffRatio*WoA + (1|Hen1) + (1|Hen2), data = dataBetween2)
null.Between = lmer(Similarity ~ 1+ (1|Hen1) + (1|Hen2), data = dataBetween2)
anova(model.Between, null.Between)
resid.Between = simulateResiduals(model.Between, 1000)
plot(resid.Between)
plotResiduals(resid.Between, form = dataBetween2$DiffRatio)
plotResiduals(resid.Between, form = dataBetween2$PenBool)
plotResiduals(resid.Between, form = dataBetween2$WoA)
summary(model.Between)
parameters(model.Between)
plot(allEffects(model.Between))
#variance explained by fixed factors and entire model
r.squaredGLMM(model.Between, null.Between)

dataBetween2[, PredictSimilarity := predict(model.Between)]

ggplot(dataBetween2[WoA %in% round(quantile(WoA)),], aes(x = DiffRatio, y = Similarity)) +
  geom_point(size=2) + 
  geom_smooth(data = dataBetween2[WoA %in% round(quantile(WoA)), mean(PredictSimilarity), by = .(WoA, DiffRatio)], 
              aes(y =V1), formula = y~x, size=1.5, colour = "red") + 
  theme_classic(base_size = 18)+ 
  facet_grid(~WoA)+
  ylab("daily between-individual similarity")

plotData = as.data.table(emmeans(model.Between, ~ pairwise ~ WoA*DiffRatio, 
                                 at =  list(DiffRatio = round(quantile(dataBetween2$DiffRatio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)

ggplot()+ 
  geom_jitter(data = dataBetween2, aes(x = WoA, y = Similarity), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = DiffRatio,  colour = as.factor(DiffRatio)), size = 1)+
  #geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

