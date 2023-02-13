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

#for repeatable results
set.seed(42)

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

##### Parameters ########

# calculate daily parameters per bird 

###### vertical travel distance ############### 
# number of vertically crossed zones during light hours, divided by the seconds of the animals spent inside
#create vector where wintergarden doesn't exist (replaced by litter)
trackingData[, distZone := Zone]
trackingData[Zone == "Wintergarten", distZone := "Litter"]
#add zone vector shifted by one
trackingData[, nextZone := c(distZone[-1], NA), by = HenID]
#calculate distance travelled by using function defineDistance
trackingData[, distVertical := apply(X = cbind(distZone, nextZone), MARGIN = 1, FUN= defineDistance), by = HenID]
#daylight selection happens when data is joined in section all parameters


###### Nestbox zone ########
#Nestbox entries per bird
#TODO: exclude hens which sleep in zone
#sift out only those in the morning, relevant for egg laying not resting (until 9)
#extract if hen was in nest zone on day or not
dailyNest = trackingData[Light == T & hour(Time) < 9, .(NestZone = ifelse(any(Zone == "Ramp_Nestbox"), 1, 0)), by = .(HenID, Date)][order(Date, HenID)]
#extract how long each hen was in the box
dailyNest[NestZone == 1, DurationNest := trackingData[Light == T & hour(Time) < 9 & Zone == "Ramp_Nestbox", sum(Duration), by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[ , NestSleep := trackingData[LightIndic == T & hour(Time) < 9, Zone == "Ramp_Nestbox", by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[NestZone == 0 | NestSleep == TRUE, DurationNest := 0]
#extract when median duration in the nest is reached
dailyNest[, MedDurNest := ifelse(NestZone == 1,  round(DurationNest/2),NA) ]
# median time point for half duration in nestbox
dailyNest[NestZone == 1, EntryNest := trackingData[Light == T & hour(Time) < 9 & Zone == "Ramp_Nestbox", Time[1], by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[NestZone == 1 & DurationNest > 0, MedTimeNest := as.ITime(EntryNest + MedDurNest)]

#switches in and out of the nestbox zone
# helper = trackingData[Light == T & hour(Time) < 9 & Zone == "Ramp_Nestbox" & LightIndic != 1, .(SwitchesNest = .N),by = .(HenID, Date)]
# dailyNest = helper[dailyNest, on = c("HenID", "Date")]
# helper = dailyNest[, .(HenID, Date, MedTimeNest)][trackingData, on = c("HenID", "Date")][Time < MedTimeNest & Time > (MedTimeNest -(3600)),]
# helper[, nextZone := c(distZone[-1], NA), by = HenID]
# #calculate distance travelled by using function defineDistance
# helper[, distVertical := apply(X = cbind(distZone, nextZone), MARGIN = 1, FUN= defineDistance), by = HenID]
# 
# dailyNest = helper[,.(preNestDist = sum(distVertical)),by = .(HenID, Date)][dailyNest, on = c("HenID", "Date")]


###### sleeping spot ############

#durations per zone per bird per day during dark hours
dailySleep = trackingData[Light == F, .(DurationNight = sum(Duration)), by = .(HenID, Zone, NightCycle)][order(NightCycle, HenID)]
dailySleep = dailySleep[, .SD[which.max(DurationNight)], by = .(HenID, NightCycle)]
#extract maximum Zone for each bird per day
dailySleep[, onTop := ifelse(Zone == "Tier_4", 1, 0)]
colnames(dailySleep)[3] = "ZoneSleep"


# ###### wintergarden use #####
# #extract if hen goes out on day or not
# dailyGarten = trackingData[, .(Out = ifelse(any(Zone == "Wintergarten"), 1, 0)), by = .(HenID, Date)][order(Date, HenID)]
# #extract how long each hen went out per day
# dailyGarten[Out == 1, DurationGarten := trackingData[Zone == "Wintergarten", .(sum(Duration)), by = .(HenID, Date)][order(Date, HenID), V1]]
# dailyGarten[Out == 0, DurationGarten := 0]
# # latency to go out
# #careful: on vaccination days garten opened later! take out vacc days
# vacc =  c(ymd("2019-11-08"), ymd("2019-12-24"), ymd("2020-01-14"),
#           ymd("2020-02-18"), ymd("2020-03-03"), ymd("2020-04-14"), 
#           ymd("2020-06-09"), ymd("2020-06-23"))
# 
# dailyGarten[Out == 1, EntryGarten := trackingData[Zone == "Wintergarten", Time[1], by = .(HenID, Date)][order(Date,HenID), V1]]
# dailyGarten[Out == 1, LatencyGarten := EntryGarten - ymd_hms(paste(Date, "10:00:00"))]
# dailyGarten[Date %in% vacc, LatencyGarten := NA]


###### Feeder reactivity #######
#Feeder runs: (ab 22.11.: 2:00), 4:00, 6:00, 8:00, 10:00, 13:00, 15:00, 16:15

dailyFeed = infeedZone(trackingData)
colnames(dailyFeed)[2] = "DurationFeed2"
colnames(dailyFeed)[3] = "DurationFeed4"
dailyFeed[, HenID := as.numeric(unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen))))]

# dailyFeed = feedReactivity(trackingData)
# 
# dailyFeedL = melt(dailyFeed,  
#                   id.vars     = c("Date"),
#                   variable.name = "Hen", 
#                   value.name = "FeedReact")
# dailyFeedL[, HenID := as.numeric(unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen))))]
# dailyFeedL = dailyFeedL[order(Date, HenID),]
# dailyFeedL[!(paste(Date, HenID) %in% trackingData[, unique(paste(Date, HenID))]), FeedReact := NA]
# dailyFeedL[!is.na(FeedReact), FeedZoneDur := trackingData[(Zone == "Tier_4" |Zone == "Tier_2") & Light == T, sum(Duration), by = .(HenID, Date)][order(Date, HenID), V1]]
# dailyFeedL = trackingData[(Zone == "Tier_4" ) & Light == T, .(FeedZone4 = sum(Duration)), by = .(HenID, Date)][dailyFeedL, on = c("HenID", "Date")]
# dailyFeedL = trackingData[(Zone == "Tier_2" ) & Light == T, .(FeedZone2 = sum(Duration)), by = .(HenID, Date)][dailyFeedL, on = c("HenID", "Date")]

###### All parameters #######
#create data.table containing all daily measures
#vertical travel distance
varOfInterest = trackingData[Light == T, .(vertTravelDist = sum(distVertical)), by = .(HenID, Date, Pen, WoA)]

#add duration for each day per zone
varOfInterest = varOfInterest[dcast(trackingData[Light == T, sum(Duration), by = .(HenID, Date, Zone)], formula = HenID + Date ~ Zone, value.var = "V1", fill = 0), on = c("HenID", "Date")]
varOfInterest[, TotalDur := rowSums(cbind(Wintergarten, Tier_2, Tier_4, Litter, Ramp_Nestbox))]
#extract max zone per day
dailyMaxZone = trackingData[Light == T, sum(Duration), by = .(HenID, Date, Zone)][, .(MaxZone = Zone[which.max(V1)]), by= .(HenID, Date)]
varOfInterest = varOfInterest[dailyMaxZone, on = c("HenID", "Date")] 

#careful with join direction not to loose the hens, for hens who don't go out for example
# -> full outer join
# sleep: max zone, duration in zone, is zone top tier?
varOfInterest = varOfInterest[dailySleep, on = c(HenID = "HenID", Date = "NightCycle")]
#delete half nights
varOfInterest = varOfInterest[!is.na(Pen),]
# wintergarden: did hen go out? duration in garten, latency to enter 
#varOfInterest = dailyGarten[varOfInterest, on = c(HenID = "HenID", Date = "Date")]
# nestbox: did hen go in nest?, total duration in nest, time point of median duration in nest 
varOfInterest = dailyNest[varOfInterest, on = c(HenID = "HenID", Date = "Date")]
# feeder zones: time spent in feeder zone during run versus outside run
varOfInterest = dailyFeed[varOfInterest, on = c(HenID = "HenID", Date = "Date")]

#add social information (dominance index)
varOfInterest = socialData[, .(HenID,Ratio,Comb)][varOfInterest, on = "HenID"] 

#varOfInterest = henDataLong[, .(HenID,Severity,WoA)][varOfInterest, on = c("HenID", "WoA")] 


##### Parameter inspection plots ####

#Highlight three most extreme
varOfInterest[, Highlight := "Any"]
varOfInterest[HenID == 82 | HenID == 97 | HenID == 33, Highlight := "Dom"]
varOfInterest[HenID == 77 | HenID == 5 | HenID == 108, Highlight := "Sub"]
#make factor out of HenID & Pen with levels sorted in ascending Ratio order
varOfInterest[,HenID := factor(HenID, levels = socialData$HenID)]
varOfInterest[,Pen := factor(Pen)]
varOfInterest[,ZoneSleep := factor(ZoneSleep, levels = c("Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
varOfInterest[,MaxZone := factor(MaxZone, levels = c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
#splitting dominance index by 0.5
varOfInterest[, RatioSplit := "Dom"]
varOfInterest[Ratio < 0.5, RatioSplit := "Sub"]


###### durations in general ####

#heatmap of durations by ratio (hens sorted)
ggplot(trackingData[Light == T, .(Duration = sum(Duration)), by = .(HenID, Zone)], aes(x = Zone, y = factor(HenID, levels = socialData$HenID))) +
  geom_tile(aes(fill = Duration), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")

#heatmap of durations by ratio (hens sorted) by WoA
plot = trackingData[, .(nDays = length(unique(Date))), by = .(HenID, WoA)][
  trackingData[Light == T, .(Duration = sum(Duration)), by = .(HenID, Zone, WoA)][WoA %in% quantile(WoA),], 
  on = .(HenID, WoA)]
plot[, dailyDuration := Duration/nDays]
ggplot(plot, aes(x = Zone, y = factor(HenID, levels = socialData$HenID))) +
  geom_tile(aes(fill = dailyDuration), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")+
  facet_grid(.~WoA)

#most common zone daily
ggplot(varOfInterest, aes(x = Date, y = factor(MaxZone, levels = c("Wintergarten","Litter", "Tier_2", "Ramp_Nestbox", "Tier_4")), colour = Ratio))+
  geom_jitter(width = 0.1)+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

ggplot(varOfInterest, aes(x = Date, y = factor(MaxZone, levels = c("Wintergarten","Litter", "Tier_2", "Ramp_Nestbox", "Tier_4")), colour = Highlight))+
  geom_jitter()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))

#Litter by Ratio
ggplot(varOfInterest, aes(x = Date, y = Litter, colour = Ratio))+
  geom_jitter(width = 0.1)+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

###### travel distance #### 
hist(varOfInterest$vertTravelDist)
#progression of travel distance over days
ggplot(varOfInterest, aes(x = Date, y = vertTravelDist, colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))
ggplot(varOfInterest, aes(x = Date, y = vertTravelDist, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

#compared to simple transition counts
ggplot(socialData[, .(HenID,Ratio,Comb)][trackingData[, .N, by = .(HenID,Date)], on = "HenID"], aes(x = Date, y = N, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

###### sleeping spot ####

hist(varOfInterest$onTop)
hist(varOfInterest[, sum(onTop), by = HenID][,V1])
#heatmap by comb
ggplot(varOfInterest[, .(Duration = sum(DurationNight)), by = .(Comb, ZoneSleep)], aes(x = ZoneSleep, y = as.factor(Comb))) +
  geom_tile(aes(fill = Duration), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")
#heatmap by ratio
ggplot(varOfInterest[, .(Duration = sum(DurationNight)), by = .(HenID, ZoneSleep)], aes(x = ZoneSleep, y = HenID)) +
  geom_tile(aes(fill = Duration), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")
#progression of most common sleeping spot on Top or not
ggplot(varOfInterest, aes(x = Date, y = onTop, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),method = "glm", 
              method.args = list(family = "binomial"),se = F)+
  scale_color_gradient(low = "blue", high = "gold")
# time spent on top during day
ggplot(varOfInterest, aes(x = Date, y = Tier_4, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")
#time spent on top and KBF
ggplot(varOfInterest, aes(x = Date, y = Tier_4, colour = Severity))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")
ggplot(varOfInterest[!is.na(Severity),], aes(x = Severity, y = Tier_4, colour = Ratio))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_grid(~WoA)
  scale_color_gradient(low = "blue", high = "gold")
#relation between time spent on top during day and sleeping spot
ggplot(varOfInterest, aes(x = Tier_4, y = onTop, colour = Ratio))+
  geom_point()+
  geom_smooth(color = "red", method = "glm", 
              method.args = list(family = "binomial"))+
  scale_color_gradient(low = "blue", high = "gold")


###### Wintergarden ####

#Went out or not?
hist(varOfInterest$Out)
hist(varOfInterest[, sum(Out), by = HenID][,V1])
ggplot(varOfInterest, aes(x = WoA, y = Out, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),method = "glm", 
              method.args = list(family = "binomial"),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

#latency to enter
hist(as.numeric(varOfInterest$LatencyGarten))
#progression of latency to enter the garden over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(LatencyGarten), colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

#duration in the wintergarten
hist(as.numeric(varOfInterest$DurationGarten))
#progression of latency to enter the garden over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(DurationGarten), colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

###### Nestbox zone ####

#Went in nest or not?
hist(varOfInterest$NestZone)
hist(varOfInterest[, sum(NestZone), by = HenID][,V1])
#heatmap by ratio
ggplot(varOfInterest[, .(Sum = .N), by = .(HenID, NestZone)], aes(x = NestZone, y = HenID)) +
  geom_tile(aes(fill = Sum), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")
ggplot(varOfInterest, aes(x = Date, y = NestZone, colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = HenID),method = "glm", 
              method.args = list(family = "binomial"),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))


hist(as.numeric(varOfInterest$MedTimeNestPure))
#progression of time to enter the nestbox over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(MedTimeNestPure), colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")


#duration in the nestbox zone
hist(as.numeric(varOfInterest$DurationNest))
#progression of duration in nestbox over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(DurationNest), colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

#switches into nestbox zone
hist(varOfInterest$preNestDist)
hist(varOfInterest$SwitchesNest)
#which variable is better fitting?
ggplot(varOfInterest, aes( x= preNestDist, y = SwitchesNest))+
  geom_point()+
  geom_smooth()
cor.test(varOfInterest$preNestDist, varOfInterest$SwitchesNest)

#progression of duration in nestbox over days
ggplot(varOfInterest, aes(x = Date, y = SwitchesNest, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

ggplot(varOfInterest, aes(x = Date, y = preNestDist, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")




###### Feed reactivity ####

hist(varOfInterest$FeedReact)
ggplot(varOfInterest[WoA %in% quantile(WoA),], aes(x = FeedReact, fill = Highlight))+
  geom_density(alpha = 0.5, position = "identity")+
  facet_grid(~ WoA)
ggplot(varOfInterest[WoA %in% quantile(WoA),], aes(x = FeedReact, fill = RatioSplit))+
  geom_density(alpha = 0.5, position = "identity")+
  facet_grid(~ WoA)  
#feed reactivity over days
ggplot(varOfInterest, aes(x = Date, y = FeedReact, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_color_gradient(low = "blue", high = "gold")
#total time in feeding zone by ratio
ggplot(varOfInterest, aes(x = Date, y = FeedZoneDur, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)
ggplot(varOfInterest, aes(x = Date, y = FeedZone4, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)
ggplot(varOfInterest, aes(x = Date, y = FeedZone2, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)
#relationship total time in zone and reactivity
ggplot(varOfInterest, aes(x = FeedZoneDur, y = FeedReact, colour = Ratio))+
  geom_point()+
  geom_smooth(colour = "red",se = T)
ggplot(varOfInterest, aes(x = FeedZone4, y = FeedReact, colour = Ratio))+
  geom_point()+
  geom_smooth(colour = "red",se = T)+
scale_color_gradient(low = "blue", high = "gold")
ggplot(varOfInterest, aes(x = FeedZone2, y = FeedReact, colour = Ratio))+
  geom_point()+
  geom_smooth(colour = "red",se = T)+
scale_color_gradient(low = "blue", high = "gold")
ggplot(varOfInterest, aes(x = FeedZone2, y = FeedZone4, colour = Ratio))+
  geom_point()+
  geom_smooth(colour = "red",se = T)+
  scale_color_gradient(low = "blue", high = "gold")


#TODO: take out outliers
cor.test(varOfInterest$FeedZoneDur, varOfInterest$FeedReact)

##### Test models ######


###### 1 daily duration in zone 4 ####

hist(varOfInterest$Tier_4)

model.Duration = glmer(Tier_4 ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest, family = poisson)
summary(model.Duration) #pen creates singularity issue
model.Duration = glmer(Tier_4 ~ Ratio*WoA + (1|HenID), data = varOfInterest, family = poisson)
resid.Duration = simulateResiduals(model.Duration, 1000)
plot(resid.Duration) #zero inflation issue most likely
plotResiduals(resid.Duration, form = varOfInterest$Ratio)
plotResiduals(resid.Duration, form = varOfInterest$WoA)
testZeroInflation(model.Duration)
model.Duration = glmmTMB(Tier_4 ~ Ratio*WoA + (1|HenID), data = varOfInterest, family = poisson)
test.Duration = glmmTMB(Tier_4 ~ Ratio*WoA + (1|HenID), data = varOfInterest, family = poisson, 
                         zi = ~1)
AIC(model.Duration, test.Duration) # including general inflation term makes fit better
resid.Duration = simulateResiduals(test.Duration, 1000)
plot(resid.Duration) #deviation okay
plotResiduals(resid.Duration, form = varOfInterest$Ratio)
plotResiduals(resid.Duration, form = varOfInterest$WoA)

null.Duration = glmmTMB(Tier_4 ~ 1 + (1|HenID), data = varOfInterest, family = poisson, 
                        zi = ~1)
AIC(test.Duration, null.Duration) #much better than null model
plot(allEffects(test.Duration)) 
summary(test.Duration)
tidy(test.Duration, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
r.squaredGLMM(test.Duration, null.Duration)

#compare against binomial and non-zero models if estimate direction is the same 
varOfInterest[, on4 := ifelse(Tier_4 >0, 1, 0)]
split1.Duration = glmer(on4 ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest, family = binomial)
resid.Duration = simulateResiduals(split1.Duration, 1000)
plot(resid.Duration) #deviation okay
plotResiduals(resid.Duration, form = varOfInterest$Ratio)
plotResiduals(resid.Duration, form = varOfInterest$WoA)
parameters(split1.Duration, exp = TRUE) # the same

split2.Duration = glmer(Tier_4 ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest[Tier_4 >0,], family = poisson)
resid.Duration = simulateResiduals(split2.Duration, 1000)
plot(resid.Duration) #deviation okay
plotResiduals(resid.Duration, form = varOfInterest[Tier_4 >0, Ratio])
plotResiduals(resid.Duration, form = varOfInterest[Tier_4 >0, WoA])
parameters(split2.Duration, exp = TRUE) # the same

plotData = as.data.table(emmeans(test.Duration, ~ pairwise ~ WoA*Ratio, type = "response",
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))$emmeans)


varOfInterest[, PredictDuration4 := predict(test.Duration)]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = Tier_4), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = rate, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = lower.CL, ymax = upper.CL, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)


###### 2 Vertical distance ####

hist(varOfInterest$vertTravelDist)
#TODO: interaction poly?
model.Travel = lmer(vertTravelDist ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest)
resid.Travel = simulateResiduals(model.Travel, 1000)
plot(resid.Travel) #deviation okay -> potential polynomial effect
plotResiduals(resid.Travel, form = varOfInterest$Ratio)
plotResiduals(resid.Travel, form = varOfInterest$WoA)#-> polynomial effect?
test.Travel = lmer(vertTravelDist ~ Ratio:WoA + Ratio + poly(WoA,2) + (1|Pen/HenID), data = varOfInterest)
resid.Travel = simulateResiduals(test.Travel, 1000)
plot(resid.Travel) #looks better
plotResiduals(resid.Travel, form = varOfInterest$Ratio)
plotResiduals(resid.Travel, form = varOfInterest$WoA)#better
AIC(model.Travel, test.Travel)# poly is better

null.Travel = lmer(vertTravelDist ~ 1 + (1|Pen/HenID), data = varOfInterest)
AIC(test.Travel, null.Travel)#test is better
summary(test.Travel)
parameters(test.Travel)
plot(allEffects(test.Travel))
r.squaredGLMM(test.Travel, null.Travel)

plotData = as.data.table(emmeans(test.Travel, ~ pairwise ~ poly(WoA,2)+WoA:Ratio+ Ratio, 
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))$emmeans)


varOfInterest[, PredictTravel := predict(test.Travel)]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = vertTravelDist), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

ggplot(varOfInterest, aes(x = Date, y = vertTravelDist, color = RatioSplit)) +
  geom_jitter(size=2) + 
  geom_line(data = varOfInterest[, mean(PredictTravel), by = .(Date, RatioSplit)], 
                              aes(x = Date, y =V1, colour = RatioSplit),size=1.5) + 
  theme_classic(base_size = 18)+ 
  ylab("daily vertical travel distance")#+ 

#plot individual variation
ggplot(data = varOfInterest, aes(x = WoA, y = PredictTravel, colour = as.factor(HenID))) + 
  geom_line(aes(group = HenID), size=1)+
  labs(y = "Predicted vertical travel distance",color = "Hen ID")+ 
  theme_classic(base_size = 18)+ theme(legend.position="none")+ 
  guides(color = guide_legend(nrow = 4))

#plot against actual data of some examples
ggplot(data = varOfInterest[Highlight != "Any",], 
       aes(x = WoA, colour = as.factor(HenID))) + 
  geom_jitter(aes(y = vertTravelDist))+
  geom_line(aes(y = PredictTravel, group = HenID), size=1)+
  labs(y = "Predicted vertical travel distance",color = "Hen ID")+ 
  theme_classic(base_size = 18)+ 
  guides(color = guide_legend(nrow = 4))



###### 3 Sleeping spot #####

hist(varOfInterest$onTop)
hist(varOfInterest[, sum(onTop), by = HenID][,V1])
model.Sleep = glmer(onTop ~ Ratio*WoA+ (1|Pen/HenID), data = varOfInterest, family = binomial)
#singularity due to Pen
model.Sleep = glmer(onTop ~ Ratio*WoA + (1|HenID),  data = varOfInterest, family = binomial)
null.Sleep = glmer(onTop ~ 1 + (1|HenID),  data = varOfInterest, family = binomial)
resid.Sleep = simulateResiduals(model.Sleep, 1000)
plot(resid.Sleep)#looks good
plotResiduals(resid.Sleep, form = varOfInterest$Ratio[!is.na(varOfInterest$onTop)])
plotResiduals(resid.Sleep, form = varOfInterest$DateID[!is.na(varOfInterest$onTop)])
AIC(model.Sleep, null.Sleep) #better than null model
summary(model.Sleep)
parameters(model.Sleep, exp = TRUE) #TODO: why increased probability with age and WOA???
plot(allEffects(model.Sleep))
r.squaredGLMM(model.Sleep, null.Sleep)

plotData = as.data.table(emmeans(model.Sleep, ~ pairwise ~ WoA*Ratio, 
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)


varOfInterest[, PredictSleep:= predict(model.Sleep)]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = onTop), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = prob, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)



###### 4 Nestbox Time #####

#only use data from when the lights went on at 02:00 otherwise too affected -> 12 days gone
varOfInterest[, fullCycle := !(day(Date) <21 & month(Date) == 11)]
nestData = varOfInterest[fullCycle == TRUE,]

hist(as.numeric(nestData[,MedTimeNest]))
model.Nest = lmer(as.numeric(MedTimeNest) ~ Ratio*WoA+ (1|Pen/HenID), data = nestData)
resid.Nest = simulateResiduals(model.Nest, 1000)
plot(resid.Nest) #polynomial necessary
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), Ratio])
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), WoA]) #poly

test.Nest = lmer(as.numeric(MedTimeNest) ~ Ratio+ Ratio:WoA +poly(WoA,2)+ (1|Pen/HenID), data = nestData)
null.Nest = lmer(as.numeric(MedTimeNest) ~ 1+ (1|Pen/HenID), data = nestData)
#model without full interaction worse
AIC(model.Nest, test.Nest) #poly better
AIC(test.Nest, null.Nest) #poly better

summary(test.Nest)
parameters(test.Nest)
plot(allEffects(test.Nest))
r.squaredGLMM(test.Nest, null.Nest)

plotData = as.data.table(emmeans(test.Nest, ~ pairwise ~ Ratio+ Ratio:WoA +poly(WoA,2), 
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)
plotData[, Time := as.ITime(response)]
plotData[, Plus := as.ITime(asymp.UCL)]
plotData[, Minus := as.ITime(asymp.LCL)]

nestData[, PredictNest:= predict(test.Nest)]

ggplot()+ 
  geom_jitter(data = nestData, aes(x = WoA, y = MedTimeNest), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = Time, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = Minus, ymax = Plus, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)


###### correlations between parameters ####

cor.test(varOfInterest$Tier_4, varOfInterest$vertTravelDist) # negative correlation -0.4
cor.test(varOfInterest$Tier_4, as.numeric(varOfInterest$MedTimeNest)) #negative correlation -0.11
cor.test(varOfInterest$Tier_4,varOfInterest$onTop)#positive correlation 0.37
cor.test(varOfInterest$vertTravelDist, as.numeric(varOfInterest$MedTimeNest))# tiny negative correlation -0.05
cor.test(varOfInterest$vertTravelDist, varOfInterest$onTop)# negative correlation -0.12
cor.test(varOfInterest$onTop, as.numeric(varOfInterest$MedTimeNest)) #no correlation -0.02


# ###### Wintergarten ####
# 
# #start model
# hist(varOfInterest$Out)
# model.Garten = glmer(Out ~ Ratio*WoA + (1|Pen/HenID), family = binomial, data = varOfInterest)
# resid.Garten = simulateResiduals(model.Garten, 1000)
# plot(resid.Garten)#good
# plotResiduals(resid.Garten, form = varOfInterest$Ratio)
# plotResiduals(resid.Garten, form = varOfInterest$WoA)
# summary(model.Garten)
# plot(allEffects(model.Garten))
# parameters(model.Garten, exp = T)
# plotData = as.data.table(emmeans(model.Garten, ~ pairwise ~ WoA*Ratio, 
#                                  at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
#                                   WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)
# 
# 
# varOfInterest[, PredictGarten:= predict(model.Garten)]
# 
# ggplot()+ 
#   geom_jitter(data = varOfInterest, aes(x = WoA, y = Out), height = 0.02, size = 1, alpha = 0.1)+
#   geom_line(data = plotData, aes(x = WoA, y = prob, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
#   geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio), alpha = 0.1)+
#   #facet_grid(.~Pen)+
#   theme_classic(base_size = 18)
# 
# 
# #model of duration very difficult to fit
# #model of latency equally bad
# # hist(as.numeric(varOfInterest$LatencyGarten))
# # model.Garten = glmer.nb(as.numeric(LatencyGarten) ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest)
# # model.Garten = glmer(as.numeric(LatencyGarten) ~ Ratio*WoA + (1|HenID), family = poisson, data = varOfInterest)
# # resid.Garten = simulateResiduals(model.Garten, 1000)
# # plot(resid.Garten) #okayish
# # plotResiduals(resid.Garten, form = varOfInterest$Ratio)
# # plotResiduals(resid.Garten, form = varOfInterest$WoA)
# # summary(model.Garten)
# # plot(allEffects(model.Garten))
# 
# 
# varOfInterest[, PredictGarten := predict(model.Garten, type="response")]
# 
# ggplot(varOfInterest, aes(x = WoA, y = Out, color = RatioSplit)) +
#   #geom_jitter(size=2) +
#   geom_line(data = varOfInterest[, mean(PredictGarten), by = .(WoA, RatioSplit)],
#             aes(x = WoA, y =V1, colour = RatioSplit),size=1.5) +
#   theme_classic(base_size = 18)+
#   ylab("Out in wintergarten")
# 
# #plot individual variation
# ggplot(data = varOfInterest, aes(x = WoA, y = PredictGarten, colour = as.factor(HenID))) +
#   geom_line(aes(group = HenID), size=1)+
#   labs(y = "Predicted probability to go in wintergarten",color = "Hen ID")+
#   theme_classic(base_size = 18)+ theme(legend.position="none")+
#   guides(color = guide_legend(nrow = 4))
# 
# #plot against actual data of some examples
# ggplot(data = varOfInterest[Highlight != "Any",],
#        aes(x = WoA, colour = as.factor(HenID))) +
#   geom_jitter(aes(y = Out))+
#   geom_line(aes(y = PredictGarten, group = HenID), size=1)+
#   labs(y = "Predicted vertical travel distance",color = "Hen ID")+
#   theme_classic(base_size = 18)+
#   guides(color = guide_legend(nrow = 4))


###### Feed reactivity #####

hist(as.numeric(varOfInterest$FeedReact))
model.Feed = lmer(FeedReact ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest)
model.Feed2 = lmer(FeedReact ~ Ratio*WoA + Ratio*MaxZone + MaxZone*WoA+ (1|Pen/HenID), data = varOfInterest)
model.Feed3 = lmer(FeedReact ~ Ratio*WoA + MaxZone + (1|Pen/HenID), data = varOfInterest)
anova(model.Feed, model.Feed2)
anova(model.Feed2, model.Feed3)
#TODO: Model with maxzone better but informative?? 
resid.Feed = simulateResiduals(model.Feed, 1000)
plot(resid.Feed) 
plotResiduals(resid.Feed, form = varOfInterest$Ratio[!is.na(varOfInterest$FeedReact)])
plotResiduals(resid.Feed, form = varOfInterest$WoA[!is.na(varOfInterest$FeedReact)])
plotResiduals(resid.Feed, form = varOfInterest$MaxZone[!is.na(varOfInterest$FeedReact)])
summary(model.Feed)
plot(allEffects(model.Feed))
parameters(model.Feed)


hist(varOfInterest$FeedZoneDur)
model.Feed = lmer(FeedZoneDur ~ Ratio*WoA+ (1|Pen/HenID), data = varOfInterest)
model.Feed = lmer(FeedZoneDur ~ Ratio*WoA+ (1|HenID), data = varOfInterest)
resid.Feed = simulateResiduals(model.Feed, 1000)
plot(resid.Feed) 
plotResiduals(resid.Feed, form = varOfInterest$Ratio[!is.na(varOfInterest$FeedReact)])
plotResiduals(resid.Feed, form = varOfInterest$WoA[!is.na(varOfInterest$FeedReact)])
summary(model.Feed)
plot(allEffects(model.Feed))


plotData = as.data.table(emmeans(model.Feed, ~ pairwise ~ WoA*Ratio, 
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)


varOfInterest[, PredictFeedReact:= predict(model.Feed)]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = FeedReact), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  #geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)




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