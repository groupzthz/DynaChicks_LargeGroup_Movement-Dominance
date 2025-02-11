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
library(Hmisc) #correlation matrix
library(corrplot)#correlation plot
library(RColorBrewer)
library(ggpubr)#plot combin

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

combDataControls = fread("ControlsComb.csv")
combDataControls[,ID := paste0(Pen,Backpack)]
combDataControls = combDataControls[,.(Comb = mean(Area)), by = .(Pen, Backpack,ID)]

#prepare WoA table 
tableWoA = data.table(Date = seq(ymd("2019-10-11"), ymd("2020-07-02"), by = "day"),
                      WoA = rep(18:55, each = 7))
#load health data
healthData = fread('HA_all.csv')
healthData[,ID := paste0(pen,backpack)]
healthData[, feathers := neck + wings +tail + cloaca+ breast]
healthData[, footproblems:= r_podo + r_bumble +r_injure + l_podo + l_bumble + l_injure]
healthData[, date := as_date(date, format = "%d.%m.%Y")]
healthData = healthData[tableWoA, on = c(date = "Date"), nomatch = NULL]

healthDataControls = healthData[grepl('[0-9]', backpack),]
healthDataControls= healthDataControls[,.(ID, pen, date, weight, feathers, wounds, comb, bare_head, footproblems)]
healthDataControls[, Pen := pen]


#healthData = healthData[(date == "28.10.2019"| date == "09.12.2019" | date == "29.06.2020" | date == "30.06.2020")
#                        & ID %in% unique(observData$ID),]
healthData = healthData[ID %in% unique(observData$ID),]
healthDataWide = dcast(healthData, formula = ID ~ WoA, value.var = list("weight", "feathers", "wounds","comb", "bare_head", "footproblems"))



KBF = fread("KBF_scores.csv")
KBF[, Date := as_date(Date, format = "%d.%m.%Y")]
KBF = KBF[tableWoA, on = "Date", nomatch = NULL]
#KBF = KBF[(Date == "28.10.2019" |Date == "09.12.2019" | Date == "29.06.2020" | Date == "30.06.2020")
#                        & HenID %in% unique(observData$ID),]

KBFControls = KBF[!is.na(no),]
KBFControls[, ID := paste0(HenID,no)]
healthDataControls = KBFControls[healthDataControls, on = "ID"]


KBF = KBF[HenID %in% unique(observData$ID),]
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
combsizeFig = ggplot(data = henData, aes(y = Ratio, x = Comb, ))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  geom_line(aes(y = PredictRatio, colour = as.factor(Pen)), size = 1.5)+
  #facet_grid(.~Pen)+
  ylab("Aggression value")+
  xlab("Comb size (cm²)")+
  scale_colour_brewer(palette = "Paired", "Pen")+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

ggsave("CombSize.tiff", combsizeFig, "tiff", width = 18, height= 14, units = "cm", dpi = 300)



###### Weight by dominance ####
#splitting dominance index by 0.5
henDataLong[, RatioSplit := "Dom"]
henDataLong[Ratio < 0.5, RatioSplit := "Sub"]


#descriptives
henDataLong[, mean(weight), by = WoA]
descrBodyM = henDataLong[RatioSplit == "Sub" , .(Low = mean(weight)), by = .(WoA)]
temp =  henDataLong[RatioSplit == "Dom" , mean(weight), by = .(WoA)]
descrBodyM[, High := temp$V1]
descrBodyM[, Diff := High-Low]
henDataLong[, .(mean(weight), sd(weight)), by = .(RatioSplit, WoA)]


#weight gain
ggplot(data = henDataLong, aes(x = WoA, y = weight))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio), size = 3)+
  geom_line(aes(group = ID, colour = Ratio))+
  geom_line(data = henDataLong[, .(meanW = mean(weight)), by = WoA], aes(y = meanW), size = 1.5)+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'WoA', y = "Weight")+
  scale_color_gradient(low = "blue", high = "gold")

#weight gain
bodyMassFig = ggplot(data = henDataLong, aes(x = WoA, y = weight/1000))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio), size = 2.2, alpha = 0.8)+
  #geom_line(aes(group = ID, colour = Ratio), size = 0.2)+
  geom_line(data = henDataLong[RatioSplit == "Dom", .(meanW = mean(weight/1000)), by = WoA], 
            aes(y = meanW, linetype = ">0.5"), size = 1.2)+
  geom_line(data = henDataLong[RatioSplit == "Sub", .(meanW = mean(weight/1000)), by = WoA], 
            aes(y = meanW, linetype = "<0.5"), size = 1.2)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Body mass (kg)")+
  scale_x_continuous(breaks = c(25, 35, 45, 55),
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


  geom_point(aes(colour = Ratio), alpha = 0.5)+
  #geom_smooth(aes(colour = Ratio,group = as.factor(HenID)),se = F)+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = ">0.5"), size = 0.9, colour = "black", se = F)+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Time on top tier (h)")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))

varOfInterest[, RatioSplit2 := "< 0.5"]
varOfInterest[RatioSplit == "Dom", RatioSplit2 := "> 0.5"]



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
model.Weight = lmer(weight ~ WoA*Ratio+(1|Pen/HenID), data = henDataLong) #singularity
model.Weight = lmer(weight ~ WoA*Ratio+(1|HenID), data = henDataLong)
null.Weight =  lmer(weight ~ 1 +(1|HenID), data = henDataLong)
resid.Weight = simulateResiduals(model.Weight, 1000)
plot(resid.Weight) #looks like polynomial effect
plotResiduals(resid.Weight, form = henDataLong$Ratio)
plotResiduals(resid.Weight, form = henDataLong$WoA) #polynomial effect

model.Weight = lmer(weight ~ poly(WoA, 2)+ WoA:Ratio + Ratio +(1|HenID), data = henDataLong)
resid.Weight = simulateResiduals(model.Weight, 1000)
plot(resid.Weight) #good
plotResiduals(resid.Weight, form = henDataLong$Ratio)
plotResiduals(resid.Weight, form = henDataLong$WoA)
summary(model.Weight)
AIC(model.Weight, null.Weight)
parameters(model.Weight)
plot(allEffects(model.Weight))
#variance explained by fixed factors and entire model
r.squaredGLMM(model.Weight, null.Weight)

red.Weight = lmer(weight ~ poly(WoA, 2) +(1|HenID), data = henDataLong)
AIC(red.Weight, model.Weight)

#TODO: starting at which WoA is does a weight difference appear?
trend.Weight = emtrends(model.Weight, "WoA", va = "Ratio", at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))
summary(trend.Weight)
emmip(model.Weight, WoA ~ Ratio, cov.reduce = range, at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))
test.Weight = lmer(weight ~ poly(WoA, 2)+ WoA:Ratio + Ratio +(1|HenID), data = henDataLong[WoA <40,])


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
henData[Ratio < 0.17 , .(mean(gain), sd(gain))]
henData[Ratio > 0.7 , .(mean(gain), sd(gain))]

#correlations with comb size
#cor.test(henData$Comb, henData$weight_20)
# correlated -> R = 0.4, p = 0.015
cor.test(henData$Comb, henData$weight_26, method = "pearson", conf.level = 0.95)
# correlated -> R = 0.55, p <0.001
cor.test(henData$Comb, henData$weight_55, method = "pearson", conf.level = 0.95)
# correlated -> R = 0.54, p <0.001
#cor.test(henData$weight_55, henData$weight_20)
# correlated -> R = 0.5, p = 0.001
cor.test(henData$weight_55, henData$weight_26)
# correlated -> R = 0.82, p < 0.001


###### Health & social data ####

henDataLong[, .(mean(Severity), sd(Severity)), by = WoA]
henDataLong[, diffKBFSev := Severity - shift(Severity), by = ID][, mean(diffKBFSev), by = WoA]

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


kbfFig = ggplot(data = henDataLong, aes(x = WoA, y = Severity))+ 
  # geom_violin()+
  geom_point(aes(color = Ratio), size = 1.5, alpha = 0.8)+
  #geom_line(aes(group = ID, colour = Ratio), size = 0.2)+
  geom_line(data = henDataLong[RatioSplit == "Dom", .(meanKBF = mean(Severity)), by = WoA], 
            aes(y = meanKBF, linetype = ">0.5"), size = 0.9)+
  geom_line(data = henDataLong[RatioSplit == "Sub", .(meanKBF = mean(Severity)), by = WoA], 
            aes(y = meanKBF, linetype = "<0.5"), size = 0.9)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of Age', y = "KBF severity")+
  scale_x_continuous(breaks = c(25, 35, 45, 55),
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


hist(henDataLong$Severity) #normal but with a lot of zeros
model.KBF = lmer(Severity ~ Ratio*WoA +(1|Pen/HenID), data = henDataLong)#singularity
model.KBF = lmer(Severity ~ Ratio*WoA +(1|HenID), data = henDataLong)
null.KBF =  lmer(Severity ~ 1 +(1|HenID), data = henDataLong)
red.KBF =  lmer(Severity ~ WoA +(1|HenID), data = henDataLong)
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

AIC(model.KBF, null.KBF)
AIC(model.KBF, red.KBF)
parameters(red.KBF)
r.squaredGLMM(red.KBF)


henDataLong[, .(mean(feathers), sd(feathers)), by = WoA]
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

plumFig = ggplot(data = henDataLong, aes(x = WoA, y = feathers))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio), size = 1.5, alpha = 0.8)+
  #geom_line(aes(group = ID, colour = Ratio),size = 0.2)+
  geom_line(data = henDataLong[RatioSplit == "Dom", .(meanKBF = mean(feathers)), by = WoA], 
            aes(y = meanKBF, linetype = ">0.5"), size = 0.9)+
  geom_line(data = henDataLong[RatioSplit == "Sub", .(meanKBF = mean(feathers)), by = WoA], 
            aes(y = meanKBF, linetype = "<0.5"), size = 0.9)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Plumage condition")+
  scale_x_continuous(breaks = c(25, 35, 45, 55),
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", "Aggress. \nvalue")+
  scale_linetype_manual(name = NULL, 
                        values = c(">0.5" = "solid", "<0.5" = "dashed"),
                        ) +
  guides(color = guide_colorbar(order = 1), 
         linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


hist(henDataLong$feathers) 
# model.Feathers = glmer(feathers ~ Ratio*WoA+(1|Pen/HenID), data = henDataLong, family = poisson)
# resid.Feathers = simulateResiduals(model.Feathers, 1000)
# plot(resid.Feathers) 
# testOverdispersion(resid.Feathers)#overdispersed -> negative binomial
# model.Feathers = glmmTMB(feathers ~ Ratio + Ratio:WoA+ WoA +(1|Pen), data = henDataLong, family = "nbinom2") 
# resid.Feathers = simulateResiduals(model.Feathers, 1000) #looks like polynomial effect
# plot(resid.Feathers)
# plotResiduals(resid.Feathers, form = henDataLong$Ratio) #good
# plotResiduals(resid.Feathers, form = henDataLong$WoA) #polynomial

model.Feathers = glmmTMB(feathers ~ Ratio*WoA +(1|Pen/HenID), data = henDataLong, family = "nbinom2") 
resid.Feathers = simulateResiduals(model.Feathers, 1000) 
plot(resid.Feathers) #polynomial effect
plotResiduals(resid.Feathers, form = henDataLong$Ratio) #good
plotResiduals(resid.Feathers, form = henDataLong$WoA)#poly

model.Feathers = glmmTMB(feathers ~ Ratio + Ratio:WoA+ poly(WoA,2, raw = T) +(1|Pen/HenID), data = henDataLong, family = "nbinom2") 
resid.Feathers = simulateResiduals(model.Feathers, 1000) 
plot(resid.Feathers) #looks okay
plotResiduals(resid.Feathers, form = henDataLong$Ratio) #good
plotResiduals(resid.Feathers, form = henDataLong$WoA)#good
null.Feathers =  glmmTMB(feathers ~ 1 +(1|Pen/HenID), data = henDataLong, family = "nbinom2")
red.Feathers = glmmTMB(feathers ~ poly(WoA,2, raw = T) +(1|Pen/HenID), data = henDataLong, family = "nbinom2") 
summary(model.Feathers)
AIC(model.Feathers, null.Feathers)
AIC(model.Feathers, red.Feathers)
plot(allEffects(model.Feathers))
tidy(red.Feathers, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
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



###joint Fig

fig = ggarrange(bodyMassFig,kbfFig, plumFig,ncol = 3, labels = c("a)", "b)", "c)"), 
                font.label=list(color="black",size=18), common.legend = TRUE, legend = "right")
ggsave("PhysicCond.tiff", fig, "tiff", width = 30, height= 10, units = "cm", dpi = 300)




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

###### Comparison with controls##########

#comb size
dataControlcomb = rbind(henData[,.(ID, Pen, Comb)], combDataControls[,.(ID, Pen, Comb)])
dataControlcomb[1:36, Cat:= "Focal"]
dataControlcomb[37:.N, Cat:= "Control"]

ggplot(data = dataControlcomb, aes(y = Comb, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

model.Comb.Control = lmer(Comb ~ Cat +(1|Pen), data = dataControlcomb)
model.Comb.Control = lm(Comb ~ Cat , data = dataControlcomb)
resid.Comb = simulateResiduals(model.Comb.Control, 1000)
plot(resid.Comb) #heteroscedacity
# estimating the variance of y for different values of x
variance = lm(abs(model.Comb.Control$residuals) ~ model.Comb.Control$fitted.values)$fitted.values^2
# calculating the weights
weights = 1 / variance
# weighted regression model
weighted.model.Comb.Control = lm(Comb ~ Cat , data = dataControlcomb, weights = weights)
resid.Comb = simulateResiduals(weighted.model.Comb.Control, 1000)
plot(resid.Comb) #good
weighted.model.Comb.Control.null = lm(Comb ~ 1 , data = dataControlcomb, weights = weights)
AIC(weighted.model.Comb.Control.null, weighted.model.Comb.Control)
anova(weighted.model.Comb.Control.null, weighted.model.Comb.Control)
summary(weighted.model.Comb.Control)
parameters(weighted.model.Comb.Control)

#health
dataControlHealth = rbind(henDataLong[WoA == 55,.(ID, Pen, weight, feathers, wounds, comb, bare_head, footproblems, Severity)], 
                          healthDataControls[,.(ID, Pen, weight, feathers, wounds, comb, bare_head, footproblems, Severity)])
dataControlHealth[1:36, Cat:= "Focal"]
dataControlHealth[37:.N, Cat:= "Control"]
dataControlHealth[Cat == "Focal", weight := weight-15.6]


ggplot(data = dataControlHealth, aes(y = weight, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

ggplot(data = dataControlHealth, aes(y = Severity, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

ggplot(data = dataControlHealth, aes(y = feathers, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

#not necesssary to compare (near zero variation)
ggplot(data = dataControlHealth, aes(y = footproblems, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

hist(dataControlHealth$weight)
model.weight.Control = lmer(weight ~ Cat +(1|Pen), data = dataControlHealth)
model.weight.Control.null = lmer(weight ~ 1 +(1|Pen), data = dataControlHealth)
resid.Weight = simulateResiduals(model.weight.Control, 1000)
plot(resid.Weight) #good

AIC(model.weight.Control, model.weight.Control.null)
anova(model.weight.Control.null, model.weight.Control)

summary(model.weight.Control)
parameters(model.weight.Control)
#variance explained by fixed factors and entire model
r.squaredGLMM(model.weight.Control, model.weight.Control.null)

hist(dataControlHealth$Severity)
model.KBF.Control = lmer(Severity ~ Cat +(1|Pen), data = dataControlHealth)
model.KBF.Control = lm(Severity ~ Cat, data = dataControlHealth)
model.KBF.Control.null = lm(Severity ~ 1 , data = dataControlHealth)
resid.KBF = simulateResiduals(model.KBF.Control, 1000)
plot(resid.KBF) #good

AIC(model.KBF.Control, model.KBF.Control.null) #not better
anova( model.KBF.Control.null,model.KBF.Control)
summary(model.KBF.Control)


hist(dataControlHealth$feathers)
model.feathers.Control = lmer(feathers ~ Cat +(1|Pen), data = dataControlHealth)
model.feathers.Control.null = lmer(feathers ~ 1+(1|Pen) , data = dataControlHealth)
resid.feathers = simulateResiduals(model.feathers.Control, 1000)
plot(resid.feathers) #good

AIC(model.feathers.Control, model.feathers.Control.null)
anova(model.feathers.Control, model.feathers.Control.null)
summary(model.feathers.Control)
parameters(model.feathers.Control)

##### Tracking data processing ####

#relevant hens
hens = sort(unique(socialData$HenID))

trackingData = prepareTrackingData(trackingData, hens)
#include Week of age
trackingData = trackingData[tableWoA, on = "Date", nomatch = NULL]

#if tracking data preparation has been performed once:
# save tracking data
#fwrite(trackingData, "trackingData.csv", sep = ";")

# can be reloaded instead of running again
#trackingData <- fread("trackingData.csv")

# can be removed from workspace (large dataset)
#rm(trackingData)

##### Movement behaviour extraction ########

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

#TODO: do we need to be careful with late entries that the duration only lasts until 9 max? -> yes

#sift out only those in the morning, relevant for egg laying not resting (until 8)
#extract if hen was in nest zone on day or not
dailyNest = trackingData[Light == T & hour(Time) < 8, .(NestZone = ifelse(any(Zone == "Ramp_Nestbox"), 1, 0)), by = .(HenID, Date)][order(Date, HenID)]
#extract how long each hen was in the box
#first extract all nest entries
temp = trackingData[Light == T & hour(Time) < 8 & Zone == "Ramp_Nestbox", .(Time, Duration), by = .(HenID, Date)][order(Date, HenID)]
#make sure last entry is not longer in duration than until 8
temp[, dayLast := 0]
temp[temp[,.(rows = .I[.N]), by = .(HenID,Date)]$rows, dayLast := 1] #boolean to filter last entry by hen and day
temp[, durationEdge := 0]
temp[dayLast == 1, durationEdge:= ymd_hms(paste(Date, "08:00:00"))-Time]
temp[dayLast == 1 & durationEdge < Duration, Duration := durationEdge]

dailyNest[NestZone == 1, DurationNest := temp[,sum(Duration), by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[ , NestSleep := trackingData[LightIndic == T & hour(Time) < 8, Zone == "Ramp_Nestbox", 
                                       by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[NestZone == 0 | NestSleep == TRUE, DurationNest := 0]
#extract when median duration in the nest is reached
dailyNest[, MedDurNest := ifelse(NestZone == 1,  round(DurationNest/2),NA) ]
# median time point for half duration in nestbox
dailyNest[NestZone == 1, EntryNest := trackingData[Light == T & hour(Time) < 8 & Zone == "Ramp_Nestbox", 
                                                   Time[1], by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[NestZone == 1 & DurationNest > 0, MedTimeNest := EntryNest + MedDurNest]

#median time from lights on that day
dailyNest[, lightsOn := trackingData[Light == T & LightIndic == T, Time, 
                                     by = .(HenID, Date)][order(Date, HenID), Time]]
dailyNest[NestZone == 1 & DurationNest > 0, MedTimeNestLights := difftime(MedTimeNest, lightsOn, units="secs")]
dailyNest[, MedTimeNest := as.ITime(MedTimeNest)]

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


###### Feeder activity #######
#Feeder runs: (ab 22.11.: 2:00), 4:00, 6:00, 8:00, 10:00, 13:00, 15:00, 16:15

#apply function to extract the daily duration in the feeding zone during feeding runs
# TODO: first 5 min after start.. makes sense? or more?
dailyFeed = infeedZone(trackingData)
colnames(dailyFeed)[2] = "DurationFeed2"
colnames(dailyFeed)[3] = "DurationFeed4"
dailyFeed[, HenID := as.numeric(unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen))))]

#old feed reactivity
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

###### All movement behaviours #######
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
# sleep: max duration on top tier?
varOfInterest = varOfInterest[dailySleep, on = c(HenID = "HenID", Date = "NightCycle")]
#delete half nights (i.e., when only data of one date is available)
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


##### Analysing movement behaviour ######

###### descriptives ####
#mean durations in each functional zone across individuals and WoA
varOfInterest[, .(mean(Tier_4)/60/60, sd(Tier_4)/60/60)]
varOfInterest[, .(mean(Tier_2)/60/60, sd(Tier_2)/60/60)]
varOfInterest[, .(mean(Litter)/60/60, sd(Litter)/60/60)]
varOfInterest[, .(mean(Wintergarten)/60, sd(Wintergarten)/60)]
varOfInterest[, .(mean(Ramp_Nestbox)/60, sd(Ramp_Nestbox)/60)]

#most common daily zone across individuals and WoA 
varOfInterest[, .(mostCommonPercent = .N/nrow(varOfInterest)), by = MaxZone]

#duration on top tier for threshold split at 0.5 of aggression value
varOfInterest[RatioSplit == "Sub" , mean(Tier_4), by = .(HenID, WoA)][WoA< 27, mean(V1)/60/60]
varOfInterest[RatioSplit  == "Sub" , mean(Tier_4), by = .(HenID, WoA)][WoA> 49, mean(V1)/60/60]
varOfInterest[RatioSplit == "Dom" , mean(Tier_4), by = .(HenID, WoA)][WoA <27, mean(V1)/60/60]
varOfInterest[RatioSplit == "Dom" , mean(Tier_4), by = .(HenID, WoA)][WoA >49, mean(V1)/60/60]

descrTopTier = varOfInterest[RatioSplit == "Sub" , .(Low = mean(Tier_4)/60/60), by = .(WoA)]
temp =  varOfInterest[RatioSplit == "Dom" , .(High = mean(Tier_4)/60/60), by = .(WoA)]
descrTopTier[, High := temp$High]
descrTopTier[, Diff := High-Low]
descrTopTier[, DiffMin := Diff*60]


#duration in the litter for three lowest and three highest aggression value animals
varOfInterest[RatioSplit == "Sub" , mean(Litter), by = .(HenID, WoA)][WoA< 27, mean(V1)/60/60]
varOfInterest[RatioSplit == "Sub" , mean(Litter), by = .(HenID, WoA)][WoA> 49, mean(V1)/60/60]
varOfInterest[RatioSplit == "Dom" , mean(Litter), by = .(HenID, WoA)][WoA <27, mean(V1)/60/60]
varOfInterest[RatioSplit == "Dom" , mean(Litter), by = .(HenID, WoA)][WoA >49, mean(V1)/60/60]

#duration in feed zones during runs
varOfInterest[, .(mean(DurationFeed2)/60, sd(DurationFeed2)/60, mean(DurationFeed4)/60, sd(DurationFeed4)/60)]
varOfInterest[RatioSplit == "Sub" , .(mean(DurationFeed2),mean(DurationFeed4)), by = .(HenID, WoA)][WoA< 27, .(mean(V1)/60, mean(V2)/60)]
varOfInterest[RatioSplit == "Sub" , .(mean(DurationFeed2),mean(DurationFeed4)), by = .(HenID, WoA)][WoA> 49, .(mean(V1)/60, mean(V2)/60)]
varOfInterest[RatioSplit == "Dom" , .(mean(DurationFeed2),mean(DurationFeed4)), by = .(HenID, WoA)][WoA< 27, .(mean(V1)/60, mean(V2)/60)]
varOfInterest[RatioSplit == "Dom" , .(mean(DurationFeed2),mean(DurationFeed4)), by = .(HenID, WoA)][WoA> 49, .(mean(V1)/60, mean(V2)/60)]

descrFeedTier = varOfInterest[RatioSplit == "Sub" , .(Low2 = mean(DurationFeed2)/60,
                                               Low4 = mean(DurationFeed4)/60), by = .(WoA)]
temp =  varOfInterest[RatioSplit == "Dom" , .(High2 = mean(DurationFeed2)/60,
                                      High4 = mean(DurationFeed4)/60), by = .(WoA)]
descrFeedTier[, High2 := temp$High2]
descrFeedTier[, High4 := temp$High4]
descrFeedTier[, Diff2 := High2-Low2]
descrFeedTier[, Diff4 := High4-Low4]



###### 1 daily duration on top tier ####


hist(varOfInterest$Tier_4) 
hist(varOfInterest$Tier_4[varOfInterest$Tier_4 != 0]) 
#looks like poisson

model.Duration = glmmTMB(Tier_4 ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest, family = poisson)
resid.Duration = simulateResiduals(model.Duration, 1000)
plot(resid.Duration) #zero inflation issue most likely
plotResiduals(resid.Duration, form = varOfInterest$Ratio)
plotResiduals(resid.Duration, form = varOfInterest$WoA)
testZeroInflation(model.Duration) #indeed zero inflation
test.Duration = glmmTMB(Tier_4 ~ Ratio*WoA + (1|HenID), data = varOfInterest, family = poisson, 
                         zi = ~1)
AIC(model.Duration, test.Duration) # including general inflation term makes fit much better
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
plot(resid.Duration) #good
plotResiduals(resid.Duration, form = varOfInterest$Ratio)
plotResiduals(resid.Duration, form = varOfInterest$WoA)
parameters(split1.Duration, exp = TRUE) # the same
split1.null = glmer(on4 ~ 1 + (1|Pen/HenID), data = varOfInterest, family = binomial)
r.squaredGLMM(split1.Duration, split1.null)

split2.Duration = glmer(Tier_4 ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest[Tier_4 >0,], family = poisson)
resid.Duration = simulateResiduals(split2.Duration, 1000)
plot(resid.Duration) #deviation okay
plotResiduals(resid.Duration, form = varOfInterest[Tier_4 >0, Ratio])
plotResiduals(resid.Duration, form = varOfInterest[Tier_4 >0, WoA])
parameters(split2.Duration, exp = TRUE) # the same
split2.null = glmer(Tier_4 ~ 1 + (1|Pen/HenID), data = varOfInterest[Tier_4 >0,], family = poisson)
r.squaredGLMM(split2.Duration, split2.null)

trend.Duration = emtrends(test.Duration, "WoA", va = "Ratio", at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")
summary(trend.Duration)
emmip(test.Duration, WoA ~ Ratio, cov.reduce = range, at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")
emmip(test.Duration, WoA ~ Ratio, cov.reduce = range, at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))


plotData = as.data.table(emmeans(test.Duration, ~ pairwise ~ WoA*Ratio, type = "response",
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))$emmeans)
plotData

varOfInterest[, PredictDuration4 := predict(test.Duration, type = "response")]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = Tier_4), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = rate, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = lower.CL, ymax = upper.CL, group = Ratio, fill = as.factor(Ratio)), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)


###### 2 Vertical distance ####
#descriptives
varOfInterest[, .(mean(vertTravelDist), sd(vertTravelDist))]
varOfInterest[RatioSplit == "Sub" , mean(vertTravelDist), by = .(HenID, WoA)][WoA< 27, mean(V1)]
varOfInterest[RatioSplit == "Sub" , mean(vertTravelDist), by = .(HenID, WoA)][WoA> 49, mean(V1)]
varOfInterest[RatioSplit == "Dom" , mean(vertTravelDist), by = .(HenID, WoA)][WoA <27, mean(V1)]
varOfInterest[RatioSplit == "Dom" , mean(vertTravelDist), by = .(HenID, WoA)][WoA >49, mean(V1)]


hist(varOfInterest$vertTravelDist) #normally distributed with some outliers maybe
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
trend.vertDistance = emtrends(test.Travel, "Ratio", va = "WoA", 
                              at =  list(Ratio = unique(varOfInterest[Ratio < 0.5 | Ratio > 0.5, Ratio])), 
                              type = "response",
                              pbkrtest.limit = 6362)
trends = as.data.table(summary(trend.vertDistance))
trends[Ratio <0.5, mean(WoA.trend)]
trends[Ratio >0.5, mean(WoA.trend)]


varOfInterest[, PredictTravel := predict(test.Travel)]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = vertTravelDist), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio, fill = as.factor(Ratio)), alpha = 0.1)+
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

#print(VarCorr(test.Travel), comp = c("Variance", "Std.Dev."))
#VarCorr(test.Travel)$"HenID:Pen"[1] / (VarCorr(test.Travel)$"HenID:Pen"[1] + 
#                                          VarCorr(test.Travel)$"Pen"[1] + 
#                                          attr(VarCorr(test.Travel), "sc")^2)



###### 3 Nestbox Time #####

#only use data from when the lights went on at 02:00 plus 3 days of adjustment
#otherwise too affected -> only 15 days gone
varOfInterest[, fullCycle := !(day(Date) <24 & month(Date) == 11)]
nestData = varOfInterest[fullCycle == TRUE,]

#descriptives
nestData[WoA< 27, .(mean(as.numeric(MedTimeNestLights), na.rm = T)/3600, sd(as.numeric(MedTimeNestLights),na.rm = T)/3600)]
nestData[WoA> 49, .(mean(as.numeric(MedTimeNestLights), na.rm = T)/3600, sd(as.numeric(MedTimeNestLights),na.rm = T)/3600)]


hist(as.numeric(nestData[,MedTimeNestLights])) #normally distributed but a bit wide

model.Nest = lmer(as.numeric(MedTimeNestLights) ~ Ratio*WoA+ (1|Pen/HenID), data = nestData)
resid.Nest = simulateResiduals(model.Nest, 1000)
plot(resid.Nest) #polynomial necessary
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), Ratio])
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), WoA]) #poly

test.Nest = lmer(as.numeric(MedTimeNestLights) ~ Ratio+ Ratio:WoA +poly(WoA,2, raw = TRUE)+ (1|Pen/HenID), data = nestData)
null.Nest = lmer(as.numeric(MedTimeNestLights) ~ 1+ (1|Pen/HenID), data = nestData)
AIC(model.Nest, test.Nest) #poly better
AIC(test.Nest, null.Nest) #poly better
test.Nest2 = lmer(as.numeric(MedTimeNestLights) ~ poly(WoA,2, raw = TRUE)+Ratio+ (1|Pen/HenID), data = nestData)
#removing interaction does not improve fit unless raw = F 
AIC(test.Nest2, test.Nest)

summary(test.Nest)
parameters(test.Nest)
plot(allEffects(test.Nest))
r.squaredGLMM(test.Nest, null.Nest)

plotData = as.data.table(emmeans(test.Nest, ~ pairwise ~ Ratio+ Ratio:WoA +poly(WoA,2),
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)
plotData[, Time := as.ITime(response)/60/60]
plotData[, Plus := as.ITime(asymp.UCL)/60/60]
plotData[, Minus := as.ITime(asymp.LCL)/60/60]

varOfInterest[fullCycle == TRUE & !is.na(MedTimeNestLights), PredictNest:= predict(test.Nest)]
nestData[!is.na(MedTimeNestLights), PredictNest:= predict(test.Nest)/60/60]

ggplot()+ 
  geom_jitter(data = nestData, aes(x = WoA, y = as.numeric(MedTimeNestLights)/60/60), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = Time, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = Minus, ymax = Plus, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)+
  ylab("median nestbox time (in h)")


###### 4 Sleeping spot #####

hist(varOfInterest$onTop)
hist(varOfInterest[, sum(onTop), by = HenID][,V1])

#descriptives
varOfInterest[RatioSplit == "Sub" , mean(onTop), by = .(HenID, WoA)][WoA< 27, mean(V1)]
varOfInterest[RatioSplit == "Sub" , mean(onTop), by = .(HenID, WoA)][WoA> 49, mean(V1)]
varOfInterest[RatioSplit == "Dom" , mean(onTop), by = .(HenID, WoA)][WoA <27, mean(V1)]
varOfInterest[RatioSplit == "Dom" , mean(onTop), by = .(HenID, WoA)][WoA >49, mean(V1)]

model.Sleep = glmer(onTop ~ Ratio*WoA+ (1|Pen/HenID), data = varOfInterest, family = binomial)
#singularity due to Pen
model.Sleep = glmer(onTop ~ WoA*Ratio + (1|HenID),  data = varOfInterest, family = binomial)
null.Sleep = glmer(onTop ~ 1 + (1|HenID),  data = varOfInterest, family = binomial)
resid.Sleep = simulateResiduals(model.Sleep, 1000)
plot(resid.Sleep)#looks good
plotResiduals(resid.Sleep, form = varOfInterest$Ratio[!is.na(varOfInterest$onTop)])
plotResiduals(resid.Sleep, form = varOfInterest$WoA[!is.na(varOfInterest$onTop)])#possible poly
AIC(model.Sleep, null.Sleep) #better than null model

summary(model.Sleep)
parameters(model.Sleep, exp = TRUE) #TODO: why increased probability with age and WOA???
plot(allEffects(model.Sleep))
r.squaredGLMM(model.Sleep, null.Sleep)

plotData = as.data.table(emmeans(model.Sleep, ~ pairwise ~ WoA*Ratio, 
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)

exampleData = as.data.table(emmeans(model.Sleep, ~ pairwise ~ WoA*Ratio, 
                                    at =  list(Ratio = unique(varOfInterest[Ratio<0.17 | Ratio >0.7, Ratio]),
                                               WoA = c(20,21,22,23,24,25,50, 51, 52, 53, 54,55)), type = "response")$emmeans)

exampleData[Ratio <0.17 & WoA <26, mean(prob)]
exampleData[Ratio >0.17 & WoA <26, mean(prob)]

exampleData[Ratio <0.17 & WoA >26, mean(prob)]
exampleData[Ratio >0.17 & WoA >26, mean(prob)]


varOfInterest[, PredictSleep:= predict(model.Sleep, type = "response")]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = onTop), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = prob, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio, fill = as.factor(Ratio)), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)+
  ylab("Probability of sleeping on top tier")




###### 5 correlations between parameters ####


relVarPred = varOfInterest[, .(Duration_top = scale(PredictDuration4), 
                               verTravel= scale(PredictTravel, center = TRUE), 
                               Dur2MedNest = scale(PredictNest, center = TRUE), 
                               SleepTop = scale(PredictSleep))]

#correlation matrix

corrMatr = cor(relVarPred,  use = "pairwise.complete.obs") 

# Reorder the correlation matrix
corrMatr <- reorder_cormat(corrMatr)
upper_tri <- get_upper_tri(corrMatr)

plotCorr = reshape2::melt(upper_tri, na.rm = TRUE)
plotCorr = plotCorr[plotCorr$Var1 != plotCorr$Var2,]
levels(plotCorr$Var1) = c("Duration on \n top tier", 
                          "Sleeping on \n top tier", 
                          "Vertical travel \ndistance",
                          "Nestbox \ntiming")
levels(plotCorr$Var2) = c("Duration on \n top tier",
                           "Sleeping on \n top tier", 
                          "Vertical travel \n distance",
                          "Nestbox \ntiming")


ggplot(plotCorr, aes(x = Var2, y = Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson's\ncorrelation") +
  theme_minimal(base_size = 16)+ # minimal theme
  theme(panel.grid.major = element_blank())+
  coord_fixed()+
  scale_x_discrete(position = "top") +
  geom_text(aes(label = round(value, 2)), size = 5)+
  labs(x = "", y = "")



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
# 
# hist(as.numeric(varOfInterest$FeedReact))
# model.Feed = lmer(FeedReact ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest)
# model.Feed2 = lmer(FeedReact ~ Ratio*WoA + Ratio*MaxZone + MaxZone*WoA+ (1|Pen/HenID), data = varOfInterest)
# model.Feed3 = lmer(FeedReact ~ Ratio*WoA + MaxZone + (1|Pen/HenID), data = varOfInterest)
# anova(model.Feed, model.Feed2)
# anova(model.Feed2, model.Feed3)
# #TODO: Model with maxzone better but informative?? 
# resid.Feed = simulateResiduals(model.Feed, 1000)
# plot(resid.Feed) 
# plotResiduals(resid.Feed, form = varOfInterest$Ratio[!is.na(varOfInterest$FeedReact)])
# plotResiduals(resid.Feed, form = varOfInterest$WoA[!is.na(varOfInterest$FeedReact)])
# plotResiduals(resid.Feed, form = varOfInterest$MaxZone[!is.na(varOfInterest$FeedReact)])
# summary(model.Feed)
# plot(allEffects(model.Feed))
# parameters(model.Feed)
# 
# 
# hist(varOfInterest$FeedZoneDur)
# model.Feed = lmer(FeedZoneDur ~ Ratio*WoA+ (1|Pen/HenID), data = varOfInterest)
# model.Feed = lmer(FeedZoneDur ~ Ratio*WoA+ (1|HenID), data = varOfInterest)
# resid.Feed = simulateResiduals(model.Feed, 1000)
# plot(resid.Feed) 
# plotResiduals(resid.Feed, form = varOfInterest$Ratio[!is.na(varOfInterest$FeedReact)])
# plotResiduals(resid.Feed, form = varOfInterest$WoA[!is.na(varOfInterest$FeedReact)])
# summary(model.Feed)
# plot(allEffects(model.Feed))
# 
# 
# plotData = as.data.table(emmeans(model.Feed, ~ pairwise ~ WoA*Ratio, 
#                                  at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
#                                             WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)
# 
# 
# varOfInterest[, PredictFeedReact:= predict(model.Feed)]
# 
# ggplot()+ 
#   geom_jitter(data = varOfInterest, aes(x = WoA, y = FeedReact), height = 0.02, size = 1, alpha = 0.1)+
#   geom_line(data = plotData, aes(x = WoA, y = emmean, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
#   #geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio), alpha = 0.1)+
#   #facet_grid(.~Pen)+
#   theme_classic(base_size = 18)


##### Movement inspection plots ####

plotbreaks = as.IDate(c("2019-11-29", "2020-02-07", "2020-04-17", "2020-06-26"))

###### Durations in general ####

#heatmap of durations by ratio (hens sorted)
ggplot(trackingData[Light == T, .(Duration = sum(Duration)), by = .(HenID, Zone)], aes(x = Zone, y = factor(HenID, levels = socialData$HenID))) +
  geom_tile(aes(fill = Duration), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")

#heatmap of durations by ratio (hens sorted) by WoA
plot = trackingData[, .(nDays = length(unique(Date))), by = .(HenID, WoA)][
  trackingData[Light == T, .(Duration = sum(Duration)/60/60), by = .(HenID, Zone, WoA)][WoA %in% quantile(WoA),], 
  on = .(HenID, WoA)]
plot[, dailyDuration := Duration/nDays]
plot[, Zone := factor(Zone, levels = c("Wintergarten","Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]

ggplot(plot, aes(x = Zone, y = factor(HenID, levels = socialData$HenID))) +
  geom_tile(aes(fill = dailyDuration), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")+
  labs(y = "Hens by increasing aggression value")+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(.~WoA)

#most common zone daily
ggplot(varOfInterest, aes(x = Date, y = factor(MaxZone, levels = c("Wintergarten","Litter", "Tier_2", "Ramp_Nestbox", "Tier_4")), colour = Ratio))+
  geom_jitter(width = 0.1)+
  #geom_smooth(aes(group = as.factor(HenID)),se = F)+
  labs(x = 'Weeks of age', y = "Zone with the highest daily duration")+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  theme_classic(base_size = 18)+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_y_discrete(labels = c("Wintergarden","Litter", "Tier 1", "Nestbox tier", "Tier 3"))+
  theme(panel.background = element_rect(color = "black", size = 1))


ggplot(varOfInterest, aes(x = Date, y = factor(MaxZone, levels = c("Wintergarten","Litter", "Tier_2", "Ramp_Nestbox", "Tier_4")), colour = Highlight))+
  geom_jitter()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))

#Litter by Ratio
ggplot(varOfInterest, aes(x = Date, y = Litter, colour = Ratio))+
  geom_jitter(width = 0.1)+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_color_gradient(low = "blue", high = "gold")


# time spent on top during day
ggplot(varOfInterest, aes(x = Date, y = Tier_4))+
  geom_point(aes(colour = Ratio))+
  geom_line(varOfInterest[RatioSplit == "Dom", .(meanDur = mean(Tier_4)), by = Date], aes(y = meanDur))+
  scale_color_gradient(low = "blue", high = "gold")



topTierFig = ggplot(data = varOfInterest, aes(x = Date, y = Tier_4/60/60))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio), alpha = 0.5)+
  #geom_smooth(aes(colour = Ratio,group = as.factor(HenID)),se = F)+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = ">0.5"), size = 0.9, colour = "black", se = F)+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Time on top tier (h)")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))

varOfInterest[, RatioSplit2 := "< 0.5"]
varOfInterest[RatioSplit == "Dom", RatioSplit2 := "> 0.5"]


#time in feed zones during feed run
ggplot(varOfInterest[WoA %in% c(24,34, 44, 54),], aes(x = DurationFeed2/60, y = DurationFeed4/60))+
  geom_density2d_filled(contour_var = "ndensity")+
  theme_classic(base_size = 18)+
  facet_grid(RatioSplit2~WoA)+
  theme_classic(base_size = 18)+
  labs(x = 'Duration on the top tier (min)', y = "Duration on the lower tier (min)", fill = "density \nrange")+
  theme(legend.position = "top")
#scale_color_gradient(low = "blue", high = "gold")

#time in feed zones during feed run
varOfInterest[, ratioFeed := (DurationFeed4/60)-(DurationFeed2/60)]
ggplot(varOfInterest, aes(x = Date, y = ratioFeed))+
  geom_point(aes(colour = Ratio), alpha = 0.5)+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = ">0.5"), size = 0.9, colour = "black", se = F)+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Ratio: lower tier/top tier")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))
  
#geom_density2d_filled(contour_var = "ndensity")+
#  theme_classic(base_size = 18)+
#  facet_grid(RatioSplit2~WoA)+
#  theme_classic(base_size = 18)+
#  labs(x = 'Duration on Tier 1 (min)', y = "Duration on Tier 3 (min)", fill = "density \nrange")+
#  theme(legend.position = "top")
#scale_color_gradient(low = "blue", high = "gold")



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



###### travel distance #### 
hist(varOfInterest$vertTravelDist)
#progression of travel distance over days
ggplot(varOfInterest, aes(x = Date, y = vertTravelDist, colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))

travelDistFig = ggplot(varOfInterest, aes(x = Date, y = vertTravelDist))+
  geom_point(aes(colour = Ratio), alpha = 0.5)+
  #geom_smooth(aes(colour = Ratio,group = as.factor(HenID)),se = F)+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanTrav = mean(vertTravelDist)), by = Date], 
              aes(y = meanTrav, linetype = ">0.5"), size = 0.9, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanTrav = mean(vertTravelDist)), by = Date], 
              aes(y = meanTrav, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Vertical travel distance")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


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

#relation between time spent on top during day and sleeping spot
ggplot(varOfInterest, aes(x = Tier_4, y = onTop, colour = Ratio))+
  geom_point()+
  geom_smooth(color = "red", method = "glm", 
              method.args = list(family = "binomial"))+
  scale_color_gradient(low = "blue", high = "gold")

sleepFig = ggplot(varOfInterest, aes(x = Date, y = onTop, colour = Ratio))+
  geom_jitter(aes(colour = Ratio), height = 0.02, alpha = 0.5)+
  #geom_smooth(aes(colour = Ratio,group = as.factor(HenID)),se = F)+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanSleep = mean(onTop, na.rm = TRUE)), by = Date], 
              aes(y = meanSleep, linetype = ">0.5"), size = 0.9, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanSleep = mean(onTop, na.rm = TRUE)), by = Date], 
              aes(y = meanSleep, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Odds to sleep on top tier")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))




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


hist(as.numeric(varOfInterest$MedTimeNest))
#progression of time to enter the nestbox over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(MedTimeNest), colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

#progression of time to enter nestbox in duration from lights on
hist(as.numeric(varOfInterest$MedTimeNestLights))
ggplot(varOfInterest, aes(x = Date, y = as.numeric(MedTimeNestLights), colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_color_gradient(low = "blue", high = "gold")

nestFig = ggplot(nestData, aes(x = Date, y = as.numeric(MedTimeNestLights)/3600, colour = Ratio))+
  geom_point(aes(colour = Ratio), alpha = 0.5)+
  #geom_smooth(aes(colour = Ratio,group = as.factor(HenID)),se = F)+
  geom_smooth(data = nestData[RatioSplit == "Dom", .(meanNest = mean(as.numeric(MedTimeNestLights), na.rm = TRUE)/3600), by = Date], 
              aes(y = meanNest, linetype = ">0.5"), size = 0.9, se = F, colour = "black")+
  geom_smooth(data = nestData[RatioSplit == "Sub", .(meanNest = mean(as.numeric(MedTimeNestLights), na.rm = TRUE)/3600), by = Date], 
              aes(y = meanNest, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Median nestbox time (h)")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


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


###joint Fig

fig = ggarrange(topTierFig,travelDistFig, nestFig,sleepFig, nrow = 4, labels = c("a)", "b)", "c)", "d)"), 
                font.label=list(color="black",size=16), common.legend = TRUE, legend = "top", align = "hv")
ggsave("DailyMove.tiff", fig, "tiff", width = 30, height= 43, units = "cm", dpi = 300)


###### Feed reactivity ####


ggplot(data = varOfInterest, aes(x = Date, y = DurationFeed2))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio))+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanDur = mean(DurationFeed2)), by = Date], aes(y = meanDur), size = 1.5, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanDur = mean(DurationFeed2)), by = Date], aes(y = meanDur), size = 1.5, linetype = "dashed", colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Date', y = "DurationFeed2")+
  scale_color_gradient(low = "blue", high = "gold")

ggplot(data = varOfInterest, aes(x = Date, y = DurationFeed4))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio))+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanDur = mean(DurationFeed4)), by = Date], aes(y = meanDur), size = 1.5, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanDur = mean(DurationFeed4)), by = Date], aes(y = meanDur), size = 1.5, linetype = "dashed", colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Date', y = "DurationFeed4")+
  scale_color_gradient(low = "blue", high = "gold")


#total duration not in feed zones during run
ggplot(data = varOfInterest, aes(x = Date, y = NotFeedZone))+ 
  # geom_violin()+
  geom_point(aes(colour = Ratio))+
  geom_smooth(data = varOfInterest[RatioSplit == "Dom", .(meanDur = mean(NotFeedZone)), by = Date], aes(y = meanDur), size = 1.5, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[RatioSplit == "Sub", .(meanDur = mean(NotFeedZone)), by = Date], aes(y = meanDur), size = 1.5, linetype = "dashed", colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Date', y = "NotFeedZone")+
  scale_color_gradient(low = "blue", high = "gold")

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


#### Time series plots #########################
#relevant times
times = list(ymd_hms(c("2019-11-10 04:00:00", "2019-11-10 17:30:00")),
             ymd_hms(c("2020-01-10 02:00:00", "2020-01-10 17:30:00")),
             ymd_hms(c("2020-03-13 02:00:00", "2020-03-13 17:30:00")),
             ymd_hms(c("2020-05-08 02:00:00", "2020-05-08 17:30:00")),
             ymd_hms(c("2020-06-26 02:00:00", "2020-06-26 17:30:00")))


#two lowest ranking hens with really low severity -> not laying???
henDataLong[WoA == 55 & Severity < 2.5,]
test1 = trackingData[Hen == "Hen_108"& WoA == 54]
test2 = trackingData[Hen == "Hen_77"& WoA == 54]
test3 = nestData[Hen == "Hen_108"]
test4 = nestData[Hen == "Hen_77"]
test3 = nestData[Hen == "Hen_84"]
test4 = nestData[Hen == "Hen_30"]

plotData = test2
plotData[, Date := as.factor(as_date(Time))]
plotData[, Time_x := as_hms(Time)]
plotData[, Zone := factor(Zone, levels= c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
ggplot(plotData, aes(x = Time_x, y = Zone)) + 
  geom_step(group = 1) + 
  geom_point() + 
  labs(x = "time", y = "Zones") +
  facet_grid(Date~.)+
  theme_bw(base_size = 18)

plotData = test3
ggplot(plotData, aes(y = MedDurNest/60, x = Date)) + 
  geom_point() + 
  theme_bw(base_size = 18)

plotData = varOfInterest
plotData[, Highlight := 0]
plotData[Hen == "Hen_108", Highlight := 5]
plotData[Hen == "Hen_77", Highlight := 5]
plotData[, Highlight := as.factor(Highlight)]
plotData = henDataLong
plotData[, Highlight := 0]
plotData[Hen == "Hen_108", Highlight := 5]
plotData[Hen == "Hen_77", Highlight := 5]

ggplot(plotData, aes(x = WoA, y = feathers, colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)


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