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
#library(igraph)
library(emmeans)
library(MuMIn) #R^2 calculations
library(glmmTMB)
library(broom.mixed)
library(Hmisc) #correlation matrix
library(corrplot)#correlation plot
library(RColorBrewer)
library(ggpubr)#plot combination
library(sjPlot) #model parameters

#for repeatable results
set.seed(42)

##### Loading and preparing Data ###########
#load tracking data
trackingData <- fread("filtered_30s_full_combined.csv")
#load social observation data
observData = fread('socialData.csv')
#load comb sizes
combData = fread('focalSelection.csv')
combData[,ID := paste0(Pen,ID)]
combData[, Date := as_date(Date, format = "%d.%m.%Y")]
combData = combData[Date < as_date("2020-01-01")]
combData = combData[, .SD, .SDcols = c("Pen", "ID", "Date", "HenID", "TagID", "Mean5_7")]
colnames(combData)[6] <- "Comb" 

#add comb sizes of control birds for comparison to focal birds
combDataControls = fread("ControlsComb.csv")
combDataControls[,ID := paste0(Pen,Backpack)]
combDataControls = combDataControls[,.(Comb = mean(Area)), by = .(Pen, Backpack,ID)]

#prepare table containing the weeks of age of the animals (WoA) 
tableWoA = data.table(Date = seq(ymd("2019-10-11"), ymd("2020-07-02"), by = "day"),
                      WoA = rep(18:55, each = 7))

#load physical assessments data
healthData = fread('HA_all.csv')
healthData[,ID := paste0(pen,backpack)]
healthData[, feathers := neck + wings +tail + cloaca+ breast]
healthData[, footproblems:= r_podo + r_bumble +r_injure + l_podo + l_bumble + l_injure]
healthData[, date := as_date(date, format = "%d.%m.%Y")]
healthData = healthData[tableWoA, on = c(date = "Date"), nomatch = NULL]

#add data of control birds
healthDataControls = healthData[grepl('[0-9]', backpack),]
healthDataControls= healthDataControls[,.(ID, pen, date, weight, feathers, wounds, comb, bare_head, footproblems)]
healthDataControls[, Pen := pen]


#select only focals' physical assessments 
healthData = healthData[ID %in% unique(observData$ID),]
#change to wide format
healthDataWide = dcast(healthData, formula = ID ~ WoA, value.var = list("weight", "feathers", "wounds","comb", "bare_head", "footproblems"))


#load radiograph data
KBF = fread("KBF_scores.csv")
KBF[, Date := as_date(Date, format = "%d.%m.%Y")]
KBF = KBF[tableWoA, on = "Date", nomatch = NULL]
KBFControls = KBF[!is.na(no),]
KBFControls[, ID := paste0(HenID,no)]

#add radiograph data for controls
healthDataControls = KBFControls[healthDataControls, on = "ID"]

#select only focals' radiographs
KBF = KBF[HenID %in% unique(observData$ID),]
#change to wide format
KBFWide = dcast(KBF, formula = HenID ~ WoA, value.var = "Severity")
#change column names
colnames(KBFWide)[2:length(colnames(KBFWide))] = paste("KBF", colnames(KBFWide)[2:length(colnames(KBFWide))], sep = "_")

#join physical assessments with radiograph data
healthData = KBF[, .(WoA, HenID, Severity)][healthData[, .(WoA, ID, weight, feathers, wounds, comb, bare_head, footproblems)], on = c(HenID = "ID",WoA = "WoA")]
#chnage to wide format
healthDataWide = KBFWide[healthDataWide, on = c(HenID = "ID")]

#delete invalid rows from observations
observData = observData[!(Exclusion == 1),]

##### Social Data #####
###### Observer reliability ##### 

reliability = observData[Reliability == 1]

Observer1 = c(reliability$Sum_Sub[reliability$Observer == 'Masha'], 
              reliability$Sum_Dom[reliability$Observer == 'Masha'])
Observer2 = c(reliability$Sum_Sub[reliability$Observer == 'Tatiana'], 
              reliability$Sum_Dom[reliability$Observer == 'Tatiana'])

#check observer agreement using CCC
res = CCC(Observer1, Observer2, ci = "z-transform", conf.level = 0.95) 
#result: CCC = 0.89 [0.81, 0.94]

###### Social Data prep & Aggression value #####
#delete duplicates of reliability
socialData = observData[!(Observer == 'Masha' & Reliability == 1),]

#Descriptives
socialData[, RatioInteractions := Sum_Actions/Hen_in_sight*60]
#mean number of interactions per bird per minute
socialData[, .(MeanInter = mean(RatioInteractions), sdInter = sd(RatioInteractions),
               MedianInter = median(RatioInteractions), IQRInter = IQR(RatioInteractions))] 

#aggregate by individual
dataIndex = socialData[, .(Dom=sum(Sum_Dom), Sub=sum(Sum_Sub), Sum_Actions = sum(Sum_Actions),
                           Affil_given = sum(Affiliative_given),
                           Affil_rec = sum(Affiliative_received)),  by = .(ID, Pen)]

#median sum of interactions per bird per observation
dataIndex[, .(MedianSum = mean(Sum_Actions), sdInter = sd(Sum_Actions), min = min(Sum_Actions), max = max(Sum_Actions))] 

#calculate relative aggression value by dividing sum of dominant behaviours by total behaviours
dataIndex[, Aggr_value := Dom/(Dom+Sub)]

dataIndex = dataIndex[order(Aggr_value),]

#add comb sizes to the information
socialData = combData[dataIndex, on = "ID"]
socialData[, Hen := paste0("Hen_", HenID)]

socialData[,ID := factor(ID, levels = ID)]

#quick overview of relative aggression values across focals
ggplot(data = socialData, aes(x = ID, y = Aggr_value, colour = as.factor(Pen)))+ 
  geom_point(size = 6)+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  theme_classic(base_size = 18)


#create a wide and a long dataset with physical assessments and social information
henData = socialData[healthDataWide, on = c(ID = "HenID")]
henDataLong = socialData[healthData, on = c(ID = "HenID")]

henData[, Pen := as.factor(Pen)]

#splitting aggression value by 0.5 
# hens with aggression value > 0.5 showed more dominant behaviours than submissive
# hens with agg. value < 0.5 showed more submissive behaviours
henData[, Aggr_valueSplit := ifelse(Aggr_value >0.5,"> 0.5" , "< 0.5")]
henDataLong[, Aggr_valueSplit := ifelse(Aggr_value >0.5,"> 0.5" , "< 0.5")]

rm(dataIndex)

##### Comb Size #######

#descriptive
henData[, .(Mean = mean(Comb), SD = sd(Comb))]
henData[, .(Mean = mean(Comb), SD = sd(Comb)), by = Aggr_valueSplit]

#does comb size predict aggression value as hypothesised?

hist(henData$Comb, breaks = seq(min(henData$Comb), max(henData$Comb), 
                                length.out = 10)) 
# near normal with missing chunk in between

model.Comb = lmer(Aggr_value ~ Comb +(1|Pen), data = henData)
resid.Comb = simulateResiduals(model.Comb, 1000)
plot(resid.Comb) #good
plotResiduals(resid.Comb, form = henData$Comb) #good

#intercept only model
null.Comb =  lmer(Aggr_value ~ 1 +(1|Pen), data = henData)
AIC(model.Comb, null.Comb) # not better than intercept only model
anova(model.Comb, null.Comb) #but significantly better here...
#interpretation more careful

summary(model.Comb)
parameters(model.Comb)
#variance explained
r.squaredGLMM(model.Comb, null.Comb)

henData[, PredictAggr_value:= predict(model.Comb)]

#plot in supplementary
# Relationship of aggression value and comb size
combsizeFig = ggplot(data = henData, aes(y = Aggr_value, x = Comb, ))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  geom_line(aes(y = PredictAggr_value, colour = as.factor(Pen)), size = 1.5)+
  #facet_grid(.~Pen)+
  ylab("Aggression value")+
  xlab("Comb size (cm²)")+
  scale_colour_brewer(palette = "Paired", "Pen")+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", linewidth = 1))

ggsave("CombSize.tiff", combsizeFig, "tiff", width = 18, height= 14, units = "cm", dpi = 300)


##### Physical condition ####

###### Weight #####

#descriptives
henDataLong[, mean(weight), by = WoA]
descrBodyM = henDataLong[Aggr_valueSplit == "< 0.5" , .(Low = mean(weight)), by = .(WoA)]
temp =  henDataLong[Aggr_valueSplit == "> 0.5" , mean(weight), by = .(WoA)]
descrBodyM[, High := temp$V1]
descrBodyM[, Diff := High-Low]
henDataLong[, .(mean(weight), sd(weight)), by = .(Aggr_valueSplit, WoA)]
#Gain descriptive
henData[, gain := weight_55 - weight_20, by = HenID]
henData[, .(mean(gain), sd(gain)), by = Aggr_valueSplit]
henDataLong[, gain := weight - shift(weight), by= (HenID)]

#how does aggression value relate to weight over time?

hist(henDataLong$weight)#near-normal with small negative skew

model.Weight = lmer(weight ~ WoA*Aggr_value+(1|Pen/HenID), data = henDataLong) #singularity
#singularity -> remove Pen
model.Weight = lmer(weight ~ WoA*Aggr_value+(1|HenID), data = henDataLong)
resid.Weight = simulateResiduals(model.Weight, 1000)
plot(resid.Weight) #looks like polynomial effect
plotResiduals(resid.Weight, form = henDataLong$Aggr_value) # okay
plotResiduals(resid.Weight, form = henDataLong$WoA) #polynomial effect
#refit with WoA as polynomial factor
model.Weight = lmer(weight ~ poly(WoA, 2)+ WoA:Aggr_value + Aggr_value +(1|HenID), data = henDataLong)
resid.Weight = simulateResiduals(model.Weight, 1000)
plot(resid.Weight) #good
plotResiduals(resid.Weight, form = henDataLong$Aggr_value) #okay
plotResiduals(resid.Weight, form = henDataLong$WoA) #good

#intercept only model
null.Weight =  lmer(weight ~ 1 +(1|HenID), data = henDataLong)
AIC(model.Weight, null.Weight) # better than intercept model
#model without aggress. value
red.Weight = lmer(weight ~ poly(WoA, 2) +(1|HenID), data = henDataLong)
AIC(red.Weight, model.Weight) # model with aggr. value better

summary(model.Weight)
parameters(model.Weight)
plot(allEffects(model.Weight))
#variance explained by fixed factors and entire model
r.squaredGLMM(model.Weight, null.Weight)

#at which WoA is does a weight difference appear?
#post-hoc comparison
trend.Weight = emtrends(model.Weight, "WoA", va = "Aggr_value", at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))
summary(trend.Weight)
emmip(model.Weight, WoA ~ Aggr_value, cov.reduce = range, at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))
test.Weight = lmer(weight ~ poly(WoA, 2)+ WoA:Aggr_value + Aggr_value +(1|HenID), data = henDataLong[WoA <40,])


#plot in paper
# body mass over time in relation to aggression value
bodyMassFig = ggplot(data = henDataLong, aes(x = WoA, y = weight/1000))+ 
  # geom_violin()+
  geom_point(aes(colour = Aggr_value), size = 2.2, alpha = 0.8)+
  #geom_line(aes(group = ID, colour = Aggr_value), size = 0.2)+
  geom_line(data = henDataLong[Aggr_valueSplit == "> 0.5", .(meanW = mean(weight/1000)), by = WoA], 
            aes(y = meanW, linetype = ">0.5"), size = 1.2)+
  geom_line(data = henDataLong[Aggr_valueSplit == "< 0.5", .(meanW = mean(weight/1000)), by = WoA], 
            aes(y = meanW, linetype = "<0.5"), size = 1.2)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Body mass (kg)")+
  scale_x_continuous(breaks = c(25, 35, 45, 55),
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", linewidth = 1))



####### Correlation: weight x comb size #########

#correlations between weight measure during social observations and comb size
# and during weight at last measure and comb size
cor.test(henData$Comb, henData$weight_26, method = "pearson", conf.level = 0.95)
# correlated -> R = 0.55, p <0.001
cor.test(henData$Comb, henData$weight_55, method = "pearson", conf.level = 0.95)
# correlated -> R = 0.54, p <0.001



#are aggression value and physical condition related?

####### KBF severity #############

#descriptives 
henDataLong[, .(mean(Severity), sd(Severity)), by = WoA]
henDataLong[, diffKBFSev := Severity - shift(Severity), by = ID][, mean(diffKBFSev), by = WoA]
henDataLong[, .(mean(Severity), sd(Severity)), by = Aggr_valueSplit]

hist(henDataLong$Severity) #normal but with a lot of zeros -> zeroinflation?

model.KBF = lmer(Severity ~ Aggr_value*WoA +(1|Pen/HenID), data = henDataLong)
#singularity -> model without Pen
model.KBF = lmer(Severity ~ Aggr_value*WoA +(1|HenID), data = henDataLong)
resid.KBF = simulateResiduals(model.KBF, 1000)
plot(resid.KBF) #problem in left bottom corner 
plotResiduals(resid.KBF, form = henDataLong[, Aggr_value]) #okay
plotResiduals(resid.KBF, form = henDataLong[, WoA]) #zero-inflation?
testZeroInflation(resid.KBF) #zero-inflation problem
#check at what age we have less zeros in data
hist(henDataLong$Severity[henDataLong$WoA >22]) #still many
hist(henDataLong$Severity[henDataLong$WoA >26]) #maybe possible

#compare full model with model with less zeros (at later age)
#to compare estimates direction (to estimate impact of zero-inflation on results)
test.KBF = lmer(Severity ~ Aggr_value*WoA +(1|Pen), data = henDataLong[WoA > 26,])
resid.KBF = simulateResiduals(test.KBF, 1000)
plot(resid.KBF) #better
plotResiduals(resid.KBF, form = henDataLong[WoA >26, Aggr_value]) #okay
plotResiduals(resid.KBF, form = henDataLong[WoA >26, WoA]) #good

summary(model.KBF) #estimate direction: Aggr_value positive, WoA positive, interaction: negative
#tested against model with > 26WoA 
summary(test.KBF) # estimate direction the same -> stick to full model

#intercept only model
null.KBF =  lmer(Severity ~ 1 +(1|HenID), data = henDataLong)
AIC(model.KBF, null.KBF) #better than intercept only
#model without aggres. value
red.KBF =  lmer(Severity ~ WoA +(1|HenID), data = henDataLong)
AIC(model.KBF, red.KBF) #better without aggres. value

summary(red.KBF)
parameters(red.KBF)
#variance explained
r.squaredGLMM(red.KBF)

#plot used in paper
#keel bone fracture severity across time in relation to aggression value
kbfFig = ggplot(data = henDataLong, aes(x = WoA, y = Severity))+ 
  # geom_violin()+
  geom_point(aes(color = Aggr_value), size = 1.5, alpha = 0.8)+
  #geom_line(aes(group = ID, colour = Aggr_value), size = 0.2)+
  geom_line(data = henDataLong[Aggr_valueSplit == "> 0.5", .(meanKBF = mean(Severity)), by = WoA], 
            aes(y = meanKBF, linetype = ">0.5"), size = 0.9)+
  geom_line(data = henDataLong[Aggr_valueSplit == "< 0.5", .(meanKBF = mean(Severity)), by = WoA], 
            aes(y = meanKBF, linetype = "<0.5"), size = 0.9)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of Age', y = "KBF severity")+
  scale_x_continuous(breaks = c(25, 35, 45, 55),
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))



####### Feather condition ########

#descriptives
henDataLong[, .(mean(feathers), sd(feathers)), by = WoA]
henDataLong[, diffFeath := feathers - shift(feathers), by = ID][, mean(diffFeath), by = WoA]
henDataLong[, .(mean(feathers), sd(feathers)), by = Aggr_valueSplit]

hist(henDataLong$feathers) # many zeros!! try poisson

model.Feathers = glmer(feathers ~ Aggr_value*WoA+(1|Pen/HenID), data = henDataLong, family = poisson)
resid.Feathers = simulateResiduals(model.Feathers, 1000)
plot(resid.Feathers) # not great fit, bad residual fit
#try negative binomial
model.Feathers = glmmTMB(feathers ~ Aggr_value*WoA +(1|Pen/HenID), data = henDataLong, family = "nbinom2") 
resid.Feathers = simulateResiduals(model.Feathers, 1000) 
plot(resid.Feathers) #polynomial effect
plotResiduals(resid.Feathers, form = henDataLong$Aggr_value) #good
plotResiduals(resid.Feathers, form = henDataLong$WoA)#poly
#fit polynomial effect for WoA
model.Feathers = glmmTMB(feathers ~ Aggr_value + Aggr_value:WoA+ poly(WoA,2, raw = TRUE) +(1|Pen/HenID), data = henDataLong, family = "nbinom2") 
resid.Feathers = simulateResiduals(model.Feathers, 1000) 
plot(resid.Feathers) #looks okay
plotResiduals(resid.Feathers, form = henDataLong$Aggr_value) #good
plotResiduals(resid.Feathers, form = henDataLong$WoA)#good

#intercept only model
null.Feathers =  glmmTMB(feathers ~ 1 +(1|Pen/HenID), data = henDataLong, family = "nbinom2")
AIC(model.Feathers, null.Feathers) #better than intercept only
#model with only WoA effect
red.Feathers = glmmTMB(feathers ~ poly(WoA,2, raw = TRUE) +(1|Pen/HenID), data = henDataLong, family = "nbinom2") 
AIC(model.Feathers, red.Feathers) #only WoA slightly better

summary(red.Feathers)
tab_model(red.Feathers)
#variance explained
r.squaredGLMM(red.Feathers, null.Feathers)

#plot used in paper
#plumage condition across time in relation to aggression value
plumFig = ggplot(data = henDataLong, aes(x = WoA, y = feathers))+ 
  # geom_violin()+
  geom_point(aes(colour = Aggr_value), size = 1.5, alpha = 0.8)+
  #geom_line(aes(group = ID, colour = Aggr_value),size = 0.2)+
  geom_line(data = henDataLong[Aggr_valueSplit == "> 0.5", .(meanKBF = mean(feathers)), by = WoA], 
            aes(y = meanKBF, linetype = ">0.5"), size = 0.9)+
  geom_line(data = henDataLong[Aggr_valueSplit == "< 0.5", .(meanKBF = mean(feathers)), by = WoA], 
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


####### Foot problems #######

#descriptives
hist(henDataLong$footproblems) # extreme 0 inflation
median(henDataLong$footproblems) # not feasible to evaluate
ggplot(data = henDataLong, aes(x = WoA, y = footproblems,colour = Aggr_value))+ 
  # geom_violin()+
  geom_point(size = 3)+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs( x = 'WoA', y = "Foot problems")+
  scale_color_gradient(low = "blue", high = "gold")

####### Wounds ####### 
hist(henDataLong$wounds) # extreme 0 inflation
median(henDataLong$wounds) # not feasible to evaluate
ggplot(data = henDataLong, aes(x = WoA, y = wounds,colour = Aggr_value))+ 
  # geom_violin()+
  geom_point(size = 3)+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs( x = 'WoA', y = "Feather loss")+
  scale_color_gradient(low = "blue", high = "gold")


###joint Figure across physical conditions

fig = ggarrange(bodyMassFig,kbfFig, plumFig,ncol = 3, labels = c("a)", "b)", "c)"), 
                font.label=list(color="black",size=18), common.legend = TRUE, legend = "right")
ggsave(path = "plots", "PhysicCond.tiff", fig, "tiff", width = 30, height= 10, units = "cm", dpi = 300)


###### Comparison with controls##########

#comb size
dataControlcomb = rbind(henData[,.(ID, Pen, Comb)], combDataControls[,.(ID, Pen, Comb)])
dataControlcomb[1:36, Cat:= "Focal"]
dataControlcomb[37:.N, Cat:= "Control"]

hist(dataControlcomb$Comb) #near normal
model.Comb.Control = lmer(Comb ~ Cat +(1|Pen), data = dataControlcomb)
#singularity -> remove Pen
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
#null model
weighted.model.Comb.Control.null = lm(Comb ~ 1 , data = dataControlcomb, weights = weights)
AIC(weighted.model.Comb.Control.null, weighted.model.Comb.Control)
#slightly better than null
anova(weighted.model.Comb.Control.null, weighted.model.Comb.Control)
#but only a trend
summary(weighted.model.Comb.Control)
parameters(weighted.model.Comb.Control)

ggplot(data = dataControlcomb, aes(y = Comb, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", linewidth = 1))


#health
dataControlHealth = rbind(henDataLong[WoA == 55,.(ID, Pen, weight, feathers, wounds, comb, bare_head, footproblems, Severity)], 
                          healthDataControls[,.(ID, Pen, weight, feathers, wounds, comb, bare_head, footproblems, Severity)])
dataControlHealth[1:36, Cat:= "Focal"]
dataControlHealth[37:.N, Cat:= "Control"]
dataControlHealth[Cat == "Focal", weight := weight-15.6]

#weight
hist(dataControlHealth$weight)
model.weight.Control = lmer(weight ~ Cat +(1|Pen), data = dataControlHealth)
resid.Weight = simulateResiduals(model.weight.Control, 1000)
plot(resid.Weight) #good

#null model
model.weight.Control.null = lmer(weight ~ 1 +(1|Pen), data = dataControlHealth)
AIC(model.weight.Control, model.weight.Control.null) #better than null
anova(model.weight.Control.null, model.weight.Control) #significantly

summary(model.weight.Control)
parameters(model.weight.Control)
#variance explained by fixed factors and entire model
r.squaredGLMM(model.weight.Control, model.weight.Control.null)

ggplot(data = dataControlHealth, aes(y = weight, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))


#KBF severity
hist(dataControlHealth$Severity)
model.KBF.Control = lmer(Severity ~ Cat +(1|Pen), data = dataControlHealth)
#singularity -> remove Pen
model.KBF.Control = lm(Severity ~ Cat, data = dataControlHealth)
resid.KBF = simulateResiduals(model.KBF.Control, 1000)
plot(resid.KBF) #good

#null model
model.KBF.Control.null = lm(Severity ~ 1 , data = dataControlHealth)
AIC(model.KBF.Control, model.KBF.Control.null) #not better
anova( model.KBF.Control.null,model.KBF.Control)
summary(model.KBF.Control)

ggplot(data = dataControlHealth, aes(y = Severity, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

#feather condition
hist(dataControlHealth$feathers)
model.feathers.Control = lmer(feathers ~ Cat +(1|Pen), data = dataControlHealth)
resid.feathers = simulateResiduals(model.feathers.Control, 1000)
plot(resid.feathers) #good

#intrecpet only model
model.feathers.Control.null = lmer(feathers ~ 1+(1|Pen) , data = dataControlHealth)
AIC(model.feathers.Control, model.feathers.Control.null) #better than intercept
anova(model.feathers.Control, model.feathers.Control.null) #but not significantly
summary(model.feathers.Control)
parameters(model.feathers.Control)

ggplot(data = dataControlHealth, aes(y = feathers, x = Cat))+ 
  geom_boxplot()+
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  theme_classic(base_size = 18)+
  theme(panel.background = element_rect(color = "black", size = 1))

#wounds and footproblems not necesssary to compare (near zero variation)


##### Tracking data processing ####

#clear sorkspace and keep only necessary data
# List all objects in the workspace
all_objects = ls()
# Specify the object you want to keep
objects_to_keep = c("trackingData", "socialData", "tableWoA")

# Remove all objects except the one you want to keep
rm(list = setdiff(all_objects, objects_to_keep))

#load functions
source('helper_functions.R')
source('prepareTracking.R')

#relevant hens
hens = sort(unique(socialData$HenID))

trackingData = prepareTrackingData(trackingData, hens)
#include Week of age
trackingData = trackingData[tableWoA, on = "Date", nomatch = NULL]

#if tracking data prepaAggr_valuen has been performed once:
# save tracking data
#fwrite(trackingData, "trackingData.csv", sep = ";")

# can be reloaded instead of running again
#trackingData = fread("trackingData.csv")

# can be removed from workspace (large dataset)
#rm(trackingData)

##### Movement behaviour extraction ########

# calculate daily parameters per bird 

###### Vertical travel distance ############### 
# number of vertically crossed zones during light hours, divided by the seconds of the animals spent inside
#create vector where wintergarden doesn't exist (replaced by litter -> same vert. distance)
trackingData[, distZone := Zone]
trackingData[Zone == "Wintergarten", distZone := "Litter"]
#add zone vector shifted by one
trackingData[, nextZone := c(distZone[-1], NA), by = HenID]
#calculate distance travelled by using function defineDistance
trackingData[, distVertical := apply(X = cbind(distZone, nextZone), MARGIN = 1, FUN= defineDistance), by = HenID]
#daylight selection happens later when data is joined in section all parameters


###### Nestbox zone ########
#Nestbox entries per bird

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

#extract duration in the zone
dailyNest[NestZone == 1, durationNest := temp[,sum(Duration), by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[ , NestSleep := trackingData[LightIndic == T & hour(Time) < 8, Zone == "Ramp_Nestbox", 
                                       by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[NestZone == 0 | NestSleep == TRUE, durationNest := 0]
#extract when median duration in the nest is reached
dailyNest[, MedDurNest := ifelse(NestZone == 1,  round(durationNest/2),NA) ]
# median time point for half duration in nestbox
dailyNest[NestZone == 1, EntryNest := trackingData[Light == T & hour(Time) < 8 & Zone == "Ramp_Nestbox", 
                                                   Time[1], by = .(HenID, Date)][order(Date, HenID), V1]]
dailyNest[NestZone == 1 & durationNest > 0, MedTimeNest := EntryNest + MedDurNest]

#median time from lights on that day
dailyNest[, lightsOn := trackingData[Light == T & LightIndic == T, Time, 
                                     by = .(HenID, Date)][order(Date, HenID), Time]]
dailyNest[NestZone == 1 & durationNest > 0, MedTimeNestLights := difftime(MedTimeNest, lightsOn, units="secs")]
dailyNest[, MedTimeNest := as.ITime(MedTimeNest)]


###### Sleeping spot ############

#durations per zone per bird per day during dark hours
dailySleep = trackingData[Light == F, .(durationNight = sum(Duration)), by = .(HenID, Zone, NightCycle)][order(NightCycle, HenID)]
#extract maximum Zone for each bird per night and define as sleep zone
dailySleep = dailySleep[, .SD[which.max(durationNight)], by = .(HenID, NightCycle)]
#create Boolean indicating whether the zone was the highest tier or not
dailySleep[, onTop := ifelse(Zone == "Tier_4", 1, 0)]
colnames(dailySleep)[3] = "ZoneSleep"


###### Feeder activity #######
#Feeder runs: (ab 22.11.: 2:00), 4:00, 6:00, 8:00, 10:00, 13:00, 15:00, 16:15

#apply function to extract the daily duration in the feeding zone during feeding runs
# evaluates the first 5 min after start of the feed chain
dailyFeed = infeedZone(trackingData)
colnames(dailyFeed)[2] = "durationFeed2"
colnames(dailyFeed)[3] = "durationFeed4"
dailyFeed[, HenID := as.numeric(unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen))))]


###### Aggregating all movement behaviours #######
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
# nestbox: did hen go in nest?, total duration in nest, time point of median duration in nest 
varOfInterest = dailyNest[varOfInterest, on = c(HenID = "HenID", Date = "Date")]
# feeder zones: time spent in feeder zone during run versus outside run
varOfInterest = dailyFeed[varOfInterest, on = c(HenID = "HenID", Date = "Date")]

#add social information (dominance index)
varOfInterest = socialData[, .(HenID,Aggr_value,Comb)][varOfInterest, on = "HenID"] 

#make factor out of HenID & Pen with levels sorted in ascending Aggr_value order
varOfInterest[,HenID := factor(HenID, levels = socialData$HenID)]
varOfInterest[,Pen := factor(Pen)]
varOfInterest[,ZoneSleep := factor(ZoneSleep, levels = c("Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
varOfInterest[,MaxZone := factor(MaxZone, levels = c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
levels(varOfInterest$MaxZone) <- c("Wintergarden", "Litter", "Lower tier", "Nestbox tier", "Top tier")
#splitting dominance index by 0.5
varOfInterest[, Aggr_valueSplit := ifelse(Aggr_value >0.5,"> 0.5" , "< 0.5")]


##### Analysing movement behaviour ######

#clear workspace and keep only necessary data
# List all objects in the workspace
all_objects = ls()
# Specify the object you want to keep
objects_to_keep = c("varOfInterest")

# Remove all objects except the one you want to keep
rm(list = setdiff(all_objects, objects_to_keep))

#for consistent plotting: fox dates for x-axis
plotbreaks = as.IDate(c("2019-11-29", "2020-02-07", "2020-04-17", "2020-06-26"))


###### descriptives ####

#mean durations in each functional zone across individuals and WoA in h and min
varOfInterest[, .(mean(Tier_4)/60/60, sd(Tier_4)/60/60)]
varOfInterest[, .(mean(Tier_2)/60/60, sd(Tier_2)/60/60)]
varOfInterest[, .(mean(Litter)/60/60, sd(Litter)/60/60)]
varOfInterest[, .(mean(Wintergarten)/60, sd(Wintergarten)/60)]
varOfInterest[, .(mean(Ramp_Nestbox)/60, sd(Ramp_Nestbox)/60)]

#most common daily zone across individuals and WoA 
varOfInterest[, .(mostCommonPercent = .N/nrow(varOfInterest)), by = MaxZone]

#plot in paper
#most common zone daily
CommonDailyZone = ggplot(varOfInterest, aes(x = Date, y = MaxZone, colour = Aggr_value))+
  geom_jitter(width = 0.1)+
  #geom_smooth(aes(group = as.factor(HenID)),se = F)+
  labs(x = 'Weeks of age', y = "Zone with the highest daily duration")+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  theme_classic(base_size = 18)+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  #scale_y_discrete(labels = c("Wintergarden","Litter", "Tier 1", "Nestbox tier", "Tier 3"))+
  theme(panel.background = element_rect(color = "black", linewidth = 1))

ggsave(path = "plots", "CommonDailyZone.tiff", CommonDailyZone, "tiff", width = 35, height= 12, units = "cm", dpi = 300)

#duration on top tier for threshold split at 0.5 of aggression value
varOfInterest[Aggr_valueSplit == "< 0.5" , mean(Tier_4), by = .(HenID, WoA)][WoA< 27, .(M = mean(V1)/60/60, SD = sd(V1)/60/60)]
varOfInterest[Aggr_valueSplit  == "< 0.5" , mean(Tier_4), by = .(HenID, WoA)][WoA> 49, .(M = mean(V1)/60/60, SD = sd(V1)/60/60)]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(Tier_4), by = .(HenID, WoA)][WoA <27, .(M = mean(V1)/60/60, SD = sd(V1)/60/60)]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(Tier_4), by = .(HenID, WoA)][WoA >49, .(M = mean(V1)/60/60, SD = sd(V1)/60/60)]

descrTopTier = varOfInterest[Aggr_valueSplit == "< 0.5" , .(Low = mean(Tier_4)/60/60), by = .(WoA)]
temp =  varOfInterest[Aggr_valueSplit == "> 0.5" , .(High = mean(Tier_4)/60/60), by = .(WoA)]
descrTopTier[, High := temp$High]
descrTopTier[, Diff := High-Low]
descrTopTier[, DiffMin := Diff*60]

#duration in the litter by aggress value and time
varOfInterest[Aggr_valueSplit == "< 0.5" , mean(Litter), by = .(HenID, WoA)][WoA< 27, mean(V1)/60/60]
varOfInterest[Aggr_valueSplit == "< 0.5" , mean(Litter), by = .(HenID, WoA)][WoA> 49, mean(V1)/60/60]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(Litter), by = .(HenID, WoA)][WoA <27, mean(V1)/60/60]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(Litter), by = .(HenID, WoA)][WoA >49, mean(V1)/60/60]

#duration in feed zones during runs
varOfInterest[, .(mean(durationFeed2)/60, sd(durationFeed2)/60, mean(durationFeed4)/60, sd(durationFeed4)/60)]
varOfInterest[Aggr_valueSplit == "< 0.5" , .(mean(durationFeed2),mean(durationFeed4)), by = .(HenID, WoA)][WoA< 27, .(mean(V1)/60, mean(V2)/60)]
varOfInterest[Aggr_valueSplit == "< 0.5" , .(mean(durationFeed2),mean(durationFeed4)), by = .(HenID, WoA)][WoA> 49, .(mean(V1)/60, mean(V2)/60)]
varOfInterest[Aggr_valueSplit == "> 0.5" , .(mean(durationFeed2),mean(durationFeed4)), by = .(HenID, WoA)][WoA< 27, .(mean(V1)/60, mean(V2)/60)]
varOfInterest[Aggr_valueSplit == "> 0.5" , .(mean(durationFeed2),mean(durationFeed4)), by = .(HenID, WoA)][WoA> 49, .(mean(V1)/60, mean(V2)/60)]

descrFeedTier = varOfInterest[Aggr_valueSplit == "< 0.5" , .(Low2 = mean(durationFeed2)/60,
                                               Low4 = mean(durationFeed4)/60), by = .(WoA)]
temp =  varOfInterest[Aggr_valueSplit == "> 0.5" , .(High2 = mean(durationFeed2)/60,
                                      High4 = mean(durationFeed4)/60), by = .(WoA)]
descrFeedTier[, High2 := temp$High2]
descrFeedTier[, High4 := temp$High4]
descrFeedTier[, Diff2 := High2-Low2]
descrFeedTier[, Diff4 := High4-Low4]

#plot in paper
#time in feed zones during feed run
FeedRun = ggplot(varOfInterest[WoA %in% c(24,34, 44, 54),], aes(x = durationFeed2/60, y = durationFeed4/60))+
  geom_density2d_filled(contour_var = "ndensity")+
  theme_classic(base_size = 18)+
  facet_grid(Aggr_valueSplit~WoA)+
  theme_classic(base_size = 18)+
  labs(x = 'Duration on the top tier (min)', y = "Duration on the lower tier (min)", fill = "density \nrange")+
  theme(legend.position = "top")
#scale_color_gradient(low = "blue", high = "gold")

ggsave(path = "plots", "FeedRun.tiff", FeedRun, "tiff", width = 32, height= 18, units = "cm", dpi = 300)

#not used in paper
#duration on lower tier during feeder runs by aggress value and time
ggplot(data = varOfInterest, aes(x = Date, y = durationFeed2))+ 
  # geom_violin()+
  geom_point(aes(colour = Aggr_value))+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "> 0.5", .(meanDur = mean(durationFeed2)), by = Date], aes(y = meanDur), size = 1.5, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "< 0.5", .(meanDur = mean(durationFeed2)), by = Date], aes(y = meanDur), size = 1.5, linetype = "dashed", colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Date', y = "durationFeed2")+
  scale_color_gradient(low = "blue", high = "gold")

#not used in paper
#duration on top tier during feeder runs by aggress value and time
ggplot(data = varOfInterest, aes(x = Date, y = durationFeed4))+ 
  # geom_violin()+
  geom_point(aes(colour = Aggr_value))+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "> 0.5", .(meanDur = mean(durationFeed4)), by = Date], aes(y = meanDur), size = 1.5, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "< 0.5", .(meanDur = mean(durationFeed4)), by = Date], aes(y = meanDur), size = 1.5, linetype = "dashed", colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Date', y = "durationFeed4")+
  scale_color_gradient(low = "blue", high = "gold")


###### 1 daily duration on top tier ####

#do hens spend different time on the top tier over time depending on aggress. value?

hist(varOfInterest$Tier_4) # zero-inflation? 
hist(varOfInterest$Tier_4[varOfInterest$Tier_4 != 0]) 
#looks a bit like poisson

model.duration = glmmTMB(Tier_4 ~ Aggr_value*WoA + (1|Pen/HenID), data = varOfInterest, family = poisson)
resid.duration = simulateResiduals(model.duration, 1000)
plot(resid.duration) #zero inflation issue most likely
plotResiduals(resid.duration, form = varOfInterest$Aggr_value)
plotResiduals(resid.duration, form = varOfInterest$WoA)
testZeroInflation(model.duration) #indeed zero inflation
#account for zero inflation
test.duration = glmmTMB(Tier_4 ~ Aggr_value*WoA + (1|Pen/HenID), data = varOfInterest, family = poisson, 
                         zi = ~1)
AIC(test2.duration, test.duration) # including general inflation term makes fit much better
resid.duration = simulateResiduals(test.duration, 1000)
plot(resid.duration) #deviation okay 
plotResiduals(resid.duration, form = varOfInterest$Aggr_value)
plotResiduals(resid.duration, form = varOfInterest$WoA)

null.duration = glmmTMB(Tier_4 ~ 1 + (1|Pen/HenID), data = varOfInterest, family = poisson, 
                        zi = ~1)
AIC(test.duration, null.duration) #much better than null model
plot(allEffects(test.duration)) 
summary(test.duration)
tidy(test.duration, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
tab_model(test.duration)
r.squaredGLMM(test.duration, null.duration)

#compare against binomial and non-zero models if estimate direction is the same 
#binomial model
varOfInterest[, on4 := ifelse(Tier_4 >0, 1, 0)]
split1.duration = glmer(on4 ~ Aggr_value*WoA + (1|Pen/HenID), data = varOfInterest, family = binomial)
resid.duration = simulateResiduals(split1.duration, 1000)
plot(resid.duration) #good
plotResiduals(resid.duration, form = varOfInterest$Aggr_value)
plotResiduals(resid.duration, form = varOfInterest$WoA)

#non-zero model
split2.duration = glmer(Tier_4 ~ Aggr_value*WoA + (1|Pen/HenID), data = varOfInterest[Tier_4 >0,], family = poisson)
resid.duration = simulateResiduals(split2.duration, 1000)
plot(resid.duration) #deviation okay
plotResiduals(resid.duration, form = varOfInterest[Tier_4 >0, Aggr_value])
plotResiduals(resid.duration, form = varOfInterest[Tier_4 >0, WoA])

parameters(test.duration, exp = TRUE) # Aggr value: IRR > 1, WoA > 1, Aggr*WoA < 1
parameters(split1.duration, exp = TRUE) # the same
parameters(split2.duration, exp = TRUE) # the same

#save predicted values for correlation
varOfInterest[,PredictDuration4 := predict(test.duration)]

#extract development over time
trend.duration = emtrends(test.duration, "WoA", va = "Aggr_value", at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")
summary(trend.duration)
emmip(test.duration, WoA ~ Aggr_value, cov.reduce = range, at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")
emmip(test.duration, WoA ~ Aggr_value, cov.reduce = range, at =  list(WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))

#plot in paper
# daily duration on top tier over time depending on aggress value 
topTierFig = ggplot(data = varOfInterest, aes(x = Date, y = Tier_4/60/60))+ 
  # geom_violin()+
  geom_point(aes(colour = Aggr_value), alpha = 0.5)+
  #geom_smooth(aes(colour = Aggr_value,group = as.factor(HenID)),se = F)+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "> 0.5", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = ">0.5"), size = 0.9, colour = "black", se = F)+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "< 0.5", .(meanDur = mean(Tier_4)), by = Date], 
              aes(y = meanDur/60/60, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Time on top tier (h)")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", linewidth = 1))


###### 2 Vertical distance ####
#descriptives
varOfInterest[, .(mean(vertTravelDist), sd(vertTravelDist))]
varOfInterest[Aggr_valueSplit == "< 0.5" , mean(vertTravelDist), by = .(HenID, WoA)][WoA< 27, mean(V1)]
varOfInterest[Aggr_valueSplit == "< 0.5" , mean(vertTravelDist), by = .(HenID, WoA)][WoA> 49, mean(V1)]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(vertTravelDist), by = .(HenID, WoA)][WoA <27, mean(V1)]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(vertTravelDist), by = .(HenID, WoA)][WoA >49, mean(V1)]


hist(varOfInterest$vertTravelDist) #normally distributed with some outliers maybe
model.Travel = lmer(vertTravelDist ~ Aggr_value*WoA + (1|Pen/HenID), data = varOfInterest)
resid.Travel = simulateResiduals(model.Travel, 1000)
plot(resid.Travel) #deviation okay -> potential polynomial effect
plotResiduals(resid.Travel, form = varOfInterest$Aggr_value)
plotResiduals(resid.Travel, form = varOfInterest$WoA)#-> polynomial effect?

#test polynomial effect
test.Travel = lmer(vertTravelDist ~ Aggr_value:WoA + Aggr_value + poly(WoA,2) + (1|Pen/HenID), data = varOfInterest)
resid.Travel = simulateResiduals(test.Travel, 1000)
plot(resid.Travel) #looks better
plotResiduals(resid.Travel, form = varOfInterest$Aggr_value)
plotResiduals(resid.Travel, form = varOfInterest$WoA)#better
AIC(model.Travel, test.Travel)# poly is better

#intercept only
null.Travel = lmer(vertTravelDist ~ 1 + (1|Pen/HenID), data = varOfInterest)
AIC(test.Travel, null.Travel)#test is better
summary(test.Travel)
parameters(test.Travel)
plot(allEffects(test.Travel))
#variance explained
r.squaredGLMM(test.Travel, null.Travel)

#at what age does the difference arise?
#post hoc trends
plotData = as.data.table(emmeans(test.Travel, ~ pairwise ~ poly(WoA,2)+WoA:Aggr_value+ Aggr_value, 
                                 at =  list(Aggr_value = round(quantile(varOfInterest$Aggr_value), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))$emmeans)
#CAREFUL: long computation time
trend.vertDistance = emtrends(test.Travel, "Aggr_value", va = "WoA", 
                              at =  list(Aggr_value = unique(varOfInterest[Aggr_value < 0.5 | Aggr_value > 0.5, Aggr_value])), 
                              type = "response",
                              pbkrtest.limit = 6362)
trends = as.data.table(summary(trend.vertDistance))
trends[Aggr_value <0.5, mean(WoA.trend)]
trends[Aggr_value >0.5, mean(WoA.trend)]

#save predicted values for correlation
varOfInterest[,PredictVertTravelDist := predict(test.Travel)]


#not used in paper
ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = vertTravelDist), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = Aggr_value,  colour = as.factor(Aggr_value)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Aggr_value, fill = as.factor(Aggr_value)), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

#used in paper
#vertical travel distance over time by aggress. value
travelDistFig = ggplot(varOfInterest, aes(x = Date, y = vertTravelDist))+
  geom_point(aes(colour = Aggr_value), alpha = 0.5)+
  #geom_smooth(aes(colour = Aggr_value,group = as.factor(HenID)),se = F)+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "> 0.5", .(meanTrav = mean(vertTravelDist)), by = Date], 
              aes(y = meanTrav, linetype = ">0.5"), size = 0.9, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "< 0.5", .(meanTrav = mean(vertTravelDist)), by = Date], 
              aes(y = meanTrav, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Vertical travel distance")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


###### 3 Nestbox Time #####

#only use data from when the lights went on at 02:00 plus 3 days of adjustment
#otherwise affects nestbox time too much -> only 15 days gone
varOfInterest[, fullCycle := !(day(Date) <24 & month(Date) == 11)]
nestData = varOfInterest[fullCycle == TRUE,]

#descriptives
nestData[WoA< 27, .(mean(as.numeric(MedTimeNestLights), na.rm = T)/3600, sd(as.numeric(MedTimeNestLights),na.rm = T)/3600)]
nestData[WoA> 49, .(mean(as.numeric(MedTimeNestLights), na.rm = T)/3600, sd(as.numeric(MedTimeNestLights),na.rm = T)/3600)]

hist(as.numeric(nestData[,MedTimeNestLights])) #normally distributed but a bit wide

model.Nest = lmer(as.numeric(MedTimeNestLights) ~ Aggr_value*WoA+ (1|Pen/HenID), data = nestData)
resid.Nest = simulateResiduals(model.Nest, 1000)
plot(resid.Nest) #polynomial necessary
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), Aggr_value])
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), WoA]) #poly

#polynomial model
test.Nest = lmer(as.numeric(MedTimeNestLights) ~ Aggr_value+ Aggr_value:WoA +poly(WoA,2, raw = TRUE)+ (1|Pen/HenID), data = nestData)
resid.Nest = simulateResiduals(test.Nest, 1000)
plot(resid.Nest) #better
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), Aggr_value])
plotResiduals(resid.Nest, form = nestData[!is.na(MedTimeNest), WoA]) #good

#intercept only model
null.Nest = lmer(as.numeric(MedTimeNestLights) ~ 1+ (1|Pen/HenID), data = nestData)
AIC(model.Nest, test.Nest) #poly better
AIC(test.Nest, null.Nest) #poly better
test.Nest2 = lmer(as.numeric(MedTimeNestLights) ~ poly(WoA,2, raw = TRUE)+Aggr_value+ (1|Pen/HenID), data = nestData)
AIC(test.Nest2, test.Nest) # #removing interaction does not improve fit unless raw = F 

summary(test.Nest)
parameters(test.Nest)
plot(allEffects(test.Nest))
#variance epxplained
r.squaredGLMM(test.Nest, null.Nest)

#save predicted values for correlation
varOfInterest[fullCycle == TRUE & !is.na(MedTimeNest),PredictNestTime := predict(test.Nest)]

#inspect trends
#post hoc
plotData = as.data.table(emmeans(test.Nest, ~ pairwise ~ Aggr_value+ Aggr_value:WoA +poly(WoA,2),
                                 at =  list(Aggr_value = round(quantile(varOfInterest$Aggr_value), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)
plotData[, Time := as.ITime(response)/60/60]
plotData[, Plus := as.ITime(asymp.UCL)/60/60]
plotData[, Minus := as.ITime(asymp.LCL)/60/60]


#not used in paper
ggplot()+ 
  geom_jitter(data = nestData, aes(x = WoA, y = as.numeric(MedTimeNestLights)/60/60), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = Time, group = Aggr_value,  colour = as.factor(Aggr_value)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = Minus, ymax = Plus, group = Aggr_value), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)+
  ylab("median nestbox time (in h)")

#used in paper
# nest box time over time by aggress value
nestFig = ggplot(nestData, aes(x = Date, y = as.numeric(MedTimeNestLights)/3600, colour = Aggr_value))+
  geom_point(aes(colour = Aggr_value), alpha = 0.5)+
  #geom_smooth(aes(colour = Aggr_value,group = as.factor(HenID)),se = F)+
  geom_smooth(data = nestData[Aggr_valueSplit == "> 0.5", .(meanNest = mean(as.numeric(MedTimeNestLights), na.rm = TRUE)/3600), by = Date], 
              aes(y = meanNest, linetype = ">0.5"), size = 0.9, se = F, colour = "black")+
  geom_smooth(data = nestData[Aggr_valueSplit == "< 0.5", .(meanNest = mean(as.numeric(MedTimeNestLights), na.rm = TRUE)/3600), by = Date], 
              aes(y = meanNest, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Median nestbox time (h)")+
  scale_x_continuous(breaks = plotbreaks, 
                     labels = c("25", "35", "45" , "55"),
                     limits = as.IDate(c("2019-11-09", NA)))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


###### 4 Sleeping spot #####

#descriptives
varOfInterest[Aggr_valueSplit == "< 0.5" , mean(onTop), by = .(HenID, WoA)][WoA< 27, mean(V1)]
varOfInterest[Aggr_valueSplit == "< 0.5" , mean(onTop), by = .(HenID, WoA)][WoA> 49, mean(V1)]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(onTop), by = .(HenID, WoA)][WoA <27, mean(V1)]
varOfInterest[Aggr_valueSplit == "> 0.5" , mean(onTop), by = .(HenID, WoA)][WoA >49, mean(V1)]


hist(varOfInterest$onTop)
hist(varOfInterest[, sum(onTop), by = HenID][,V1])

model.Sleep = glmer(onTop ~ Aggr_value*WoA+ (1|Pen/HenID), data = varOfInterest, family = binomial)
#singularity due to Pen
model.Sleep = glmer(onTop ~ WoA*Aggr_value + (1|HenID),  data = varOfInterest, family = binomial)
resid.Sleep = simulateResiduals(model.Sleep, 1000)
plot(resid.Sleep)#looks good
plotResiduals(resid.Sleep, form = varOfInterest$Aggr_value[!is.na(varOfInterest$onTop)])
plotResiduals(resid.Sleep, form = varOfInterest$WoA[!is.na(varOfInterest$onTop)])#possible poly

#intercept only model
null.Sleep = glmer(onTop ~ 1 + (1|HenID),  data = varOfInterest, family = binomial)
AIC(model.Sleep, null.Sleep) #better than null model

summary(model.Sleep)
parameters(model.Sleep, exp = TRUE)
tab_model(model.Sleep)
plot(allEffects(model.Sleep))
#variance explained
r.squaredGLMM(model.Sleep, null.Sleep)

#save predicted values for correlation
varOfInterest[,PredictSleep := predict(model.Sleep)]

#extract means and trends post hoc
plotData = as.data.table(emmeans(model.Sleep, ~ pairwise ~ WoA*Aggr_value, 
                                 at =  list(Aggr_value = round(quantile(varOfInterest$Aggr_value), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)), type = "response")$emmeans)
# not used for paper
ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = onTop), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = prob, group = Aggr_value,  colour = as.factor(Aggr_value)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Aggr_value, fill = as.factor(Aggr_value)), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)+
  ylab("Probability of sleeping on top tier")

#used in the paper
# probability to sleep on the top tier over time y aggress. value 
sleepFig = ggplot(varOfInterest, aes(x = Date, y = onTop, colour = Aggr_value))+
  geom_jitter(aes(colour = Aggr_value), height = 0.02, alpha = 0.5)+
  #geom_smooth(aes(colour = Aggr_value,group = as.factor(HenID)),se = F)+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "> 0.5", .(meanSleep = mean(onTop, na.rm = TRUE)), by = Date], 
              aes(y = meanSleep, linetype = ">0.5"), size = 0.9, se = F, colour = "black")+
  geom_smooth(data = varOfInterest[Aggr_valueSplit == "< 0.5", .(meanSleep = mean(onTop, na.rm = TRUE)), by = Date], 
              aes(y = meanSleep, linetype = "<0.5"), size = 0.9, colour = "black", se = F)+
  theme_classic(base_size = 18)+
  labs(x = 'Weeks of age', y = "Odds to sleep on top tier")+
  scale_x_continuous(breaks = plotbreaks, #which(!duplicated(varOfInterest[,WoA]))
                     labels = c("25", "35", "45" , "55"))+
  scale_color_gradient(low = "blue", high = "gold", name = "Agress. \nvalue")+
  scale_linetype_manual(name = NULL, values = c(">0.5" = "solid", "<0.5" = "dashed")) +
  guides(color = guide_colorbar(order = 1), linetype = guide_legend(order = 2)) +
  theme(panel.background = element_rect(color = "black", size = 1))


###joint Fig

fig = ggarrange(topTierFig,travelDistFig, nestFig,sleepFig, nrow = 4, labels = c("a)", "b)", "c)", "d)"), 
                font.label=list(color="black",size=16), common.legend = TRUE, legend = "top", align = "hv")
ggsave(path = "plots", "DailyMove.tiff", fig, "tiff", width = 20, height= 43, units = "cm", dpi = 300)


###### 5 correlations between parameters ####

#scale values
relVarPred = varOfInterest[, .(duration_top = scale(PredictDuration4), 
                               verTravel= scale(PredictVertTravelDist, center = TRUE), 
                               Dur2MedNest = scale(PredictNestTime, center = TRUE), 
                               SleepTop = scale(PredictSleep))]

#correlation matrix
corrMatr = cor(relVarPred,  use = "pairwise.complete.obs") 

#load functions
source('helper_functions.R')

# Reorder the correlation matrix
corrMatr <- reorder_cormat(corrMatr)
upper_tri <- get_upper_tri(corrMatr)

#reshape matrix for plotting with ggplot
plotCorr = reshape2::melt(upper_tri, na.rm = TRUE)
plotCorr = plotCorr[plotCorr$Var1 != plotCorr$Var2,]
#renae levels for plot
levels(plotCorr$Var1) = c("Duration on \n top tier", 
                          "Sleeping on \n top tier", 
                          "Vertical travel \ndistance",
                          "Nestbox \ntiming")
levels(plotCorr$Var2) = c("Duration on \n top tier",
                           "Sleeping on \n top tier", 
                          "Vertical travel \n distance",
                          "Nestbox \ntiming")

#plot used in supplementary
# correlation matrix 
CorrPlot = ggplot(plotCorr, aes(x = Var2, y = Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson's\ncorrelation") +
  theme_minimal(base_size = 9)+ # minimal theme
  theme(panel.grid.major = element_blank())+
  coord_fixed()+
  scale_x_discrete(position = "top") +
  geom_text(aes(label = round(value, 2)), size = 3)+
  labs(x = "", y = "")

ggsave(path = "plots", "CorrPlot.tiff", CorrPlot, "tiff", width = 10, height= 8, units = "cm", dpi = 300)
