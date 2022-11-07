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

#load functions
source("functions.R")

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
#load health data
healthData = fread('HA_all.csv')
healthData[,ID := paste0(pen,backpack)]
healthData = healthData[(date == "09.12.2019" | date == "29.06.2020" | date == "30.06.2020")
                        & ID %in% unique(observData$ID),]
healthData[, feathers := neck + wings +tail + cloaca+ breast]
healthData[, injuries:= wounds + r_podo + r_bumble +r_injure + l_podo + l_bumble + l_injure]
healthData[, WoA := ifelse(date == "09.12.2019", 26, 55)]
healthDataWide = dcast(healthData, formula = ID ~ WoA, value.var = list("weight", "feathers", "injuries","comb", "bare_head"))

KBF = fread("KBF_scores.csv")
KBF = KBF[(Date == "09.12.2019" | Date == "29.06.2020" | Date == "30.06.2020")
                        & HenID %in% unique(observData$ID),]
KBF[, WoA := ifelse(Date == "09.12.2019", 26, 55)]
KBFWide = dcast(KBF, formula = HenID ~ WoA, value.var = "Severity")
colnames(KBFWide)[colnames(KBFWide) == "26"] = "KBF_26"
colnames(KBFWide)[colnames(KBFWide) == "55"] = "KBF_55"
healthData = KBF[, .(WoA, HenID, Severity)][healthData[, .(WoA, ID, WoA, weight, feathers, injuries, comb, bare_head)], on = c(HenID = "ID",WoA = "WoA")]
healthDataWide = KBFWide[healthDataWide, on = c(HenID = "ID")]

#delete invalid rows
observData = observData[!(Exclusion == 1),]


##### Observer reliability ##### 

reliability = observData[Reliability == 1]

Observer1 = c(reliability$Sum_Sub[reliability$Observer == 'Masha'], 
              reliability$Sum_Dom[reliability$Observer == 'Masha'])
Observer2 = c(reliability$Sum_Sub[reliability$Observer == 'Tatiana'], 
              reliability$Sum_Dom[reliability$Observer == 'Tatiana'])

res = CCC(Observer1, Observer2, ci = "z-transform", conf.level = 0.95) 
#result: CCC = 0.89 [0.8, 0.94]

lab <- paste("CCC: ", round(res$rho.c[,1], digits = 2), " (95% CI ", 
             round(res$rho.c[,2], digits = 2), " - ",
             round(res$rho.c[,3], digits = 2), ")", sep = "")
z <- lm(Observer1 ~ Observer2)

par(pty = "s")
plot(Observer1, Observer2, xlab = "Observer 1", 
     ylab = "Observer 2", pch = 16)
abline(a = 0, b = 1, lty = 2)
abline(z, lty = 1)
legend(x = "topleft", legend = c("Line of perfect concordance", 
                                 "Reduced major axis"), lty = c(2,1), lwd = c(1,1), bty = "n")
text(x = 1.55, y = 3.8, labels = lab)

tmp.mean <- mean(res$blalt$delta)
tmp.sd <- sqrt(var(res$blalt$delta))

plot(res$blalt$mean, res$blalt$delta, pch = 16, 
     xlab = "Average PEFR ", 
     ylab = "Difference in PEFR") 
abline(h = tmp.mean, lty = 1, col = "gray")
abline(h = tmp.mean - (2 * tmp.sd), lty = 2, col = "gray")
abline(h = tmp.mean + (2 * tmp.sd), lty = 2, col = "gray")
legend(x = "topleft", legend = c("Mean difference", 
                                 "Mean difference +/ 2SD"), lty = c(1,2), bty = "n")
legend(x = 0, y = 125, legend = c("Difference"), pch = 16, 
       bty = "n")


##### Social Data prep & Index #####
socialData = observData[!(Observer == 'Masha' & Reliability == 1),]

#Descriptives
socialData[, ratioInteractions := Sum_Actions/Hen_in_sight]
#mean number of interactions per bird per minute
socialData[, .(MeanInter = mean(ratioInteractions), sdInter = sd(ratioInteractions))] 

#aggregate by individual
dataIndex = socialData[, .(Dom=sum(Sum_Dom), Sub=sum(Sum_Sub), Sum_Actions = sum(Sum_Actions),
                           Affil_given = sum(Affiliative_given),
                           Affil_rec = sum(Affiliative_received)),  by = .(ID, Pen)]

#median sum of interactions per bird per observation
dataIndex[, .(MedianSum = mean(Sum_Actions), sdInter = sd(Sum_Actions), min = min(Sum_Actions), max = max(Sum_Actions))] 

dataIndex[, Ratio := Dom/(Dom+Sub)]

dataIndex = dataIndex[order(Ratio),]

socialData = combData[dataIndex, on = "ID"]

ggplot(data = socialData, aes(x = ID, y = Ratio, colour = as.factor(Pen)))+ 
  geom_point(size = 6)+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  theme_classic(base_size = 18)


henData = socialData[healthDataWide, on = c(ID = "HenID")]
henDataLong = socialData[healthData, on = c(ID = "HenID")]

##### Badges of Status #######

#comb size by dominance
ggplot(data = henData, aes(x = Comb, y = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'Comb Size', y = "Dominance ratio")

#weight at 26 by dominance
ggplot(data = henData, aes(x = weight_26, y = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  geom_point(aes(x = weight_55), size = 3, colour = "grey")+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  geom_smooth(aes(x = weight_55), method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'Weight', y = "Dominance ratio")

#weight at 55 by dominance
ggplot(data = henData, aes(x = Comb, y = weight_55))+ 
  # geom_violin()+
  geom_point(size = 6, aes(colour = Ratio))+
  #geom_hline(yintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'Comb Size', y = "Weight")

cor.test(henData$Comb, henData$weight_26)
# correlated -> R = 0.56, p <0.001
cor.test(henData$Comb, henData$weight_55)
# correlated -> R = 0.54, p <0.001
cor.test(henData$weight_55, henData$weight_26)
# correlated -> R = 0.82, p <0.001

model.Comb = lmer(Ratio ~ Comb +(1|Pen), data = henData)
null.Comb =  lmer(Ratio ~ 1 +(1|Pen), data = socialData)
resid.Comb = simulateResiduals(model.Comb, 1000)
plot(resid.Comb)
plotResiduals(resid.Comb, form = socialData$Mean5_7)
summary(model.Comb)
anova(model.Comb, null.Comb)
parameters(model.Comb)

henData[, PredictRatio:= predict(model.Comb)]
ggplot(data = henData, aes(y = Ratio, x = Comb, ))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 2)+
  geom_line(aes(y = PredictRatio, colour = as.factor(Pen)))+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

##### Weight gain by dominance ####

#weight gain
ggplot(data = henDataLong, aes(x = as.factor(WoA), y = weight))+ 
  # geom_violin()+
  geom_point(size = 6, aes(colour = Ratio))+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'WoA', y = "Weight")

henData[, gain := weight_55-weight_26]

#gain by dominance
ggplot(data = henData, aes(y = gain, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'weight gain', y = 'Weight gain', x = "Dominance ratio")


model.Gain = lmer(gain ~ Ratio +(1|Pen), data = henData)
null.Gain =  lmer(gain ~ 1 +(1|Pen), data = henData)
resid.Gain = simulateResiduals(model.Gain, 1000)
plot(resid.Gain)
plotResiduals(resid.Gain, form = henData$Ratio)
summary(model.Gain)
anova(model.Gain, null.Gain)
parameters(model.Gain)

henData[, PredictGain:= predict(model.Gain)]
ggplot(data = henData, aes(x = Ratio, y = gain, ))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 2)+
  geom_line(aes(y = PredictGain, colour = as.factor(Pen)), size = 1.5)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)


##### Health & social data ####

#KBF and dominance
ggplot(data = henData, aes(y = KBF_55, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Health & Dominance', y = 'KBF Severity', x = "Dominance ratio")

#KBF development
ggplot(data = henDataLong, aes(x = as.factor(WoA), y = Severity))+ 
  # geom_violin()+
  geom_point(size = 6, aes(colour = Ratio))+
  geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'WoA', y = "KBF Severity")
#two lowest ranking hens with really low severity -> not laying???

model.KBF = lmer(KBF_55 ~ Ratio +(1|Pen), data = henData)
model.KBF = lm(KBF_55 ~ Ratio, data = henData)
null.KBF =  lm(KBF_55 ~ 1, data = henData)
resid.KBF = simulateResiduals(model.KBF, 1000)
plot(resid.KBF)
plotResiduals(resid.KBF, form = henData$Ratio)
summary(model.KBF)
anova(model.KBF, null.KBF)
parameters(model.KBF)

#feather cover
ggplot(data = henData, aes(y = feathers_55, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Health & Dominance', y = 'Feather cover', x = "Dominance ratio")

#feather development
ggplot(data = henDataLong, aes(x = as.factor(WoA), y = feathers))+ 
  # geom_violin()+
  geom_point(size = 6, aes(colour = Ratio))+
  #geom_text(aes(label = ID),hjust=1, vjust=1)+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'WoA', y = "Feather loss")

model.Feathers = lmer(feathers_55 ~ Ratio +(1|Pen), data = henData)
null.Feathers =  lmer(feathers_55 ~ 1 +(1|Pen), data = henData)
resid.Feathers = simulateResiduals(model.Feathers, 1000)
plot(resid.Feathers)
plotResiduals(resid.Feathers, form = henData$Ratio)
summary(model.Feathers)
anova(model.Feathers, null.Feathers)
parameters(model.Feathers)


#injuries 
ggplot(data = henData, aes(y = injuries_55, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Health & Dominance', y = 'injuries', x = "Dominance ratio")

#comb pecks
ggplot(data = henData, aes(y = comb_55, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Health & Dominance', y = 'injuries', x = "Dominance ratio")

#bare head
ggplot(data = henData, aes(y = bare_head_55, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Health & Dominance', y = 'injuries', x = "Dominance ratio")

##### Gentle beak pecking behaviour##########

ggplot(data = henData, aes(y = Affil_given, x = Ratio))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 2)+
  geom_smooth(method = lm)+
  theme_classic(base_size = 18)

ggplot(data = henData, aes(y = Affil_rec, x = Ratio))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 2)+
  geom_smooth(method = lm)+
  theme_classic(base_size = 18)

##### Extraction tracking data #########

#add numerical hen identifier to tracking data
trackingData[, HenID := unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen)))]
trackingData[, HenID := as.numeric(HenID)]
trackingData[, Time := ymd_hms(Time)]

#relevant hens
hens = sort(unique(socialData$HenID))

#splitting data into Hens
splitHen = splitHenData(trackingData)

hen_list <- vector(mode='list', length=length(hens))
start = ymd_hms("2019-11-11 00:00:00")
end = ymd_hms("2019-12-20 23:59:59")
i = 1
for (hen in hens) {
    hen_list[[i]] = extractInterval(splitHen[[hen]], start, end)
    i = i+1
} 
trackingRel = mergeHenData(hen_list)
#add Pen
trackingRel[, Pen := as.numeric((unlist(regmatches(PackID, gregexpr('\\(?[0-9,.]+', PackID)))))]
#add date grouping
trackingRel[, Date:= as_date(Time)]
#fwrite(trackingRel, file = "relTrackingData.csv", sep = ";")
rm(trackingData)

########### data checks #############################

#there is data for every day?
length(unique(trackingRel$Date)) == 40
#which day is missing?
unique(trackingRel$Date)
#-> 05.12. is missing
#check start and end time of each day by pen
dayPenEntries = trackingRel[, .(min = min(Time), max = max(Time)), by = .(Date, Pen)]
#-> missing data on 04.12. -06.12. -> exclude
#-> health assessment on 09.12. -> exclude
trackingRel = trackingRel[!(Date == as.Date("2019-12-04")|
                              Date == as.Date("2019-12-06")|
                            Date == as.Date("2019-12-09")),]
length(unique(trackingRel$Date)) #-> now we have 36 days in total of data


#all individuals have data?
length(unique(trackingRel$HenID)) == 36
#all individuals have data every day?
entriesPerDay = trackingRel[, length(unique(Date)), by = HenID]
entriesPerDay[V1 != 36,] #Hen 1 (12), 65(32) and 84(34) don't have every day
#which days are missing?
days = unique(trackingRel$Date)

trackingRel[HenID == 1,unique(Date)] #no data between the 22.11. and 19.12.
trackingRel[HenID == 65,unique(Date)] #no data between 2.12. and 08.12.
trackingRel[HenID == 84,unique(Date)] #no data between 07.12. and 08.12.


#exclude HenID 1 and exclude December days before battery replacement
trackingRel = trackingRel[!(Date == as.Date("2019-12-02")|
                              Date == as.Date("2019-12-03")|
                              Date == as.Date("2019-12-07")|
                            Date == as.Date("2019-12-08")|
                              HenID == 1),]
entriesPerDay = trackingRel[, length(unique(Date)), by = HenID]
length(unique(trackingRel$Date)) #-> now we have 32 days in total of data


########### time series plot examples #########################
#relevant times
times = list(ymd_hms(c("2019-11-12 03:30:00", "2019-11-12 17:30:00")),
             ymd_hms(c("2019-11-21 03:30:00", "2019-11-21 17:30:00")),
             ymd_hms(c("2019-11-30 02:30:00", "2019-11-30 17:30:00")),
             ymd_hms(c("2019-12-08 02:30:00", "2019-12-08 17:30:00")),
             ymd_hms(c("2019-12-20 02:30:00", "2019-12-20 17:30:00")))

hen_list <- vector(mode='list', length=length(hens))
j = 1
for (hen in hens) {
  hen_list[[j]] <- vector(mode='list', length=5)
  for (i in 1:5){
    hen_list[[j]][[i]] = extractInterval(splitHen[[hen]], times[[i]][1], times[[i]][2])
  }
  j= j+1
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

plotData = trackingFull
plotData[, Date := as.factor(Date)]
plotData[, Time_x := as_hms(Time)]
plotData[, Zone := factor(Zone, levels= c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]

#ideas for plots: for each individual an average step function plot
#requires: full sequence for every day

splitHen = splitHenData(trackingFull)
hens = unique(trackingFull$HenID)
average_hen <- vector(mode='list', length=length(hens))
i = 1
for (i in 1:length(hens)) {
  cat("Hen", i)
  list = fillSeqHenData(splitHen[[i]])
  #cat("done")
  average_hen[[i]] = averageDailyZone(list)
} 


#make overview plots per hen
for (i in 1:length(hens)){
  plotData = average_hen[[i]]
  plotData[, Zone := factor(Zone, levels= c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]
  ggplot(plotData, aes(x = Time, y = Zone)) + 
    geom_step(group = 1) + 
    geom_point() + 
    labs(x = "time", y = "Zones") +
    ggtitle(hens[i])+
    theme_bw(base_size = 18)
  
  
  ggsave(filename = paste0(hens[i],".png"), plot = last_plot(), 
         width = 40, height = 18, dpi = 300, units = "cm")
}

########## Data preparations ################

#data cleaning step
#remove Wintergarten entries outside time that Wintergarten can be open
#TODO: ask Laura for dates open and times
#delete all Wintergarten entries at night (not possible, tracking system errors)
trackingRel = trackingRel[!(hour(Time) > 16 & Zone == "Wintergarten"), ]
#delete all Wintergarten entries before 9 in the morning (not possible, tracking system error)
trackingRel = trackingRel[!(hour(Time) < 9 & Zone == "Wintergarten"), ]

#add beginning and end of each day by looping through every hen and every day and adding the stamps
# (from original data set (splitHen))
hen_list <- vector(mode='list', length=length(unique(trackingRel$hens)))
splitHen = splitHenData(trackingRel)
hens = unique(trackingRel$HenID)
days = as.character(unique(trackingRel$Date))
j = 1
for (hen in hens) {
  hen_list[[j]] <- vector(mode='list', length=length(days))
  i = 1
  for (day in days){
    start = ymd_hms(paste(day, "00:00:00"))
    end = ymd_hms(paste(day, "23:59:59"))
    hen_list[[j]][[i]] = extractInterval(splitHen[[j]], start, end)
    i = i+1
  }
  hen_list[[j]] = mergeHenData(hen_list[[j]])
  j= j+1
} 
trackingFull = mergeHenData(hen_list)

#remove unnecessary data from workspace
rm(splitHen)

#make sure date is correct (now with inserting of end and start of day)
trackingFull[, Date:= as_date(Time)]
# mark true transitions
trackingFull[, TrueTransition := TRUE]
trackingFull[, DayIndic := F]
trackingFull[hour(Time) == 0 & minute(Time) == 0 & second(Time) == 0, TrueTransition := FALSE] 
trackingFull[hour(Time) == 23 & minute(Time) == 59 & second(Time) == 59, TrueTransition := FALSE] 
trackingFull[hour(Time) == 0 & minute(Time) == 0 & second(Time) == 0, DayIndic := TRUE] 
trackingFull[hour(Time) == 23 & minute(Time) == 59 & second(Time) == 59, DayIndic := TRUE] 


#include light-dark cycle
trackingFull[, Light := TRUE]
trackingFull[hour(Time) > 17, Light := FALSE]
trackingFull[month(Time)== 11 & day(Time) < 15 & hour(Time) < 4, Light := FALSE]
trackingFull[month(Time)== 11 & day(Time) < 22 & day(Time) > 14 & hour(Time) < 3, Light := FALSE]
trackingFull[(month(Time)== 12| (month(Time)== 11 & day(Time) > 21)) & hour(Time) < 2, Light := FALSE]

#add beginning and end of the night
#split series into hens
hen_list = splitHenData(trackingFull)
#create empty list
full_list = vector(mode='list', length= length(hen_list))
for (i in 1:length(unique(trackingFull$HenID))){
  full_list[[i]] <- vector(mode='list', length=length(days))
  
  j = 1
  for (day in days) {
    relDat = hen_list[[i]][Date == day,]

    #to keep track if an entry already exists at lights on
    start = F
    
    if( month(day) == 11 & day(day) < 15){
      if(any(hour(relDat$Time) == 4 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
        start = T
      } else{
        dayStart = "04:00:00"
      }
      
    } else if (month(day)== 11 & day(day) < 22 & day(day) >14){
              if(any(hour(relDat$Time) == 3 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
                  start = T
                } else{
                  dayStart = "03:00:00"
                }
    } else{
      if(any(hour(relDat$Time) == 2 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
        start = T
      } else{
      dayStart = "02:00:00"
      }
    }
    if(start == F){
      #use first indication of lights as index
      indexStart = which(relDat$Light == T)[1]
      startEntry = relDat[indexStart-1,]
      startEntry[, Time := ymd_hms(paste(day, dayStart))]
      startEntry[, Light := T]
      startEntry[, TrueTransition := F]
      startEntry[, DayIndic := F]
      relDat = insertRow(relDat, indexStart, startEntry)
    }
    if(!any(hour(relDat$Time) == 17 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
      dayEnd = "17:00:00"
      #use last entry before lights out time
      indexEnd = which(hour(relDat$Time) > 17)[1]
      endEntry = relDat[indexEnd-1,]
      endEntry[, Time := ymd_hms(paste(day, dayEnd))]
      endEntry[, Light := F]
      endEntry[, TrueTransition := F]
      endEntry[, DayIndic := F]
      relDat = insertRow(relDat, indexEnd, endEntry)
    }
    full_list[[i]][[j]] = relDat
    j =j+1
  }
  full_list[[i]] = mergeHenData(full_list[[i]])
}
trackingFull = mergeHenData(full_list)

#add duration
trackingFull[, Duration := (shift(Time, type="lead") - Time), by = HenID]
#set all last day entries to 1 sec so that missing days don't count into the duration
trackingFull[hour(Time)== 23 & minute(Time)== 59, Duration := 1]


####### Parameters ########

#transitions per bird
transitions = trackingFull[TrueTransition == TRUE, .(Transitions = .N), by = .(HenID)]
#transitions per bird per day
transDaily = trackingFull[TrueTransition == TRUE, .(Transitions = .N), by = .(HenID, Date)]
#transitions per bird per day (light hours)
transDailyL = trackingFull[TrueTransition == TRUE & Light == T, .(Transitions = .N), by = .(HenID, Date)]
#transitions per bird per day (dark hours)
transDailyD = trackingFull[TrueTransition == TRUE & Light == F, .(Transitions = .N), by = .(HenID, Date)]


#durations per zone per bird
durations = trackingFull[, .(Duration = sum(Duration)), by = .(HenID, Zone)]
#durations per zone per bird per day
durDaily = trackingFull[, .(Duration = sum(Duration)), by = .(HenID, Zone, Date)]



###### vertical travel distance ############### 
# number of vertically crossed zones during light hours, divided by the seconds of the animals spent inside
#create vector where wintergarden doesn't exist (replaced by litter)
trackingFull[, distZone := Zone]
trackingFull[Zone == "Wintergarten", distZone := "Litter"]
#add zone vector shifted by one
trackingFull[, nextZone := c(distZone[-1], NA), by = HenID]
#calculate distance travelled by using function defineDistance
trackingFull[, distVertical := apply(X = cbind(distZone, nextZone), MARGIN = 1, FUN= defineDistance), by = HenID]

###### sleeping spot ############
#indicate which night times belong together
trackingFull[Light == T, NightCycle := Date]
trackingFull[hour(Time) > 17|hour(Time) == 17 , NightCycle := Date]
trackingFull[hour(Time) < 4 & Light == F, NightCycle := Date-1]
#durations per zone per bird per day during dark hours
durDailyD = trackingFull[Light == F, .(DurationNight = sum(Duration)), by = .(HenID, Zone, NightCycle)]
durDailyD = durDailyD[order(HenID),]
#extract maximum Zone for each bird per day
mainDurDailyD = durDailyD[, .SD[which.max(DurationNight)], by = .(HenID, NightCycle)]
mainDurDailyD[, onTop := ifelse(Zone == "Tier_4", 1, 0)]


###### wintergarden use #####
#all Wintergarten entries per bird
dailyGarten = trackingFull[Zone == "Wintergarten", .(Duration = Duration), by = .(HenID, Time, Date)]
#extract if hen goes out on day or not
inGarten = trackingFull[, .(Out = ifelse(any(Zone == "Wintergarten"), 1, 0)), by = .(HenID, Date)]
#extract how long each hen went out per day
durDailyGarten = dailyGarten[, .(DurationGarten = sum(Duration)), by = .(HenID, Date)]
# latency to go out
latGarten = dailyGarten[, .(LatencyGarten = Time[1] - ymd_hms(paste(as_date(Time[1]), "09:00:00"))), by = .(HenID, Date)]

###### Time in nestbox zone ########
#all Nestbox entries per bird
dailyNest = trackingFull[Light == T & Zone == "Ramp_Nestbox", .(Duration = Duration), by = .(HenID, Time, Date)]
#extract if hen was in nest zone on day or not
inNest = trackingFull[Light == T, .(NestZone = ifelse(any(Zone == "Ramp_Nestbox"), 1, 0)), by = .(HenID, Date)]
#extract how long each hen went out per day
durDailyNest = dailyNest[, .(DurationNest = sum(Duration), MedDurNest = sum(Duration)/2), by = .(HenID, Date)]
# add duration and medDur to daily nest zone entries for latency calculations
dailyNest = durDailyNest[dailyNest, on = c("HenID", "Date")]
# median time point for half duration in nestbox
timeNest = dailyNest[, .(TimeNest = Time[1] + MedDurNest), by = .(HenID, Date)]
timeNest[, del := duplicated(TimeNest), by = .(HenID, Date)]
timeNest = timeNest[del == F, .(HenID, Date, TimeNest)]


###### Feeder reactivity #######
#TODO: ask Laura for times of feeder run


#create data.table containing all daily measures
varOfInterest = trackingFull[Light == T, .(vertTravelDist = sum(distVertical)), by = .(HenID, Date, Pen)]
#careful with join not to loose the hens for hens who don't go out for example
#full outer join
varOfInterest = varOfInterest[mainDurDailyD, on = c(HenID = "HenID", Date = "NightCycle")]
varOfInterest = varOfInterest[!is.na(Pen),]
varOfInterest = inGarten[varOfInterest, on = c(HenID = "HenID", Date = "Date")]
varOfInterest = durDailyGarten[varOfInterest, on = c(HenID = "HenID", Date = "Date")]
varOfInterest[is.na(DurationGarten), DurationGarten := 0]
varOfInterest = latGarten[varOfInterest, on = c(HenID = "HenID", Date = "Date")]
varOfInterest = durDailyNest[varOfInterest, on = c(HenID = "HenID", Date = "Date")]
varOfInterest[is.na(DurationNest), DurationNest := 0]
varOfInterest[is.na(MedDurNest), MedDurNest := 0]
varOfInterest = timeNest[varOfInterest, on = c(HenID = "HenID", Date = "Date")]

#add socialInformation
varOfInterest = socialData[, .(HenID,Ratio)][varOfInterest, on = "HenID"] 


##### Parameter inspection plots ####

#Highlight three most extreme
varOfInterest[, Highlight := "Any"]
varOfInterest[HenID == 82 | HenID == 97 | HenID == 33, Highlight := "Dom"]
varOfInterest[HenID == 77 | HenID == 5 | HenID == 108, Highlight := "Sub"]
#make factor out of HenID with levels sorted in ascending Ratio order
varOfInterest[,HenID := factor(HenID, levels = socialData$HenID)]


###### sleeping spot ####
hist(varOfInterest$onTop)
hist(varOfInterest[, sum(onTop), by = HenID][,V1])
#heatmap
ggplot(varOfInterest, aes(x = Zone, y = as.factor(HenID))) +
  geom_tile(aes(fill = as.numeric(DurationNight)), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")
#progression of most common sleeping spot on Top or not
ggplot(varOfInterest, aes(x = Date, y = onTop, colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = HenID),method = "glm", 
              method.args = list(family = "binomial"),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))


###### travel distance #### 
hist(varOfInterest$vertTravelDist)
#progression of travel distance over days
ggplot(varOfInterest, aes(x = Date, y = vertTravelDist, colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))

###### Wintergarten ####

#Went out or not?
hist(varOfInterest$Out)
hist(varOfInterest[, sum(Out), by = HenID][,V1])
ggplot(varOfInterest, aes(x = Date, y = Out, colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = HenID),method = "glm", 
              method.args = list(family = "binomial"),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))

#latency to enter
hist(as.numeric(varOfInterest$LatencyGarten))
#progression of latency to enter the garden over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(LatencyGarten), colour = Highlight))+
  geom_point()+
  geom_smooth(data = varOfInterest[Highlight != "Any",], aes(group = HenID),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))

#duration in the wintergarten
hist(as.numeric(varOfInterest$DurationGarten))
#progression of latency to enter the garden over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(DurationGarten), colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))

###### Nestbox zone ####

#latency to enter
varOfInterest[, TimeNestPure := as.ITime(TimeNest)]
hist(as.numeric(varOfInterest$TimeNestPure))
#progression of time to enter the nestbox over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(TimeNestPure), colour = Ratio))+
  geom_point()+
  geom_smooth(data = varOfInterest[Highlight != "Any",], aes(group = HenID),se = F)#+
  #scale_colour_manual(values = c("grey", "red", "blue"))


#duration in the nestbox zone
hist(as.numeric(varOfInterest$DurationNest))
#progression of duration in nestbox over days
ggplot(varOfInterest, aes(x = Date, y = as.numeric(DurationNest), colour = Highlight))+
  geom_point()+
  geom_smooth(aes(group = HenID),se = F)+
  scale_colour_manual(values = c("grey", "red", "blue"))



##### Test models ######

varOfInterest[, RatioSplit := "Dom"]
varOfInterest[Ratio < 0.5, RatioSplit := "Sub"]
varOfInterest = varOfInterest[order(Date),]
varOfInterest[,  DateID := as.numeric(Date)]

###### Vertical distance ####
model.Travel = lmer(vertTravelDist ~ Ratio + DateID + (1|Pen/HenID), data = varOfInterest)
resid.Travel = simulateResiduals(model.Travel, 1000)
plot(resid.Travel)
plotResiduals(resid.Travel, form = varOfInterest$Ratio)
plotResiduals(resid.Travel, form = varOfInterest$DateID)
summary(model.Travel)

varOfInterest[, PredictTravel := predict(model.Travel)]

ggplot(varOfInterest, aes(x = Date, y = vertTravelDist, color = RatioSplit)) +
  geom_jitter(size=2) + 
  geom_line(data = varOfInterest[, mean(PredictTravel), by = .(Date, RatioSplit)], 
                              aes(x = Date, y =V1, colour = RatioSplit),size=1.5) + 
  theme_classic(base_size = 18)+ 
  ylab("daily vertical travel distance")#+ 

#plot individual variation
ggplot(data = varOfInterest, aes(x = Date, y = PredictTravel, colour = as.factor(HenID))) + 
  geom_line(aes(group = HenID), size=1)+
  labs(y = "Predicted vertical travel distance",color = "Hen ID")+ 
  theme_classic(base_size = 18)+ theme(legend.position="none")+ 
  guides(color = guide_legend(nrow = 4))

#plot against actual data of some examples
ggplot(data = varOfInterest[Highlight != "Any",], 
       aes(x = Date, colour = as.factor(HenID))) + 
  geom_jitter(aes(y = vertTravelDist))+
  geom_line(aes(y = PredictTravel, group = HenID), size=1)+
  labs(y = "Predicted vertical travel distance",color = "Hen ID")+ 
  theme_classic(base_size = 18)+ 
  guides(color = guide_legend(nrow = 4))

#divide the variance explained by animal_id by the total phenotypic variance 
#(animal_id + month:year + year + residual variance)
print(VarCorr(model.Travel), comp = c("Variance", "Std.Dev."))
VarCorr(model.Travel)$"HenID:Pen"[1] / (VarCorr(model.Travel)$"HenID:Pen"[1] + 
                                VarCorr(model.Travel)$"Pen"[1] + 
                                attr(VarCorr(model.Travel), "sc")^2)
#-> within -individual repeatability of 0.31 (within a pen)

#to get uncertainity estimate: simulate model 1000 times
set.seed(1) 
library(arm)
simulated <- sim(model.Travel, n.sim = 1000)
posterior_HenID <- apply(simulated@ranef$"HenID:Pen"[ , , 1],1,var) 
posterior_Pen <- apply(simulated@ranef$"Pen"[ , , 1],1,var) 
posterior_residual <- simulated@sigma^2
quantile(posterior_HenID / (posterior_HenID + posterior_Pen + posterior_residual), prob=c(0.025, 0.5, 0.975))

#coeffeicient of variation for beteween individual variance
CVi <- sqrt(posterior_HenID) / summary(model.Travel)$coefficients[1] 
quantile(CVi,prob=c(0.025, 0.5, 0.975)) #TODO: how can this be negative? what went wrong?

#behavioural type
#the behavioral type is the best linear unbiased prediction (BLUP) of the 
#random effect, i.e. the prediction for mean behavioral expression for each individual.
library(merTools)
randomSims <- as.data.table(REsim(model.Travel, n.sims = 1000)) 
head(randomSims[randomSims$groupFctr=="HenID:Pen",])
# add the dominance ratio of the individual 
randomSims[, HenID := as.factor(unlist(strsplit(groupID,':'))[2*(1:dim(randomSims)[1])-1])] 
randomSims <- merge(randomSims[randomSims$groupFctr=="HenID:Pen",], 
                    varOfInterest[!duplicated(HenID),c("HenID","Ratio", "Pen", "RatioSplit")])
# add identifier to color individuals uniquely 
randomSims[, ID := ifelse(HenID %in% c(77, 5, 108, 82, 97, 33), as.character(HenID), "Other individuals")]

#add population intercept for easier interpretation on transition scale
#TODO: why all negative????
randomSims[, meanTravel := mean + fixef(model.Travel)["(Intercept)"]]

# order by transitions
randomSims = randomSims[order(mean),]
randomSims[,HenID := factor(HenID, levels = as.character(HenID))]

#plot
ggplot(data = randomSims, aes(x = as.factor(HenID), y = mean))+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = RatioSplit), size = 2)+
  geom_point(aes(fill = ID), shape = 21, size = 5)+
  theme_classic(base_size = 18)+
  scale_fill_manual(values = c("red","blue","orange","yellow", "lightblue","white", "gray"))+
  scale_color_manual(values = c("grey", "black"))

###### Wintergarten ####

#start model
model.Garten = glmer(Out ~ Ratio + DateID + (1|Pen/HenID), family = binomial, data = varOfInterest)
model.Garten = glmer(Out ~ Ratio + DateID + (1|HenID), family = binomial, data = varOfInterest)
resid.Garten = simulateResiduals(model.Garten, 1000)
plot(resid.Garten)
plotResiduals(resid.Garten, form = varOfInterest$Ratio)
plotResiduals(resid.Garten, form = varOfInterest$DateID)
summary(model.Garten)

#TODO: Something wrong here, values weird
varOfInterest[, PredictGarten := predict(model.Garten)]

# ggplot(varOfInterest, aes(x = Date, y = Out, color = RatioSplit)) +
#   geom_jitter(size=2) + 
#   geom_line(data = varOfInterest[, mean(PredictGarten), by = .(Date, RatioSplit)], 
#             aes(x = Date, y =V1, colour = RatioSplit),size=1.5) + 
#   theme_classic(base_size = 18)+ 
#   ylab("Out in wintergarten") 
# 
# #plot individual variation
# ggplot(data = varOfInterest, aes(x = Date, y = PredictTravel, colour = as.factor(HenID))) + 
#   geom_line(aes(group = HenID), size=1)+
#   labs(y = "Predicted vertical travel distance",color = "Hen ID")+ 
#   theme_classic(base_size = 18)+ theme(legend.position="none")+ 
#   guides(color = guide_legend(nrow = 4))
# 
# #plot against actual data of some examples
# ggplot(data = varOfInterest[Highlight != "Any",], 
#        aes(x = Date, colour = as.factor(HenID))) + 
#   geom_jitter(aes(y = vertTravelDist))+
#   geom_line(aes(y = PredictTravel, group = HenID), size=1)+
#   labs(y = "Predicted vertical travel distance",color = "Hen ID")+ 
#   theme_classic(base_size = 18)+ 
#   guides(color = guide_legend(nrow = 4))

###### Nestbox Time #####
model.Nest = lmer(as.numeric(DurationNest) ~ Ratio + DateID + (1|Pen/HenID), data = varOfInterest)
model.Nest = lmer(as.numeric(TimeNest) ~ Ratio + DateID + (1|HenID),  data = varOfInterest)
resid.Nest = simulateResiduals(model.Nest, 1000)
plot(resid.Nest)
plotResiduals(resid.Nest, form = varOfInterest$Ratio[!is.na(varOfInterest$TimeNest)])
plotResiduals(resid.Nest, form = varOfInterest$DateID[!is.na(varOfInterest$TimeNest)])
summary(model.Nest)
