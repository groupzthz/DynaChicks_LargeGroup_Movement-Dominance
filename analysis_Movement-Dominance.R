#clear workspace
rm(list = ls())

#load packages
library(data.table)
library('ggplot2')
library('DescTools') #for CCC
library(lubridate)
library(lme4)
library(DHARMa)
library(effects)
source("functions.R")
library(hms)

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


##### Social Index #####

socialData = observData[!(Observer == 'Masha' & Reliability == 1),]

dataIndex = observData[, .(Dom=sum(Sum_Dom), Sub=sum(Sum_Sub)),  by = ID]

dataIndex[, Ratio := Dom/(Dom+Sub)]

dataIndex[, Pen := factor(c(rep(3,6), rep(4,6), rep(5,6), rep(10,6), rep(11,6), rep(12,6)))]

dataIndex = dataIndex[order(Ratio),]

socialData = combData[dataIndex, on = "ID"]

ggplot(data = socialData, aes(x = ID, y = Ratio, colour = as.factor(Pen)))+ 
  geom_point(size = 6)+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  theme_classic(base_size = 18)

##### Extraction tracking data #########

#add numerical hen identifier to tracking data
trackingData[, HenID := unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen)))]
trackingData[, HenID := as.numeric(HenID)]
trackingData[, Time := ymd_hms(Time)]

#relevant hens
hens = sort(unique(socialData$HenID))
#relevant times
times = list(ymd_hms(c("2019-11-12 03:30:00", "2019-11-12 17:30:00")),
             ymd_hms(c("2019-11-21 03:30:00", "2019-11-21 17:30:00")),
             ymd_hms(c("2019-11-30 02:30:00", "2019-11-30 17:30:00")),
             ymd_hms(c("2019-12-08 02:30:00", "2019-12-08 17:30:00")),
             ymd_hms(c("2019-12-20 02:30:00", "2019-12-20 17:30:00")))

#splitting data into Hens
splitHen = splitHenData(trackingData)


#TODO: change this to 23:59:59
hen_list <- vector(mode='list', length=length(hens))
start = ymd_hms("2019-11-11 02:00:00")
end = ymd_hms("2019-12-20 18:00:00")
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
#-> helath assessment on 09.12. -> exclude
trackingRel = trackingRel[!(Date == as.Date("2019-12-04")|
                              Date == as.Date("2019-12-06")|
                            Date == as.Date("2019-12-09")),]
length(unique(trackingRel$Date)) #-> now we have 37 days in total of data


#all individuals have data?
length(unique(trackingRel$HenID)) == 36
#all individuals have data every day?
entriesPerDay = trackingRel[, length(unique(Date)), by = HenID]
entriesPerDay[V1 != 37,] #Hen 1 (12), 65(33) and 84(35) don't have every day
#which days are missing?
days = unique(trackingRel$Date)

trackingRel[HenID == 1,unique(Date)] #no data between the 22.11. and 19.12.
trackingRel[HenID == 65,unique(Date)] #no data between 2.12. and 08.12.
trackingRel[HenID == 84,unique(Date)] #no data between 07.12. and 08.12.


#exclude HenID 1 and exclude december days before battery replacement
trackingRel = trackingRel[!(Date == as.Date("2019-12-02")|
                              Date == as.Date("2019-12-03")|
                              Date == as.Date("2019-12-07")|
                            Date == as.Date("2019-12-08")|
                              HenID == 1),]
entriesPerDay = trackingRel[, length(unique(Date)), by = HenID]


########### plot examples #########################
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



########## Parameters ################


#mark true transitions
trackingRel[, TrueTrans := TRUE]
#mark the first stamps as false transitions
trackingRel[1:length(unique(HenID)), TrueTrans := FALSE]

#add duration for all entries -> careful need to first insert begin and end of day stamps if daily durations want to be calculated!

#add end of each day
trackingFull = rbind(trackingRel, 
                     trackingRel[, tail(.SD, 1) , by = .(HenID, Date)][, c("Time", "TrueTrans") := .(ymd_hms(paste(Date, "23:59:59")), FALSE)])[order(HenID,Date)]

#add beginning of each day
trackingFull = rbind(trackingRel[, head(.SD, 1) , by = .(HenID, Date)][, c("Time", "TrueTrans") := .(ymd_hms(paste(Date, "00:00:00")), FALSE)], 
                     trackingFull)[order(HenID,Date)]

#include light-dark cycle
trackingFull[, Light := TRUE]
trackingFull[hour(Time) > 17, Light := FALSE]
trackingFull[month(Time)== 11 & day(Time) < 15 & hour(Time) < 4, Light := FALSE]
trackingFull[month(Time)== 11 & day(Time) < 22 & day(Time) > 14 & hour(Time) < 3, Light := FALSE]
trackingFull[(month(Time)== 12| (month(Time)== 11 & day(Time) > 21)) & hour(Time) < 2, Light := FALSE]

#add beginning and end of the night


#add duration
trackingFull[, Duration := (shift(Time, type="lead") - Time), by = HenID]
# Set duration for last day entry as 1 sec
trackingFull[hour(Time)== 23 & minute(Time)== 59, Duration := 1]

#transitions per bird
transitions = trackingFull[TrueTrans == TRUE, .(Transitions = .N), by = .(HenID)]
#transitions per bird per day
transDaily = trackingFull[TrueTrans == TRUE, .(Transitions = .N), by = .(HenID, Date)]

#durations per zone per bird
durations = trackingFull[, .(Duration = sum(Duration)), by = .(HenID, Zone)]
#durations per zone per bird per day
durDaily = trackingFull[, .(Duration = sum(Duration)), by = .(HenID, Zone, Date)]

plotDataDur = durations[socialData, on = "HenID"]

plotDataTrans = transitions[socialData, on = "HenID"]

ggplot(plotDataDur, aes(x = Ratio, y = Duration, colour = Zone))+
  geom_point()+
  geom_smooth()

ggplot(plotDataTrans, aes(x = Ratio, y = Transitions))+
  geom_point()+
  geom_smooth()

plotDataTransDaily = transDaily[socialData, on = "HenID"]

plotDataDurDaily = durDaily[socialData, on = "HenID"]

dataTransVar = plotDataTrans[, .(mean = mean(Transitions), sd = sd(Transitions), Ratio = mean(Ratio)), by = HenID]
dataDurVar = plotDataDur[, .(mean = mean(Duration), sd = sd(Duration), Ratio = mean(Ratio)), by = .(HenID, Zone)]

ggplot(plotDataTransDaily, aes(x = Ratio, y = Transitions))+
  geom_point(aes(colour = as.factor(HenID)))+
  geom_point(data = plotDataTransDaily[order(Ratio), .(mean = mean(Transitions)), by = Ratio], aes(x = Ratio, y = mean), 
             colour = "red", size = 3)

ggplot(plotDataTransDaily, aes(x = as.factor(Date), y = Transitions))+
  geom_line(aes(group = HenID), colour = "darkgrey")+
  geom_point(aes(colour = Ratio), size = 3)


ggplot(na.omit(plotDataDurDaily), aes(x = as.factor(Date), y = Duration))+
  geom_line(aes(group = HenID), colour = "darkgrey")+
  geom_point(aes(colour = Ratio))+
  facet_grid(.~Zone)

ggplot(na.omit(plotDataDurDaily[Zone == "Ramp_Nestbox",]), aes(x = as.factor(Date), y = Duration))+
  geom_line(aes(group = HenID), colour = "darkgrey")+
  geom_point(aes(colour = Ratio))

########## Test plot ###############

plotData = trackingFull
plotData[, Date := as.factor(Date)]
plotData[, Time_x := as_hms(Time)]
plotData[, Zone := factor(Zone, levels= c("Wintergarten", "Litter", "Tier_2", "Ramp_Nestbox", "Tier_4"))]

#ideas for plots: for each individual an average step function plot
#requires: average movement data: 

############ Comb Size ##############


ggplot(data = socialData, aes(x = Mean5_7, y = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_hline(yintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'Comb Size', y = "Dominance rank")


comb.model = lmer(Ratio ~ Mean5_7 + (1|Pen), data = socialData)
resid.comb = simulateResiduals(comb.model)
plot(resid.comb)
plotResiduals(resid.comb, form = socialData$Mean5_7)

summary(comb.model)
comb.model.null = lmer(Ratio ~ 1 + (1|Pen), data = socialData)
anova(comb.model, comb.model.null, type = "chisq")
