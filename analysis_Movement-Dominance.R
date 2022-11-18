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
healthData[, footproblems:= r_podo + r_bumble +r_injure + l_podo + l_bumble + l_injure]
healthData[, WoA := ifelse(date == "09.12.2019", 26, 55)]
healthDataWide = dcast(healthData, formula = ID ~ WoA, value.var = list("weight", "feathers", "wounds","comb", "bare_head", "footproblems"))

KBF = fread("KBF_scores.csv")
KBF = KBF[(Date == "09.12.2019" | Date == "29.06.2020" | Date == "30.06.2020")
                        & HenID %in% unique(observData$ID),]
KBF[, WoA := ifelse(Date == "09.12.2019", 26, 55)]
KBFWide = dcast(KBF, formula = HenID ~ WoA, value.var = "Severity")
colnames(KBFWide)[colnames(KBFWide) == "26"] = "KBF_26"
colnames(KBFWide)[colnames(KBFWide) == "55"] = "KBF_55"

healthData = KBF[, .(WoA, HenID, Severity)][healthData[, .(WoA, ID, WoA, weight, feathers, wounds, comb, bare_head, footproblems)], on = c(HenID = "ID",WoA = "WoA")]
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


###### Social Data prep & Index #####
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

henData[, Pen := as.factor(Pen)]

rm(dataIndex)
###### Badges of Status #######

#comb size by dominance
ggplot(data = henData, aes(x = Comb, y = Ratio))+ 
  # geom_violin()+
  geom_point(aes(colour = as.factor(Pen)),size = 6)+
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

#weight at 55 by comb size
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
null.Comb =  lmer(Ratio ~ 1 +(1|Pen), data = henData)
resid.Comb = simulateResiduals(model.Comb, 1000)
plot(resid.Comb)
plotResiduals(resid.Comb, form = henData$Comb)
summary(model.Comb)
anova(model.Comb, null.Comb)
parameters(model.Comb)

henData[, PredictRatio:= predict(model.Comb)]
ggplot(data = henData, aes(y = Ratio, x = Comb, ))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  geom_line(aes(y = PredictRatio, colour = as.factor(Pen)))+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

###### Weight gain by dominance ####

#TODO: check out gain ab 20WoA


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
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  geom_line(aes(y = PredictGain, colour = as.factor(Pen)), size = 1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)


###### Health & social data ####

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
  labs(x = 'WoA', y = "KBF Severity")
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
ggplot(data = henData, aes(y = feathers_55, x = Ratio, colour = Pen))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, se = F)+
  theme_classic(base_size = 18)+
  labs(y = 'Feather cover', x = "Dominance ratio")

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


#foot problems
ggplot(data = henData, aes(y = footproblems_55, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 6)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'Health & Dominance', y = 'injuries', x = "Dominance ratio")

#injuries
ggplot(data = henData, aes(y = wounds_55, x = Ratio))+ 
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


#health pca
health.pca <- prcomp(henData[,.(feathers_55, footproblems_55, 
                                KBF_55, gain)], scale = TRUE, center = T)
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
###### Extraction tracking data #########

#add numerical hen identifier to tracking data
trackingData[, HenID := unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen)))]
trackingData[, HenID := as.numeric(HenID)]
trackingData[, Time := ymd_hms(Time)]

#relevant hens
hens = sort(unique(socialData$HenID))

#splitting data into Hens
splitHen = splitHenData(trackingData)

hen_list <- vector(mode='list', length=length(hens))
start = ymd_hms("2019-11-09 00:00:00")
end = ymd_hms("2020-06-28 23:59:59")
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

###### data checks #############################

#there is data for every day?
allDays = seq(ymd("2019-11-09"), to = ymd("2020-06-28"), by = "day")
length(unique(trackingRel$Date)) == length(allDays)
#which days are missing?
allDays[which(!(allDays %in% unique(trackingRel$Date)))]
#missing days:
#"2019-12-05" "2019-12-22" "2019-12-27" "2019-12-28" "2019-12-29" "2019-12-30" 
#"2019-12-31" "2020-01-01" "2020-01-02" "2020-01-03" "2020-01-04" "2020-05-19" "2020-05-24"

#days to exclude
excl = c(ymd("2019-12-09"), ymd("2020-01-06"), ymd("2020-02-10"),
         ymd("2020-03-16"), ymd("2020-05-04"), ymd("2020-06-02"), #health assessments
         ymd("2019-12-04"), ymd("2019-12-06"), ymd("2020-01-12"), ymd("2020-06-20"), #missing data
         seq(ymd("2020-01-21"), to = ymd("2020-02-03"), by = "day"), #different configuration
         seq(ymd("2020-05-16"), to = ymd("2020-05-26"), by = "day"), # experimental light change
         ymd("2020-01-07"),ymd("2020-02-25"),ymd("2020-03-11") #other disturbances
         )

trackingRel = trackingRel[!(Date %in% excl),]

length(unique(trackingRel$Date)) #-> now we have 184 days in total of data

#do all individuals have data?
length(unique(trackingRel$HenID)) == 36
#all individuals have data every day?
entriesPerDay = trackingRel[, length(unique(Date)), by = HenID][V1 != 184,][, miss := 184-V1]
entriesPerDay

#which days are missing for those hens?
days = unique(trackingRel$Date)

trackingRel[HenID == 1, days[which(!(days %in% Date))]] #42 days miss, no data between the 22.11. - 26.12. + 05.01. + 09.- 20.01.
trackingRel[HenID == 17, days[which(!(days %in% Date))]] # 2 days miss, no data on the 23.12. & 05.01.
trackingRel[HenID == 33, days[which(!(days %in% Date))]] # 1 day miss, no data on the 17.03.
trackingRel[HenID == 35, days[which(!(days %in% Date))]] # 1 day miss, no data on the 23.12.
trackingRel[HenID == 39, days[which(!(days %in% Date))]] # 1 day miss, no data on the 05.01.
trackingRel[HenID == 65,days[which(!(days %in% Date))]] #8 day miss, no data between 2.12. - 08.12. + 23.12. - 26.12.
trackingRel[HenID == 84,days[which(!(days %in% Date))]] #2 days miss, no data between 07.12. and 08.12.


#exclude identified problem days
trackingRel = trackingRel[!(Date == as.Date("2020-01-05")|
                              Date == as.Date("2019-12-23")),]
entriesPerDay = trackingRel[, length(unique(Date)), by = HenID][V1 != 182,][, miss := 182-V1]
length(unique(trackingRel$Date)) #-> now we have 172 days in total of data
#TODO: what to do with hens with missing days? (especially HenID 1 (40 miss) & 65 (7 miss))

###### time series plot examples #########################
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

# #ideas for plots: for each individual an average step function plot
# #requires: full sequence for every day
# 
# splitHen = splitHenData(trackingFull)
# hens = unique(trackingFull$HenID)
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

###### Data preparations ################

#data cleaning step
#remove Wintergarten entries outside time that Wintergarten can be open
#delete all Wintergarten entries at night (not possible, tracking system errors)
trackingRel = trackingRel[!(hour(Time) > 16 & minute(Time) > 45 & Zone == "Wintergarten"), ]
#delete all Wintergarten entries before 9 in the morning (not possible, tracking system error)
trackingRel = trackingRel[!(hour(Time) < 10 & Zone == "Wintergarten"), ]

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
rm(splitHen, trackingRel)

#make sure date is correct (now with inserting of end and start of day)
trackingFull[, Date:= as_date(Time)]

#include Week of age
tableWoA = data.table(Date = seq(ymd("2019-10-11"), ymd("2020-07-02"), by = "day"),
                      WoA = rep(18:55, each = 7))
trackingFull = trackingFull[tableWoA, on = "Date", nomatch = NULL]

# mark true transitions
trackingFull[, TrueTransition := TRUE]
trackingFull[, DayIndic := F]
trackingFull[hour(Time) == 0 & minute(Time) == 0 & second(Time) == 0, TrueTransition := FALSE] 
trackingFull[hour(Time) == 23 & minute(Time) == 59 & second(Time) == 59, TrueTransition := FALSE] 
trackingFull[hour(Time) == 0 & minute(Time) == 0 & second(Time) == 0, DayIndic := TRUE] 
trackingFull[hour(Time) == 23 & minute(Time) == 59 & second(Time) == 59, DayIndic := TRUE] 


#include light-dark cycle
#TODO: how did it work on the 29.03.2020?? (Zeitumstellung)
trackingFull[, Light := TRUE]
trackingFull[hour(Time) > 17, Light := FALSE]
trackingFull[month(Time)== 11 & day(Time) < 15 & hour(Time) < 4, Light := FALSE]
trackingFull[month(Time)== 11 & day(Time) < 22 & day(Time) > 14 & hour(Time) < 3, Light := FALSE]
trackingFull[(month(Time)!= 11| (month(Time)== 11 & day(Time) > 21)) & hour(Time) < 2, Light := FALSE]

#catching Wintergarten errors of inserting time stamps
trackingFull[Light == F & Zone == "Wintergarten", Zone := "Litter"]

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
    if (nrow(relDat)){
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
    }
    full_list[[i]][[j]] = relDat
    j =j+1
  }
  full_list[[i]] = mergeHenData(full_list[[i]])
}
trackingFull = mergeHenData(full_list)

rm(hen_list, full_list)

#add duration
trackingFull[, Duration := (shift(Time, type="lead") - Time), by = HenID]
#set all last day entries to 1 sec so that missing days don't count into the duration
trackingFull[hour(Time)== 23 & minute(Time)== 59, Duration := 1]

#indicate which night times belong together
trackingFull[Light == T, NightCycle := Date]
trackingFull[hour(Time) > 17|hour(Time) == 17 , NightCycle := Date]
trackingFull[hour(Time) < 4 & Light == F, NightCycle := Date-1]
#delete dupliacted entries
trackingFull = trackingFull[Duration > 0,]

##### Sequence similarity ####

#sequence similarity calculated by how many seconds two timelines agree and how many they disagree
#calculated within birds -> comparing consecutive (or further away?) days
#Question: split by day and night hours?
#Problem: need to deal with time shift? I think not will resolve themselves  

###### Between individuals ####
betweenIndividuals = similarityBetween(trackingFull)


mean(as.matrix(betweenIndividuals[, -c(1,2)]), na.rm = T)

ggplot(data = betweenIndividuals[, c(1,2,4)], aes(x = V1, y = V2))+
  geom_tile(aes(fill = as.numeric(`2019-11-12`)), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")

betweenIndividuals[, pair := paste(V1, "-", V2)]
betweenIndividualsL = melt(betweenIndividuals,                          # Reshape data from wide to long format
                              id.vars     = c("pair", "V1", "V2"),
                           variable.name = "Date", 
                           value.name = "Similarity")
betweenIndividualsL[, Date := as_date(Date)]


ggplot(data = betweenIndividualsL[V1 == "Hen_108",], aes(x = Date, y = Similarity))+
  geom_line(aes(group = pair))

agrBetween = betweenIndividualsL[, .(Mean = mean(Similarity, na.rm = T), 
                                     SD = sd(Similarity, na.rm = T)), by = .(pair, V1, V2)]
ggplot(data = agrBetween, aes(x = V1, y = V2))+
  geom_tile(aes(fill = Mean), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")


#Assortativity coefficients close to 1 indicate that there is very high 
#likelihood of two vertices with the same property being connected.
assort = data.table(Date = as_date(colnames(betweenIndividuals)[-c(1,2, .N)]))
assort[, Assort := 0]
# high degree assortativity is a measure of preferential attachment in organizations, 
# where highly connected vertices are connected with each other and a large number 
# of vertices with low degree make up the remainder of the network.
assort[, Degree := 0]
dailySimil = vector(mode='list', length= length(colnames(betweenIndividuals)[-c(1,2, .N)]))
i = 1
for(day in colnames(betweenIndividuals)[-c(1,2, .N)]){
  cat(day)
  dailySimil[[i]] = as.matrix(rbind(cbind(data.table(Hen_108 = rep(NA, 34)), 
                                  dcast(betweenIndividuals[,.SD, .SDcols = c("V1", "V2", (day))], V1 ~ V2)[,-1]),
                            data.table(Hen_108 = NA), fill = T))
  rownames(dailySimil[[i]]) = colnames(dailySimil[[i]])
  simMatrix = dailySimil[[i]]
  simMatrix = simMatrix *10
  #maybe set threshold smarter by checking 
  simMatrix[simMatrix < 3.5] = 0
  g  <- graph.adjacency(simMatrix,mode = "undirected", weighted = T, diag = F)
  #add network node attribute
  V(g)$domIndex <- socialData$Ratio[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
  assort[Date == (day), Assort := assortativity(g, V(g)$domIndex)]
  assort[Date == (day), Degree := assortativity_degree(g)]
  i = i+1
}

similSimply = lapply(dailySimil, function(x){ x[x < 0.6] = 0; x[(x > 0.6)|(x == 0.6)] = 1; x})

similSimply = Reduce('+', similSimply)
g  <- graph.adjacency(similSimply,mode = "undirected", diag = F)
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

### Isolate by Nestbox behaviour

btwIndivNest = similarityBetween(trackingFull, zone = "Ramp_Nestbox", interval = "morning")

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

withinIndividuals = similarityWithin(trackingFull)
withinIndividualsL = melt(withinIndividuals,                          # Reshape data from wide to long format
                           id.vars     = c("Date"),
                           variable.name = "Hen", 
                           value.name = "Similarity")
withinIndividualsL[, HenID := as.numeric(unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen))))]
withinIndividualsL = withinIndividualsL[socialData[,c(1,4,13)], on = "HenID"]


ggplot(data = na.omit(withinIndividualsL), aes(x = Date, y = Similarity))+
  geom_line(aes(group = Hen, colour = Ratio))

##### Parameters ########

# calculate daily parameters per bird 

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

#durations per zone per bird per day during dark hours
durDailyD = trackingFull[Light == F, .(DurationNight = sum(Duration)), by = .(HenID, Zone, NightCycle)]
durDailyD = durDailyD[order(HenID),]
#extract maximum Zone for each bird per day
mainDurDailyD = durDailyD[, .SD[which.max(DurationNight)], by = .(HenID, NightCycle)]
mainDurDailyD[, onTop := ifelse(Zone == "Tier_4", 1, 0)]


###### wintergarden use #####
#all Wintergarten entries per bird
#careful: on vaccination days garten opened later!
dailyGarten = trackingFull[Zone == "Wintergarten", .(Duration = Duration), by = .(HenID, Time, Date)]
#extract if hen goes out on day or not
inGarten = trackingFull[, .(Out = ifelse(any(Zone == "Wintergarten"), 1, 0)), by = .(HenID, Date)]
#extract how long each hen went out per day
durDailyGarten = dailyGarten[, .(DurationGarten = sum(Duration)), by = .(HenID, Date)]
# latency to go out
latGarten = dailyGarten[, .(LatencyGarten = Time[1] - ymd_hms(paste(as_date(Time[1]), "10:00:00"))), by = .(HenID, Date)]

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
#TODO: 
#Feeder runs: (ab 22.11.: 2:00), 4:00, 6:00, 8:00, 10:00, 13:00, 15:00, 16:15


#create data.table containing all daily measures
varOfInterest = trackingFull[Light == T, .(vertTravelDist = sum(distVertical)), by = .(HenID, Date, Pen, WoA)]
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
varOfInterest = socialData[, .(HenID,Ratio,Comb)][varOfInterest, on = "HenID"] 


##### Parameter inspection plots ####

#Highlight three most extreme
varOfInterest[, Highlight := "Any"]
varOfInterest[HenID == 82 | HenID == 97 | HenID == 33, Highlight := "Dom"]
varOfInterest[HenID == 77 | HenID == 5 | HenID == 108, Highlight := "Sub"]
#make factor out of HenID with levels sorted in ascending Ratio order
varOfInterest[,HenID := factor(HenID, levels = socialData$HenID)]


###### sleeping spot ####

#TODO: try plot with comb size
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
ggplot(varOfInterest, aes(x = Date, y = vertTravelDist, colour = Ratio))+
  geom_point()+
  geom_smooth(aes(group = as.factor(HenID)),se = F)
  #scale_colour_manual(values = c("grey", "red", "blue"))

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

###### Vertical distance ####

model.Travel = lmer(vertTravelDist ~ Ratio*WoA + (1|Pen/HenID), data = varOfInterest)
resid.Travel = simulateResiduals(model.Travel, 1000)
plot(resid.Travel)
plotResiduals(resid.Travel, form = varOfInterest$Ratio)
plotResiduals(resid.Travel, form = varOfInterest$WoA)
summary(model.Travel)
plot(allEffects(model.Travel))

varOfInterest[, PredictTravel := predict(model.Travel)]

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

#divide the variance explained by animal_id by the total phenotypic variance 
#(animal_id + month:year + year + residual variance)
print(VarCorr(model.Travel), comp = c("Variance", "Std.Dev."))
VarCorr(model.Travel)$"HenID:Pen"[1] / (VarCorr(model.Travel)$"HenID:Pen"[1] + 
                                VarCorr(model.Travel)$"Pen"[1] + 
                                attr(VarCorr(model.Travel), "sc")^2)
#-> within -individual repeatability of 0.35 (within a pen)

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
quantile(CVi,prob=c(0.025, 0.5, 0.975))

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
randomSims[, meanTravel := mean + fixef(model.Travel)["(Intercept)"]]

# order by transitions
randomSims = randomSims[order(meanTravel),]
randomSims[,HenID := factor(HenID, levels = as.character(HenID))]

#plot
ggplot(data = randomSims, aes(x = as.factor(HenID), y = meanTravel))+
  geom_errorbar(aes(ymin = meanTravel-sd, ymax = meanTravel+sd, color = RatioSplit), size = 2)+
  geom_point(aes(fill = ID), shape = 21, size = 5)+
  theme_classic(base_size = 18)+
  scale_fill_manual(values = c("red","blue","orange","yellow", "lightblue","white", "gray"))+
  scale_color_manual(values = c("grey", "black"))

###### Wintergarten ####

#start model
model.Garten = glmer(Out ~ Ratio*WoA + (1|Pen/HenID), family = binomial, data = varOfInterest)
resid.Garten = simulateResiduals(model.Garten, 1000)
plot(resid.Garten)
plotResiduals(resid.Garten, form = varOfInterest$Ratio)
plotResiduals(resid.Garten, form = varOfInterest$WoA)
summary(model.Garten)
plot(allEffects(model.Garten))

varOfInterest[, PredictGarten := predict(model.Garten, type="response")]

ggplot(varOfInterest, aes(x = WoA, y = Out, color = RatioSplit)) +
  #geom_jitter(size=2) +
  geom_line(data = varOfInterest[, mean(PredictGarten), by = .(WoA, RatioSplit)],
            aes(x = WoA, y =V1, colour = RatioSplit),size=1.5) +
  theme_classic(base_size = 18)+
  ylab("Out in wintergarten")

#plot individual variation
ggplot(data = varOfInterest, aes(x = WoA, y = PredictGarten, colour = as.factor(HenID))) +
  geom_line(aes(group = HenID), size=1)+
  labs(y = "Predicted probability to go in wintergarten",color = "Hen ID")+
  theme_classic(base_size = 18)+ theme(legend.position="none")+
  guides(color = guide_legend(nrow = 4))

#plot against actual data of some examples
ggplot(data = varOfInterest[Highlight != "Any",],
       aes(x = WoA, colour = as.factor(HenID))) +
  geom_jitter(aes(y = Out))+
  geom_line(aes(y = PredictGarten, group = HenID), size=1)+
  labs(y = "Predicted vertical travel distance",color = "Hen ID")+
  theme_classic(base_size = 18)+
  guides(color = guide_legend(nrow = 4))

###### Nestbox Time #####
model.Nest = lmer(as.numeric(TimeNest) ~ Ratio*WoA+ (1|Pen/HenID), data = varOfInterest)
model.Nest = lmer(as.numeric(TimeNest) ~ Ratio*WoA + (1|HenID),  data = varOfInterest)
resid.Nest = simulateResiduals(model.Nest, 1000)
plot(resid.Nest)
plotResiduals(resid.Nest, form = varOfInterest$Ratio[!is.na(varOfInterest$TimeNest)])
plotResiduals(resid.Nest, form = varOfInterest$DateID[!is.na(varOfInterest$TimeNest)])
summary(model.Nest)
