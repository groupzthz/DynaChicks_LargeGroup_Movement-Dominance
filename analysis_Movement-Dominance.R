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
socialData[, Hen := paste0("Hen_", HenID)]

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

henData[, PredictRatio:= predict(model.Comb)]
ggplot(data = henData, aes(y = Ratio, x = Comb, ))+ 
  geom_point(aes(colour = as.factor(Pen)), size = 3)+
  geom_line(aes(y = PredictRatio, colour = as.factor(Pen)))+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

###### Weight gain by dominance ####


#weight gain
ggplot(data = henDataLong, aes(x = as.factor(WoA), y = weight))+ 
  # geom_violin()+
  geom_point(size = 6, aes(colour = Ratio))+
  geom_line(aes(group = ID))+
  theme_classic(base_size = 18)+
  labs(title = 'Visual signalling', x = 'WoA', y = "Weight")

henData[, gain := weight_55-weight_20]

#gain by dominance
ggplot(data = henData, aes(y = gain, x = Ratio))+ 
  # geom_violin()+
  geom_point(size = 3)+
  geom_vline(xintercept= 0.5, linetype='dashed')+
  geom_smooth(method = 'lm', size = 1.5, colour= 'red')+
  theme_classic(base_size = 18)+
  labs(title = 'weight gain', y = 'Weight gain', x = "Dominance ratio")


# model.Gain = lmer(gain ~ Ratio +(1|Pen), data = henData)
# null.Gain =  lmer(gain ~ 1 +(1|Pen), data = henData)
# resid.Gain = simulateResiduals(model.Gain, 1000)
# plot(resid.Gain)
# plotResiduals(resid.Gain, form = henData$Ratio)
# summary(model.Gain)
# anova(model.Gain, null.Gain)
# parameters(model.Gain)
# 
# henData[, PredictGain:= predict(model.Gain)]
# ggplot(data = henData, aes(x = Ratio, y = gain, ))+ 
#   geom_point(aes(colour = as.factor(Pen)), size = 3)+
#   geom_line(aes(y = PredictGain, colour = as.factor(Pen)), size = 1)+
#   #facet_grid(.~Pen)+
#   theme_classic(base_size = 18)

#TODO: take out second order interaction
model.Weight = lmer(weight ~ poly(WoA, 2)+ WoA:Ratio + Ratio +(1|Pen), data = henDataLong)
null.Weight =  lmer(weight ~ 1 +(1|Pen), data = henDataLong)
resid.Weight = simulateResiduals(model.Weight, 1000)
plot(resid.Weight)
plotResiduals(resid.Weight, form = henDataLong$Ratio)
plotResiduals(resid.Weight, form = henDataLong$WoA)
summary(model.Weight)
anova(model.Weight, null.Weight)
parameters(model.Weight)
plot(allEffects(model.Weight))

library(emmeans)
plotData = as.data.table(emmeans(model.Weight, ~ pairwise ~ WoA*Ratio, at =  list(Ratio = quantile(henDataLong$Ratio),
                                                         WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))$emmeans)

henDataLong[, PredictWeight:= predict(model.Weight)]
henDataLong[, RatioSplit := ifelse(Ratio >0.5,"Dom" , "Sub")]

ggplot()+ 
  geom_point(data = henDataLong, aes(x = WoA, y = weight), size = 3)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = lower.CL, ymax = upper.CL, group = Ratio), alpha = 0.1)+
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

model.KBF = lmer(Severity ~ Ratio*WoA +(1|Pen), data = henDataLong)
null.KBF =  lmer(Severity ~ 1 +(1|Pen), data = henDataLong)
resid.KBF = simulateResiduals(model.KBF, 1000)
plot(resid.KBF)
testZeroInflation(resid.KBF) 
#zero-inflation problem:
#tested by running model without 20 WoA and 26 WoA -> estimates do not change   
plotResiduals(resid.KBF, form = henDataLong$Ratio)
plotResiduals(resid.KBF, form = henDataLong$WoA)
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

model.Feathers = glmer.nb(feathers ~ Ratio + Ratio:WoA+ poly(WoA,2) +(1|Pen), data = henDataLong) 
null.Feathers =  glmer.nb(feathers ~ 1 +(1|Pen), data = henDataLong)
resid.Feathers = simulateResiduals(model.Feathers, 1000)
plot(resid.Feathers)
plotResiduals(resid.Feathers, form = henDataLong$Ratio)
plotResiduals(resid.Feathers, form = henDataLong$WoA)
summary(model.Feathers)
anova(model.Feathers, null.Feathers)
parameters(model.Feathers) 
plot(allEffects(model.Feathers))


#foot problems
ggplot(data = henDataLong, aes(x = as.factor(WoA), y = footproblems))+ 
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

#relevant hens
hens = sort(unique(socialData$HenID))

trackingData = prepareTrackingData(trackingData, hens)
#include Week of age
trackingData = trackingData[tableWoA, on = "Date", nomatch = NULL]

##### Sequence similarity ####

#sequence similarity calculated by how many seconds two timelines agree and how many they disagree
#calculated within birds -> comparing consecutive (or further away?) days
#Question: split by day and night hours?
#Problem: need to deal with time shift? I think not will resolve themselves  

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
                                unlist(lapply(betweenIndividuals2, function(x){ sum(x<0.2, na.rm = TRUE)})))
)

ggplot(meanData, aes(x = Date))+
  geom_pointrange(aes(y = Mean, ymin = Mean-SD, ymax = Mean+SD))+
  geom_point(aes(y = Median), color = "red")+
  #geom_point(aes(y = Max))+
  theme_classic(base_size = 18)
ggplot(meanData, aes(x = Date))+
  #geom_point(aes(y = Below))+
  #geom_smooth(aes(y = Below))+
  geom_point(aes(y = Above))+
  geom_smooth(aes(y = Above))+
  theme_classic(base_size = 18)

similSimply1 = lapply(betweenIndividuals, function(x){ x[x < 0.7|is.na(x)] = 0; x[(x > 0.7)|(x == 0.7)] = 1; x})
similSimply1 = Reduce('+', similSimply1)
similSimply2 = lapply(betweenIndividuals2, function(x){ x[x < 0.7|is.na(x)] = 0; x[(x > 0.7)|(x == 0.7)] = 1; x})
similSimply2 = Reduce('+', similSimply2)
similSimply = similSimply1 + similSimply2

similSimply1 = lapply(betweenIndividuals, function(x){ x[is.na(x)] = 0; x})
similSimply1 = Reduce('+', similSimply1)
similSimply2 = lapply(betweenIndividuals2, function(x){ x[is.na(x)] = 0; x})
similSimply2 = Reduce('+', similSimply2)
similSimply = similSimply1 + similSimply2
noData1 = lapply(betweenIndividuals, function(x){ x[!is.na(x)] = 0; x[is.na(x)] = 1; x})
noData1 = Reduce('+', noData1)
noData2 = lapply(betweenIndividuals2, function(x){ x[!is.na(x)] = 0; x[is.na(x)] = 1; x})
noData2 = Reduce('+', noData2)
noData = (noData1 + noData2 - 179)*-1

similSimply = similSimply/noData

mean(similSimply, na.rm = TRUE)
max(similSimply, na.rm = TRUE)

similSimply[is.na(similSimply)] = 0
similSimply = similSimply + t(similSimply)


#[similSimply <0.45]= 0

g  <- graph.adjacency(similSimply,mode = "undirected", weighted = TRUE, diag = F)
bc <- edge.betweenness.community(g)
#par(mfrow=c(1,2))
plot(as.dendrogram(bc))
#network vertex names
V(g)$name
#inspect network edge and node attributes
edge_attr(g)
vertex_attr(g)
V(g)$domIndex <- socialData$Ratio[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
V(g)$Pen <- socialData$Pen[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
plot(g, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     vertex.color=c( "pink", "skyblue")[1+(V(g)$domIndex>0.5)] ) 
assortativity(g, V(g)$domIndex, types1 = graph.strength(g))
E(g)$width <- E(g)$weight*3
plot(g, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     vertex.color=V(g)$Pen) 
assortativity(g, V(g)$Pen)
library(assortnet)
assortment.continuous(similSimply, V(g)$domIndex, weighted = TRUE, SE = FALSE, M = 1)

#old first tests
# mean(as.matrix(betweenIndividuals[, -c(1,2)]), na.rm = T)
# 
# ggplot(data = betweenIndividuals[, c(1,2,4)], aes(x = V1, y = V2))+
#   geom_tile(aes(fill = as.numeric(`2019-11-12`)), colour = "white") +
#   scale_fill_gradient(low = "white", high = "red")
# 
# betweenIndividuals[, pair := paste(V1, "-", V2)]
# betweenIndividualsL = melt(betweenIndividuals,                          # Reshape data from wide to long format
#                               id.vars     = c("pair", "V1", "V2"),
#                            variable.name = "Date", 
#                            value.name = "Similarity")
# betweenIndividualsL[, Date := as_date(Date)]
# 
# 
# ggplot(data = betweenIndividualsL[V1 == "Hen_108",], aes(x = Date, y = Similarity))+
#   geom_line(aes(group = pair))
# 
# agrBetween = betweenIndividualsL[, .(Mean = mean(Similarity, na.rm = T), 
#                                      SD = sd(Similarity, na.rm = T)), by = .(pair, V1, V2)]
# ggplot(data = agrBetween, aes(x = V1, y = V2))+
#   geom_tile(aes(fill = Mean), colour = "white") +
#   scale_fill_gradient(low = "white", high = "red")
# 
# 
# #Assortativity coefficients close to 1 indicate that there is very high 
# #likelihood of two vertices with the same property being connected.
# assort = data.table(Date = as_date(colnames(betweenIndividuals)[-c(1,2, .N)]))
# assort[, Assort := 0]
# # high degree assortativity is a measure of preferential attachment in organizations, 
# # where highly connected vertices are connected with each other and a large number 
# # of vertices with low degree make up the remainder of the network.
# assort[, Degree := 0]
# dailySimil = vector(mode='list', length= length(colnames(betweenIndividuals)[-c(1,2, .N)]))
# i = 1
# for(day in colnames(betweenIndividuals)[-c(1,2, .N)]){
#   cat(day)
#   dailySimil[[i]] = as.matrix(rbind(cbind(data.table(Hen_108 = rep(NA, 34)), 
#                                   dcast(betweenIndividuals[,.SD, .SDcols = c("V1", "V2", (day))], V1 ~ V2)[,-1]),
#                             data.table(Hen_108 = NA), fill = T))
#   rownames(dailySimil[[i]]) = colnames(dailySimil[[i]])
#   simMatrix = dailySimil[[i]]
#   simMatrix = simMatrix *10
#   #maybe set threshold smarter by checking 
#   simMatrix[simMatrix < 3.5] = 0
#   g  <- graph.adjacency(simMatrix,mode = "undirected", weighted = T, diag = F)
#   #add network node attribute
#   V(g)$domIndex <- socialData$Ratio[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
#   assort[Date == (day), Assort := assortativity(g, V(g)$domIndex)]
#   assort[Date == (day), Degree := assortativity_degree(g)]
#   i = i+1
# }
# 
# similSimply = lapply(dailySimil, function(x){ x[x < 0.6] = 0; x[(x > 0.6)|(x == 0.6)] = 1; x})
# 
# similSimply = Reduce('+', similSimply)
# g  <- graph.adjacency(similSimply,mode = "undirected", diag = F)
# bc <- edge.betweenness.community(g)
# #par(mfrow=c(1,2))
# plot(as.dendrogram(bc))
# #network vertex names
# V(g)$name
# #inspect network edge and node attributes
# edge_attr(g)
# vertex_attr(g)
# V(g)$domIndex <- socialData$Ratio[match(V(g)$name, paste0("Hen_",as.character(socialData$HenID)))]
# plot(g, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
#      vertex.color=c( "pink", "skyblue")[1+(V(g)$domIndex>0.5)] ) 
# assortativity(g, V(g)$domIndex)


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

#build Model
hist(dataBetween$Similarity)

model.Between = lmer(Similarity ~ DiffRatio*PenBool+ DiffRatio*WoA + (1|Hen1) + (1|Hen2), data = dataBetween)
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

dataWithin[, PredictWithin := predict(model.Within)]


#test random Ratio allocation
#TODO: or randomise similarity within day to take out day und ratio effect

betaTable = data.table(Run = "Orig", t(fixef(model.Between)))
colnames(betaTable)[3] = "DiffRatio_rand"
colnames(betaTable)[6] = "DiffRatio_rand:PenBoolTRUE"
colnames(betaTable)[7] = "DiffRatio_rand:WoA"

set.seed(42)
for (i in 1:1000){
  cat("Run:",i, "\n")
  socialData[, randRatio := sample(Ratio)]
  dataBetween[, Hen1_Ratio_rand := socialData$randRatio[match(Hen1, socialData$Hen)]]
  dataBetween[, Hen2_Ratio_rand := socialData$randRatio[match(Hen2, socialData$Hen)]]
  dataBetween[, DiffRatio_rand := abs(Hen1_Ratio_rand - Hen2_Ratio_rand)]
  model.Between.rand = lmer(Similarity ~ DiffRatio_rand*PenBool+ DiffRatio_rand*WoA + (1|Hen1) + (1|Hen2), data = dataBetween)
  entry = data.table(Run = as.character(i), t(fixef(model.Between.rand)))
  betaTable= rbind(betaTable, entry)
}

ggplot(data = betaTable, aes(x= `DiffRatio_rand:WoA`))+
  geom_histogram()+
  geom_vline(xintercept = betaTable[Run == "Orig", `DiffRatio_rand:WoA`], colour = "red")

p = (sum(betaTable[Run != "Orig", `DiffRatio_rand:WoA`] <= betaTable[Run == "Orig", `DiffRatio_rand:WoA`]) + 1)/dim(betaTable)[1]

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

#TODO: null-model für p-value

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

model.Within = lmer(Similarity ~ Ratio*WoA + (1|Pen/HenID), data = dataWithin)
#model.Within = lmer(Similarity ~ Ratio + (WoA|HenID), data = dataWithin)
resid.Within = simulateResiduals(model.Within, 1000)
plot(resid.Within) #TODO: okay? -> bootstrapping?
plotResiduals(resid.Within, form = dataWithin$Ratio)
plotResiduals(resid.Within, form = dataWithin$WoA)
summary(model.Within)
parameters(model.Within)
plot(allEffects(model.Within))

dataWithin[, PredictWithin := predict(model.Within)]

ggplot(dataWithin, aes(x = Date, y = Similarity, color = RatioSplit)) +
  geom_jitter(size=2) + 
  geom_line(data = dataWithin[, mean(PredictWithin), by = .(Date, RatioSplit)], 
            aes(x = Date, y =V1, colour = RatioSplit),size=1.5) + 
  theme_classic(base_size = 18)+ 
  ylab("daily within-individual similarity")#+ 

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

#coeffeicient of variation for beteween individual variance
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

###### sleeping spot ############

#durations per zone per bird per day during dark hours
durDailyD = trackingData[Light == F, .(DurationNight = sum(Duration)), by = .(HenID, Zone, NightCycle)]
durDailyD = durDailyD[order(HenID),]
#extract maximum Zone for each bird per day
mainDurDailyD = durDailyD[, .SD[which.max(DurationNight)], by = .(HenID, NightCycle)]
mainDurDailyD[, onTop := ifelse(Zone == "Tier_4", 1, 0)]


###### wintergarden use #####
#all Wintergarten entries per bird
#careful: on vaccination days garten opened later!
dailyGarten = trackingData[Zone == "Wintergarten", .(Duration = Duration), by = .(HenID, Time, Date)]
#extract if hen goes out on day or not
inGarten = trackingData[, .(Out = ifelse(any(Zone == "Wintergarten"), 1, 0)), by = .(HenID, Date)]
#extract how long each hen went out per day
durDailyGarten = dailyGarten[, .(DurationGarten = sum(Duration)), by = .(HenID, Date)]
# latency to go out
latGarten = dailyGarten[, .(LatencyGarten = Time[1] - ymd_hms(paste(as_date(Time[1]), "10:00:00")), 
                            Time = Time[1]), by = .(HenID, Date)]

###### Time in nestbox zone ########
#all Nestbox entries per bird
dailyNest = trackingData[Light == T & Zone == "Ramp_Nestbox", .(Duration = Duration), by = .(HenID, Time, Date)]
#extract if hen was in nest zone on day or not
inNest = trackingData[Light == T, .(NestZone = ifelse(any(Zone == "Ramp_Nestbox"), 1, 0)), by = .(HenID, Date)]
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
varOfInterest = trackingData[Light == T, .(vertTravelDist = sum(distVertical)), by = .(HenID, Date, Pen, WoA)]
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
hist(as.numeric(varOfInterest$TimeNest))
model.Nest = lmer(as.numeric(TimeNest) ~ Ratio*WoA+ (1|Pen/HenID), data = varOfInterest)
model.Nest = lmer(as.numeric(TimeNest) ~ Ratio*WoA + (1|HenID),  data = varOfInterest)
resid.Nest = simulateResiduals(model.Nest, 1000)
plot(resid.Nest)
plotResiduals(resid.Nest, form = varOfInterest$Ratio[!is.na(varOfInterest$TimeNest)])
plotResiduals(resid.Nest, form = varOfInterest$DateID[!is.na(varOfInterest$TimeNest)])
summary(model.Nest)


#### Time series plots #########################
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