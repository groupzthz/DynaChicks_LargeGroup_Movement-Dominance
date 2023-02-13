###### Vertical distance ####

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

plotData = as.data.table(emmeans(model.Travel, ~ pairwise ~ poly(WoA,2)*Ratio, 
                                 at =  list(Ratio = round(quantile(varOfInterest$Ratio), digits = 2),
                                            WoA = c(20, 26, 30, 35, 40, 47, 51, 55)))$emmeans)


varOfInterest[, PredictTravel := predict(model.Travel)]

ggplot()+ 
  geom_jitter(data = varOfInterest, aes(x = WoA, y = vertTravelDist), height = 0.02, size = 1, alpha = 0.1)+
  geom_line(data = plotData, aes(x = WoA, y = emmean, group = Ratio,  colour = as.factor(Ratio)), size = 1)+
  geom_ribbon(data = plotData, aes(x = WoA, ymin = asymp.LCL, ymax = asymp.UCL, group = Ratio), alpha = 0.1)+
  #facet_grid(.~Pen)+
  theme_classic(base_size = 18)

ggplot(varOfInterest[WoA %in% quantile(WoA),], aes(x = Ratio, y = Similarity)) +
  geom_point(size=2) + 
  geom_smooth(data = dataWithin[WoA %in% quantile(WoA), mean(PredictWithin), by = .(WoA, Ratio)], 
              aes(y =V1), formula = y~x, method = lm, size=1.5, colour = "red") + 
  theme_classic(base_size = 18)+ 
  facet_grid(~WoA)+
  ylab("daily within-individual similarity")

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

