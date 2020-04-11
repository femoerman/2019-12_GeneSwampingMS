#1) Clear memory
rm(list=ls())

#2) Load packages
{
  library(tidyverse)
  library(MuMIn)
  library(nlme)
  library(car)
}

#3) Set working directory and load in data
{
  #setwd("F:/Documenten/PhD/08_Big_Experiment_1")
  setwd("/media/felix/DataDrive2/Documenten/PhD/08_Gene_Swamping")
  load("2_data/6_EvolvedStrainsPosteriors/SummarisedPosteriorsEvolved.RData")
  load("2_data/5_EvolvedStrains/1_RawDensityData.RData")
  load("2_data/6_EvolvedStrainsPosteriors/SummarisedPosteriorsAncestorRedoneAll.RData")
}
sumdata <- sumoutputall

#4) Combine posterior data with raw data and un-log transform posteriors
{
  dd <- dd %>% filter(hours==0) %>% arrange(curveID) %>% select(curveID, strain, Treatment, testpH, Gradient, Sex, Gene.Flow, ID)
  plotdata <- cbind(dd, sumoutput)
  
  plotdata$r0.mean <- ifelse(is.na(plotdata$logr0.mean), 0, exp(plotdata$logr0.mean))
  plotdata$K.mean <- ifelse(is.na(plotdata$logK.mean), NA, exp(plotdata$logK.mean))
  plotdata$alpha.mean <- ifelse(is.na(plotdata$logalpha.mean), NA, exp(plotdata$logalpha.mean))
  plotdata$d.mean <- ifelse(is.na(plotdata$logd.mean), NA, exp(plotdata$logd.mean))
  plotdata$Sex <- factor(ifelse(plotdata$Sex=="y", "Sexual", "Asexual"), levels = c("Asexual", "Sexual"))
  plotdata$Gradient <- factor(ifelse(plotdata$Gradient=="y", "Gradient", "Uniform"), levels = c("Uniform", "Gradient"))
  plotdata$Gene.Flow <- factor(ifelse(plotdata$Gene.Flow=="y", "Present", "Absent"), levels = c("Absent", "Present"))
  
  sumdata$testpH <- paste(sumdata$pH, "j")
}

#5) Further prepare data
{
  
  #Filter for evolution data and pH 4 or 6.5
  dd.evo <- filter(plotdata, strain=="mix") %>% filter(testpH %in% c(4, 6.5)) %>% filter((testpH==4 & Gradient == "Gradient") | (testpH!=4 & Gradient != "Gradient"))
  sumdata.evo <- filter(sumdata, strain=="mix") %>% filter(pH %in% c(4, 6.5))
}




#Do analysis with all r0 data?
{
  dd.all <- filter(plotdata, strain=="mix")
  ggplot(dd.all, aes(x= exp(-testpH), y = exp(logr0.mean), colour = Sex, group = ID)) + geom_point() + geom_line() + facet_grid(Gradient~Gene.Flow)
  ggplot(dd.all, aes(x=testpH, y = exp(logr0.mean), colour = Sex, group = ID)) + geom_point() + geom_line() + facet_grid(Gradient~Gene.Flow)
  
  #Standardize to ancestor
  {
    dd.anc.sum <- filter(sumdata, strain=="mix") %>% group_by(pH) %>% summarize(meangrowth = mean(r0mean))
    
    #Calculate log/log ratio of growth for evolved lines
    dd.all$logratio <- NA
    for (i in 1:nrow(dd.all)){
      dd.all[i, ]$logratio <- log2(exp(dd.all[i, ]$logr0.mean)/filter(dd.anc.sum, pH==dd.all[i, ]$testpH)$meangrowth[1])
    }
    dd.all <- filter(dd.all, !is.na(logratio))
    
  }
}

#Create factorial variable for pH of the assay medium
dd.all$pHfact <- as.factor(dd.all$testpH)

#Fit model, tetspH
full.model <- lme(data = dd.all, logratio ~ Gradient*Sex*Gene.Flow*pHfact, random =  ~ 1|ID, method = "ML")
t.all <- dredge(full.model, rank = BIC)
best.model <- lme(data = dd.all, logratio ~ Gradient + pHfact, random =  ~ 1|ID, method = "REML")
anova(best.model)
summary(best.model)
Anova(best.model, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)

#Create prediction data
prd.data2 <- expand.grid(Gradient = unique(dd.all$Gradient),pHfact =unique(dd.all$pHfact),
                         Sex = unique(dd.all$Sex), Gene.Flow = unique(dd.all$Gene.Flow))
predictions2 <- predict(best.model, newdata = prd.data2, se.fit = T, level = 0)
prd.data2$meansize <-predictions2$fit
prd.data2$uppersize <- predictions2$fit + 1.96*predictions2$se.fit
prd.data2$lowersize <- predictions2$fit - 1.96*predictions2$se.fit
prd.data2$pHfact <- as.numeric(as.character(prd.data2$pHfact))

#Plot model predictions with data
ggplot(dd.all, aes(x = testpH, y = logratio, colour = Gradient)) + geom_point(size=2) + 
  geom_ribbon(data = prd.data2, inherit.aes = F, mapping = aes(x = pHfact, ymax = uppersize, ymin = lowersize, fill = Gradient), alpha = 0.3) + 
  geom_line(data = prd.data2, inherit.aes = F, mapping = aes(x = pHfact, y = meansize, colour = Gradient))+
  ylab(expression("Intrinsic rate of increase, r"[0]*" (log-ratio response) ")) + xlab("Assay medium pH") + 
  theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5), legend.position = "none")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black')) +
  scale_color_manual(values=c("#2d33f9", "#e41a1c"), breaks=c("Uniform", "Gradient"), name="Abiotic conditions", labels=c("Uniform", "Gradient")) +
  scale_fill_manual(values=c("#2d33f9","#e41a1c"), breaks=c("Uniform", "Gradient"), name="Abiotic conditions", labels=c("Uniform", "Gradient"))+
  facet_grid(Gene.Flow~Sex)
ggsave(filename = "4_results/Figures/06_growthRateAllpH", device = "png", width = 8.5, height = 5.88, dpi = 300)
