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
  setwd("/media/felix/DataDrive2/Documenten/PhD/08_Gene_Swamping/")
  load("2_data/6_EvolvedStrainsPosteriors/SummarisedPosteriorsEvolved.RData")
  load("2_data/5_EvolvedStrains/1_RawDensityData-Evolvedlines.RData")
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


#6) evolution model test without ancestor
{
  fullmodel <- lm(data = dd.evo, logr0.mean ~ Sex*Gene.Flow*Gradient, na.action = "na.fail")
  t <- dredge(fullmodel)
  best.model <- lm(data = dd.evo, logr0.mean ~ Sex + Gene.Flow + Gradient + Sex*Gene.Flow, na.action = "na.fail")
  anova(best.model)
  summary(best.model)
  Anova(best.model, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
}

#7) Make model predictions and plot evolution model without ancestor
{
  #Create prediction data
  prd.data <- expand.grid(Sex = unique(dd.evo$Sex), Gene.Flow = unique(dd.evo$Gene.Flow), Gradient = unique(dd.evo$Gradient))
  predictions <- predict(best.model, newdata = prd.data, se.fit = T)
  prd.data$meanr0 <-exp(predictions$fit)
  prd.data$upperr0 <- exp(predictions$fit + 1.96*predictions$se.fit)
  prd.data$lowerr0 <- exp(predictions$fit - 1.96*predictions$se.fit)
  
  #Plot data + predictions
  ggplot(dd.evo, aes(y=exp(logr0.mean), x=Sex, color = Sex)) + geom_point(size=2) + facet_grid(Gene.Flow~Gradient) + 
    geom_boxplot(inherit.aes = F, data = prd.data, mapping = aes(fill = Sex, colour =NA, middle = meanr0, ymax = upperr0, ymin = lowerr0, upper = upperr0, lower = lowerr0, x = Sex), stat = "identity", alpha = 0.3) + 
    geom_boxplot(inherit.aes = F, data = prd.data, mapping = aes(fill = NA, middle = meanr0, ymax = meanr0, ymin = meanr0, upper = meanr0, lower = meanr0, x = Sex), stat = "identity", alpha = 0.3) + 
    ylab(expression("Intrinsic rate of increase, r"[0]*" (1/h) ")) + xlab("Reproductive strategy") + 
    theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), axis.title.x = element_blank(),
                            axis.text.x=element_text(size=16, color = "black"), legend.position="none")+
    theme(strip.background =element_rect(fill="grey"))+
    theme(strip.text = element_text(colour = 'black')) +
    scale_color_manual(values=c("#8c510a", "#01665e"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual")) +
    scale_fill_manual(values=c("#8c510a", "#01665e"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual"))
  ggsave(filename = "4_results/Figures/03_growthRateEvolution_rawdata", device = "png", width = 8.5, height = 5.88, dpi = 300)
}

#8) Do model including ancestor data
{
  
  dd.evo <- select(dd.evo, Treatment, Gradient, Sex, Gene.Flow, logr0.mean, logalpha.mean, logK.mean, testpH)
  
  #Calculate mean growth rate per pH value
  sumdata.anc <- filter(sumdata, strain=="mix") %>% filter(pH %in% c(4, 6.5))
  sumdata.anc <- mutate(sumdata.anc, Treatment = "anc", Sex = "No sex", Gradient = ifelse(pH==6.5, "No gradient", "Gradient"), Gene.Flow = "No gene flow", testpH=pH)
  dd.anc <- select(sumdata.anc, Treatment, Gradient, Sex, Gene.Flow, logr0.mean, logalpha.mean, logK.mean, testpH)
  dd.anc.sum <- group_by(dd.anc, testpH) %>% summarise(mean.log.r0 = mean(logr0.mean))
  
  #Calculate log/log ratio of growth for evolved lines
  dd.evo$logratio <- log2(exp(dd.evo$logr0.mean) / exp(ifelse(dd.evo$testpH==4.0, dd.anc.sum[which(dd.anc.sum$testpH==4.0), ]$mean.log.r0[1], dd.anc.sum[which(dd.anc.sum$testpH==6.5), ]$mean.log.r0[1])))
  
  #Do model selection
  fullmodel2 <- lm(data = dd.evo, logratio ~ Sex*Gene.Flow*Gradient, na.action = "na.fail")
  t2 <- dredge(fullmodel2)
  model.avg(t2)
  best.model2 <- lm(data = dd.evo, logratio ~ Sex + Gene.Flow + Gradient + Sex*Gene.Flow, na.action = "na.fail")
  best.model3 <- lm(data = dd.evo, logratio ~ Sex + Gene.Flow + Gradient + Sex*Gene.Flow + Gradient*Gene.Flow, na.action = "na.fail")
  best.model4 <- lm(data = dd.evo, logratio ~ Sex + Gene.Flow + Gradient + Sex*Gene.Flow + Gradient*Sex, na.action = "na.fail")
  anova(best.model2)
  summary(best.model2)
  Anova(best.model2, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
  
  
}

#9) Make model predictions and plot evolution model with ancestor
{
  #Create prediction data
  prd.data2 <- expand.grid(Sex = unique(dd.evo$Sex), Gene.Flow = unique(dd.evo$Gene.Flow), Gradient = unique(dd.evo$Gradient))
  predictions2 <- predict(best.model2, newdata = prd.data2, se.fit = T)
  prd.data2$meanr0 <-predictions2$fit
  prd.data2$upperr0 <- predictions2$fit + 1.96*predictions2$se.fit
  prd.data2$lowerr0 <- predictions2$fit - 1.96*predictions2$se.fit

  
  #Plot data + predictions
  ggplot(dd.evo, aes(y=logratio, x=Sex, color = Sex)) + geom_point(size=2) + facet_grid(Gene.Flow~Gradient) +
    geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(fill = Sex, colour = NA, middle = meanr0, ymax = upperr0, ymin = lowerr0, upper = upperr0, lower = lowerr0, x = Sex), stat = "identity", alpha = 0.3) + 
    geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(fill = NA, middle = meanr0, ymax = meanr0, ymin = meanr0, upper = meanr0, lower = meanr0, x = Sex), stat = "identity", alpha = 0.3) + 
    ylab(expression("Intrinsic rate of increase, r"[0]*" (log-ratio response) ")) + xlab("Reproductive strategy") + 
    theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), axis.title.x = element_blank(),
                            axis.text.x=element_text(size=16, color = "black"))+
    theme(strip.background =element_rect(fill="grey"))+
    theme(strip.text = element_text(colour = 'black')) +
    scale_color_manual(values=c("#8c510a", "#01665e"), breaks=c("No sex (Asex)", "Sex (Sex)"), name="Reproduction", labels=c("Asexual", "Sexual")) +
    scale_fill_manual(values=c("#8c510a", "#01665e"), breaks=c("No sex (Asex)", "Sex (Sex)"), name="Reproduction", labels=c("Asexual", "Sexual"))
  ggsave(filename = "4_results/Figures/03_growthRateEvolution.png", device = "png", width = 8.5, height = 5.88, dpi = 300)
}

#Plot data + predictions
ggplot(dd.evo, aes(y=logratio, x=Sex, color = Sex)) + geom_point(size=2) + facet_grid(Gene.Flow~Gradient) +
  geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(fill = Sex, colour = NA, middle = meanr0, ymax = upperr0, ymin = lowerr0, upper = upperr0, lower = lowerr0, x = Sex), stat = "identity", alpha = 0.3) + 
  geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(fill = NA, middle = meanr0, ymax = meanr0, ymin = meanr0, upper = meanr0, lower = meanr0, x = Sex), stat = "identity", alpha = 0.3) + 
  ylab(expression("Intrinsic rate of increase, r"[0]*" (log-ratio response) ")) + xlab("Reproductive strategy") + 
  theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), axis.title.x = element_blank(),
                          axis.text.x=element_text(size=16, color = "black"), legend.position = "bottom")+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black')) +
  scale_color_manual(values=c("#8c510a", "#01665e"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual")) +
  scale_fill_manual(values=c("#8c510a", "#01665e"), breaks=c("Asexual", "Sexual"), name="Reproduction", labels=c("Asexual", "Sexual"))
ggsave(filename = "4_results/Figures/03_growthRateEvolution_legend.png", device = "png", width = 8.5, height = 5.88, dpi = 300)
