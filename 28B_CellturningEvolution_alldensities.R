#1) Clear memory
{
  rm(list=ls())
}

#2) load libraries
{
  library(tidyverse)
  library(readr)
  library(MuMIn)
  library(nlme)
}

#3) Set working directory
{
  setwd("/media/felix/DataDrive2/Documenten/PhD/08_Gene_Swamping/")
  # setwd("F:/Documenten/PhD/08_Big_Experiment_1")
  # setwd("E:/PhD/08_Big_Experiment_1")
}

#4) read in data
{
  #For evolved lines
  timepoints <- c( "t1", "t2", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12")
  load("2_data/5_EvolvedStrains/t0/Population_Data.RData")
  dd <- rbind(pop_output, pop_output, pop_output, pop_output, pop_output, pop_output, pop_output, pop_output)
  dd$hours <- 0
  dd$timepoint <- 't0'
  dd$testpH <- c(rep(6.5, 49), rep(6.0, 49), rep(5.5, 49), rep(5, 49), rep(4.5, 49), rep(4, 49), rep(3.5, 49), rep(3, 49))
  for (p in timepoints){
    path <- paste("2_data/5_EvolvedStrains/", p, "/Population_Data.RData", sep="")
    load(path)
    if (p=="t5"){
      pop_output$time <- pop_output$hours
    }
    if (p=="t1"){
      pop_output$time <- c(rep(8, 36), rep(9, 46), 10, rep(9, 40), rep(10, 73), rep(11, 36),rep(12, 13), rep(11, 49), rep(9, 12), rep(10, 25), rep(11, 12), rep(8, 25), rep(9, 24))
    }
    pop_output$timepoint <- p
    dd <- rbind(dd, pop_output)
    
    #Do the same for the ancestral lines
    timepoints <- c("t1", "t2", "t4", "t5", "t6", "t7", "t8", "t9", "t10", "t11", "t12")
    load("2_data/1_Ancestors/t0/Population_Data.RData")
    ddanc <- pop_output
    ddanc$timepoint <- 't0'
    for (p in timepoints){
      path <- paste("2_data/1_Ancestors/", p, "/Population_Data.RData", sep="")
      load(path)
      pop_output$timepoint <- p
      ddanc <- rbind(ddanc, pop_output)
    }
    
    #Read in summarized data and filter
    #Read in summarized ancestor data
    load("2_data/summarizedAncestorData.RData")
    sumdata.anc <-sumdata %>% filter(strain=="mix")
    #Read in summarized evolved data
    load("2_data/6_EvolvedStrainsPosteriors/SummarisedPosteriorsEvolved.RData")
    sumdata <-sumoutput
    
    ddanc$ID <- paste("A", as.numeric(ddanc$file), sep="")
    ddanc$testpH <- ddanc$pH
    ddanc$Treatment <- "Anc"
    ddanc$Gradient <- "No gradient"
    ddanc$Sex <- "No sex"
    ddanc$Gene.Flow <- "No gene flow"
    ddanc$date <- 0
    ddanc$time <- 0
  }
}

#5) Filter data for Evolution lines and paste together
{
  #Filter data and select correct variables
  ddanc <- filter(ddanc, strain=="mix") %>% select(hours, testpH, indiv_per_volume, major_mean, gross_speed_mean, sd_turning_mean, ID, testpH, dilution, Treatment,
                                                   Gradient, Sex, Gene.Flow, date, time, timepoint, minor_mean)
  ddevo <- filter(dd, strain=="mix") %>% select(hours, testpH, indiv_per_volume, major_mean, gross_speed_mean, sd_turning_mean, ID, testpH, dilution, Treatment,
                                                Gradient, Sex, Gene.Flow, date, time, timepoint, minor_mean)
  
  #Bind together
  dd.all <- rbind(ddanc, ddevo)
}

#6) Data managing
{
  #6.1) Define hours variable
  dd.all$curveID <- paste(dd.all$ID, dd.all$testpH)
  dd.all$days<-  ifelse(dd.all$Treatment=="EVO", as.integer(as.Date(dd.all$date, "%d.%m.%y"))-17868, 0)
  dd.all$hours <- ifelse(dd.all$Treatment=="EVO", dd.all$days*24+dd.all$time - ifelse(dd.all$ID<23, 17, 18), dd.all$hours)
  
  #6.2) create indiv_per_ml variable
  dd.all$indiv_per_ml <- dd.all$indiv_per_volume*dd.all$dilution*1000
  
  #6.3) Remove missing datapoints
  dd.all <- dd.all[complete.cases(dd.all), ]
  
  #6.4) Rename variables properly
  dd.all$Sex <- factor(ifelse(dd.all$Sex=="y", "Sex", "No sex"), levels = c("No sex", "Sex"))
  dd.all$Gradient <- factor(ifelse(dd.all$Gradient=="y", "Gradient", "No gradient"), levels = c("No gradient", "Gradient"))
  dd.all$Gene.Flow <- factor(ifelse(dd.all$Gene.Flow=="y", "Gene flow", "No gene flow"), levels = c("No gene flow", "Gene flow"))
  
  
  #6.4) calculate t0, inflection point and K point
  #Make dataset for t1
  dd.t0 <- filter(dd.all, timepoint=="t1")
  
  #6.5) Redefine ID for sumdata.anc, select variables and bind together
  #Rename variables so they correspond
  sumdata.anc$ident <- unique(paste(ddanc$ID, ddanc$testpH))
  sumdata.anc$logr0.mean <- log(sumdata.anc$r0mean)
  sumdata.anc$logalpha.mean <- log(sumdata.anc$alfamean)
  sumdata.anc$logK.mean <- sumdata.anc$log_Kmean
  sumdata.anc$logd.mean <- sumdata.anc$log_dmean
  
  #6.6) Select the necessary variables
  sumdata.anc <- select(sumdata.anc, ident, logr0.mean, logd.mean, logalpha.mean, logK.mean)
  sumdata.evo <- filter(sumdata, ident %in% unique(dd.all$curveID))  %>% select(ident, logr0.mean, logd.mean, logalpha.mean, logK.mean)
  
  #6.7) Bind together
  sumdata <- rbind(sumdata.anc, sumdata.evo)
  
  #6.8) Calculate mid log phase and K for all populations
  {
    sumdata <- drop_na(sumdata)
    dd.K <- data.frame()
    dd.inf <- data.frame()
    
    
    row.names(sumdata) <- sumdata$ident
    for (i in sumdata$ident){
      tempdata <- filter(dd.all, curveID==i)
      timeK <- min(filter(tempdata, indiv_per_ml>=0.99*exp(sumdata[i, ]$logK.mean))$hours)
      timeK <- ifelse(timeK==Inf, tempdata[which.max(tempdata$indiv_per_ml), ]$hours, timeK)
      if(timeK> 0){
        timeInflect <- tempdata[which.min((tempdata[which(tempdata$hours<timeK), ]$indiv_per_ml - exp(sumdata[i, ]$logK.mean)/2)^2), ]$hours
        dd.inf <- rbind(dd.inf, filter(tempdata, hours==timeInflect))
      }
      dd.K <- rbind(dd.K, filter(tempdata, hours==timeK))
    }
  }
  
  #6.9) remove datapoints where no good datapoint was available (i.e. where population crashed and never reached K or mid log phase)
  dd.K <- na.omit(dd.K)
  dd.inf <- na.omit(dd.inf)
  
  #6.10) Create %K variable
  dd.evo <-filter(dd.all, Treatment=="EVO")%>% na.omit()
  {
    dd.evo$percK <- NA
    for (i in 1:nrow(dd.evo)){
      K <- exp(filter(sumdata.evo, ident == dd.evo[i, "curveID"])$logK.mean)
      dd.evo[i, "percK"] <- dd.evo[i, "indiv_per_ml"]/K
    }
  }
  dd.evo <- na.omit(dd.evo)
  dd.evo$percK <- as.numeric(dd.evo$percK)
  
  
}

ggplot(dd.evo, aes(x = percK, y = sd_turning_mean, group = Gradient, colour = Gradient)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(~testpH)



#7) Fit model over all data for size - evolved only
{
  full.model.speed.evo <- lme(data = dd.evo, sd_turning_mean ~ percK*testpH*Gradient*Sex +percK*testpH*Gradient*Gene.Flow +percK*testpH*Sex*Gene.Flow + percK*Gradient*Sex*Gene.Flow + testpH*Gradient*Sex*Gene.Flow, random =  list(ID=~1, curveID=~1), method = "ML")
  
  #Do manual comparison
  {
    full.model.speed.evo <- lme(data = dd.evo, sd_turning_mean ~ percK*testpH*Gradient*Sex*Gene.Flow, random =  list(ID=~1, curveID=~1), method = "ML")
    t1 <-  dredge(lme(data = dd.evo, sd_turning_mean ~ testpH*Gradient*Sex*Gene.Flow, random =  list(ID=~1, curveID=~1), method = "ML"), rank = BIC)
    t2 <-  dredge(lme(data = dd.evo, sd_turning_mean ~ percK*Gradient*Sex*Gene.Flow, random =  list(ID=~1, curveID=~1), method = "ML"), rank = BIC)
    t3 <-  dredge(lme(data = dd.evo, sd_turning_mean ~ percK*testpH*Sex*Gene.Flow, random =  list(ID=~1, curveID=~1), method = "ML"), rank = BIC)
    t4 <-  dredge(lme(data = dd.evo, sd_turning_mean ~ percK*testpH*Gradient*Sex, random =  list(ID=~1, curveID=~1), method = "ML"), rank = BIC)
    t5 <-  dredge(lme(data = dd.evo, sd_turning_mean ~ percK*testpH*Gradient*Gene.Flow, random =  list(ID=~1, curveID=~1), method = "ML"), rank = BIC)
    BIC(full.model.speed.evo)
    t1[1, "BIC"]
    t2[1, "BIC"]
    t3[1, "BIC"]
    t4[1, "BIC"]
    t5[1, "BIC"]
  }
  
  #best model
  best.model.speed.evo <- lme(data = dd.evo, sd_turning_mean ~ testpH*percK , random =  list(ID=~1, curveID=~1))
  anova(best.model.speed.evo)
  summary(best.model.speed.evo)
  Anova(best.model.speed.evo, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
}

#8) Make predictions and plot + raw data
{
  #Create prediction data
  prd.data1 <- expand.grid(Gradient = unique(dd.evo$Gradient), testpH = seq(from = min(dd.evo$testpH), to = max(dd.evo$testpH), length.out = 7),
                           Sex = unique(dd.evo$Sex), Gene.Flow = unique(dd.evo$Gene.Flow), percK = seq(from=0, to=2, length.out = 100))
  predictions1 <- predict(best.model.speed.evo, newdata = prd.data1, se.fit = T, level = 0)
  prd.data1$meanspeed <-predictions1$fit
  prd.data1$upperspeed <- predictions1$fit + 1.96*predictions1$se.fit
  prd.data1$lowerspeed <- predictions1$fit - 1.96*predictions1$se.fit
  
  #Plot data and predictions
  ggplot(dd.evo, aes(x = percK, y = sd_turning_mean, colour = Gradient)) + geom_point(size=1) + 
    geom_ribbon(data = prd.data1, inherit.aes = F, mapping = aes(x = percK, ymax = upperspeed, ymin = lowerspeed, fill = Gradient), alpha = 0.3) + 
    geom_line(data = prd.data1, inherit.aes = F, mapping = aes(x = percK, y = meanspeed, colour = Gradient)) + 
    ylab("Turning speed (rad/s)") + xlab("Density (proportion of K)") +
    theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5), legend.position = "none", axis.text.x = element_text(angle = 90))+
    theme(strip.background =element_rect(fill="grey"))+
    theme(strip.text = element_text(colour = 'black')) +
    scale_color_manual(values=c("#2d33f9", "#e41a1c"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient", "Gradient")) +
    scale_fill_manual(values=c("#2d33f9", "#e41a1c"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient", "Gradient")) + 
    facet_grid(~testpH)
  ggsave(filename = "4_results/Figures/09_TurningEvolutionAllDensities", device = "png", width = 9.53, height = 5.88, dpi = 300)
}
