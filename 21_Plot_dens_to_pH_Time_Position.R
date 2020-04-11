#Clear memory
rm(list=ls())

#load libraries
library(ggplot2)
library(tidyverse)
library(scales)
library(nlme)
library(MuMIn)
library(ggpubr)
library(car)

#Set working directory
#setwd("E:/PhD/08_Big_Experiment_1")
setwd("/media/felix/DataDrive2/Documenten/PhD/08_Gene_Swamping/")

#read in data
load(file = "2_data/EvolutionTimeSeries.RData")

#Define indiv_per_ml
dd$indiv_per_ml <- dd$indiv_per_volume*dd$dilution*1000
#Fix other error in dates
dd$date <- ifelse(dd$date=="25.10.18", "24.10.18", as.character(dd$date))
#Calculate number of days since start
dd$days <- as.numeric(as.Date(dd$date, format="%d.%m.%y")-as.Date("28.09.18", format="%d.%m.%y")) + 2

#Create dataset for % dispersers and densities
dd2 <- dd[seq(1, length(dd$"file"), 2),c(2:9, 11, 26)]
dd2$perc.disp <- dd[seq(1, length(dd$file), 2), ]$indiv_per_ml/dd[seq(2, length(dd$file), 2), ]$indiv_per_ml
dd2$dens.front <- dd[seq(1, length(dd$file), 2), ]$indiv_per_ml+dd[seq(2, length(dd$file), 2), ]$indiv_per_ml

#Replace infinites by NA
dd2[dd2 == Inf] <- NA
dd.evo <- filter(dd2, Treatment == "EVO")

#3) read in expansion data
{
  dd3 <- read_tsv("2_data/0_b_PositionDuringExperiment.txt")
  dd4 <- read_tsv("2_data/6_Survival_Data.csv")
  dd5 <- read_tsv("2_data/0_c_pHDuringExperiment.txt")
}

#4) Data handling
{
  #4.1) Create dataframe that is easy to handle
  days <- as.integer(as.Date(colnames(dd3)[2:24], "%m/%d/%y"))
  days <- days - days[1] - 2
  colnames(dd3) <- c("Number", 2:24)
  
  #4.2) Get all positions in one row
  positions <- data.frame(row.names = 1:(104*23))
  positions$Culture <- rep(1:104, 23)
  positions$Day <- NA
  positions$pos <- NA
  positions$pH <- NA
  for (i in 1:23){
    start <- (i-1)*104+1
    positions[start:(start+104-1), ]$Day <- rep(days[i], 104)
    positions[start:(start+104-1), ]$pos <- pull(dd3[, i+1])
    positions[start:(start+104-1), ]$pH <- pull(dd5[, i+1])
  }
  
  #5.3) copy other variables
  positions$Sex <- dd4$sex
  positions$Strain <- dd4$Strain
  positions$Gradient <- dd4$Gradient
  positions$GeneFlow <- dd4$`Gene flow`
  positions$Sex <- ifelse(positions$Sex == "y", "y", "n")
  
  #5.4) Create treatment variable
  positions$Treatment <- as.factor(paste(positions$Strain, positions$Gradient, positions$`Gene flow`, positions$Sex))
  positions$Strain <- ifelse(positions$Strain=="CU427.4", "clone 1", ifelse(positions$Strain=="CU428.2", "clone 2", 
                                                                            ifelse(positions$Strain=="SB3539", "clone 3", ifelse(positions$Strain=="B2086.2", "clone 4", "mix"))))
}

#6) Filter for evolution and create common identifier
{
  dens.evo <- filter(dd2, Treatment == "EVO")
  pos.evo <- filter(positions, Strain == "mix")
  dens.evo$ident <- paste(dens.evo$ID, dens.evo$days)
  pos.evo$ident <- paste(pos.evo$Culture, pos.evo$Day)
  
  #Filter pos.evo based on dens.evo entries
  pos.evo2 <- filter(pos.evo, ident %in% dens.evo$ident)
  dens2.evo <- filter(dens.evo, ident %in%pos.evo2$ident)
}


#7)Combine in one nice dataset
{
  data.all <- dens2.evo
  data.all$pos <- pos.evo2$pos
  
  #Change varables to nicer names
  data.all$Sex <- factor(ifelse(data.all$Sex=="y", "Sex (Sex)", "No sex (Asex)"), levels = c("No sex (Asex)", "Sex (Sex)"))
  data.all$Gradient <- factor(ifelse(data.all$Gradient=="y", "Gradient", "No gradient"), levels = c("No gradient", "Gradient"))
  data.all$Gene.Flow <- factor(ifelse(data.all$Gene.Flow=="y", "Gene flow (GF)", "No gene flow (NoGF)"), levels = c("No gene flow (NoGF)", "Gene flow (GF)"))
  data.all$pHcalc <- pos.evo2$pH
}

#9) Fit model for density in terms of Position, sex, gradient and gene flow
{
  full.model <- lme(data=data.all, dens.front~Gradient*Sex*Gene.Flow*pos, random =  ~ 1|ID, method = "ML")
  comp <- dredge(full.model, rank = BIC)
  best.model <-  lme(data=data.all, dens.front~Gradient*pos, random =  ~ 1|ID, method = "REML")
  r.squaredLR(best.model)
  anova(best.model)
  summary(best.model)
  Anova(best.model, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
  
}

#Redo figure with plotted model predictions
#Define prediction data
d.Grad <- filter(data.all, Gradient==unique(data.all$Gradient)[1])
d.const <- filter(data.all, Gradient==unique(data.all$Gradient)[2])
prd.data1 <- expand.grid(pos = seq(min(d.Grad$pos), max(d.Grad$pos), by = 1), Gradient = unique(d.Grad$Gradient))
prd.data2 <- expand.grid(pos = seq(min(d.const$pos), max(d.const$pos), by = 1), Gradient = unique(d.const$Gradient))
prd.data.all <- rbind(prd.data1, prd.data2)
predictions1 <- predict(newdata = prd.data.all, best.model, se.fit = T, level = 0)
prd.data.all$mean <- predictions1$fit
prd.data.all$upper <- predictions1$fit + 1.96*predictions1$se.fit
prd.data.all$lower <- predictions1$fit - 1.96*predictions1$se.fit

#Colour Sex events black
data.all$SexEvent <- ifelse(data.all$date %in% c("10.10.18", "24.10.18", "07.11.18", "21.11.18"), "Centrifugation event", "No centrifugation event")
data.all$colour = ifelse(data.all$SexEvent=="Centrifugation event", "Centrifugation event", as.character(data.all$Gradient))

library(ggstance)

ggplot(data.all, aes(y=dens.front, x=pos, colour=as.character(Gradient) )) + geom_point(size=2) + facet_grid(Sex~Gene.Flow) + geom_line(aes(group=as.factor(ID))) + 
  xlab("Space (patches)") + geom_ribbon(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, ymin = lower, ymax = upper, fill=as.character(Gradient)), alpha = 0.3) +
  geom_line(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, y = mean, colour=Gradient), size=2) + 
  ylab("Population density, N (indiv/mL)")+ scale_y_continuous(labels = scientific) + 
  theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5))+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black')) +
  scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)")) +
  scale_fill_manual(values=c("#e41a1c","#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)"))

ggsave(filename = "4_results/Figures/02_DensityDuringEvolutionB", device = "png", width = 9, height = 5.88, dpi = 300)

ggplot(data.all, aes(x = pos, y = Gradient)) + geom_boxploth() + facet_grid(Sex~Gene.Flow) 
ggplot(data.all, aes(x = pos, y = perc.disp, colour = Gradient)) + geom_point() + facet_grid(Sex~Gene.Flow)  + geom_smooth(method = "lm")

#10) Remake figure expansion distance boxplot
{
  #Read in survival data
  surv <- read_tsv("2_data/6_Survival_Data.csv")
  
  #Do a linear model to assess expansion distance for surviving populations
  {
    #Get expansion distances for surviving populations, on the last day
    dd.dist <- filter(positions, Strain=="mix", Day==max(positions$Day)) %>% filter(Culture%in% filter(surv, Survival==1)$Number)
    dd.dist$Sex <- factor(ifelse(dd.dist$Sex=="y", "Sex (Sex)", "No sex (Asex)"), levels = c("No sex (Asex)", "Sex (Sex)"))
    dd.dist$Gradient <- factor(ifelse(dd.dist$Gradient=="y", "Gradient", "No gradient"), levels = c("No gradient", "Gradient"))
    dd.dist$GeneFlow <- factor(ifelse(dd.dist$GeneFlow=="y", "Gene flow (GF)", "No gene flow (NoGF)"), levels = c("No gene flow (NoGF)", "Gene flow (GF)"))
    
    #Perform the linear model
    fullmodel <- lm(pos~Gradient*GeneFlow*Sex, dd.dist, na.action = "na.fail")
    t <- dredge(fullmodel)
    bestmodel <- lm(pos ~Gradient, dd.dist)
    anova(bestmodel)
    summary(bestmodel)
    Anova(bestmodel, contrasts=list(topic=contr.sum, sys=contr.sum), type=3)
    
    #Prepare the prediction data
    prd.datadist <- expand.grid(Gradient = unique(dd.dist$Gradient))
    predictions <- predict(bestmodel, newdata = prd.datadist, se.fit = T)
    prd.datadist$mean <- predictions$fit
    prd.datadist$upper <- predictions$fit + 1.96*predictions$se.fit
    prd.datadist$lower <- predictions$fit - 1.96*predictions$se.fit
  }
  
  #Remake the original figure
  ggplot(dd.dist, aes(x=Gradient, y=pos, colour=Gradient)) + xlab("Gradient treatment") + ylab("Expansion distance") + facet_grid(Sex~GeneFlow) +
    geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(x = Gradient, lower = lower, upper = upper, middle = mean, ymin = lower, ymax = upper, fill = Gradient, colour=NA), stat = "identity", alpha = 0.3) +
    geom_boxplot(inherit.aes = F, data = prd.data2, mapping = aes(x = Gradient, lower = mean, upper = mean, middle = mean, ymin = mean, ymax = mean, fill = NA), stat = "identity", alpha = 0.3) +
    geom_point(position = position_dodge2(width = 0.5), size=3) + 
    theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5))+
    theme(strip.background =element_rect(fill="grey"))+
    theme(strip.text = element_text(colour = 'black')) +
    scale_color_manual(values=c("#2d33f9", "#e41a1c"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)")) +
    scale_fill_manual(values=c("#2d33f9","#e41a1c"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)"))
  
}

#11) Combine both plots in one nice figure
{
  #Create colour variable
  data.all$colourscheme <- as.character(ifelse(data.all$ID %in% c(3, 12, 13, 15, 16, 17, 19), "Gradient (extinct)", ifelse(data.all$Gradient == "No gradient", "Uniform (surviving)", "Gradient (surviving)")))
  
  prd.data.all$colourscheme<- as.character(ifelse(prd.data.all$Gradient=="Gradient", "Gradient (surviving)", "Uniform (surviving)"))
  
  #Add for the extinct populations density=0 for the last timepoint to show final density in plot
  temp <- filter(data.all, days==min(days), ID %in% c(3, 12, 15, 16, 17, 19))
  temp$days <- 61
  temp$dens.front<-0
  for (i in temp$ID){
    temp[which(temp$ID==i), "pos"] <- max(filter(data.all, ID==i)$pos)
  }
  data.all <- rbind(data.all, temp)
  
  #Create all figures separately for densiy during range expansion
  {
    #1) Sex and gene flow
    fig1.1 <- ggplot(filter(data.all, Sex == "Sex (Sex)", Gene.Flow == "Gene flow (GF)"), aes(y=dens.front, x=pos, colour=colourscheme )) + geom_point(size=2, alpha = 0.3) + geom_line(aes(group=as.factor(ID)), alpha = 0.3) + 
      geom_point(inherit.aes = F, data = filter(data.all, days==61, Sex == "Sex (Sex)", Gene.Flow == "Gene flow (GF)"), mapping = aes(y=dens.front, x=pos, colour=colourscheme ), size=3)+
      xlab("Space (patches)") + geom_ribbon(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, ymin = lower, ymax = upper, fill=colourscheme), alpha = 0.3) +
      geom_line(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, y = mean, colour=colourscheme), size=2) + 
      ylab("Population density, N\n (indiv/mL)")+ scale_y_continuous(labels = scientific, limits = c(-0.2e5, 3.3e5)) + guides(fill=F)+ xlim(3, 24)+ 
      theme_light() +   theme(legend.position="none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5),
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())+
      scale_color_manual(values=c("#000000", "#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Gradient (extinct)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Gradient (extinct)", "Uniform (surviving)")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Uniform (surviving)"))
    
    #2) Sex and no gene flow
    fig2.1 <- ggplot(filter(data.all, Sex == "Sex (Sex)", Gene.Flow != "Gene flow (GF)"), aes(y=dens.front, x=pos, colour=colourscheme )) + geom_point(size=2, alpha = 0.3) + geom_line(aes(group=as.factor(ID)), alpha = 0.3) + 
      geom_point(inherit.aes = F, data = filter(data.all, days==61, Sex == "Sex (Sex)", Gene.Flow != "Gene flow (GF)"), mapping = aes(y=dens.front, x=pos, colour=colourscheme ), size=3)+
      xlab("Space (patches)") + geom_ribbon(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, ymin = lower, ymax = upper, fill=colourscheme), alpha = 0.3) +
      geom_line(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, y = mean, colour=colourscheme), size=2) + 
      ylab("Population density, N\n (indiv/mL)")+ scale_y_continuous(labels = scientific, limits = c(-0.2e5, 3.3e5)) + guides(fill=F)+ xlim(3, 24)+ 
      theme_light() +   theme(legend.position="none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5),
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())+
      scale_color_manual(values=c("#000000", "#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Gradient (extinct)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Gradient (extinct)", "Uniform (surviving)")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Uniform (surviving)"))
    
    #3) No sex and gene flow
    fig3.1 <- ggplot(filter(data.all, Sex != "Sex (Sex)", Gene.Flow == "Gene flow (GF)"), aes(y=dens.front, x=pos, colour=colourscheme )) + geom_point(size=2, alpha = 0.3) + geom_line(aes(group=as.factor(ID)), alpha = 0.3) + 
      geom_point(inherit.aes = F, data= filter(data.all, days==61, Sex != "Sex (Sex)", Gene.Flow == "Gene flow (GF)"), mapping = aes(y=dens.front, x=pos, colour=colourscheme ), size=3)+
      xlab("Space (patches)") + geom_ribbon(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, ymin = lower, ymax = upper, fill=colourscheme), alpha = 0.3) +
      geom_line(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, y = mean, colour=colourscheme), size=2) + 
      ylab("Population density, N\n (indiv/mL)")+ scale_y_continuous(labels = scientific, limits = c(-0.2e5, 3.3e5)) + guides(fill=F)+ xlim(3, 24)+ 
      theme_light() +   theme(legend.position="none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5),
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())+
      scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Uniform (surviving)")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Uniform (surviving)"))
    
    #4) No sex and no gene flow
    fig4.1 <- ggplot(filter(data.all, Sex != "Sex (Sex)", Gene.Flow != "Gene flow (GF)"), aes(y=dens.front, x=pos, colour=colourscheme )) + geom_point(size=2, alpha = 0.3) + geom_line(aes(group=as.factor(ID)), alpha = 0.3) + 
      geom_point(inherit.aes = F, data = filter(data.all, days==61, Sex != "Sex (Sex)", Gene.Flow != "Gene flow (GF)"), mapping = aes(y=dens.front, x=pos, colour=colourscheme ), size=3)+
      xlab("Space (patches)") + geom_ribbon(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, ymin = lower, ymax = upper, fill=colourscheme), alpha = 0.3) +
      geom_line(inherit.aes=F, data = prd.data.all, mapping = aes(x = pos, y = mean, colour=colourscheme), size=2) + 
      ylab("Population density, N\n (indiv/mL)")+ scale_y_continuous(labels = scientific, limits = c(-0.2e5, 3.3e5)) + guides(fill=F)+ xlim(3, 24)+ 
      theme_light() +   theme(legend.position="none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5),
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())+
      scale_color_manual(values=c("#000000", "#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Gradient (extinct)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Gradient (extinct)", "Uniform (surviving)")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("Gradient (surviving)", "Uniform (surviving)"), name="Abiotic conditions", labels=c("Gradient (surviving)", "Uniform (surviving)"))
    
  }
  
  #Create all figures separately for expansion distance
  {
    dd.dist$col <- ifelse(dd.dist$Gradient=="Gradient", "Gradient", "Uniform")
    prd.datadist$col <- ifelse(prd.datadist$Gradient=="Gradient", "Gradient", "Uniform")
    #1) Sex and gene flow
    fig1.2 <- ggplot(filter(dd.dist, Sex == "Sex (Sex)", GeneFlow == "Gene flow (GF)"), aes(y=col, x=pos, colour=col)) + xlab("Space (patches)") + ylab("Abiotic conditions") +
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = lower, xupper=upper, xmiddle=mean, xmin = lower, xmax = upper, fill = Gradient, colour = NA), stat = "identity", alpha = 0.3) + 
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = mean, xupper=mean, xmiddle=mean, xmin = mean, xmax = mean, fill = NA), stat = "identity", alpha = 0.3) + 
      geom_point(position = position_dodge2(width = 0.5), size=3) +xlim(3, 24)+
      theme_light() +   theme(legend.position = "none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5),
                              axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())+
      scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) 
    
    #2) Sex and no gene flow
    fig2.2 <- ggplot(filter(dd.dist, Sex == "Sex (Sex)", GeneFlow != "Gene flow (GF)"), aes(y=col, x=pos, colour=col)) + xlab("Space (patches)") + ylab("Abiotic conditions") +
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = lower, xupper=upper, xmiddle=mean, xmin = lower, xmax = upper, fill = Gradient, colour = NA), stat = "identity", alpha = 0.3) + 
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = mean, xupper=mean, xmiddle=mean, xmin = mean, xmax = mean, fill = NA), stat = "identity", alpha = 0.3) + 
      geom_point(position = position_dodge2(width = 0.5), size=3) + xlim(3, 24)+
      theme_light() +   theme(legend.position = "none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5),
                              axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())+
      scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) 
    
    #3) No sex and gene flow
    fig3.2 <- ggplot(filter(dd.dist, Sex != "Sex (Sex)", GeneFlow == "Gene flow (GF)"), aes(y=col, x=pos, colour=col)) + xlab("Space (patches)") + ylab("Abiotic conditions") +
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = lower, xupper=upper, xmiddle=mean, xmin = lower, xmax = upper, fill = Gradient, colour = NA), stat = "identity", alpha = 0.3) + 
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = mean, xupper=mean, xmiddle=mean, xmin = mean, xmax = mean, fill = NA), stat = "identity", alpha = 0.3) + 
      geom_point(position = position_dodge2(width = 0.5), size=3) + xlim(3, 24)+
      theme_light() +   theme(legend.position = "none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5))+
      scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) 
    
    #4) No sex and no gene flow
    fig4.2 <- ggplot(filter(dd.dist, Sex != "Sex (Sex)", GeneFlow != "Gene flow (GF)"), aes(y=col, x=pos, colour=col)) + xlab("Space (patches)") + ylab("Abiotic conditions") +
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = lower, xupper=upper, xmiddle=mean, xmin = lower, xmax = upper, fill = Gradient, colour = NA), stat = "identity", alpha = 0.3) + 
      geom_boxploth(data = prd.datadist, inherit.aes = F, mapping = aes(y = col, xlower = mean, xupper=mean, xmiddle=mean, xmin = mean, xmax = mean, fill = NA), stat = "identity", alpha = 0.3) + 
      geom_point(position = position_dodge2(width = 0.5), size=3) + xlim(3, 24)+
      theme_light() +   theme(legend.position = "none", axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                              axis.title=element_text(size=14), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5))+
      scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) +
      scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) 
    
  }
  
  #Combine 2 by 2 the figures for expansion distance and density during evolution
  {
    fig4 <- ggarrange(fig1.1, fig1.2, nrow=2, align = "v", heights = c(2, 1), labels = c("F", "H"), label.x = 0.04, label.y = 0.96) #Sex and gene flow
    fig2 <- ggarrange(fig2.1, fig2.2, nrow=2, align = "v", heights = c(2, 1), labels = c("B", "D"), label.x = 0.04, label.y = 0.96) #Sex and no gene flow
    fig3 <- ggarrange(fig3.1, fig3.2, nrow=2, align = "v", heights = c(2, 1), labels = c("E", "G"), label.x = 0.20, label.y = 0.96) #No sex and gene flow
    fig1 <- ggarrange(fig4.1, fig4.2, nrow=2, align = "v", heights = c(2, 1), labels = c("A", "C"), label.x = 0.20, label.y = 0.96) #No sex and no gene flow
    
    fig.final <- ggarrange(fig1, fig2, fig3, fig4, nrow=2, ncol=2, widths = c(1.2, 1))
    ggsave(fig.final, filename = "4_results/Figures/FigCombined.pdf", device = "pdf", units = "cm", width = 30, height = 24, dpi=250)
    ggsave(fig.final, filename = "4_results/Figures/FigCombined.png", device = "png", units = "cm", width = 30, height = 24, dpi=250)
  }
  
  
  
  
}

#12) Create figure showing perc of dispersers over time
{
  dd.disp <- filter(dd2, Treatment=="EVO")
  ggplot(dd.disp, aes(x = days, y = perc.disp, colour = Gradient)) + geom_point() + facet_grid(Sex~Gene.Flow)
  ggplot(dd.disp, aes(x = days, y = perc.disp, colour = Gradient, group = ID)) + geom_point() + facet_grid(Sex~Gene.Flow) + geom_line()
  ggplot(dd.disp, aes(x = dens.front, y = perc.disp, colour = Gradient)) + geom_point() + facet_grid(Sex~Gene.Flow) + geom_smooth()
}

#13)
{
  summary.disp <- data.all %>% group_by(Gradient, Sex, Gene.Flow) %>% summarise(meandisp = mean(perc.disp, na.rm=T))
  dat.disp.all <- filter(data.all, !is.na(data.all$perc.disp))
  
  ggplot(dat.disp.all, aes(y=perc.disp*100, x=pos, colour=as.character(Gradient) )) + geom_point(size=2) + facet_grid(Gene.Flow~Sex) + geom_line(aes(group=as.factor(ID))) +
    xlab("Space (patches)") +
    ylab("Percentage of dispersers")+ scale_y_continuous() +
    theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5), legend.position = "none")+
    theme(strip.background =element_rect(fill="grey"))+
    theme(strip.text = element_text(colour = 'black')) +
    scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient")) +
    scale_fill_manual(values=c("#e41a1c","#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("Uniform", "Gradient"))
  
  ggsave(filename = "4_results/Figures/DispersalRates", device = "png", width = 9, height = 5.88, dpi = 300)

}