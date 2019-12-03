#Clear memory
rm(list=ls())

#load libraries
library(ggplot2)
library(tidyverse)
library(scales)
library(nlme)
library(MuMIn)

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

}

#Redo figure with plotted model predictions
#Define prediction data
prd.data <- expand.grid(pos = seq(min(data.all$pos), max(data.all$pos), by = 1), Gradient = unique(data.all$Gradient))
predictions <- predict(newdata = prd.data, best.model, se.fit = T, level = 0)
prd.data$mean <- predictions$fit
prd.data$upper <- predictions$fit + 1.96*predictions$se.fit
prd.data$lower <- predictions$fit - 1.96*predictions$se.fit

#Colour Sex events black
data.all$SexEvent <- ifelse(data.all$date %in% c("10.10.18", "24.10.18", "07.11.18", "21.11.18"), "Centrifugation event", "No centrifugation event")
data.all$colour = ifelse(data.all$SexEvent=="Centrifugation event", "Centrifugation event", as.character(data.all$Gradient))


ggplot(data.all, aes(y=dens.front, x=pos, colour=as.character(Gradient) )) + geom_point(size=2) + facet_grid(Sex~Gene.Flow) + geom_line(aes(group=as.factor(ID))) + 
  xlab("Space (patches)") + geom_ribbon(inherit.aes=F, data = prd.data, mapping = aes(x = pos, ymin = lower, ymax = upper, fill=as.character(Gradient)), alpha = 0.3) +
  geom_line(inherit.aes=F, data = prd.data, mapping = aes(x = pos, y = mean, colour=Gradient), size=2) + 
  ylab("Population density, N (indiv/mL)")+ scale_y_continuous(labels = scientific) + 
  theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5))+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black')) +
  scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)")) +
  scale_fill_manual(values=c("#e41a1c","#2d33f9"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)"))
ggsave(filename = "4_results/Figures/02_DensityDuringEvolutionB", device = "png", width = 9, height = 5.88, dpi = 300)

