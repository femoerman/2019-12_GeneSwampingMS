#1) Clear memory
{
  rm(list=ls())
}

#2) Read packages 
{
  library(readr)
  library(tidyverse)
  library(rethinking)
  library(MuMIn)
}

#3) Set working directory and read in data
{
  # setwd("E:/PhD/08_Big_Experiment_1")
  setwd("/media/felix/DataDrive2/Documenten/PhD/08_Gene_Swamping/")
  dd <- read_tsv("2_data/0_b_PositionDuringExperiment.txt")
  dd2 <- read_tsv("2_data/6_Survival_Data.csv")
}

#4) Data handling
{
  #4.1) Create dataframe that is easy to handle
  days <- as.integer(as.Date(colnames(dd)[2:24], "%m/%d/%y"))
  days <- days - days[1]
  colnames(dd) <- c("Number", 1:23)
  
  #4.2) Get all positions in one row
  positions <- data.frame(row.names = 1:(104*23))
  positions$Culture <- rep(1:104, 23)
  positions$Day <- NA
  positions$pos <- NA
  for (i in 1:23){
    start <- (i-1)*104+1
    positions[start:(start+104-1), ]$Day <- rep(days[i], 104)
    positions[start:(start+104-1), ]$pos <- pull(dd[, i+1])
  }
  
  #5.3) copy other variables
  positions$Sex <- dd2$sex
  positions$Strain <- dd2$Strain
  positions$Gradient <- dd2$Gradient
  positions$GeneFlow <- dd2$`Gene flow`
  positions$Sex <- ifelse(positions$Sex == "y", "y", "n")
  
  #5.4) Create treatment variable
  positions$Treatment <- as.factor(paste(positions$Strain, positions$Gradient, positions$`Gene flow`, positions$Sex))
  positions$Strain <- ifelse(positions$Strain=="CU427.4", "clone 1", ifelse(positions$Strain=="CU428.2", "clone 2", 
                      ifelse(positions$Strain=="SB3539", "clone 3", ifelse(positions$Strain=="B2086.2", "clone 4", "mix"))))
}

#5) Plot raw expanxion rate data
#Test effect expansion rate evolved populations
dd <- filter(positions, Strain=="mix", Day==max(positions$Day)) %>% filter(Culture%in% filter(dd2, Survival==1)$Number)
dd$Sex <- factor(ifelse(dd$Sex=="y", "Sex (Sex)", "No sex (Asex)"), levels = c("No sex (Asex)", "Sex (Sex)"))
dd$Gradient <- factor(ifelse(dd$Gradient=="y", "Gradient", "No gradient"), levels = c("No gradient", "Gradient"))
dd$GeneFlow <- factor(ifelse(dd$GeneFlow=="y", "Gene flow (GF)", "No gene flow (NoGF)"), levels = c("No gene flow (NoGF)", "Gene flow (GF)"))

#Fit statstical model
fullmodel <- lm(pos~Gradient*GeneFlow*Sex, dd, na.action = "na.fail")
t <- dredge(fullmodel)
bestmodel <- lm(pos ~Gradient, dd)
anova(bestmodel)
summary(bestmodel)

#Create data showing prediction
prd.data <- expand.grid(Gradient = unique(dd$Gradient))
predictions <- predict(bestmodel, newdata = prd.data, se.fit = T)
prd.data$mean <- predictions$fit
prd.data$upper <- predictions$fit + 1.96*predictions$se.fit
prd.data$lower <- predictions$fit - 1.96*predictions$se.fit

#Plot data + model prediction

ggplot(dd, aes(x=Gradient, y=pos, colour=Gradient)) + xlab("Gradient treatment") + ylab("Expansion distance") + facet_grid(Sex~GeneFlow) +
  geom_boxplot(inherit.aes = F, data = prd.data, mapping = aes(x = Gradient, lower = lower, upper = upper, middle = mean, ymin = lower, ymax = upper, fill = Gradient, colour=NA), stat = "identity", alpha = 0.3) +
  geom_boxplot(inherit.aes = F, data = prd.data, mapping = aes(x = Gradient, lower = mean, upper = mean, middle = mean, ymin = mean, ymax = mean, fill = NA), stat = "identity", alpha = 0.3) +
  geom_point(position = position_dodge2(width = 0.5), size=3) + 
  theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                          axis.title=element_text(size=16), strip.text = element_text(size=16), plot.title = element_text(size=16, hjust = 0.5))+
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'black')) +
  scale_color_manual(values=c("#2d33f9", "#e41a1c"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)")) +
  scale_fill_manual(values=c("#2d33f9","#e41a1c"), breaks=c("No gradient", "Gradient"), name="Evolution treatment", labels=c("No gradient (Const)", "Gradient (Grad)"))
ggsave(filename = "4_results/Figures/01_expansionRate", device = "png", width = 8.5, height = 5.88, dpi = 300)
