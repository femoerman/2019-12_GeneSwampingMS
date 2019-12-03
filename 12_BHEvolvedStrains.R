#1) Clear memory
rm(list=ls())

#2) Load packages
library(PBPGM)
library(tidyverse)

#3) Check prior distributions
hist(exp(rnorm(1e4, log(10e4), 0.5)), breaks=50)    #K-prior
hist(exp(rnorm(1e4, -2.3, 0.5)), breaks=50)         #r0-prior
hist(exp(rnorm(1e4, -2.3, 1.5)), breaks=10)         #d-prior
hist(exp(rnorm(1e4, log(1e3), 0.5)), breaks=50)     #N0-prior

#4) Set working directory and load data evolved strains
setwd("/media/felix/DataDrive2/Documenten/PhD/08_Gene_Swamping")
load("2_data/5_EvolvedStrains/1_RawDensityData-Evolvedlines.RData")

#5) Prepare data for growth curve fitting
fitdata <- select(dd, indiv_per_ml, hours, curveID)
colnames(fitdata) <- c("popsize", "time", "ident")

#6) Fit Beverton-Holt model
fulloutput <- BevertonHolt(fitdata, log(1e5), 0.5, -2.3, 0.5, -2.3, 1.5, log(1e3),0.5, graphname="2_data/6_EvolvedStrainsPosteriors/DataFit", output="both")
sumoutput <- fulloutput$sumoutput
posteriorsevolved <- fulloutput$fulloutput

#7) Save output files
save(sumoutput, file="2_data/6_EvolvedStrainsPosteriors/SummarisedPosteriorsEvolved.RData")
save(posteriorsevolved, file="2_data/6_EvolvedStrainsPosteriors/RawPosteriorsEvolved.RData")
