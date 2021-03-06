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
  library(car)
}

#3) Set working directory and read in data
{
  #setwd("F:/Documenten/PhD/08_Gene_Swamping")
  setwd("/media/felix/DataDrive2/Documenten/PhD/08_Gene_Swamping/")
  dd <- read_tsv("2_data/6_Survival_Data.csv")
}

#4) Data handling
{
  #4.1) Create treatment variable
  dd$Treatment <- as.factor(paste(dd$Strain, dd$Gradient, dd$`Gene flow`, dd$sex))
  dd$Strain <- ifelse(dd$Strain=="CU427.4", "clone 1", ifelse(dd$Strain=="CU428.2", "clone 2", 
                                                                            ifelse(dd$Strain=="SB3539", "clone 3", ifelse(dd$Strain=="B2086.2", "clone 4", "mix"))))
 
  #Filter data for evolution part and rename variables properly
  dd.evo <- filter(dd, Strain=="mix")
  dd.evo$sex <- factor(ifelse(dd.evo$sex=="y", "Sexual", "Asexual"), levels = c("Asexual", "Sexual"))
  dd.evo$Gradient <- factor(ifelse(dd.evo$Gradient=="y", "Gradient", "Uniform"), levels = c("Uniform", "Gradient"))
  dd.evo$'Gene flow' <- factor(ifelse(dd.evo$'Gene flow'=="y", "Present", "Absent"), levels = c("Absent", "Present"))
  
  #Make survival class
  dd.evo$Surv <- ifelse(dd.evo$Survival == 1, "yes", "no")
  
  #Summarize to probabilities
  dd.evo$class <- paste(dd.evo$Gradient, dd.evo$'Gene flow', dd.evo$sex)
  dd.evo.sum <- group_by(dd.evo, class)%>% summarize(probSurv = sum(Survival)/5)
  dd.evo.sum$Gradient <-ifelse("Uniform" %in% dd.evo.sum$class, "Uniform", "Gradient")
  dd.evo.sum$Sex <-ifelse("Asexual" %in% dd.evo.sum$class, "Asexual", "Sexual")
  dd.evo.sum$GeneFlow <-ifelse("Absent" %in% dd.evo.sum$class, "Absent", "Present")
  
}


#redo model using Bayesian framework
{
  library(rethinking)
  dd <- dd.evo
  ddstats <- list(surv = dd$Survival, strain = dd$Strain, grad = ifelse(dd$Gradient=="Uniform", 0, 1), geneflow = ifelse(dd$`Gene flow`=="Present", 1, 0), sex = ifelse(dd$sex=="Sexual", 1, 0))

  
  #6.1) Survival ~ 0
  {
    m01 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a ,
        a ~dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m01)
  }
  
  #6.2) Survival ~ Gradient
  {
    m02 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + grad*bgrad,
        a ~dnorm(0, 10),
        bgrad ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m02)
  }
  
  #6.3) Survival ~ geneflow
  {
    m03 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow,
        a ~dnorm(0, 10),
        bgeneflow ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m03)
  }
  
  #6.4) Survival ~ sex
  {
    m04 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + sex*bsex,
        a ~dnorm(0, 1),
        bsex ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m04)
  }
  
  #6.5) Survival ~ gradient + geneflow
  {
    m05 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + grad*bgrad + geneflow*bgeneflow,
        a ~dnorm(0, 1),
        bgrad ~ dnorm(0, 10),
        bgeneflow ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m05)
  }
  
  #6.6) Survival ~ gradient * geneflow
  {
    m06 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + grad*bgrad + geneflow*bgeneflow + bint1*grad*geneflow,
        a ~dnorm(0, 1),
        bgrad ~ dnorm(0, 10),
        bgeneflow ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m06)
  }
  
  #6.7) Survival ~ gradient + sex
  {
    m07 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + grad*bgrad + sex*bsex,
        a ~dnorm(0, 1),
        bgrad ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m07)
  }
  
  #6.8) Survival ~ gradient * sex
  {
    m08 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + grad*bgrad + sex*bsex + bint1*grad*sex,
        a ~dnorm(0, 1),
        bgrad ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m08)
  }
  
  #6.9) Survival ~ geneflow + sex
  {
    m09 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m09)
  }
  
  #6.10) Survival ~ geneflow * sex
  {
    m10 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bint1*geneflow*sex,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m10)
  }
  
  #6.11) Survival ~ gradient + geneflow + sex
  {
    m11 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m11)
  }
  
  #6.12) Survival ~ gradient * geneflow + sex
  {
    m12 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*grad*geneflow,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m12)
  }
  
  #6.13) Survival ~ gradient + geneflow * sex
  {
    m13 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*sex*geneflow,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m13)
  }
  
  #6.14) Survival ~ gradient * sex + geneflow
  {
    m14 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*sex*grad,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m14)
  }
  
  #6.15) Survival ~ gradient * sex * geneflow
  {
    m15 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*sex*geneflow+ bint2*grad*geneflow + bint3*grad*sex + bint4*sex*grad*geneflow,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10),
        bint2 ~ dnorm(0, 10),
        bint3 ~ dnorm(0, 10),
        bint4 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m15)
  }
  
  #6.16) Survival ~ gradient * sex + geneflow*sex
  {
    m16 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*sex*grad+ bint2*sex*geneflow,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10),
        bint2 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m16)
  }
  
  #6.17) Survival ~ gradient * sex + geneflow*gradient
  {
    m17 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*sex*grad+ bint2*grad*geneflow,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10),
        bint2 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m17)
  }
  
  #6.18) Survival ~ gradient * geneflow + geneflow*sex
  {
    m18 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*sex*geneflow+ bint2*grad*geneflow,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10),
        bint2 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m18)
  }
  
  #6.19) Survival ~ gradient * geneflow + geneflow*sex
  {
    m19 <- map2stan(
      alist(
        surv ~ dbinom(1, p),
        logit(p) <- a + geneflow*bgeneflow + sex*bsex + bgrad*grad + bint1*sex*geneflow+ bint2*grad*geneflow + bint3*grad*sex,
        a ~dnorm(0, 1),
        bgeneflow ~ dnorm(0, 10),
        bsex ~ dnorm(0, 10),
        bgrad ~ dnorm(0, 10),
        bint1 ~ dnorm(0, 10),
        bint2 ~ dnorm(0, 10),
        bint3 ~ dnorm(0, 10)
      ) ,data=ddstats ,
      iter=1e4, warmup=2e3,
      control=list(adapt_delta=0.95), WAIC=FALSE )
    precis(m19)
  }
  
  compare(m01, m02, m03, m04, m05, m06, m07, m08, m09, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19)
  precis(m05)
  
  #Predict based on weighted model
  dd.predict2 <- expand.grid(geneflow = unique(ddstats$geneflow), sex = unique(ddstats$sex), grad = unique(ddstats$grad))
  pred <- ensemble(m01, m02, m03, m04, m05, m06, m07, m08, m09, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, data = dd.predict2)
  
  #combine predictins with variables
  dd.predict2$meansurv <- apply(pred$link, 2, mean)
  dd.predict2$uppersurv <- apply(pred$link, 2, PI, prob = 0.95)[2, ]
  dd.predict2$lowersurv <- apply(pred$link, 2, PI, prob = 0.95)[1, ]
  dd.predict2$Gradient <- factor(ifelse(dd.predict2$grad==0, "Uniform", "Gradient"), levels = c("Uniform", "Gradient"))
  dd.predict2$Sex <- factor(ifelse(dd.predict2$sex==1, "Sexual", "Asexual"), levels = c("Asexual", "Sexual"))
  dd.predict2$GeneFlow <- factor(ifelse(dd.predict2$geneflow==1, "Present", "Absent"), levels = c("Absent", "Present"))
  
  #Plot predictions and data
  dd.evo <- mutate(dd.evo, Sex=sex)
  dd.evo$GeneFlow <- dd.evo$`Gene flow`
  ggplot(dd.evo, aes(x=Gradient, y=as.numeric(Survival), colour=Gradient)) + ylab("Survival probability")  + facet_grid(GeneFlow~Sex) +
    geom_boxplot(inherit.aes = F, data = dd.predict2, mapping = aes(x = Gradient, lower = lowersurv, upper = uppersurv, middle = meansurv, ymin = lowersurv, ymax = uppersurv, fill = Gradient, colour=NA), stat = "identity", alpha = 0.3) +
    geom_boxplot(inherit.aes = F, data = dd.predict2, mapping = aes(x = Gradient, lower = meansurv, upper = meansurv, middle = meansurv, ymin = meansurv, ymax = meansurv, fill = NA), stat = "identity", alpha = 0.3) +
    geom_point(position = position_dodge2(width = 0.5), size=3) + 
    theme_light() +   theme(axis.text=element_text(size=12), legend.text=element_text(size=16),legend.title=element_text(size=16), strip.text.x=element_text(20),
                            axis.title=element_text(size=16), strip.text = element_text(size=16), axis.title.x = element_blank(),
                            axis.text.x=element_text(size=16, color = "black"), legend.position="none")+
    theme(strip.background =element_rect(fill="grey"))+
    theme(strip.text = element_text(colour = 'black')) +
    scale_color_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("Uniform", "Gradient"), name="Abiotic conditions", labels=c("Uniform", "Gradient")) +
    scale_fill_manual(values=c("#e41a1c", "#2d33f9"), breaks=c("Uniform", "Gradient"), name="Abiotic conditions", labels=c("Uniform", "Gradient"))

  ggsave(filename = "4_results/Figures/05_ExtinctionProbability.png", device = "png", width = 8.5, height = 5.88, dpi = 300)

}
