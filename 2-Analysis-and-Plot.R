setwd("~/Unckless lab/ExperimentalEvolution/Pseud.spp/")

summarize <- function(kin_file, concentration, treatment, class) {
  
  library(tidyverse)
  library(growthcurver)
  library(reshape2)
  
  dod.func <- function(well) {
    dOD <- c()
    for (i in c(1:length(well))) {
      dOD[i] <- well[i+1] - well[i]
    }
    return(dOD)
  }
  maxv.func <- function(dOD, int) {
    
    unsort.dod <- dod.func(dOD) %>% na.omit() %>% as.vector()
    sort.dOD <- unsort.dod %>% sort()
    
    maxVs <- sort.dOD[(length(sort.dOD)-4):(length(sort.dOD))]
    outlier <- max(abs(mean(maxVs) - maxVs))
    outlier.index <- which(abs(mean(maxVs) - maxVs) == outlier)
    
    maxV <- maxVs[-outlier.index] %>% mean()
    tmaxv <- which(unsort.dod %in% maxVs[-outlier.index]) %>% mean()
    tmaxv <- tmaxv/(60/int)
    
    out <- c(maxV, tmaxv)
    
    return(out)
  }
  lt.func <- function(well, int){
    blanked <- well - well[1]
    lagtime <- which(blanked >= 0.01) %>% min() / (60/int)
    return(lagtime)
  }
  
  #vector of all the inner wells that have growing cultures. Outer wells filled with water.
  wells <- c("Time", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11",
             "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11",
             "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11",
             "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11",
             "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11",
             "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11")
  x <- read.delim(file = kin_file, sep = ",", header = T) %>% #read file
    select(all_of(wells)) %>% #select relevant wells
    as.data.frame() #convert to dataframe
  
  kin <- SummarizeGrowthByPlate(x) %>% #calculate growth characteristics
    select(all_of(c("sample", "r", "k", "auc_e", "t_gen"))) #select relevant variables
  
  dOD <- x[nrow(x),-1] - x[1,-1] #calculate difference between first and last read
  dOD <- t(dOD) #transpose
  colnames(dOD) <- c("deltaOD") #rename
  
  wells <- wells[-1]
  maxVs <- matrix(ncol=2, nrow = length(wells))
  rownames(maxVs) <- wells
  lagtime <- c()
  for(well in wells) {
    maxVs[well,] <- maxv.func(x[,well], int = 30)
    lagtime[well] <-  lt.func(x[,well], int = 30)
    if (lagtime[well] == Inf) {
      lagtime[well] <-  NA
    }
  }
  calc.vals <- cbind(maxVs, lagtime)
  colnames(calc.vals) <- c("maxV", "t_maxV", "lagtime")
  
  
  conc <- read.delim(file = concentration, header = F, sep = "\t") %>% 
    t() %>% melt()
  conc <- conc$value
  
  treatment <- read.delim(file = treatment, header = F, sep = "\t") %>% 
    t() %>% melt()
  treatment <- treatment$value
  
  class <- read.delim(file = class, header = F, sep = "\t") %>% 
    t() %>% melt()
  class <- class$value
  
  
  out <- cbind(class, treatment, conc, kin, dOD, calc.vals) #combine all varialbles into a single dataframe
  
  return(out)
}
pa01 <- summarize(kin_file = "PA01_MIC_ABDOPT.csv",
                       concentration = "PA01_MIC_ABDOPT_conc.txt",
                       treatment = "PA01_MIC_ABDOPT_treat.txt",
                       class = "PA01_MIC_ABDOPT_class.txt")
ggplot(data = pa01, aes(x = conc, y = auc_e, color = class)) +
  geom_point() +
  geom_smooth(se=F) +
  facet_wrap(~factor(treatment, levels=c("API", "BAC", "DRO", "ONC", "PYR", "TUR"))) +
  ggtitle("PA01 MIC") +
  theme_bw() + 
  xlab(expression(paste("Conc (", mu,"g/ml)"))) + ylab("Area Under the Curve") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")

summarize2 <- function(kin_file, concentration, treatment, class) {
  
  library(tidyverse)
  library(growthcurver)
  library(reshape2)
  
  dod.func <- function(well) {
    dOD <- c()
    for (i in c(1:length(well))) {
      dOD[i] <- well[i+1] - well[i]
    }
    return(dOD)
  }
  maxv.func <- function(dOD, int) {
    
    unsort.dod <- dod.func(dOD) %>% na.omit() %>% as.vector()
    sort.dOD <- unsort.dod %>% sort()
    
    maxVs <- sort.dOD[(length(sort.dOD)-4):(length(sort.dOD))]
    outlier <- max(abs(mean(maxVs) - maxVs))
    outlier.index <- which(abs(mean(maxVs) - maxVs) == outlier)
    
    maxV <- maxVs[-outlier.index] %>% mean()
    tmaxv <- which(unsort.dod %in% maxVs[-outlier.index]) %>% mean()
    tmaxv <- tmaxv/(60/int)
    
    out <- c(maxV, tmaxv)
    
    return(out)
  }
  lt.func <- function(well, int){
    blanked <- well - well[1]
    lagtime <- which(blanked >= 0.01) %>% min() / (60/int)
    return(lagtime)
  }
  
  #vector of all the inner wells that have growing cultures. Outer wells filled with water.
  wells <- c("Time", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9",
                     "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", 
                     "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", 
                     "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", 
                     "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9",
                     "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9")
  x <- read.delim(file = kin_file, sep = ",", header = T) %>% #read file
    select(all_of(wells)) %>% #select relevant wells
    as.data.frame() #convert to dataframe
  
  kin <- SummarizeGrowthByPlate(x) %>% #calculate growth characteristics
    select(all_of(c("sample", "r", "k", "auc_e", "t_gen"))) #select relevant variables
  
  dOD <- x[nrow(x),-1] - x[1,-1] #calculate difference between first and last read
  dOD <- t(dOD) #transpose
  colnames(dOD) <- c("deltaOD") #rename
  
  wells <- wells[-1]
  maxVs <- matrix(ncol=2, nrow = length(wells))
  rownames(maxVs) <- wells
  lagtime <- c()
  for(well in wells) {
    maxVs[well,] <- maxv.func(x[,well], int = 30)
    lagtime[well] <-  lt.func(x[,well], int = 30)
    if (lagtime[well] == Inf) {
      lagtime[well] <-  NA
    }
  }
  calc.vals <- cbind(maxVs, lagtime)
  colnames(calc.vals) <- c("maxV", "t_maxV", "lagtime")
  
  
  conc <- read.delim(file = concentration, header = F, sep = "\t") %>% 
    t() %>% melt()
  conc <- conc$value
  
  treatment <- read.delim(file = treatment, header = F, sep = "\t") %>% 
    t() %>% melt()
  treatment <- treatment$value
  
  class <- read.delim(file = class, header = F, sep = "\t") %>% 
    t() %>% melt()
  class <- class$value
  
  
  out <- cbind(class, treatment, conc, kin, dOD, calc.vals) #combine all varialbles into a single dataframe
  
  return(out)
}
Pspp <- summarize2(kin_file = "Pseud-spp_MIC_BOT_1.csv",
                      concentration = "Pseud-spp_MIC_BOT_1_conc.txt",
                      treatment = "Pseud-spp_MIC_BOT_1_treat.txt",
                      class = "Pseud-spp_MIC_BOT_1_class.txt")

ggplot(Pspp, aes(x=as.factor(conc), y=auc_e, fill = treatment)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
  facet_wrap(~factor(class, 
                     levels=c("P.aeruginosa(PA01)", "P.putida(KL2480?)", 
                              "P.fluorescens(Dmel)","P.fluorescens(Pf0-1)")))+
  ggtitle("Pseudomonas spp. AMP Susceptibility") +
  theme_classic() + 
  xlab(expression(paste("Conc (", mu,"g/ml)"))) + ylab("Area Under the Curve") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


P.put.flu1 <- summarize(kin_file = "P_put_flu_MIC_BOT_1.csv",
                   concentration = "P_put_flu_MIC_BOT_1_conc.txt",
                   treatment = "P_put_flu_MIC_BOT_1_treat.txt",
                   class = "P_put_flu_MIC_BOT_1_class.txt")
ggplot(P.put.flu1, aes(x=conc, y=auc_e)) + 
  geom_point()+
  geom_line(se=F)+
  facet_wrap(~class*treatment)+
  ggtitle("Pseudomonas putida + fluorescens") +
  theme_classic() + 
  xlab(expression(paste("Conc (", mu,"g/ml)"))) + ylab("Area Under the Curve") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

P.put.flu2 <- summarize(kin_file = "P_put_flu_MIC_BOT_2.csv",
                        concentration = "P_put_flu_MIC_BOT_2_conc.txt",
                        treatment = "P_put_flu_MIC_BOT_2_treat.txt",
                        class = "P_put_flu_MIC_BOT_2_class.txt")

ggplot(P.put.flu2, aes(x=conc, y=auc_e)) + 
  geom_point()+
  geom_line(se=F)+
  facet_wrap(~class*treatment, scales = "free_x")+
  ggtitle("Pseudomonas putida + fluorescens") +
  theme_classic() + 
  xlab(expression(paste("Conc (", mu,"g/ml)"))) + ylab("Area Under the Curve") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(rbind(P.put.flu1, P.put.flu2), aes(x=conc, y=auc_e)) + 
  geom_point()+
  geom_line()+
  facet_wrap(~class*treatment, scales = "free_x")+
  ggtitle("Pseudomonas putida + fluorescens") +
  theme_classic() + 
  xlab(expression(paste("Conc (", mu,"g/ml)"))) + ylab("Area Under the Curve") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

P.put.flu.abx <- summarize(kin_file = "P-put.flu_MIC_GKSpStT_1.csv",
                        concentration = "P-put.flu_MIC_GKSpStT_1_conc.txt",
                        treatment = "P-put.flu_MIC_GKSpStT_1_treat.txt",
                        class = "P-put.flu_MIC_GKSpStT_1_class.txt")


ggplot(P.put.flu.abx, aes(x=conc, y=auc_e, color=class)) + 
  geom_point()+
  geom_line()+
  facet_wrap(~treatment, scales="free_x")+
  ggtitle("Pseudomonas putida + fluorescens") +
  theme_classic() + 
  xlab(expression(paste("Conc (", mu,"g/ml)"))) + ylab("Area Under the Curve") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")  +
  labs(color=NULL)

P.put.mic <- summarize(kin_file = "P_put_MIC_BOTGKT.csv",
                           concentration = "P_put_MIC_BOTGKT_conc.txt",
                           treatment = "P_put_MIC_BOTGKT_treat.txt",
                           class = "P_put_MIC_BOTGKT_class.txt")


ggplot(P.put.mic, aes(x=conc, y=auc_e, color=class)) + 
  geom_point()+
  geom_line()+
  facet_wrap(~factor(treatment, levels = c("BAC", "ONC", "TUR", "GEN", "KAN", "TOB")),
             scales="free_x")+
  ggtitle("Pseudomonas putida MICs") +
  theme_classic() + 
  xlab(expression(paste("Conc (", mu,"g/ml)"))) + ylab("Area Under the Curve") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")  +
  labs(color=NULL)



inf <- read.delim(file = "P_put_flu_Dmel_infections_data.txt", header = T, sep = "\t")
inf$rep=as.factor(inf$rep)
inf$perc.live=100*(inf$live/10)



ggplot(inf, aes(x=day, y = live))  +
  geom_point()+
  geom_line(aes(color=rep)) +
  facet_wrap(~bacteria)+
  theme_classic()
ggplot(inf, aes(x=day, y = perc.live))  +
  geom_point()+
  geom_line(aes(color=rep)) +
  facet_wrap(~bacteria)+
  theme_classic()
ggplot(na.omit(inf), aes(x=day, y = mean.live))  +
  geom_point()+
  geom_line() +
  facet_wrap(~bacteria)+
  theme_classic()














