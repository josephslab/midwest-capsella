# Reanalysis of Pi Values for Chrissy's manuscript

# Pi at 4 fold degenerate sites across midwestern cities in 10kb windows

repo_location <- "/mnt/home/wils1582/midwest-capsella/"
# load data
pi <- read.delim(paste0(repo_location, "data/Midwest_4-fold_10kb_pi.txt"))

# quick violin of all data
library(tidyverse)
ggplot() + geom_boxplot(data=pi, 
                        aes(x=pop, y = avg_pi))

# just the means
pi2 <- pi %>%
  summarize(avg_pi_genome = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm=T),
            .by=pop)

# mean plot
ggplot() + geom_point(data=pi2, 
                        aes(x=pop, y = avg_pi_genome, color=pop),
                      size = 5) +
  theme_classic()

###### Adding the permutation for adding confidence intervals to the means
### The following code is left exactly as Chrissy implemented originally.
#Subsetting data sets
Detroit<-subset(pi, pop =="D")

Pi_Det<-sum(Detroit$count_diffs) / sum(Detroit$count_comparisons) #Checking if avg Pi is the same and it is. 
Pi_Det

GR<-subset(pi, pop =="GR")

Pi_GR<-sum(GR$count_diffs) / sum(GR$count_comparisons)

Kal<-subset(pi, pop =="Kal")

Mad<-subset(pi, pop =="Mad")

Mil<-subset(pi, pop =="Mil")
# Bootstrapping to get CI for each population for the mean -----------------------------


# Detroit ----------------------------
n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Detroit$window_pos_1)# Getting unique genimic sites
Detroit$window_pos_1<-as.numeric(Detroit$window_pos_1)

library(dplyr)
Det_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Detroit %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    mean(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Det_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Det_bootstrap_pi<-Det_bootstrap_pi[!is.na(Det_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI) for the mean
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Det_bootstrap_pi) > 0) {
  pi_CI <- quantile(Det_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Mean Pi:", mean(Det_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Mean Pi:", mean(Det_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Det_bootstrap_pi, breaks=30, col="lightblue", xlab= "Mean Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Det_Pi_summary <- data.frame(Detroit %>%
                               group_by(pop) %>%
                               summarise(
                                 mean_pi = mean(avg_pi, na.rm = TRUE),
                                 median_pi = median(avg_pi, na.rm = TRUE),
                                 lower_CI = quantile(Det_bootstrap_pi, 0.025, na.rm = TRUE),
                                 upper_CI = quantile(Det_bootstrap_pi, 0.975, na.rm = TRUE)
                               ))

print(Det_Pi_summary)
#POP MEAN_PI LOWER CI  UPPER CI Median_pi 
#D    0.0966    0.0955    0.0975      0.0711

# Grand Rapids ------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(GR$window_pos_1)# Getting unique genimic sites
GR$window_pos_1<-as.numeric(GR$window_pos_1)

library(dplyr)
GR_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- GR %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    mean(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(GR_bootstrap_pi)
#Remove NA values before calculating confidence intervals
GR_bootstrap_pi<-GR_bootstrap_pi[!is.na(GR_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(GR_bootstrap_pi) > 0) {
  pi_CI <- quantile(GR_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Mean Pi:", mean(GR_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Mean Pi:", mean(GR_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(GR_bootstrap_pi, breaks=30, col="lightblue", xlab= "Mean Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


GR_Pi_summary <- data.frame(GR %>%
                              group_by(pop) %>%
                              summarise(
                                mean_pi = mean(avg_pi, na.rm = TRUE),
                                median_pi = median(avg_pi, na.rm = TRUE),
                                lower_CI = quantile(GR_bootstrap_pi, 0.025, na.rm = TRUE),
                                upper_CI = quantile(GR_bootstrap_pi, 0.975, na.rm = TRUE)
                              ))

print(GR_Pi_summary)
#POP  MEAN PI  LOWER CI  UPPER CI Median_pi 
#GR    0.08282   0.0819   0.0836      0.0612

# Kalamazoo ---------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Kal$window_pos_1)# Getting unique genimic sites
Kal$window_pos_1<-as.numeric(Kal$window_pos_1)

library(dplyr)
Kal_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Kal %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    mean(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Kal_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Kal_bootstrap_pi<-Kal_bootstrap_pi[!is.na(Kal_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Kal_bootstrap_pi) > 0) {
  pi_CI <- quantile(Kal_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Mean Pi:", mean(Kal_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Mean Pi:", mean(Kal_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Kal_bootstrap_pi, breaks=30, col="lightblue", xlab= "Mean Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Kal_Pi_summary <- data.frame(Kal %>%
                               group_by(pop) %>%
                               summarise(
                                 mean_pi = mean(avg_pi, na.rm = TRUE),
                                 median_pi = median(avg_pi, na.rm = TRUE),
                                 lower_CI = quantile(Kal_bootstrap_pi, 0.025, na.rm = TRUE),
                                 upper_CI = quantile(Kal_bootstrap_pi, 0.975, na.rm = TRUE)
                               ))

print(Kal_Pi_summary)

#pop  mean pi lower ci  upper ci
#Kal   0.0995   0.0985    0.1004

# Milwaukee ---------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Mil$window_pos_1)# Getting unique genimic sites
Mil$window_pos_1<-as.numeric(Mil$window_pos_1)

library(dplyr)
Mil_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Mil %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    mean(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Mil_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Mil_bootstrap_pi<-Mil_bootstrap_pi[!is.na(Mil_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Mil_bootstrap_pi) > 0) {
  pi_CI <- quantile(Mil_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Mean Pi Mil:", mean(Mil_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Mean Pi Mil:", mean(Mil_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Mil_bootstrap_pi, breaks=30, col="lightblue", xlab= "Mean Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Mil_Pi_summary <-data.frame(Mil %>%
                              group_by(pop) %>%
                              summarise(
                                mean_pi = mean(avg_pi, na.rm = TRUE),
                                median_pi = median(avg_pi, na.rm = TRUE),
                                lower_CI = quantile(Mil_bootstrap_pi, 0.025, na.rm = TRUE),
                                upper_CI = quantile(Mil_bootstrap_pi, 0.975, na.rm = TRUE)
                              ))

print(Mil_Pi_summary)
#pop   mean pi   lower ci   upper ci  median
#Mil   0.0848     0.0840     0.0857    0.05936

# Madison -----------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Mad$window_pos_1)# Getting unique genimic sites
Mad$window_pos_1<-as.numeric(Mad$window_pos_1)

library(dplyr)
Mad_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Mad %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    mean(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Mad_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Mad_bootstrap_pi<-Mad_bootstrap_pi[!is.na(Mad_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Mad_bootstrap_pi) > 0) {
  pi_CI <- quantile(Mad_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Mean Pi:", mean(Mad_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Mean Pi:", mean(Mad_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Mad_bootstrap_pi, breaks=30, col="lightblue", xlab= "Mean Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Mad_Pi_summary <-  data.frame(Mad %>%
                                group_by(pop) %>%
                                summarise(
                                  mean_pi = mean(avg_pi, na.rm = TRUE),
                                  median_pi = median(avg_pi, na.rm = TRUE),
                                  lower_CI = quantile(Mad_bootstrap_pi, 0.025, na.rm = TRUE),
                                  upper_CI = quantile(Mad_bootstrap_pi, 0.975, na.rm = TRUE)
                                ))

print(Mad_Pi_summary)

# Merging the data together for a plot and analysis -----------------------

Cities_Pi_CI <- bind_rows(Det_Pi_summary, GR_Pi_summary, Kal_Pi_summary, Mil_Pi_summary, Mad_Pi_summary)

head(Cities_Pi_CI)

str(Cities_Pi_CI)
#write.csv(Cities_Pi_CI, "C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Genotypic data/Cities_Pi_CI_20March_MWB.csv")


# Mean plot ---------------------------------------------------------------
#Cities_Pi_CI<-read.csv("C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Genotypic data/Cities_Pi_CI_20March_MWB.csv")


Cities_Pi_CI$pop<-c("Detroit", "Grand Rapids", "Kalamazoo", "Milwaukee", "Madison")
head(Cities_Pi_CI)

Cities_Pi_CI_reorderd<-Cities_Pi_CI %>% arrange(pop)
head(Cities_Pi_CI_reorderd)

Cities_Pi_CI_reorderd$pop<-as.factor(Cities_Pi_CI_reorderd$pop)
#install.packages("ggplot2")
# library(ggplot2)
# windowsFonts(A = windowsFont("Times New Roman"))
#tiff ("Pi Plot_20March_MWBcode.tiff", width = 12, height = 9, units = 'in', res = 80)
MeanCI_plot<-ggplot(Cities_Pi_CI_reorderd, aes(x=pop, y=mean_pi, fill = pop)) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, size = 1)+
  geom_point(aes(shape = pop, colour=pop, fill=pop), size = 8) +
  scale_colour_manual(values = c("black", "black", "black", "black", "black"), 
                      labels = c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee"))+
  scale_fill_manual(values = c("blue", "#06d6a0", "purple","gold","#ff6f61"),  
                    labels = c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee")) + 
  scale_shape_manual(values=c(21,24, 22, 23,25),
                     labels = c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee"))+ 
  theme(
    panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.background = element_rect(colour = NA),
    legend.key = element_rect(colour = "white", fill = NA),
    axis.line=element_line(size = 0.5, colour = "black"),
    axis.text=element_text(colour="black", size = 25, family="A"),
    axis.title=element_text(size=30, family = "A")
  )+
  theme(legend.text = element_text(size = 25, family="A"), #Text was originally 20 and 25 but increased it due to looking inconsistent with other plots
        legend.title = element_text(size = 30, family="A")) +
  labs (y = expression("Mean " * pi), x= "City")
MeanCI_plot
#while (!is.null(dev.list()))  dev.off()

ggsave(plot=MeanCI_plot, 
       filename="~/mean_pi_ci.jpeg",
       height = 10, width = 15, units = "in")

# Bootstrapping to get CI for each population for the median!! -----------------------------


# Detroit ----------------------------
n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Detroit$window_pos_1)# Getting unique genimic sites
Detroit$window_pos_1<-as.numeric(Detroit$window_pos_1)

library(dplyr)
Det_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Detroit %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    median(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Det_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Det_bootstrap_pi<-Det_bootstrap_pi[!is.na(Det_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI) for the median
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Det_bootstrap_pi) > 0) {
  pi_CI <- quantile(Det_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Median Pi:", median(Det_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Median Pi:", median(Det_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Det_bootstrap_pi, breaks=30, col="lightblue", xlab= "Median Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Det_Pi_summary <- data.frame(Detroit %>%
                               group_by(pop) %>%
                               summarise(
                                 mean_pi = mean(avg_pi, na.rm = TRUE),
                                 median_pi = median(avg_pi, na.rm = TRUE),
                                 lower_CI = quantile(Det_bootstrap_pi, 0.025, na.rm = TRUE),
                                 upper_CI = quantile(Det_bootstrap_pi, 0.975, na.rm = TRUE)
                               ))

print(Det_Pi_summary)
#POP MEAN_PI LOWER CI  UPPER CI Median_pi 
#D    0.0966    0.0955    0.0975      0.0711

# Grand Rapids ------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(GR$window_pos_1)# Getting unique genimic sites
GR$window_pos_1<-as.numeric(GR$window_pos_1)

library(dplyr)
GR_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- GR %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    median(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(GR_bootstrap_pi)
#Remove NA values before calculating confidence intervals
GR_bootstrap_pi<-GR_bootstrap_pi[!is.na(GR_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(GR_bootstrap_pi) > 0) {
  pi_CI <- quantile(GR_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Median Pi:", median(GR_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Median Pi:", median(GR_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(GR_bootstrap_pi, breaks=30, col="lightblue", xlab= "Median Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


GR_Pi_summary <- data.frame(GR %>%
                              group_by(pop) %>%
                              summarise(
                                mean_pi = mean(avg_pi, na.rm = TRUE),
                                median_pi = median(avg_pi, na.rm = TRUE),
                                lower_CI = quantile(GR_bootstrap_pi, 0.025, na.rm = TRUE),
                                upper_CI = quantile(GR_bootstrap_pi, 0.975, na.rm = TRUE)
                              ))

print(GR_Pi_summary)
#POP  MEAN PI  LOWER CI  UPPER CI Median_pi 
#GR    0.08282   0.0819   0.0836      0.0612

# Kalamazoo ---------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Kal$window_pos_1)# Getting unique genimic sites
Kal$window_pos_1<-as.numeric(Kal$window_pos_1)

library(dplyr)
Kal_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Kal %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    median(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Kal_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Kal_bootstrap_pi<-Kal_bootstrap_pi[!is.na(Kal_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Kal_bootstrap_pi) > 0) {
  pi_CI <- quantile(Kal_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Median Pi:", median(Kal_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Mean Pi:", median(Kal_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Kal_bootstrap_pi, breaks=30, col="lightblue", xlab= "Median Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Kal_Pi_summary <- data.frame(Kal %>%
                               group_by(pop) %>%
                               summarise(
                                 mean_pi = mean(avg_pi, na.rm = TRUE),
                                 median_pi = median(avg_pi, na.rm = TRUE),
                                 lower_CI = quantile(Kal_bootstrap_pi, 0.025, na.rm = TRUE),
                                 upper_CI = quantile(Kal_bootstrap_pi, 0.975, na.rm = TRUE)
                               ))

print(Kal_Pi_summary)

#pop  mean pi lower ci  upper ci
#Kal   0.0995   0.0985    0.1004

# Milwaukee ---------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Mil$window_pos_1)# Getting unique genimic sites
Mil$window_pos_1<-as.numeric(Mil$window_pos_1)

library(dplyr)
Mil_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Mil %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    median(boot_sample$avg_pi, na.rm=TRUE)#Computing the mean Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Mil_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Mil_bootstrap_pi<-Mil_bootstrap_pi[!is.na(Mil_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Mil_bootstrap_pi) > 0) {
  pi_CI <- quantile(Mil_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Median Pi Mil:", median(Mil_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Median Pi Mil:", median(Mil_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Mil_bootstrap_pi, breaks=30, col="lightblue", xlab= "Median Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Mil_Pi_summary <-data.frame(Mil %>%
                              group_by(pop) %>%
                              summarise(
                                mean_pi = mean(avg_pi, na.rm = TRUE),
                                median_pi = median(avg_pi, na.rm = TRUE),
                                lower_CI = quantile(Mil_bootstrap_pi, 0.025, na.rm = TRUE),
                                upper_CI = quantile(Mil_bootstrap_pi, 0.975, na.rm = TRUE)
                              ))

print(Mil_Pi_summary)
#pop   mean pi   lower ci   upper ci  median
#Mil   0.0848     0.0840     0.0857    0.05936

# Madison -----------------------------------------------------------------

n_boot<-1000 #Defining the number of bootstraps
mysites <- unique(Mad$window_pos_1)# Getting unique genimic sites
Mad$window_pos_1<-as.numeric(Mad$window_pos_1)

library(dplyr)
Mad_bootstrap_pi<-replicate(n_boot, { #bootstrap function
  sampled_sites<-sample(mysites, replace=T, size=length(mysites)) #sample sites with replacement
  boot_sample<- Mad %>% filter(window_pos_1 %in% sampled_sites)#Filter the dataset to include only the sampled sites
  if (nrow(boot_sample)>0){
    median(boot_sample$avg_pi, na.rm=TRUE)#Computing the median Pi for the bootstrap sample
  }else {
    NA #return NA if no valid data
  }
})

#warnings()
head(Mad_bootstrap_pi)
#Remove NA values before calculating confidence intervals
Mad_bootstrap_pi<-Mad_bootstrap_pi[!is.na(Mad_bootstrap_pi)]


#Compute bootstrap confidence intervals (95% CI)
#pi_CI<-quantile(Det_bootstrap_pi, probs = c(0.025, 0.975))
if (length(Mad_bootstrap_pi) > 0) {
  pi_CI <- quantile(Mad_bootstrap_pi, probs = c(0.025, 0.975), na.rm = TRUE)
  cat("Bootstrap Median Pi:", median(Mad_bootstrap_pi), "\n")
  cat("95% CI:", pi_CI[1], "-", pi_CI[2], "\n")
  #Results
  cat("Bootstrap Median Pi:", median(Mad_bootstrap_pi), "\n")
  #
  cat("95% CI:", pi_CI[1], "_", pi_CI[2], "\n")
  
  #plot bootstrap distribution- inside previous function
  hist(Mad_bootstrap_pi, breaks=30, col="lightblue", xlab= "Median Pi")
  abline(v=pi_CI, col="red", lwd=2, lty=2)}


Mad_Pi_summary <-  data.frame(Mad %>%
                                group_by(pop) %>%
                                summarise(
                                  mean_pi = mean(avg_pi, na.rm = TRUE),
                                  median_pi = median(avg_pi, na.rm = TRUE),
                                  lower_CI = quantile(Mad_bootstrap_pi, 0.025, na.rm = TRUE),
                                  upper_CI = quantile(Mad_bootstrap_pi, 0.975, na.rm = TRUE)
                                ))

print(Mad_Pi_summary)

#pop   mean pi   lower ci  upper ci    median
#Mad    0.113     0.11189     0.1142    0.0816



# Merging the data together for a plot and analysis -----------------------

Cities_Pi_CI <- bind_rows(Det_Pi_summary, GR_Pi_summary, Kal_Pi_summary, Mil_Pi_summary, Mad_Pi_summary)

head(Cities_Pi_CI)

str(Cities_Pi_CI)
#write.csv(Cities_Pi_CI, "C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Genotypic data/Cities_Pi_CI_MWB_MedianCI_MAY.csv")


# Creating plot -----------------------------------------------------------



# Median plot -------------------------------------------------------------

Cities_Pi_CI$pop<-c("Detroit", "Grand Rapids", "Kalamazoo", "Milwaukee", "Madison")
head(Cities_Pi_CI)

Cities_Pi_CI_reorderd<-Cities_Pi_CI %>% arrange(pop)
head(Cities_Pi_CI_reorderd)

Cities_Pi_CI_reorderd$pop<-as.factor(Cities_Pi_CI_reorderd$pop)
#Plotting the Median with Mayas windowed analysis file

head(Cities_Pi_CI_reorderd)

#Median plot

#windowsFonts(A = windowsFont("Times New Roman"))
#tiff ("Pi Plot_May_MWBcode_Median_Errorbars.tiff", width = 12, height = 9, units = 'in', res = 80)
Median_plot<-ggplot(Cities_Pi_CI_reorderd, aes(x=pop, y=median_pi, fill = pop)) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, size = 1)+
  geom_point(aes(shape = pop, colour=pop, fill=pop), size = 8) +
  scale_colour_manual(values = c("black", "black", "black", "black", "black"), 
                      labels = c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee"))+
  scale_fill_manual(values = c("blue", "#06d6a0", "purple","gold","#ff6f61"),  
                    labels = c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee")) + 
  scale_shape_manual(values=c(21,24, 22, 23,25),
                     labels = c("Detroit", "Grand Rapids", "Kalamazoo", "Madison", "Milwaukee"))+ 
  theme(
    panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.background = element_rect(colour = NA),
    legend.key = element_rect(colour = "white", fill = NA),
    axis.line=element_line(size = 0.5, colour = "black"),
    axis.text=element_text(colour="black", size = 25, family="A"),
    axis.title=element_text(size=30, family = "A")
  )+
  theme(legend.text = element_text(size = 25, family="A"), #Text was originally 20 and 25 but increased it due to looking inconsistent with other plots
        legend.title = element_text(size = 30, family="A"))+
  labs (y = expression("Median " * pi), x= "City")
Median_plot
#while (!is.null(dev.list()))  dev.off()


# Histograms of data normality --------------------------------------------



library(car)

#is data normally distributed?
hist(pi$avg_pi, probability = TRUE, col = "lightblue", main = "Histogram with Normal Curve")
#no due to zero inflation.

hist(Detroit$avg_pi, probability = TRUE, col = "lightblue", main = "Histogram with Normal Curve")
#no due to zero inflation.

hist(GR$avg_pi, probability = TRUE, col = "lightblue", main = "Histogram with Normal Curve")
#no due to zero inflation.


# Permutation ANOVA  ------------------------------------------------------



# Observed ANOVA F-statistic
pi$pop<-as.factor(pi$pop)

anova_result <- aov(avg_pi ~ pop, data = pi)
summary(anova_result)
obs_F <- summary(anova_result)[[1]]["pop", "F value"]
obs_F
# Permutation test function
perm_test <- function(data, n_perm = 10000) {
  perm_F_values <- replicate(n_perm, {
    shuffled_pop <- sample(data$pop)  # Shuffle population labels
    perm_anova <- aov(avg_pi ~ shuffled_pop, data = data)  # Recalculate ANOVA
    return(summary(perm_anova)[[1]]["shuffled_pop", "F value"])
  })
  
  # Compute p-value
  p_value <- mean(perm_F_values >= obs_F)
  return(p_value)
}

# Run permutation test
set.seed(123)
p_value <- perm_test(pi)

# Print results
cat("Observed ANOVA F-statistic:", round(obs_F, 4), "\n")
#1417.209
#With Mayas windowed data- 340.7789 
cat("Permutation test p-value:", round(p_value, 6), "\n")
#0


# Pairwise Comparison test ------------------------------------------------


# pairwise test

# install.packages("RVAideMemoire")
# citation("RVAideMemoire")
# library(RVAideMemoire)
# 
# # Assume your data is stored in a data frame called `data`
# pairwise_result <- pairwise.perm.t.test(Midwest_pi$avg_pi, Midwest_pi$pop)
# 
# # View the results
# print(pairwise_result)
# 
# #Output shows pvalues for each comparison.
# 
# #D      GR     Kal    Mad   
# #GR  0.0025 -      -      -     
# # Kal 0.0044 0.0025 -      -     
# # Mad 0.0025 0.0025 0.0025 -     
# # Mil 0.0025 0.0160 0.0025 0.0025
# 
# #P value adjustment method: fdr 


