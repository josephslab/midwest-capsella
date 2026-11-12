#Creating sumarise dataset from final data used in rest of paper.
#Then perge with PCA's and catagorise, and then conduct phenotypic data analysis

Eigenvecs<-read.csv("C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Genotypic data/PCA Eigenvectors_March.csv")
head(Eigenvecs)
Phenotypes<-read.csv ("C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Phenotypic data to analyse/FinalDataExclOutliers_Nov27.csv")

Phenotypes<- Phenotypes[Phenotypes$GrowingConditions == 'GC', ]
Phenotypes$FlowerTime<-Phenotypes$DaysPlantToBolt + Phenotypes$DaysBoltToFlower

Phenotypes$DateCulled<-NULL
Phenotypes$FiveFruitCount<-NULL
Phenotypes$Date_Flower<-NULL
Phenotypes$Date_Bolt<-NULL
Phenotypes$DateTransplanted<-NULL
Phenotypes$Date_1stGerm<-NULL
Phenotypes$DateOutFridge<-NULL
Phenotypes$DatePlatedandFridged<-NULL
Phenotypes$Sequenced<-NULL
Phenotypes$Long<-NULL
Phenotypes$Lat<-NULL
Phenotypes$Unique.ID<-NULL
Phenotypes$F2_Plant_ID <-NULL
Phenotypes$Notes <-NULL

#install.packages("plyr")
library(plyr)

#install.packages("ggplot2")
library(ggplot2)



SumarisedPD  <- ddply(Phenotypes, c("Genotype", "City","ParentGrowingConditions","Population","Individual","GrowingConditions"), summarise,
                      N    = length(F2_Rep),
                      Rosette_D = round(mean(Rosette_D, na.rm=TRUE),2),
                      Lobed_level = round(mean(Lobed_level, na.rm=TRUE),2),
                      Height_Flower = round(mean(Height_Flower, na.rm=TRUE),2),
                      Days_Fower_Cull = round(mean(Days_Fower_Cull, na.rm=TRUE),2),
                      Height_Harvest = round(mean(Height_Harvest, na.rm=TRUE),2),
                      N_Fruits = round(mean(N_Fruits, na.rm=TRUE),2),
                      DaysGermToBolt = round(mean(DaysGermToBolt, na.rm=TRUE),2),
                      Rosette_branches = round(mean(Rosette_branches, na.rm=TRUE),2),
                      Cauline_branches = round(mean(Cauline_branches, na.rm=TRUE),2),
                      DaysBoltToFlower = round(mean(DaysBoltToFlower, na.rm=TRUE),2),
                      DaysPlantToBolt = round(mean(DaysPlantToBolt, na.rm=TRUE),2),
                      Seed_count = round(mean(Seed_count, na.rm=TRUE),2),
                      Fecundity = round(mean(Fecundity, na.rm=TRUE),2),
                      BushyScore = round(mean(BushyScore, na.rm=TRUE),2),
                      FlowerTime = round(mean(FlowerTime, na.rm=TRUE),2))


head(SumarisedPD)

#write.csv(SumarisedPD, "C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Phenotypic data to analyse/SumDataForPCA_May.csv", row.names=FALSE)


PCA_Phenotypes <- merge(SumarisedPD, Eigenvecs, by = "Genotype")
#View(PCA_Phenotypes)

#write.csv(PCA_Phenotypes, "C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Phenotypic data to analyse/PCA_Sum.csv", row.names=FALSE)

#Labeled column which loads on PC1 and PC2 and excluded outlier

#Includes flowering time
PCA_Sum_CatPCs<-read.csv ("C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Phenotypic data to analyse/PCA_Sum_CatPCs.csv")
head(PCA_Sum_CatPCs)

str(PCA_Sum_CatPCs)
PCA_Sum_CatPCs$PCsCAT<-as.factor(PCA_Sum_CatPCs$PCsCAT)
PCA_Sum_CatPCs$Genotype <-as.factor(PCA_Sum_CatPCs$Genotype)

library(ggplot2)
#Fecundity
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Fecundity.png", width = 20, height = 6, units = 'in', res = 300)
Fecundity <- ggplot(data = PCA_Sum_CatPCs, aes(y = Fecundity, x = Genotype, colour = PCsCAT, shape = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  scale_colour_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(1, 16))+
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Fecundity", x = "Genotype") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
Fecundity
#dev.off()


#Rosette_D
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Rosette_D.png", width = 20, height = 6, units = 'in', res = 300)
Rosette_D <- ggplot(data = PCA_Sum_CatPCs, aes(y = Rosette_D, x = Genotype, colour = PCsCAT, shape = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  scale_colour_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(1, 16))+
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Rosette diameter (cm)", x = "Genotype") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
Rosette_D
#dev.off()


#Height_Flower
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Height_Flower.png", width = 20, height = 6, units = 'in', res = 300)
Height_Flower <- ggplot(data = PCA_Sum_CatPCs, aes(y = Height_Flower, x = Genotype, colour = PCsCAT, shape = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  scale_colour_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(1, 16))+
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Height at flowering (cm)", x = "Genotype") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
Height_Flower
#dev.off()



#Height_Harvest
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Height_Harvest.png", width = 20, height = 6, units = 'in', res = 300)
Height_Harvest <- ggplot(data = PCA_Sum_CatPCs, aes(y = Height_Harvest, x = Genotype, colour = PCsCAT, shape = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  scale_colour_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(1, 16))+
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Height at harvest (cm)", x = "Genotype") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
Height_Harvest
#dev.off()


#DaysPlantToBolt
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs DaysPlantToBolt.png", width = 20, height = 6, units = 'in', res = 300)
DaysPlantToBolt <- ggplot(data = PCA_Sum_CatPCs, aes(y = DaysPlantToBolt, x = Genotype, colour = PCsCAT, shape = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  scale_colour_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(1, 16))+
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Bolting time (days)", x = "Genotype") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
DaysPlantToBolt
#dev.off()

#Plot points in the background


#Bushiness
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs BushyScore.png", width = 20, height = 6, units = 'in', res = 300)
BushyScore <- ggplot(data = PCA_Sum_CatPCs, aes(y = BushyScore, x = Genotype, colour = PCsCAT, shape = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  scale_colour_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(1, 16))+
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Bushiness score", x = "Genotype") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
BushyScore
#dev.off()

#FlowerTime
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs FlowerTime.png", width = 20, height = 6, units = 'in', res = 300)
FlowerTime <- ggplot(data = PCA_Sum_CatPCs, aes(y = FlowerTime, x = Genotype, colour = PCsCAT, shape = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 6) +
  scale_colour_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(1, 16))+
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Flowering time (days)", x = "Genotype") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
FlowerTime
#dev.off()


# GGarrange figures -------------------------------------------------------

#install.packages("gridExtra") 
library(gridExtra)

combined_PCplots_supp<-grid.arrange(DaysPlantToBolt, Rosette_D, FlowerTime, Height_Flower, Height_Harvest, BushyScore, Fecundity, ncol = 1, padding = unit(2, "cm"))

width <- 21 # in inches 
height <-30 # in inches 
dpi <- 600 # resolution in dots per inch 

ggsave("combined_PCPLOTS_Easiertosee.tiff", combined_PCplots_supp, width = width, height = height, dpi = dpi)


# Point Plots of PCs rather than genotypes --------------------------------

cluster_n <- sapply(
  split(PCA_Sum_CatPCs$Genotype, PCA_Sum_CatPCs$PCsCAT),
  function(x) length(unique(x))
)

library(ggplot2)
#Fecundity
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Fecundity.png", width = 20, height = 6, units = 'in', res = 300)
FecundityP <- ggplot(data = PCA_Sum_CatPCs, aes(y = Fecundity, x = PCsCAT, colour = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size=1.2, colour="black") +
  stat_summary(fun = mean, geom = "point", size = 7) +
  scale_colour_manual(values = c("purple", "darkblue")) +
  scale_x_discrete(labels = c("PC1" = "Cluster 1",
                              "PC2" = "Cluster 2")) +
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Fecundity", x = "Cluster") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
FecundityP

#dev.off()


#Rosette_D
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Rosette_D.png", width = 20, height = 6, units = 'in', res = 300)
Rosette_DP <- ggplot(data = PCA_Sum_CatPCs, aes(y = Rosette_D, x = PCsCAT, colour = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size=1.2, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 7) +
  scale_colour_manual(values = c("purple", "darkblue")) +
  scale_x_discrete(labels = c("PC1" = "Cluster 1",
                              "PC2" = "Cluster 2")) +
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Rosette diameter (cm)", x = "Cluster") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
Rosette_DP
#dev.off()


#Height_Flower
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Height_Flower.png", width = 20, height = 6, units = 'in', res = 300)
Height_FlowerP <- ggplot(data = PCA_Sum_CatPCs, aes(y = Height_Flower, x = PCsCAT, colour = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size=1.2, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 7) +
  scale_colour_manual(values = c("purple", "darkblue")) +
  scale_x_discrete(labels = c("PC1" = "Cluster 1",
                              "PC2" = "Cluster 2")) +
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Height at flowering (cm)", x = "Cluster") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
Height_FlowerP
#dev.off()



#Height_Harvest
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs Height_Harvest.png", width = 20, height = 6, units = 'in', res = 300)
Height_HarvestP <- ggplot(data = PCA_Sum_CatPCs, aes(y = Height_Harvest, x = PCsCAT, colour = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size=1.2, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 7) +
  scale_colour_manual(values = c("purple", "darkblue")) +
  scale_x_discrete(labels = c("PC1" = "Cluster 1",
                              "PC2" = "Cluster 2")) +
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Height at harvest (cm)", x = "Cluster") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
Height_HarvestP
#dev.off()


#DaysPlantToBolt
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs DaysPlantToBolt.png", width = 20, height = 6, units = 'in', res = 300)
DaysPlantToBoltP <- ggplot(data = PCA_Sum_CatPCs, aes(y = DaysPlantToBolt, x = PCsCAT, colour = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size=1.2, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 7) +
  scale_colour_manual(values = c("purple", "darkblue")) +
  scale_x_discrete(labels = c("PC1" = "Cluster 1",
                              "PC2" = "Cluster 2")) +
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Bolting time (days)", x = "Cluster") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
DaysPlantToBoltP
#dev.off()

#Plot points in the background


#Bushiness
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs BushyScore.png", width = 20, height = 6, units = 'in', res = 300)
BushyScoreP <- ggplot(data = PCA_Sum_CatPCs, aes(y = BushyScore, x = PCsCAT, colour = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size=1.2, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 7) +
  scale_colour_manual(values = c("purple", "darkblue")) +
  scale_x_discrete(labels = c("PC1" = "Cluster 1",
                              "PC2" = "Cluster 2")) +
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Bushiness score", x = "Cluster") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
BushyScoreP
#dev.off()

#FlowerTime
windowsFonts(A = windowsFont("Times New Roman"))
#png ("PCs FlowerTime.png", width = 20, height = 6, units = 'in', res = 300)
FlowerTimeP <- ggplot(data = PCA_Sum_CatPCs, aes(y = FlowerTime, x = PCsCAT, colour = PCsCAT)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size=1.2, color = "black") +
  stat_summary(fun = mean, geom = "point", size = 7) +
  scale_colour_manual(values = c("purple", "darkblue")) +
  scale_x_discrete(labels = c("PC1" = "Cluster 1",
                              "PC2" = "Cluster 2")) +
  theme(axis.line = element_line(size = 0.9, colour = "black", linetype = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 20, family = "A"),  # Changes only x-axis text angle
        axis.text.y = element_text(color = "black", size = 20, family = "A"),  # Keeps y-axis text unchanged
        axis.title = element_text(size = 30, face = "plain"), 
        axis.title.x = element_text(vjust = -3), 
        axis.title.y = element_text(vjust = 3),
        panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA), 
        text = element_text(colour = "black", size = 20, family = "A"), 
        plot.margin = margin(1, 1, 1.5, 1.2, "cm")) +
  labs(y = "Flowering time (days)", x = "Cluster") +
  theme(legend.position = "top right", legend.justification = 'right', legend.text = NULL)
FlowerTimeP
#dev.off()


# GGarrange figures point plots -------------------------------------------------------

#install.packages("gridExtra") 
library(gridExtra)

combined_Point_PCplots_supp<-grid.arrange(DaysPlantToBoltP, Rosette_DP, FlowerTimeP, Height_FlowerP, Height_HarvestP, BushyScoreP, FecundityP, ncol = 2, padding = unit(2, "cm"))

width <- 21 # in inches 
height <-30 # in inches 
dpi <- 600 # resolution in dots per inch 

ggsave("combined_Point_PCPLOTS_Cluster26.tiff", combined_Point_PCplots_supp, width = width, height = height, dpi = dpi)



# Checking normality for analysis -----------------------------------------

#Rosette_D
qqnorm(PCA_Sum_CatPCs$Rosette_D, main='Normal')
qqline(PCA_Sum_CatPCs$Rosette_D)
#Not that normal

#BushyScore
qqnorm(PCA_Sum_CatPCs$BushyScore, main='Normal')
qqline(PCA_Sum_CatPCs$BushyScore)
#Not bad but not great


#DaysPlantToBolt
qqnorm(PCA_Sum_CatPCs$DaysPlantToBolt, main='Normal')
qqline(PCA_Sum_CatPCs$DaysPlantToBolt)
#Terrible

#Height_Harvest
qqnorm(PCA_Sum_CatPCs$Height_Harvest, main='Normal')
qqline(PCA_Sum_CatPCs$Height_Harvest)
#Not really that normal

#Height_Flower
qqnorm(PCA_Sum_CatPCs$Height_Flower, main='Normal')
qqline(PCA_Sum_CatPCs$Height_Flower)
#Awful

#Fecundity
qqnorm(PCA_Sum_CatPCs$Fecundity, main='Normal')
qqline(PCA_Sum_CatPCs$Fecundity)
#Not bad but not great.


# Kruskal Wallis tests ----------------------------------------------------

#Rosette_D
Rosette_D_KW <- kruskal.test(Rosette_D ~ PCsCAT, data = PCA_Sum_CatPCs)
print(Rosette_D_KW)
#x2= 0.0006, df=1, p=0.9811

#BushyScore
BushyScore_KW <- kruskal.test(BushyScore ~ PCsCAT, data = PCA_Sum_CatPCs)
print(BushyScore_KW)
#Kruskal-Wallis chi-squared = 0.063514, df = 1, p-value = 0.801


#DaysPlantToBolt
DaysPlantToBolt_KW <- kruskal.test(DaysPlantToBolt ~ PCsCAT, data = PCA_Sum_CatPCs)
print(DaysPlantToBolt_KW)
#Kruskal-Wallis chi-squared = 0.015047, df = 1, p-value = 0.9024

#Height_Harvest
Height_Harvest_KW <- kruskal.test(Height_Harvest ~ PCsCAT, data = PCA_Sum_CatPCs)
print(Height_Harvest_KW)
#Kruskal-Wallis chi-squared = 0.84156, df = 1, p-value = 0.359

#Height_Flower
Height_Flower_KW <- kruskal.test(Height_Flower ~ PCsCAT, data = PCA_Sum_CatPCs)
print(Height_Flower_KW)
#Kruskal-Wallis chi-squared = 0.32149, df = 1, p-value = 0.5707

#Fecundity
Fecundity_KW <- kruskal.test(Fecundity ~ PCsCAT, data = PCA_Sum_CatPCs)
print(Fecundity_KW)
#Kruskal-Wallis chi-squared = 0.60902, df = 1, p-value = 0.4352

#FLOWERTIME
FlowerTime_KW <- kruskal.test(FlowerTime ~ PCsCAT, data = PCA_Sum_CatPCs)
print(FlowerTime_KW)
#Kruskal-Wallis chi-squared = 0.45893, df = 1, p-value = 0.4981

