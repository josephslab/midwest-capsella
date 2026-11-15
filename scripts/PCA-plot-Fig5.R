#PCA plot

#The sum of the eigenvalues is the same as the sum of the diagonal in the covariance matrix that is created when you use the --make-rel flag in plink. This sum of the diagonal is the overall variability.

#To get the proportion of variance, divide the eigenvalues by this "sum of the diagonal". Alternatively, divide by the sum of the eigenvalues - BUT ONLY if you are sure you have returned all of the principal components (the default in PLINK is just to show the top 20).

#A quick way to get just the diagonals from the plink.rel file is to use the awk code:
  
#  awk '{print $NR}' plink.rel > plink.rel.diag





CityPCA<-read.csv("C:/Wondering Researcher/Academics/Michigan State University/Chrissy/UrbanPlantProject/Genotypic data/MidwestEigenvec_Feb.csv")


PCACityExcludingGROutlier <- CityPCA[-c(14), ]

PCACityExcludingGROutlier$City<-as.factor(PCACityExcludingGROutlier$City)

library(ggplot2)

windowsFonts(A = windowsFont("Times New Roman"))
tiff ("PCA_WithoutOutlier_VarianceExplained.tiff", width = 8, height = 8, units = 'in', res = 80)
scatterV3V4<-ggplot(PCACityExcludingGROutlier, aes(x=PCA1, y=PCA2)) + 
  geom_point(aes(shape = City, colour=City, fill=City), size = 7) +
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
    legend.position = c(0.80, 0.80),
    legend.background = element_rect(colour = NA),
    legend.key = element_rect(colour = "white", fill = NA),
    axis.line=element_line(size = 0.5, colour = "black"),
    axis.text=element_text(colour="black", size = 20, family="A"),
    axis.title=element_text(size=25, family = "A")
  )+
  theme(legend.text = element_text(size = 20, family="A"),
        legend.title = element_text(size = 25, family="A"))+
  labs (y = "PC 2 (explains 9.64% of varience)", x= "PC 1 (explains 11.96% of variation)")
scatterV3V4
#while (!is.null(dev.list()))  dev.off()





