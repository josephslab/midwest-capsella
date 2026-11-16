#PCA plot

library(dplyr)
library(ggplot2)

#load in data
midwestEv<-read.table('../data/midwestCapsella-pca.eigenvec')
names(midwestEv)[2] <- "id"
midwestPopInfo <- read.table('../data/midwest_pixy_pops.txt')
names(midwestPopInfo) <- c("id",'pop')
midwestPopInfo$pop = factor(midwestPopInfo$pop)

#add pop info
midwestEv <- dplyr::left_join(midwestEv, midwestPopInfo, by="id")


#get perc var explained
system("awk '{print $NR}' ../data/midwestCapsella-pca.rel > ../data/midwestCapsella-pca.rel.diag")
mydiag = read.table('../data/midwestCapsella-pca.rel.diag')
totalVar = sum(mydiag$V1)
evals = read.table('../data/midwestCapsella-pca.eigenval')
evals$perc = evals$V1/totalVar

##make pca plot colored by city
scatterV3V4<-ggplot(midwestEv, aes(x=V3, y=V4)) + 
  geom_point(aes(shape = pop, colour=pop, fill=pop), size = 7) +
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
    axis.line=element_line(linewidth = 0.5, colour = "black"),
    axis.text=element_text(colour="black", size = 20, family="A"),
    axis.title=element_text(size=25, family = "A")
  )+
  #theme(legend.text = element_text(size = 20, family="A"),
   #     legend.title = element_text(size = 25, family="A"))+
  labs (y = paste("PC 2 (",round(evals$per[2],digits=4)*100,"%)",sep=""), 
        x= paste("PC 1 (",round(evals$per[1],digits=4)*100,"%)",sep=""))






scatterV3V4
#while (!is.null(dev.list()))  dev.off()





