#Project: Endotypes (Scrip 2) 
#
# Purpose:Clustering of the biomarkers - Random-Walk
# Version: 2  
# Date: 25/10/22
# Author:HPF
#
# Input: Dataset - cluster bio markers 
#       
# Output: minimal clustering validation -by random walker 
#
# Add: N/A
#
# Dependencies: NRHI 10 Random Walk 251022
#
# Notes: jump from NHRI 7 

# =    1 working space =========================================================

setwd("F:/")
list.files("F:/")
#if (! dir.exists("Project_IPF")) dir.create("Project_IPF")
setwd("F:/NRHI 10 RandomWalk 251022")
list.files("F:/NRHI 10 RandomWalk 251022")

# = 2 packages need ============================================================
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")  
#BiocManager::install(c("limma"))

library(egg) # for ggarrange
library(inspectdf) # check entire data.frame for variable types, etc.
library(visdat) # visualize missing data
library(corrr)# correlations
library(tidyverse) # all tidyverse packages
library(tidymodels) # meta package for modeling
library(psych) # for skewness and kurtosis
library(socviz) # for %nin%
library(WGCNA)
library(viridis)

# = 3 open database ============================================================

MWD<-read.csv("PROFILEEndoClustOutPut300622_455x29.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))


#= 4 Data selection ============================================================

MWD<-MWD%>%select(patientid,CYFRA211,PROC3,C6M,PROC6,C3M,REC1M,PROC28,XFIB,PROC4,SPD,CA199,CA125,MMP7)

rownames(MWD)<-MWD$patientid
MWD$patientid=NULL

#= 5 data processing  ==========================================================

A2 = adjacency((MWD),power=12 ,type='distance')   #distance matrix by adjacency
prob_mat = TOMdist(A2)
prob_matdf<-data.frame(prob_mat)
rownames(A2)->rownames(prob_matdf)
rownames(A2)->colnames(prob_matdf)
A2df<-data.frame(A2)

# = 6 Graphic networks  ========================================================

rplot_labs_x1 <- ggplot2::labs( title = "Distances between Biomarkers",
                                subtitle = " ",
                                caption = "PROFILE") 

A2df %>% corrr::rplot(rdf = ., shape = 19,
                      colors = c("firebrick", "dodgerblue")) +  rplot_labs_x1

network_labs_x1 <- ggplot2::labs( title = "Biomarker network",
                                  caption = "PROFILE")

network_plot_x1 <- A2df %>% dplyr::select_if(is.numeric) %>% 
                            corrr::correlate(x = ., 
                            method = "spearman", quiet = TRUE) %>%  corrr::network_plot(rdf = ., 
                            colors = c("firebrick", "white", "dodgerblue"),
                            min_cor = .50) + 
network_labs_x1
network_plot_x1<-network_plot_x1+theme(plot.title = element_text(size=15, hjust = 0.5, face="bold"))
network_plot_x1    



# = 7 Random Walker ============================================================

set.seed(2021)

find_pos_prob <- function(prob_mat, no_of_steps){
  
  x <- c(1:nrow(prob_mat))               # index for nodes
  position <- 1                          # initiating from 1st Node
  occured <- rep(0,nrow(prob_mat))       # initiating occured count
  
  for (i in 1:no_of_steps)   {
    # update position at each step and increment occurence
    position  <-  sample(x, 1, prob = prob_mat[position,])      
    occured[position] <- occured[position] + 1
  }
  return (occured/no_of_steps)
}

RW<-data.frame(find_pos_prob(prob_mat, 1000)) 
RW<-RW%>%mutate(name=colnames(MWD))
names(RW)[1:2]<-paste(c("Probability","Feature"))
RW<-RW%>%arrange(-Probability)
str(RW)
RW<-RW%>%mutate(Colour=case_when(Probability>0.094 ~ 1,Probability<0.063 ~ 3, Probability<0.095 & Probability>0.063 ~ 2))

IMPPredictors<-ggplot(RW[1:13,])+
               geom_bar(aes(x=reorder(Feature,Probability), y=Probability, fill=as.factor(Colour)) ,stat='identity')+
               xlab(label = "Top Probabilities")+
               theme_bw() +
               scale_fill_manual("legend", values = c("1" = "darkgreen", "2" ="royalblue", "3" = "firebrick"))+  
               geom_hline(yintercept= 0.081, linetype="dashed", color = "red", size=0.5)+
               geom_hline(yintercept= 0.062, linetype="dashed", color = "red", size=0.5)+
               labs(x="Top Probabilities",y="Probabilities", title="Top Probabilities for 3 clusters ",caption = " ")+
  theme(plot.title = element_text(size=15, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=12, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=12, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=10, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=10, angle=0))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  coord_flip()

IMPPredictors 


IMPPredictors2 <-ggplot(RW) + geom_bar(aes(reorder(Feature,Probability),Probability, fill = (-Probability)), stat = 'identity') + 
  scale_fill_viridis_c(option = 'magma') +
  theme_bw() +
  theme(legend.position = "none")+
  geom_hline(yintercept= 0.075, linetype="dashed", color = "red", size=0.5)+
  geom_hline(yintercept= 0.062, linetype="dashed", color = "red", size=0.5)+
  labs(x="Top Probabilities",y="Probabilities", title="Top Probabilities for 3 clusters ",caption = " ")+
  theme(plot.title = element_text(size=15, hjust = 0.5, face="bold"),
        axis.title.x = element_text(face="bold", size=12, hjust = 0.5),
        axis.title.y = element_text(face="bold", size=12, hjust = 0.5)) +
  theme(axis.text.x = element_text(face="bold", color="#993333",size=10, angle=0),
        axis.text.y = element_text(face="bold", color="#993333",size=10, angle=0))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  coord_flip()

IMPPredictors2

dev.off()
dpi=100
tiff("Network factors.tiff",  res=dpi, height=12*dpi, width=20*dpi)
grid.arrange(IMPPredictors,network_plot_x1,ncol=2)
dev.off()

write.csv(RW, file="TopPredictors 251022.csv", na="")