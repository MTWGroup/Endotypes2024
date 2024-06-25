#Project: Endotypes 
#
# Purpose: Clustering of the biomarkers
# Version: 1 
# Date: 21/06/22
# Author:HPF
#
# Input: DataSet
#       
# Output: Fully validated clusters of the biomarkers - to endotype
#         Cluster ensembles from a diverse set of algorithms
#         algorithms =  k-means,Hierarchical Clustering,GMM
#         Distance Measures = "euclidean"
#
# Add: N/A
#
# Dependencies:
#
# Notes: jump to script 3
#
# =    1 working space =========================================================
setwd("G:/")
list.files("G:/")
setwd("G:/")
list.files("G:/")
# = 2 packages need ============================================================

#Grafics  
library(ggplot2)
library(cowplot)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(plotly)

# Edting 
library(plyr)
library(doBy)
library(dplyr)
library(tidyverse)

# Analytic packages
library(factoextra)
library(NbClust)
library(cluster)
library(pvclust)
library(dendextend)
library(pheatmap)
library(survminer)
library(FactoMineR)
library(factoextra)
library(survival)
library(cowplot)
library(flextable)
library(finalfit)
library(WGCNA)
library(umap)
library(shipunov)
library(diceR)

# = 3 open databases ===========================================================

MWD<- read.csv("MWDataBIO.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
MWD$X=NULL

MWD_HC<-read.csv("PROFILEE.csv",header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("","NA"))
MWD_HC<-MWD_HC%>%select(patientid,colourHC1,ClusterH1)

# = 4 select data ==============================================================

MTData<-MWD%>%select(patientid,ProgY,Status,Dead_3Y,Dead_5Y)
MWB<-MWD%>%select(patientid,CYFRA211,PROC3,C6M,PROC6,C3M,REC1M,PROC28,XFIB,PROC4,SPD,CA199,CA125,MMP7)
rownames(MWB)<-MWB$patientid
MWB$patientid=NULL

# = 5 Normalize the customer data into the same scale ==========================

MWDx1D50Tr <- scale(t(MWB))
MWDx1D50<-t(MWDx1D50Tr)
sapply(MWDx1D50, function(x) sum(is.na(x)))

# = 6 detection of outliers for HC =============================================

A <- adjacency(t(MWDx1D50), type = "distance") # adjency calcutaion 
k  <- as.numeric(apply(A, 2, sum)) - 1 # calculates the whole network connectivity
Z.k <- scale(k) # standardized connectivity
thresholdZ.k = -3.5 # Z.k value is below the threshold == from -2.5<x<-5 
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")
sampleTree <- hclust(as.dist(1 - A), method = "average")
datColors <- data.frame(outlierC = outlierColor)

#dev.off()
dpi=300
tiff("Detection of outliers.tiff", res=dpi, height=12*dpi, width=12*dpi)

plotDendroAndColors(sampleTree, groupLabels = names(datColors), colors = datColors,      
                    main = "detection of unclustering datapoints")

dev.off()

remove.samples = Z.k < thresholdZ.k | is.na(Z.k)
MWDx1D50 = MWDx1D50[!remove.samples, ]

# = 7 Assessing clustering tendency ============================================
#To assess the clustering tendency, the Hopkins' statistic and a visual approach can be used. 
#This can be performed using the function get_clust_tendency() [factoextra package], 
#which creates an ordered dissimilarity image (ODI).
#Hopkins statistic: If the value of Hopkins statistic is close to 1 (far above 0.5), 
#then we can conclude that the dataset is significantly clusterable

OptCl50<-fviz_nbclust(MWDx1D50, FUN= hcut, method = "wss", k.max= 20)+
  labs(subtitle = "", tag = " ", title = " ")
OptCl50

# Find best number of clusterw
set.seed(1990)
# list of potential indices (note: I removed some which took longer to run)
indices <- c("kl", "ch", "hartigan", "ccc", "scott",
             "marriot", "trcovw", "tracew", "friedman", "rubin",
             "cindex", "db", "silhouette", "duda", "pseudot2", 
             "beale", "ratkowsky", "ball", "ptbiserial", "gap",
             "frey", "mcclain", "dunn", "sdindex", "sdbw")

# initialize var to collect clustering results
results <- list()
# loop over indices with try function (which continues running even with errors)
for (i in 1:length(indices)) {
  print(paste0("Trying ", indices[i], " index..."))
  results[[i]] <- try(NbClust(data=MWDx1D50,min.nc=2,max.nc=10, index=indices[i], method="kmeans")) 
}

results[[23]]

num_clust <- list()
for (i in 1:length(results)){
  num_clust[[i]] <- try(as.numeric(results[[23]]$Best.nc[1]))
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

paste0("Based on a number of criteria, we will select ", getmode(num_clust), " clusters.")

Index1<-data.frame(t(data.frame(results[[1]]$Best.nc,
                                results[[2]]$Best.nc,
                                results[[3]]$Best.nc,
                                results[[11]]$Best.nc,
                                results[[12]]$Best.nc,
                                results[[13]]$Best.nc,
                                results[[14]]$Best.nc,
                                results[[15]]$Best.nc,
                                results[[16]]$Best.nc,
                                results[[17]]$Best.nc,
                                results[[18]]$Best.nc,
                                results[[19]]$Best.nc,
                                results[[20]]$Best.nc,
                                results[[21]]$Best.nc,
                                results[[22]]$Best.nc,
                                results[[23]]$Best.nc,
                                results[[24]]$Best.nc,
                                results[[25]]$Best.nc)))



ClsrTD50 <- get_clust_tendency(MWDx1D50, n = 3, gradient = list(low = "steelblue",high = "white"))    
ClsrTD50[1]


DisClD50 <- dist(MWDx1D50, method = "euclidean")


DispPlotD50<-fviz_dist(DisClD50, lab_size = 7, gradient = list(low = "steelblue",high = "white"))
DispPlotD50<-DispPlotD50+labs(title=" ",caption = "hopkins_stat=0.8408")+
  theme(plot.title = element_text(size=15, hjust = 0.5, face="bold"))
DispPlotD50

dev.off()
dpi=300
tiff("Determining the number of clusters Kmeans.tiff",  res=dpi, height=6*dpi, width=13*dpi)
grid.arrange(DispPlotD50,OptCl50,ncol=2)
dev.off()


# = 8  Custom distance function ================================================

CC <- consensus_cluster(MWDx1D50, 
                        nk = 3, 
                        p.item = 0.8, 
                        reps = 5,
                        hc.method = "ward.D2",
                        algorithms = c("km","hc","som","gmm"),
                        distance = c("euclidean"),
                        scale = F,
                        seed.data = 3084,
                        progress = FALSE)

co <- capture.output(str(CC))
strwrap(co, width = 80)

DF<-data.frame(CC)

# = 8 Compute consensus matrix   ===============================================

#dev.off()
dpi=300
tiff("Cluster Consessus Cl3.tiff",  res=dpi, height=10*dpi, width=14*dpi)

par(mfrow=c(2,2))

HC3<- CC[, , "HC_Euclidean", "3", drop = FALSE]
cm <- consensus_matrix(HC3)
dim(cm)

HC1 <- graph_heatmap(HC3)


KM3<- CC[, , "KM", "3", drop = FALSE]
cm <- consensus_matrix(KM3)
dim(cm)

KM1 <- graph_heatmap(KM3)


SOM3<- CC[, , "SOM", "3", drop = FALSE]
cm <- consensus_matrix(SOM3)
dim(cm)

SOM1 <- graph_heatmap(SOM3)

GMM3<- CC[, , "GMM", "3", drop = FALSE]
cm <- consensus_matrix(GMM3)
dim(cm)

GMM1 <- graph_heatmap(GMM3)

dev.off()


# = 9 Combine consensus summaries   ============================================

ccomb_matrix <- consensus_combine(CC, element = "matrix")
ccomb_class <- data.frame(consensus_combine(CC, element = "class"))

ccomp <- consensus_evaluate(MWDx1D50, CC, plot = F)
ccomp$pac

ctrim <- consensus_evaluate(MWDx1D50, CC, trim = TRUE, reweigh = FALSE, n = 2)
ctrim$trim.obj$top.list

rownames(ccomb_class)<-rownames(MWDx1D50)

ccomb_classBest<-data.frame(ctrim$trim.obj$E.new)

ccomb_class<-ccomb_class%>%mutate(X3.GMM=case_when(X3.GMM==1~ 2,X3.GMM==2 ~1,X3.GMM==3~3))

Vote2<-data.frame(apply(ccomb_class,1,function(x) names(which.max(table(x)))))
names(Vote2)[1]<-paste("ClusterC2")

NClusterCon2<-data.frame(Vote2%>%group_by(ClusterC2)%>%dplyr:: summarise(n = n()))

Vote2<-Vote2%>%mutate(colourCon2=case_when(ClusterC2==3 ~ "darkgreen", ClusterC2==2 ~ "darkred", ClusterC2==1 ~ "royalblue4"))
Vote2$patientid<-rownames(Vote2)
Vote2<-Vote2%>%select(patientid,ClusterC2,colourCon2)


ccomb_class$patientid<-rownames(ccomb_class)
names(ccomb_class)[1:4]<-paste(c("H_Cluster","KMEANS","SOM","GMM"))
ccomb_class<-ccomb_class%>%select(patientid,KMEANS,H_Cluster,SOM,GMM)

MTData<-join(MTData,Vote2,by="patientid")
MTData<-join(MTData,ccomb_class,by="patientid")
MTData<-na.omit(MTData)


datColors2B <- data.frame(MTData$patientid,MWD_HC$ClusterH1,MTData$H_Cluster,MTData$KMEANS,MTData$SOM,MTData$GMM,MTData$colourCon2,MWD_HC$colourHC1)
names(datColors2B)[1:8]<-paste (c("patientid","H_Cluster1","H_Cluster","KMEANS","SOM","GMM","colourCon2","colourHC1"))

datColors2B<-datColors2B%>%mutate(colourHC=case_when(H_Cluster==3 ~ "darkgreen", H_Cluster==2 ~ "darkred", H_Cluster==1 ~ "royalblue4"))%>%
  mutate(colourKM=case_when(KMEANS==3 ~ "darkgreen", KMEANS==2 ~ "darkred", KMEANS==1 ~ "royalblue4"))%>%
  mutate(colourSOM=case_when(SOM==3 ~ "darkgreen", SOM==2 ~ "darkred", SOM==1 ~ "royalblue4"))%>%
  mutate(colourGMM=case_when(GMM==3 ~ "darkgreen", GMM==2 ~ "darkred", GMM==1 ~ "royalblue4"))

rownames(datColors2B)<-datColors2B$patientid
datColors2B$patientid=NULL


# = 10 consensus evaluation ====================================================

set.seed(2022)
HC_3 <- ccomb_class$`3`[, "HC_Euclidean"]
sig_objHC <- sigclust(MWDx1D50, k = 3, nsim = 100, labflag = 0, label = DIA_3)
co <- capture.output(str(sig_objHC))
strwrap(co, width = 80)


set.seed(2022)
KM_3 <- ccomb_class$`3`[, "KM"]
sig_objKM <- sigclust(MWDx1D50, k = 3, nsim = 100, labflag = 0, label = KM_3)
coKm <- capture.output(str(sig_objKM))
strwrap(coKm, width = 80)


set.seed(2022)
GMM_3 <- ccomb_class$`3`[, "GMM"]
sig_objGMM <- sigclust(MWDx1D50, k = 3, nsim = 100, labflag = 0, label = GMM_3)
coGmm <- capture.output(str(sig_objGMM))
strwrap(coGmm, width = 80)


set.seed(2022)
SOM_3 <- ccomb_class$`3`[, "SOM"]
sig_objSOM <- sigclust(MWDx1D50, k = 3, nsim = 100, labflag = 0, label = GMM_3)
coSOM <- capture.output(str(sig_objSOM))
strwrap(coSOM, width = 80)



# = 11 grafical evaluation of all clustering ===================================

A2 = adjacency(t(MWDx1D50), power = 12, type='unsigned')
dissTOM = TOMdist(A2)
geneTreeA = hclust(as.dist(dissTOM), method = "average")
plot(geneTreeA)

dev.off()
hc50 <- hclust(dist(MWDx1D50, method = "euclidean"), method = "ward.D2")
plot(hc50)

datColors2BS<-data.frame(datColors2B[ ,c("colourHC1","colourHC","colourKM","colourSOM","colourGMM","colourCon2")])


dev.off()
dpi=300
tiff("Cluster Dendogram Cl3All.tiff",  res=dpi, height=10*dpi, width=17*dpi)

plotDendroAndColors(hc50, colors = datColors2BS, groupLabels = c("H Cluster 1","H Cluster","KMEANS","SOM","GMM","Consensus"), dendroLabels = FALSE, hang = 0.03,      
                    addGuide = TRUE, guideHang = 0.05, main ="Algorithmic comparison" )


dev.off()

# = 11 eigenvalues and clusters =================================================

MEs0 = moduleEigengenes(t(MWDx1D50), MTData$colourCon2)$eigengenes
MEsC = orderMEs(MEs0)

nGenes = ncol(t(MWDx1D50)) 
nSamples = nrow(t(MWDx1D50))
modTraitCor = cor(MEsC ,use = "p") 
modTraitP = corPvalueStudent(modTraitCor, nSamples)

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")",      
                   sep = "") 
dim(textMatrix) = dim(modTraitCor)


cols <- colorRampPalette(brewer.pal(5, "RdBu"))(256)



#dev.off()
dpi=100
tiff("Cluster Intra Module relationships Cl3.tiff",  res=dpi, height=10*dpi, width=17*dpi)

labeledHeatmap(Matrix = modTraitCor, xLabels = names(MEsC), 
               yLabels = names(MEsC),ySymbols = names(MEsC), 
               colorLabels = FALSE, colors = rev(cols),
               textMatrix = textMatrix, setStdMargins =FALSE, cex.text = 1, 
               zlim = c(-1,1), main = paste("Intra Module relationships"))
dev.off()



# = 12  calculate the module membership values (module connectivity kME) =======

moduleColors <- MTData$colourCon2
datKME <- signedKME(t(MWDx1D50), MEsC)
colorOfColumn <- substring(names(datKME),4)

cmd1 = cmdscale(as.dist(dissTOM), 2)
par(mfrow = c(2,2))
plot(cmd1, col = MTData$colourCon2, main = "MDS plot Patient Vs Biomarkers", xlab = "Scaling Dimension 1", 
     ylab = "Scaling Dimension 2")


cmd2=cmdscale(as.dist(dissTOM),4)
par(mfrow=c(2,3))
plot(cmd2[,c(1,2)], col= as.character(MTData$colourCon2) )
plot(cmd2[,c(1,3)], col= as.character(MTData$colourCon2) )
plot(cmd2[,c(1,4)], col= as.character(MTData$colourCon2) )
plot(cmd2[,c(2,3)], col= as.character(MTData$colourCon2) )
plot(cmd2[,c(2,4)], col= as.character(MTData$colourCon2) )
plot(cmd2[,c(3,4)], col= as.character(MTData$colourCon2) )


library(scatterplot3d)


dev.off()
dpi=300
tiff("MDS plot Cl3.tiff",  res=dpi, height=8*dpi, width=10*dpi)
par(mfrow=c(1,1))

scatterplot3d(cmd2[,1:3], color=as.character(MTData$colourCon2),  
              xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", zlab="Scaling Dimension 3"
              ,cex.axis=1,angle=320,pch=16, main = " ")


dev.off()




# = 11 UMAP visualization biomarker patients ===================================

MTData<-cbind(MTData,datKME,MWDx1D50)
#names(MTData)[2:13]<-paste (c("CYFRA211_NN","PROC3_NN","C6M_NN","PRC6_NN","C3M_NN","REC1M_NN","PROC28_NN",
#                              "XFIB_NN","PROC4_NN","SPD_NN","CA199_NN","CA125_NN"))


MTData<-join(MTData,MWD_HC,by="patientid")

write.table(MTData, "PROFILEEndoClustOutPut300622_455x29.csv", row.names = F, sep = ",")


UMAPdata = data.frame(MTData[15:27])
UMAP.labels = MTData[ ,"ClusterC2"]


MDumap = umap(UMAPdata)
MDumap$layout

colors<-MTData%>%group_by(colourCon2)%>%dplyr:: summarise(n = n())

plot.UMAP<-function(x, labels,
                    main=" ",
                    colors=c("royalblue4","darkred","darkgreen"),
                    pad=0.1, cex=1.5, pch=20, add=FALSE, legend.suffix="",
                    cex.main=1, cex.legend=1.5) {
  layout = x
  if (is(x, "umap")) {
    layout = x$layout   }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)  
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topleft"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

#dev.off()
dpi=300
tiff("Cluster UMAP vissaulization 226022.tiff",  res=dpi, height=8*dpi, width=8*dpi)
plot.UMAP(MDumap,UMAP.labels, cex = 1.5 )
dev.off()

