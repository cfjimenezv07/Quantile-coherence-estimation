# Clustering FTS based on QC
# load the required packages
packages <- c("QZ", "matrixcalc", "quantreg","fields","colorRamps", "cluster", 
              "factoextra","graphics","AR","TSA","MTS","tsDyn","qfa","parallel")
## Now load or install&load all
package_check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

color.qper<-topo.colors(1024)
color2.qper<-rev(color.qper)
color.cqper<-rainbow(1024)
color.qper<-matlab.like2(1024)

# Set the directory
setwd("~/My Drive/Spring 2023/STAT 397/PhD project 2/Codes_paper/")
dir.d<- "./"

################################################################################
# Clustering Financial time series
################################################################################
#Load the QC maps and compute the features
Smoothed_financial_stocks <- readRDS("Stocks_1.rds")
nf<-dim(Smoothed_financial_stocks[[1]])[1]
ff<-c(0:(nf-1))/nf
sel.f <- which(ff < 0.5 & ff > 0)
frequ<-ff[sel.f]

# Quantile levels
tau.min <- 0.04
tau.max <- 0.96
tau.del<- 0.01
tau<-seq(tau.min,tau.max,tau.del)
ntau<-length(tau)
sel.tau<-which(tau >= 0.05 & tau <= 0.95)

# Compute the features for the distance matrix for clustering
len_vec=dim(Smoothed_financial_stocks[[1]])[1]*dim(Smoothed_financial_stocks[[1]])[2]
features_extraction<-matrix(0,nrow=length(Smoothed_financial_stocks),ncol = len_vec)
for (i in 1:length(Smoothed_financial_stocks)) {
  mat=Smoothed_financial_stocks[[i]]
  features_extraction[i,]=fBasics::vec(mat)
}


#Compute the distance
dist_norm<-as.matrix(stats::dist(features_extraction))
row.names(dist_norm)<- stocks_names
colnames(dist_norm)<-row.names(dist_norm)
clusters <- hclust(d=as.dist(dist_norm))
plot(clusters)

# Elbow plot for selecting the optimal number of clusters
source("fviz_nbclust.R")
fviz_nbclust1(x=features_extraction, FUNcluster =hcut, method = "wss")


clusterCut <- cutree(clusters,3)
which(clusterCut==2)
which(clusterCut==1)
which(clusterCut==3)
Clusters<-list(which(clusterCut==2),which(clusterCut==1),which(clusterCut==3))



################################################################################
# Some plots by clusters
################################################################################


Max=max(sapply(Smoothed_financial_stocks, max))
for (ww in 1:length(Clusters)) {
  Cluster<-Smoothed_financial_stocks[Clusters[[ww]]]
  Stocks_names<-stocks_names[Clusters[[ww]]]
  Beta=beta[Clusters[[ww]]]
  dir.l<-"~/My Drive/Spring 2023/STAT 397/PhD project 2/Plots_QC_od1/"
  pdf(paste0(dir.l, "QC_maps_cluster_",ww,".pdf"))
  
  for (i in 1:length(Cluster)) {
    qper1=Cluster[[i]]
    tlab6<-paste("Quantile coherence:", Stocks_names[i],",","Beta=",Beta[i] )
    qfa.plot(ff[sel.f],tau[sel.tau],qper1,rg.qper=range(c(0,Max)),rg.tau=range(tau[sel.tau]),color=color.qper,tlab=tlab6)
    
  }
  
  dev.off()
}
