# Financial stocks clustering based on ordinary VAR coherence 

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

# ###############################################################################
# Define frequencies 
# ###############################################################################
#Fourier frequencies
nf<-5739 #2010-2019
ff<-c(0:(nf-1))/nf
sel.f <- which(ff < 0.5 & ff > 0)
frequ<-ff[sel.f]

# ###############################################################################
# load the data. Financial time series
# ###############################################################################
stocks_names <- readRDS("./stocks_names.rds")
file<-"Financial_data.csv" # All 52 financial time series. Time period: 2010-2019
x0<-read.csv(file,header=T) #Column 1 is date, column 10 the S&P 500 index.
sp500<-x0[,10] # Financial index S&P500
Stocks<-x0[,-c(1,10)] #52 Financial stocks 
y2<-diff(log(sp500)) #Compute the logarithmic differences for the S&P financial index.

#Function for the estimation of the ordinary VAR coherence
source("VAR_coh.R")

# ###############################################################################
# Estimation of the Ordinary VAR coherence
# ###############################################################################
# The order of the model is automatically selected by the AIC criterion. 
ord_coherence<-matrix(0,nrow=(length(sp500)-1)/2,ncol=dim(Stocks)[2])
Order_selected<-c()
for (w in 1:dim(Stocks)[2]) {
  y1<-diff(log(Stocks[,w])) #Compute the logarithmic differences for the stock w
  VAR2coh<-var_coh(x=cbind(y1,y2),order.max=10,aic=TRUE)
  ord_coherence[,w] <- VAR2coh$var_coh[sel.f]
  Order_selected[w] <- VAR2coh$order
}

##############################################################################
# Compute the clusters based on OC for comparisons
##############################################################################


# Compute the distance matrix
norm_coh<-as.matrix(stats::dist(t(ord_coherence)))
rownames(norm_coh)<- stocks_names


# Hierarchical clustering based on OC
clusters_ord <- hclust(d=as.dist(norm_coh))
plot(clusters_ord) 

# Elbow plot for selecting the optimal number of clusters
source("fviz_nbclust.R")
fviz_nbclust2(ord_coherence, hcut, method = "wss") 


clusterCut_ord <- cutree(clusters_ord,3)
which(clusterCut_ord==1)
which(clusterCut_ord==2)
which(clusterCut_ord==3)


# ###############################################################################
# # Exercise 2. Estimation of the Ordinary VAR coherence
# ###############################################################################
# The order of the model is 1. AIC=FALSE 
ord_coherence<-matrix(0,nrow=(length(sp500)-1)/2,ncol=dim(Stocks)[2])
for (w in 1:dim(Stocks)[2]) {
  y1<-diff(log(Stocks[,w])) #Compute the logarithmic differences for the stock w
  VAR2coh<-var_coh(x=cbind(y1,y2),order.max=1,aic=FALSE)
  ord_coherence[,w] <- VAR2coh$var_coh[sel.f]
}

##############################################################################
# Compute the clusters based on OC for comparisons
##############################################################################


# Compute the distance matrix
norm_coh<-as.matrix(stats::dist(t(ord_coherence)))
rownames(norm_coh)<- stocks_names


# Hierarchical clustering based on OC
clusters_ord <- hclust(d=as.dist(norm_coh))
plot(clusters_ord) 

# Elbow plot for selecting the optimal number of clusters
source("fviz_nbclust.R")
fviz_nbclust2(ord_coherence, hcut, method = "wss") 


clusterCut_ord <- cutree(clusters_ord,3)
which(clusterCut_ord==1)
which(clusterCut_ord==2)
which(clusterCut_ord==3)


# ###############################################################################
# # Exercise 3. Estimation of the Ordinary VAR coherence
# ###############################################################################
# The order of the model is automatically selected by AIC. maximum.order=2.AIC=TRUE 
ord_coherence<-matrix(0,nrow=(length(sp500)-1)/2,ncol=dim(Stocks)[2])
for (w in 1:dim(Stocks)[2]) {
  y1<-diff(log(Stocks[,w])) #Compute the logarithmic differences for the stock w
  VAR2coh<-var_coh(x=cbind(y1,y2),order.max=2,aic=TRUE)
  ord_coherence[,w] <- VAR2coh$var_coh[sel.f]
}

##############################################################################
# Compute the clusters based on OC for comparisons
##############################################################################


# Compute the distance matrix
norm_coh<-as.matrix(stats::dist(t(ord_coherence)))
rownames(norm_coh)<- stocks_names


# Hierarchical clustering based on OC
clusters_ord <- hclust(d=as.dist(norm_coh))
plot(clusters_ord) 

# Elbow plot for selecting the optimal number of clusters
source("fviz_nbclust.R")
fviz_nbclust2(ord_coherence, hcut, method = "wss") 


clusterCut_ord <- cutree(clusters_ord,3)
which(clusterCut_ord==1)
which(clusterCut_ord==2)
which(clusterCut_ord==3)


