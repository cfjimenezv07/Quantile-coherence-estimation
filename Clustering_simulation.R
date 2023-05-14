#Clustering_simulation with smoothing
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
## Load the simulated bivariate time series and their respective QC estimations
################################################################################
#They are 100 samples from each of the 4 models defined in the paper
Simulated_time_series <- readRDS("Simulated_time_series.rds")
# To each of the bivariate time series we compute the QC based on the proposed 
# semi parametric estimator.
QC_simulated_time_series_smoothed <- readRDS("./QC_simulated_time_series_smoothed.rds")
ns=100 #Number of time series per model
nm=4   # Number of simulated models

###############################################################################
# Clustering simulation with the estimation of the QC
###############################################################################
# Get the features from the QC estimation
Features_extraction<-list()
for (ii in 1:length(QC_simulated_time_series_smoothed)) {
  QC_map<-QC_simulated_time_series_smoothed[[ii]]
  len_vec=dim( QC_map[[1]])[1]*dim(QC_map[[1]])[2]
  
  features_extraction<-matrix(0,nrow=length(QC_map),ncol = len_vec)
  for (jj in 1:length(QC_map)) {
    mat<-QC_map[[jj]]
    features_extraction[jj,]<-fBasics::vec(mat)
  }
  
  Features_extraction[[ii]]<-features_extraction
}


Final_extracted_features_smoothed<-matrix(0,nrow = ns*nm,ncol = dim(Features_extraction[[3]])[2])
for (ij in 1:length(Features_extraction)) {
  features_model<-Features_extraction[[ij]]
  Final_extracted_features_smoothed[(ns*ij-(ns-1)):(ns*ij),1:dim(features_model)[2]]<-features_model
}

#Compute the distance for hierarchical Clustering
dist_norm<-as.matrix(stats::dist(Final_extracted_features_smoothed))
rownames(dist_norm)<-c(rep(1,ns),rep(2,ns),rep(3,ns),
                       rep(nm,ns))
colnames(dist_norm)<-rownames(dist_norm)
clusters <- hclust(d=as.dist(dist_norm))
#plot(clusters)
clusterCut <- cutree(clusters,4)


# Elbow plot for selecting the optimal number of clusters
source("fviz_nbclust.R")
fviz_nbclust1(Final_extracted_features_smoothed, hcut, method = "wss") 


#####################################################################################
# Simulation clustering based on VAR ordinary coherence
######################################################################################

#(a)fit a vector AR model to the original vector time series using the function ar;
VAR_pred<-list()
AR.coeff <- list()
AR.order<-list()
for (i in 1:length(Simulated_time_series)) {
  X_1<-Simulated_time_series[[i]][[1]]
  X_2<-Simulated_time_series[[i]][[2]]
  var.pred<-list() 
  ar.coeff<-list()
  aic_order<-list()
  for (j in 1:dim(X_1)[2]) {
    coeff=ar(x=cbind(X_1[,j],X_2[,j]))
    ar.coeff[[j]]=aperm(coeff$ar,c(2,3,1))
    var.pred[[j]]=coeff$var.pred
    order <- coeff$order
    if(order>=10){
      aic_order[[j]]=10
    }else{
      aic_order[[j]]=order
    }
    
  }
  VAR_pred[[i]]<- var.pred
  AR.coeff[[i]]<- ar.coeff
  AR.order[[i]]<- aic_order
}

#(b) compute the spectral matrix using the parameters from the fitted vector 
# AR model.
var2coh<-function(x,ar.order,var.pred,ar.coeff){
  nc=2
  nf=dim(x)[1]
  freq <- c(0:(nf - 1)) / nf
  varspec=array(data=NA,dim=c(nc,nc,nf))
  order_=ar.order
  varcoeff=ar.coeff
  Sigma=var.pred
  for (v in 1:nf ) {
    u=matrix(0,2,2)
    for (r in 1:order_) {
      u=u+varcoeff[,,r]*complex(real = cos(2*pi*r*freq[v]),imaginary =-sin(2*pi*r*freq[v]))
    }
    U=diag(nc)-u
    varspec[,,v]=solve(U)%*%Sigma%*%QZ::H(solve(U))
  }
  
  list(spec=varspec)
}

# Compute the spectral matrix for all time series
VAR_spec<-list()
for (i in 1:length(Simulated_time_series)) {
  X_1<-Simulated_time_series[[i]][[1]]
  X_2<-Simulated_time_series[[i]][[2]]
  var.pred  <-  VAR_pred[[i]]
  ar.coeff  <-  AR.coeff[[i]]
  aic_order <-  AR.order[[i]]
  var_spec   <-   list()
  for (j in 1:dim(X_1)[2]) {
    var_spec[[j]] <- var2coh(x=cbind(X_1[,j],X_2[,j]),ar.order=aic_order[[j]],
                             var.pred=var.pred[[j]],ar.coeff=ar.coeff[[j]])$spec
  }
  VAR_spec[[i]]<-var_spec
}

#(c) compute the coherence from the AR spectral matrix.

k<-1
kk<-2
VAR_Coh<-list()
for (i in 1:length(Simulated_time_series)) {
  n <-dim(VAR_spec[[i]][[i]])[3]
  ff<-c(0:(n-1))/n
  nf<-length(ff)
  sel.f <- which(ff < 0.5 & ff > 0)
  var.spec <- VAR_spec[[i]]
  VAR_coh<-matrix(0,nrow=n,ncol = length(var.spec))
  for (jj in 1:dim(VAR_coh)[2]) {
    var.coh_ <- var.spec[[jj]]
    coh<- c()
    for (ij in 1:n) {
      var.coh <-var.coh_[,,ij]
      num<-Mod(var.coh[k,kk])^2
      den<-(Re(var.coh[kk,kk])*Re(var.coh[k,k]))
      coh[ij]<-num/den
    }
    VAR_coh[,jj]<-coh
  }
  VAR_Coh[[i]]<-VAR_coh
}

##############################################################################
# Simulation clustering based on VAR coherence
###############################################################################

# Feature extraction from VAR coherence
Final_extracted_features_VAR_coh<-matrix(0,nrow = ns*nm,ncol = 512)
for (ij in 1:length(VAR_Coh)) {
  features_model<-t(VAR_Coh[[ij]])
  Final_extracted_features_VAR_coh[(ns*ij-99):(ns*ij),1:dim(features_model)[2]]<-features_model
}

norm_coh<-as.matrix(dist(Final_extracted_features_VAR_coh))

rownames(norm_coh)<-c(rep("Sim 1",ns),rep("Sim 2",ns),rep("Sim 3",ns),
                      rep("Sim 4",ns))
colnames(norm_coh)<-c(rep("Sim 1",ns),rep("Sim 2",ns),rep("Sim 3",ns),
                      rep("Sim 4",ns))

clusters_VAR_coh <- hclust(d=as.dist(norm_coh))
plot(clusters_VAR_coh)
clusterCut_VAR_coh <- cutree(clusters_VAR_coh,4)
c_1<-which(clusterCut_VAR_coh==1)
c_2<-which(clusterCut_VAR_coh==2)
c_3<-which(clusterCut_VAR_coh==3)
c_4<-which(clusterCut_VAR_coh==4)

# Elbow plot for selecting the optimal number of clusters
fviz_nbclust2(Final_extracted_features_VAR_coh, hcut, method = "wss") 



