# Alternative Methods for the simulation settings

# load the required packages
packages <- c("QZ", "matrixcalc", "quantreg","fields","colorRamps", "cluster", 
              "factoextra","graphics","AR","TSA","MTS","tsDyn","qfa")
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
# Define frequencies and quantile levels
# ###############################################################################
#Fourier frequencies
nf<-512
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

# ###############################################################################
# Simulate a dataset. For example data from model 3
# ###############################################################################
n<-512
m<-5000
set.seed(2000)
xx<-matrix(NA,nrow=n,ncol=m)
xx2<-xx
for(i in c(1:m)) {
  
  # - component 1 (lowpass) -
  a1<- 0.80
  tmp<-ar.sim(20*n,a1,1)[-c(1:(19*n))]
  tmp<-(tmp-mean(tmp))/sd(tmp)
  x1<-tmp[1:n]
  
  # - component 2 (highPASS) -
  a2<- -0.70
  tmp<-ar.sim(20*n,a2,1)[-c(1:(19*n))]
  tmp<-(tmp-mean(tmp))/sd(tmp)
  x2<-tmp[1:n]
  
  # - component 3 (bandpass) -
  r<-0.9         # bandwidth of AR(2)
  f.ar<-0.20     # peak frequency of AR(2)
  aa1<- 2*r*cos(2*pi*f.ar)
  aa2<- -r*r
  v<-(1-aa2)/(1+aa2)/((1-aa2)^2-aa1*aa1)
  tmp0<-ar.sim(21*n,c(aa1,aa2),1)[-c(1:(19*n))]
  tmp0<-(tmp0-mean(tmp0))/sd(tmp0)
  x3<-tmp0[1:n]
  # - delayed copy of component 3
  d<-10
  x4<-tmp0[c(1:n)+d]
  
  # - nonlinear mixing -
  w1<-0.8
  w2<-0.1
  y1<-0.8
  y2<--0.8
  y12<-0.4
  y22<--0.4
  w12<-0.0
  w22<-0.5
  wt<-psi.wts(x1,y1,y2,w1,w2)
  x0<-wt*x2+(1-wt)*x1
  wt2<-psi.wts(x0,y12,y22,w12,w22)
  x<-wt2*x3+(1-wt2)*x0
  
  xx[,i]<-x
  xx2[,i]<-x4
}

# ###############################################################################
# 0. Simulate the ground truth for the quantile spectral matrix
# ###############################################################################
n.cores <- detectCores()-2
if(n.cores>1) {
  library(foreach) 
  library(doParallel) # libraries required for parallel computing
  cl <- makeCluster(n.cores)
  clusterExport(cl, c("rq"))
  registerDoParallel(cl) 
}
qcper.sim<-list()
qacf.sim<-list()
for(i in c(1:m)) {
  y<-cbind(xx[,i],xx2[,i])
  y.qdft<-qdft(y,tau,n.cores = n.cores)
  qcper.sim[[i]]<-qdft2qcper(y.qdft)
  qacf.sim[[i]]<-qdft2qacf(y.qdft)
}
if(n.cores>1) stopCluster(cl)

# simulated quantile spectra and coherence
nc<-2
qcper.mean<-array(0,dim=dim(qcper.sim[[1]]))
for(i in c(1:m)) {
  for(k in c(1:nc)) {
    for(kk in c(1:nc)) {
      qcper.mean[k,kk,,]<-qcper.mean[k,kk,,]+qcper.sim[[i]][k,kk,,]
    }
  }
}
qcper.mean<-qcper.mean/m # The truth quantile spectral matrix

# ###############################################################################
# 1. Smoothing splines independently across frequencies and quantiles
# ###############################################################################

 # data: is a nf by nq matrix representing the number of considered frequencies and nq number of quantiles 
 #vetsmo: auxiliary function to perform smoothing splines across the rows of data. We use the usual GCV criterion.
vetsmo= function(data){
  m = dim(data)[1]
  n = dim(data)[2]
  x = matrix(0, m, n)
  for(i in 1:m){
    x[i,] = smooth.spline(data[i,])$y
  }
  return(x)
}

#Smooth_splines perform smoothing spline over the frequencies and then over the quantiles. 
smooth_splines<-function(data){
  sel.f=dim(data)[1]
  sel.tau=dim(data)[2]
  quant_peri_splines=matrix(0,sel.f,sel.tau)
  for(i in 1:sel.tau){
    quant_peri_splines[,i]      = smooth.spline(data[,i])$y}
  smooth_qcper    = vetsmo(quant_peri_splines)
  return(smooth_qcper)
}

# ###############################################################################
# 2. 2D kernel smoothing
# ###############################################################################
# here we use the 2D kernel smoothing for images. 
smooth_2Dkernel<-function(data){
  smooth_qcper_2DKernel = image.smooth(data)$z
  return(smooth_qcper_2DKernel)
}

spec.normalize<-function(qper) {
  if(is.matrix(qper)) {
    return( apply(qper,2,FUN=function(x) { x/sum(x) }) )
  } else {
    return( qper/sum(qper) )
  }
}


# ###############################################################################
# Semi-parametric estimator of the quantile coherence
# ###############################################################################

p_max <- 10 # Maximum order to be selected.
nc <- 2 #number of time series
nq <- 93 # number of quantiles
source("aux_semi_parametric_estimation_QC.R")
source("joint_smoothing.R")

# ###############################################################################
# True coherence estimate
# ###############################################################################

tmp0<-Mod(qcper.mean[k,kk,sel.f,sel.tau])^2/(Re(qcper.mean[k,k,sel.f,sel.tau])*Re(qcper.mean[kk,kk,sel.f,sel.tau]))
np<-length(c(tmp0))

# ####################################################################################################
# Compute the alternative methods for comparisons to the results with the semi parametric estimator
# ####################################################################################################
#RMSE for all the competitors
eps<-1e-6
m=200
rmse_parspec          <-rep(NA,m)
rmse_semiparspec      <-rep(NA,m)
rmse_smooth_spline    <-rep(NA,m)
rmse_smooth_2Dkernel  <-rep(NA,m)
ptm <- proc.time()
for (i in c(1:m)) {
  k=1  # First component of the bivariate time series
  kk=2 # Second component of the bivariate time series
  # qacf at each of the simulated time series
  qacf <- qacf.sim[[i]]
  #1.  Parametric coherence without smoothing
  qcper.par<-var2spec(p_max=p_max,qacf=qacf,order.form="aic",custom_order=NULL)$spec
  qcoh.par1<-Mod(qcper.par[k,kk,sel.f,sel.tau])^2/(Re(qcper.par[k,k,sel.f,sel.tau])*Re(qcper.par[kk,kk,sel.f,sel.tau]))
  qcoh.par1[qcoh.par1>1-eps]<-1-eps
  rmse_semiparspec[i]<-sum(Mod(qcoh.par1-tmp0)^2)/np
  
  #2.  Semi-parametric coherence with smoothing
  smoothed_qcoh.par1 <- Smoothing_across_quantiles(data=qcoh.par1,K=5,n.cores=detectCores()-2)$coh_smoothed
  rmse_parspec[i]<-sum(Mod(smoothed_qcoh.par1-tmp0)^2)/np
  
  # normalization for the computations in smoothing splines and 2D kernel smoothing
  qper12_a <- Mod(qcper[k,kk,,])
  qper12 <-spec.normalize(qper12_a)
  qper11_a <- Re(qcper[k,kk,,])
  qper11 <-spec.normalize(qper11_a)
  qper22_a <- Re(qcper[k,kk,,])
  qper22 <-spec.normalize(qper22_a)
  
  #3.Smoothing the raw coherence with smoothing splines across quantiles
  # and then across frequencies.
  smooth.qper12<-smooth_splines(qper12)
  smooth.qper11<-smooth_splines(qper11)
  smooth.qper22<-smooth_splines(qper22)
  qcoh.par2<-(smooth.qper12[sel.f,sel.tau])^2/(smooth.qper11[sel.f,sel.tau]*smooth.qper22[sel.f,sel.tau])
  qcoh.par2[qcoh.par2>1-eps]<-1-eps
  rmse_smooth_spline[i]<-sum(Mod(qcoh.par2-tmp0)^2)/np

  #4.Smoothing the raw coherence with 2D kernel smoothing
  KSsmooth.qper12<-smooth_2Dkernel(qper12)
  KSsmooth.qper11<-smooth_2Dkernel(qper11)
  KSsmooth.qper22<-smooth_2Dkernel(qper22)
  qcoh.par3<-(KSsmooth.qper12[sel.f,sel.tau])^2/(KSsmooth.qper11[sel.f,sel.tau]*KSsmooth.qper22[sel.f,sel.tau])
  qcoh.par3[qcoh.par3>1-eps]<-1-eps
  rmse_smooth_2Dkernel[i]<-sum(Mod(qcoh.par3-tmp0)^2)/np

}


proc.time()-ptm

#Results
#1. Parametric spectrum without smoothing
mean(sqrt(rmse_parspec))
sd(sqrt(rmse_parspec))

#2. Semi-parametric spectrum 
mean(sqrt(rmse_semiparspec))
sd(sqrt(rmse_semiparspec))

#3.Smoothing the raw coherence with smoothing splines across quantiles 
# and then across frequencies.
mean(sqrt(rmse_smooth_spline))
sd(sqrt(rmse_smooth_spline))

#4.Smoothing the raw coherence with 2D kernel smoothing 
mean(sqrt(rmse_smooth_2Dkernel))
sd(sqrt(rmse_smooth_2Dkernel))
