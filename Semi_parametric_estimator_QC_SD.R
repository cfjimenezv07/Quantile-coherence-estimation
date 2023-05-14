#Semiparametric estimator of quantile coherence for simulated data

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

# Load auxiliary files 
# source("qfa_lib_qdft_v2.txt") #Necessary functions for the qcper and qacf computations
source("qfa_lib.txt") #Necessary functions for simulated time series.

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

################################################################################
##Simulated data
################################################################################
# As an example consider one realization of the mixture model 3 in the paper.
n<-512
m<-1
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
# compute the qacf
# ###############################################################################
y<-cbind(xx,xx2)
y.qdft <- qfa::qdft(y,tau,n.cores=12)
qcper <- qfa::qdft2qper(y.qdft)
qacf  <- qfa::qdft2qacf(y.qdft)

# ###############################################################################
# Parametric estimator of the quantile coherence
# ###############################################################################

p_max <- 10 # Maximum order to be selected.
nc <- 2 #number of time series
nq <- 93 # number of quantiles
source("aux_semi_parametric_estimation_QC.R")
varspec1=var2spec(p_max=p_max,qacf=qacf,order.form="aic",custom_order=NULL)$spec

# ###############################################################################
# Define the components of the bivariate time series
# ###############################################################################
k=1  # First component of the bivariate time series
kk=2 # Second component of the bivariate time series

#Plot the parametric estimator.
num<-Mod(varspec1[k,kk,sel.f,sel.tau])^2
den<-(Re(varspec1[kk,kk,sel.f,sel.tau])*Re(varspec1[k,k,sel.f,sel.tau]))
qcoh<-num/den
eps<-1e-6
qcoh[qcoh>1-eps]<-1-eps
QCpar<-qcoh # parametric quantile coherence estimate
lab<-c("Series 1","Series 2")
tlab6<-paste0("Quantile Coherence: parametric spectrum\n(",lab[k],", ",lab[kk],"): ")
qfa::qfa.plot(ff[sel.f],tau[sel.tau],QCpar,rg.qper=range(c(0,QCpar)),rg.tau=range(tau[sel.tau]),color=color.qper,tlab=tlab6)


# ###############################################################################
# Semi-parametric estimator for the quantile coherence
# ###############################################################################
source("joint_smoothing.R")
smoothed_varspec1 <- Smoothing_across_quantiles(data=QCper,K=5)$coh_smoothed

#Plot the semi-parametric estimator.
smoothed_QCpar<-smoothed_varspec1 # Semi-parametric quantile coherence estimate
lab<-c("Series 1","Series 2")
tlab6<-paste0("Quantile Coherence: Semi-parametric spectrum\n(",lab[k],", ",lab[kk],"): ")
qfa::qfa.plot(ff[sel.f],tau[sel.tau],smoothed_QCpar,rg.qper=range(c(0,smoothed_QCpar)),rg.tau=range(tau[sel.tau]),color=color.qper,tlab=tlab6)


