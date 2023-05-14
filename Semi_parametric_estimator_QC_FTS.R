#Semiparametric estimator of quantile coherence

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
# Define frequencies and quantile levels
# ###############################################################################
#Fourier frequencies
nf<-5739 #2010-2019
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
# load the data. Financial time series
# ###############################################################################
stocks_names <- readRDS("./stocks_names.rds")
file<-"Financial_data.csv" # All 52 financial time series. Time period: 2010-2019
x0<-read.csv(file,header=T) #Column 1 is date, column 10 the S&P 500 index.
sp500<-x0[,10] # Financial index S&P500
Stocks<-x0[,-c(1,10)] #52 Financial stocks 
y2<-diff(log(sp500)) #Compute the logarithmic differences for the S&P financial index.
source("aux_semi_parametric_estimation_QC.R") # for the parametric estimator
source("joint_smoothing.R") #smoothing the parametric estimator of QC
# ###############################################################################
# Parametric estimator of the quantile coherence
# ###############################################################################

  p_max <- 10 # Maximum order to be selected.
  nc <- 2 #number of time series
  nq <- 93 # number of quantiles

# #####################################################################################
# Compute the Semi-parametric quantile coherence estimator for financial time series
# #####################################################################################
  # This function computes the semi-parametric estimator of the proposed QC
qaff <- function(w,order.form=c("aic","bic","aicc","None","custom_order"),custom_order=NULL,ncores=1){
    # ncores refers to the number of cores to be used for the qdft computation
    ##############################################################################
    # compute the qacf
    ##############################################################################
    y1<-diff(log(Stocks[,w])) #Compute the logarithmic differences for the stock w
    y=cbind(y1,y2)
    y.qdft <- qfa::qdft(y,tau,n.cores=ncores)
    qcper <- qfa::qdft2qper(y.qdft)
    qacf  <- qfa::qdft2qacf(y.qdft)
    #Compute the parametric estimator
    varspec1=var2spec(p_max=p_max,qacf=qacf,order.form=order.form,custom_order=custom_order)$spec
    k=1  # First component of the bivariate time series
    kk=2 # Second component of the bivariate time series
    
    #Compute the QC
    num<-Mod(varspec1[k,kk,sel.f,sel.tau])^2
    den<-(Re(varspec1[kk,kk,sel.f,sel.tau])*Re(varspec1[k,k,sel.f,sel.tau]))
    qcoh<-num/den
    eps<-1e-6
    qcoh[qcoh>1-eps]<-1-eps
    QCpar<-qcoh # parametric quantile coherence estimate
    #Smooth the parametric estimator across quantiles
    smoothed_varspec1 <- Smoothing_across_quantiles(data=QCpar,K=5)$coh_smoothed
    return(smoothed_varspec1)
}

# Estimate the QC maps for each of the 52 stocks in parallel

# Estimation with order selected based on AIC
A <- mclapply(1:dim(Stocks)[2], qaff, mc.preschedule = TRUE, mc.cores = 15,order.form="aic")
# Estimation with order equals 1
B <- mclapply(1:dim(Stocks)[2], qaff, mc.preschedule = TRUE, mc.cores = 15,order.form="custom_order",custom_order=1)
# Estimation with max.order=2  on the VAR model
p_max <- 2 # Maximum order to be selected.
C <- mclapply(1:dim(Stocks)[2], qaff, mc.preschedule = TRUE, mc.cores = 15,order.form="None")


################################################################################
# Get the plots for all the 52 stocks
################################################################################

Max=max(sapply(Smoothed_financial_stocks, max))
beta <- readRDS("./beta.rds")

# Set the directory for plots
dir.l<-"~/My Drive/Spring 2023/STAT 397/PhD project 2/Plots_paper/"
pdf(paste0(dir.l, "All_QC_maps.pdf"))

for (i in 1:length(Smoothed_financial_stocks)) {
  qper1=Smoothed_financial_stocks[[i]][sel.f,]
  tlab6<-paste("Quantile coherence:", stocks_names[i],",", "Beta=",beta[i] )
  qfa.plot(ff[sel.f],tau[sel.tau],qper1,rg.qper=range(c(0,Max)),rg.tau=range(tau[sel.tau]),color=color.qper,tlab=tlab6)
  
}

dev.off()




