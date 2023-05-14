# Auxiliary functions required for the separametric estimator of quantile coherence

# load the required packages
packages <- c("QZ", "matrixcalc", "quantreg","fields","colorRamps", "cluster", 
              "factoextra","graphics","AR","TSA","MTS","tsDyn")
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
# Durbin-Levinson algorithm for the estimation of the VAR coefficients
# ###############################################################################

# Function for the delta computation in the Durbin-Levinson Algorithm

Delta_n_func <- function (n,phi,qacf_fix_quantile) {
  #n-> order of the model
  partial <- qacf_fix_quantile[, , (n + 2)]
  for (i in (n - 1):0) {
    partial <- partial - (phi[, , (n - i)] %*% qacf_fix_quantile[, , (i + 2)])
  }
  partial
}

# ###############################################################################
# For a particular quantile compute the VAR coefficients,
#and the aic, bic,aicc for order selection if desired.
# ###############################################################################

compute_coeff <- function (p_max,qacf_fix_quantile) {
  
  # qacf_fix_quantile:  it is the qacf at a particular quantile level.
  nc  <- dim(qacf_fix_quantile)[1] 
  nn  <- dim(qacf_fix_quantile)[3]/2
  n   <- dim(qacf_fix_quantile)[3]
  V           <- array(data = 0, dim = c(nc, nc, (p_max + 1)))
  V_tilde     <- array(data = 0, dim = c(nc, nc, (p_max + 1)))
  Delta       <- array(data = 0, dim = c(nc, nc, (p_max + 1)))
  phi         <- array(data = 0, dim = c(nc, nc, p_max))
  phi_tilde   <- array(data = 0, dim = c(nc, nc,p_max))
  Psi         <- array(data = 0, dim = c(nc, nc, p_max))
  
  V[, , 1]           <- qacf_fix_quantile[, , 1] 
  V_tilde[, , 1]     <- qacf_fix_quantile[, , 1]
  Delta[, , 1]       <- qacf_fix_quantile[, , 2]
  
  aic=rep(NA,p_max+1)
  bic=rep(NA,p_max+1)
  aicc=rep(NA,p_max+1)
  
  aic[1]=n*log(det(V[,,1]))
  bic[1]=n*log(det(V[,,1]))
  aicc[1]=n*log(det(V[,,1])) 
  
  for (i in 1:p_max) {
    V_sr             <- MTS::msqrt(V[, ,i])$mtxsqrt
    V_sr_tilde       <- MTS::msqrt(V_tilde[, ,i])$mtxsqrt
    Psi[, , i]       <- solve(V_sr)%*%Delta[, , i]%*%solve(V_sr_tilde)
    phi[, , i]       <- V_sr%*%Psi[, , i]%*%solve(V_sr_tilde)
    phi_tilde[, , i] <- V_sr_tilde%*%t(Psi[, , i])%*%solve(V_sr)
    
    
    if (i > 1) {
      for(j in 1:(i-1)){
        phi[,,j] <- phi_temp[,,j] - phi[, , i] %*% phi_tilde_temp[,,(i-j)]
        phi_tilde[,,j] <- phi_tilde_temp[,,j] - phi_tilde[, , i] %*% phi_temp[,,(i-j)]
      }
    }
    
    
    V[, , (i + 1)]           <- V_sr%*%(diag(nc)-Psi[, , i]%*%t(Psi[, , i]))%*%t(V_sr)
    V_tilde[, , (i + 1)]     <- V_sr_tilde%*%(diag(nc)-t(Psi[, , i])%*%Psi[, , i])%*%t(V_sr_tilde)
    
    #Assure symmetry  of the V matrix
    
    if(!isSymmetric(V[, , (i + 1)])){
      V[, , (i + 1)]=(V[, , (i + 1)]+t(V[, , (i + 1)]))/2
    }
    
    if(!isSymmetric(V_tilde[, , (i + 1)])){
      V_tilde[, , (i + 1)]=(V_tilde[, , (i + 1)]+t(V_tilde[, , (i + 1)]))/2
    }
    
    Delta[, , (i + 1)]       <- Delta_n_func(i, phi,qacf_fix_quantile)
    
    
    aic[i+1]=n*log(abs(det(V[,,(i+1)])))+2*nc^2*i
    bic[i+1]=n*log(abs(det(V[,,(i+1)])))+log(n)*i*nc^2
    aicc[i+1]=n*log(abs(det(V[,,(i+1)])))+(2*i*n*nc^2)/(n-i*nc^2-1)
    
    phi_temp         <- array(data = 0, dim = c(nc, nc, i))
    phi_tilde_temp   <- array(data = 0, dim = c(nc, nc, i))
    phi_temp[,,1:i]         <- phi[,,1:i]
    phi_tilde_temp[,,1:i]   <- phi_tilde[,,1:i]
  
    
  }
  
  
  list(V = V[,,p_max+1], phi = phi,aic=aic,bic=bic,aicc=aicc,
       V_tilde = V_tilde[,,p_max+1],phi_tilde = phi_tilde,all_V=V
       ,all_delta=Delta,all_V_tilde=V_tilde,Psi=Psi)
}


# ###############################################################################
# Compute the VAR coefficients for all the quantile levels
# ###############################################################################

# The user can choose how to perform the order of the VAR model
# AIC, BIC, AICC are automatic order selection criteria
# "None" uses the maximum order as the selected order
# custom_order fits a VAR model with the desired customized order

qacf2VAR <- function (p_max,qacf,order.form=c("aic","bic","aicc","None","custom_order"),custom_order=NULL) {
  nc=dim(qacf)[1] # from the qcper. Quantile Fourier cross-periodogram
  nq=dim(qacf)[4] # from the qcper. Quantile Fourier cross-periodogram
  
  if(order.form=="None"){
    phi         <- array(data = 0, dim = c(nc, nc, p_max,nq)) 
    V           <- array(data = 0, dim = c(nc, nc,nq))
    for (k in 1:nq) {
      result_coeff <- compute_coeff(p_max,qacf_fix_quantile=qacf[,,,k])
      V[,,k] <- result_coeff$V
      phi[,,,k] <- result_coeff$phi
    }
    selec_order=p_max
    aicm<-NULL
    bicm<-NULL
    aiccm<-NULL
  }
   else if(order.form=="custom_order"){
     phi         <- array(data = 0, dim = c(nc, nc, p_max,nq)) 
     V           <- array(data = 0, dim = c(nc, nc,nq))
     for (k in 1:nq) {
       result_coeff <- compute_coeff(p_max,qacf_fix_quantile=qacf[,,,k])
       V[,,k] <- result_coeff$V
       phi[,,,k] <- result_coeff$phi
     }
     selec_order=custom_order
     aicm<-NULL
     bicm<-NULL
     aiccm<-NULL
     
   }else{
    aic<-matrix(NA,nrow=p_max+1,ncol=nq)
    bic<-matrix(NA,nrow=p_max+1,ncol=nq)
    aicc<-matrix(NA,nrow=p_max+1,ncol=nq)
    for (k in 1:nq) {
      result_coeff <-compute_coeff(p_max=p_max,qacf_fix_quantile=qacf[,,,k]) 
      aic[,k]=result_coeff$aic
      bic[,k]=result_coeff$bic
      aicc[,k]=result_coeff$aicc
    }
    aicm<-apply(aic,1,mean,na.rm=T)
    bicm<-apply(bic,1,mean,na.rm=T)
    aiccm<-apply(aicc,1,mean,na.rm=T)
    if(order.form=="aic") selec_order=which(aicm==min(aicm))
    if(order.form=="bic") selec_order=which(bicm==min(bicm)) 
    if(order.form=="aicc") selec_order=which(aiccm==min(aiccm))
    selec_order=(selec_order-1)
    phi         <- array(data = 0, dim = c(nc, nc, selec_order,nq)) 
    V           <- array(data = 0, dim = c(nc, nc,nq))
    
    #Compute the VAR coefficients for all the quantile levels
    
    for (k in 1:nq) {
      result_coeff <- compute_coeff(p_max=selec_order,qacf_fix_quantile=qacf[,,,k])
      V[,,k] <- result_coeff$V
      phi[,,,k] <- result_coeff$phi
    }
    
    
  }
  
  
  list(phi = phi,V=V,selec_order=selec_order, aic=aicm,bic=bicm,aicc=aiccm)
}


# ###############################################################################
# Compute the parametric spectrum coming from the VAR coefficients
# ###############################################################################

var2spec<-function(p_max,qacf,order.form=c("aic","bic","aicc","None","custom_order"),custom_order=NULL){
  
  nc=dim(qacf)[2] # from the qcper. Quantile Fourier cross-periodogram
  nq=dim(qacf)[4] # from the qcper. Quantile Fourier cross-periodogram
  nf=dim(qacf)[3] # from the qcper. Quantile Fourier cross-periodogram
  
  freq <- ff # Fourier frequencies.
  
  varspec=array(data=NA,dim=c(nc,nc,nf,nq))
  varspecfreq=array(data=NA,dim=c(nc,nc,nf))
  if(order.form=="None"){
    result=qacf2VAR(p_max=p_max,qacf=qacf,order.form="None",custom_order=NULL)
    for (k in 1:nq) {
      varcoeff=result$phi[,,,k]
      Sigma=result$V[,,k]
      for (v in 1:nf ) {
        u=matrix(0,nc,nc)
        for (r in 1:p_max) {
          u=u+varcoeff[,,r]*base::complex(real = cos(2*pi*r*freq[v]),imaginary =-sin(2*pi*r*freq[v]))
        }
        U=diag(nc)-u
        varspecfreq[,,v]=solve(U)%*%Sigma%*%QZ::H(solve(U))
      }
      varspec[,,,k]=varspecfreq
    }
    
  } else if(order.form=="custom_order"){
    result=qacf2VAR(p_max=p_max,qacf=qacf,order.form="custom_order",custom_order=custom_order)
    order_=result$selec_order
    if(order_==0){
      for (k in 1:nq) {
        Sigma= result$V[,,k]
        varspec[,,,k]=Sigma
      }
    }else {
      for (k in 1:nq) {
        varcoeff=result$phi[,,,k]
        Sigma=result$V[,,k]
        for (v in 1:nf ) {
          u=matrix(0,nc,nc)
          for (r in 1:order_) {
            u=u+varcoeff[,,r]*complex(real = cos(2*pi*r*freq[v]),imaginary =-sin(2*pi*r*freq[v]))
          }
          U=diag(nc)-u
          varspecfreq[,,v]=solve(U)%*%Sigma%*%QZ::H(solve(U))
        }
        varspec[,,,k]=varspecfreq
      }
    }
    
  } else if(order.form=="aic"){
    result=qacf2VAR(p_max=p_max,qacf=qacf,order.form="aic",custom_order=NULL)
    order_=result$selec_order
    if(order_==0){
      for (k in 1:nq) {
        Sigma= result$V[,,k]
        varspec[,,,k]=Sigma
      }
    }else {
      for (k in 1:nq) {
        varcoeff=result$phi[,,,k]
        Sigma=result$V[,,k]
        for (v in 1:nf ) {
          u=matrix(0,nc,nc)
          for (r in 1:order_) {
            u=u+varcoeff[,,r]*complex(real = cos(2*pi*r*freq[v]),imaginary =-sin(2*pi*r*freq[v]))
          }
          U=diag(nc)-u
          varspecfreq[,,v]=solve(U)%*%Sigma%*%QZ::H(solve(U))
        }
        varspec[,,,k]=varspecfreq
      }
    }
    
  
  }else if(order.form=="bic"){
    result=qacf2VAR(p_max=p_max,qacf=qacf,order.form="bic",custom_order=NULL)
    order_=result$selec_order
    if(order_==0){
      for (k in 1:nq) {
        Sigma= result$V[,,k]
        varspec[,,,k]=Sigma
      }
    }else {
      for (k in 1:nq) {
        varcoeff=result$phi[,,,k]
        Sigma=result$V[,,k]
        for (v in 1:nf ) {
          u=matrix(0,nc,nc)
          for (r in 1:order_) {
            u=u+varcoeff[,,r]*complex(real = cos(2*pi*r*freq[v]),imaginary =-sin(2*pi*r*freq[v]))
          }
          U=diag(nc)-u
          varspecfreq[,,v]=solve(U)%*%Sigma%*%QZ::H(solve(U))
        }
        varspec[,,,k]=varspecfreq
      }
    }
  }else{
    result=qacf2VAR(p_max=p_max,qacf=qacf,order.form="aicc",custom_order=NULL)
    order_=result$selec_order
    if(order_==0){
      for (k in 1:nq) {
        Sigma= result$V[,,k]
        varspec[,,,k]=Sigma
      }
    }else {
      for (k in 1:nq) {
        varcoeff=result$phi[,,,k]
        Sigma=result$V[,,k]
        for (v in 1:nf ) {
          u=matrix(0,nc,nc)
          for (r in 1:order_) {
            u=u+varcoeff[,,r]*complex(real = cos(2*pi*r*freq[v]),imaginary =-sin(2*pi*r*freq[v]))
          }
          U=diag(nc)-u
          varspecfreq[,,v]=solve(U)%*%Sigma%*%QZ::H(solve(U))
        }
        varspec[,,,k]=varspecfreq
      }
    }
  }
  
  list(spec=varspec)
}


