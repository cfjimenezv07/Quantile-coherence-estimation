# Function to compute the ordinary VAR coherence

var_coh<-function(x,order.max=10,aic=TRUE){
  #x: It's a bivariate time series. It's a matrix of two columns.
  k=1
  kk=2
  #x:= vector of two time series
  coeff    <-  ar(x=x, order.max = order.max ,aic = aic)
  ar.coeff <-  aperm(coeff$ar,c(2,3,1))
  var.pred <-  coeff$var.pred
  ar.order <-  coeff$order
  
  if(ar.order==0){
    nc=2
    nf=dim(x)[1]
    freq <- c(0:(nf - 1)) / nf
    varspec=array(data=NA,dim=c(nc,nc,nf))
    order_=ar.order
    varcoeff=ar.coeff
    Sigma=var.pred
    for (v in 1:nf ) {
      u=matrix(0,2,2)
      U=diag(nc)-u
      varspec[,,v]=Sigma
    }
    num<-Mod(varspec[k,kk,])^2
    den<-(Re(varspec[kk,kk,])*Re(varspec[k,k,]))
    VAR_coh<-num/den
  }else{
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
    num<-Mod(varspec[k,kk,])^2
    den<-(Re(varspec[kk,kk,])*Re(varspec[k,k,]))
    VAR_coh<-num/den
  }
  
  
  
  list(order=ar.order,var_coh=VAR_coh)
}
