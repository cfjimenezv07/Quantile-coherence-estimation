# Code that implements the joint smoothing splines. 



######################################################################################
# Smooth the spectrum functions
################################################################################
library(parallel)

# Reproducing kernel and matries for the cubic splines
R_kernel<-function(x,z){
  R=1/4*((z-1/2)^2-1/12)*((x-1/2)^2-12)-1/24*((abs(x-z)-1/2)^4-1/2*(abs(x-z)-1/2)^2+7/240)
  return(R)
}

#We plan to smooth the sequences of quantile coherence at each frequency
#Define the quantile levels,  the  independent variable x
tau.min <- 0.04
tau.max <- 0.96
tau.del<- 0.01
tau<-seq(tau.min,tau.max,tau.del)
ntau<-length(tau)
sel.tau<-which(tau >= 0.05 & tau <= 0.95)
x<-tau[sel.tau]

# The solution to the smoothing spline problem comes in the form 
#\hat{\beta}(\lambda)=(X^/topX+\lambdaS)^{-1}X^\topY. 
# Please see. Wood,Simon N..Generalized Additive Models: An Introduction with R, Second Edition.United Kingdom:CRC Press,2017.

#Compute the R matrix for the smoothing spline solution
matrix_R<-outer(x,x,FUN = R_kernel)
nq<- length(sel.tau) ##number of total quantiles
# Construction of the design matrix in the solution of the smoothing spline problem
ones<-rep(1,nq)
X<-cbind(ones,x,matrix_R) # Design matrix in the solution of smoothing spline problem
b_1<-matrix(0,nrow = (nq+nc),ncol = nc) # submatrices for the matrix S in the solution of the smoothing spline problem
b_2<-matrix(0,nrow=nc,ncol = nq) # submatrices for the matrix S in the solution of the smoothing spline problem
m_1<-rbind(b_2,matrix_R) # submatrices for the matrix S in the solution of the smoothing spline problem
S<-cbind(b_1,m_1) # Matrix S in the solution of the smoothing spline problem

#basis functions
s1<-psych::tr(t(X)%*%X)
s2<-psych::tr(S)
r=s1/s2 # Computation of the r parameter for the computational implementation in R. The solution to the
#smoothing spline problem depends on the hyperparameter lambda. 
# As with previous implementations we use r to compute the spar parameter. 
# Please see the documentation of the smooth.spline in R

#Implementation in C of the matrix operation needed for the solution of the smoothing spline problem.
library(Rcpp)
library(RcppArmadillo)
cppFunction(depends = "RcppArmadillo", 
            ' 
            arma::vec GLS_cpp( arma::mat X, arma::mat X_training, arma::mat training, arma::mat S, double Lambda) {
            arma::vec beta_hat = inv(X_training.t() * X_training + Lambda * S) * X_training.t() * training;
            return X * beta_hat;
            }
            ')
cppFunction(depends = "RcppArmadillo", 
            ' 
            arma::vec GLS_cpp2( arma::mat X, arma::mat Y, arma::mat S, double Lambda) {
            arma::vec beta_hat = inv(X.t() * X + Lambda * S) * X.t() * Y;
            return X * beta_hat;
            }
            ')
Smoothing_across_quantiles <- function(data,K){
  # data:  it should be a nf times nq matrix. nf: number of frequencies, nq: number of quantiles
  # K: the number of partitions for the Cross-validation
  Y <- data
  nfreq=dim(Y)[1]
  cross.v <- function(s){
    Lambda <- r * 256 ^ (3 * s - 1)
    ECM <- function(m){
      testing_index <-(as.integer((nq/K) * m) - as.integer((nq/K)-1)):as.integer((nq/K)*m)
      PE <- function(h){
        Y_h <- Y[h, ]
        training <- matrix(Y_h[-testing_index], ncol = 1)
        testing <- mean(Y_h[testing_index])
        X_training <- X[-testing_index, ]
        Y_hat <- GLS_cpp(X, X_training, training, S, Lambda)
        Y_hat <- mean(Y_hat[testing_index])
        return(sum((Y_hat - testing) ^ 2))
      }
      return(sum(sapply(1:nfreq, PE)))
    }
    return(sum(sapply(1:K, ECM)))
  }
  pe2 <- optimize(cross.v, c(-1.5, 1.5))
  sel_lambda <- r * 256 ^ (3 * pe2$minimum - 1) # see smooth.spline function in R
  vetsmo <- function(data){
    m <- dim(data)[1]
    n = dim(data)[2]
    x = matrix(0, m, n)
    for(i in 1:m){
      Y=as.matrix(data[i,])
      Lambda=  sel_lambda
      x[i,] =  t(GLS_cpp2(X, Y, S, Lambda))
    }
    return(x)
  }
  coh_smoothed <- vetsmo(data)
  return(list(coh_smoothed = coh_smoothed, sel_lambda = sel_lambda))
}
