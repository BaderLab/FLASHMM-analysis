##################################################
##Restricted maximum likelihood (REML)  for LMM estimation 
##Iterative algorithms for REML: 
##(1) Gradient methods: Newton–Raphson (NR), Fisher scoring (FS), and average information (AI)
##(2) Expectation-maximization (EM) algorithm, and 
##(3) Iterated MINQE (minimum norm quadratic estimation) 
##
##lmmfit: (reml.c v5.2)
##- REML algorithm working on columns of data matrix
##- Assuming the number of columns is not large
##- REML with FS
##
##The new version of lmmfit (lmmfit.nt)
##- No transpose of observed measurements, Y, a responses-by-samples matrix

##Inputs
##Y: a features-by-samples matrix of observed measurements, genes-by-cells matrix for scRNA-seq.
##X:  a design matrix for fixed effects, with row names identical to the column names of Y.
##Z = [Z1, ..., Zk],  a design matrix for different types (groups) of random effects
##  Zi, i=1,...,k, the design matrix for the i-th type (grouping) random effects
##  Every Zi is associated with a grouping factor. 
##  k: number of the types of the random effects
##d = (m1,...,mk), mi = ncol(Zi), number of columns in Zi
##  m1 + ... + mk = ncol(Z), number of columns in Z	
##theta0 = (s1, ...,sk, s_{k+1}), a vector of initial values of the variance components
##  si = sigma_i^2, the variance component of the i-th type random effects
##  s_{k+1} = sigma^2, the variance component of model residual error
##output.cov: If TRUE, lmm output the covariance matrices for the estimated coefficients, 
##  which are needed for testing contrasts.
##output.RE: If TRUE, output the best linear unbiased prediction (BLUP) of the random effects.

##Outputs
##theta: a matrix of the variance component estimates, 
##  each column corresponds to a sample and each row one variance component. 
##  The last row is the variance component of the residual error.
##se: standard errors of the estimated theta
##coef: a matrix of the fixed effects (coefficients)
##cov: a array of covariance matrices of the estimated coefficients (fixed effects)
##RE: a matrix of random effects
##################################################

lmmfit <- function(Y, X, Z, d, theta0 = NULL, method = "REML-FS", max.iter = 50, epsilon = 1e-5, output.cov = TRUE, output.RE = FALSE)
{
stopifnot(!any(is.na(Y)), !any(is.na(X)), !any(is.na(Z)))
stopifnot(ncol(Y) == nrow(X), ncol(Y) == nrow(Z))

##summary statistics
n <- nrow(X)
XY <- t(Y%*%X)
ZX <- t(Z)%*%X #as.matrix(t(Z)%*%X)
ZY <- t(Y%*%Z) #as.matrix(t(Z)%*%Y)
ZZ <- t(Z)%*%Z #as.matrix(t(Z)%*%Z)

#XXinv <- ginv(t(X)%*%X)
XXinv <- try(chol2inv(chol(t(X)%*%X)), silent = TRUE)
if (inherits(XXinv, "try-error")) XXinv <- ginv(t(X)%*%X)
Ynorm <- rowSums(Y*Y) #colSums(Y*Y)

rm(X, Y, Z)

##LMM fitting
lmm(XY, ZX, ZY, ZZ, XXinv = XXinv, Ynorm = Ynorm, n = n, d = d, theta0 = theta0, method = method, max.iter = max.iter, epsilon = epsilon, output.cov = output.cov, output.RE = output.RE)
}
