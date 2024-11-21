##################################################
##Restricted maximum likelihood (REML)  for LMM estimation 
##Iterative algorithms for REML: 
##(1) Gradient methods: Newton–Raphson (NR), Fisher scoring (FS), and average information (AI)
##(2) Expectation-maximization (EM) algorithm, and 
##(3) Iterated MINQE (minimum norm quadratic estimation) 
##
##reml.c: (v5.2)
##- REML algorithm working on columns of data matrix
##- Assuming the number of columns is not large
##- REML with FS
##
##lmm was named as lmmfitSS before.
##lmmfitSS: using correlation-related summary statistics (SS): 
##XY, ZX and ZY, as inputs to fit LMM.

##Inputs
##XXinv = (X'X)^{-1}
##ZZ: ZZ <- t(Z)%*%Z
##Ynorm: Ynorm <- colSums(Y*Y) 
##XXinv: XXinv <- ginv(t(X)%*%X)
##XY: XY <- t(X)%*%Y
##ZX: ZX <- t(Z)%*%X
##ZY: ZY <- t(Z)%*%Y
#n <- nrow(Y) 
##d = (m1,...,mk), mi = ncol(Zi), number of columns in Zi
##  m1 + ... + mk = ncol(Z), number of columns in Z	
##sigma2 = (s1, ...,sk, s_{k+1}), a vector of initial values of the variance components
##  si = sigma_i^2, the variance component of the i-th type random effects
##  s_{k+1} = sigma^2, the variance component of model residual error
##output.cov: If TRUE, lmm output the covariance matrices for the estimated coefficients, 
##  which are needed for testing contrasts.

##Outputs
##theta: a matrix of the variance component estimates, 
##  each column corresponds to a sample and each row one variance component. 
##  The last row is the variance component of the residual error.
##se: standard errors of the estimated theta
##coef: a matrix of the fixed effects (coefficients)
##cov: a array of covariance matrices of the estimated coefficients (fixed effects)
##################################################

lmm <- function(XY, ZX, ZY, ZZ, XXinv, Ynorm, n, d, sigma2 = NULL, 
		method = "REML-FS", max.iter = 50, epsilon = 1e-5, output.cov = TRUE)
{
stopifnot(!any(is.na(XY)), !any(is.na(ZX)), !any(is.na(ZY)))
p <- ncol(ZX)
k <- length(d)  

stopifnot(sum(d) == ncol(ZZ))

##xxz = (X'X)^{-1}X'Z
##zrz = Z'RZ
##zry = Z'Ry
##yry = [y1'Ry1,...,ym'Rym]
xxz <- XXinv%*%t(ZX)
zrz <- ZZ - ZX%*%(XXinv%*%t(ZX))
zry <- ZY - ZX%*%(XXinv%*%XY)
yry <- Ynorm - colSums(XY*(XXinv%*%XY))


##########
##estimated variance components in iterations
##dl: first derivatives of log likelihood

niter <- NULL ##number of iterations
dlogL <- NULL ##derivatives of log-likelihoods at the last iteration
theta <- matrix(nrow = k + 1, ncol = ncol(XY), dimname = list(c(1:k, 0), colnames(XY)))
se.theta <- theta
beta <- matrix(nrow = nrow(XY), ncol = ncol(XY), dimnames = dimnames(XY))
covbeta <- array(dim = c(nrow(XY), nrow(XY), ncol(XY)), 
	dimnames = list(rownames(XY), rownames(XY), colnames(XY)))
tval <- beta

for (jy in 1:ncol(ZY)) {
	if (is.null(sigma2)) {
		s <- c(rep(0, k), yry[jy]/(n-p))
	} else s <- sigma2

dl <- 100
iter <- 0
while ((max(abs(dl)) > epsilon)	& (iter < max.iter)){
	iter <- iter + 1
	
	fs <- matrix(NA, k+1, k+1)	##Fisher scoring matrix
	dl <- rep(NA, k+1) ##dl: derivatives of log-likelihood
		
	sr <- s[1:k]/s[k+1]
	##M = (SZ'RZ + I)^{-1}
	M <- solve(sweep(zrz, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
	ZRZ <- zrz%*%M
	ZR2Z <- ZRZ%*%M
	yRZ <- t(zry[, jy])%*%M
	
	mi <- 0
	for (i in 1:k){	
		ik <- (mi+1):(mi+d[i])
		dl[i] <- (sum((yRZ[ik])^2)/s[k+1]^2 - sum(diag(ZRZ[ik, ik, drop = FALSE]))/s[k+1])/2

	mj <- 0
	for (j in 1:i){
		ji <- (mj+1):(mj+d[j])
		fs[i, j] <- sum((ZRZ[ji, ik])^2)/s[k+1]^2/2
		fs[j, i] <- fs[i, j]
		mj <- mj + d[j]
		}
		
	j <- k+1		
	fs[i, j] <- sum(diag(ZR2Z[ik, ik, drop = FALSE]))/s[k+1]^2/2
	fs[j, i] <- fs[i, j]
	mi <- mi + d[i]
	}
	
	i <- k+1
	fs[i, i] <- (n - p - sum(d) + sum(t(M)*M))/s[k+1]^2/2
	
	yR2y <- yry[jy] - sum(((t(M) + diag(sum(d)))%*%zry[, jy])*(M%*%(rep(sr, times = d)*zry[, jy])))
	dl[i] <-  (yR2y/s[k+1]^2 - (n-p-sum(d)+sum(diag(M)))/s[k+1])/2

	##The H, FS, and AI matrices can be singular.
	#Minv <- solve(M)
	Minv <- ginv(fs)
	s <- s + Minv%*%dl
	}
	

if (max(abs(dl)) > epsilon) {
	warningText <- paste0("The first derivatives of log likelihood for Y", jy)
	dlText <- paste0(ifelse(abs(dl) > 1e-3, round(dl, 4), format(dl, digits = 3, scientific = TRUE)), collapse = ", ")
	warning(paste0(warningText, ": ", dlText, ", doesn't reach the zero, epsilon ", epsilon))
	}

##
sr <- s[1:k]/s[k+1]
M <- solve(sweep(ZZ, 1, STATS = rep(sr, times = d), FUN = "*") + diag(sum(d)))
M <- sweep(M, 2, STATS = rep(sr, times = d), FUN = "*") 
xvx <- XXinv + xxz%*%(ginv(diag(sum(d)) - M%*%(ZX%*%xxz))%*%(M%*%t(xxz)))
xvy <- XY[, jy] - t(ZX)%*%(M%*%ZY[, jy])
b <- xvx%*%xvy ##beta, fixed effects
covbeta[,,jy] <- (xvx + t(xvx))*(s[k+1]/2)

##outputs
niter <- c(niter, iter)
theta[, jy] <- s 
se.theta[, jy] <- sqrt(diag(Minv))
beta[, jy] <- b 
dlogL <- cbind(dlogL, dl)
tval[, jy] <- b/sqrt(diag(as.matrix(covbeta[,,jy])))
}
pval <- 2 * pt(-abs(tval), df = n-p)

if (!output.cov) covbeta <- NULL

list(method = method, dlogL = dlogL, niter = niter, coef = beta, t = tval, p = pval, 
	cov = covbeta, df = n-p, theta = theta, se = se.theta)
}
