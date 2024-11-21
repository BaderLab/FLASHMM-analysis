##Variation (variance) explained by LMM
##
##Arguments:
##X:  a design matrix for fixed effects.
##Z = [Z1, ..., Zk],  a design matrix for different types (groups) of random effects
##  Zi, i=1,...,k, the design matrix for the i-th type (grouping) random effects
##  k: number of the types of the random effects
##d = (m1,...,mk), mi = ncol(Zi), number of columns in Zi
##lmm.fit: output of lmm or lmmfit for fitting LMM.
##beta: fixed effects, beta = lmm.fit$coef.
##theta: variance components, theta = lmm.fit$theta.
##  The last row of theta must be the variance components of residual errors.
##covariates: character vector for partitioning (groupping) fixed effects. 
##  Covariates[i] has a design matrix Xi, i=1,...,L.
##  X = [X1,...,XL].
##
##Value
##total: the variations explained by all fixed effects (X), all random effects (Z), and residual errors (R).
##Vxs: sample variances for covariates X
##Vx: the variations explained by covariates based on orthogonal decomposition.
##Vz: the variations explained by random effects.
##fxs: fraction of variation explained by covariates based on proportion of the sample variances.
##fx: fraction of variation explained by covariates based on orthogonal decomposition.
##fz: fraction of variation explained by random effects.
##fr: fraction of variation explained by residual errors.

varexp <- function(X, Z, d, lmm.fit, beta, theta, covariates = colnames(X), output.fraction = TRUE)
{
if (!missing(lmm.fit)){
	beta <- lmm.fit$coef
	theta <- lmm.fit$theta
} else {
	if (missing(beta) | missing(theta)) stop("Miss beta and/or theta.")
}
beta <- as.matrix(beta)
theta <- as.matrix(theta)
X <- as.matrix(X)
Z <- as.matrix(Z)

##variance components = 0 if negative
theta <- theta*(theta > 0)
n <- nrow(X)

stopifnot(ncol(X) == nrow(beta), ncol(Z) == sum(d), 
	length(d) == nrow(theta) - 1, ncol(beta) == ncol(theta))

##centering: PnX, PnZ
X <- sweep(X, MARGIN = 2, STATS = colMeans(X), FUN = "-")
Z <- sweep(Z, MARGIN = 2, STATS = colMeans(Z), FUN = "-")

##remove intercept
j <- (colSums(abs(X)) < 1e-10)
nocovariate <- (sum(j) >= ncol(X))
if (any(j)){
	covariates <- setdiff(covariates, colnames(X)[j])
	X <- X[, !j, drop = FALSE]
	beta <- beta[!j, , drop = FALSE]
}

##variances explained by residual errors
if (rownames(theta)[length(d)+1] != "var0") message("Warning: the last row of theta should be variance components of residual errors.")
Vr <- theta[length(d) + 1, ]

##variances explained by random effects
Vz <- matrix(NA, nrow = ncol(theta), ncol = nrow(theta) - 1)
rownames(Vz) <- colnames(theta)
if (is.null(names(d))){
	if (length(d) == 1) colnames(Vz) <- "Z"
	else colnames(Vz) <- paste0("Z", 1:length(d))
} else colnames(Vz) <- names(d)

d0 <- 0
for (i in 1:length(d)){
	j <- (d0+1):(d0+d[i])
	zz <- sum(Z[,j]^2)/(n-1)
	Vz[,i] <- zz*theta[i,]
	d0 <- d[i]
}
Vztot <- rowSums(Vz)

##variances explained by fixed effects
if (nocovariate) {
	Vxtot <- rep(0, ncol(beta))
	Vx <- matrix(0, nrow = ncol(beta), ncol = 1)
	rownames(Vx) <- colnames(beta)
	colnames(Vx) <- "intercept"	
	Vxsv <- Vx
} else {
##For all fixed effects
Vxtot <- colSums((X%*%beta)^2)/(n-1)

##sample variance
Vxsv <- matrix(NA, nrow = ncol(beta), ncol = length(covariates))
rownames(Vxsv) <- colnames(beta)
colnames(Vxsv) <- covariates
for (i in 1:length(covariates)){
	#j <- grep(covariates[i], colnames(X))
	j <- which(colnames(X) == covariates[i])
	Vxsv[, i] <- colSums((X[, j]%*%beta[j, , drop = FALSE])^2)/(n-1)
	}
	
##For each of fixed effects (covariates)
Vx <- matrix(NA, nrow = ncol(beta), ncol = length(covariates))
	rownames(Vx) <- colnames(beta)
	colnames(Vx) <- covariates

for (i in 1:length(covariates)){
	#j <- grep(covariates[i], colnames(X))
	j <- which(colnames(X) == covariates[i])
	if (length(j) == 0) stop(paste0(covariates[i], " is not in colnames(X)."))
	XXi <- ginv(t(X[, j])%*%X[, j])
	HXbi <- X[, j]%*%(XXi%*%((t(X[, j])%*%X))%*%beta)
	Vx[, i] <- colSums((HXbi)^2)/(n-1)		
	}
}

##fraction of variation explained
if (output.fraction){
	##total variance
	Vt <- Vxtot + Vztot + Vr
	##fractions
	fz <- Vz/Vt
	fr <- Vr/Vt
	##orthogonal decomposition
	fx <- Vx/Vt
	##proportional partitioning
	fxsv <- (Vxsv/rowSums(Vxsv))*(Vxtot/Vt)
} else {
	fxsv <- NULL
	fx <- NULL
	fz <- NULL
	fr <- NULL
}
	
list(total = data.frame(X = Vxtot, Z = Vztot, R = Vr), Vxs = Vxsv, Vx = Vx, Vz = Vz, fxs = fxsv, fx = fx, fz = fz, fr = fr)
}