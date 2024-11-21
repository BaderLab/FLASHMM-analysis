	
##################################################
##Simulation analysis I
##null hypothesis of no treatment effect. 
##H0: beta = 0

##########	
##Simulating data
##trt: treatments
##cls: sample clusters
##Y: responses
##model: y = trt*beta + cls*b + e

##null hypothesis of no treatment effect. 
##random effects b ~ Norm(0, sigma^2_b)
##Null model: beta = 0

##variance components
sigma2e <- 1    ##for residual errors
sigma2b <- 0.16 ##for random effects

##number of treatments
ntrt <- 2
##number of clusters
ncls <- 16
##number of simulations
nSimu <- 1000 
##number of fixed effects != 0, 
##the first 100 out of 1000 simulations
#nFixedEffects <- 100 
nFixedEffects <- 0
sizeSamples <- c(5e3, 1e4, 1.5e4, 2e4, 2.5e4, 3e4)

##########
##Run time for LMM fitting only
##not including LMM test
##
##lmer
##lme4::lmer for LMM fitting
##lmerTest::lmer for LMM fixed effect testing 
##(lme4 doesn't provide LMM testing)

##different sample sizes
#nSamples <- 10000
for (nSamples in sizeSamples){
	
set.seed(231212)

##sample clusters and treatments
dat <- data.frame(
	cls = as.factor(sample.int(ncls, nSamples, replace = T)),
	trt = LETTERS[sample.int(ntrt, nSamples, replace = T)])
	
##responses
Y <- matrix(NA, nrow = nrow(dat), ncol = nSimu)
colnames(Y) <- paste0("Y", 1:nSimu)
for (j in 1:ncol(Y)) {
for (cls in unique(dat$cls)){
	index <- (dat$cls == cls)
	Y[index, j] <- rnorm(1, sd = sqrt(sigma2b)) + rnorm(sum(index), sd = sqrt(sigma2e))
	}##cls
	}##i

##treatment effect
#indexB <- (dat$trt == "B")
#for (j in 1:nFixedEffects){
#	betaj <- runif(1, 0.25, 1)
#	#Y[indexB, j] <- Y[indexB, j] + rnorm(sum(indexB), mean = betaj, sd = 0.1)
#	}
	
##########
##LMM fit and test

##outputs
##lmmfit
felmm <- NULL ##fe: fixed effects
tvlmm <- NULL ##tv: t-value
plmm <- NULL ##p: p-value
slmm <- NULL ##s: sigma^2
rtlmm <- NULL ##rt: running time

##lmer with default setting
felme <- NULL
tvlme <- NULL
plme <- NULL
slme <- NULL
rtlme <- NULL

##lmer with control changed
felme2 <- NULL
tvlme2 <- NULL
plme2 <- NULL
slme2 <- NULL
rtlme2 <- NULL


##
for (j in 1:ncol(Y)) {
	
#####	
##lmmfit
##design matrix for fixed effects
X <- model.matrix(~ trt, data = dat)
##design matrix for random effects
Z <- model.matrix(~ 0 + cls, data = dat)
d <- ncol(Z)

##LMM fit
t1 <- Sys.time()
fit <- lmmfit(Y = Y[, j, drop = F], X, Z, d=d, epsilon = 1e-8)
	t2 <- Sys.time()
	#rtlmm <- rtlmm + difftime(t2, t1, units = "secs")
	rtlmm <- c(rtlmm, difftime(t2, t1, units = "secs"))
	rm(t1, t2)

test <- lmmtest(fit)

felmm <- cbind(felmm, fit$coef)
tvlmm <- cbind(tvlmm, test[, grep("_t", colnames(test))])
plmm <- cbind(plmm, test[, grep("_pvalue", colnames(test))])
slmm <- cbind(slmm, fit$theta)

rm(fit)

#####
##lmer
##the observations
dat$y <- Y[, j]

##lmer with default setting
##LMM fitting: lme4::lmer 
t1 <- Sys.time()
lmerFit <- lme4::lmer(y ~ trt + (1 | cls), data = dat)
	t2 <- Sys.time()
	#rtlme <- rtlme + difftime(t2, t1, units = "secs") 
	rtlme <- c(rtlme, difftime(t2, t1, units = "secs"))
	rm(t1, t2)

sfit <- summary(lmerFit)$coefficients
felme <- cbind(felme, sfit[, "Estimate"])
tvlme <- cbind(tvlme, sfit[, "t value"])
slme <- cbind(slme, as.data.frame(VarCorr(lmerFit))$vcov)

##LMM testing: lmerTest::lmer
lmerFit <- lmerTest::lmer(y ~ trt + (1 | cls), data = dat)
sfit <- summary(lmerFit)$coefficients
plme <- cbind(plme, sfit[, "Pr(>|t|)"])

rm(sfit, lmerFit)

#summary(lmerFit)
#summary(lmerFit, ddf = "Satterthwaite")
#summary(lmerFit, ddf = "Kenward-Roger") #taking a long time!

#####
##lmer with control changed
control <- lmerControl(
		check.conv.grad = .makeCC("warning", tol = 2e-5, relTol = NULL),
		optCtrl=list(xtol_rel=1e-8, xtol_abs=1e-11, ftol_abs=1e-11))

##LMM fitting: lme4::lmer 
t1 <- Sys.time()
lmerFit <- lme4::lmer(y ~ trt + (1 | cls), data = dat, control = control)
	t2 <- Sys.time()
	#rtlme <- rtlme + difftime(t2, t1, units = "secs") 
	rtlme2 <- c(rtlme2, difftime(t2, t1, units = "secs"))
	rm(t1, t2)

sfit <- summary(lmerFit)$coefficients
felme2 <- cbind(felme2, sfit[, "Estimate"])
tvlme2 <- cbind(tvlme2, sfit[, "t value"])
slme2 <- cbind(slme2, as.data.frame(VarCorr(lmerFit))$vcov)

##LMM testing: lmerTest::lmer
lmerFit <- lmerTest::lmer(y ~ trt + (1 | cls), data = dat, control = control)
sfit <- summary(lmerFit)$coefficients
plme2 <- cbind(plme2, sfit[, "Pr(>|t|)"])

##VarCorr: variance = (std.dev)^2!!!
rm(sfit, lmerFit)	

}##nSimu


#####
##fit LMM using all Y
t1 <- Sys.time()
fit <- lmmfit(Y, X, Z, d=d, epsilon = 1e-8)
	t2 <- Sys.time()
	rtlmm_mult <- difftime(t2, t1, units = "secs")
	
	rm(t1, t2)

test <- lmmtest(fit)
	
range(fit$niter)
range(fit$dlogL)

##compare with the single column outputs
cat(paste0("nSamples = ", nSamples, "\n"))
print(range(fit$coef - felmm))	
print(range(fit$theta - slmm))	
print(range(test[, grep("_t", colnames(test))] - t(tvlmm)))
print(range(test[, grep("_pvalue", colnames(test))] - t(plmm)))

#####
rtlmm <- as.difftime(rtlmm, units = "secs")
rtlme <- as.difftime(rtlme, units = "secs")

save(felmm, tvlmm, plmm, slmm, rtlmm, felme, tvlme, plme, slme, rtlme, rtlmm_mult, 
	felme2, tvlme2, plme2, slme2, rtlme2, 
	nSimu, nFixedEffects, nSamples, sigma2e, sigma2b, 
	file = paste0("simu2_beta_lme4_n", nSamples, ".RData"))

}##nSamples

Warning messages:
1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge with max|grad| = 2.14479e-05 (tol = 2e-05, component 1)
2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge with max|grad| = 2.14479e-05 (tol = 2e-05, component 1)
3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge with max|grad| = 2.17682e-05 (tol = 2e-05, component 1)
4: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
  Model failed to converge with max|grad| = 2.17682e-05 (tol = 2e-05, component 1)


#####
##plots

##QQ-plots of p-values
par(mfrow = c(3,2), mar = c(4.1,4.1,1.1,1.1))
maxdiff <- NULL
maxdiff2 <- NULL
maxratio <- NULL
maxratio2 <- NULL
pvdiff <- NULL

runtime <- NULL
runtime4 <- NULL
timeratio <- NULL
timeratio4 <- NULL
runtime.lmerTest <- NULL

#sizeSamples <- c(5e3, 1e4, 1.5e4, 2e4, 2.5e4, 3e4)
for (nSamples in sizeSamples){

load(file = paste0(dirOut, "/simu2_beta_lme4_n", nSamples, ".RData"))

maxdiff <- cbind(maxdiff, 
	c(betas = max(abs(felmm - felme)), 
	't-values' = max(abs(tvlmm - tvlme)), 
	#'p-values' = max(abs(plmm - plme)), 
	sigma2s = max(abs(slmm - slme))))

maxdiff2 <- cbind(maxdiff2, 
	c(betas = max(abs(felmm - felme2)), 
	't-values' = max(abs(tvlmm - tvlme2)), 
	#'p-values' = max(abs(plmm - plme2)), 
	sigma2s = max(abs(slmm - slme2))))

maxratio <- cbind(maxratio, 
	c(betas = max(abs(felmm/felme - 1)), 
	't-values' = max(abs(tvlmm/tvlme - 1)), 
	#'p-values' = max(abs(plmm/plme - 1)), 
	sigma2s = max(abs(slmm/slme - 1))))

maxratio2 <- cbind(maxratio2, 
	c(betas = max(abs(felmm/felme2 - 1)), 
	't-values' = max(abs(tvlmm/tvlme2 - 1)), 
	#'p-values' = max(abs(plmm/plme2 - 1)), 
	sigma2s = max(abs(slmm/slme2 - 1))))

pvdiff <- cbind(pvdiff,
	c(apply(abs(plmm - plme), 1, max), apply(abs(plmm - plme2), 1, max))) 

tunits <- units(rtlmm_mult)

timeratio <- cbind(timeratio, 
	c('lmer-accurate/lmmfit' = sum(rtlme2)/rtlmm_mult[[1]], 
	'lmer-default/lmmfit' = sum(rtlme)/rtlmm_mult[[1]]))
timeratio4 <- cbind(timeratio4, 
	c('lmer-accurate/lmmfit' = sum(rtlme2)/rtlmm_mult[[1]], 
	'lmer-default/lmmfit' = sum(rtlme)/rtlmm_mult[[1]],
	'lmmfit-one/lmmfit' = sum(rtlmm)/rtlmm_mult[[1]]
	))
	
runtime <- cbind(runtime, 
	c('lmer-accurate' = sum(rtlme2), 'lmer-default' = sum(rtlme), lmmfit = rtlmm_mult))
runtime4 <- cbind(runtime4, 
	c('lmer-accurate' = sum(rtlme2), 'lmer-default' = sum(rtlme), 
	'lmmfit-one' = sum(rtlmm), lmmfit = rtlmm_mult))


p1 <- plmm[2, ]
p2 <- plme2[2, ]

q1 <- qqpvalue(p1, plot.it = F)
q2 <- qqpvalue(p2, plot.CI.only = T, add.grid = F)
points(q1$x, q1$y, col = "red", pch = 16, cex = 0.9)
lines(q2$x, q2$y, col = "blue", lwd = 1.5)
legend("topleft", c("lmer", "lmmfit"), title = paste0("n = ", nSamples),
	pch = c(NA, 16), lty = 1, col = c("blue", "red"), bty = "n")
	

rm(felmm, tvlmm, plmm, slmm, rtlmm, rtlmm_mult)
rm(felme, tvlme, plme, slme, rtlme, felme2, tvlme2, plme2, slme2, rtlme2)

##lmerTest:lmer
##both LMM fitting and testing
load(file = paste0(dirOut, "/simu1_beta_n", nSamples, ".RData"))
runtime.lmerTest <- c(runtime.lmerTest, sum(rtlme))
}


##differences between estimates obtained by lmer and lmmfit
colnames(maxdiff) <- paste0("n=", sizeSamples)
colset <- 1 + 1:nrow(maxdiff)
barplot(maxdiff, ylim = c(0, 1.2*max(maxdiff)), 
	ylab = "Maximum absolute differences of the estimates", beside = T, col = colset)
#rownames(maxdiff)
#txt <- c(expression(beta), "t-values", expression(paste(sigma^2, "s")))
txt <- c(expression(beta), "t-values", expression(sigma[k]^2))
legend("topleft", txt, pch = 15, col = colset, bty = "n")
	
	
##lmer (accurate) with control changed
colnames(maxdiff2) <- paste0("n=", sizeSamples)
colset <- 1 + 1:nrow(maxdiff2)
barplot(maxdiff2, ylim = c(0, 1.2*max(maxdiff2)), 
	ylab = "Maximum absolute differences of the estimates", beside = T, col = colset)
#rownames(maxdiff2)
#txt <- c(expression(beta), "t-values", expression(paste(sigma^2, "s")))
txt <- c(expression(beta), "t-values", expression(sigma[k]^2))
legend("topleft", txt, pch = 15, col = colset, bty = "n")
	

##ratios between estimates obtained by lmer and lmmfit
colnames(maxratio) <- paste0("n=", sizeSamples)
colset <- 1 + 1:nrow(maxratio)
barplot(maxratio, ylim = c(0, 1.2*max(maxratio)), 
	ylab = "Max |ratios  - 1|", beside = T, col = colset)
#rownames(maxdiff)
#txt <- c(expression(beta), "t-values", expression(paste(sigma^2, "s")))
txt <- c(expression(beta), "t-values", expression(sigma[k]^2))
legend("topleft", txt, pch = 15, col = colset, bty = "n")
	
	
##
colnames(maxratio2) <- paste0("n=", sizeSamples)
colset <- 1 + 1:nrow(maxratio2)
barplot(maxratio2, ylim = c(0, 1.2*max(maxratio2)), 
	ylab = "Max |ratios  - 1|", beside = T, col = colset)
#rownames(maxdiff)
#txt <- c(expression(beta), "t-values", expression(paste(sigma^2, "s")))
txt <- c(expression(beta), "t-values", expression(sigma[k]^2))
legend("topleft", txt, pch = 15, col = colset, bty = "n")
	
	

##run times
colnames(runtime) <- paste0("n=", sizeSamples)
#par(mfrow = c(1, 2), mar = c(5.1, 4.1, 4.1, 1.1))

plot(runtime[1,], ylim = range(runtime), axes = F, type = "n", 
	xlab = NA, ylab = paste0("Run time in ", tunits))
#grid()
box()
axis(1, 1:ncol(runtime), labels = colnames(runtime))
axis(2)
colset <- c("red", "green", "blue") #1 + 1:nrow(runtime)
for (i in 1:nrow(runtime)){
	lines(runtime[i,], col = colset[i], lty = 3)
	points(runtime[i, ], col = colset[i], pch = 16, cex = 1.2)
}
legend("topleft", rownames(runtime), pch = 16, lty = 3, col = colset, bty = "n")
	

##ratio
colnames(timeratio) <- paste0("n=", sizeSamples)

plot(timeratio[1,], ylim = range(timeratio), axes = F, type = "n", 
	xlab = NA, ylab = "Ratio of run time")
#grid()
box()
axis(1, 1:ncol(timeratio), labels = colnames(timeratio))
axis(2)
colset <- c("red", "green") #1 + 1:nrow(timeratio)
for (i in 1:nrow(timeratio)){
	lines(timeratio[i,], col = colset[i], lty = 3)
	points(timeratio[i, ], col = colset[i], pch = 16, cex = 1.2)
}
legend("topleft", rownames(timeratio), pch = 16, lty = 3, col = colset, bty = "n")

	
##run times 4
rtime <- runtime4
colnames(rtime) <- paste0("n=", sizeSamples)

plot(rtime[1,], ylim = range(rtime), axes = F, type = "n", 
	xlab = NA, ylab = paste0("Time in ", tunits))
#grid()
box()
axis(1, 1:ncol(rtime), labels = colnames(rtime))
axis(2)
colset <- c("red", "green", "purple", "blue") #1 + 1:nrow(rtime)
for (i in 1:nrow(rtime)){
	lines(rtime[i,], col = colset[i], lty = 3)
	points(rtime[i, ], col = colset[i], pch = 16, cex = 1.2)
}
legend("topleft", rownames(rtime), pch = 16, lty = 3, col = colset, bty = "n")
	
	
##ratio
colnames(timeratio4) <- paste0("n=", sizeSamples)

plot(timeratio4[1,], ylim = range(timeratio4), axes = F, type = "n", 
	xlab = NA, ylab = "Ratio of run time")
#grid()
box()
axis(1, 1:ncol(timeratio4), labels = colnames(timeratio4))
axis(2)
colset <- c("red", "green", "purple") #1 + 1:nrow(timeratio4)
for (i in 1:nrow(timeratio4)){
	lines(timeratio4[i,], col = colset[i], lty = 3)
	points(timeratio4[i, ], col = colset[i], pch = 16, cex = 1.2)
}
legend("topleft", rownames(timeratio4), pch = 16, lty = 3, col = colset, bty = "n")

	
##
##run times for lmerTest:lmer
rtime <- rbind('lmerTest::lmer' = runtime.lmerTest, runtime[-1, ])
colnames(rtime) <- paste0("n=", sizeSamples)
rownames(rtime)[grep("default", rownames(rtime))] <- 'lme4::lmer'
rtime

plot(rtime[1,], ylim = range(rtime), axes = F, type = "n", 
	xlab = NA, ylab = paste0("Time in ", tunits))
#grid()
box()
axis(1, 1:ncol(rtime), labels = colnames(rtime))
axis(2)
colset <- c("red", "green", "blue") #1 + 1:nrow(rtime)
for (i in 1:nrow(rtime)){
	lines(rtime[i,], col = colset[i], lty = 3)
	points(rtime[i, ], col = colset[i], pch = 16, cex = 1.2)
}
legend("topleft", rownames(rtime), pch = 16, lty = 3, col = colset, bty = "n")
	

##################################################



	
##################################################
##Simulation analysis II
##null hypothesis: sigma2b = 0 
##
##simulating data with 900 sigma2b = 0 and 100 sigma2b ~ uniform(0.1, 0.5) 

##########	
##Simulating data
##trt: treatments
##cls: sample clusters
##Y: responses
##model: y = trt*beta + cls*b + e

##null hypothesis: sigma2b = 0 
##random effects: b ~ Norm(0, sigma2b)

##variance components
sigma2e <- 1  ##for residual errors
beta <- 0.5

##number of treatments
ntrt <- 2
##number of clusters
ncls <- 16
##sample sizes
nSamples <- 10000
##number of simulations
nSimu <- 1000 
#nRandomEffects <- 50
nRandomEffects <- 100
	
		
set.seed(231212)
##sample clusters and treatments
dat <- data.frame(
	cls = as.factor(sample.int(ncls, nSamples, replace = T)),
	trt = LETTERS[sample.int(ntrt, nSamples, replace = T)])

#####	
##responses
##no random effects
Y <- matrix(rnorm(nrow(dat)*nSimu, mean = 0, sd = sqrt(sigma2e)), nrow = nrow(dat), ncol = nSimu)
colnames(Y) <- paste0("Y", 1:nSimu)

##treatment effect
index <- (dat$trt == "B")
for (j in 1:ncol(Y)){
	betaj <- rnorm(1, m = beta, sd = 1)
	Y[index, j] <- Y[index, j] + rnorm(sum(index), mean = betaj, sd = 0.1)
	}
	
##random effects
for (i in 1:nRandomEffects) {
for (cls in unique(dat$cls)){
	index <- (dat$cls == cls)
	#Y[index, i] <- Y[index, i] + rnorm(1, mean = 0, sd = sqrt(sigma2b)) 
	Y[index, i] <- Y[index, i] + rnorm(1, mean = 0, sd = runif(1, 0.1, 0.5)) 
	}##cls
	}##i


##########
##LMM fit and test
##lmmfit
##design matrix for fixed effects
X <- model.matrix(~ trt, data = dat)
##design matrix for random effects
Z <- model.matrix(~ 0 + cls, data = dat)
d <- ncol(Z)

##LMM fit
t1 <- Sys.time()
fit <- lmmfit(Y = Y, X, Z, d=d, epsilon = 1e-8)
test <- lmmtest(fit)
	t2 <- Sys.time()
	difftime(t2, t1)


##
range(fit$niter)
range(fit$dlogL)
range(fit$theta[1,])
range(fit$theta[2,])

lambdaZ <- svd(Z)$d
-1/max(lambdaZ)
gamma <- fit$theta[1, ]/fit$theta[2,]
all(gamma > - 1/max(lambdaZ))
range(gamma)


##hypothesis tests for variance components
fit$theta[, 1:5]
stats0 <- fit$theta[1, ]/fit$se[1, ]
p <- pnorm(stats0, lower.tail = F)

##Normal or CHi-square approximation provide 
##an excessively conservative approximation to the null distribution. 
##This is safe when the evidence against the null is overwhelming.
qqpvalue(p[-(1:nRandomEffects)])
qqpvalue(p)
hist(p)


##########
##permutation
nPerm <- 10000

##design matrix for fixed effects
X <- model.matrix(~ trt, data = dat)
##design matrix for random effects
Z <- model.matrix(~ 0 + cls, data = dat)
d <- ncol(Z)
trtset <- unique(dat$trt)

set.seed(231218)

theta <- NULL
stats <- NULL
t1 <- Sys.time()
for (k in 1:nPerm){
	index <- rep(NA, nrow(dat))
	for (trt in trtset){
		itrt <- which(dat$trt == trt)
		index[itrt] <- sample(itrt, replace = F)
	}
	##LMM fit
	fitk <- lmmfit(Y = Y[index, ], X, Z, d=d, epsilon = 1e-8)
	theta <- rbind(theta, fitk$theta[1, ])
	stats <- rbind(stats, fitk$theta[1, ]/fitk$se[1, ])
}

t2 <- Sys.time()
difftime(t2, t1)
##Time difference of 8.418366 hours

#thetaPerm <- theta
save(fit, theta, stats, nPerm, nSamples, nRandomEffects,
	file = paste0(dirOut, "/simu_theta_permutation.RData"))


##########
##p-values and QQ-plot

##lmmfit$theta
##permutation beta
load(file = paste0(dirOut, "/simu_theta_permutation.RData"))

##p-values for testing the zero variance components
p <- pnorm(fit$theta[1, ]/fit$se[1, ], lower.tail = F)

##Permutation p-values
pp <- (1+ colSums(theta >= rep(fit$theta[1, ], each = nrow(theta))))/(1+nrow(theta))

#nRandomEffects, 100

qqp <- qqpvalue(pp[-(1:nRandomEffects)], plot.it = F)
qq <- qqpvalue(p[-(1:nRandomEffects)], plot.CI.only = T, add.grid = F)
colset <- c("green", "blue")
#lines(qq$x, qq$y, col = colset[1], lwd = 1.5)
points(qq$x, qq$y, col = colset[1], pch = 16, cex = 0.6)
points(qqp$x, qqp$y, col = colset[2], pch = 16, cex = 0.6)
legend("topleft", c("lmmfit", "permutation"), pch = 16, col = colset, bty = "n")

	figFile <- paste0(dirFig, "/simu_theta_qqplot")
	#save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)


##################################################

	
	
##################################################
##Simulation analysis III
##Generate counts data based on negative binomial distribution.
##Simulating data: 
##4900 genes with beta = 0 and 
##100 genes with beta ~ runiform(0.25, 1)
 

##########	
##Simulating data
##trt: treatments
##cls: sample clusters
##Y: responses
##model: 
##mu = trt*beta + cls*b + e
##Y ~ rnbinom(n, size = 1, mu = mu)

##variance components
sigma2b <- 0.16 ##for random effects

##number of treatments
ntrt <- 2
##number of clusters (subjects)
ncls <- 16
##number of simulations
nSimu <- 5000 
##number of fixed effects != 0, 
##the first 100 out of 5000 simulations
nFixedEffects <- 100 
##sample sizes
nSamples <-  2e4

##
set.seed(240108)

##sample clusters and treatments
dat <- data.frame(
	cls = as.factor(sample.int(ncls, nSamples, replace = T)),
	trt = LETTERS[sample.int(ntrt, nSamples, replace = T)])
	
##counts
##generated by rnbinom
dispersion <- 0.1

Y <- matrix(NA, nrow = nrow(dat), ncol = nSimu)
colnames(Y) <- paste0("Y", 1:nSimu)

##
indexB <- (dat$trt == "B")

for (j in 1:nFixedEffects) {
betaj <- runif(1, 0.25, 1)
for (cls in unique(dat$cls)){
	index <- (dat$cls == cls)
	m <- rnorm(1, sd = sqrt(sigma2b)) ##mean of log-Y
	Y[index, j] <- rnbinom(n = sum(index), mu = exp(m), size = dispersion)
	i <- (indexB & index)
	if (length(i) > 0)
	Y[i, j] <- rnbinom(n = sum(i), mu = exp(m + betaj), size = dispersion)
	}##cls
	}##i

for (j in (1+nFixedEffects):ncol(Y)) {
for (cls in unique(dat$cls)){
	index <- (dat$cls == cls)
	m <- rnorm(1, sd = sqrt(sigma2b)) ##mean of log-Y
	Y[index, j] <- rnbinom(n = sum(index), mu = exp(m), size = dispersion)
	}##cls
	}##i

##
range(Y)
#[1]   0 331
sum(Y > 0)
#[1] 21361377
sum(Y == 0)/length(Y)
#[1] 0.7863862
sum(Y <= 1)/length(Y)
#[1] 0.857389
sum(Y <= 2)/length(Y)
#[1] 0.8926895


##########
##lmmfit
##design matrix for fixed effects
nGenes <- rowSums(Y)
#hist(log(nGenes))

X <- model.matrix(~ log(nGenes) + trt, data = dat)

##design matrix for random effects
Z <- model.matrix(~ 0 + cls, data = dat)
d <- ncol(Z)

##LMM fit
t1 <- Sys.time()
fit <- lmmfit(Y = log2(1+Y), X, Z, d=d, epsilon = 1e-8)
	t2 <- Sys.time()
	difftime(t2, t1)

timeUnits <- "secs"
rtlmm <- difftime(t2, t1, units = timeUnits)

##LMM test	
test <- lmmtest(fit)
head(test)

p <- test[, "trtB_pvalue"]
hist(p)
sum(p <= 0.05)
range(p[(1:nFixedEffects)])
range(p[-(1:nFixedEffects)])

qqpvalue(p)
qqpvalue(p[-(1:nFixedEffects)])

fdr <- p.adjust(p, method = "BH")
sum(fdr <= 0.05)
sum(which(fdr <= 0.05) <= nFixedEffects)


##########
##nebula

##random effect (sample groups)
Z <- dat$trt
table(Z)
length(Z) 

##
##NBGMM' is for fitting a negative binomial gamma mixed model. 
##'PMM' is for fitting a Poisson gamma mixed model. 
##'NBLMM' is for fitting a negative binomial lognormal mixed model (the same model as that in the lme4 package).

#model <- "NBGMM"
model <- "NBLMM"

t1 <- Sys.time()
negbn <- NULL
#The cells in the count matrix need to be grouped by the subjects
o <- order(Z)
negbn <- nebula(count = t(Y[o, ]), id = as.factor(Z[o]), pred = X[o, ], 
	model = model, covariance = T, output_re = T)
	t2 <- Sys.time()

difftime(t2, t1)

timeUnits <- "secs"
rtnebula <- difftime(t2, t1, units = timeUnits)


##########
##comparison

##LMM test
test <- lmmtest(fit)
head(test)
plmm <- test[, "trtB_pvalue"]
hist(plmm)
tvlmm <- test[, "trtB_t"]

##nebula test
st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)
pnbn <- st[, "p_trtB"]
tvnbn <- st[, "logFC_trtB"]/st[, "se_trtB"]
range(pnbn - 2*pnorm(-abs(tvnbn)))

##t-values
plot(tvlmm, tvnbn, xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.65, col = "green")
abline(0, 1, col = "green")
i <- 1:nFixedEffects
points(tvlmm[i], tvnbn[i], col = "blue", cex = 0.65)
	
##histograms
hist(plmm, xlab = "lmmfit p-values", main = NA)
hist(pnbn, xlab = "nebula p-values", main = NA)

##QQplots
qq1 <- qqpvalue(plmm, plot.it = F)
qq2 <- qqpvalue(pnbn, plot.CI.only = T, add.grid = F)
colset <- c("blue", "green")
lines(qq1$x, qq1$y, col = colset[1])
lines(qq2$x, qq2$y, col = colset[2])
points(qq1$x, qq1$y, col = colset[1], cex = 0.6, pch = 16)
points(qq2$x, qq2$y, col = colset[2], cex = 0.6, pch = 16)
legend("topleft", c("lmmfit", "nebula"), pch = 16, lty = 1, col = colset, bty = "n")

##no effects
j <- (1+nFixedEffects):length(plmm)
qq1 <- qqpvalue(plmm[j], plot.it = F)
qq2 <- qqpvalue(pnbn[j], plot.CI.only = T, add.grid = F)
colset <- c("blue", "green")
lines(qq1$x, qq1$y, col = colset[1])
lines(qq2$x, qq2$y, col = colset[2])
points(qq1$x, qq1$y, col = colset[1], cex = 0.6, pch = 16)
points(qq2$x, qq2$y, col = colset[2], cex = 0.6, pch = 16)
legend("topleft", c("lmmfit", "nebula"), pch = 16, lty = 1, col = colset, bty = "n")
	


##AUCs

cutoff <- 10^(seq(-20, 0, by = 0.01))
cutoff <- cutoff[cutoff < max(plmm, pnbn)]
i <- 1:nFixedEffects
cutoff <- cutoff[cutoff >= max(min(plmm[-i]), min(pnbn[-i]))]

##lmmfit
pv <- plmm
p0 <- pv[-(1:nFixedEffects)]
p1 <- pv[1:nFixedEffects]
auclmm <- AUC(p0 = p0, p1 = p1, cutoff = cutoff, plot.it = F)

##nebula
pv <- pnbn
p0 <- pv[-(1:nFixedEffects)]
p1 <- pv[1:nFixedEffects]
aucnbn <- AUC(p0 = p0, p1 = p1, cutoff = cutoff, plot.it = F)

##plot
plot(auclmm$FPR, auclmm$TPR, log = "x", type = "n", 
	xlab = "False positive rate", ylab = "True positive rate")
lines(auclmm$FPR, auclmm$TPR, col = colset[1], lty = 1)
lines(aucnbn$FPR, aucnbn$TPR, col = colset[2], lty = 2)
txt <- paste0(c("lmmfit AUC = ", "nebula AUC = "), round(c(auclmm$AUC, aucnbn$AUC), 4))
legend("bottomright", txt, col = colset, lty = 1:2, bty = "n")


##################################################


	