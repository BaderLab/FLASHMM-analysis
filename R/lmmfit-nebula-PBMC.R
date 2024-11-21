
##################################################
##load data
##PBMC data
datafile <- paste0(dirData, "/PBCMdata_ExperimentHubKang_muscat.RData")
load(file = datafile)

rm(logcounts)

##'sample_id': unique sample identifiers
##'cluster_id': subpopulation (cell cluster) assignments 
##'group_id': experimental group/condition (control/treatment, healthy/diseased)
table(coldata$group_id)
# ctrl  stim 
#13463 13357 

table(coldata$sample_id)
coldata$sample_id <- gsub("stim", "", gsub("ctrl", "", coldata$sample_id))
table(coldata$sample_id)
# 101 1015 1016 1039  107 1244 1256 1488 
#2109 5374 3963 1075 1164 3733 4394 5008 

##counts
dim(counts)
#[1]  7118 26820

##Filtering genes
##Mixed models in muscat package - 3.4 Cell-level analysis:
##(1) subpopulations with at least 10 cells in at least 2 samples (not necessary?)
##(2) genes with a count >= 1 in at least 20 cells

all(colnames(counts) == rownames(coldata))

##number of celss
nCells <- rowSums(counts > 0)
#hist(log2(nCells))

minCells <- 2^5

##number of cells in a group_id
nCellsgrp <- do.call(cbind, 
		tapply(1:ncol(counts), as.factor(coldata$group_id), 
		function(j) rowSums(counts[, j, drop = F] > 0))
		)

minCellsgrp <- 10

##number of counts
ncts <- rowSums(counts)
#hist(log2(ncts)) #outliers in the upper tail

maxCounts <- 2^20
sum(ncts > maxCounts)

##number of counts in a group_id
nctsgrp <- do.call(cbind, 
		tapply(1:ncol(counts), as.factor(coldata$group_id), 
		function(j) rowSums(counts[, j, drop = F]))
		)


##nebula filtering:
##Filtering out low-expressed genes can be specified by cpc=0.005 (i.e., counts per cell<0.5%). 
##The argument cpc is defined by the ratio between the total count of the gene and the number of cells.

cpc <- rowSums(counts)/ncol(counts)
sum(cpc <= 0.005) 
#[1] 31


##Filtering
minCells <- 32 #2^5
minCellsgrp <- 15
maxCounts <- 2^20
minCountsgrp <- 2*minCellsgrp
minCountsgrp
mincpc <- 0.005

index <- (nCells >= minCells) & (rowSums(nCellsgrp >= minCellsgrp) >= 2)
index <- index & (ncts <= maxCounts) & (rowSums(nctsgrp >= minCountsgrp) >= 2)
sum(index)
index <- index & (cpc > mincpc)
sum(index)

counts <- counts[index, ]
rm(index)
dim(counts)
#[1]  7017 26820



##################################################
#NEBULA
##https://github.com/lhe17/nebula
##Checking convergence for the summary statistics and quality control
##  1: The convergence is reached due to a sufficiently small improvement of the function value.
##-10: The convergence is reached because the gradients are close to zero 
##     (i.e., the critical point) and no improvement of the function value can be found.
##
##Depending on the concrete application, 
##the estimated gene-specific overdispersions can also be taken into consideration in quality control. 
##For example, when testing differential expression for a variable, 
##genes with a very large estimated cell-level overdispersion should be filtered out because such genes have huge unexplained noises.
##
##If the variable of interest is subject-level, 
##genes with a very large subject-level overdispersion (>1) 
##should be removed or interpreted cautiously as well.
##
##The NBLMM is the same model as that adopted in the glmer.nb function in the lme4 R package, 
##but is computationally much more efficient by setting method='LN'. 


##raw counts
Y <- counts
dim(Y) 
#[1]  7017 26820

nGenes <- colSums(Y)

##
rm(counts)

##nebula
##fixed effect desigm matrix
X <- model.matrix(~ log(nGenes) + cluster_id + cluster_id:group_id, data = coldata)
colnames(X) <- gsub("cluster_id", "", colnames(X))
colnames(X) <- gsub("group_id", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))

head(X)
dim(X)

##random effect (sample groups)
Z <- coldata$sample_id
table(Z)
length(Z) #[1] 26820


##model:
##NBGMM' is for fitting a negative binomial gamma mixed model. 
##'PMM' is for fitting a Poisson gamma mixed model. 
##'NBLMM' is for fitting a negative binomial lognormal mixed model (the same model as that in the lme4 package).

model <- "NBLMM"

t1 <- Sys.time()
negbn <- NULL
##The cells in the count matrix need to be grouped by the subjects
o <- order(Z)
negbn <- nebula(count = Y[, o], id = as.factor(Z[o]), pred = X[o, ], 
	model = model, covariance = T, output_re = T)
t2 <- Sys.time()

difftime(t2, t1)
#Time difference of 53.49769 mins

timeUnits <- "secs"
rtnebula <- difftime(t2, t1, units = timeUnits)

##nebula outputs
##summary (statistics): 
##The estimated coefficient, standard error and p-value for each predictor.
str(negbn)

plot(negbn$overdispersion$Subject)
plot(negbn$overdispersion$Cell)

##
##https://github.com/lhe17/nebula
##  1: The convergence is reached due to a sufficiently small improvement of the function value.
##-10: The convergence is reached because the gradients are close to zero 
##     (i.e., the critical point) and no improvement of the function value can be found.
table(negbn$convergence)
# -40  -10    1 
#   5  101 6911 


##################################################
##LMM fit

dim(Y)
#[1]  7017 26820

##log-transformation
##log2(1 + counts)
##log2(1+Y) by groupping to reduce data size
ngrp <- 3
sizegrp <- round(nrow(Y)/ngrp)
for (i in 1:ngrp){
	j <- (1+(i-1)*sizegrp):(min(nrow(Y), i*sizegrp))
	print(range(j))
	Y[j, ] <- log2(1 + Y[j, ])
}

dim(Y) 
#[1]  7017 26820

##fixed effect desigm matrix
X <- model.matrix(~ log(nGenes) + cluster_id + cluster_id:group_id, data = coldata)
colnames(X) <- gsub("cluster_id", "", colnames(X))
colnames(X) <- gsub("group_id", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub(" ", "_", colnames(X))

head(X)
dim(X)

##random effects
##sample groups
Z <- model.matrix(~ 0 + sample_id, data = coldata)
colnames(Z) <- gsub(".*_id", "", colnames(Z))
dim(Z) #[1] 26820    8
head(Z)

d <- ncol(Z)
n <- nrow(X)

##########
##lmmfitSS
class(Y)
class(X)

t1 <- Sys.time()
XY <- t(Y%*%X)
ZX <- t(Z)%*%X #as.matrix(t(Z)%*%X)
ZY <- t(Y%*%Z) #as.matrix(t(Z)%*%Y)
ZZ <- t(Z)%*%Z #as.matrix(t(Z)%*%Z)

XXinv <- ginv(t(X)%*%X)
Ynorm <- rowSums(Y*Y) #colSums(Y*Y)
t2 <- Sys.time()
	difftime(t2, t1)
	
	class(XY)
	class(ZY)

##dgeMatrix class
##takes more run time than matrix class!
t1 <- Sys.time()
fit <- lmmfitSS(XY, ZX, ZY, ZZ = ZZ, XXinv = XXinv, Ynorm = Ynorm, n = n, d = d, 
	max.iter = 100, epsilon = 1e-5)
		
t2 <- Sys.time()
difftime(t2, t1) 
#Time difference of 2.434709 mins

##Matrix class
XY <- as.matrix(XY)
ZY <- as.matrix(ZY)
t1 <- Sys.time()
fitss <- lmmfitSS(XY, ZX, ZY, ZZ = ZZ, XXinv = XXinv, 
	Ynorm = Ynorm, n = n, d = d, max.iter = 100, epsilon = 1e-5)
		
t2 <- Sys.time()
difftime(t2, t1)
#Time difference of 32.95579 secs

##comparison
##approximately identical
identical(fit, fitss)
#[1] FALSE
range(fit$dlogL - fitss$dlogL) 
#[1] -8.381903e-09  1.810491e-05
range(fit$niter - fitss$niter) 
#[1] 0 0
range(fit$theta - fitss$theta) 
#[1] -4.329870e-15  7.702172e-15
range(fit$se - fitss$se)       
#[1] -1.457168e-15  1.092876e-16
range(fit$coef - fitss$coef)   
#[1] -2.580158e-13  2.744471e-13
range(fit$cov - fitss$cov)     
#[1] -1.239141e-16  1.294104e-15


##########
##Fit LMM by lmmfit.nt
##Operating on "matrix" or "array"  is faster than "dgCMatrix"!!!
Y <- as.matrix(Y)

timeUnits <- "secs"
maxIter <- 100 
epsilon <- 1e-5 

t1 <- Sys.time()
fitnt <- lmmfit.nt(Y = Y, X = X, Z = Z, d = d, max.iter = maxIter, epsilon = epsilon)
t2 <- Sys.time()
	difftime(t2, t1)
	#Time difference of 39.17241 secs
	rtlmm <- difftime(t2, t1, units = timeUnits) 

table(fitnt$niter)

##comparison between lmmfitSS and lmmfit.nt
fit <- fitnt
identical(fit, fitss)
#[1] FALSE
range(fit$dlogL - fitss$dlogL) 
#[1] -0.002343416  0.002813816
range(fit$niter - fitss$niter) 
#[1] 0 0
range(fit$theta - fitss$theta) 
#[1] -8.995027e-13  1.172396e-13
range(fit$se - fitss$se)       
#[1] -7.765490e-15  1.484923e-15
range(fit$coef - fitss$coef)   
#[1] -5.328182e-12  1.014300e-12
range(fit$cov - fitss$cov)     
#[1] -1.692743e-14  1.001282e-14


##########
##Fit LMM by lmmfit
##transpose
Y <- t(Y)
dim(Y)

#timeUnits <- "secs"
#maxIter <- 100 
#epsilon <- 1e-5 

t1 <- Sys.time()
fit <- lmmfit(Y = Y, X = X, Z = Z, d = d, max.iter = maxIter, epsilon = epsilon)
t2 <- Sys.time()
	difftime(t2, t1) 
	#Time difference of 37.8269 secs
	
##comparison between lmmfit.nt and lmmfit
identical(fit, fitnt)
#[1] TRUE


##################################################
##comparison of results

##########
##lmmfit

rtlmm; timeUnits; maxIter; epsilon
#Time difference of 49.96773 secs
#[1] "secs"
#[1] 200
#[1] 1e-08

##fixed effects
felmm <- fit$coef

##variance components
slmm <- fit$theta

##LMM tests
test <- lmmtest(fit)
dim(test)
head(test)

##t-values
tvlmm <- test[, grep("_t", colnames(test)), drop = F]
dim(tvlmm)
head(tvlmm)

##p-values
plmm <- test[, grep("_pvalue", colnames(test)), drop = F]
dim(plmm)
head(plmm)


##########
##nebula

rtnebula
#Time difference of 3209.861 secs
rtnebula[[1]]/rtlmm[[1]]
#[1] 64.23869
str(negbn)

table(negbn$convergence)
# -40  -10    1 
#   5  101 6911 

st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)

any(is.na(st))
sum(is.na(st))
st[rowSums(is.na(st)) > 0, ]
st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)

##fixed effects, se, and p-values
iFC <- grep("logFC_", colnames(st))
ise <- grep("se_", colnames(st))
ipv <- setdiff(grep("p_", colnames(st)), c(iFC, ise))

##fixed effects
b <- st[, iFC]

##se
se <- st[, ise]

##t-values
tv <- b/se

##p-values
pv <- st[, ipv]

range(pv - 2*pnorm(-abs(tv)), na.rm = T)
#[1] -2.331468e-15  2.386980e-15


##########
##comparison of lmmfit and nebula
##for the genes that lmmfit and nebula are convergent.

all(rownames(tv) == rownames(tvlmm))
any(is.na(tv))
any(is.na(tvlmm))

##genes with convergence
i <- (negbn$overdispersion$Subject < 1) 
i <- i & (negbn$overdispersion$Cell < 100) & (negbn$convergence >= -10)
sum(i)
#[1] 7000
i <- i & (colSums(abs(fit$dlogL) < 1e-5) == nrow(fit$dlogL))
sum(i)
#[1] 6967
i <- i & (apply(!is.na(tv), 1, all)) & (apply(!is.na(tvlmm), 1, all))
sum(i)
#[1] 6958


##t-values
j <- 2:ncol(tv)
plot(as.matrix(tvlmm[i, j]), as.matrix(tv[i, j]), 
	xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1)
	
##p-values
j <- 2:ncol(pv)
#j <- grep(":stim", colnames(pv))
plot(as.matrix(-log10(plmm[i, j])), as.matrix(-log10(pv[i, j])), 
	xlab = "lmmfit -log10(p-values)", ylab = "nebula -log10(p-values)", cex = 0.6)
abline(0, 1)
	
	
##cell-type specific
jset <- grep(":stim", colnames(pv))
nc <- round(sqrt(length(jset)))
nr <- ceiling(length(jset)/nc)
par(mfrow = c(nr, nc), mar = c(5.1,4.1,2.1,1.1))
for (j in jset){
	nm <- gsub("p", "+", gsub("_t", "", colnames(tvlmm)[j]))
	plot(tvlmm[i,j], tv[i,j], cex = 0.6,
	xlab = "lmmfit t-values", ylab = "nebula t-values", main = nm, cex.main = 0.8)
	abline(0,1)
	}
	

##histograms of p-values
jset <- grep(":stim", colnames(pv))
par(mfrow = c(4, 4))
for (j in jset){
	h1 <- hist(pv[i,j], plot = F)
	h2 <- hist(plmm[i,j], plot = F)
	ylim <- c(0, 1.1*max(h1$counts, h2$counts))
	nm <- gsub("p", "+", gsub("_pvalue", "", colnames(plmm)[j]))

	par(mar = c(5.1,4.1,2.1,0))
	#nm <- gsub("p_", "", colnames(pv)[j])
	hist(pv[i,j], ylim = ylim, xlab = "nebula p-values", 
		main = nm, cex.main = 0.8, cex.axis = 0.8)

	par(mar = c(5.1,2.1,2.1,2.1))	
	#nm <- gsub("_pvalue", "", colnames(plmm)[j])
	hist(plmm[i, j], ylim = ylim, xlab = "lmmfit p-value", ylab = NULL, 
		main = nm, cex.main = 0.8, cex.axis = 0.8, yaxt = "n")
	}



##barplot of common genes in the top
all(rownames(plmm) == rownames(pv))
glist <- rownames(plmm)[i]

jset <- grep(":stim", colnames(pv))
ntops <- c(250, 500)

commonGenes <- NULL
for (nTop in ntops){
cgenes <- NULL	
#par(mfrow = c(4, 4))
for (j in jset){
	p <- NULL
	p <- pv[i,j]
	g1 <- glist[order(p)[1:nTop]]
	
	p <- NULL
	p <- plmm[i,j]
	g2 <- glist[order(p)[1:nTop]]

	##proportion of the common genes in the top
	cgenes <- c(cgenes, length(intersect(g1, g2))/nTop)
	}
	commonGenes <- cbind(commonGenes, cgenes)
}
rownames(commonGenes) <- gsub("p_", "+_", gsub(":stim_pvalue", "", colnames(plmm)[jset]))
colnames(commonGenes) <- ntops
commonGenes

colset <- rainbow(nrow(commonGenes))

par(mar = c(5.1, 4.1, 4.1, 8.1))
barplot(commonGenes, beside = T, ylim = c(0, 1), col = colset,
	xlab = "Number of top genes", ylab = "Proportion of the common genes in the top")
legend("topleft", rownames(commonGenes), inset = c(1, 0), xpd = T, pch = 15, 
	col = colset, bty = "n", cex = 0.7)

##################################################


	