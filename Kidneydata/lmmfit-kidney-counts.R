##########
library(MASS)
library(Matrix)
#library(lme4)
#library(lmerTest)
library(nebula)
#library(dplyr)
library(Seurat)

#library(patchwork)
#library(lme4)
#library(nnls)
#library(nlme)

#library(monocle)
#packageVersion('monocle')

##the defined R functions
source(paste0(dirwork, "/ut/Bader/MM/R/lmm.R"))
source(paste0(dirwork, "/ut/Bader/MM/R/lmmfit.R"))
#source(paste0(dirwork, "/ut/Bader/MM/R/lmmtest.R"))

##################################################


##################################################
#options(Seurat.object.assay.version = "v5")

##load data
#datafile <- "/scData/Human_Kidney_data.rds"
datafile <- "/scData/Kidney_raw_counts_decontaminated.rds"
dat <- readRDS(file = datafile)
#str(dat)
dat

#####
#counts <- GetAssayData(object = dat, slot = "counts")
counts <- GetAssayData(dat, layer = 'counts', assay = 'RNA')

class(counts)
dim(counts)
#[1] 22484 27677
min(counts)
min(counts[counts > 0])
#[1] 5.097617e-05
max(counts)
#[1] 25332.21


i <- 1
table(counts[i,])
table(round(counts[i,]))


#####
coldata <- dat@meta.data
dim(coldata)
head(coldata)
all(rownames(coldata) == colnames(counts))
#[1] TRUE


##########
##Filter cells

##the cells included
dim(counts)
##number of features
nFeature <- colSums(counts > 0)
minnFeature <- 100 
##library size
libsize <- colSums(counts)
minlibsize <- 2^9
maxlibsize <- 2^16
##number of cells in a cell-type (cluster)
nCellsType <- table(coldata$Cell_Types_Broad)
minnCellsType <- 20

##filtering
all(colnames(counts) == rownames(coldata))

j <- (nFeature >= minnFeature) & (libsize >= minlibsize) & (libsize <= maxlibsize)
sum(j)
#[1] 27566
##remove "Podocyte"
j <- j & (coldata$Cell_Types_Broad %in% names(nCellsType)[nCellsType >= minnCellsType])
sum(j)
#[1] 27550

counts <- counts[, j]
coldata <- coldata[j, ]
rm(j)

all(colnames(counts) == rownames(coldata))
sort(table(coldata$Cell_Types_Broad))

##########
##Filtering genes
##Mixed models in muscat package - 3.4 Cell-level analysis:
##(1) subpopulations with at least 10 cells in at least 2 samples (not necessary?)
##(2) genes with a count >= 1 in at least 20 cells

all(colnames(counts) == rownames(coldata))

##(1) number of celss
#nCells <- rowSums(counts > 0)
nCells <- rowSums(counts >= 1)
#hist(log2(nCells))

minCells <- 2^4 #2^5 #2^6
sum(nCells >= minCells)

##(2) number of cells in a group_id (sex)
nCellsgrp <- do.call(cbind, 
		tapply(1:ncol(counts), as.factor(coldata$sex), 
		function(j) rowSums(counts[, j, drop = F] >= 1))
		#function(j) rowSums(counts[, j, drop = F] > 0))
		)

dim(nCellsgrp)
head(nCellsgrp)

minCellsgrp <- 10 #50#used to be 10
sum(rowSums(nCellsgrp >= minCellsgrp) == ncol(nCellsgrp))


##(3) number of counts
nCounts <- rowSums(counts)
#hist(log2(nCounts)) #outliers in the upper tail
#boxplot(log2(nCounts))

maxCounts <- 2^20
sum(nCounts > maxCounts)
minCounts <- 2^6
sum(nCounts >= minCounts)



##(4) nebula filtering:
##Filtering out low-expressed genes can be specified by cpc=0.005 (i.e., counts per cell<0.5%). 
##The argument cpc is defined by the ratio between the total count of the gene and the number of cells.

cpc <- rowSums(counts)/ncol(counts)
range(cpc[cpc > 0])
sum(cpc <= 0.005) 
sum(cpc > 0.005) 
sum(cpc < 0.001) 
sum(cpc >= 0.001) 

#####
##Filtering
minCells <- 16 #2^6#32 #2^5#
minCellsgrp <- 10 #25 #15#50
minCounts <- 2^6
maxCounts <- 2^20
#minCountsgrp <- 2*minCellsgrp
#minCountsgrp
mincpc <- 0.005

length(nCells)
dim(nCellsgrp)

index <- (nCells >= minCells) & (rowSums(nCellsgrp >= minCellsgrp) >= ncol(nCellsgrp))
sum(index)
#[1] 15733
#index <- index & (nCounts >= minCounts) & (nCounts <= maxCounts)
#sum(index)
#[1] 15033
index <- index & (nCounts >= minCounts)
sum(index)
#[1] 15050
index <- index & (cpc > mincpc)
sum(index)
#[1] 15033
#[1] 14175

counts <- counts[index, ] 
dim(counts)
#[1] 14727 27650
#[1] 11968 26202
#[1] 14175 27550
dim(coldata)

rm(index)
rm(dat)

all(colnames(counts) == rownames(coldata))

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

#Y <- counts[gset, ]
##raw counts
table(counts[1, ])
table(round(counts[1, ]))

Y <- round(counts)
dim(Y) 
table(Y[1,])

##log library size
loglib <- log(colSums(Y))
#hist(loglib)


##nebula
##fixed effect desigm matrix
#X <- model.matrix(~ Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata)
X <- model.matrix(~ loglib + Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata)
dim(X)
#[1] 26202    35
#[1] 27550    37
colnames(X) <- gsub("Cell_Types_Broad", "", colnames(X))
colnames(X) <- gsub("sex", "", colnames(X))
#colnames(X) <- gsub("Phase", "", colnames(X))

colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub("\\-", "_", colnames(X))
colnames(X) <- gsub(" ", ".", colnames(X))

head(X)
dim(X)

Z <- coldata$sampleID
table(Z)
length(Z) 


##model:
##NBGMM' is for fitting a negative binomial gamma mixed model. 
##'PMM' is for fitting a Poisson gamma mixed model. 
##'NBLMM' is for fitting a negative binomial lognormal mixed model (the same model as that in the lme4 package).

#model <- "NBLMM" #
model <- "NBGMM"

t1 <- Sys.time()
##The cells in the count matrix need to be grouped by the subjects
o <- order(Z)
negbn <- nebula(count = Y[, o], id = as.factor(Z[o]), pred = X[o, ], 
                model = model, cpc = mincpc, covariance = T, output_re = T)
	t2 <- Sys.time()
	difftime(t2, t1)
	#Time difference of 2.542573 hours
	
	
	timeUnits <- "secs"
	rtnebula <- difftime(t2, t1, units = timeUnits)
	rtnebula
	#Time difference of 9153.262 secs
	save(negbn, t1, t2, rtnebula, file = paste0(dirOut, "/kidney-counts-nebula.RData"))


table(negbn$convergence) 
#  -40   -25   -10     1 
#   55   516    75 13529 
      
# 1:   The convergence is reached due to a sufficiently small improvement of the function value.
#-10: The convergence is reached because the gradients are close to zero (i.e., the critical point) and no improvement of the function value can be found.
# -25: The Hessian matrix is either almost singular or not positive definite.
# -30: The convergence fails because the likelihood function returns NaN.
# -40: The convergence fails because the critical point is not reached and no improvement of the function value can be found.
   	


##################################################
##lmmfit
Y <- round(counts)

##log library size
loglib <- log(colSums(Y))

rm(counts)

##log-transformation
##log2(1 + counts)
##log2(1+Y) by groupping to reduce data size
ngrp <- 6
sizegrp <- round(nrow(Y)/ngrp)
if (ngrp*sizegrp < nrow(Y)) sizegrp <- round(nrow(Y)/ngrp) + 1
for (i in 1:ngrp){
  j <- (1+(i-1)*sizegrp):(min(nrow(Y), i*sizegrp))
  print(range(j))
  Y[j, ] <- log2(1 + Y[j, ])
}

table(Y[nrow(Y), ])
#               0                1 1.58496250072116                2 2.32192809488736 
#           26719              702              114               14                1 


##fixed effect desigm matrix
#X <- model.matrix(~ Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata)
X <- model.matrix(~ loglib + Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata)
colnames(X) <- gsub("Cell_Types_Broad", "", colnames(X))
colnames(X) <- gsub("sex", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub("\\-", "_", colnames(X))
colnames(X) <- gsub(" ", ".", colnames(X))

head(X)
dim(X)

##random effects
##sample groups
Z <- model.matrix(~ 0 + sampleID, data = coldata)
colnames(Z) <- gsub(".*sampleID", "", colnames(Z))
dim(Z)
head(Z)

d <- ncol(Z)

##########
##Fit LMM by lmmfit
epsilon <- 1e-5
max.iter <- 100

t1 <- Sys.time()
fit <- lmmfit(Y, X, Z, d=d, max.iter = max.iter, epsilon = epsilon)
	t2 <- Sys.time()
	rtlmm <- difftime(t2, t1, units = "secs")
	rtlmm
	#Time difference of 155.5957 secs

	#save(fit, epsilon, max.iter, rtlmm, file = paste0(dirOut, "/kidney-counts-lmmfit.RData"))
		
table(fit$niter)

#There were 23 warnings (use warnings() to see them)
warnings()

##########
n <- nrow(X)
XY <- t(Y%*%X)
ZX <- t(Z)%*%X #as.matrix(t(Z)%*%X)
ZY <- t(Y%*%Z) #as.matrix(t(Z)%*%Y)
ZZ <- t(Z)%*%Z #as.matrix(t(Z)%*%Z)

#XXinv <- ginv(t(X)%*%X)
XXinv <- chol2inv(chol(t(X)%*%X))
Ynorm <- rowSums(Y*Y) 


t1 <- Sys.time()
fitss <- lmm(XY, ZX, ZY, ZZ = ZZ, XXinv = XXinv, Ynorm = Ynorm, 
		n = n, d = d, max.iter = max.iter, epsilon = epsilon)
	t2 <- Sys.time()
	difftime(t2, t1)
	#Time difference of 14.9471 mins
	#rtlmm <- difftime(t2, t1, units = timeUnits) 

table(fitss$niter)
identical(fit, fitss)

rm(XY, ZX, ZY, ZZ, XXinv, Ynorm)


##################################################
##comparison of results

##########
##lmmfit

##convergence
table(fit$niter)
sum(fit$niter >= max.iter)
#[1] 4

##fixed effects
felmm <- t(fit$coef)
##variance components
slmm <- t(fit$theta)

##negative variance components
i <- apply(fit$theta < 0, 2, any)
sum(i)
#[1] 69
##p-values for testing variance components
pv.theta <- pnorm(fit$theta/fit$se.theta, lower.tail = F)
pv.theta[1, i]
sum(pv.theta[1,] < 0.05)
#[1] 12967
sum(pv.theta[1,] > 0.05)
#[1] 1208

hist(pv.theta[1, ])

#####
##LMM tests
#test <- lmmtest(fit)
##t-values
#tvlmm <- test[, grep("_t", colnames(test)), drop = F]
tvlmm <- t(fit$t)

dim(tvlmm)
head(tvlmm)
sum(is.na(tvlmm))
#[1] 28
sum(apply(is.na(tvlmm), 1, any))
#[1] 19
tvlmm[apply(is.na(tvlmm), 1, any), ]
hist(tvlmm)
hist(tvlmm[abs(tvlmm) < 5])

##p-values
#pvlmm <- test[, grep("_pvalue", colnames(test)), drop = F]
pvlmm <- t(fit$p)

dim(pvlmm)
head(pvlmm)
hist(pvlmm)
hist(pvlmm[, grep(":Male", colnames(pvlmm))])


##########
##nebula

rtnebula[[1]]/rtlmm[[1]]
#[1] 58.82723

#str(negbn)
table(negbn$convergence)
#  -40   -25   -10     1 
#   55   516    75 13529 
   
st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)

any(is.na(st))
sum(is.na(st))
#[1] 38184

##fixed effects, se, and p-values
iFC <- grep("logFC_", colnames(st))
ise <- grep("se_", colnames(st))
ipv <- grep("p_", colnames(st))

##fixed effects
b <- as.matrix(st[, iFC])
##se
se <- as.matrix(st[, ise])
##t-values
tv <- b/se
##p-values
pv <- as.matrix(st[, ipv])

range(pv - 2*pnorm(-abs(tv)), na.rm = T)
#[1] -2.720046e-15  3.136380e-15

##
hist(tv)
hist(tv[abs(tv) < 5])

hist(pv)
hist(pv[, grep(":Male", colnames(pv))])




##########
##comparison of lmmfit and nebula

all(rownames(tv) == rownames(tvlmm))

#####
##number of significant genes
##by FDR
colnames(pvlmm)
colnames(pv)
all(colnames(pvlmm) == gsub("p_", "", colnames(pv)))
j <- grep("Male", colnames(pv))

fdrlmm <- apply(pvlmm[, j], 2, p.adjust, method = "fdr")
head(fdrlmm)
fdrneb <- apply(pv[, j], 2, p.adjust, method = "fdr")
head(fdrneb)

fdrcut <- 0.05
fdrcut <- 0.01

ng <- rbind(lmmfit = colSums(fdrlmm < fdrcut, na.rm = T),
	nebula = colSums(fdrneb < fdrcut, na.rm = T))
o <- order(ng["lmmfit", ], decreasing = T)
ng <- ng[, o]

bp <- barplot(ng, beside = T, axes = T, xaxt = "n", col = c("red", "blue"), border = NA,
	ylim = c(0, 1.2*max(ng)), ylab = paste0("Number of genes with FDR < ", fdrcut))
text(colMeans(bp), par("usr")[3], labels = colnames(ng), srt = 45, adj = c(1.1, 1.1), xpd = TRUE, cex = 0.7)
legend("topright", rownames(ng), col = c("red", "blue"), pch = 15, bty = "n")
  
table(coldata$Cell_Types_Broad)


#####
##p-values
plot(as.matrix(-log10(pvlmm[, j])), as.matrix(-log10(pv[, j])), 
     xlab = "lmmfit -log10(p-values)", ylab = "nebula -log10(p-values)", cex = 0.6)
abline(0, 1)
	
	

#####
##histograms of p-values for male vs female
index <- grep("Male", colnames(pvlmm))
length(index)
par(mfrow = c(6, 3), mar = c(2.1,2.1,2.1,1.1))
for (j in index[o]){
  nm <- gsub("p", "+", gsub("_pvalue", "", colnames(pvlmm)[j]))
  hist(pvlmm[,j], xlab = "nebula p-values", main = nm, cex.main = 0.9, cex.axis = 0.8)
  }
 
 	
##histograms of p-values
index <- grep("Male", colnames(pv))
length(index)
par(mfrow = c(6, 3), mar = c(2.1,2.1,2.1,1.1))
for (j in index[o]){
  nm <- gsub("p", "+", gsub("p_", "", colnames(pv)[j]))
  hist(pv[,j], xlab = "nebula p-values", main = nm, cex.main = 0.9, cex.axis = 0.8)
  }
  
  

##histogram of p-values
par(mfrow = c(2, 1), mar = c(5.1, 4.1, 1.1, 2.1))
j <- grep("Male", colnames(tv))
hist(as.matrix(pvlmm[, j]), xlab = "lmmfit p-values for Male vs Female within cell-types", main = NULL)
hist(as.matrix(pv[, j]), xlab = "nebula p-values for Male vs Female within cell-types", main = NULL)


##histograms of p-values
index <- grep(":Male", colnames(pv))
length(index)
for (k in 1:3){
	
kj <- index[((k-1)*6+1):(k*6)]	
par(mfrow = c(3, 4))
for (j in kj){
  h1 <- hist(pv[,j], plot = F)
  h2 <- hist(pvlmm[,j], plot = F)
  ylim <- c(0, 1.1*max(h1$counts, h2$counts))
  nm <- gsub("p", "+", gsub("_pvalue", "", colnames(pvlmm)[j]))
  
  par(mar = c(5.1,4.1,3.1,0))
  #nm <- gsub("p_", "", colnames(pv)[j])
  hist(pv[,j], ylim = ylim, xlab = "nebula p-values", 
       main = nm, cex.main = 0.9, cex.axis = 0.8)
  
  par(mar = c(5.1,2.1,3.1,2.1))	
  #nm <- gsub("_pvalue", "", colnames(pvlmm)[j])
  hist(pvlmm[, j], ylim = ylim, xlab = "lmmfit p-value", ylab = NULL, 
       main = nm, cex.main = 0.9, cex.axis = 0.8, yaxt = "n")
}
	
	
}


#####
##all t-values 
##t-values
j <- 2:ncol(tv)
plot(as.matrix(tvlmm[, j]), as.matrix(tv[, j]),
     xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1, col = "gray")


##Male vs Female
j <- grep("Male", colnames(tv))

##t-values
plot(as.matrix(tvlmm[, j]), as.matrix(tv[, j]),
     xlab = "lmmfit t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1, col = "gray")
	


##t-values
##cell-type specific
index <- grep(":Male", colnames(pv))
length(index)
nc <- round(sqrt(length(index)))
nr <- ceiling(length(index)/nc)
par(mfrow = c(nr, nc), mar = c(2.1,2.1,1.1,1.1))
plot(0,0, type = "n", axes = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
text(0, 0, "t-values\nnebula vs lmmfit", cex = 1.2)
for (j in index){
  nm <- gsub("p", "+", gsub("_t", "", colnames(tvlmm)[j]))
  plot(tvlmm[,j], tv[,j], cex = 0.6, xlab = NA, ylab = NA, main = nm, cex.main = 0.8, cex.axis = 0.8)
  abline(0,1)
}
	

##################################################


