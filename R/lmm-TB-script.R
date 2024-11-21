##########
##Directories
##Outputs
today <- Sys.Date()
dirOut <- paste0(dirWork, "/Outputs/TB/")
dirFig <- paste0(dirOut, "/Figures")

##
library(MASS)
library(Matrix)
#library(dplyr)
library(Seurat)

library(lme4)
#library(lmerTest)
#library(nebula)

##the R functions
source("R/varest.R")
source("R/sstat.R")
source("R/lmm.R")
source("R/lmmtest.R")


##########
##data
datafile <- "scData/Nathan_NatImm_2021.rds"
dat <- readRDS(file = datafile)
#str(dat)
dat

An object of class Seurat 
33569 features across 500089 samples within 1 assay 
Active assay: originalexp (33569 features, 0 variable features)
 2 layers present: counts, data

##
#counts <- GetAssayData(object = dat, slot = "counts")
counts <- GetAssayData(dat, layer = 'counts', assay = 'originalexp')

dim(counts)
#[1]  33569 500089

##
coldata <- dat@meta.data

dim(coldata)
head(coldata)
all(rownames(coldata) == colnames(counts))
#[1] TRUE

table(coldata$donor)
length(unique(coldata$donor))
#[1] 259
table(coldata$sex)
#     F      M 
#280735 219354 
table(coldata$batch)
length(unique(coldata$batch))
#[1] 46

table(coldata$cluster_name)
length(table(coldata$cluster_name))
#[1] 29
#table(coldata$TB_status)
#   CASE CONTROL 
# 241814  258275 
          
          
##Filter cells
##by removing the cells beyond the extreme of lower and upper whisker.
##library size
libsize <- colSums(counts)
range(libsize)
#[1]    835 163673

libsbp <- boxplot(libsize)
libsbp$stats

sum(libsize > libsbp$stats[5])
minlibsize <- libsbp$stats[1]
maxlibsize <- 2^15 #libsbp$stats[5]
sum(libsize > maxlibsize)

j <- (libsize >= minlibsize) & (libsize <= maxlibsize)
sum(j)
#[1] 499973

##
counts <- counts[, j]
coldata <- coldata[j, ]
rm(j)

all(colnames(counts) == rownames(coldata))
rm(libsize)

##Filtering genes
##(1) number of celss
nCells <- rowSums(counts > 0)
#hist(log2(nCells))
minCells <- 2^9
sum(nCells >= minCells)

##(2) Filtering out low-expressed genes can be specified by cpc=0.005 (i.e., counts per cell<0.5%). 
##The argument cpc is defined by the ratio between the total count of the gene and the number of cells.

cpc <- rowSums(counts)/ncol(counts)
range(cpc[cpc > 0])
#minCountsgrp
mincpc <- 0.005

index <- (nCells >= minCells) & (cpc > mincpc)

counts <- counts[index, ] 
dim(counts)
#[1]  11596 499973
dim(coldata)

coldata$loglib <- log(colSums(counts))
hist(coldata$loglib)

#save(coldata, file = "/scData/Nathan_NatImm_coldata.RData")
##################################################


##################################################
##LMM

##log-transformation
##log2(1 + counts)
##Y = log2(1+Y) 

Y <- counts
dim(Y)
rm(counts)


##by groupping to reduce data size
ngrp <- 50
sizegrp <- round(nrow(Y)/ngrp)
sizegrp
if (ngrp*sizegrp < nrow(Y)) sizegrp <- round(nrow(Y)/ngrp) + 1
for (i in 1:ngrp){
  j <- (1+(i-1)*sizegrp):(min(nrow(Y), i*sizegrp))
  Yj <- log2(1 + Y[j, ])
  print(range(j))
  #Y[j, ] <- log2(1 + Y[j, ])

save(j, Yj, file = "/scData/Nathan_NatImm_Y", i, ".RData")
rm(Yj)
}


##########
##design matrices
##load meta data
load(file = "/scData/Nathan_NatImm_coldata.RData")
	dim(coldata)
	#[1] 499973     19
	sum(is.na(coldata))
	#[1] 520

##remove NAs
if (any(is.na(coldata$cluster_name))){
	#Yj <- Yj[, !(is.na(coldata$cluster_name))]
	iNA <- (is.na(coldata$cluster_name))
	coldata <- coldata[!iNA,]
}

sum(iNA)

table(coldata$season)
#Spring Summer Winter 
# 69324 186527 243862 
table(coldata[, c("season", "TB_status")]) 
table(coldata[, c("sex", "TB_status")]) 
table(coldata[, c("season", "sex")]) 
 
#####
##random effects 

if (length(grep("batch", coldata$batch)) == 0) coldata$batch <- paste0("batch", coldata$batch)
table(coldata$batch)

Zd <- model.matrix(~ 0 + donor, data = coldata)
#Zb <- model.matrix(~ 0 + batch, data = coldata)

Z <- Zd
d <- ncol(Zd)
dim(Z)
d

rm(Zd)

#####
##fixed effects
table(coldata$cluster_name)
coldata$cluster_name <- gsub(" ", ".", coldata$cluster_name)
coldata$cluster_name <- gsub("\\+", "p", coldata$cluster_name)
coldata$cluster_name <- gsub("\\-", "_", coldata$cluster_name)
coldata$cluster_name <- gsub("\\/", "s", coldata$cluster_name)
table(coldata$cluster_name)
length(table(coldata$cluster_name))

table(coldata$TB_status)
##change levels (orders)
coldata$TB_status[grep("CASE", coldata$TB_status)] <- "trt"
coldata$TB_status[grep("CONTROL", coldata$TB_status)] <- "ctr"
	table(coldata$TB_status)
	#   ctr    trt 
	#258225 241748 
	#258084 241629 

##normalized the age
hist(coldata$age)
rk <- rank(coldata$age)
age <- qnorm(rk/(max(rk)+1))
hist(age)
#plot(age, coldata$age)
coldata$age <- age
rm(age)

##age, sex and season
X <- model.matrix(~ loglib + age + sex + season + cluster_name + cluster_name:TB_status, data = coldata)
	
colnames(X) <- gsub("cluster_name", "", colnames(X))
colnames(X) <- gsub("TB_status", "", colnames(X))
#head(X)
colnames(X)

rm(coldata)

#####
##computing summary statistics

XY <- NULL
ZY <- NULL
Ynorm <- NULL
gindex <- NULL

t1 <- Sys.time()
for (i in 1:50){
##load the log-transformed data
load(file = paste0("scData/Nathan_NatImm_Y", i, ".RData"))
	gindex <- c(gindex, j)
	#dim(Yj)
	#[1]    232 499973
	Yj <- Yj[, !iNA]
	#stopifnot(all(colnames(Yj) == rownames(coldata)))
	
	XY <- cbind(XY, t(Yj%*%X))
	ZY <- cbind(ZY, t(Yj%*%Z))
	Ynorm <- c(Ynorm, rowSums(Yj*Yj))
	
	rm(j, Yj)
}
##
n <- nrow(X)
XX <- t(X)%*%X
ZX <- t(Z)%*%X #as.matrix(t(Z)%*%X)
ZZ <- t(Z)%*%Z #as.matrix(t(Z)%*%Z)

t2 <- Sys.time()
rtSS <- difftime(t2, t1)
rtSS

rm(X, Z)


#####
epsilon <- 1e-5
max.iter <- 50

t1 <- Sys.time()
fit <- lmm(XX, XY, ZX, ZY, ZZ, Ynorm = Ynorm, n = n, d = d, max.iter = max.iter, epsilon = epsilon)
#fitss <- lmm(summary.stats = ss, d = d, max.iter = max.iter, epsilon = epsilon)
	t2 <- Sys.time()
	rtlmm <- difftime(t2, t1)
	rtlmm

##################################################


##################################################
##lmer 

##design matrices
##load meta data
load(file = paste0("scData/Nathan_NatImm_coldata.RData")

##Relabel the batches as a charactor instead of an integer.
if (length(grep("batch", coldata$batch)) == 0) coldata$batch <- paste0("batch", coldata$batch)
table(coldata$batch)

##remove NAs
dim(coldata)
if (any(is.na(coldata$cluster_name))){
	iNA <- (is.na(coldata$cluster_name))
	coldata <- coldata[!iNA,]
	}

sum(iNA)
dim(coldata)

##fixed effects
table(coldata$cluster_name)
coldata$cluster_name <- gsub(" ", ".", coldata$cluster_name)
coldata$cluster_name <- gsub("\\+", "p", coldata$cluster_name)
coldata$cluster_name <- gsub("\\-", "_", coldata$cluster_name)
coldata$cluster_name <- gsub("\\/", "s", coldata$cluster_name)
table(coldata$cluster_name)
length(table(coldata$cluster_name))

table(coldata$TB_status)
##change levels (orders)
coldata$TB_status[grep("CASE", coldata$TB_status)] <- "trt"
coldata$TB_status[grep("CONTROL", coldata$TB_status)] <- "ctr"
	table(coldata$TB_status)
	#   ctr    trt 
	#258084 241629 


##lmer with default setting
#modelformula <- formula(y ~ loglib + cluster_name + cluster_name:TB_status + (1 | donor) + (1 | batch))
modelformula <- formula(y ~ loglib + cluster_name + cluster_name:TB_status + (1 | donor))

gnm <- NULL
felme <- NULL
tvlme <- NULL
plme <- NULL
slme <- NULL
rtlme <- 0
timeUnits <- "secs"

convlme <- NULL
singlme <- NULL

#####
##loop for subsets of genes
for (i in 1:50){
##load the log-transformed data
load(file = paste0("scData/Nathan_NatImm_Y", i, ".RData"))
Y <- Yj
rm(Yj)
Y <- Y[, !(iNA)]

stopifnot(all(colnames(Y) == rownames(coldata)))

gnm <- c(gnm, rownames(Y))

print(i)
print(Sys.time())

##loop for each gene
for (j in 1:nrow(Y)){

##the observations
coldata$y <- Y[j, ]

##lmer with default setting
##LMM fitting: lme4::lmer 
t1 <- Sys.time()
lmerFit <- lme4::lmer(modelformula, data = coldata)
	t2 <- Sys.time()
	rtlme <- rtlme + difftime(t2, t1, units = timeUnits) 
	rm(t1, t2)
convlme <- c(convlme, lmerFit@optinfo$conv$opt)
#convlme <- c(convlme, length(lmerFit@optinfo$conv$lme4))
singlme <- c(singlme, isSingular(lmerFit))

sfit <- summary(lmerFit)$coefficients
felme <- cbind(felme, sfit[, "Estimate"])
tvlme <- cbind(tvlme, sfit[, "t value"])
slme <- cbind(slme, as.data.frame(VarCorr(lmerFit))$vcov)

}##j

print(Sys.time())

save(felme, tvlme, plme, slme, gnm, rtlme, convlme, singlme, file = paste0(dirOut, "/TB_lmer_default.RData"))
}##i

colnames(felme) <- gnm
colnames(tvlme) <- gnm
#colnames(plme) <- gnm
colnames(slme) <- gnm
save(felme, tvlme, plme, slme, rtlme, convlme, singlme, file = paste0(dirOut, "/TB_lmer_default.RData"))

##################################################


##################################################
##comparison of results


#####
##lmmfit
load(paste0(dirOut, "/TB_lmm_intercept.RData"))

rtlmm; max.iter; epsilon

##convergence
table(fit$niter)
sum(fit$niter >= max.iter)

##fixed effects
felmm <- t(fit$coef)
##variance components
slmm <- t(fit$theta)

##negative variance components
i <- apply(fit$theta < 0, 2, any)
sum(i)

##LMM tests
#test <- lmmtest(fit)
##t-values
#tvlmm <- test[, grep("_t", colnames(test)), drop = F]
tvlmm <- t(fit$t)

##p-values
#pvlmm <- test[, grep("_pvalue", colnames(test)), drop = F]
pvlmm <- t(fit$p)


#####
##lmer
load(file = paste0(dirOut, "/TB_lmer_default_intercept.RData"))

rtlmm
#Time difference of 1.395676 hours
rtlme
#Time difference of 200331.7 secs
rtlme[[1]]/(3600*rtlmm[[1]])
#[1] 39.8715

rtime <- cbind("lmm" = rtlmm[[1]], "lmmfit" = rtlmm[[1]] + rtSS[[1]]/60, "lmer" = rtlme[[1]]/3600)
rtime

#####
ba <- barplot(rtime, beside = T, ylim = c(0, 1.05*max(rtime)),
	ylab = "Computation time (hours)", col = "gray")
text(ba, rtime, round(rtime[1,], 1), pos = 3)


#####
##the genes that 
##(1) lmer fitting is not boundary (singular).
##(2) lmer fitting is convergent.

##not boundary (singular)
i <- (!singlme)
sum(i)
sum(!i)

##convergence
table(convlme)
sum(convlme == 0)

iconv <- ((!singlme) & (convlme == 0))
sum(iconv)

##lmmfit convergence
ipos <- (apply(slmm > 0, 1, all) & (fit$niter < max.iter))
sum(ipos)

#####
##differences between lmm and lmer

all(rownames(felmm) == colnames(felme))
range(felmm - t(felme))
#[1] -9.211326e-07  2.303989e-07

range(t(tvlmm) - tvlme)
#[1] -0.001447692  0.000317034
nx <- (nrow(felme) - 1)/2
range(t(tvlmm[, -(1:(nx+1))]) - tvlme[-(1:(nx+1)), ])
#[1] -2.837501e-05  3.123577e-05

range(t(slmm) - slme)
#[1] -2.486041e-06  4.338225e-07


##
d <- list('Coefficients' = c(t(felmm) - felme), 
	"t-values" = t(tvlmm) -tvlme, 
	'Variance components' = c(t(slmm) - slme))
boxplot(d, ylab = "Differences between lmm and lmer", cex.axis = 0.9, cex.lab = 0.9)

##
d <- list('Coefficients' = c(t(felmm) - felme), 
	'Variance components' = c(t(slmm) - slme))
boxplot(d, ylab = "Differences between lmm and lmer", cex.axis = 0.9, cex.lab = 0.9)


##########
##DE genes
dim(pvlmm)
dim(tvlmm)
dim(felmm)
identical(dimnames(pvlmm), dimnames(felmm))

icell <- grep(":", colnames(pvlmm))
length(icell)
pv <- pvlmm[, icell]
fe <- felmm[, icell]

hist(pv)

par(mfrow = c(6,5), mar = c(0,0,0,0))
for (i in 1:29){
	hist(pv[,i], axes = F, xlab = NULL, ylab = NULL, main = NA)
	mtext(colnames(pv)[i], line = -2, cex = 0.6, col = "blue")
}
	
	
##FDR
fdr <- apply(pv, 2, p.adjust, method = "fdr")
head(fdr)

##Numbers of DE genes
#colSums(fdr <= 0.05 & fe > 0.1)
colSums(fdr <= 0.05 & fe > 0)

##################################################


