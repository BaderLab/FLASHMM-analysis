

##################################################
##load data
#datafile <- "scData/Human_Kidney_data.rds"
datafile <- "scData/Kidney_raw_counts_decontaminated.rds"
dat <- readRDS(file = datafile)
#str(dat)

#counts <- GetAssayData(object = dat, slot = "counts")
counts <- GetAssayData(dat, layer = 'counts', assay = 'RNA')

dim(counts)
#[1] 22484 27677
max(counts)
#[1] 25332.21

#####
coldata <- dat@meta.data
dim(coldata)
head(coldata)
all(rownames(coldata) == colnames(counts))
#[1] TRUE


##'sampleID': unique sample identifiers
##'Cell_Types_Broad': subpopulation (cell cluster) assignments 
##'sex': make vs female experimental group/condition (control/treatment, healthy/diseased)

table(coldata$sampleID)

CD45_1 CD45_10  CD45_2  CD45_3  CD45_4  CD45_5  CD45_6  CD45_7  CD45_8  CD45_9 
    584     714     315     938     126     344    1102    1182     553    1041 
 Total1  Total2  Total3  Total4  Total5  Total6  Total7  Total8  Total9 
   1080    4738    2366    3739    3280     606    1811    1894    1264 

table(coldata$sex)
# Female   Male 
# 15945  11732 
table(coldata$sampletype)
#  LD 
#6899 
sort(table(coldata$Cell_Types_Broad))

       Podocyte        CCD-like             PEC          B cell              U2             DCT 
             16              45              64              68             107             109 
            CNT            IC-B       Mesangial              U1            IC-A              PC 
            122             130             130             214             304             455 
        NK cell     Endothelial          T cell             MNP        LOH-like            cTAL 
            466             572             938             999            1011            1156 
Proximal Tubule 
          20771 
         


#####
##Filter cells
##by removing the cells beyond the extreme of lower and upper whisker.

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


##Filtering
minCells <- 16 
minCellsgrp <- 10 
minCounts <- 2^6
maxCounts <- 2^20
mincpc <- 0.005

length(nCells)
dim(nCellsgrp)

index <- (nCells >= minCells) & (rowSums(nCellsgrp >= minCellsgrp) >= ncol(nCellsgrp))
index <- index & (nCounts >= minCounts)
index <- index & (cpc > mincpc)
sum(index)

counts <- counts[index, ] 
dim(counts)
dim(coldata)

all(colnames(counts) == rownames(coldata))
##################################################



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

Y <- round(counts)
dim(Y) 
table(Y[1,])

##log library size
loglib <- log(colSums(Y))
#hist(loglib)

##
rm(counts)

##nebula
##fixed effect desigm matrix
X <- model.matrix(~ loglib + Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata)
dim(X)

colnames(X) <- gsub("Cell_Types_Broad", "", colnames(X))
colnames(X) <- gsub("sex", "", colnames(X))
#colnames(X) <- gsub("Phase", "", colnames(X))

colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub("\\-", "_", colnames(X))
colnames(X) <- gsub(" ", ".", colnames(X))


Z <- coldata$sampleID
table(Z)
length(Z) 


##model:
##NBGMM' is for fitting a negative binomial gamma mixed model. 
##'PMM' is for fitting a Poisson gamma mixed model. 
##'NBLMM' is for fitting a negative binomial lognormal mixed model (the same model as that in the lme4 package).

#model <- "NBLMM"
model <- "NBGMM"

t1 <- Sys.time()
##The cells in the count matrix need to be grouped by the subjects
o <- order(Z)
negbn <- nebula(count = Y[, o], id = as.factor(Z[o]), pred = X[o, ], 
                model = model, cpc = mincpc, covariance = T, output_re = T)
	t2 <- Sys.time()
	difftime(t2, t1)
	
	timeUnits <- "secs"
	rtnebula <- difftime(t2, t1, units = timeUnits)
	rtnebula
	##Time difference of 8986.294 secs
	save(negbn, t1, t2, rtnebula, model, file = paste0(dirOut, "/kidney-counts-nebula.", model, ".RData"))


table(negbn$convergence) 
          
# 1:   The convergence is reached due to a sufficiently small improvement of the function value.
#-10: The convergence is reached because the gradients are close to zero (i.e., the critical point) and no improvement of the function value can be found.
# -25: The Hessian matrix is either almost singular or not positive definite.
# -30: The convergence fails because the likelihood function returns NaN.
# -40: The convergence fails because the critical point is not reached and no improvement of the function value can be found.

##################################################


##################################################
##lmmfit
Y <- round(counts)
dim(Y)

##log library size
loglib <- log(colSums(Y))

#####
##log-transformation
##log2(1 + counts)
##Y = log2(1+Y) 
##by groupping to reduce data size
ngrp <- 6
sizegrp <- round(nrow(Y)/ngrp)
if (ngrp*sizegrp < nrow(Y)) sizegrp <- round(nrow(Y)/ngrp) + 1
for (i in 1:ngrp){
  j <- (1+(i-1)*sizegrp):(min(nrow(Y), i*sizegrp))
  print(range(j))
  Y[j, ] <- log2(1 + Y[j, ])
}


#####
X <- model.matrix(~ loglib + Cell_Types_Broad + Cell_Types_Broad:sex, data = coldata)

dim(X)

colnames(X) <- gsub("Cell_Types_Broad", "", colnames(X))
colnames(X) <- gsub("sex", "", colnames(X))
#colnames(X) <- gsub("Phase", "", colnames(X))
colnames(X) <- gsub("\\+", "p", colnames(X))
colnames(X) <- gsub("\\-", "_", colnames(X))
colnames(X) <- gsub(" ", ".", colnames(X))

##random effects
##sample groups
Z <- model.matrix(~ 0 + sampleID, data = coldata)
colnames(Z) <- gsub(".*sampleID", "", colnames(Z))
dim(Z)

d <- ncol(Z)

#####
epsilon <- 1e-5
max.iter <- 100

t1 <- Sys.time()
fit <- lmmfit(Y, X, Z, d=d, max.iter = max.iter, epsilon = epsilon)
	t2 <- Sys.time()
	difftime(t2, t1)

	rtlmm <- difftime(t2, t1, units = "secs")
	rtlmm
	#Time difference of 92.03347 secs
	save(fit, epsilon, max.iter, rtlmm, file = paste0(dirOut, "/kidney-counts-lmm.beta.RData"))
		
table(fit$niter)

#####
##lmm

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
	#Time difference of 1.111575 mins

identical(fit, fitss)

##################################################



##################################################
##lmer 
##
dat <- coldata[, c("Cell_Types_Broad", "sex", "sampleID")]
table(dat$Cell_Types_Broad)
dat$Cell_Types_Broad <- gsub("\\+", "p", dat$Cell_Types_Broad)
dat$Cell_Types_Broad <- gsub("\\-", "_", dat$Cell_Types_Broad)
dat$Cell_Types_Broad <- gsub(" ", ".", dat$Cell_Types_Broad)

modelformula <- formula(y ~ loglib + Cell_Types_Broad + Cell_Types_Broad:sex + (1 | sampleID))

##lmer with default setting
felme <- NULL
tvlme <- NULL
plme <- NULL
slme <- NULL
rtlme <- 0
timeUnits <- "secs"

convlme <- NULL
singlme <- NULL

for (j in 1:nrow(Y)){

##the observations
dat$y <- Y[j, ]

##lmer with default setting
##LMM fitting: lme4::lmer 
t1 <- Sys.time()
lmerFit <- lme4::lmer(modelformula, data = dat)
#lmerFit <- lme4::lmer(modelformula, data = dat, control = control)
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

##LMM testing: lmerTest::lmer
#lmerFit <- lmerTest::lmer(modelformula, data = dat)
#lmerFit <- lmerTest::lmer(modelformula, data = dat, control = control)
#sfit <- summary(lmerFit)$coefficients
#plme <- cbind(plme, sfit[, "Pr(>|t|)"])

#rm(sfit, lmerFit)
}

##
colnames(felme) <- colnames(Y)
colnames(tvlme) <- colnames(Y)
colnames(plme) <- colnames(Y)
colnames(slme) <- colnames(Y)

save(felme, tvlme, plme, slme, rtlme, convlme, singlme, 
	file = paste0(dirOut, "/kidney-counts_lmer_default.RData"))

##
rtlme
#Time difference of 7166.913 secs

##################################################




##################################################
##comparison of results

##########
##lmmfit
load(paste0(dirOut, "/kidney-counts-lmm.beta.RData"))

rtlmm; max.iter; epsilon

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


##p-values
#pvlmm <- test[, grep("_pvalue", colnames(test)), drop = F]
pvlmm <- t(fit$p)

dim(pvlmm)
head(pvlmm)


##########
##lmer
load(file = paste0(dirOut, "/kidney-counts_lmer_default.RData"))

rtlmm
Time difference of 92.03347 secs
rtlme
Time difference of 7166.913 secs
rtlme[[1]]/rtlmm[[1]]
[1] 77.87289


##the genes that 
##(1) lmer fitting is not boundary (singular).
##(2) lmer fitting is convergent.

##not boundary (singular)
i <- (!singlme)
sum(i)
#[1] 14106
sum(!i)
#[1] 69

##convergence
table(convlme)
sum(convlme == 0)

iconv <- ((!singlme) & (convlme == 0))
sum(iconv)
#[1] 14106

##lmmfit convergence
ipos <- (apply(slmm > 0, 1, all) & (fit$niter < max.iter))
sum(ipos)
#[1] 14106
sum(ipos & iconv)


##
i <- (iconv & ipos)
sum(i)

range(felmm[i, ] - t(felme[, i]))
#[1] -1.570977e-05  5.202275e-06
range(slmm[i, ] - t(slme[, i]))
#[1] -2.252868e-06  3.987087e-06
range(t(tvlmm[i, ]) - tvlme[, i])
#[1] -0.0013624379  0.0006275755

#range(t(pvlmm[i,]) - plme[, i])  ##lmerTest::lmer

#####
i <- (iconv & ipos)
sum(i)
d <- list('Coefficients' = c(t(felmm[i,]) - felme[, i]), 
	"t-values" = t(tvlmm[i,]) - tvlme[, i], 
	'Variance components' = c(t(slmm[i, ]) - slme[, i]))
boxplot(d, ylab = "Differences between lmm and lmer", cex.axis = 0.9, cex.lab = 0.9)
	

##########
##nebula
#load(paste0(dirOut, "/kidney-counts-nebula.NBLMM.RData"))
#nblmm <- negbn
#rtnblmm <- rtnebula

load(paste0(dirOut, "/kidney-counts-nebula.NBGMM.RData"))
nbgmm <- negbn

rtlmm
Time difference of 92.03347 secs
rtlme
Time difference of 7166.913 secs
rtlme[[1]]/rtlmm[[1]]
[1] 77.87289

rtlme[[1]]/rtlmm[[1]]
#[1] 77.87289
rtnebula[[1]]/rtlmm[[1]]
#[1] 99.45579
rtnblmm[[1]]/rtlmm[[1]]
#[1] 412.11

##lmm:
#Time difference of 1.111575 mins

#####
rtime <- cbind("lmm" = 1.111575, "lmmfit" = 92.03347/60, "lmer" = 7166.913/60)
rtime

ba <- barplot(rtime, beside = T, ylim = c(0, 1.05*max(rtime)),
	ylab = "Computation time (minutes)", col = "gray")
text(ba, rtime, round(rtime[1,], 1), pos = 3)


#####
#rtime <- cbind("lmm" = 1.111575, "lmmfit" = rtlmm/60, "lmer" = rtlme/60, 
#	"nebula\nNBGMM" = rtnebula/60, "nebula\nNBLMM" = rtnblmm/60)

rtime <- cbind("lmm" = 1.111575, "lmmfit" = rtlmm/60, "lmer" = rtlme/60, "nebula" = rtnebula/60)
rtime <- as.matrix(rtime)
rtime

ba <- barplot(rtime, beside = T, ylim = c(0, 1.05*max(rtime)),
	ylab = "Computation time (minutes)", col = "gray")
text(ba, rtime, round(rtime[1,], 1), pos = 3)


#####
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
#ipv <- setdiff(grep("p_", colnames(st)), c(iFC, ise))
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

##NBLMM
#st <- nblmm$summary
#rownames(st) <- st$gene
#dim(st)
##NBGMM vs NBLMM
#i <- (negbn$convergence >= -10) & (nblmm$convergence >= -10)
#plot(b[i, ], as.matrix(st[i, iFC])); abline(0, 1)
#plot(se[i, ], as.matrix(st[i, ise])); abline(0, 1)




##########
##comparison of lmmfit and nebula
##for the genes that lmmfit and nebula are convergent.
##https://github.com/lhe17/nebula
##Based on the empirical simulation study in (He et al. 2021), genes with an estimated cell-level overdispersion >100 should be removed for a data set with at least 50 cells per subject.
##If the variable of interest is subject-level, genes with a very large subject-level overdispersion (>1) should be removed or interpreted cautiously as well.

all(rownames(tv) == rownames(tvlmm))

##genes with convergence
iGenes <- (negbn$convergence >= -10)
sum(iGenes)
#[1] 13604
##overdispersion
iGenes <- iGenes & (negbn$overdispersion$Subject <= 1) & (negbn$overdispersion$Cell < 100) 
sum(iGenes)
#[1] 13586
##gens without NA
sum(apply(!is.na(tv[iGenes, ]), 1, all)) 

##negative variance components
sum(apply(slmm[iGenes, ] < 0, 1, any))
ipos <- (apply(slmm > 0, 1, all))
sum(ipos)
#[1] 14106

#####
##all t-values 
##t-values
j <- 2:ncol(tv)
j <- 3:ncol(tv)
i <- (iGenes & ipos)

plot(as.matrix(tvlmm[i, j]), as.matrix(tv[i, j]),
     xlab = "lmm t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1, col = "gray")
	
##Male vs Female
j <- grep("Male", colnames(tv))
j <- grep(":Male", colnames(tv))
j


#xylim <- range(tvlmm[i, j], tv[i, j], na.rm = T)
plot(as.matrix(tvlmm[i, j]), as.matrix(tv[i, j]), xlab = "lmm t-values", ylab = "nebula t-values", cex = 0.6)
abline(0, 1, col = "gray")
	
	
##p-values
hist(pvlmm[i,j], xlab = "lmm p-values", main = "")

hist(pv[i,j], xlab = "nebula p-values", main = "")


#####
##number of significant genes
##by Bonferroni Correction 
sum(i)

colSums(pvlmm[i,j] <= 0.05/sum(i))
         B.cell:Male        CCD_like:Male             CNT:Male            cTAL:Male 
                  48                   25                  646                  523 
            DCT:Male     Endothelial:Male            IC_A:Male            IC_B:Male 
                 119                  283                  235                  148 
       LOH_like:Male       Mesangial:Male             MNP:Male         NK.cell:Male 
                  10                  175                  522                   66 
             PC:Male             PEC:Male Proximal.Tubule:Male          T.cell:Male 
                 536                  255                   19                   58 
             U1:Male              U2:Male 
                 316                  203 

colSums(pv[i,j] <= 0.05/sum(i))
         p_B.cell:Male        p_CCD_like:Male             p_CNT:Male            p_cTAL:Male 
                     1                      0                      1                     31 
            p_DCT:Male     p_Endothelial:Male            p_IC_A:Male            p_IC_B:Male 
                     1                     11                      4                      2 
       p_LOH_like:Male       p_Mesangial:Male             p_MNP:Male         p_NK.cell:Male 
                     3                      4                     38                      1 
             p_PC:Male             p_PEC:Male p_Proximal.Tubule:Male          p_T.cell:Male 
                    13                      0                     71                     11 
             p_U1:Male              p_U2:Male 
                     6                      6 

##effects
i <- ipos
colSums((pvlmm[i,j] <= 0.05/sum(i)) & (abs(felmm[i,j]) > 0.25))
         B.cell:Male        CCD_like:Male             CNT:Male            cTAL:Male 
                  31                   24                  456                   97 
            DCT:Male     Endothelial:Male            IC_A:Male            IC_B:Male 
                  42                   49                   87                   63 
       LOH_like:Male       Mesangial:Male             MNP:Male         NK.cell:Male 
                   2                   68                  192                   29 
             PC:Male             PEC:Male Proximal.Tubule:Male          T.cell:Male 
                 145                  153                    2                   13 
             U1:Male              U2:Male 
                 123                  117 

i <- iGenes
colSums((pv[i,j] <= 0.05/sum(i)) & (abs(b[i,j]) > 0.25))
         p_B.cell:Male        p_CCD_like:Male             p_CNT:Male            p_cTAL:Male 
                     1                      0                      1                     31 
            p_DCT:Male     p_Endothelial:Male            p_IC_A:Male            p_IC_B:Male 
                     1                     11                      4                      2 
       p_LOH_like:Male       p_Mesangial:Male             p_MNP:Male         p_NK.cell:Male 
                     3                      4                     38                      1 
             p_PC:Male             p_PEC:Male p_Proximal.Tubule:Male          p_T.cell:Male 
                    13                      0                     71                     11 
             p_U1:Male              p_U2:Male 
                     6                      6 


#####
##number of significant genes
##by FDR
i <- ipos
fdr <- apply(pvlmm[i, j], 2, p.adjust, method = "BH")
colSums((fdr <= 0.05) & (abs(felmm[i,j]) > 0.25))
         B.cell:Male        CCD_like:Male             CNT:Male            cTAL:Male 
                  69                   68                 1391                  131 
            DCT:Male     Endothelial:Male            IC_A:Male            IC_B:Male 
                 147                   79                  147                  149 
       LOH_like:Male       Mesangial:Male             MNP:Male         NK.cell:Male 
                   2                  103                  272                   36 
             PC:Male             PEC:Male Proximal.Tubule:Male          T.cell:Male 
                 239                  448                   10                   19 
             U1:Male              U2:Male 
                 326                  264 
colSums((fdr <= 0.05) & (abs(felmm[i,j]) > 0.5))
         B.cell:Male        CCD_like:Male             CNT:Male            cTAL:Male 
                  20                   19                  200                   20 
            DCT:Male     Endothelial:Male            IC_A:Male            IC_B:Male 
                   9                    8                   10                   14 
       LOH_like:Male       Mesangial:Male             MNP:Male         NK.cell:Male 
                   0                   18                   59                    5 
             PC:Male             PEC:Male Proximal.Tubule:Male          T.cell:Male 
                  31                   53                    4                    5 
             U1:Male              U2:Male 
                  46                   62 
                  
colSums((fdr <= 0.05) & (abs(felmm[i,j]) >= 1))

##
i <- iGenes
fdr <- apply(pv[i, j], 2, p.adjust, method = "BH")
colSums((fdr <= 0.05) & (abs(b[i,j]) > 0.25))
         p_B.cell:Male        p_CCD_like:Male             p_CNT:Male            p_cTAL:Male 
                     1                      0                      1                     90 
            p_DCT:Male     p_Endothelial:Male            p_IC_A:Male            p_IC_B:Male 
                     1                     36                      4                      5 
       p_LOH_like:Male       p_Mesangial:Male             p_MNP:Male         p_NK.cell:Male 
                     5                      7                    166                      1 
             p_PC:Male             p_PEC:Male p_Proximal.Tubule:Male          p_T.cell:Male 
                    28                      0                    428                     36 
             p_U1:Male              p_U2:Male 
                    12                      7 
colSums((fdr <= 0.05) & (abs(b[i,j]) > 0.5))
         p_B.cell:Male        p_CCD_like:Male             p_CNT:Male            p_cTAL:Male 
                     1                      0                      1                     79 
            p_DCT:Male     p_Endothelial:Male            p_IC_A:Male            p_IC_B:Male 
                     1                     34                      4                      5 
       p_LOH_like:Male       p_Mesangial:Male             p_MNP:Male         p_NK.cell:Male 
                     5                      7                    140                      1 
             p_PC:Male             p_PEC:Male p_Proximal.Tubule:Male          p_T.cell:Male 
                    25                      0                    261                     34 
             p_U1:Male              p_U2:Male 
                    11                      7 
                    
                    
                    
#####
##DE genes in Proximal.Tubule
cls <- "Proximal.Tubule"
fdrcut <- 0.01

##lmmfit
p <- pvlmm
p <- p[, grep("Male", colnames(p))]
dim(p)
fdr <- p.adjust(p, method = "fdr")
fdr <- matrix(fdr, nrow = nrow(p), dimnames = dimnames(p))
glmm <- rownames(fdr)[which(fdr[, grep(cls, colnames(fdr))] <= fdrcut)]
length(glmm)

##nebula
p <- pv
p <- p[, grep("Male", colnames(p))]
dim(p)
fdr <- p.adjust(p, method = "fdr")
fdr <- matrix(fdr, nrow = nrow(p), dimnames = dimnames(p))
gneb <- rownames(fdr)[which(fdr[, grep(cls, colnames(fdr))] <= fdrcut)]
length(gneb)

##Venn diagram of top genes
#library(limma)

glist <- rownames(fdr) #union(glmm, gneb)
length(glist)
indx <- matrix(0, nrow = length(glist), ncol = 2)
	rownames(indx) <- glist
	colnames(indx) <- c("nebula", "lmm")
	indx[gneb, 1] <- 1
	indx[glmm, 2] <- 1		
	cgset <- vennCounts(indx)
	#main <- gsub("p_", "", colnames(pv)[j])
	vennDiagram(cgset, counts.col = "blue", cex = c(0.9, 0.9, 0.8), lwd = 1.5)
	##vennDiagram(cgset, counts.col = "blue", cex = c(0.9, 0.9, 0.8), lwd = 1.5, mar = rep(0, 4))
	#mtext(main, line = -5)


##################################################


