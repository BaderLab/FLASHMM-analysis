##################################################
source("R/qqpvalue.R")
source("/R/AUC.R")
source("/R/save.figure.R")

resFig <- 720
scaFig <- 7
widFig <- scaFig*resFig
heiFig <- scaFig*resFig

##folder of simulation results
dirOut <- "Simulations"


##################################################
##scRNA-seq simulation results

##differences between lmm and lmer
vd <- NULL
fd <- NULL
td <- NULL
pd <- NULL
logp12d <- NULL

runtime <- NULL

##lmm and nebula (nbn)
##t-values, p-values, and AUC
tlmm <- NULL
tnbn <- NULL
plmm <- NULL
pnbn <- NULL
almm <- NULL #AUC
anbn <- NULL #AUC

##FPR and TPR
FPRlmm <- NULL
FPRnbn <- NULL
TPRlmm <- NULL
TPRnbn <- NULL


#####
Ng <- 6000
#snm <- "beta0.5"
snm <- ""
direct <- "" #"dn" #"up"

#NcList <- c(4e4, 6e4, 8e4, 1e5)
NcList <- c(2e4, 4e4, 6e4, 8e4, 10e4, 12e4)

QQp <- T
if (QQp) par(mfrow = c(3, 2), mar = c(4.5, 4.5, 1, 1.1))

options(scipen=0)

for (Nc in NcList){

##plme
load(file = paste0(dirOut, "/simuNBMM", direct, "_Ng", Ng, "Nc", Nc, snm, "_lmerpv.RData"))
rtlme1 <- rtlme
rm(rtlme, convlme, singlme)

load(file = paste0(dirOut, "/simuNBMM", direct, "_Ng", Ng, "Nc", Nc, snm, "_lmer.RData"))
a <- NBMM.par

load(file = paste0(dirOut, "/simuNBMM", direct, "_Ng", Ng, "Nc", Nc, snm, "_lmm.RData"))
b <- NBMM.par
b$lmm.epsilon <- NULL
print(identical(a, b))
#[1] TRUE

load(file = paste0(dirOut, "/simuNBMM", direct, "_Ng", Ng, "Nc", Nc, snm, "_nebula.RData"))
#str(NBMM.par)
d <- NBMM.par
d$lmm.epsilon <- NULL
d$nebula.model <- NULL
d$nebula.method <- NULL
print(identical(b, d))
#[1] TRUE

rm(a, b, d)

print(identical(fit, fitss))
#[1] TRUE

runtime <- cbind(runtime, c(lmm = rtlmmSS, lmmfit = rtlmm, lmer = rtlme, nebula = rtnebula))

				
#####
##boundary (singular) 
#i <- (convlme == 0) & (!singlme)
i <- (!singlme)
print(sum(i))
vd <- c(vd, fit$theta[, i] - slme[, i])
fd <- c(fd, fit$coef[, i] - flme[, i])
td <- c(td, fit$t[, i] - tlme[, i])
pd <- c(pd, fit$p[, i] - plme[, i])


i <- (!singlme)
p1 <- c(fit$p[, i])
p2 <- c(plme[, i])
i <- (p1 < 1e-12)&(p2 < 1e-12)
#sum(i) #including all covariates, such as log(libsize) and cell-tpes
#plot(p1[i], p2[i])
#plot(-log10(p1[i]), -log10(p2[i]))
logp12d <- c(logp12d, log10(p1[!i]) - log10(p2[!i]))


range(td)
range(vd)
range(fd)
range(pd)
range(logp12d)


##t-values and p-values
range(fit$p - 2*pt(-abs(fit$t), df = fit$df))
#[1] 0 0

##t-values
tv <- t(fit$t)
dim(tv)
tv <- tv[, grep("trtB", colnames(tv))]
dim(tv)
head(tv)

tlmm <- cbind(tlmm, c(tv))

##p-values
pv <- fit$p
if (direct == "up") pv <- pt(fit$t, df = fit$df, lower.tail = F)
if (direct == "dn") pv <- pt(fit$t, df = fit$df, lower.tail = T)
pv <- t(pv)

dim(pv)
pv <- pv[, grep("trtB", colnames(pv))]
dim(pv)
head(pv)

plmm <- cbind(plmm, c(pv))

#hist(pv)
#qqpvalue(pv)

##p0: p-values for non-DE genes
##p1: p-values for DE genes
DEgenes <- NBMM.par$DEgenes
head(DEgenes)

p0 <- NULL
p1 <- NULL
for (cls in unique(DEgenes$cluster)){
	DEg <- DEgenes$gene[DEgenes$cluster == cls]
	indexDE <- (rownames(pv) %in% DEg)
	j <- grep(paste0("cls", cls, ":trt"), colnames(pv))
	p0 <- c(p0, pv[!indexDE, j])
	p1 <- c(p1, pv[indexDE, j])	
}

length(p0)
length(p1)
stopifnot(length(p1) + length(p0) == nrow(pv)*ncol(pv))

p0lmm <- p0
p1lmm <- p1

##
#qqpvalue(p0)

##type-I-error rate
##false positive rate (FPR)
cutoff <- 0.05 #0.04937081
FPR <- sum(p0 <= cutoff)/length(p0)
FPR
pcutoff <- c(0.05, 0.01)
FPR <- sapply(pcutoff, function(cutoff) sum(p0 < cutoff)/length(p0))
names(FPR) <- paste0("p", pcutoff)
FPR

##true positive rate(TPR or power)
TPR <- sum(p1 <= cutoff)/length(p1)
TPR

##false discovery rate (FDR)
FDR <- sum(p0 <= cutoff)/(sum(p0 <= cutoff) + sum(p1 <= cutoff))
FDR


#####
##nebula test
st <- negbn$summary
rownames(st) <- st$gene
dim(st)
#head(st)

st <- st[, grep("trtB", colnames(st))]
dim(st)
#head(st)


##t-values
tv <- st[, grep("logFC", colnames(st))]/st[, grep("se", colnames(st))]
tv <- as.matrix(tv)

tnbn <- cbind(tnbn, c(tv))


##p-values
pv <- st[, grep("p_", colnames(st))]
range(pv - 2*pnorm(-abs(tv)))
#[1] -3.053113e-15  3.275158e-15
if (direct == "up") pv <- pnorm(tv, lower.tail = F)
if (direct == "dn") pv <- pnorm(tv, lower.tail = T)

dim(pv)
head(pv)
pnbn <- cbind(pnbn, c(as.matrix(pv)))

#range(pnbn - 2*pnorm(-abs(tnbn)))

##p0: p-values for non-DE genes
##p1: p-values for DE genes
head(DEgenes)
p0 <- NULL
p1 <- NULL
for (cls in unique(DEgenes$cluster)){
	DEg <- DEgenes$gene[DEgenes$cluster == cls]
	indexDE <- (rownames(pv) %in% DEg)
	j <- grep(paste0("cls", cls, ":trt"), colnames(pv))
	p0 <- c(p0, pv[!indexDE, j])
	p1 <- c(p1, pv[indexDE, j])	
}

length(p1)
length(p1) + length(p0) == nrow(pv)*ncol(pv)
stopifnot(length(p1) + length(p0) == nrow(pv)*ncol(pv))

p0nbn <- p0
p1nbn <- p1

##
#qqpvalue(p0)

##type-I-error rate
##false positive rate (FPR)
cutoff <- 0.05
FPR <- sum(p0 <= cutoff)/length(p0)
FPR
pcutoff <- c(0.05, 0.01)
FPR <- sapply(pcutoff, function(cutoff) sum(p0 < cutoff)/length(p0))
names(FPR) <- paste0("p", pcutoff)
FPR

##true positive rate(TPR or power)
TPR <- sum(p1 <= cutoff)/length(p1)
TPR

##false discovery rate (FDR)
FDR <- sum(p0 <= cutoff)/(sum(p0 <= cutoff) + sum(p1 <= cutoff))
FDR

##
rm(fit, fitss, negbn)

##QQplots
if (QQp){

options(scipen=999)

qq1 <- qqpvalue(p0lmm, plot.it = F)
qq2 <- qqpvalue(p0nbn, plot.CI.only = T, shadow = T, add.grid = F)
colset <- c("red", "blue")
lines(qq2$x, qq2$y, col = colset[2])
lines(qq1$x, qq1$y, col = colset[1])
#points(qq1$x, qq1$y, col = colset[1], cex = 0.6, pch = 16)
#points(qq2$x, qq2$y, col = colset[2], cex = 0.6, pch = 16)
legend("topleft", c("lmm", "nebula"), title = paste0("\nn=",Nc), 
	pch = 16, lty = 1, col = colset, bty = "n")

options(scipen=0)
}##QQP
		
#####
##AUCs

lowcut <- log10(quantile(c(p1lmm, p1nbn), prob = 0.1))
lowcut
#cutoff <- 10^(seq(lowcut, 0, by = 0.01))
cutoff <- 10^(seq(lowcut, 0, by = 0.0025))
cutoff <- cutoff[cutoff < max(p0lmm, p1lmm, p0nbn, p1nbn)]
cutoff <- cutoff[cutoff > max(min(p0lmm), min(p0nbn))]
range(cutoff)
length(cutoff)

##lmm
auc <- AUC(p0 = p0lmm, p1 = p1lmm, cutoff = cutoff, plot.it = F)
auc$FDR <- auc$FPR/(auc$FPR + auc$TPR)

#auc$cutoff[which.min(abs(auc$cutoff - 0.05))]
#auc$FPR[which.min(abs(auc$cutoff - 0.05))]
#max(auc$FPR[auc$cutoff < 0.05])
#auc$cutoff[which.min(abs(auc$cutoff - 0.01))]
#auc$FPR[which.min(abs(auc$cutoff - 0.01))]
FPRlmm <- cbind(FPRlmm, 
	c(p0.05 = max(auc$FPR[auc$cutoff < 0.05]), p0.01 = max(auc$FPR[auc$cutoff < 0.01])))

#auc$FDR[which.min(abs(auc$FDR - 0.05))]
#auc$TPR[which.min(abs(auc$FDR - 0.05))]
#max(auc$TPR[auc$FDR < 0.05])
#auc$FDR[which.min(abs(auc$FDR - 0.1))]
#auc$TPR[which.min(abs(auc$FDR - 0.1))]
TPRlmm <- cbind(TPRlmm, 
	c(FDR0.05 = max(auc$TPR[auc$FDR < 0.05]), 
	FDR0.1 = max(auc$TPR[auc$FDR < 0.1]), 
	AUC = auc$AUC))

almm <- c(almm, list(auc))
rm(auc)

##nebula
auc <- AUC(p0 = p0nbn, p1 = p1nbn, cutoff = cutoff, plot.it = F)
auc$FDR <- auc$FPR/(auc$FPR + auc$TPR)

#auc$cutoff[which.min(abs(auc$cutoff - 0.05))]
#auc$FPR[which.min(abs(auc$cutoff - 0.05))]
#max(auc$FPR[auc$cutoff < 0.05])
#auc$cutoff[which.min(abs(auc$cutoff - 0.01))]
#auc$FPR[which.min(abs(auc$cutoff - 0.01))]
FPRnbn <- cbind(FPRnbn, 
	c(p0.05 = max(auc$FPR[auc$cutoff < 0.05]), p0.01 = max(auc$FPR[auc$cutoff < 0.01])))

#auc$FDR[which.min(abs(auc$FDR - 0.05))]
#auc$TPR[which.min(abs(auc$FDR - 0.05))]
#max(auc$TPR[auc$FDR < 0.05])
#auc$FDR[which.min(abs(auc$FDR - 0.1))]
#auc$TPR[which.min(abs(auc$FDR - 0.1))]
TPRnbn <- cbind(TPRnbn, 
	c(FDR0.05 = max(auc$TPR[auc$FDR < 0.05]), 
	FDR0.1 = max(auc$TPR[auc$FDR < 0.1]),
	AUC = auc$AUC))

anbn <- c(anbn, list(auc))
rm(auc)
}##Nc

figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_qqplot_null")
save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)
save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)


dev.off()

#####
##names
runtime
options(scipen=999)
colnames(runtime) <- paste0("n=", NcList)
runtime
#          n=20000    n=40000    n=60000     n=80000    n=100000    n=120000
#lmm      20.58146   18.32584   17.44987    17.16134    17.78659    17.40549
#lmmfit   27.57673   30.31667   39.99761    40.26968    52.55419    62.38063
#lmer   1508.92989 2982.24092 4493.96598  5640.97453  7346.69779  8947.09545
#nebula 2382.54303 5186.87765 8457.04869 12358.65186 16492.39994 19913.46079


##FPR
colnames(FPRlmm) <- colnames(runtime)
colnames(FPRnbn) <- colnames(runtime)
colnames(TPRlmm) <- colnames(runtime)
colnames(TPRnbn) <- colnames(runtime)
FPRlmm
FPRnbn
TPRlmm
TPRnbn


##AUC
str(almm)
names(almm) <- colnames(runtime)
names(anbn) <- colnames(runtime)
str(anbn)

options(scipen=0)


#####
##running time
library(xtable)

tunits <- "minutes"
rtime <- runtime/60
rtime
capt <- paste0("lmm, lmmSS, lmer and nebula computation time in ", 
		tunits, " in different sample sizes $n$.")
xtable(rtime, caption = capt, label = "tab:nbinomtime")

##lmmfit and lmmSS
tunits <- "minutes"
rtime <- runtime/60
rtime
colnames(rtime) <- gsub("n=", "", colnames(rtime))

plot(rtime[1,], ylim = range(rtime), axes = F, type = "n", 
	xlab = "Number of cells", ylab = paste0("Computation time in ", tunits))
#grid()
box()
axis(1, 1:ncol(rtime), labels = colnames(rtime))
axis(2)
colset <- c("red", "red", "green", "blue") 
pchset <- c(4, 1, 17, 15)
cexset <- c(1, 1, 1, 1)
for (i in 1:nrow(rtime)){
	if (i != 1) lines(rtime[i,], col = colset[i], lty = 3)
	points(rtime[i, ], col = colset[i], pch = pchset[i], cex = cexset[i])
}
legend("topleft", rownames(rtime), pch = pchset, lty = 3, col = colset, bty = "n")
	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_runtime")
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)

##ratio
timeratio <- rtime/matrix(rtime[rownames(rtime) == "lmm", ], nrow = nrow(rtime), ncol = ncol(rtime), byrow = T)
timeratio
timeratio <- timeratio[rownames(timeratio) != "lmm",]
timeratio
rownames(timeratio) <- paste0(rownames(timeratio), " vs lmm")

plot(timeratio[1,], ylim = range(timeratio), axes = F, type = "n", 
	xlab = "Number of cells", ylab = "Ratio of computation time")
#grid()
box()
axis(1, 1:ncol(rtime), labels = colnames(rtime))
axis(2)
for (i in 1:nrow(timeratio)){
	lines(timeratio[i, ], lty = 3, col = colset[i+1])
	points(timeratio[i, ], pch = pchset[i+1], col = colset[i+1], cex = 1.1)
	}
legend("topleft", rownames(timeratio), pch = pchset[-1], lty = 3, col = colset[-1], bty = "n")
	
	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_timeratio")
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)



##lmmfit
rtime <- runtime[-1, ]/60
rtime
colnames(rtime) <- gsub("n=", "", colnames(rtime))

plot(rtime[1,], ylim = range(rtime), axes = F, type = "n", 
	xlab = "Number of cells", ylab = paste0("Computation time in ", tunits))
#grid()
box()
axis(1, 1:ncol(rtime), labels = colnames(rtime))
axis(2)
colset <- c("red", "green", "blue") 
pchset <- c(16, 17, 15)
cexset <- c(1, 1, 1)
for (i in 1:nrow(rtime)){
	lines(rtime[i,], col = colset[i], lty = 3)
	points(rtime[i, ], col = colset[i], pch = pchset[i], cex = cexset[i])
}
legend("topleft", rownames(rtime), pch = pchset, lty = 3, col = colset, bty = "n")

	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_runtime_lmmfit")
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)


#####
##FPR
nm <- "p0.01"; p0 <- 0.01
nm <- "p0.05"; p0 <- 0.05

cnm <- gsub("n=", "", colnames(FPRlmm))

a <- rbind(lmm = FPRlmm[nm,], nebula = FPRnbn[nm,])
#nm <- "Type I error rate"
ylim <- range(a) #c(0.98*min(a), 1.02*max(a))
#par(mar = c(0, 4.1, 5.1, 2.1))
plot(a[1,], ylim = ylim, axes = F, type = "n", xlab = "Number of cells", ylab = "FPR")
box()
#mtext(nm, 3, line = 0.5, cex = 0.9)
axis(1, 1:ncol(a), labels = cnm)
axis(2)
abline(h = p0, col = "gray")
colset <- c("red", "blue") 
pchset <- c(16, 17)
for (i in 1:nrow(a)){
	lines(a[i,], col = colset[i], lty = 3)
	points(a[i, ], col = colset[i], pch = pchset[i])
}
legend("bottomright", rownames(a), pch = pchset, lty = 3, col = colset, bty = "n")

	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_FPR", nm)
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)

dev.off()

#####
##TPR
nm <- "FDR0.05"
#nm <- "FDR0.1"
#nm <- "AUC"

cnm <- gsub("n=", "", colnames(TPRlmm))

a <- rbind(lmm = TPRlmm[nm,], nebula = TPRnbn[nm,])
#nm <- "Power at FDR = 0.1"
ylim <- range(a) #c(0.98*min(a), 1.02*max(a))
#par(mar = c(0, 4.1, 5.1, 2.1))
plot(a[1,], ylim = ylim, axes = F, type = "n", xlab = "Number of cells", ylab = "FPR")
box()
#mtext(nm, 3, line = 0.5, cex = 0.9)
axis(1, 1:ncol(a), labels = cnm)
axis(2)
abline(h = 0.05, col = "gray")
colset <- c("red", "blue") 
pchset <- c(16, 17)
for (i in 1:nrow(a)){
	lines(a[i,], col = colset[i], lty = 3)
	points(a[i, ], col = colset[i], pch = pchset[i])
}
legend("bottomright", rownames(a), pch = pchset, lty = 3, col = colset, bty = "n")

	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_TPR_", nm)
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)

dev.off()

#####
##ROC
str(almm)
str(anbn)

colset <- c("red", "blue")
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1.1))
for (i in 1:length(almm)){
	auclmm <- almm[[i]]
	aucnbn <- anbn[[i]]
ylim <- range(auclmm$TPR, aucnbn$TPR)
plot(auclmm$FPR, auclmm$TPR, ylim = ylim, type = "n", cex.axis = 0.9, cex.lab = 0.9,
	xlab = "False positive rate", ylab = "True positive rate")
lines(auclmm$FPR, auclmm$TPR, col = colset[1], lty = 1)
lines(aucnbn$FPR, aucnbn$TPR, col = colset[2], lty = 2)
txt <- paste0(c("lmm AUC = ", "nebula AUC = "), round(c(auclmm$AUC, aucnbn$AUC), 2))
legend("bottomright", txt, col = colset, lty = 1:2, bty = "n", title = names(almm)[i], cex = 0.9)
}
	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_AUC")
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)

dev.off()

#####
##t-values
plot(tlmm, tnbn, xlab = "lmm t-values", ylab = "nebula t-values", cex = 0.65)
abline(0, 1)
	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_t")
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)

##p-values
range(pnbn - 2*pnorm(-abs(tnbn)))
#[1] -3.608225e-15  3.247402e-15


##
i <- (plmm == 0) | (pnbn == 0)
plmm[i]
pnbn[i]

plot(-log10(plmm[!i]), -log10(pnbn[!i]), 
	xlab = "lmm: -log10(p-values)", ylab = "nebula: -log10(p-values)", cex = 0.65)
abline(0, 1)
	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_p")
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)


##########
##boxplot of differences of variance componenets, coefficients, and t-values 
##between lmm and lmer with default setting

options(scipen=0)

d <- list('Variance components' = c(vd), 'Coefficients' = c(fd), "t-values" = c(td), 
	"p-values" = pd)
boxplot(d, ylab = "Differences between lmm and lmer", cex.axis = 0.9, cex.lab = 0.9)
	figFile <- paste0(dirFig, "/simuNBMM", direct, snm, "_lmmfitdiff")
	save.figure(figFile, "png", widthscale = scaFig, heightscale = scaFig, res = resFig)
	save.figure(figFile, "pdf", widthscale = scaFig, heightscale = scaFig, res = resFig)

range(vd)
#[1] -9.020614e-07  7.116696e-06
range(fd)
#[1] -5.198020e-08  1.400356e-07
range(td)
#[1] -0.0003256032  0.0020077457
range(pd)
#[1] -4.447217e-06  4.578032e-06



##################################################
##################################################



	