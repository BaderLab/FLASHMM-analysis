
##########
AUC <- function(FPR, TPR, p0, p1, y, yhat, cutoff, plot.it = FALSE, type = "l", add.AUC = TRUE, ...){
##AUC: area under the ROC (Receiver Operating Characteristic) curve
##FPR: false positive rate
##TPR: true positive rate
##The area under the ROC curve:
##consisting of rectangles and triangles (half of rectangles) 
##p0: p-values of the negative controls (under H0)
##p1: p-values of the positive controls (under H1)
##y: a vector of binary variable 0 and 1
##yhat: predicted y
##cutoff default:
##cutoff = sort(unique(c(p0, p1))) for p-values
##cutoff = seq(0, 1, by = 0.01) for y

if (missing(y) | missing(yhat)){
	if (missing(p0) | missing(p1)){
	##sorted 
	o <- order(FPR, TPR)
	FPR <- FPR[o]
	TPR <- TPR[o]
	} else {
	if (missing(cutoff)) cutoff <- sort(unique(c(p0, p1)))
	FPR <- rep(NA, length(cutoff))
	TPR <- FPR
	for (i in 1:length(cutoff)){
		FPR[i] <- sum(p0 <= cutoff[i])/length(p0)
		TPR[i] <- sum(p1 <= cutoff[i])/length(p1)
	}	
	}
} else {
if (missing(cutoff)) cutoff <- seq(0, 1, by = 0.01)
cutoff <- sort(cutoff, decreasing = TRUE)
FPR <- NULL
TPR <- NULL
for (p0 in cutoff) {
	tab <- table(y, yhat >= p0)
	if (all(colnames(tab) == "FALSE")) tab <- cbind(tab, "TRUE" = 0)
	FPR <- c(FPR, tab["0", "TRUE"]/sum(tab["0", ]))
	TPR <- c(TPR, tab["1", "TRUE"]/sum(tab["1", ]))
}
}

FPR0 <- c(0, FPR, 1)
TPR0 <- c(0, TPR, 1)
dFPR <- c(diff(FPR0), 0)
dTPR <- c(diff(TPR0), 0)

auc <- sum(TPR0 * dFPR) + sum(dTPR * dFPR)/2

if (plot.it) {
	plot(FPR0, TPR0, xlab = "False positive rate", ylab = "True positive rate", type = type, ...)
	abline(0, 1, col = "gray")
	if (add.AUC) text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), paste0("AUC=", round(auc, 4)))
	}
#return(auc)
list(AUC = auc, cutoff = cutoff, FPR = FPR, TPR = TPR)
}
