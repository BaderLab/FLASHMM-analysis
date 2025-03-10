
qqpvalue <- function(p, CI = 0.95, ylim, xlab = "Expected: -log10(uniform)", ylab = "Observed: -log10(p-values)", shadow = TRUE, shadow.col = "gray", diagonal = TRUE, diagonal.col = "gray90", zero.rm = TRUE, add.grid = FALSE, log.qq = TRUE, plot.it = TRUE, plot.CI.only = FALSE, trend = FALSE, ...)
{
##p: a vector of pvalues
##CI: if numerical, add confidence intervals
##trend: 
##- If TRUE, add a straight line to plot. 
##- The line is a fit of the quantiles within [0, 0.95].


stopifnot(CI > 0.5)

##remove NA
p <- p[!is.na(p)]

stopifnot(p >= 0)
orderp <- order(p)
#p <- sort(p)
p <- p[orderp]

if (zero.rm) {
	orderp <- orderp[p != 0]
	p <- p[p != 0]
} else if (any(p == 0)) {
	log.qq = FALSE
	message("Warning: cannot plot log-QQ because some p = 0.")
	}

##null distribution
n <- length(p)
p0 <- (1:n - 0.5)/n

##confidence intervals
##- the jth order statistic of uniform(0,1) has a beta(j,n-j+1) distribution
	j <- 1:n
	CI <- (1 - CI)/2
	ciu <- qbeta(1 - CI, j, n - j + 1)
	cil <- qbeta(CI, j, n - j + 1)

##output CIs
pci <- data.frame(p = p, lb = cil, ub = ciu)

##log-QQ
if (log.qq){
	p0 <- -log10(p0)
	p <- -log10(p)
	ciu <- -log10(ciu)
	cil <- -log10(cil)
	}


##########
if (plot.it){

if (missing(ylim)) ylim <- range(c(p, ciu, cil))
else ylim <- range(c(ylim, p, ciu, cil))

##plot
#qqplot(p0, p, ylim = ylim, xlab = xlab, ylab = ylab, type = "n", ...)
dots <- list(...)
if (any(names(dots) == "main")) main = dots$main else main = NULL
if (any(names(dots) == "cex.main")) cex.main = dots$cex.main else cex.main = 1 
if (any(names(dots) == "cex.lab")) cex.lab = dots$cex.lab else cex.lab = 1 
if (any(names(dots) == "cex.axis")) cex.axis = dots$cex.axis else cex.axis = 1 

plot(p0, p, ylim = ylim, xlab = xlab, ylab = ylab, cex.lab = cex.lab, cex.axis = cex.axis, 
type = "n", main = main, cex.main = cex.main)
if (add.grid) {grid(); box()}

##add shadows
	if (shadow){
	x <- seq(min(p0), max(p0), length = min(1e3, length(p0)*20))
	y1 <- approx(p0, cil, n = length(x))$y
	y2 <- approx(p0, ciu, n = length(x))$y
	polygon(c(x, rev(x)), c(y2,rev(y1)), border = shadow.col, col = shadow.col)
	} 
	
	if (diagonal) abline(0, 1, col = diagonal.col)

	if (!plot.CI.only){
	points(p0, p, ...)
	if (trend) {
		i <- (p0 <= min(1, quantile(p0, probs = 0.25))) 
		fit <- lm(p[i] ~ p0[i])
		abline(coef = coefficients(fit), col = "gray", lwd = 0.6, lty = "dashed")
		}#
	}#!plot.CI.only
}#plot.it

##########
invisible(list(x = p0, y = p, p.ci = pci, order = orderp))
}

##example
##p <- pnorm(rnorm(100))
##qqpvalue(p, CI = NULL, col = "blue", diagonal.col = "red", cex = 0.5, cex.lab = 0.5, cex.main = 0.5)
##qqpvalue(p, CI = F, col = "blue", cex = 0.5, cex.lab = 0.85, cex.main = 0.85)
##qqpvalue(p, CI = T, col = "blue", cex = 0.5, cex.lab = 0.85, cex.main = 0.85)
##qqpvalue(p, col = "blue", cex = 0.5, cex.lab = 0.85, main = "QQ-plot", cex.main = 0.85)
##qqpvalue(p, col = "blue", cex = 0.5, cex.lab = 0.85, main = "QQ-plot", cex.main = 0.85, log.qq = F)



