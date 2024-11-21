
##########
##lmm

load(paste0(dirOut, "/kidney-counts-lmm.RData"))

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
##p-values for testing variance components

##t-values
tvlmm <- t(fit$t)

dim(tvlmm)
head(tvlmm)
sum(is.na(tvlmm))
#[1] 28
sum(apply(is.na(tvlmm), 1, any))
#[1] 19

##p-values
pvlmm <- t(fit$p)

dim(pvlmm)
head(pvlmm)

#####
##fixed effects (logFC), t-values and p-values 
##for Male vs Female specific to a cell-type
index <- grep(":Male", colnames(pvlmm))

fe <- felmm[, index]
tv <- tvlmm[, index]
pv <- pvlmm[, index]


##########
##nebula
#load(paste0(dirOut, "/kidney-counts-nebula.NBLMM.RData"))
##Using default setting
load(paste0(dirOut, "/kidney-counts-nebula.default.RData"))

#####
#str(negbn)
table(negbn$convergence)
 
st <- negbn$summary
rownames(st) <- st$gene
dim(st)
head(st)

any(is.na(st))
sum(is.na(st))

##fixed effects, se, and p-values
iFC <- grep("logFC_", colnames(st))
ise <- grep("se_", colnames(st))
ipv <- grep("p_", colnames(st))

##fixed effects
feneb <- as.matrix(st[, iFC])
##se
se <- as.matrix(st[, ise])
##t-values
tvneb <- feneb/se
##p-values
pvneb <- as.matrix(st[, ipv])

range(pvneb - 2*pnorm(-abs(tvneb)), na.rm = T)
#[1] -2.720046e-15  3.136380e-15

#####
##fixed effects (logFC), t-values and p-values 
##for Male vs Female specific to a cell-type
index <- grep(":Male", colnames(pvneb))

fe <- feneb[, index]
tv <- tvneb[, index]
pv <- pvneb[, index]
