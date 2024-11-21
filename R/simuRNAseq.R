##########	
##Simulating scRNA-seq counts data by negative binomial distribution (NB):
##(1) NB means and dispersions for genes estimated by the method-of-moments estimate (MME) 
##    using a reference dataset as background or control.
##(2) Multiple subjects (individuals or samples), multiple cell-types (clusters), 
##    and two treatments (conditions) for differential expressed (DE) genes.
##(3) The meta data consisting of samples, clusters of cell types, and treatments 
##    will be randomly generated if it is not provided.
##(4) The DE gene efects: beta ~ +/- uniform(minbeta, maxbeta), up and down-regulated DE genes

##Arguments
##- counts: genes-by-cells matrix of reference counts 
##  used for estimating the mean and dispersion of NB distribution.
##- nGenes: number of genes
##- nCells: number of cells
##- metadata: meta data consisting of 4 columns: 
##  sam (sample labels), 
##  cls (cluster lables of cell types), 
##  trt (treatments or conditions) and
##  libsize (library size or total counts per cell),
##  which is randomly generated if it is not provided.
##- ntrt: number of treatments (only one condition is treated.)
##- ncls: number of clusters (cell-types)
##- nsam: number of samples (subjects or individuals)
##- nDEgenes: total number of DE genes
##- nDEgenesType: number of DE genes specific to a cell type, named by cell cluster labels.
##- pDEgenesType, proportion of DE genes in a cell-type. Default NULL means equal proportion.
##- adjust.library.size, if TRUE, adjust library sizes using the reference counts.
##- minbeta and maxbeta: lower and upper bounds of DE gene effects, 
##  beta ~ Uniform(minbeta, maxbeta) for the DE genes.
##- one.direction: 
##  TRUE (only up- or down-regulated DE genes) or 
##  FALSE (including both up- and down-regulated DE genes).
##- trt: one of metadata$trt (treatments), specifying which condition is treated.
##- var.randomeffects: variance of random effects

##Value
##- mean.dispersion: data frame of the means and dispersions
##- metadata: meta data consisting of 4 columns: 
##  sam: sample labels
##  cls: cluster lables of cell types
##  trt: two treatments (conditions)
##  libsize: library sizes
##- counts: genes-by-cells matrix of the simulated counts
##- DEgenes: data frame of DE genes consisting of 3 columns:
##  gene, beta (effect), and cluster to which the gene is specific.
##- treatment: condition treated.

##Note:
##A random seed is recommended to be specified by set.seed before simulating data.
##Filtering the genes with Var < Mean or Mean = 0.

##########
simuRNAseq <- function(counts, nGenes = nrow(counts), nCells = ncol(counts), metadata = NULL, nsam = 25, ncls = 10, ntrt = 2, nDEgenes = ncls, nDEgenesType, pDEgenesType = NULL, adjust.library.size = TRUE, minbeta = 0.25, maxbeta = 1, one.direction = FALSE, trt = NULL, var.randomeffects = 0.01)
{
stopifnot(minbeta < maxbeta)

##MME dispersions
mu <- rowMeans(counts)
if (any(mu <= 0)){
	index <- (mu > 0)
	message(paste0("Warnings: removing ", sum(!index), " row(s) with zero mean."))
	if (nGenes == nrow(counts)) nGenes <- sum(index)
	counts <- counts[index, ]
	mu <- mu[index]
	}
	
n <- ncol(counts)
v <- rowSums(counts^2)/(n-1) - n*mu^2/(n-1)
size <- mu^2/(v - mu)

if (any(v <= mu)){
	index <- (v > mu)
	message(paste0("Message: removing ", sum(!index), " row(s) with var <= mu."))
	if (nGenes == nrow(counts)) nGenes <- sum(index)
	mu <- mu[index]
	size <- size[index]
	counts <- counts[index, ]
	}

##mean and dispersion for the genes in the reference data
meandisp <- data.frame(mu = mu, dispersion = size)
rownames(meandisp) <- rownames(counts)

##mean and dispersion for the genes to be simulated
if (nrow(meandisp) == nGenes){
	musize <- meandisp
} else musize <- meandisp[sample.int(nrow(meandisp), nGenes, replace = (nrow(meandisp) < nGenes)), ]

mu <- musize$mu      
size <- musize$dispersion 

##subjects (samples)
##clusters (cell-types)
##treatments (conditions)
if (is.null(metadata)){
	metadata <- data.frame(
	sam = as.factor(sample.int(nsam, nCells, replace = TRUE)),
	cls = as.factor(sample.int(ncls, nCells, replace = TRUE)),
	trt = LETTERS[sample.int(ntrt, nCells, replace = TRUE)])
	if (adjust.library.size){
		if (ncol(counts) != nCells){
		metadata$libsize <- sample(colSums(counts), nCells, replace = (ncol(counts) < nCells))
		} else metadata$libsize <- colSums(counts)
	} else metadata$libsize <- 1
	rownames(metadata) <- paste0("cell", 1:nCells)
} else {
	stopifnot(all(c("sam", "cls", "trt", "libsize") %in% colnames(metadata)))
	if (!is.null(trt)) stopifnot(trt %in% metadata$trt)
	if (nrow(metadata) != nCells) 
	metadata <- metadata[sample.int(nrow(metadata), nCells, replace = (nrow(metadata) < nCells)), ]	
	ntrt <- length(unique(metadata$trt))
	}

##Number of DE genes in each cell-type (cluster)
clusters <- sort(unique(metadata$cls))
if (missing(nDEgenesType)){
	ncls <- length(clusters)
	if (is.null(pDEgenesType)) pDEgenesType <- rep(1/ncls, ncls)
	stopifnot(length(pDEgenesType) == ncls)
	nDEgenesType <- c(rmultinom(1, nDEgenes, prob = pDEgenesType))
	names(nDEgenesType) <- clusters
} else {
	stopifnot(length(nDEgenesType) == length(clusters))
	if (is.null(names(nDEgenesType))){
	message("Messages: 'nDEgenesType' is named by cell-types (clusters)")
	names(nDEgenesType) <- clusters
	}
}
nDEgenesType <- nDEgenesType[nDEgenesType > 0]	
clusters <- NULL

##nDEgenes: numbers of DE genes
##Ng0: number of non-DE genes
##Ng: total genes
if (sum(nDEgenesType) != nDEgenes){
	nDEgenes <- sum(nDEgenesType)
	message(paste0("Warning: number of DE genes is changed as ", nDEgenes, "."))
	}
Ng <- nrow(musize)
stopifnot(Ng >= nDEgenes)
Ng0 <- Ng - nDEgenes

##Counts
counts <- matrix(NA, nrow = Ng, ncol = nCells)
rownames(counts) <- rownames(musize)
colnames(counts) <- rownames(metadata) #paste0("cell", 1:nCells)

##(1) Constant effects (non-DE or no change): diff(betas) = 0 or 
##    betas = log(mu) = log(mu0) + centered-library-size
##(2) Add random samples effects: log(mu_s) = Xs*beta + Zs*Bs = log(mu) + Zs*Bs
##(3) Treatment (fixed) effects in one subpopulation (one cluster) with trt = "B".

NgList <- c("0"=Ng0, nDEgenesType)	

##DE genes
DEgenes <- data.frame(gene = rownames(musize)[(Ng0+1):Ng], 
		beta = runif(nDEgenes, minbeta, maxbeta), 
		cluster = NA)
if (!one.direction) DEgenes$beta <- (2*rbinom(nDEgenes, 1, prob = 0.5)-1)*DEgenes$beta
rownames(DEgenes) <- (Ng0+1):Ng
for (i in 1:length(nDEgenesType)){
	j <- (1+sum(NgList[1:i])):sum(NgList[1:(1+i)]) - Ng0
	DEgenes$cluster[j] <- names(nDEgenesType)[i]
	}

##center log-library-size
libsize <- log(metadata$libsize)
libsize <- libsize - mean(libsize)

##trt: the treatment considered as the last condition if not specified.
if (is.null(trt)) {
	trt <- sort(unique(metadata$trt))
	trt <- trt[length(trt)]
	}
	
##two conditions (control and treatment)
if (ntrt == 1) message("Message: only one condition, no treatment effects.")
if (ntrt > 2) message("Warning: more than 2 conditions. Only one condition is considered as treatment.")

##counts generated by sampling
for (sam in unique(metadata$sam)){
	index <- (metadata$sam == sam)
	Ns <- sum(index)	
	##contant effects + random effects
	lfc <- log(mu) + rnorm(length(mu), sd = sqrt(var.randomeffects))
	lfcsam <- lfc + matrix(libsize[index], nrow = Ng, ncol = Ns, byrow = TRUE)
	counts[, index] <- matrix(rnbinom(Ng*Ns, size = size, mu = exp(lfcsam)), nrow = Ng, ncol = Ns)
	
	##treatment effects
	for (i in 1:length(nDEgenesType)){
		indexB <- index & (metadata$trt == trt) & (metadata$cls == names(nDEgenesType)[i])
		Nsi <- sum(indexB)	
		if (Nsi > 0){
		Ngi <- nDEgenesType[i]		
		j <- (1+sum(NgList[1:i])):sum(NgList[1:(1+i)])
		lfci <- lfc[j] + DEgenes$beta[j - Ng0]
		lfci <- lfci + matrix(libsize[indexB], nrow = Ngi, ncol = Nsi, byrow = TRUE)
		counts[j, indexB] <- matrix(rnbinom(Ngi*Nsi, size = size[j], mu = exp(lfci)), nrow = Ngi, ncol = Nsi)
		}
	}
	}

list(mean.dispersion = meandisp, metadata = metadata, counts = counts, DEgenes = DEgenes, treatment = trt)
}
