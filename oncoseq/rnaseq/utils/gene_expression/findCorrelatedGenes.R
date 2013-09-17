# packages
library(qvalue)

exprmat <- function(mdata, cdata, pdata) {
	# make exprmat object
	m <- list()
	class(m) <- "exprmat"	
	m$pdata <- pdata
	m$mdata <- mdata
	m$cdata <- cdata
	return(m)
}

read.exprmat <- function(phenoFile, countFile, key="library_id") {
	# read pheno and count data
	pdata <- read.table(phenoFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names=key, check.names=FALSE)
	firstid <- rownames(pdata)[1]
	cdata <- read.table(countFile, header=TRUE, sep="\t", row.names="tracking_id", check.names=FALSE)
	firstcol <- which(colnames(cdata) == firstid)[1]
	# make exprmat object
	m <- list()
	class(m) <- "exprmat"	
	m$pdata <- pdata
	m$mdata <- cdata[,1:firstcol-1]
	m$cdata <- cdata[,firstcol:ncol(cdata)]
	return(m)
}

getCorrelatedGenes <- function(m, gene_id, dolog=FALSE) {
	# log transform
	if (dolog) {
		exprmat <- log2(m$cdata + 1.0)
	} else {
		exprmat <- m$cdata
	}
	# observed correlation
	cors <- cor(matrix(exprmat[geneid,],ncol=1), y=t(exprmat))
	r <- t(cors)
	# return
	res <- cbind(m$mdata, r)
	res
}

args <- commandArgs(trailingOnly = TRUE)
phenoFile <- args[1]
exprFile <- args[2]
dolog <- args[3]
outputFile <- args[4]

dolog <- ifelse(dolog == "y", TRUE, FALSE)
m <- read.exprmat(phenoFile, exprFile, key="library_id")

res <- getCorrelatedGenes(m, gene_id, dolog)
write.table(res, file=outputFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

