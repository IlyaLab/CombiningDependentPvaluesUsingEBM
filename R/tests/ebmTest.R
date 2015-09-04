#
# TEST SCRIPT
# Empirical Browns Method of Combining P-values.
# Created on Mon Jun 22, 2015
#
# author: William Poole: wpoole@systemsbiology.org
# ported to R: David L Gibbs: dgibbs@systemsbiology.org
#


source("ebm.R")
options(digits=16)

# Test with random data

print("****************************************************************")
print("RANDOM DATA")
randData <- read.table("RandomData.tsv", sep="\t", header=F, stringsAsFactors=F)
a <- as.numeric(randData[1,-1])
rd <- randData[-1,-1]
pvals <- sapply(1:10, function(i) cor.test(a, as.numeric(rd[i,]))$p.value)
print(empiricalBrownsMethod(data_matrix=rd, p_values=pvals, extra_info=T))
print("****************************************************************")

# Test with cancer data #

pathways <- read.table("pathways.tsv", sep="\t", header=T, stringsAsFactors=F)
allPvals <- read.table("CDH4_Pvalues.tsv", sep="\t", stringsAsFactors=F, header=T )
dat <- read.table("ReducedFeatureMatrix.tsv", sep="\t", stringsAsFactors=F, header=F)
allPvals <- (unique(allPvals))

print("Glypican 3 Network")
glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"]
glypPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% glypGenes]
glypDat <- dat[dat$V1 %in% glypGenes, 2:ncol(dat)]
print(empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=T))
print("****************************************************************")

print("SUMO Reg")
sumoGenes <- pathways$gene[pathways$pathway == "SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION"]
sumoPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% sumoGenes]
sumoDat <- dat[dat$V1 %in% sumoGenes, 2:ncol(dat)]
print(empiricalBrownsMethod(data_matrix=sumoDat, p_values=sumoPvals, extra_info=T))
print("****************************************************************")

print("FOXA1 Network")
foxGenes <- pathways$gene[pathways$pathway == "FOXA1 TRANSCRIPTION FACTOR NETWORK"]
foxPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% foxGenes]
foxDat <- dat[dat$V1 %in% foxGenes, 2:ncol(dat)]
print(empiricalBrownsMethod(data_matrix=foxDat, p_values=foxPvals, extra_info=T))
print("****************************************************************")
