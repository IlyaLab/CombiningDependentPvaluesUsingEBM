#
# TEST SCRIPT
# Empirical Browns Method of Combining P-values.
# Created on Mon Jun 22, 2015
#
# author: William Poole: wpoole@systemsbiology.org
# ported to R: David L Gibbs: dgibbs@systemsbiology.org
#

options(digits=16)

context("Main Function Testing")

# Test with random data

data("testData")
#randData <- read.table("inst/extdata/RandomData.tsv", sep="\t", header=F, stringsAsFactors=F)
#pathways <- read.table("inst/extdata/pathways.tsv", sep="\t", header=T, stringsAsFactors=F)
#allPvals <- read.table("inst/extdata/CDH4_Pvalues.tsv", sep="\t", stringsAsFactors=F, header=T )
#dat <- read.table("inst/extdata/ReducedFeatureMatrix.tsv", sep="\t", stringsAsFactors=F, header=F)
#allPvals <- (unique(allPvals))

test_that("random results are correct", {
  a <- as.numeric(randData[1,-1])
  rd <- randData[-1,-1]
  pvals <- sapply(1:10, function(i) cor.test(a, as.numeric(rd[i,]))$p.value)
  res0 <- empiricalBrownsMethod(data_matrix=rd, p_values=pvals, extra_info=T)
  expected = list(P_Brown=0.7228817373295435, P_Fisher=0.8613842570343421, Scale_Factor_C=2.45800963585645, DF_Brown=8.136664603851868)
  expect_equal(res0, expected)
})

# Test the pathways #

#Should give:
#pathway: FOXA1 TRANSCRIPTION FACTOR NETWORK
#(7.7778969794178595e-53, 4.043406925735029e-139, 2.7193665607584965, 21.328496436251836)
#pathway: SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION
#(1.6980563950404756e-41, 6.4438388244313223e-45, 1.0877310573657077, 18.386897997043924)
#pathway: GLYPICAN 3 NETWORK
#(4.8216794064099692e-07, 1.4387321406058163e-08, 1.2976927497874169, 10.788378067376447)
#(Pbrown,Pfisger,Scale_Factor C,DFbrown)


test_that("the glypican 3 network example works", {
  glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"]
  glypPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% glypGenes]
  glypDat <- dat[dat$V1 %in% glypGenes, 2:ncol(dat)]
  res0 <- empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=T)
  expected = list(P_Brown=4.821679406409943e-07, P_Fisher=1.438732140605804e-08, Scale_Factor_C=1.297692749787417, DF_Brown=10.78837806737645)
  expect_equal(res0, expected)
})

test_that("the SUMO pathway example works", {
  sumoGenes <- pathways$gene[pathways$pathway == "SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION"]
  sumoPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% sumoGenes]
  sumoDat <- dat[dat$V1 %in% sumoGenes, 2:ncol(dat)]
  res0 <- empiricalBrownsMethod(data_matrix=sumoDat, p_values=sumoPvals, extra_info=T)
  expected = list(P_Brown=1.698056395040471e-41, P_Fisher=6.443838824431309e-45, Scale_Factor_C=1.087731057365708, DF_Brown=18.38689799704392)
  expect_equal(res0, expected)
})

test_that("the FOXA1 network example works", {
  foxGenes <- pathways$gene[pathways$pathway == "FOXA1 TRANSCRIPTION FACTOR NETWORK"]
  foxPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% foxGenes]
  foxDat <- dat[dat$V1 %in% foxGenes, 2:ncol(dat)]
  res0 <- empiricalBrownsMethod(data_matrix=foxDat, p_values=foxPvals, extra_info=T)
  expected = list(P_Brown=7.777896979417795e-53, P_Fisher=4.043406925735124e-139, Scale_Factor_C=2.719366560758496, DF_Brown=21.32849643625184)
  expect_equal(res0, expected)
})
