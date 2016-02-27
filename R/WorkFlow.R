# Copyright 2015, Institute for Systems Biology.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# Author: William Poole
# Ported to R: David L Gibbs
# Email: dgibbs@systemsbiology.org / william.poole@systemsbiology.org / tknijnen@systemsbiology.org
# Created: June 2015
# Updated: February 2016

source("EmpiricalBrownsMethod/R/ebm.R")
options(digits=16)

# Test with artificical data

print("****************************************************************")
print("RANDOM DATA")
randData <- read.table("../Data/RandomData.tsv", sep="\t", header=F, stringsAsFactors=F)
a <- as.numeric(randData[1,-1])
rd <- randData[-1,-1]
pvals <- sapply(1:10, function(i) cor.test(a, as.numeric(rd[i,]))$p.value)
print(empiricalBrownsMethod(data_matrix=rd, p_values=pvals, extra_info=T))
print("****************************************************************")

#Should give:
#$P_Brown
#[1] 0.7228817373295435
#
#$P_Fisher
#[1] 0.8613842570343421
#
#$Scale_Factor_C
#[1] 2.45800963585645
#
#$DF_Brown
#[1] 8.136664603851868

print(kostsMethod(data_matrix=as.matrix(rd), p_values=pvals, extra_info=T))
print("****************************************************************")

#Should give:
# $P_test
# [1] 0.701752883272515
# 
# $P_Fisher
# [1] 0.8613842570343421
# 
# $Scale_Factor_C
# [1] 2.814405567447344
# 
# $DF
# [1] 7.106296345959808


# Test with cancer data #

pathways <- read.table("../Data/pathways.tsv", sep="\t", header=T, stringsAsFactors=F)
allPvals <- read.table("../Data/CDH4_Pvalues.tsv", sep="\t", stringsAsFactors=F, header=T )
dat <- read.table("../Data/ReducedFeatureMatrix.tsv", sep="\t", stringsAsFactors=F, header=F)
allPvals <- (unique(allPvals))

print("Glypican 3 Network")
glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"]
glypPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% glypGenes]
glypDat <- dat[dat$V1 %in% glypGenes, 2:ncol(dat)]
print(empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=T))
print("****************************************************************")

#Should give:
#$P_Brown
#[1] 4.821679406409943e-07
#
#$P_Fisher
#[1] 1.438732140605804e-08
#
#$Scale_Factor_C
#[1] 1.297692749787417
#
#$DF_Brown
#[1] 10.78837806737645

print(kostsMethod(data_matrix=as.matrix(glypDat), p_values=glypPvals, extra_info=T))
print("****************************************************************")

#Should give:
# $P_test
# [1] 7.570776008807138e-07
# 
# $P_Fisher
# [1] 1.438732140605804e-08
# 
# $Scale_Factor_C
# [1] 1.349048766471012
# 
# $DF
# [1] 10.3776826664485

print("SUMO Reg")
sumoGenes <- pathways$gene[pathways$pathway == "SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION"]
sumoPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% sumoGenes]
sumoDat <- dat[dat$V1 %in% sumoGenes, 2:ncol(dat)]
print(empiricalBrownsMethod(data_matrix=sumoDat, p_values=sumoPvals, extra_info=T))
print("****************************************************************")

#Should give:
#$P_Brown
#[1] 1.698056395040471e-41
#
#$P_Fisher
#[1] 6.443838824431309e-45
#
#$Scale_Factor_C
#[1] 1.087731057365708
#
#$DF_Brown
#[1] 18.38689799704392

print("FOXA1 Network")
foxGenes <- pathways$gene[pathways$pathway == "FOXA1 TRANSCRIPTION FACTOR NETWORK"]
foxPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% foxGenes]
foxDat <- dat[dat$V1 %in% foxGenes, 2:ncol(dat)]

print(empiricalBrownsMethod(data_matrix=foxDat, p_values=foxPvals, extra_info=T))
print("****************************************************************")

#Should give:
#$P_Brown
#[1] 7.777896979417795e-53
#
#$P_Fisher
#[1] 4.043406925735124e-139
#
#$Scale_Factor_C
#[1] 2.719366560758496
#
#$DF_Brown
#[1] 21.32849643625184



