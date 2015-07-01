# -*- coding: utf-8 -*-
"""
Copyright 2015, Institute for Systems Biology.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Author: William Poole
Email: william.poole@systemsbiology.org / tknijnen@systemsbiology.org
Created: June 2015
"""

from numpy import *
from EmpiricalBrownsMethod import *
from scipy.stats import pearsonr
import PathwayParser as PP




# ARTIFICIAL DATASET
#RandomData.tsv contains gaussian random data.
#--Independent Var [line 1] are 25 samples from a unit normal distribution.
#--Depedent Var 1-10 [line 2-11] are each 25 samples drawn from a 10 dimensional normal distribution centered at the origin with off diagonal terms a=0.25.
#--The P values from a pearson correlation between the independent var and each dependent var are combined

raw_data = open("../Data/RandomData.tsv")
data = []
for line in raw_data:
    L = line.replace("\n", "").split("\t")
    if "Independent Var" in L[0]:
        indV = array([float(l) for l in L[1:]])
    else:
        data.append([float(l) for l in L[1:]])

raw_data.close()
data = np.array(data)
pvals = [pearsonr(indV, data[i])[1] for i in range(data.shape[0])]
transformed_data1 = TransformData(data[0, :])
print EmpiricalBrownsMethod(data, pvals, extra_info = True)

#Should give:
#(0.72288173732954353, 0.86138425703434118, 2.4580096358564503, 8.1366646038518677)
#(Pbrown,Pfisger,Scale_Factor C,DFbrown)




# TCGA dataset
#Pathways.tsv contains a list of 45 genes that belong to 3 pathways:
#--'FOXA1 TRANSCRIPTION FACTOR NETWORK', 'SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION', 'GLYPICAN 3 NETWORK'
#CDH4_Pvalues.tsv contains P values from the spearman correlation between CHD4 and each of the 45 genes from TCGA GBM [data from feb 12th 2013].
#The P-values for each set of genes in each pathway are combined using our method or fishers method

pathways = ['FOXA1 TRANSCRIPTION FACTOR NETWORK', 'SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION', 'GLYPICAN 3 NETWORK']
gene_sublist = []
for p in pathways:
    gene_sublist+=PP.AllPathwayGeneDict[p]

gene_sublist = list(set(gene_sublist))
PathwayGeneDict = {p:[] for p in pathways}
#load pathways:
f = open("../Data/pathways.tsv")
f.readline()
gene_list = []
for line in f:
    L = line.replace("\n", "").split("\t")
    PathwayGeneDict[L[0]].append(L[1])
    gene_list.append(L[1])

f.close()
gene_list = list(set(gene_list))

PValueDict = {}
f = open("../Data/CDH4_Pvalues.tsv")
f.readline()
for line in f:
    L = line.replace("\n", "").split("\t")
    PValueDict[L[0]] = float(L[1])
f.close()
GeneData = {}
FM = open("../Data/ReducedFeatureMatrix.tsv")
for line in FM:
    L = line.replace("\n", "").split("\t")
    GeneData[L[0]] = [float(l) for l in L[1:]]
FM.close()

for p in pathways:
    print "pathway:", p
    DataMatrix = array([GeneData[g] for g in PathwayGeneDict[p] if g in GeneData])
    Pvalues = array([PValueDict[g] for g in PathwayGeneDict[p] if g in PValueDict])
    print EmpiricalBrownsMethod(DataMatrix, Pvalues, extra_info = True)


#Should give:
#pathway: FOXA1 TRANSCRIPTION FACTOR NETWORK
#(7.7778969794178595e-53, 4.043406925735029e-139, 2.7193665607584965, 21.328496436251836)
#pathway: SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION
#(1.6980563950404756e-41, 6.4438388244313223e-45, 1.0877310573657077, 18.386897997043924)
#pathway: GLYPICAN 3 NETWORK
#(4.8216794064099692e-07, 1.4387321406058163e-08, 1.2976927497874169, 10.788378067376447)
#(Pbrown,Pfisger,Scale_Factor C,DFbrown)



