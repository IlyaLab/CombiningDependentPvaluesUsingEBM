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
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import leaves_list
#import RandomPathwayParser as RPP

print "parsing pathways"
PathwayGeneDict = {}
GenePathwayDict = {}
f = open( "../Data/only_NCI_Nature_ver4.tab")

for line in f:
    ind1 = line.index("\t")
    pathway = line[:ind1]
    ind2 = line[ind1+1:].index("\t")+ind1+1
    genes = []
    iterate = True
    ind1 = ind2
    while (iterate):
        ind1 = ind2
        try:
            ind2 = line[ind1+1:].index(",")+ind1+1
            gene = line[ind1+1:ind2]
            genes.append(gene)
            if gene not in GenePathwayDict:
                GenePathwayDict[gene] = [pathway]
            else:
                GenePathwayDict[gene].append(pathway)
            if ind2 == ind1:
                iterate = False 
        except ValueError:
            iterate = False
    PathwayGeneDict[pathway]=genes
f.close()


#Sort pathway list using hierarchical clustering

global_pathway_list = [p for p in PathwayGeneDict.keys() if len(PathwayGeneDict[p])>=1]
L = len(global_pathway_list)
DM = zeros((L, L))
for i in range(L):
    for j in range(L):
        p1 = global_pathway_list[i]
        p2 = global_pathway_list[j]
        intersection = len([g for g in PathwayGeneDict[p1] if g in PathwayGeneDict[p2]])
        union = len(set(PathwayGeneDict[p1]+PathwayGeneDict[p2]))
        DM[i, j] = 1.0-1.0*intersection/union
DM_c = squareform(DM)
linkageMatrix = linkage(DM_c)

pathway_list_clustering = leaves_list(linkageMatrix)
sorted_pathway_list = array(global_pathway_list)[pathway_list_clustering]


pw_ontology = open( "../Data/pathway_ontology_annotations.txt")
PathwayParentDict = {}
PathwayChildDict = {}
parent_list = []
child_list = []
leaf_list =  []
all_pathway_list = []
for line in pw_ontology:
    L = line.replace("\n", "").replace("\r", "").split("\t")
    pathway = str.upper(L[3])
    all_pathway_list.append(pathway)
    parent = str.upper(L[4])
    if pathway in PathwayParentDict:
        PathwayParentDict[pathway].append(parent)
    else:
        PathwayParentDict[pathway]=[parent]

    if parent != '':
        child_list.append(pathway)
    if pathway in parent_list and parent != '':
        parent_list.remove(pathway)
    if parent not in child_list and parent not in parent_list:
        parent_list.append(parent)
    if parent not in PathwayChildDict:
        PathwayChildDict[parent] = [pathway]
    else:
        PathwayChildDict[parent].append(pathway)

AllPathwayGeneDict = {}
AllGenePathwayDict = {}

def recursive_gene_pathway_merger(pathway_list):
    genes = []
    for pathway in pathway_list:
        if pathway in PathwayGeneDict:
            sub_genes = list(set(PathwayGeneDict[pathway]))
            if pathway in PathwayChildDict:
                recursive_gene_pathway_merger(PathwayChildDict[pathway])
        elif pathway in PathwayChildDict:
            sub_genes = recursive_gene_pathway_merger(PathwayChildDict[pathway])
        else:
            sub_genes = []
        genes += sub_genes
        AllPathwayGeneDict[pathway]=sub_genes
        
    return list(set(genes))

leaf_list = [p for p in PathwayParentDict if p not in PathwayChildDict]
recursive_gene_pathway_merger([''])

for pathway in [p for p in AllPathwayGeneDict if p != ""]:
    genes = AllPathwayGeneDict[pathway]
    for gene in genes:
        if gene in AllGenePathwayDict:
            #print "can't be?"
            AllGenePathwayDict[gene].append(pathway)
        else:
            AllGenePathwayDict[gene]=[pathway]
            #print "could be?"

original_leaves = [p for p in PathwayGeneDict if p in leaf_list and len(PathwayGeneDict[p])>0]



def pathway_gene_overlap(pathway_list):
    overlap_matrix = zeros((len(pathway_list), len(pathway_list)))
    for i in range(len(pathway_list)):
        for j in range(i+1, len(pathway_list)):
            p1 = pathway_list[i]
            p2 = pathway_list[j]
            gl1 = AllPathwayGeneDict[p1]
            gl2 = AllPathwayGeneDict[p2]
            if len(gl1+gl2)>0:
                overlap_dist = 1.0*len([g for g in gl1 if g in gl2])/(len(gl1)+len(gl2))
                overlap_matrix[i, j] = overlap_dist
    print pathway_list
    return overlap_matrix
    

#Create Random Pathway original pathway mapping
"""
PRL = [(len(RPP.PathwayGeneDict[p]), p) for p in RPP.PathwayGeneDict if len(RPP.PathwayGeneDict[p])>0]
PL = [(len(PathwayGeneDict[p]), p) for p in PathwayGeneDict if len(PathwayGeneDict[p]) > 0]
PRL.sort()
PL.sort()
RandomPathwayMapping = {}
for i in range(len(PL)):
    pathway = PL[i][1]
    random_pathway = PRL[i][1]
    RandomPathwayMapping[pathway] = random_pathway
    RandomPathwayMapping[random_pathway] = pathway
"""


