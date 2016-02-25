---
title: "Empirical Browns Method"
author: "David L Gibbs"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Empirical Browns Method}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


## Abstract

This package provides an empirical adaptation of Brown’s Method, which is an extension of Fisher’s Method, for combining dependent P-values. This is appropriate for highly correlated data sets commonly found in high-throughput biological experiments.

In our corresponding paper (submitted) we show that Fisher’s Method is biased given dependent sets of P-values. When applied on the same data sets, the Empirical Brown’s Method provides a better null distribution and a more conservative result.

## Introduction

In order to integrate the large and diverse datasets found in systems biology, it is common to combine P-values from multiple statistical tests. The earliest method to combine independent P-values is seen in the work of Fisher (1). Brown later extended Fisher’s Method to the case where P-values are assumed to be drawn from a multivariate normal distribution with a known covariance matrix (2). Of all methods for combining P-values, Brown’s most simply combines equally weighted dependent P-values. However, instead of using numerical integration, our adaptation of Brown’s Method uses the empirical cumulative distribution function derived from the data, making our method dramatically more efficient and suitable for large omics data.

## Using the function

This function is used to combine the p-values of non-independent tests into one statistic.  The toy example here is to combine gene-gene correlation tests given a pathway, which is described as a set of genes. The CDH4 gene is thought to be important in the development of Glioblastoma cancer, so a correlation test was performed comparing the gene expression of a list of genes to the expression of CDH4. Then, these correlation p-values were combined on pathways, indicating that the pathway could be aberrant.

The method is highly general, and can be used flexibly.

```r
data(ebmTestData)

# first get the genes for the specified pathway
glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"]

# then gather the correlation test p-values
glypPvals <- allPvals$pvalue.with.CHD4[match(glypGenes, allPvals$gene)];

# the gene expression for genes in the given pathway
glypDat   <- dat[match(glypGenes, dat$V1), 2:ncol(dat)];

# the call to combine p-values
empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=TRUE);

# Alternatively, to use Kost's methods
kostsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=TRUE);
```

## References

[1] Fisher, R.A. (1948) Questions and answers /#14, The American Statistician, 2, 30-31.
[2] Brown, M. (1975) A method for combining non-independent, one-sided tests of significance, Biometrics, 31, 987992.
[3] Kost J.T. et al. (2002) Combining dependent P-values, Statistics & Probability Letters, 60, 183190.
[4] Whitlock M.C. (2005) Combining probability from independent tests: the weighted Z-method is superior to Fishers approach, Journal of Evolutionary Biology, 18, 1368-1373.
[5] Aerts S. et al. (2006) Gene prioritization through genomic data fusion, Nature Biotechnology, 24, 537 - 544.
[6] Zaykin D.V. et al. (2002) Truncated product method for combining P-values, Genet Epidemiol. 22, 170-85.
[7] The Cancer Genome Atlas Research Network. (2008). Comprehensive genomic characterization defines human glioblastoma genes and core pathways, Nature 455, 1061-1068.
[8] Schaefer C.F. et al. (2009) The Pathway Interaction Database, Nucleic Acids Res. 37, 674-679.
