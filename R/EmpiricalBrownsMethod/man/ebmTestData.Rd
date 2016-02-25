\name{ebmTestData}
\docType{data}
\alias{ebmTestData}
\alias{randData}
\alias{dat}
\alias{allPvals}
\alias{pathways}

\title{Data used in tests and examples.}

\description{
  This data is used in the unit tests and usage examples. There are four items:

  allPvals, dat, pathways, and randData. allPvals is a data.frame of p-values
  for the spearman correlation between CHD4 and each of the 45 genes.

  dat is the gene expression data corresponding to genes in allPvals.

  pathways is a data.frame listing gene membership for 3 biochemical pathways.

  randData is a gaussian generated data set,
  emphasizing dependence among variables. Independent Var [line 1] are 25
  samples from a unit normal distribution. Depedent Var 1-10 [line 2-11] are
  each 25 samples drawn from a 10 dimensional normal distribution centered at
  the origin with off diagonal terms a=0.25. The P values from a pearson
  correlation between the independent var and each dependent var are combined.
}
\value{
  data objects in the environment
}
\usage{data(ebmTestData)}
\format{Rdata object}
\source{GEO and generated.}
\keyword{datasets}
