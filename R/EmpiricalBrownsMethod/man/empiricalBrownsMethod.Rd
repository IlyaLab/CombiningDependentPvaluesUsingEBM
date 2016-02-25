% File inst/doc/empiricalBrownsMethod.Rd
\name{empiricalBrownsMethod}
\alias{empiricalBrownsMethod}
\title{The Empirical Browns Method For Combining P-values}
\description{Combining P-values from multiple statistical tests is common in
bioinformatics. However, this procedure is non-trivial for dependent P-values.
This package provides an empirical adaptation of Brown’s Method (an extension
of Fisher’s Method) for combining dependent P-values which is appropriate for
highly correlated data sets found in high-throughput biological experiments.}
\usage{
  empiricalBrownsMethod(data_matrix, p_values, extra_info)
}
\arguments{
  \item{data_matrix}{An m x n numeric matrix with m variables in rows and n
                      samples in columns.}
  \item{p_values}{A numeric vector of p-values with length m.}
  \item{extra_info}{boolean, TRUE additionally returns the p-value from
                    Fisher's method, the scale factor c, and the new degrees of
                    freedom from Brown's Method }
}
\value{
  The output is a list containing
  list(P_Brown=p_brown, P_Fisher=p_fisher, Scale_Factor_C=c, DF_Brown=df_brown)
  \item{P_Brown}{ p-value for Brown's method}
  \item{P_Fisher}{p-value for Fisher's method}
  \item{Scale_Factor}{the scale factor c}
  \item{DF_Brown}{the degrees of freedom used in Brown's method}
}
\examples{
## restore the saved values to the current environment
  data(ebmTestData)
  glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"]
  glypPvals <- allPvals$pvalue.with.CHD4[match(glypGenes, allPvals$gene)];
  glypDat   <- dat[match(glypGenes, dat$V1), 2:ncol(dat)];
  empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=TRUE);
}
\keyword{multivariate}
