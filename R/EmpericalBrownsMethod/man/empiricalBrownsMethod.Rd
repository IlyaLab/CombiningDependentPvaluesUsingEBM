% File inst/doc/empiricalBrownsMethod.Rd
\name{empiricalBrownsMethod}
\alias{empiricalBrownsMethod}
\title{EmpiricalBrownsMethod}
\description{Combine p-values}
\usage{
  empiricalBrownsMethod(data_matrix, p_values, extra_info)
}
\arguments{
  \item{data_matrix}{a}
  \item{p_values}{b}
  \item{extra_info}{c}
}
\examples{
## restore the saved values to the current environment
  data(testData)
#  glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"];
#  glypPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% glypGenes];
#  glypDat   <- dat[dat$V1 %in% glypGenes, 2:ncol(dat)];
#  print(empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=T));
}
\keyword{file}
