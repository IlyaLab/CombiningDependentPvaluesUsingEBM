\name{calculateCovariances}
\alias{calculateCovariances}
\title{Calculate pairwise covariances between transformed raw data vectors}
\description{Calculate pairwise covariances between transformed raw data vectors}
\usage{
  calculateCovariances(data_matrix)
}
\arguments{
  \item{data_matrix}{An m x n numeric matrix with m variables in rows and n
                      samples in columns.}
}
\value{
  The output is an m x m matrix of calculated covariances
}
\examples{
  data(ebmTestdata)
  glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"]
  glypDat   <- as.matrix(dat[match(glypGenes, dat$V1), 2:ncol(dat)])
  calculateCovariances(glypDat)
}