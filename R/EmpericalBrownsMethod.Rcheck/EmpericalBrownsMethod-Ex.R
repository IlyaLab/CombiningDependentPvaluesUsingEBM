pkgname <- "EmpericalBrownsMethod"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('EmpericalBrownsMethod')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("empiricalBrownsMethod")
### * empiricalBrownsMethod

flush(stderr()); flush(stdout())

### Name: empiricalBrownsMethod
### Title: The Empirical Browns Method For Combining P-values
### Aliases: empiricalBrownsMethod
### Keywords: multivariate

### ** Examples

## restore the saved values to the current environment
  data(ebmTestData)
  glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"]
  glypPvals <- allPvals$pvalue.with.CHD4[match(glypGenes, allPvals$gene)];
  glypDat   <- dat[match(glypGenes, dat$V1), 2:ncol(dat)];
  empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=TRUE);



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
