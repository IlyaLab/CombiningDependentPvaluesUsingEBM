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
### Title: EmpiricalBrownsMethod
### Aliases: empiricalBrownsMethod
### Keywords: file

### ** Examples

## restore the saved values to the current environment
  data(testData)
#  glypGenes <- pathways$gene[pathways$pathway == "GLYPICAN 3 NETWORK"];
#  glypPvals <- allPvals$pvalue.with.CHD4[allPvals$gene %in% glypGenes];
#  glypDat   <- dat[dat$V1 %in% glypGenes, 2:ncol(dat)];
#  print(empiricalBrownsMethod(data_matrix=glypDat, p_values=glypPvals, extra_info=T));



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
