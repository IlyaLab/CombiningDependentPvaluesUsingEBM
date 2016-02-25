#Copyright (c) 2015, Institute for Systems Biology
#
#Permission is hereby granted, free of charge, to any person obtaining
#a copy of this software and associated documentation files (the
#"Software"), to deal in the Software without restriction, including
#without limitation the rights to use, copy, modify, merge, publish,
#distribute, sublicense, and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to
#the following conditions:

#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Author: William Poole
# Ported to R: David L Gibbs
# Email: dgibbs@systemsbiology.org / william.poole@systemsbiology.org / tknijnen@systemsbiology.org
# Created: June 2015

pop.var <- function(x) var(x) * (length(x)-1) / length(x)
pop.sd <- function(x) sqrt(pop.var(x))


#Input: raw data vector (of one variable) with no missing samples. May be a list or an array.
#Output Transforemd data vector w.
transformData <- function(data_vector) {
    dvm = mean(data_vector)
    dvsd = pop.sd(data_vector)
    s = (data_vector-dvm)/dvsd
    distr = ecdf(s)
    sapply(s, function(a) -2*log(distr(a)))
}


#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array.
#       Note: Method does not deal with missing values within the data.
#Output: An m x m matrix of pairwise covariances between transformed raw data vectors
calculateCovariances <- function(data_matrix){
    transformed_data_matrix = apply(data_matrix, MARGIN=1, FUN=transformData)
    covar_matrix = cov(transformed_data_matrix)
    covar_matrix
  }


#Input: A m x m numpy array of covariances between transformed data vectors and a vector of m p-values to combine.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
combinePValues <- function(covar_matrix, p_values, extra_info = FALSE){
    N = ncol(covar_matrix) # number of samples
    df_fisher = 2.0*N
    Expected  = 2.0*N
    cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=FALSE)]))
    Var = 4.0*N+cov_sum
    c = Var/(2.0*Expected)
    df_brown = (2.0*Expected^2)/Var
    if (df_brown > df_fisher) {
        df_brown = df_fisher
        c = 1.0
    }
    x = 2.0*sum( -log(p_values) )

    p_brown = pchisq(df=df_brown, q=x/c, lower.tail=FALSE)
    p_fisher = pchisq(df=df_fisher, q=x, lower.tail=FALSE)

    if (extra_info) {
        return(list(P_test=p_brown, P_Fisher=p_fisher, Scale_Factor_C=c, DF=df_brown))
    }
    else {
        return(p_brown)
    }
}

#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
#       A vector of m P-values to combine. May be a list or of type numpy.array.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
empiricalBrownsMethod <- function(data_matrix, p_values, extra_info = FALSE) {
  # inputs must be numeric
    covar_matrix = calculateCovariances(data_matrix)
    return(combinePValues(covar_matrix, p_values, extra_info))
}

#

#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numeric matrix
#       A numeric vector of m P-values to combine.
#Output: A combined P-value using Kost's Method.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
kostsMethod <- function(data_matrix, p_values, extra_info = FALSE) {
    covar_matrix <- calculateKostCovariance(data_matrix)
    combinePValues(covar_matrix, p_values, extra_info = extra_info)
}

#Input correlation between two n x n data vectors.
#Output: Kost's approximation of the covariance between the -log cumulative distributions. This is calculated with a cubic polynomial fit.
kostPolyFit <- function(cor) {
    a1 <- 3.263
    a2 <- 0.710
    a3 <- 0.027 #Kost cubic coeficients
    (a1*cor + a2*cor^2 + a3*cor^3)
}

#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numeric matrix.
#       Note: Method does not deal with missing values within the data.
#Output: An m x m matrix of pairwise covariances between the data vectors calculated using Kost's polynomial fit and the pearson correlation function.
calculateKostCovariance <- function(data_matrix) {
    m = nrow(data_matrix)
    covar_matrix = mat.or.vec(m, m)
    for (i in 1:m) {
        for (j in i:m) {
            res0 <- cor.test(data_matrix[i,], data_matrix[j,])
            cor <- res0$estimate
            p_val <- res0$p.value
            covar = kostPolyFit(cor)
            covar_matrix[i, j] = covar
            covar_matrix[j, i] = covar
          }
        }
    return(covar_matrix)
}
