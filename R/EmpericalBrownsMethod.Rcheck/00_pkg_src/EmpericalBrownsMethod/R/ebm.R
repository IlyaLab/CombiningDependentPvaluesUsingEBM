# Copyright 2015, Institute for Systems Biology.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
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
combinePValues <- function(covar_matrix, p_values, extra_info = F){
    N = ncol(covar_matrix) # number of samples
    df_fisher = 2.0*N
    Expected  = 2.0*N
    cov_sum <- (2*sum(covar_matrix[lower.tri(covar_matrix, diag=F)]))
    Var = 4.0*N+cov_sum
    c = Var/(2.0*Expected)
    df_brown = (2.0*Expected^2)/Var
    if (df_brown > df_fisher) {
        df_brown = df_fisher
        c = 1.0
    }
    x = 2.0*sum( -log(p_values) )

    p_brown = pchisq(df=df_brown, q=x/c, lower.tail=F)
    p_fisher = pchisq(df=df_fisher, q=x, lower.tail=F)

    if (extra_info) {
        return(list(P_Brown=p_brown, P_Fisher=p_fisher, Scale_Factor_C=c, DF_Brown=df_brown))
    }
    else {
        return(p_brown)
    }
}

#Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
#       A vector of m P-values to combine. May be a list or of type numpy.array.
#Output: A combined P-value.
#        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
empiricalBrownsMethod <- function(data_matrix, p_values, extra_info = F) {
  # inputs must be numeric
    covar_matrix = calculateCovariances(data_matrix)
    return(combinePValues(covar_matrix, p_values, extra_info))
}

#
