% Copyright 2015, Institute for Systems Biology.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
% http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
% 
% Author: William Poole
% Ported to Matlab: Theo Knijnenburg
% Email: william.poole@systemsbiology.org / tknijnen@systemsbiology.org
% Created: February 2016

function [Pkost,Pfisher,Ckost,DFkost] = KostsMethod(data_matrix, p_values)

%license, author, date

%% Inputes and outputs

% Input:  An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
%         A vector of m P-values to combine. May be a list or of type numpy.array.
% Output: A combined P-value.
%         If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method

%% Compute P-values
covar_matrix = CalculateCovariances(data_matrix);
[Pkost,Pfisher,Ckost,DFkost] = CombinePValues(covar_matrix, p_values);

%% Functions

%Calculate Covariances
    function covar_matrix = CalculateCovariances(data_matrix);
        [m,n] = size(data_matrix);
        cor = corr(data_matrix');
        %Kost's polynomial fit
        a1 = 3.263;
        a2 = 0.710; 
        a3 = 0.027;
        covar_matrix = a1.*cor+a2.*cor.^2+a3.*cor.^3;
    end

% Combining P-values
    function [p_kost,p_fisher,c,df_kost] = CombinePValues(covar_matrix, p_values);
        
        m = size(covar_matrix,1);
        df_fisher = 2.0*m;
        Expected = 2.0*m;
        cov_sum = sum(sum(covar_matrix))-sum(diag(covar_matrix));
        Var = 4.0*m+cov_sum;
        c = Var/(2*Expected);
        df_kost = 2*(Expected^2)/Var;
        if df_kost > df_fisher;
            df_kost = df_fisher;
            c = 1;
        end
        
        x = 2*sum(-log(p_values));
        p_kost = chi2cdf(1.0*x/c,df_kost,'upper');
        p_fisher = chi2cdf(1.0*x,df_fisher,'upper');
    end

end

