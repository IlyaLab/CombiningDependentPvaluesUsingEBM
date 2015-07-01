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
% Created: June 2015

function [Pbrown,Pfisher,Cbrown,DFbrown] = EmpiricalBrownsMethod(data_matrix, p_values)

%license, author, date

%% Inputes and outputs

% Input:  An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
%         A vector of m P-values to combine. May be a list or of type numpy.array.
% Output: A combined P-value.
%         If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method

%% Compute P-values
covar_matrix = CalculateCovariances(data_matrix);
[Pbrown,Pfisher,Cbrown,DFbrown] = CombinePValues(covar_matrix, p_values);

%% Functions

%Calculate Covariances
    function covar_matrix = CalculateCovariances(data_matrix);
        [m,n] = size(data_matrix);
        transformed_data_matrix = NaN(size(data_matrix));
        for f = 1:m;
            transformed_data_matrix(f,:) = TransformData(data_matrix(f,:));
        end
        covar_matrix = cov(transformed_data_matrix',0);
    end

% Transform data
    function transformed_data_vector = TransformData(data_vector);
        m = mean(data_vector);
        sd = std(data_vector,1);
        s = (data_vector-m)./sd;
        [~,~,i] = unique(s);
        z = ecdf(s); z(1) = []; %ecdf
        transformed_data_vector = -2*log(z(i));
    end

% Combining P-values
    function [p_brown,p_fisher,c,df_brown] = CombinePValues(covar_matrix, p_values);
        
        m = size(covar_matrix,1);
        df_fisher = 2.0*m;
        Expected = 2.0*m;
        cov_sum = sum(sum(covar_matrix))-sum(diag(covar_matrix));
        Var = 4.0*m+cov_sum;
        c = Var/(2*Expected);
        df_brown = 2*(Expected^2)/Var;
        if df_brown > df_fisher;
            df_brown = df_fisher;
            c = 1;
        end
        
        x = 2*sum(-log(p_values));
        p_brown = chi2cdf(1.0*x/c,df_brown,'upper');
        p_fisher = chi2cdf(1.0*x,df_fisher,'upper');
    end

end

