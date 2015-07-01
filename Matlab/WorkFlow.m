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

%% Initialize
Init

%% ARTIFICIAL DATASET

%% Load artificial dataset
D = importdata('../Data/RandomData_withoutheaders.tsv');

%% Pearson correlation
[R,P] = corr(D(1,:)',D(2:end,:)');

%% Emperical Browns Methods
data_matrix = D(2:end,:);
p_values = P;
[Pbrown,Pfisher,Cbrown,DFbrown] = EmpiricalBrownsMethod(data_matrix, p_values)

%% Should give...

% Pbrown =
% 
%     0.7229
% 
% 
% Pfisher =
% 
%     0.8614
% 
% 
% Cbrown =
% 
%     2.4580
% 
% 
% DFbrown =
% 
%     8.1367


%% Load TCGA dataset
% parseTCGAdata
load('../Data/TCGAdata.mat','Pathways','CDH4pvalues','ExpressionData','Genes');

for p = 1:size(Pathways,1);
    genes_in_pathway = Pathways{p,2};
    [check,pos] = ismember(genes_in_pathway,CDH4pvalues{1});
    if ~all(check==1);error('error');end
    p_values = CDH4pvalues{2}(pos);
    [check,pos] = ismember(genes_in_pathway,Genes);
    if ~all(check==1);error('error');end
    data_matrix = ExpressionData(pos,:);
    
    display(Pathways{p,1});
    [Pbrown,Pfisher,Cbrown,DFbrown] = EmpiricalBrownsMethod(data_matrix, p_values)
    
end


%% Should give...
% FOXA1 TRANSCRIPTION FACTOR NETWORK
% 
% Pbrown =
% 
%    7.7779e-53
% 
% 
% Pfisher =
% 
%   4.0434e-139
% 
% 
% Cbrown =
% 
%     2.7194
% 
% 
% DFbrown =
% 
%    21.3285
% 
% GLYPICAN 3 NETWORK
% 
% Pbrown =
% 
%    4.8217e-07
% 
% 
% Pfisher =
% 
%    1.4387e-08
% 
% 
% Cbrown =
% 
%     1.2977
% 
% 
% DFbrown =
% 
%    10.7884
% 
% SUMOYLATION BY RANBP2 REGULATES TRANSCRIPTIONAL REPRESSION
% 
% Pbrown =
% 
%    1.6981e-41
% 
% 
% Pfisher =
% 
%    6.4438e-45
% 
% 
% Cbrown =
% 
%     1.0877
% 
% 
% DFbrown =
% 
%    18.3869














