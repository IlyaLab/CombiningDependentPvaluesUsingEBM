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

%% Load TCGA dataset
%pathways
fid = fopen('../Data/pathways.tsv');
PW = textscan(fid,'%s%s','Headerlines',1,'Delimiter','\t');
fclose(fid);
Pathways = unique(PW{1});
for p = 1:length(Pathways);
    pos = find(strcmp(Pathways{p},PW{1}));
    Pathways{p,2} = unique(PW{2}(pos));
end

%genes and p-values
fid = fopen('../Data/CDH4_Pvalues.tsv');
CDH4pvalues = textscan(fid,'%s%n','Headerlines',1,'Delimiter','\t');
fclose(fid);

%expression data
fid = fopen('../Data/ReducedFeatureMatrix.tsv');
NL = 10;
tline = cell(NL,1);
not = zeros(NL,1);
for n = 1:10;
    tline{n} = fgetl(fid);
    display(tline{n});
    not(n) = sum(double(tline{n})==9);
end
fclose(fid);

unot = unique(not);
if length(unot)~=1;
    error(num2str(unot));
else
    rs = [];
    for n = 1:unot+1;
        if n>1;
            rs = cat(2,rs,'%n');
        else
            rs = cat(2,rs,'%s');
        end
    end
end

%read in file
fid = fopen('../Data/ReducedFeatureMatrix.tsv');
E = textscan(fid,rs,'Delimiter','\t','Headerlines',0);
fclose(fid);

Genes = E{1};
ExpressionData = NaN(length(Genes),unot);
for n = 1:unot;
    ExpressionData(:,n) = E{n+1};
end

save('../Data/TCGAdata.mat','Pathways','CDH4pvalues','ExpressionData','Genes');