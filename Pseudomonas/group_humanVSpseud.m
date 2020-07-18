function [human_idx, PA_idx, other_idx] = group_humanVSpseud(data)
%MS returned some proteins that are not pseudomonas
%Filter data to determine which proteins are: pseudomonas, human, other

% Inputs
%   data -- data strucrure

% Output
%   human_idx -- indices of human proteins
%   PA_idx -- indices of pseudomonas proteins
%   other_idx -- indices of proteins that were neither taged as neither human nor pseudomonas

%% 

%extract headers that contrain "human" or "pseudomona
isHuman = contains(data.Fasta_headers,'human','IgnoreCase',true);
isPA = contains(data.Fasta_headers,'Pseudomonas','IgnoreCase',true);

% indices of stritcly human or pseudomonas genes 
human_idx = find(isHuman==1);
PA_idx = find(isPA==1);

% indices of proteins that are neigher strictly human nor strictly pseudomonas
human_PA = cat(1,find(isHuman==1),find(isPA==1)); %human plus pseudo indices
allIDX = 1:1:length(data.Fasta_headers); %array of all possible indices
other_idx = setdiff(allIDX,human_PA);

end