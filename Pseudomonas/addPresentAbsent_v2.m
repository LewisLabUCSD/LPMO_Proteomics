function [otherList,sigList] = addPresentAbsent_v2(data, otherList,structFields,sigList,human_idx)
% Add genes that are present/absent to the sig lists

% Inputs
%   data -- data structure
%   otherList -- data structure with significant gene information NOT separated by neg vs pos FC
%   structFields -- number of samples (excluding replicates)
%   sigList -- data structure with significant gene information separated by neg vs pos FC
%   human_idx -- Indices of human genes to filter out

% Output
%   otherList -- updated otherList with present/absent sig proteins
%   sigList -- updated otherlist with present/absent sig proteins

%%

%For each sample, sum # NON-nans per row/protein across the three replictes
for i = 1:6
    nanMatrix(:,i) = sum(~isnan(data.(strcat('sample_',num2str(i)))),2);
end

strainNumber = 2; % 2 strains (WT and KO) per media 
for i=0:2 %for all 3 media types
    WTloc = i*strainNumber+1; %index of WT in nanMatrix: 1,3,5
    WT = nanMatrix(:,WTloc);
    KO = nanMatrix(:,WTloc+1); 
    WT_detect = find(WT>=2); % indices of proteins that were detected in 2 or more of the 3 replicates in WT
    WT_nan = find(WT==0); % indices of proteins that were detected in none of the 3 replicates in WT
    KO_detect = find(KO>=2); % indices of proteins that were detected in 2 or more of the 3 replicates in KO
    KO_nan = find(KO==0); % indices of proteins that were detected in none of the 3 replicates in KO
    dwn_KO = WT_detect(ismember(WT_detect,KO_nan)); % indices of proteins that were detected in WT and absent in KO
    up_KO = KO_detect(ismember(KO_detect,WT_nan)); % indices of proteins that were detected in KO and absent in WT

    % trim proteins that are human
    dwn_KO = setdiff(dwn_KO,human_idx);
    up_KO=setdiff(up_KO,human_idx);

    % create vectors indicating wich proteins are "up" vs "down" regulated in KO compared to WT
    % note these proteins will not have an associated q value, therefore create a vector of "NaN" for their associated sig statistic
    FC_up = repmat({'up'},length(up_KO),1);
    FC_dwn = repmat({'down'},length(dwn_KO),1);
    q_nan = repmat({'NaN'},length(up_KO)+length(dwn_KO),1);

    %Add up/down info to sig lists
    otherList.(structFields{i+1}).GeneNames = [otherList.(structFields{i+1}).GeneNames;data.geneNames_merge(up_KO);data.geneNames_merge(dwn_KO)];
    otherList.(structFields{i+1}).ProteinIDs = [otherList.(structFields{i+1}).ProteinIDs;data.Protein_IDs(up_KO);data.Protein_IDs(dwn_KO)];
    otherList.(structFields{i+1}).ProteinNames = [otherList.(structFields{i+1}).ProteinNames;data.proteinNames_merge(up_KO);data.proteinNames_merge(dwn_KO)];
    otherList.(structFields{i+1}).logFC = [cellstr(num2str(otherList.(structFields{i+1}).logFC));FC_up;FC_dwn];
    otherList.(structFields{i+1}).q = [cellstr(num2str(otherList.(structFields{i+1}).q));q_nan];
    otherList.(structFields{i+1}).Fasta_headers = [otherList.(structFields{i+1}).Fasta_headers;data.Fasta_headers(up_KO);data.Fasta_headers(dwn_KO)];
    otherList.(structFields{i+1}).virulenceFactor = [otherList.(structFields{i+1}).virulenceFactor;data.virulenceFactor(up_KO);data.virulenceFactor(dwn_KO)];

    sigList.(structFields{i+1}).KEGG_pos = [sigList.(structFields{i+1}).KEGG_pos;data.KEGG_singular(up_KO)];
    sigList.(structFields{i+1}).KEGG_neg = [sigList.(structFields{i+1}).KEGG_neg;data.KEGG_singular(dwn_KO)];
end
end
    