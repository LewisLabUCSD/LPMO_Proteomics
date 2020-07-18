function [sigPath_sort,sort_score,MpathSig,pathway] = pathwayEnrichment_v3(fileName,data,analysisID,sigList,titles,pos_OR_neg,cuttOff,q_thresh,meanMatrix)
% Hypergeometric enrichment of KEGG pathways

% Inputs
%   fileName -- file containing KEGGid-Pathway mapping
%   data -- data structure
%   analysisID -- id of specific condition (from sigList fields, corresponds to media)
%   sigList -- data structure with significant gene information separated by neg vs pos FC
%   titles -- condition specific titles
%   pos_OR_neg -- analysis of up or down regulated pathways
%   cuttOff -- only include pathways have >= cuttOFF proteins involved
%   q_thresh -- FDR threshold for significance
%   meanMatrix -- matrix with mean sample values


% Output
%   sigPath_sort -- significant pathways sorted based on enrichment scores
%   sort_score -- sorted pathway enrichment scores
%   MpathSig -- binary matrix of sig proteins vs pathways with applied cuttoff (cuttOff)
%   pathway -- list of all pathways included in the analysis


%% Read in ID conversion files
[~,txt_KEGGpathwayID,~] = xlsread(strcat(fileName,'.xlsx'),'KEGG-to-pathwayID'); %KEGGid to pathway map
pathwayID = txt_KEGGpathwayID(2:end,1);
KEGGid = txt_KEGGpathwayID(2:end,2);

[~,txt_KEGGpathway,~] = xlsread(strcat(fileName,'.xlsx'),'pathwayID-to-pathway'); %pathwayID to pathway
pathwayID2 = txt_KEGGpathway(2:end,1);
pathway = txt_KEGGpathway(2:end,2);

%% For each condition (analysisID), create binary matrices of proteins vs pathways (1 if gene is present in pathway)

%%%%%%%%%% Binary matrix for all background proteins
%extract WT sample from meanMatrix
if analysisID ==1
    WT = meanMatrix(:,1);
elseif  analysisID == 2
    WT = meanMatrix(:,3);
else
    WT = meanMatrix(:,5);
end

%mean matrix of just KO samples (remove WT from meanMatrix)
meanMatrix_trim = meanMatrix;
meanMatrix_trim(:,[1,3,5])=[];

WT_prot = find(~isnan(WT)); % protein indices in WT sample
sample_prot = find(~isnan(meanMatrix_trim(:,analysisID))); % protein indices in KO sample
pre_proteinList = data.KEGG_singular(unique([WT_prot;sample_prot])); %background, all proteins that are present in either WT or KO
proteinList = pre_proteinList(~cellfun(@isempty,pre_proteinList)); %remove empty cells (those that were not able to map to KEGG ids)
%proteinList = intersect(proteinList,KEGGid); %only include proteins for which we have pathway info on

%initialize binary matrix (proteins vs pathways)
Mpathway = zeros(length(proteinList), length(pathwayID2)); 

%fill in matrix
for i = 1:length(proteinList)
    idx = find(strcmp(proteinList{i},KEGGid)); % indices of KEGGid-pathway pairs (protein may map to multiple pathways)
    for j=1:length(idx)
        genePathway = pathwayID{idx(j)};
        idx2 = find(strcmp(genePathway,pathwayID2)); % find corresponding pathway indices in map of pathway id to pathway
        Mpathway(i,idx2) = 1; % add 1 to binary matrix
    end
end



%%%%%%%%%% Binary matrix for significant proteins
sigList_fields = fields(sigList);
sigProteins = sigList.(sigList_fields{analysisID}).(pos_OR_neg); %list of sig proteins 
sigProteins_trimed = sigProteins(~cellfun(@isempty,sigProteins)); %remove empty (protiens that weren't able to map to KEGGid)
countEmpty_sig = length(sigProteins)-length(sigProteins_trimed); %allows you to count you losses due to mapping
fprintf('%d/%d sig genes  were not mapped to a KEGGids\n',countEmpty_sig,length(sigProteins))
% sigProteinList = intersect(sigProteins_trimed,KEGGid);
sigProteinList = sigProteins_trimed; %only use proteins that were able to map to KEGGids

%subset binary matrix of all background proteins to only include sig proteins
[~, idx_sig] = ismember(sigProteinList,proteinList); 
MpathSig = Mpathway(idx_sig,:);

%define variables for hypergeometric enrichment
X = sum(MpathSig); %vector with sum of sig proteins that belong to each pathway
M = size(Mpathway,1); %background proteins that were successfully convereted to KEGGids
K = sum(Mpathway); %vector with sum of background proteins that belong to each pathway
N = size(MpathSig,1); %sig proteins minus genes that were not converted to KEGGids

%hypergeometric enrichment with FDR correction
hyge_p_b4 = 1-hygecdf(X,M,K,N);
hyge_p = hyge_p_b4+hygepdf(X,M,K,N);
[hyge_q]=mafdr(hyge_p,'BHFDR',true);

figure()
subplot(2,1,1);
histogram(hyge_p);
title({titles{analysisID}{1};'Distribution of p-values'});
subplot(2,1,2);
histogram(hyge_q)
title({titles{analysisID}{1};'Distribution of q-values'});


if contains(pos_OR_neg,'pos')
    direction = 'Positively enriched: ';
elseif contains(pos_OR_neg,'neg')
    direction = 'Negatively enriched: ';
end
subtitle = strcat({direction},titles{analysisID}{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Redo analyisis with min protein cutoff for pathways

%% trim matrices to only include pathways that have >= 8 proteins involved
pw_idx = find(sum(Mpathway)>=cuttOff); %indices of patwhays that have 8 or more protein involved
Mpathway = Mpathway(:,pw_idx);
MpathSig = MpathSig(:,pw_idx);
pathway = pathway(pw_idx); %list of all pathways included in analysis

%% Hypergeometric with trimmed pathways
X = sum(MpathSig);
M = size(Mpathway,1); %all genes that were successfully convereted to uniprotIDs minus those that were not present in the ID-to-GOterm map (PAmap)
K = sum(Mpathway);
N = size(MpathSig,1); %sig proteins minus genes that were not converted to uniprotIDs minus those that were not present in the PAmap

hyge_p_b4 = 1-hygecdf(X,M,K,N);
hyge_p = hyge_p_b4+hygepdf(X,M,K,N);
[~, ~, ~, hyge_q]=fdr_bh(hyge_p);

figure()
subplot(2,1,1);
histogram(hyge_p);
title({titles{analysisID}{1};'Distribution of p-values'; 'Only including pathways with 8 or more proteins'});
subplot(2,1,2);
histogram(hyge_q)
title({titles{analysisID}{1};'Distribution of q-values'; 'Only including pathways with 8 or more proteins'});

%% Significantly enriched pathways

path_idx = find(hyge_q>0 & hyge_q<q_thresh); %apply FDR threshold to find sig pathways
sigPath_term = pathway(path_idx); %convert pathway ids to pathways

enrichmentScore=-log10(hyge_q(path_idx)); %calculate enrichment score (-log10(FDR))

%sort significant pathways based on enrichment score
[sort_score,sort_idx] = sort(enrichmentScore);
sigPath_sort = sigPath_term(sort_idx);
sigPath_sort = erase(sigPath_sort,'- Pseudomonas aeruginosa UCBPP-PA14'); %trim pathway labels


figure('Position', get(0, 'Screensize'));
figure()
fig=barh(sort_score);
set(gca,'ytick',1:length(sigPath_sort));
set(gca,'yticklabel',sigPath_sort)
xlabel({'Enrichment Score'; '-log10(q-value)'})
title({'Pathway Enrichment';subtitle{1}});

end
