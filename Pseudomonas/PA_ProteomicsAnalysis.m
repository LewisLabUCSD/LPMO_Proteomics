%% Create data strucure
fileName = 'PA_proteome';
sheetName = 'Proteome'; 
sampleNumber = 6;
Other=["Protein IDs"; "Fasta headers"; "Gene names";"Protein names"]; %other column information for import
norm='none'; %type of normalization: z_genes, z_samples, quantile,none
[data] = DataStructure(fileName, sheetName, sampleNumber, Other, norm);

figure()
boxplot(data.allData)
title('Data Distribution');
ylabel(' Log2(LFQ)')
xticklabels(data.sampleNames_txt);
set(gca,'XTickLabelRotation',45)
% Note: data looks good without any further normalization

%% Add geneID conversions
fileName = 'refSeq_uniprot';
[data] = addIDs(fileName,data);

%Fill empty KEGG entries with original protein ID with prefix "noKEGG"
data.KEGG_singular(cellfun(@isempty,data.KEGG_singular)) = strcat('noKEGG_',data.Protein_IDs(cellfun(@isempty,data.KEGG_singular)));

%%%%% Note that we now have two verions of gene names and protein names
%%%%% Determine if these overlap, and combine them to create comprehensive gene and protein name lists

%Combine the gene names in the original data and the conversion gene names for a more comprehensive list of gene names
nonempty_Gene_names = find(~cellfun(@isempty,data.Gene_names));
nonempty_geneNames = find(~cellfun(@isempty,data.geneNames));
intersect(nonempty_Gene_names,nonempty_geneNames) %Note:there is no overlap
data.geneNames_merge = data.geneNames;
for i =1:length(nonempty_Gene_names)
    data.geneNames_merge{nonempty_Gene_names(i)} = data.Gene_names{nonempty_Gene_names(i)};
end

%Combine the protein names in the original data and the conversion protein names for a more comprehensive list of protein names
nonempty_Protein_names = find(~cellfun(@isempty,data.Protein_names));
nonempty_proteinNames = find(~cellfun(@isempty,data.proteinNames));
intersect(nonempty_Protein_names,nonempty_proteinNames) %Note:there is no overlap
data.proteinNames_merge = data.proteinNames;
for i =1:length(nonempty_Protein_names)
    data.proteinNames_merge{nonempty_Protein_names(i)} = data.Protein_names{nonempty_Protein_names(i)};
end

%% add virulence annotations
[data] = addVirulence(data);
%% Create comparison titles for plots and structure names
Media = {'LB'; 'RPMIserum'; 'RPMI'}; % three conditions
Strains = unique(data.sampleNames_txt,'stable');
titles = {}; %initialize cell for plot titles
structFields = {}; %initialize cell for structure field names
for i = 1:length(Media)
    for j = 2:length(Strains)
        titles{end+1} = strcat(Strains{j},{' vs. '},Strains{1},{' in '},Media{i});
        structFields{end+1} = strcat(Media{i},'_',Strains{j});
    end
end
% Structure field names can't handle delta symbol
structFields = strrep(structFields,'?','del');

%% human vs pseudomonas
% some of the MS results mapped to human genes (and/or genes that are annotated as neither PA nor human
% retrieve indices of human vs PA vs other for future filtering
[human_idx, PA_idx, other_idx] = group_humanVSpseud(data);

%% Pairwise ttest and FDR adjustment
[Ttest] = F_Ttest_PA(data);

%% Volcano plots and significant differential expression lists

q_thresh = 0.05; %FDR threshold for significance
[sigList, otherList,meanMatrix] = newVolcanoPlots_PA(data,titles, structFields,Ttest,human_idx,q_thresh);

% In addition to significant proteins that meet the threshold criteria, we
% also want to define proteins that are absent in one sample (not produced
% or so low that they are not detected) and present in the other as a
% differnt class of significant and add them to the sig lists
[otherList,sigList] = addPresentAbsent_v2(data, otherList,structFields,sigList,human_idx);

%% Pathway enrichment
PosNeg = {'KEGG_pos','KEGG_neg'};
cuttOff = 8;
for i= 1:length(fields(sigList)) % for each media condition
    p_score = [];
    p_lab = [];
    p_dir = [];
    for j = 1:2 % perform enrichment for upregulated and downregulated pathways (pos/neg)
       pos_OR_neg = PosNeg{j};
       [sigPath_sort,sort_score,MpathSig,pathway]=pathwayEnrichment_v3('refSeq_uniprot',data,i,sigList,titles,pos_OR_neg,cuttOff,q_thresh,meanMatrix);
       if j == 1
           dir = repmat('pos',length(sort_score),1);
       else
           dir = repmat('neg',length(sort_score),1);
       end
    p_score = [p_score;sort_score'];
    p_lab = [p_lab;sigPath_sort];
    p_dir = [p_dir;dir];
    end
    %save enrichment results in structures
    pEnrichment.(structFields{i}).pathway = p_lab;
    pEnrichment.(structFields{i}).direction = p_dir;
    pEnrichment.(structFields{i}).score = p_score;
    pathwayMatrix.(structFields{i}).M = MpathSig; %sig protein binary matrix
    pathwayMatrix.(structFields{i}).pathway = pathway; %pathways corresponding to binary matrix columns
end
