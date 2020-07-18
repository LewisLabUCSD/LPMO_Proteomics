%% Create data strucure
fileName = 'mouseProteome';
totSamples = 19; %count each individual replicate as a sinlge sample
Other=["Protein IDs"; "Fasta headers";"Protein names";"Gene names"];
norm='none'; %z_genes, z_samples, quantile,none
norm2='none';
[data] = DataStructure_mouse(fileName, totSamples, Other, norm, norm2);

figure()
boxplot(data.allData)
title('Data Distribution');
ylabel(' Log2(LFQ)')
xticklabels(data.sampleNames_txt);
set(gca,'XTickLabelRotation',45)
% Note: data looks good without any further normalization

%% Pairwise ttest and FDR adjustment comparing WT and KO vs CTR
% F-test to determine equal variance
[Ttest] = F_Ttest_v2(data);

%% Volcano plots and significant differential expression lists
q_thresh = 0.05; %FDR threshold for significance
[sigList, otherList,foldChange] = volcanoPlots(data,q_thresh,Ttest);

%% Identify DE proteins that are unique and common to both WT and KO comparisons against control

commonSpleen = intersect(otherList.SpleenKO.GeneNames,otherList.SpleenWT.GeneNames);

unique_spleen = setxor(otherList.SpleenKO.ProteinIDs,otherList.SpleenWT.ProteinIDs);

[~,ia_spleen,~] = intersect(data.Protein_IDs,unique_spleen);

non_common.Spleen.genes = data.Gene_names(ia_spleen);
non_common.Spleen.proteinIDS = data.Protein_IDs(ia_spleen);
non_common.Spleen.KO_log2FC = foldChange.SpleenKO(ia_spleen);
non_common.Spleen.KO_q = Ttest.Spleen.KO.q(ia_spleen);
non_common.Spleen.WT_log2FC = foldChange.SpleenWT(ia_spleen);
non_common.Spleen.WT_q = Ttest.Spleen.WT.q(ia_spleen);



