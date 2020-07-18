load('proteomeAnalysis')

%% Create new data element that groups together data by sample
for i=1:2 % separate WT vs KO replicates
    data.LB_bySample.(strcat('s',num2str(data.sampleNames_num(i)))) = data.LB(:,(((i-1)*3)+1):(i*3));
    data.RPMIs_bySample.(strcat('s',num2str(data.sampleNames_num(i)))) = data.RPMIs(:,(((i-1)*3)+1):(i*3));
    data.RPMI_bySample.(strcat('s',num2str(data.sampleNames_num(i)))) = data.RPMI(:,(((i-1)*3)+1):(i*3));
end

%% Impute
newFields = fields(data.LB_bySample);
%initialize matrices to hold imputed values
impMatrix_LB = [];
impMatrix_RPMIs = [];
impMatrix_RPMI = [];
for k = 1:2 % for WT and KO
    %impute across replicates in LB
    LBimp = data.LB_bySample.(newFields{k});
    LBimp = knnimpute(LBimp);
    impMatrix_LB = [impMatrix_LB,LBimp];
    %impute across replicates in RPMIs
    RPMIsimp = data.RPMIs_bySample.(newFields{k});
    RPMIsimp = knnimpute(RPMIsimp);
    impMatrix_RPMIs = [impMatrix_RPMIs,RPMIsimp];
    %impute across replicates in RPMI
    RPMIimp = data.RPMI_bySample.(newFields{k});
    RPMIimp = knnimpute(RPMIimp);
    impMatrix_RPMI = [impMatrix_RPMI,RPMIimp];
end

%% Trim human proteins
impMatrix_LB(human_idx,:)=[];
impMatrix_RPMIs(human_idx,:)=[];
impMatrix_RPMI(human_idx,:)=[];
impMatrix_all = [impMatrix_LB impMatrix_RPMIs impMatrix_RPMI];
%vector corresponding to indices of proteins left
gene_indices=[1:1:length(data.Protein_IDs)]';
gene_indices(human_idx)=[];

%% If some genes are NaN across all samples, remove them from the analysis
nanSum = sum(isnan(impMatrix_all),2);
keepRow = find(nanSum ~= size(impMatrix_all,2));
impMatrix_trim = impMatrix_all(keepRow,:);
gene_indices = gene_indices(keepRow);

%clean up by creating a single structure of the imputed matrices for each media condition
impMatrix.LB = impMatrix_trim(:,1:6);
impMatrix.RPMIs = impMatrix_trim(:,7:12);
impMatrix.RPMI = impMatrix_trim(:,13:18);

%% Replace remaining Nan with min value from each column
minImp.LB = fillmissing(impMatrix.LB,'constant',min([data.LB_bySample.s1 , data.LB_bySample.s2]));
minImp.RPMIs = fillmissing(impMatrix.RPMIs,'constant',min([data.RPMIs_bySample.s1 , data.RPMIs_bySample.s2]));
minImp.RPMI = fillmissing(impMatrix.RPMI,'constant',min([data.RPMI_bySample.s1 , data.RPMI_bySample.s2]));
minImp_all = [minImp.LB,minImp.RPMIs,minImp.RPMI];


%% PCA and Clustering
titleLoop = [{'LB'} {'RPMI-serum'} {'RPMI'}];
scatterLegend = unique(data.sampleNames_txt,'stable');
rep = kron([1,2],ones(1,3));
heatmapLegend = scatterLegend(rep);
mediaLoop = {'LB','RPMIs','RPMI'};

%%%%%%%%%%%%%%%% PCA Cluster
[coeff_all,score_pca_all,latent_all,tsquared_all,explained_all,mu_all] = pca(minImp_all');
%// WT LB
group1_all = score_pca_all(1:3,:);
%// KO LB
group2_all = score_pca_all(4:6,:);
%// WT RPMIs
group3_all = score_pca_all(7:9,:);
%// KO RPMIs
group4_all = score_pca_all(10:12,:);
%// WT RPMI
group5_all = score_pca_all(13:15,:);
%// KO RPMI
group6_all = score_pca_all(16:18,:);

fig_all = figure ();
scatter(group1_all(:,1), group1_all(:,2),140,'red', 'filled','o');
hold on
scatter(group2_all(:,1), group2_all(:,2),140,'red', 'filled','s');
scatter(group3_all(:,1), group3_all(:,2),140,'green', 'filled','o');
scatter(group4_all(:,1), group4_all(:,2),140,'green', 'filled','s');
scatter(group5_all(:,1), group5_all(:,2),140,'blue', 'filled','o');
scatter(group6_all(:,1), group6_all(:,2),140,'blue', 'filled','s');
title('Clustering of KO vs WT in All Media Conditions');
xlabel(strcat('Principal Component 1 (',num2str(explained_all(1)),'%)'))
ylabel(strcat('Principal Component 2 (',num2str(explained_all(2)),'%)'))
legend([{'WT in LB'},{'KO in LB'},{'WT in RPMIs'},{'KO in RPMIs'},{'WT in RPMI'},{'KO in RPMI'}], 'Location', 'northwest');


%%%%%%%%%%%%%%%% Hierarchical clustering, standardize across rows/genes
leg_LB = strcat('LB::',heatmapLegend);
leg_RPMIs = strcat('RPMIs::',heatmapLegend);
leg_RPMI = strcat('RPMI::',heatmapLegend);
leg_all = [leg_LB, leg_RPMIs, leg_RPMI];

colormap('redblue')
cgo_all = clustergram(minImp_all,'Standardize','Row');
set(cgo_all,'ColumnLabels',leg_all,'ColumnLabelsRotate',45);
set(cgo_all,'RowLabels',gene_indices);
cgo_all.Colormap = redblue;
cgo_all.Dendrogram = [5 0];
%Use plot to determine cluster IDs
Clusters={1986,1996, 1995,1994,1993};

% Retrieve cluster information
for i =1:length(Clusters)
    clusterID = strcat('cluster_',(num2str(Clusters{i})));
    clustStruct.(clusterID) = clusterGroup(cgo_all,Clusters{i},'row','InfoOnly','true');
    figure(i)
    plot(mean(clustStruct.(clusterID).ExprValues),'LineWidth',2) %plot average of gene values
    title(strcat('Cluster ',num2str(i)))
    xlim([1 18])
    ylim([24,30])
    set(gca, 'XTick', [1:1:18], 'XTickLabel', clustStruct.(clusterID).ColumnNodeNames)
    xtickangle(45)
    yl = ylim;
    patch([1 3 3 1],[yl(1) yl(1) yl(2) yl(2)],[0.8 0.8 0.8])
    patch([4 6 6 4],[yl(1) yl(1) yl(2) yl(2)],[0.94 0.94 0.94])
    patch([7 9 9 7],[yl(1) yl(1) yl(2) yl(2)],[0.8 0.8 0.8])
    patch([10 12 12 10],[yl(1) yl(1) yl(2) yl(2)],[0.94 0.94 0.94])
    patch([13 15 15 13],[yl(1) yl(1) yl(2) yl(2)],[0.8 0.8 0.8])
    patch([16 18 18 16],[yl(1) yl(1) yl(2) yl(2)],[0.94 0.94 0.94])
    set(gca,'children',flipud(get(gca,'children')))
end

addTitle(cgo_all,'Hierarchical clustering of all samples');
plot(cgo_all)

%% Determine which genes belong to each cluster
clusterIDs = fields(clustStruct);
allIDs = table(data.Protein_IDs, data.geneNames, data.KEGG_singular, data.uniprotID_singular, ...
    'VariableNames',{'proteinID','geneName', 'KEGG', 'uniprot'});
clusterGenes = allIDs(str2double(clustStruct.(clusterIDs{2}).RowNodeNames),:);




