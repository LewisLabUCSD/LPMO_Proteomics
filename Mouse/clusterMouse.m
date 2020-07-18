load('proteomeAnalysis')

%% Impute
sampleField = {'KO','WT','CTR'};
ImpMatrix.All = [];
for j = 1:length(sampleField)
    tissueSamples = data.(strcat('Spleen',sampleField{j}));
    tissueSamples = knnimpute(tissueSamples);
    ImpMatrix.(strcat('Spleen',sampleField{j})) = tissueSamples;
    ImpMatrix.All = [ImpMatrix.All tissueSamples];
end

%% Check if some genes are NaN across all samples, remove them from the analysis
nanSum = sum(isnan(ImpMatrix.All),2);
keepRow = find(nanSum ~= size(ImpMatrix.All,2)); 
proteins_spleen = data.Protein_IDs(keepRow);
spleenImp_trim = ImpMatrix.All(keepRow,:);

%% Replace remaining Nan with min from that column
spleenImp_trim_min = fillmissing(spleenImp_trim,'constant',min(spleenImp_trim));

%% Hierarchical clustering
spleen_lab = data.sampleNames_txt(contains(data.sampleNames_txt,'Spleen')); %define labels for plots
colormap('redblue')   
cgo_spleen = clustergram(spleenImp_trim_min,'Standardize','Row');
set(cgo_spleen,'ColumnLabels',spleen_lab,'ColumnLabelsRotate',45);
cgo_spleen.Colormap = redblue;
addTitle(cgo_spleen,'Hierarchical clustering spleen samples');

%% PCA
[~,score_spleen,~,~,explained_spleen,~] = pca(spleenImp_trim_min');
%// KO
group1 = score_spleen(1:7,:);
%// WT
group2 = score_spleen(8:15,:);
%// CTR
group3 = score_spleen(16:19,:);
%// Plot as separate colours
fig = figure();
plot(group1(:,1), group1(:,2), 'b.', group2(:,1), group2(:,2), 'r.',...
    group3(:,1), group3(:,2),'g.','MarkerSize',50); 
xlabel(strcat('Principal Component 1 (',num2str(explained_spleen(1)),'%)'))
ylabel(strcat('Principal Component 2 (',num2str(explained_spleen(2)),'%)'))
title('PCA clustering of Spleen samples')
legend({'KO','WT','CTR'}, 'Location', 'northeast');

%% Polar dendrogram
Z= linkage(spleenImp_trim_min');
% Z= linkage(minImp.All(:,1:19)');

figure()
[h,T,perm] = dendrogram(Z,0,'Orientation','left','colorthreshold','default','Labels',data.sampleNames_txt(:,1:19));

figure()
[h2,T2,perm2] = polardendrogram_HM(Z,data.sampleNames_txt(1:19),0,'colorthreshold','default');
zoom(0.8);
view(2);
