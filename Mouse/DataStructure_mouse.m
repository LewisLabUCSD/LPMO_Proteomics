function [data] = DataStructure_mouse(fileName, totSamples, Other,norm,norm2)
%Create data structure with all data and metadata for downstream analysis

% Inputs
%   fileName -- name of excel file containing your omics data
%   sheetName -- name of sheet in excel file containing your data
%   sampleNumber -- number of samples (do not count replicates as individual samples)
%   Other -- other column information you would like to include 
%   norm -- type of normalization to use:
    % z-genes: calculates z scores across rows/genes
    % z_samples: calculates z scores across columns/samples
    % quantile: quantile normalization
%   norm2 -- option to quantile normalize z-scores    

% Output
%   data -- data structure containing data and metadata

%% Read in data
[num,txt,~] = xlsread(strcat(fileName,'.xlsx')); %read in data
data.sampleNames_num = num(1,1:(totSamples)); %each sample is associated with a numerical value
data.sampleNames_txt = txt(2,1:(totSamples)); %sample names
data.allData = num(3:end,1:totSamples); %omics data

%% Normalization
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan')); %ignore nan

figure(1)
boxplot(data.allData)
title('Pre-norm');
ylabel('proteins (LFQ)')
xticklabels(data.sampleNames_txt);
set(gca,'XTickLabelRotation',45)
ylabel('proteins (LFQ)')


if strcmp(norm,'z_genes')
    %normalize across rows/genes
    zScores_genes = zscor_xnan(data.allData'); 
    data.allData = zScores_genes';
elseif strcmp(norm,'z_samples')
    %normalize across columns/samples
    data.allData = zscor_xnan(data.allData);
elseif strcmp(norm,'quantile')
    data.allData = quantilenorm(data.allData);
end

% option to quantile normalize z-scores
if strcmp(norm2,'quantile')
    data.allData = quantilenorm(data.allData);
end

figure(2)
boxplot(data.allData)
title('Post-norm')
xticklabels(data.sampleNames_txt);
set(gca,'XTickLabelRotation',45)
ylabel('protein z-scores')

%% Group samples
% separate based on sample: WT, KO, CTR
uniqueSampleNames = unique(data.sampleNames_txt,'stable');
for i = 1:length(uniqueSampleNames)
    ind = find(strcmp(data.sampleNames_txt,uniqueSampleNames{i}));
    data.(uniqueSampleNames{i}) = data.allData(:,ind);
end

% add additional information as indicated by "Other"
for j=1:length(Other)
    [r,c]=find(strcmp(txt,Other{j}));
    name = strrep(Other{j},' ','_');
    data.(name) = txt(3:end,c);
end

end

