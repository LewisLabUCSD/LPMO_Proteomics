function [data] = DataStructure(fileName,sheetName, sampleNumber, Other,norm)
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

% Output
%   data -- data structure containing data and metadata

%% Read in data
[num,txt,~] = xlsread(strcat(fileName,'.xlsx'),sheetName); %read in data
data.sampleNames_num = unique(num(1,1:(3*sampleNumber)),'stable'); %each sample is associated with a numerical value
data.sampleNames_txt = txt(1,1:(3*sampleNumber)); %sample names
data.allData = num(2:end,1:3*sampleNumber); %omics data
samples_per_media = sampleNumber/3; %number of samples per media type

%% Normalization
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan')); %ignore nan


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

%% Group samples 
%separate data based on media condition
data.LB = data.allData(:,1:(samples_per_media*3));
data.RPMIs = data.allData(:,((samples_per_media*3)+1):((samples_per_media*2)*3));
data.RPMI = data.allData(:,(((samples_per_media*2)*3)+1):((samples_per_media*3)*3));

data.groups = num(1,1:(3*sampleNumber)); %vector of numerical sample names, 3 replicates per sample 

% separated based on sample (3 replicates per sample)
for i=1:length(data.sampleNames_num)
    data.(strcat('sample_',num2str(data.sampleNames_num(i))))= data.allData(:,(((i-1)*3)+1):(i*3));
end

% add additional information as indicated by "Other"
for j=1:length(Other)
    [r,c]=find(strcmp(txt,Other{j}));
    name = strrep(Other{j},' ','_');
    data.(name) = txt(3:end,c);
end

end

