function [data] = addVirulence(data)
%Add virulence annotations to data structure

% Inputs
%   data -- data strucrure

% Output
%   data -- data structure with virulence annotations

%% Read in data
[~,txt,~] = xlsread('PA virulence.xlsx'); %document containing PA virulence annotations
virulenceGene = txt(2:end,3);
virulenceProduct = txt(2:end,4);
virulenceSource = txt(2:end,5);
virulenceCategory = txt(2:end,9);

% Gene vector can't contain blanks, fill in blank entries with 'NA'
genes2compare = data.geneNames_merge;
idx = find(cellfun(@isempty,genes2compare) == 1);
genes2compare(idx) = {'NA'};
genes2compare = lower(genes2compare);

% Add virulence annotations to the data structure
virulenceFactor = cell(length(data.Protein_IDs),1); %initialize dummy vector for virulance annotations
for i =1:length(virulenceGene)
    vir = virulenceGene{i};
    if ~isempty(vir)
        [ia,ib] = ismember(lower(vir),genes2compare);
        if ia == 1
        virulenceFactor(ib) = {'virulence factor'};
        end
    end
end
        
data.virulenceFactor = virulenceFactor;
end
        
    
    