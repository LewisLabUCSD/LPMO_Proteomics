function [data] = addIDs(fileName,data)
% Add additional protein/gene IDs to data structure

% Inputs
%   data -- data strucrure
%   fileName -- file containing mapping between gene/protein IDs

% Output
%   data -- data structure with additinal gene/protein IDs (uniprot, gene names, protein names, and KEGG)

%% Read in data
[~,txt_uniprot,~] = xlsread(strcat(fileName,'.xlsx'),'refseq-to-uniprot'); %from uniprot ID mapping tool
[~,txt_geneName,~] = xlsread(strcat(fileName,'.xlsx'),'uniprot-to-geneName'); %from uniprot ID mapping tool
[~,txt_KEGG,~] = xlsread(strcat(fileName,'.xlsx'),'uniprot-to-KEGG'); %from uniprot ID mapping tool

%%  To convert from refSeq to uniprot
refSeq_id = txt_uniprot(2:end,1); %ref seq ids
uniprotID_refSeq = txt_uniprot(2:end,2); %uniprot ids
% create data field for uniprot IDs
ID = cell(length(data.Protein_IDs),1,1); %dummy for uniprot
IDsingular = ID; %dummy for uniprot singlular IDs

for i = 1:length(refSeq_id)
    proteinID = refSeq_id{i};
    idx = find(contains(data.Protein_IDs,proteinID)==1);
    for j=1:length(idx)
        if isempty(ID{idx(j)})
            ID{idx(j)}=uniprotID_refSeq{i};
            IDsingular{idx(j)} = uniprotID_refSeq{i};
        else
            ID{idx(j)}=strcat(ID{idx(j)},';',uniprotID_refSeq{i});
        end
    end
end
% Note: there is some wierd strings of ";", remove
data.uniprotID_singular = strrep(IDsingular,';;;;;','');
data.uniprotIDs = strrep(ID,';;;;;','');

%% To convert from uniprot to geneNames
uniprotID_geneName = txt_geneName(2:end,1);
geneNames = txt_geneName(2:end,2);
% create data field for geneNames 
ID2 = cell(length(data.Protein_IDs),1,1); %dummy for gene names

for i = 1:length(uniprotID_geneName)
    proteinID2 = uniprotID_geneName{i};
    idx2 = find(contains(data.uniprotIDs,proteinID2)==1);
    for j=1:length(idx2)
        if isempty(ID2{idx2(j)})
            ID2{idx2(j)}=geneNames{i};
        else
            ID2{idx2(j)}=strcat(ID2{idx2(j)},';',geneNames{i});
        end
    end
end
data.geneNames = ID2;

%% To convert to protein names
protName_refSeq = txt_uniprot(2:end,5);
% create data field for proteinNames 
IDprotName = cell(length(data.Protein_IDs),1,1);%dummy for protein names

for i = 1:length(refSeq_id)
    refSeq_prot = refSeq_id{i};
    idxProt = find(contains(data.Protein_IDs,refSeq_prot)==1);
    for j=1:length(idxProt)
        if isempty(IDprotName{idxProt(j)})
            IDprotName{idxProt(j)}=protName_refSeq{i};
        else
            IDprotName{idxProt(j)}=strcat(IDprotName{idxProt(j)},';',protName_refSeq{i});
        end
    end
end
% Note: there is some wierd strings of ";", remove
data.proteinNames = strrep(IDprotName,';;;;;',''); %proteome

%% To convert from uniprot to KEGG
uniprotID_KEGG = txt_KEGG(2:end,1);
KEGG = txt_KEGG(2:end,2);
% create data field for KEGG ids
ID3 = cell(length(data.Protein_IDs),1,1);%dummy for KEGG
ID3singular = cell(length(data.Protein_IDs),1,1);%dummy for KEGG

for i = 1:length(uniprotID_KEGG)
    proteinID3 = uniprotID_KEGG{i};
    idx3 = find(contains(data.uniprotIDs,proteinID3)==1);
    for j=1:length(idx3)
        if isempty(ID3{idx3(j)})
            ID3{idx3(j)}=KEGG{i};
            ID3singular{idx3(j)} = KEGG{i};
        else
            ID3{idx3(j)}=strcat(ID3{idx3(j)},';',KEGG{i});
        end
    end
end

data.KEGG = ID3;
data.KEGG_singular = ID3singular;
end