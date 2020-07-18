function [sigList, otherList,meanMatrix] = newVolcanoPlots_PA(data,titles, structFields,Ttest,human_idx,q_thresh)
%Calculate significantly differentially expressed proteins and generate volcano plots

% Inputs
%   data -- data structure
%   titles -- condition specific titles
%   structFields -- number of samples (excluding replicates)
%   Ttest -- structure containing t/F test statistics 
%   human_idx -- Indices of human genes to filter out
%   q_thresh -- FDR threshold for significance

% Output
%   sigList -- data structure with significant gene information separated by neg vs pos FC
%   otherList -- data structure with significant gene information NOT separated by neg vs pos FC
%   meanMatrix -- matrix with mean sample values

%%
media = {'LB', 'RPMIs', 'RPMI'};
strains = {'Isogenic_mutant'};

meanMatrix = []; %initialize matrix
for j = 1:3

    % Set media-specific title plot and structure fields
    media_title = titles{j};
    media_field = structFields{j};
    % Extact media specific data
    WT=data.(media{j})(:,1:3); %WT data
    KO = data.(media{j})(:,4:6); %KO data
    % Calculate mean and fold change
    meanMatrix = [meanMatrix, nanmean(WT,2)]; %WT mean
    meanMatrix = [meanMatrix, nanmean(KO,2)]; %KO mean
    FC = nanmean(KO,2)-nanmean(WT,2); %Fold change
    % Apply significant threshold: FC>= 1.5 & q <= 0.05
    qvals = Ttest.(media{j}).Isogenic_mutant.q;
    sigPos = find(qvals<=0.05 & FC >= log2(1.5)); %indices of significant with pos FC compared to WT
    sigPos = setdiff(sigPos,human_idx); %remove human proteins
    sigNeg = find(qvals<=0.05 & FC <= -log2(1.5)); %indices of significant with neg FC compared to WT
    sigNeg = setdiff(sigNeg,human_idx); %remove human proteins

    % Volcano plot
    fig = figure();
    scatter(FC,-log10(qvals),25,'k');
    hold on
    xlim([-8 4])
    ylim([0 6])
    scatter(FC(sigPos),-log10(qvals(sigPos)),25,'r','filled');
    scatter(FC(sigNeg),-log10(qvals(sigNeg)),25,'r','filled');
    line(xlim,[-log10(q_thresh),-log10(q_thresh)],'Color','red','LineStyle','--');
    line([log2(1.5),log2(1.5)],ylim,'Color','red','LineStyle','--');
    line([-log2(1.5),-log2(1.5)],ylim,'Color','red','LineStyle','--');
    title(media_title);
    xlabel('log2 Fold Change')
    ylabel('-log10 FDR')
    hold off
        
    % Save lists of significant genes in structures for future analysis
    sigList.(media_field).logFC_pos = FC(sigPos);
    sigList.(media_field).logFC_neg = FC(sigNeg);
    sigList.(media_field).q_pos = qvals(sigPos);
    sigList.(media_field).q_neg = qvals(sigNeg);
    sigList.(media_field).ProteinIDs_pos = data.Protein_IDs(sigPos);
    sigList.(media_field).ProteinIDs_neg = data.Protein_IDs(sigNeg);
    sigList.(media_field).GeneNames_pos = data.geneNames_merge(sigPos);
    sigList.(media_field).GeneNames_neg = data.geneNames_merge(sigNeg);
    sigList.(media_field).ProteinNames_pos = data.proteinNames_merge(sigPos);
    sigList.(media_field).ProteinNames_neg = data.proteinNames_merge(sigNeg);
    sigList.(media_field).KEGG_pos = data.KEGG_singular(sigPos);
    sigList.(media_field).KEGG_neg = data.KEGG_singular(sigNeg);

    otherList.(media_field).logFC = [FC(sigPos);FC(sigNeg)];
    otherList.(media_field).q = [qvals(sigPos);qvals(sigNeg)];
    otherList.(media_field).ProteinIDs = [data.Protein_IDs(sigPos);data.Protein_IDs(sigNeg)];
    otherList.(media_field).GeneNames = [data.geneNames_merge(sigPos);data.geneNames_merge(sigNeg)];
    otherList.(media_field).ProteinNames = [data.proteinNames_merge(sigPos);data.proteinNames_merge(sigNeg)];
    otherList.(media_field).virulenceFactor = [data.virulenceFactor(sigPos);data.virulenceFactor(sigNeg)];
    otherList.(media_field).Fasta_headers = [data.Fasta_headers(sigPos);data.Fasta_headers(sigNeg)];
    
end
end