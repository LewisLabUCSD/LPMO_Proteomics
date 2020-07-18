function [sigList, otherList,foldChange]= volcanoPlots(data,q_thresh,Ttest)
%Calculate significantly differentially expressed proteins and generate volcano plots

% Inputs
%   data -- data structure
%   Ttest -- structure containing t/F test statistics 
%   q_thresh -- FDR threshold for significance

% Output
%   sigList -- data structure with significant gene information separated by neg vs pos FC
%   otherList -- data structure with significant gene information NOT separated by neg vs pos FC
%   foldChange -- structure containing sample FCs compared to control

%%
strain = {'KO','WT'}; % samples to compare against control (CTR)

for j = 1:2
    %----------- calculate fold change, filter out those with NA test stats
    FC = nanmean(data.(strcat('Spleen',(strain{j}))),2)-nanmean(data.SpleenCTR,2);
    plotIDX = ~isnan(Ttest.Spleen.(strain{j}).q);
    plotIDX_idx = find(plotIDX == 1);
    plotQ = Ttest.Spleen.(strain{j}).q(plotIDX);
    plotQ_log = -log10(plotQ);
    plotFC = FC(plotIDX);
    %----------- find those that meet significant thresholds
    sig_pos = find(plotQ_log >= -log10(q_thresh) & plotFC >= log2(1.5));
    sig_neg = find(plotQ_log >= -log10(q_thresh) & plotFC <= -log2(1.5));
    fullIDX_pos = plotIDX_idx(sig_pos);
    fullIDX_neg = plotIDX_idx(sig_neg);
    %----------- plot
    figure()
    scatter(plotFC,plotQ_log,25,'k');
    hold on
    xlim([-10 5])
    ylim([0 4])
    scatter(plotFC(sig_pos),plotQ_log(sig_pos),25,'r','filled');
    scatter(plotFC(sig_neg),plotQ_log(sig_neg),25,'r','filled');
    line(xlim,[-log10(q_thresh),-log10(q_thresh)],'Color','red','LineStyle','--');
    line([log2(1.5),log2(1.5)],ylim,'Color','red','LineStyle','--');
    line([-log2(1.5),-log2(1.5)],ylim,'Color','red','LineStyle','--');
    title({'Spleen' ; strcat(strain{j},' vs CTR')});
    xlabel('log2 Fold Change')
    ylabel('-log10 FDR')
    hold off
    %----------- sig list (for enrichment)
    sigList.(strcat('Spleen',(strain{j}))).logFC_pos = plotFC(sig_pos);
    sigList.(strcat('Spleen',(strain{j}))).logFC_neg = plotFC(sig_neg);
    sigList.(strcat('Spleen',(strain{j}))).q_pos = plotQ(sig_pos);
    sigList.(strcat('Spleen',(strain{j}))).q_neg = plotQ(sig_neg);
    sigList.(strcat('Spleen',(strain{j}))).ProteinIDs_pos = data.Protein_IDs(fullIDX_pos);
    sigList.(strcat('Spleen',(strain{j}))).ProteinIDs_neg = data.Protein_IDs(fullIDX_neg);
    sigList.(strcat('Spleen',(strain{j}))).GeneNames_pos = data.Gene_names(fullIDX_pos);
    sigList.(strcat('Spleen',(strain{j}))).GeneNames_neg = data.Gene_names(fullIDX_neg);
    sigList.(strcat('Spleen',(strain{j}))).ProteinNames_pos = data.Protein_names(fullIDX_pos);
    sigList.(strcat('Spleen',(strain{j}))).ProteinNames_neg = data.Protein_names(fullIDX_neg);
    %----------- other list (for print out)
    otherList.(strcat('Spleen',(strain{j}))).logFC = [plotFC(sig_pos);plotFC(sig_neg)];
    otherList.(strcat('Spleen',(strain{j}))).q = [plotQ(sig_pos);plotQ(sig_neg)];
    otherList.(strcat('Spleen',(strain{j}))).ProteinIDs = [data.Protein_IDs(fullIDX_pos); data.Protein_IDs(fullIDX_neg)];
    otherList.(strcat('Spleen',(strain{j}))).GeneNames = [data.Gene_names(fullIDX_pos); data.Gene_names(fullIDX_neg)];
    otherList.(strcat('Spleen',(strain{j}))).ProteinNames = [data.Protein_names(fullIDX_pos); data.Protein_names(fullIDX_neg)];
    %----------- save FC in structure
    foldChange.(strcat('Spleen',(strain{j}))) = FC;
end
        
    