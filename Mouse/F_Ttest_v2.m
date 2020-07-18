function [Ttest] = F_Ttest_v2(data)
%Create data structure with F and T test statistics of 1) WT vs CTR, 2) KO vs CTR

% Inputs
%   data -- data structure

% Output
%   Ttest -- structure containing 1)F test statistic, 2)p-value from t-test, 3) FDR adjusted p-value

%% 
strain = {'KO','WT'}; % compate WT and KO vs control (CTR)

for j = 1:length(strain)
    h = vartest2(data.(strcat('Spleen',strain{j})),data.SpleenCTR, 'Dim',2); % Two sample F test for equal variance
    pArray = NaN(size(h)); % initilize dummy vector of sig values
    % Run ttest with equal and unequal variance based on result from F test
    for k = 1:length(h)
        F = h(k);
        if F == 1
            [~,p_unequal] = ttest2(data.(strcat('Spleen',strain{j}))(k,:),data.SpleenCTR(k,:),'Dim',2,'Vartype','unequal');
            pArray(k) = p_unequal;
        elseif F == 0
            [~,p_equal] = ttest2(data.(strcat('Spleen',strain{j}))(k,:),data.SpleenCTR(k,:),'Dim',2);
            pArray(k) = p_equal;
        else
            pArray(k) = pArray(k);
        end
    Ttest.Spleen.(strain{j}).Ftest = h;
    Ttest.Spleen.(strain{j}).pVal = pArray; 
    % BH adjustment
    [FDR] = mafdr(pArray,'BHFDR',true);
    Ttest.Spleen.(strain{j}).q = FDR;
    end
end
    
end