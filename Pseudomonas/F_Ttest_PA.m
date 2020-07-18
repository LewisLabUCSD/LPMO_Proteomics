function [Ttest] = F_Ttest_PA(data)
%Create data structure with F and T test statistics

% Inputs
%   data -- data structure

% Output
%   Ttest -- structure containing 1)F test statistic, 2)p-value from t-test, 3) FDR adjusted p-value

%% LB
WT_LB=data.LB(:,1:3); %WT samples
KO_LB=data.LB(:,4:6); %KO samples
h = vartest2(WT_LB,KO_LB, 'Dim',2); % Two sample F test for equal variance
pArray = NaN(size(h)); % initilize dummy vector of sig values 

% Run ttest with equal and unequal variance based on result from F test
for k = 1:length(h)
        F = h(k);
        if F == 1
            [~,p_unequal] = ttest2(WT_LB(k,:),KO_LB(k,:),'Dim',2,'Vartype','unequal');
            pArray(k) = p_unequal;
        elseif F == 0
            [~,p_equal] = ttest2(WT_LB(k,:),KO_LB(k,:),'Dim',2);
            pArray(k) = p_equal;
        else
            pArray(k) = pArray(k);
        end
end
Ttest.LB.Isogenic_mutant.Ftest = h;
Ttest.LB.Isogenic_mutant.pVal = pArray; 

% BH adjustment
[FDR] = mafdr(pArray,'BHFDR',true);
Ttest.LB.Isogenic_mutant.q = FDR;

%% RPMIs
WT_RPMIs=data.RPMIs(:,1:3);
KO_RPMIs=data.RPMIs(:,4:6); %KO samples
h = vartest2(WT_RPMIs,KO_RPMIs, 'Dim',2); % Two sample F test for equal variance
pArray = NaN(size(h)); % initilize dummy vector of sig values 

% Run ttest with equal and unequal variance based on result from F test
for k = 1:length(h)
        F = h(k);
        if F == 1
            [~,p_unequal] = ttest2(WT_RPMIs(k,:),KO_RPMIs(k,:),'Dim',2,'Vartype','unequal');
            pArray(k) = p_unequal;
        elseif F == 0
            [~,p_equal] = ttest2(WT_RPMIs(k,:),KO_RPMIs(k,:),'Dim',2);
            pArray(k) = p_equal;
        else
            pArray(k) = pArray(k);
        end
end
Ttest.RPMIs.Isogenic_mutant.Ftest = h;
Ttest.RPMIs.Isogenic_mutant.pVal = pArray; 

% BH adjustment
[FDR] = mafdr(pArray,'BHFDR',true);
Ttest.RPMIs.Isogenic_mutant.q = FDR;

%% RPMI
WT_RPMI=data.RPMI(:,1:3);
KO_RPMI=data.RPMI(:,4:6); %KO samples
h = vartest2(WT_RPMI,KO_RPMI, 'Dim',2); % Two sample F test for equal variance
pArray = NaN(size(h)); % initilize dummy vector of sig values 

% Run ttest with equal and unequal variance based on result from F test
for k = 1:length(h)
        F = h(k);
        if F == 1
            [~,p_unequal] = ttest2(WT_RPMI(k,:),KO_RPMI(k,:),'Dim',2,'Vartype','unequal');
            pArray(k) = p_unequal;
        elseif F == 0
            [~,p_equal] = ttest2(WT_RPMI(k,:),KO_RPMI(k,:),'Dim',2);
            pArray(k) = p_equal;
        else
            pArray(k) = pArray(k);
        end
end
Ttest.RPMI.Isogenic_mutant.Ftest = h;
Ttest.RPMI.Isogenic_mutant.pVal = pArray; 

% BH adjustment
[FDR] = mafdr(pArray,'BHFDR',true);
Ttest.RPMI.Isogenic_mutant.q = FDR;

end