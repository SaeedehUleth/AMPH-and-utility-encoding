function [Rew,Eff] = SlctiveNeurons(trial_events,RampHeight,bad_trials,nCells,nBins,varargin)
% categorizing units to Reward- and Effortiselective ones in each bin
% ANOVA(p<0.05)
% 
% @Feb 2018-SH
Plim = 0.05;
assignopts(who, varargin);

[vec_R2,vec_R1,vec_E1,vec_E2] = deal(zeros(nCells,nBins));
[Pval_R1,Pval_R2,Pval_E1,Pval_E2] = deal(zeros(nCells,nBins));
[efct_sz_R1,efct_sz_R2,efct_sz_E1,efct_sz_E2] = deal(zeros(nCells,nBins));
[avg_FR_R1,avg_FR_R2,avg_FR_E1,avg_FR_E2] = deal(zeros(nCells,2*nBins));

for i = 1:nCells
    [R,E] = Categorizing_Cell_FR(i,trial_events,'CellFR_PerBin','RampHeight',RampHeight,'bad_trials',bad_trials);
    
    avg_FR_R1(i,:) = [mean(R{1,1},1),mean(R{2,1},1)];
    avg_FR_R2(i,:) = [mean(R{1,2},1),mean(R{2,2},1)];
    
    avg_FR_E1(i,:) = [mean(E{1,1},1),mean(E{2,1},1)];
    avg_FR_E2(i,:) = [mean(E{1,2},1),mean(E{2,2},1)];
    
    [Pval_R1(i,:),vec_R1(i,:),efct_sz_R1(i,:),~,~] = Categorizing_Cell(R(1:2,1),'Plim',Plim);
    [Pval_R2(i,:),vec_R2(i,:),efct_sz_R2(i,:),~,~] = Categorizing_Cell(R(1:2,2),'Plim',Plim);
    [Pval_E1(i,:),vec_E1(i,:),efct_sz_E1(i,:),~,~] = Categorizing_Cell(E(1:2,1),'Plim',Plim);
    [Pval_E2(i,:),vec_E2(i,:),efct_sz_E2(i,:),~,~] = Categorizing_Cell(E(1:2,2),'Plim',Plim);
end
vec_R = vec_R1 + vec_R2;
vec_R(vec_R> 1) = 1;
vec_E = vec_E1 + vec_E2;
vec_E(vec_E> 1) = 1;

Rew = struct('vec1',vec_R1,'vec2',vec_R2,'vec',vec_R,'Pval1',Pval_R1,'Pval2',Pval_R2,...
    'avg_FR1',avg_FR_R1,'avg_FR2',avg_FR_R2,'efct_sz1',efct_sz_R1,'efct_sz2',efct_sz_R2);
Eff = struct('vec1',vec_E1,'vec2',vec_E2,'vec',vec_E,'Pval1',Pval_E1,'Pval2',Pval_E2,...
    'avg_FR1',avg_FR_E1,'avg_FR2',avg_FR_E2,'efct_sz1',efct_sz_E1,'efct_sz2',efct_sz_E2);
return
