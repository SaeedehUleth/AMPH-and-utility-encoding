function [Pval,vec,effect_sz,mvec,groupCompare] = Categorizing_Cell(I,varargin)
% Categorizing a cell based on anova and t-test of its firing rate
% Categories: Reward, Effort
% 
% INPUTS:   - I: either R or E which are the outputs of the function "Categorizing_Cell_FR"
%            R or E is a MATLAB cell format containing FR matrix of the 
%            cell in each subcategory:
%            R = {LRewFR}{HRewFR}
%            E = {LEffFR}{MEffFR}{HEffFR}
% 
%            - Plim(Optional): limit of P value to be accepted as a significant diff.
%             Its defult value is 0.05.
% 
% OUTPUT:   - vec: A logical vector of size nBins showing in which bin Pval
%            is less than Plim
% 
%           - mvec: is same as vec but based on pairwise comparison. 
%            It is a binary vector of size nBins.
%            (or a matrix of size 3 x nBins if I=E).
%            1 means the cell is showing significant difference in
%            subcategories of I. If size of I is 2, then mvec = vec;
% 
%           - effect_sz: effect size for anova1 which is partial eta squred
%           obtaining by deviding SS(effect)/SS(total);
%           Remeber: SS(effect) = SS(between) and SS(total) = SS(effect)+SS(erroe)
%           eta_squared~ .0099 small, .0588 medium, and .1379 large effect size
%           In a one-way ANOVA, Eta Squared and Partial Eta Squared will be equal, 
%           but this isn’t true in models with more than one independent variable.
% 
%           - groupCompare: is a matrix of pairwise comparison results.
%           e.g. For example, suppose one row of it contains the following entries.
%           1.0000  2.0000  1.9442  8.2206  14.4971 0.0432
%           These numbers indicate that the mean of group 1 minus the mean 
%           of group 2 is estimated to be 8.2206, and a 95% confidence interval
%           for the true difference of the means is [1.9442, 14.4971]. 
%           The p-value for the corresponding hypothesis test that the difference
%           of the means of groups 1 and 2 is significantly different from zero is 0.0432.
%            
%  
% @Feb 2018 - SH

Plim = 0.05;
assignopts(who, varargin);

nBins = size(I{1},2);
groupCompare = cell(1,nBins);
[Pval,effect_sz] = deal(zeros(1,nBins));
mvec = zeros(1,nBins);
samplesz = [size(I{1},1) size(I{2},1)];
for i = 1:nBins
    sample_vec = [I{1}(:,i);I{2}(:,i)];
    membership_vec = [ones(samplesz(1),1);2*ones(samplesz(2),1)];
    [Pval(i),tbl,multicomp] = anova1(sample_vec,membership_vec,'off');
    effect_sz(i) = tbl{2,2}/tbl{4,2};
    groupCompare{i} = multcompare(multicomp,'CType','bonferroni','Display','off');
    if groupCompare{i}(1,6)<=Plim
        mvec(1,i) = 1;
    end
end
vec = Pval<Plim;
return