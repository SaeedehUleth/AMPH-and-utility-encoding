function [R,E] = Categorizing_Cell_FR(cellnum,trial_events,fld,varargin)
% Categorizing FR of a cell in reward and effort trials
% 
%            Left       |    Right
% ------------------------------------
%           LR(0inch)   |    LR(0inch)   -> 20 force trials
% ------------------------------------
% Rew:    R1:LR(0inch)  |   R2:HR(0inch)   -> 20 force trials
%         R1:HR(0inch)  |   R2:LR(0inch)   -> 20 force trials
% ------------------------------------
% Eff:    E1:HR(9inch)  |   E2:HR(0inch)   -> 20 force trials
%         E1:HR(0inch)  |   E2:HR(9inch)   -> 20 force trials
% 
% 
% INPUTS:   -cellnum: is the cell number(row number) 
%           in the matrix of trial_events().fld
%           -trial_events is a structure of trial events.
%           -fld(string) is the name of the filed which is to be
%           analyzed for each category. e.g. FirRatePerBin or SortCellFR_PerBin
%           -(Optional)Reward: LRew/HRew are the time duration in which 
%           low/high reward feeders are open; unit:'ms'
%           -(Optional)RampHeight: LHeight/HHeight are heights of the ramp 
%           for low/med/high effort; unit:'inch'
%           -(Optional)bad_trials: is a vector of trial numbers you don't 
%           want to be involved in analysis.
% 
% OUTPUT:   -R(Reward),E(Effort): are three cells, each of which
%            containing FR of the cell(determined by 'cellnum') in that categgory.
%           R = {R1:LRewFR}{R2:HRewFR}   
%               {R1:HRewFR}{R2:LRewFR}
%           E = {E1:HEffFR}{E2:LEffFR}
%               {E1:LEffFR}{E2:HEffFR}
%           Assuming left trials as R1(or E1) and right trials as R2(or E2)
% 
% Example: [S,R,E] = ctgcellFR(1,trial_events,'CellFR_PerBin','RampHeight',[0 9],'bad_trials',[1 2 3]);
% 
% 
% @Feb 2018 - Saeedeh Hashemnia

%% Inputs
Reward = [200;600];
RampHeight = [0;9];
bad_trials = [];
assignopts(who, varargin);
LRew = Reward(1); HRew = Reward(2);
LHeight = RampHeight(1); HHeight = RampHeight(2);
trial_events(bad_trials(:)) = [];
nTrials = length(trial_events);
R = cell(2,2);    % R = {R1:LRewFR}{R2:HRewFR}   
%                       {R1:HRewFR}{R2:LRewFR}
E = cell(2,2);    % E = {E1:HEffFR}{E2:LEffFR}
%                       {E1:LEffFR}{E2:HEffFR}

%% Reward category
R_trials = trialn_inEachGroup('R',trial_events,'RampHeight',RampHeight,'Reward',Reward,'bad_trials',bad_trials);
for ii = 1:2
    for jj = 1:2
        for tr = R_trials{ii,jj}
            R{ii,jj}(end+1,:) = trial_events(tr).(fld)(cellnum,:);
        end
    end
end

%% Effort category
E_trials = trialn_inEachGroup('E',trial_events,'RampHeight',RampHeight,'Reward',Reward,'bad_trials',bad_trials);
for ii = 1:2
    for jj = 1:2
        for tr = E_trials{ii,jj}
            E{ii,jj}(end+1,:) = trial_events(tr).(fld)(cellnum,:);
        end
    end
end
return