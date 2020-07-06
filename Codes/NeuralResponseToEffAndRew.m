% This code calculate neural responses to effort and reward.
% A cell is considered eff- (rew-) selective if its neural firing 
% descriminates the level of effort (reward) significantly. (ANOVA p<0.05, effect-size>0.1379)
% 
% @Feb 2018_SH

clc; clear; close all;
MatlabRoot = '/Volumes/Seagate Backup Plus Drive/My CCBN PC/Saeedeh_MATLAB/Ramp_task_Drug/codes/Elife_Revisions/Repository_Github/';
addpath(genpath([MatlabRoot 'Lib']));
addpath(genpath(MatlabRoot));
load([MatlabRoot 'Data/Data_Info.mat'])
DatatRoot = [MatlabRoot 'Data/'];
%% Variables
Plim = 0.05; % criterion of p value
nBins = 36;  % #spatial bins across the whole maze
%% Main
s = 18; % chosen sessions: 18th session in Sessions structure
[Rew,Eff] = deal(cell(2,1)); %1st row for pre-injection and 2nd row for post-inj.
for D = 0:1  % D = 0 for pre-injection and D = 1 for post-inj.
    cd(DatatRoot)
    load(['Data_' Sessions(s).ID '_' num2str(Sessions(s).Session) '_D' num2str(D)])
    RampHeight = unique([trial_events.heights]);
    nCells = Sessions(s).nCells;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % selectivity: Finding selective neurons in each spatial bin
    [Rew{D+1,1},Eff{D+1,1}] = SlctiveNeurons(trial_events,RampHeight,[],nCells,nBins,'Plim',Plim);
end

