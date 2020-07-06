% Behavioural analysis for AMPH- and AMPH+ comparison on Ramp-task
% Comparison of median time spent on reward zones (based on t_CF and
% t_SF)
% 
% @Jan 2018_SH

clc; clear; close all;
MatlabRoot = '/Volumes/Seagate Backup Plus Drive/My CCBN PC/Saeedeh_MATLAB/Ramp_task_Drug/codes/Elife_Revisions/Repository_Github/';
addpath(genpath([MatlabRoot 'Lib']));
addpath(genpath(MatlabRoot));
load([MatlabRoot 'Data/Data_Info.mat'])
DatatRoot = [MatlabRoot 'Data/'];
%% Inputs
% There are 22 sessions in total, we here run the analysis for 1 session
% as an example.
[dur_CF,dur_SF_LR,dur_SF_HR] = deal(cell(2,22));
[Pre_TORZ,Post_TORZ] = deal(cell(1,22));
[Pre_m_TORZ,Pre_SEM_TORZ,Post_m_TORZ,Post_SEM_TORZ] = deal(zeros(2,22));
%% Main
s = 18; % chosen sessions: 18th session in Sessions structure
a = load([DatatRoot 'Data_' Sessions(s).ID '_' num2str(Sessions(s).Session) '_D0']);
Pre_trial_events = a.trial_events;

b = load([DatatRoot 'Data_' Sessions(s).ID '_' num2str(Sessions(s).Session) '_D1']);
Post_trial_events = b.trial_events;

% Time occupancy on Reward Zones
[Pre_TORZ{s},dur_CF{1,s},dur_SF_LR{1,s},dur_SF_HR{1,s},Pre_m_TORZ(:,s),Pre_SEM_TORZ(:,s)] = TORZ(Pre_trial_events);
[Post_TORZ{s},dur_CF{2,s},dur_SF_LR{2,s},dur_SF_HR{2,s},Post_m_TORZ(:,s),Post_SEM_TORZ(:,s)] = TORZ(Post_trial_events);