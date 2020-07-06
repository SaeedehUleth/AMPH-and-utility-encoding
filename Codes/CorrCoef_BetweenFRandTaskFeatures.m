% This code calculate correlation coefficients between neural responses 
% and task features (effort and reward).
% 
% @Feb 2018_SH

clc; clear; close all;
MatlabRoot = '/Volumes/Seagate Backup Plus Drive/My CCBN PC/Saeedeh_MATLAB/Ramp_task_Drug/codes/Elife_Revisions/Repository_Github/';
addpath(genpath([MatlabRoot 'Lib']));
addpath(genpath(MatlabRoot));
load([MatlabRoot 'Data/Data_Info.mat'])
DatatRoot = [MatlabRoot 'Data/'];
%% Main
% There are 22 sessions in total, we here run the analysis for 1 session
% as an example.
s = 18; % chosen sessions: 18th session in Sessions structure
% D = 0 for pre-injection and D = 1 for post-inj.
[rou_E,rou_R,beta_E,beta_R] = deal(cell(4,22));  % 4 columns as {D0(left);D0(right);D1(left);D1(Right)}
for D = 0:1
    cd(DatatRoot)
    load(['Data_' Sessions(s).ID '_' num2str(Sessions(s).Session) '_D' num2str(D)],'trial_events');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % categorizing firing rate
    nCell = Sessions(s).nCells;     % # recorded cells
    nTrials = length(trial_events); % # trials
    
    % R1/E1:left and R2/E2:Right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = trialn_inEachGroup('R',trial_events);
    % if type is R: T{1:2,1} = R1:{LR;HR} and T{1:2,2} = R2:{HR;LR}
    % LR: left side, low hight     HR: right side, low hight
    % HR: left side, low hight     LR: right side, low hight
    LR1 = [T{1,1}]; HR1 = [T{2,1}]; LR2 = [T{1,2}]; HR2 = [T{2,2}];
    T = trialn_inEachGroup('E',trial_events);
    % if type is E: T{1:2,1} = E1:{HE;LE} and T{1:2,2} = E2:{LE;HE}
    % HE: left side, high rew     LE: right side, high rew
    % LE: left side, high rew     HE: right side, high rew
    LE1 = [T{2,1}]; HE1 = [T{1,1}]; LE2 = [T{1,2}]; HE2 = [T{2,2}];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preparing x values for the regresion
    x1_Eff = [ones(length(LE1),1);2*ones(length(HE1),1)];
    x2_Eff = [ones(length(LE2),1);2*ones(length(HE2),1)];
    x1_Rew = [ones(length(LR1),1);2*ones(length(HR1),1)];
    x2_Rew = [ones(length(LR2),1);2*ones(length(HR2),1)];
    
    %% preparing y values for the regresion
    [FR_trial_event,LE1_FR,LE2_FR,HE1_FR,HE2_FR,LR1_FR,LR2_FR,HR1_FR,HR2_FR] = deal(cell(nCell,1));
    for celln = 1:nCell
        for i = 1:nTrials
            FR_trial_event{celln}(i,:) = trial_events(i).CellFR_PerBin(celln,:);
        end
        % for some sessions some of the followings are empty
        LE1_FR{celln} = FR_trial_event{celln}(LE1,:);
        LE2_FR{celln} = FR_trial_event{celln}(LE2,:);
        HE1_FR{celln} = FR_trial_event{celln}(HE1,:);
        HE2_FR{celln} = FR_trial_event{celln}(HE2,:);
        LR1_FR{celln} = FR_trial_event{celln}(LR1,:);
        LR2_FR{celln} = FR_trial_event{celln}(LR2,:);
        HR1_FR{celln} = FR_trial_event{celln}(HR1,:);
        HR2_FR{celln} = FR_trial_event{celln}(HR2,:);
    end
    %% Calculating Rou and Beta values
    % Preparing x-values for the regression
    x1_ux_Eff = x1_Eff - mean(x1_Eff); x1_ux_Rew = x1_Rew - mean(x1_Rew);
    x2_ux_Eff = x2_Eff - mean(x2_Eff); x2_ux_Rew = x2_Rew - mean(x2_Rew);
    
    % Pre-allocations
    [rou_E{(2*D)+1,s},rou_R{(2*D)+1,s}] = deal(zeros(nCell,36));
    [beta_E{(2*D)+1,s},beta_R{(2*D)+1,s}] = deal(zeros(nCell,36));
    [rou_E{(2*D)+2,s},rou_R{(2*D)+2,s}] = deal(zeros(nCell,36));
    [beta_E{(2*D)+2,s},beta_R{(2*D)+2,s}] = deal(zeros(nCell,36));
    
    for BIN = 1:36
        for celln = 1:nCell
            % Preparing y-values for the regression
            y1_Eff = [LE1_FR{celln}(:,BIN);HE1_FR{celln}(:,BIN)];
            y2_Eff = [LE2_FR{celln}(:,BIN);HE2_FR{celln}(:,BIN)];
            y1_uy_Eff = y1_Eff - mean(y1_Eff); y2_uy_Eff = y2_Eff - mean(y2_Eff);
            
            y1_Rew = [LR1_FR{celln}(:,BIN);HR1_FR{celln}(:,BIN)];
            y2_Rew = [LR2_FR{celln}(:,BIN);HR2_FR{celln}(:,BIN)];
            y1_uy_Rew = y1_Rew - mean(y1_Rew); y2_uy_Rew = y2_Rew - mean(y2_Rew);
            
            if ~isempty(x1_Eff)
                % calculating rou values based on the formulla
                rou_E{(2*D)+1,s}(celln,BIN) = sum(x1_ux_Eff.*y1_uy_Eff)/(sqrt(sum(x1_ux_Eff.^2))*sqrt(sum(y1_uy_Eff.^2))); %refer to wikipedia(Pearson correlation coefficient) if any question,rou = cov(X,Y)/sigma(X).sigma(Y)
                rou_R{(2*D)+1,s}(celln,BIN) = sum(x1_ux_Rew.*y1_uy_Rew)/(sqrt(sum(x1_ux_Rew.^2))*sqrt(sum(y1_uy_Rew.^2)));  %rou = cov(X,Y)/sigma(X).sigma(Y)
                % fitting a linear model and calculating beta values
                mdl_Eff = fitlm(x1_Eff,y1_Eff,'linear');  % in which R-squared equals to rou^2
                beta_E{(2*D)+1,s}(celln,BIN) = mdl_Eff.Coefficients.Estimate(2); % y ~ beta*x+alpha where alpha is mdl.Coefficients.Estimate(1);
                mdl_Rew = fitlm(x1_Rew,y1_Rew,'linear');  % in which R-squared equals to rou^2
                beta_R{(2*D)+1,s}(celln,BIN) = mdl_Rew.Coefficients.Estimate(2); % y ~ beta*x+alpha where alpha is mdl.Coefficients.Estimate(1);
            end
            if ~isempty(x2_Eff)
                % calculating rou values based on the formulla
                rou_E{(2*D)+2,s}(celln,BIN) = sum(x2_ux_Eff.*y2_uy_Eff)/(sqrt(sum(x2_ux_Eff.^2))*sqrt(sum(y2_uy_Eff.^2))); %refer to wikipedia(Pearson correlation coefficient) if any question,rou = cov(X,Y)/sigma(X).sigma(Y)
                rou_R{(2*D)+2,s}(celln,BIN) = sum(x2_ux_Rew.*y2_uy_Rew)/(sqrt(sum(x2_ux_Rew.^2))*sqrt(sum(y2_uy_Rew.^2)));  %rou = cov(X,Y)/sigma(X).sigma(Y)
                % fitting a linear model and calculating beta values
                mdl_Eff = fitlm(x2_Eff,y2_Eff,'linear');  % in which R-squared equals to rou^2
                beta_E{(2*D)+2,s}(celln,BIN) = mdl_Eff.Coefficients.Estimate(2); % y ~ beta*x+alpha where alpha is mdl.Coefficients.Estimate(1);
                mdl_Rew = fitlm(x2_Rew,y2_Rew,'linear');  % in which R-squared equals to rou^2
                beta_R{(2*D)+2,s}(celln,BIN) = mdl_Rew.Coefficients.Estimate(2); % y ~ beta*x+alpha where alpha is mdl.Coefficients.Estimate(1);
            end
        end
    end
end
