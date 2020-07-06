% This code is producing a cache file for each session
% 
% Cache contains: -'spikes': All cell spikes of the esession
%                 -'trial_events': trial_event of all trial in which the first
%                   trial and Sessions(s).rmTrials have been removed.
%                 -'oVideo': The original Video structure before
%                   interpolation. Position are mapped to the reference
%                   (P1958_26) and timestamps are multiplied by 1e-6, so the
%                   unit is sec. oViideo contains only two fields:
%                   timestamps and corresponding position. rest of analysis
%                   has been applied to the 'Video' not 'oVideo'.
%                 -'Video': is the interpolated Video structure.
% 
% @Jan 2018-SH

clc; clear; close all;
MatlabRoot = 'C:\Users\saeedeh\Documents\MATLAB\';
addpath(genpath([MatlabRoot 'lib']));
addpath('C:\Users\saeedeh\Documents\MATLAB\Neuralynx');
addpath(genpath(MatlabRoot));
load([MatlabRoot 'Data\Data_Info.mat'])

D = 0;   % D must be 0 for PreDrug and 1 for PstDrug
%% Inputs
for s = 20
    if D == 0
        rmTrials = Sessions(s).PreD_rmTrials;
        t_end = Sessions(s).PreD_t_end;
        b = Sessions(s).PreD_LabviewLogFile;
    else
        rmTrials = Sessions(s).PstD_rmTrials;
        t_end = Sessions(s).PstD_t_end;
        b = Sessions(s).PstD_LabviewLogFile;
    end
    data_root = ['D:\Data\OrgClustered\Amph\phase2\' Sessions(s).ID '_' num2str(Sessions(s).Session) 'p\'];

    %% Loading Spikes and time series
    tFiles = dir([data_root 'tfiles\TT*.t']);
    nCells = length(tFiles);spikes = cell(1,nCells);
    for ii = 1:nCells
        spikes(ii) = load_spikes({[data_root 'tfiles\' tFiles(ii).name]});
    end
%     spikes is a cell containing timeseries of spikes of each neuron. the
%     unit is sec.

    %% fixing the timestamp offset between cheetah and NI
    % Assuming NI timestamp as the reference
    [epoch_times, trial_events ] = read_events_ramp([data_root 'Events.txt'], [data_root b]);
    % trial_events time stamp is based on cheetah
    % e.g. timestamp offset:  cheetah event file ts - NI ts = -6365651.8 microsec
    %% *************************************************************************
    %%*************************************************************************
    % Fixing problems in trial_events(not for P1958) if there is any:
    choice = zeros(1,length(trial_events));
    for i = 1:length(trial_events)
        choice(i) = find(trial_events(i).gates == 1);
        trial_events(i).choice = choice(i)-1;
        trial_events(i).feeder_ts(3) = trial_events(i).feeder_off_ts(3)-200002;
        a = trial_events(i).feeder_off_ts(choice(i));
        k = (trial_events(i).feeder_dur(choice(i))*1e3)+2;
        trial_events(i).feeder_ts(choice(i)) = a-double(k);
    end
    %%*************************************************************************
    %% *************************************************************************
    SF = zeros(length(trial_events),1);
    for i = 1:length(trial_events)
        SF(i) = trial_events(i).feeder_off_ts(trial_events(i).choice+1);
    end
    trial_events(SF>t_end) = [];
    if D == 1
        trial_events(SF<Sessions(s).PreD_t_end) = [];
    end
    % Some trials need to be removed: those have been mentioned in
    % Sessions.rmTrials and also those happned after Sessions.t_end
    trial_events(rmTrials) = []; 
    
    nTrials = length(trial_events);
    %% Video data
    % Load video information, map the position based on origin and corner,
    % clean the Video from redundant timestamps of before and after the task
    V = load([data_root 'vtm1_' num2str(D+1) '.pvd']);
    % Map the Video positions (limitation and registration)
    Video = D_vid_map(V,Sessions,s);

    timestamp = [Video.timestamp];
    Video((timestamp<trial_events(1).feeder_ts(3)*1e-6) | (timestamp>t_end*1e-6)) = [];

    % Interpolation
    oVideo = Video; clear Video;
    
    % Interpolation of position by upsampling
    Video = InterpVid(oVideo,10);
    
    %% Completing Video and trial_events structures
    % For each cell: matching spike time series to the position time seri
    Video = nSpk_Pos(Video,spikes);

    % Matching trial number to each video sample point
    % matching number of spikes to each trial in trial_events
    % each trial is starting from centre feeder and ending to it too.
    [trial_events,Video] = nSpikes_trialnum(Video,trial_events);

    % Assigning Path choice to the video data points
    Video = Video_choice(Video,trial_events);

    % Assigning bin number to the video data points
    bin_num = binnum(Video,'org', Sessions(1).Origin);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correction of bin number 36 which somtimes happen at the beginning of the
    % trial(change it to bin # 1)
    trialnum = [Video.trialnum];
    for i = 1:nTrials-1
        f36 = find((bin_num == 36)&(trialnum == i));
        f1 = find((bin_num == 1)&(trialnum == i),1);
        bin_num(f36<f1) = 1;
    end
%     Putting bin_num into Video
    for ii = 1:length(Video)
        Video(ii).binnum = bin_num(ii);
    end
    
    % find time occupancy for each bin and also number of spikes
    % per bin for each trial
    trial_events = Occup_nSpkPerBin(Video,trial_events);

    %% rate of the spikes in each bin and each trial
    % dividing of nSpikesPerBin to the TimeOccupPerBin
    % If there is only one data point of an specific bin, corresponding
    % TimeOccupPerBin will be 0. To avoid deviding by 0 and because minimum
    for n1 = 1:nTrials
        trial_events(n1).TimeOccupPerBin(trial_events(n1).TimeOccupPerBin == 0) = 0.003;
        FR = trial_events(n1).nSpikesPerBin./repmat(trial_events(n1).TimeOccupPerBin,nCells,1);
        trial_events(n1).FirRatePerBin = FR;
    end
    cd('C:\Users\saeedeh\Documents\MATLAB\Ramp_task_Drug\Results\cache\')
    save(['cache_' Sessions(s).ID '_' num2str(Sessions(s).Session) '_D' num2str(D)],'spikes','trial_events','oVideo','Video')
    clearvars -except Sessions s D
end