function [TORZ,dur_CF,dur_SF_LR,dur_SF_HR,m_TORZ,SEM_TORZ] = TORZ(trial_events)
% TORZ stands for Time Occupancy on Reward Zones.
% This function gets trial_events and uses its t_CF and t_SF fields to
% give you matrix of TORZ and its mean and SEM over trials.
% 
% INPUT:     -trial_events: an structure containing t_CF and t_SF fields.
% 
% OUTPUTS:   -TORZ: 2 row matrix of spent time on CF and SF respectively in
%             each trial
%            -m_TORZ: average of TORZ over trials.
%            -SEM_TORZ: SEM of m_TORZ over trials.
% 
% @Jan 2018_SH

t = [trial_events.t_CF;trial_events.t_SF;trial_events.t_ramp];
dur_CF = t(2,:)-t(1,:); dur_CF(isnan(dur_CF)) = 0;
dur_SF_LR = t(4,:)-t(3,:); dur_SF_LR(isnan(dur_SF_LR)) = 0;

TORZ = [dur_CF;dur_SF_LR];

% ===========================
% modified the code as below:
% ===========================
% Only keep those SF trails that the feeder was set on LR
feeder = reshape([trial_events.feeder_dur],3,[]); feeder = feeder(1:2,:)';
LR_value = min(min(feeder));
HR_value = max(max(feeder));
ch = [trial_events.choice]'+1; SR = zeros(1,length(ch));
for ii = 1:length(ch)
    SR(ii) = feeder(ii,ch(ii));  %Side-Reward value
end
dur_SF_HR = dur_SF_LR;
dur_SF_LR = dur_SF_LR(SR == LR_value);
dur_SF_HR = dur_SF_HR(SR == HR_value);

m_TORZ = [mean(dur_CF),mean(dur_SF_LR)];
SEM_TORZ = [std(dur_CF,[],2)/sqrt(length(dur_CF)),std(dur_SF_LR,[],2)/sqrt(length(dur_SF_LR))];

return