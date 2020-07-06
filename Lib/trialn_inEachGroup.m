function T = trialn_inEachGroup(type,trial_events,varargin)
% Finding trial numbers in each group of reward and effort.
% 
%            Left       |    Right
% ------------------------------------
%           LR(0inch)   |   LR(0inch)   -> 20 force trials
% ------------------------------------
% Rew:    R1:LR(0inch)  |   R2:HR(0inch)   -> 20 force trials
%         R1:HR(0inch)  |   R2:LR(0inch)   -> 20 force trials
% ------------------------------------
% Eff:    E1:HR(9inch)  |   E2:HR(0inch)   -> 20 force trials
%         E1:HR(0inch)  |   E2:HR(9inch)   -> 20 force trials
% 
% INPUTS:   -type: categorization will be based on type.
%                 -R:Reward, E:Effort
%           -trial_events

% OUTPUT:   -T: is a cell containing the trial number in each group
% 
%               if type is R: T{1:2,1} = R1:{LR;HR} and T{1:2,2} = R2:{HR;LR}                  
%               LR: left side, low hight     HR: right side, low hight
%               HR: left side, low hight     LR: right side, low hight
% 
%               if type is E: T{1:2,1} = E1:{HE;LE} and T{1:2,2} = E2:{LE;HE}                  
%               HE: left side, high rew     LE: right side, high rew
%               LE: left side, high rew     HE: right side, high rew
%
% @Feb 2018- SH


%% Inputs
Reward = [200;600];
RampHeight = [0;9];
bad_trials = [];
assignopts(who, varargin);

LRew = Reward(1); HRew = Reward(2);
LHeight = RampHeight(1); HHeight = RampHeight(2);
trial_events(bad_trials(:)) = [];
nTrials = length(trial_events);
[LR1,HR1,LR2,HR2,HE1,LE1,HE2,LE2] = deal([]);
choice = [trial_events.choice];
%% Reward
if strcmp (type,'R')
    for i = 1:nTrials
        if trial_events(i).heights == [LHeight,LHeight]
            if trial_events(i).feeder_dur == [LRew,HRew,LRew]
                if choice(i) == 0
                    LR1(1,end+1) = i;
                else
                    HR2(1,end+1) = i;
                end
            elseif trial_events(i).feeder_dur == [HRew,LRew,LRew]
                if choice(i) == 0
                    HR1(1,end+1) = i;
                else
                    LR2(1,end+1) = i;
                end
            end
        end
    end
    T = {LR1,HR2;HR1,LR2};
end
%% Effort
if strcmp (type,'E')
    for i = 1:nTrials     
        if trial_events(i).feeder_dur == [HRew,HRew,LRew]
            if trial_events(i).heights == [HHeight,LHeight]
                if choice(i) == 0
                    HE1(1,end+1) = i;
                else
                    LE2(1,end+1) = i;
                end
            elseif trial_events(i).heights == [LHeight,HHeight]
                if choice(i) == 0
                    LE1(1,end+1) = i;
                else
                    HE2(1,end+1) = i;
                end
            end
        end
    end
    T = {HE1,LE2;LE1,HE2};
end

return