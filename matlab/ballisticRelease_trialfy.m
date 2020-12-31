% trialfy ballistic task
% origianl date: 2020-02-05
% author: Cleave
function ballisticRelease_trialfy(fname, dstname)
load(fname);
% these are task states
ST_BGN = 1; 
ST_PRT = 2;
ST_FCR = 3; % force ramp
ST_MOV = 4;
ST_HLD = 5;
ST_END = 6;
ST_RST = 7;

avoidVals.fth = 0.0220; % don't know why this value, but it is useless -cg
% to be sure: the trial method
trial_count = sum(diff(Data.TrialNo) ~= 0); 

% find indexes
idx_BGN = find([0 diff(Data.TaskStateCodes.Values)] ~= 0 & Data.TaskStateCodes.Values(1:end) == ST_BGN);
idx_PRT = find([0 diff(Data.TaskStateCodes.Values)] ~= 0 & Data.TaskStateCodes.Values(1:end) == ST_PRT);
idx_FCR = find([0 diff(Data.TaskStateCodes.Values)] ~= 0 & Data.TaskStateCodes.Values(1:end) == ST_FCR);
idx_MOV = find([0 diff(Data.TaskStateCodes.Values)] ~= 0 & Data.TaskStateCodes.Values(1:end) == ST_MOV);
idx_HLD = find([0 diff(Data.TaskStateCodes.Values)] ~= 0 & Data.TaskStateCodes.Values(1:end) == ST_HLD);
idx_END = find([0 diff(Data.TaskStateCodes.Values)] ~= 0 & Data.TaskStateCodes.Values(1:end) == ST_END);
idx_RST = find([0 diff(Data.TaskStateCodes.Values)] ~= 0 & Data.TaskStateCodes.Values(1:end) == ST_RST);

% assign value in each trial
for trial_i = 1:trial_count-1 % remove last trial in case of unfinished
    % display 
    if mod(trial_i/trial_count*100,10) < 1
        fprintf('%03d %% \n',floor(trial_i/trial_count*100));
    end
    trial_bgn = idx_BGN(trial_i);
    trial_end = idx_RST(trial_i);
    
    % **1. Trial time and index
    trial(trial_i).trialno = Data.TrialNo(trial_bgn);
    trial(trial_i).block = Data.BlockNo(trial_bgn);
%    trial(trial_i).outcome = ismember(ST_HLD, [trial(trial_i).states]); %%%!!! question here! do hold represents a stop?
    trial(trial_i).outcome = unique(Data.OutcomeMasks.Success(trial_bgn:trial_end));
    trial(trial_i).combo = Data.ComboNo(trial_end);
    [trial(trial_i).states, trial(trial_i).states_idx, ~] = unique(Data.TaskStateCodes.Values(trial_bgn:trial_end));
     
    % **2. Task:    Target
    target = sort(unique(Data.TaskJudging.Target(5, trial_bgn:trial_end)));     % target_theta
    targetl = sort(unique(Data.TaskJudging.Target(6, trial_bgn:trial_end)));    % target_r
    trial(trial_i).target = target(end); % choose bigger one
    trial(trial_i).targetl = targetl(end);
    %               ForceThreshold % why?
    fth = unique(nonzeros(Data.TaskJudging.Target(4, trial_bgn:trial_end))); 
    fth = setdiff(fth, avoidVals.fth);
    if isempty(fth), fth = nan; end 
    trial(trial_i).fth = fth;
    
    % **3. Trial:   index
    trial(trial_i).bgn = trial_bgn;
    trial(trial_i).end = trial_end; % choose ST_RST rather than ST_END because reset represents a true end
    %               time
    trial(trial_i).time = Data.Time(:, trial_bgn:trial_end);
    %               positions
    trial(trial_i).positionAct = Data.Position.Actual(trial_bgn:trial_end,:)';
    trial(trial_i).positionCnt = Data.Position.Center(:,trial_bgn:trial_end);
    trial(trial_i).jointTorque = Data.Position.JointTorque(:,trial_bgn:trial_end);
    %               velocity
    trial(trial_i).velocity = Data.Velocity.Actual(:,trial_bgn:trial_end);
    %               forces
    trial(trial_i).force = Data.Force.Sensor(:,trial_bgn:trial_end);
    trial(trial_i).force_FtSeq = Data.Force.FtSeq(:,trial_bgn:trial_end);
    

    
end
save(dstname, 'trial', 'ST_*');
end