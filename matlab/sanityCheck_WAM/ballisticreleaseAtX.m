%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is aimed at testing different directions of perturbation, it
% contains following parts: 
%       1. Find the current problem that the x direction is over-dampped
%       2. Trying some wam configurations so that x direction is less
%       damped;
%       3. Find the 'sweet spot' so that both x and y direction are less
%       damped. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ballisticreleaseAtX
% according to Andy's instruction, try to do x direction ballistic release 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART1, THAT X DIRECTION IS OVERLY DAMPED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpDatarg2(3727);
sstmp = SessionScan(3727);
sstmp.plotTrialfyPositionh_all()
sgtitle('x direction springs, x')
sstmp.plotTrialfyForceh_all()
sgtitle('x direction springs, F')


% use another session as a comparation
sstmp1 = SessionScan(3722);
sstmp1.plotTrialfyPositionh_all();
sgtitle('y direction springs, x')
sstmp1.plotTrialfyForceh_all()
sgtitle('y direction springs, F')

% use another force level, but still x to do the release
sstmp2 = SessionScan(3728);
 sstmp2.plotTrialfyPositionh_all();
 sstmp2.plotTrialfyVelocityh_all();
sgtitle('x direction springs, x')
 sstmp2.plotTrialfyForceh_all();
sgtitle('x direction springs, F')

figure(); 
for i = 1:sstmp2.trials_num
    axh(1) = subplot(3,1,1); hold on;
    plot(sstmp2.trials(i).data.t_shift, sstmp2.trials(i).data.v(1,:));
        axh(2) = subplot(3,1,2); hold on;
    plot(sstmp2.trials(i).data.t_shift, sstmp2.trials(i).data.v(2,:));
        axh(3) = subplot(3,1,3); hold on;
    plot(sstmp2.trials(i).data.t_shift, sstmp2.trials(i).data.v(3,:));
end
linkaxes(axh, 'x');
% acc
figure(); 
for i = 1:sstmp2.trials_num
    axh(1) = subplot(3,1,1); hold on;
    plot(sstmp2.trials(i).data.t_shift(2:end), diff(smooth(sstmp2.trials(i).data.v(1,:))));
        axh(2) = subplot(3,1,2); hold on;
    plot(sstmp2.trials(i).data.t_shift(2:end), diff(smooth(sstmp2.trials(i).data.v(2,:))));
        axh(3) = subplot(3,1,3); hold on;
    plot(sstmp2.trials(i).data.t_shift(2:end), diff(smooth(sstmp2.trials(i).data.v(3,:))));
end
linkaxes(axh, 'x');

% plot out these plot to to show in next time! 

%% 
%%%%%%%%%%%%%%%%%%%
% still in x direction, but decrease the damping (even use the negative),
% to see if the over-damping problem has been solved.

%%% try robot with different perturbation magnitude
%sstmp = SessionScan(3731);
%sstmp = SessionScan(3732); % up to 16N perturbation, not working yet (1 dip at the perturbation)
% sstmp = SessionScan(3733); % up to 20N perturbation, up to -10 Nm/s^2 damping. 
% sstmp = SessionScan(3734); % up to 18N perturbation, up to -5 Nm/s^2 damping. 
%sstmp = SessionScan(3735); % up to 12N/15N perturbation, up to -5 Nm/s^2 damping. 
% sstmp = SessionScan(3736); % up to 12N perturbation, up to -5 Nm/s^2 damping. 
sstmp = SessionScan(3737); % 12N perturbation of 30 trials of step&stoc
% sstmp = SessionScan(3740); % 18N perturbation of 30 trials of step&stoc
% 20N -5Nm/s2 perturbation might work (as it has two dip), (need to feel)
%sstmp = SessionScan(3722);
sstmp.plotTrialfyPositionh_all()
sstmp.plotTrialfyVelocityh_all()
figure(); 
% buttorworth filter
freq = 500;
 cf = 20; % cutoff freqnency
[b,a] = butter(4,cf/(freq/2)); % make filter
% x = filtfilt(b,a,x); % apply fitler

for i = 1:sstmp.trials_num
   
    v1 = filtfilt(b,a,sstmp.trials(i).data.v(1,:));
    v2 = filtfilt(b,a,sstmp.trials(i).data.v(2,:));
    v3 = filtfilt(b,a,sstmp.trials(i).data.v(3,:));
    axh(1) = subplot(3,1,1); hold on;
    plot(sstmp.trials(i).data.t_shift(2:end), diff(v1));
    axh(2) = subplot(3,1,2); hold on;
    plot(sstmp.trials(i).data.t_shift(2:end), diff(v2));
    axh(3) = subplot(3,1,3); hold on;
    plot(sstmp.trials(i).data.t_shift(2:end), diff(v3));
end
title(axh(1), 'a-x');title(axh(2), 'a-y');title(axh(3), 'a-z');
linkaxes(axh, 'x');


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART2, THAT TRYING NEW ROBOT CONFIGURATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try roobt at a differnt configuration: 
% configuration at the following: 
%           jp = [-1.57, -0.785, 0, 1.57];
%           tp = [0, 0.7707, 0]; 

% cpDatarg2(3730);
sstmp = SessionScan(3730);

sstmp.plotTaskEndpointPosition_all()
sstmp.plotTaskJointPosition_all()

%% data from yqr
sstmp = SessionScan(3739); 
sstmp.plotTrialfyPositionh_all()
sstmp.plotTrialfyForceh_all()

%% tidy up ss3737 & 3740. 
data = cell(2,1,1,15,3);
sstmp1 = SessionScan(3737); % 12N perturbation, Br=-5
celltmp1 = sstmp1.export_as_formatted_hybridss(1);
data(1,1,1,:,:) = celltmp1(1:15,:);
% save('data/processedData/ss3737.mat', 'data');



sstmp2 = SessionScan(3740); % 18N perturbation, Br=-5
celltmp2 = sstmp2.export_as_formatted_hybridss(1);
data(2,1,1,:,:) = celltmp2(1:15,:);
save('data/processedData/ss3737_3740.mat', 'data');

%% see the joint positions in +-10 cm at a new configuration: 
% jp = [-1.567,     -0.7614,    0,  1.281]
% cp = [0.07657,    0.7798,     0];
% cpDatarg2(3742)
wamtmp = SessionScanWam(3742);
figure();
axh(1) = subplot(2,1,1); grid on; hold on;
plot(wamtmp.time, wamtmp.jp, 'Marker', '.');
legend('j1', 'j2', 'j3', 'j4');
axh(2) = subplot(2,1,2); grid on; hold on;
plot(wamtmp.time, wamtmp.tp, 'Marker', '.');
legend('x', 'y', 'z');
linkaxes(axh,'x');

%% try the spring test in the new configuration: 
sstmp(1) = SessionScan(3743);          % right
sstmp(1).plotTrialfyPositionh_all()
sstmp(2) = SessionScan(3744);          % left
sstmp(2).plotTrialfyPositionh_all()

sstmp(3) = SessionScan(3746);          % front
sstmp(3).plotTrialfyPositionh_all()
sstmp(4) = SessionScan(3747);          % back
sstmp(4).plotTrialfyPositionh_all()

% TODO: 
%  show the time-x relationship in the 4 directions
% show the trajectories in the 4 directions. 

% position
sstmp(1).plotTrialfyPositionh_all()
sstmp(2).plotTrialfyPositionh_all()
sstmp(3).plotTrialfyPositionh_all()
sstmp(4).plotTrialfyPositionh_all()

% velocity

sstmp(1).plotTrialfyVelocityh_all()
sstmp(2).plotTrialfyVelocityh_all()
sstmp(3).plotTrialfyVelocityh_all()
sstmp(4).plotTrialfyVelocityh_all()

%% plot trajectory over here
figure(); hold on; grid on;
for ss_i = 1:4 
    for trial_i = 1:sstmp(ss_i).trials_num
        idx = sstmp(ss_i).trials(trial_i).data.ts == 5 | ...
                sstmp(ss_i).trials(trial_i).data.ts == 6;
        idxf = find(idx);
        plot(sstmp(ss_i).trials(trial_i).data.x(1,idx),...
            sstmp(ss_i).trials(trial_i).data.x(2,idx));
        plot(sstmp(ss_i).trials(trial_i).data.x(1,idxf(end)),...
            sstmp(ss_i).trials(trial_i).data.x(2,idxf(end)),...
            'Marker', '*');
    end
end
title('trajectories in release');
xlim([-0.05 0.25]);
ylim([[0.65 0.90]]);

%% plot perturbation informations

% buttorworth filter
freq = 500;
 cf = 20; % cutoff freqnency
[b,a] = butter(4,cf/(freq/2)); % make filter
% x = filtfilt(b,a,x); % apply fitler
direction_str = {'right', 'left', 'front', 'back'};
for ssi = 1:4
    fh(ssi) = figure(); 
    sgtitle(direction_str{ssi});
for i = 1:sstmp(ssi).trials_num
    axh(1) = subplot(3,1,1); hold on;
    v1 = filtfilt(b,a,sstmp(ssi).trials(i).data.v(1,:));
    v2 = filtfilt(b,a,sstmp(ssi).trials(i).data.v(2,:));
    v3 = filtfilt(b,a,sstmp(ssi).trials(i).data.v(3,:));
    plot(sstmp(ssi).trials(i).data.t_shift(2:end), diff(v1));
        axh(2) = subplot(3,1,2); hold on;
    plot(sstmp(ssi).trials(i).data.t_shift(2:end), diff(v2));
        axh(3) = subplot(3,1,3); hold on;
    plot(sstmp(ssi).trials(i).data.t_shift(2:end), diff(v3));

%     plot(sstmp(ssi).trials(i).data.t_shift(2:end), diff(smooth(sstmp(ssi).trials(i).data.v(1,:))));
%         axh(2) = subplot(3,1,2); hold on;
%     plot(sstmp(ssi).trials(i).data.t_shift(2:end), diff(smooth(sstmp(ssi).trials(i).data.v(2,:))));
%         axh(3) = subplot(3,1,3); hold on;
%     plot(sstmp(ssi).trials(i).data.t_shift(2:end), diff(smooth(sstmp(ssi).trials(i).data.v(3,:))));

%     axh(1) = subplot(3,1,1); hold on;
%     plot(sstmp(ssi).trials(i).data.t_shift, sstmp(ssi).trials(i).data.v(1,:));
%         axh(2) = subplot(3,1,2); hold on;
%     plot(sstmp(ssi).trials(i).data.t_shift, sstmp(ssi).trials(i).data.v(2,:));
%         axh(3) = subplot(3,1,3); hold on;
%     plot(sstmp(ssi).trials(i).data.t_shift, sstmp(ssi).trials(i).data.v(3,:));
end
linkaxes(axh, 'xy');

end



%% run two other positions 
% cpDatarg2(3748);

% first, scan the positions 
wamtmp = SessionScanWam(3748);
figure(); 
axh(1) = subplot(2,1,1); 
plot(wamtmp.time, wamtmp.tp); 
legend('x', 'y','z');
axh(2) = subplot(2,1,2); 
plot(wamtmp.time, wamtmp.jp); 
legend('1','2','3','4');
linkaxes(axh,'x');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stack all the new configuration data in one matrix 
% do the following spring tests: 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read all the data 
ss_num = [...
            3336, 3724, 3727, 3727
            3743, 3744, 3746, 3747; ...
            3752, 3750, 3751, 3749; ...
            3753, 3755, 3754, 3756;...
            3759, 3761, 3763, 3764]';       % This is the 'sweet point' 
        % ... sequence: R, L, U, D
ssnum = ss_num(:);
for i = 1:length(ssnum)
    sstmp(i) = SessionScan(ssnum(i));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% show all the data 

figure(); hold on; grid on;
for ss_i = 1:length(ssnum)
    for trial_i = 1:sstmp(ss_i).trials_num
        idx = sstmp(ss_i).trials(trial_i).data.ts == 5 | ...
                sstmp(ss_i).trials(trial_i).data.ts == 6;
        idxf = find(idx);
        plot(sstmp(ss_i).trials(trial_i).data.x(1,idx),...
            sstmp(ss_i).trials(trial_i).data.x(2,idx));
        plot(sstmp(ss_i).trials(trial_i).data.x(1,idxf(end)),...
            sstmp(ss_i).trials(trial_i).data.x(2,idxf(end)),...
            'Marker', '*');
    end
end
title('trajectories in release');
    
for i = 1:4
    sstmp(0+i).plotTrialfyVelocityh_all();
end

%% see the 3rd position
cpDatarg2(3757)
wamtmp = SessionScanWam(3757);
figure(); 
axh(1) = subplot(2,1,1); 
plot(wamtmp.time, wamtmp.tp); 
legend('x', 'y','z');
axh(2) = subplot(2,1,2); 
plot(wamtmp.time, wamtmp.jp); 
legend('1','2','3','4');
linkaxes(axh,'x');

% Is the perturbation wrong? 
% cpDatarg2(3759)
wamtmp = SessionScanWam(3760);
axh(1) = subplot(3,1,1); 
plot(wamtmp.time, wamtmp.tp); 
legend('x', 'y','z');
axh(2) = subplot(3,1,2); 
plot(wamtmp.time, wamtmp.jp); 
legend('1','2','3','4');
axh(3) = subplot(3,1,3); 
plot(wamtmp.time, wamtmp.cf); 
legend('1','2','3');
linkaxes(axh,'x');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART3, ANALYSE THE 'SWEET SPOT'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. Export the data into standard spring format
% format: 5D (direction-force-stiffness-trial-perturbation

% note: the data feed into James' code should be 3-by-3 (force vs springs),
% here I just use 1-by-1, and copy and paste it if nessesary 

ss_num = [3759 3761 3763 3764]; % right, left, front, back

for dir_i = 1:4
    ss_tmp = SessionScan(ss_num(dir_i));
    celltmp = ss_tmp.export_as_formatted_hybridss(1); % stochastic data not recorded 
    % TO ACCOMODATE FOR JAMES' CODE (ABLE TO RUN X DIRECTION), SAVE X
    % DIRECTION DATA INTO Y DIRECTION HERE 
    % 
    if (ss_tmp.ssnum == 3759 || ss_tmp.ssnum == 3761)
        if ss_tmp.ssnum == 3761 % move to left direction 
            k = -1;
        else 
            k = 1;
        end
        for i = 1:size(celltmp,2)         % types of pert
            for j = 1:size(celltmp,1)     % trials
                if isempty(celltmp{j,i})
                    continue;
                else
                    % change y to x
                    % x
                    tmp = celltmp{j,i}.x(2,:);
                    celltmp{j,i}.x(2,:) = celltmp{j,i}.x(1,:) * k;
                    celltmp{j,i}.x(1,:) = tmp * k;
                    % v
                    tmp = celltmp{j,i}.v(2,:);
                    celltmp{j,i}.v(2,:) = celltmp{j,i}.v(1,:) * k;
                    celltmp{j,i}.v(1,:) = tmp * k;
                    % f
                    tmp = celltmp{j,i}.f(2,:);
                    celltmp{j,i}.f(2,:) = celltmp{j,i}.f(1,:) * k;
                    celltmp{j,i}.f(1,:) = tmp * k;
                    % Fp
                    tmp = celltmp{j,i}.Fp(2,:);
                    celltmp{j,i}.Fp(2,:) = celltmp{j,i}.Fp(1,:) * k;
                    celltmp{j,i}.Fp(1,:) = tmp * k;
                end
            end
        end
    end
    
    for fce_i = 1:3
        for tar_i = 1:3 % step perts
%         ss_tmp = SessionScan(ss_num(fce_i, tar_i));
%         celltmp = ss_tmp.export_as_formatted_hybridss(1); % stochastic data not recorded 
        trials_num = size(celltmp,1);
        if trials_num>15
            data(dir_i,fce_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(dir_i,fce_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
        end
    end
end
save('data/processedData/ss3759_3764.mat', 'data')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1. Feed the data into James' code 

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3759_3764.mat');
data1 = reshape(data(1,1,:,:,:,:),1,3,3,15,3); % 3, 4, up-down

this = crossConditionAnalysis(data1, 1, 1, 1,'spring');

this.k_hat_pulse
% crossConditionAnalysis(data,dexSubject,dexDirection,dexDistance,subjectType)
% dexSubject: (spring) stiffness levels and (human) target distances
% dexDirection: 

% NEED TO CHECK THE FITTING! 
%% try 
%% % tmptmp: try the code to see if it works 
% Test from Chenguang Full data (11/19/2021)
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3486_3534.mat');
data_human = reshape(data(1,3,:,:,:,:),1,3,3,15,3); % originally 3
clear data
dexSubject = 1; % [Chenguang]
dexForce = 1; % [15N, 20N, 25N]
dexDistance = 1; % [2.5cm, 5cm, 7.5cm]
depMeasures_human2 = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');


%% 
% %%%%%%%%%%%%%%%%%%%%%%% here we try a new posture %%%%%%%%%%%%%%%%%%%%%%%
wamtmp = SessionScanWam(3829);
axh(1) = subplot(3,1,1); 
plot(wamtmp.cf')
axh(2) = subplot(3,1,2); 
plot(wamtmp.tp');
axh(3) = subplot(3,1,3); 
plot(wamtmp.tv');
linkaxes(axh, 'x')

%%
sstmp = SessionScan(3830)
axh(1) = subplot(3,1,1); 
plot(sstmp.wam.cf')
axh(2) = subplot(3,1,2); 
plot(sstmp.wam.tp');
axh(3) = subplot(3,1,3); 
plot(sstmp.wam.tv');
linkaxes(axh, 'x')
%%  visualize joint position vs tool position 
clear axh
figure(); 
axh(1) = subplot(2,1,1); 
plot(sstmp.wam.time, sstmp.wam.jp, 'Marker', '.'); 
legend('1','2','3','4');
axh(2) = subplot(2,1,2); 
plot(sstmp.wam.time, sstmp.wam.tp, 'Marker', '.');
legend('x', 'y', 'z');


%%
% sstmp = SessionScan(3831); %  6N pulse
sstmp = SessionScan(3834); %  6N pulse, B_hold = -5
% sstmp = SessionScan(3832); % 12N pulse
% sstmp = SessionScan(3833); % 8N pulse
sstmp.plotTrialfyVelocityh_all()
axh(1) = subplot(3,1,1); 
plot(sstmp.wam.cf')
axh(2) = subplot(3,1,2); 
plot(sstmp.wam.tp');
axh(3) = subplot(3,1,3); 
plot(sstmp.wam.tv');
linkaxes(axh, 'x')
linkaxes(axh, 'x');
celltmp = sstmp.export_as_formatted_hybridss(1);
for i = 1:size(celltmp,2)         % types of pert
    for j = 1:size(celltmp,1)     % trials
        if isempty(celltmp{j,i})
            continue;
        else
            % change y to x
            % x
            tmp = celltmp{j,i}.x(2,:);
            celltmp{j,i}.x(2,:) = celltmp{j,i}.x(1,:) * k;
            celltmp{j,i}.x(1,:) = tmp * k;
            % v
            tmp = celltmp{j,i}.v(2,:);
            celltmp{j,i}.v(2,:) = celltmp{j,i}.v(1,:) * k;
            celltmp{j,i}.v(1,:) = tmp * k;
            % f
            tmp = celltmp{j,i}.f(2,:);
            celltmp{j,i}.f(2,:) = celltmp{j,i}.f(1,:) * k;
            celltmp{j,i}.f(1,:) = tmp * k;
            % Fp
            tmp = celltmp{j,i}.Fp(2,:);
            celltmp{j,i}.Fp(2,:) = celltmp{j,i}.Fp(1,:) * k;
            celltmp{j,i}.Fp(1,:) = tmp * k;
        end
    end
end

for i = 1:3
    for j = 1:3
        data(1,1,i,j,:,:) = celltmp(1:15,:);
    end
end
% pertend to save data here 
data1 = reshape(data(1,1,:,:,:,:),1,3,3,15,3); % 3, 4, up-down

this = crossConditionAnalysis(data1, 1, 1, 1,'spring');

this.k_hat_pulse


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try another robot configurations q = [0 0 0 1.57] %%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss_num = {  [3854,3855 ], [ 3840,3841 ], [ 3842,3843,3844];
            [3853      ], [ 3838,3839 ], [ 3848,3849];       
            [3852      ], [ 3836,3837 ], [ 3850,3851];      };
      
% try if all sessions work 
for i = 1:3 
    for j = 1:3 
        for k = length(ss_num{i,j})
            try 
            sstmp = SessionScan(ss_num{i,j}(k));
            catch 
                disp(['Session' num2str(ss_num{i,j}(k)) 'dont work!']);
            end
        end
    end
end

%% overlap all the celltmp data and sort depending on the magnitude of the perturbation 

ss_num = {  [3854,3855 ], [ 3840,3841 ], [ 3842,3843,3844];
            [3853      ], [ 3838,3839 ], [ 3848,3849];       
            [3852      ], [ 3836,3837 ], [ 3850,3851];      };
pert_f = 12; % only use 12N pert
pertM_num = 1 + 5;     % 1 without pert, 5 pert force combinations
data = cell(1, size(ss_num,1), size(ss_num,2), 5, pertM_num+1); % The stochastic ones are attached at the end 
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        paramtmp= cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
            [celltmptmp, paramtmptmp] = ss_tmp.export_as_formatted_hybridss(1);
            % .... TODO: wirte another function that return the K, B so
            % that I can sort...
            
            % check the size of celltmptmp
            cell_avail_num = zeros(1,3);
            for i = 1:size(celltmptmp,1)
                for j = 1:size(celltmptmp,2)
                    if ~isempty(celltmptmp{i,j})
                        cell_avail_num(j) = cell_avail_num(j) + 1;
                    end
                end
            end
            
            for j = 1:size(celltmptmp,2)
                celltmp(cell_idx_from(j)+(1:cell_avail_num(j)),j) = ...
                    celltmptmp(1:cell_avail_num(j),j);
                paramtmp(cell_idx_from(j)+(1:cell_avail_num(j)),j) = ...
                    paramtmptmp(1:cell_avail_num(j),j);
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varParam = celltmp(:,2);
        paramtmp_varParam= paramtmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varParam)
            if isempty(celltmp_varParam{ti})
                continue;
            end
            pertMag(ti) = 0;
            WamBp(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varParam{ti}.t, celltmp_varParam{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varParam{ti}.t, celltmp_varParam{ti}.Fp);
            end
            if max(abs(celltmp_varParam{ti}.Fp(1,:))) ~= 0
                pertMag(ti) = max(abs(celltmp_varParam{ti}.Fp(1,:))); % 500Hz
                WamBp(ti) = paramtmp_varParam{ti}.wamBp;
            end
        end
        % classify the pulse time into cells that have different Magnitude
        % and wamBp
        
        % get the index of each delay interval 
        [pertmag_unq, ia, ic] = unique(pertMag + WamBp*10);
        pertmag_unq
        idx_trialsPertT = cell(1, length(pertmag_unq));
        trials_num_max = 0;
        for pi = 1:length(pertmag_unq) 
            idx_trialsPertT{pi} = find((pertMag + WamBp*10) == pertmag_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertmag_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertmag_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end 
        
        % now save the celltmp in a x-direction manner 
        for i = 1:size(celltmp,2)         % types of pert
            for j = 1:size(celltmp,1)     % trials
                if isempty(celltmp{j,i})
                    continue;
                else
                    % change y to x
                    % x
                    tmp = celltmp{j,i}.x(2,:);
                    celltmp{j,i}.x(2,:) = celltmp{j,i}.x(1,:) * k;
                    celltmp{j,i}.x(1,:) = tmp * k;
                    % v
                    tmp = celltmp{j,i}.v(2,:);
                    celltmp{j,i}.v(2,:) = celltmp{j,i}.v(1,:) * k;
                    celltmp{j,i}.v(1,:) = tmp * k;
                    % f
                    tmp = celltmp{j,i}.f(2,:);
                    celltmp{j,i}.f(2,:) = celltmp{j,i}.f(1,:) * k;
                    celltmp{j,i}.f(1,:) = tmp * k;
                    % Fp
                    tmp = celltmp{j,i}.Fp(2,:);
                    celltmp{j,i}.Fp(2,:) = celltmp{j,i}.Fp(1,:) * k;
                    celltmp{j,i}.Fp(1,:) = tmp * k;
                end
            end
        end
        

        data(1,frc_i,dist_i,1:15,1:pertM_num) = celltmp1(1:15,1:pertM_num);
        data(1,frc_i,dist_i,1:15,pertM_num+1) = celltmp(1:15,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3836_3855.mat', 'data'); % 12N perturbation, various time

% 
%% tmptmp %% now save the data so that x-direction overwrite at y position 
load('data/processedData/ss3836_3855.mat', 'data');
% 
k = -1;
for frc_i = 1:3
    for dist_i = 1:3
        for pert_i = 1:7
            celltmp = reshape(data(1,frc_i,dist_i,1:15,pert_i), 15, 1);
            for j = 1:size(celltmp,1)     % trials
                if isempty(celltmp{j})
                    continue;
                else
                    % change y to x
                    % x
                    tmp = celltmp{j}.x(2,:);
                    celltmp{j}.x(2,:) = celltmp{j}.x(1,:) * k;
                    celltmp{j}.x(1,:) = tmp * k;
                    % v
                    tmp = celltmp{j}.v(2,:);
                    celltmp{j}.v(2,:) = celltmp{j}.v(1,:) * k;
                    celltmp{j}.v(1,:) = tmp * k;
                    % f
                    tmp = celltmp{j}.f(2,:);
                    celltmp{j}.f(2,:) = celltmp{j}.f(1,:) * k;
                    celltmp{j}.f(1,:) = tmp * k;
                    % Fp
                    tmp = celltmp{j}.Fp(2,:);
                    celltmp{j}.Fp(2,:) = celltmp{j}.Fp(1,:) * k;
                    celltmp{j}.Fp(1,:) = tmp * k;
                end
            end
            data(1,frc_i,dist_i,1:15, pert_i) = celltmp;
        end
    end
end
save('data/processedData/ss3836_3855_cvtxy.mat', 'data');

%% show how the data look like 

% data: 1: no-pert; 2: 6N pulse B=-5; 3: 6N pulse B=0; 4: 8N pulse B=0; 
%       5: 10N pulse B=0; 6: 12N pulse B=0;     7: stoc
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3836_3855_cvtxy.mat');
% clear data
pert_result = zeros(5,15);
for pert_i = 1:5
    pert_type = [1,pert_i+1,7];
    data_human = reshape(data(:,:,:,:,pert_type),1,3,3,15,3); % originally 3
    dexSubject = 1; % [Chenguang]
    dexForce = 1; % [15N, 20N, 25N]
    dexDistance = 1; % [2.5cm, 5cm, 7.5cm]
    depMeasures_human2 = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'spring');
    pert_result(pert_i,:) = reshape([depMeasures_human2.k_hat_pulse],1,15);
end


fh = figure();
subplot(1,2,1);
plot(pert_result', 'o'); legend('6N B-5','6N B0', '8N B0', '10N B0', '12N B0');
xlabel('trial');
ylabel('Estimated stiffness');
title('Stiffness Estimation on different perturbation');

subplot(1,2,2);
errorbar(nanmean(pert_result'), nanstd(pert_result'))
xlim([0,6]);
xticks([1:5])
xticklabels({'6N B-5','6N B0', '8N B0', '10N B0', '12N B0'});
ylabel('Estimated stiffness');
title('Stiffness Estimation mean and std');


saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_WAM/tryingNewWamConfiguration/StiffnessEstOnConf5.png');
% Cannot tell which one is better though.

%% show the trajectories 
sstmp = SessionScan(3838);
sstmp.plotTrialfyVelocityh_all()
xlim([-1 1])

figure(); 
for i = 1:sstmp.trials_num
    Fpeak = max(abs(sstmp.trials(i).data.Fp(1,:))); 
    wamBp = sstmp.trials(i).wamBp;
    if (Fpeak == 6 && wamBp == 0)
        axh(1) = subplot(3,1,1); hold on;
    plot(sstmp.trials(i).data.t_shift(2:end), diff(smooth(sstmp.trials(i).data.v(1,:))));
        axh(2) = subplot(3,1,2); hold on;
    plot(sstmp.trials(i).data.t_shift(2:end), diff(smooth(sstmp.trials(i).data.v(2,:))));
        axh(3) = subplot(3,1,3); hold on;
    plot(sstmp.trials(i).data.t_shift(2:end), diff(smooth(sstmp.trials(i).data.v(3,:))));
    end
end
linkaxes(axh, 'x');
xlim([-4 -2]);