% According to the meeting at 2021-12-21, I start to try the pulse during
% the movement. 
%
% In this script, I'll collect(tidy up) the data in the perturbation with
% different magnitude (to find out the better one), and maybe different
% perturbation time when the handle was released during motion. 
% The current thought is to fix the duration of the perturbation. 

%% Find the possible perturbation magnitude using a short perturbation 
% pert duration: 010s, pert magnitude: 6,12,18,24,30
% both perturbation at hold and at release 

% 1. Check if release has been perturbed, whether it influence the movement
clear all; close all; clc;
wamtmp = SessionScanWam(3766);

figure();
axh(1) = subplot(3,1,1);
plot(wamtmp.time, wamtmp.state);
axh(2) = subplot(3,1,2);
plot(wamtmp.time, wamtmp.cf(2,:));
axh(3) = subplot(3,1,3);
plot(wamtmp.time, wamtmp.tv(2,:));

linkaxes(axh, 'x');

% did not see change here. What about overlay the velocities?  

sstmp = SessionScan(3766);
sstmp.plotTrialfyVelocityh()
% is the force change according to the perturbation? 
sstmp.plotTrialfyForceh()

%% A session with both perturbed and non-perturbed trials ... 
clear all; close all; clc;

% cpDatarg2(3767);
sstmp = SessionScan(3767); 
sstmp.plotTrialfyVelocityh()

%% %%%%%%%%%%%%%%%% The first try of variate spring stiffness %%%%%%%%%%%%
% The same parameter with the other spring tests, and have perturb during
% release (at the 100ms of the release). 
% Export to the format so that James is able to work 
% 
%  ss_num = [  3774        3768        3775
%              3773        3770        3777
%              3772        3771        3778];     % 6N cases
%  ss_num = [  3774        3769        3775                   
%              3773        3770        3777
%              3772        3771        3778];     % 12N cases % tocopy: 3768, 3769
 ss_num = [  3787        3782        3781
             3786        3783        3780
             3785        3784        3779];     % 18N,20N cases

pert_f = 24; % only use 6N perturb
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        ss_tmp = SessionScan(ss_num(frc_i, dist_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
        
        if (ss_num(frc_i, dist_i) == 3770)
            celltmp(12:15,1) = celltmp(1:4,1);
        end
        
        trials_num = size(celltmp,1);
        % check the celltmp, only use 6N perturb 
        trials_Pf = zeros(1,(trials_num)); 
        for ti = 1:trials_num
            if isempty(celltmp{ti,2})
                trials_Pf(ti) = nan;
                continue
            end
            trials_Pf(ti) = max(abs(celltmp{ti,2}.Fp(2,:))); % could only be 6 or 12
            plot(celltmp{ti,2}.t,celltmp{ti,2}.Fp(:,:));
        end
        celltmp1 = cell(size(celltmp));
        idx = find(trials_Pf == pert_f);
        celltmp1(:,1) = celltmp(:,1);
        celltmp1(1:length(idx),2) = celltmp(idx,2);
        celltmp1(:,3) = celltmp(:,3);
        
        if trials_num>15
            data(1,frc_i,dist_i,:,:) = celltmp1(1:15,:);
        else
            data(1,frc_i,dist_i,1:trials_num,:) = celltmp1(:,:);
        end
    end
end
%save('data/processedData/ss3768_3770.mat', 'data'); % 6N perturbation
%save('data/processedData/ss3769_3770.mat', 'data');   % 12N perturbation

% save('data/processedData/ss3779_3783.mat', 'data'); % 18N perturbation
 save('data/processedData/ss3784_3787.mat', 'data'); % 24N perturbation

 %% tidy 4 force levels in one file 
 %...TO BE FINISHED
% load('data/processedData/ss3784_3787.mat', 'data'); % 24N perturbation
 
%% display them ...
load('data/processedData/ss3768_3770.mat', 'data');  % 6N perturbation on x, 
% load('data/processedData/ss3769_3770.mat', 'data');  % 12N perturbation on x, 
% load('data/processedData/ss3779_3783.mat', 'data'); % 18N perturbation
% load('data/processedData/ss3784_3787.mat', 'data'); % 24N perturbation
Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 
colors = colormap('lines');
r = size(Data, 1); % subj
c = size(Data, 2); % direction
d = size(Data, 3); % target
l = size(Data, 4); % trials
p = size(Data, 5); % perturbation type
idx_last = 200;
if_subtract = 0;
epoc_type = 2;
plot_type = 4; % 1displacement
for ri = 1:r % subj
    for ci = 1:c % direction
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        axh(ri, ci) = subplot(c,r,r*(ci-1) + ri);grid on;hold on;
%         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
        for di = 1:d % target distance
        %for di = 3 % target distance
            for li = 2%1:p % perturbation
                trial_num = length(Data(ri,ci,di,:,li));
                for ti = 1:trial_num % each trial
                    if (isempty(Data{ri,ci,di,ti,li}))
                        continue;
                    end
                    
                    switch epoc_type
                        case 1
                            idx = find(Data{ri,ci,di,ti,li}.Fp(2,:)~=0 & ...  
                                Data{ri,ci,di,ti,li}.ts==4);  % pert at y
                            idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                            if li == 1
                                display('ERROR: should use li == 2!!!');
                            end
                        case 2
                            idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                            idx = (idx(1)-100):idx(end);
                            %idx = (idx(1)):idx(end);
                    end
                    %plot(Data{ri,ci,di,ti,li}.Fp(2,:));
                    %idx = find(Data{ri,ci,di,ti,li}.Fp(2,:)~=0);
                    %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                    %idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                    %idx = (idx-9):idx(end);
                    %idx = (idx(1)):(idx(end)+100);
                    time = t_step*(idx-idx(1));
                    %time = idx-idx(1);
                    switch plot_type
                        case 2
                            dat = Data{ri,ci,di,ti,li}.f(2,idx);
                            titlestr = 'force';
                        case 1
                            dat = Data{ri,ci,di,ti,li}.x(2,idx);
                            %dat = dat - dat(1);
                            titlestr = 'displacement';
                        case 3
                            dat = Data{ri,ci,di,ti,li}.Fp(2,idx);
                            titlestr = 'Fp';
                        case 4
                            dat = Data{ri,ci,di,ti,li}.v(2,idx);
                            titlestr = 'velocity';
                        case 5
                            dat = Data{ri,ci,di,ti,li}.tq(3,idx);
                            titlestr = 'torque3';
                        case 6
                            dat = Data{ri,ci,di,ti,li}.x(:,idx);
                            dat_submean = dat - mean(dat(:,1:50),2);
                            dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                            dat = dat_norm .* sign(dat_submean(2,:));
                            titlestr = 'norm displacement';
                        case 7 % the force mode
                            dat = Data{ri,ci,di,ti,li}.f(:,idx);
                            dat_submean = dat - mean(dat(:,1:50),2);
                            dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                            dat = dat_norm .* sign(dat_submean(2,:));
                            titlestr = 'norm force';
                        
                    end
                    if (if_subtract)
                        dat = dat - mean(dat(1:50));
                    end
                    plot(time, dat, 'Color', colors(4*(li-1)+di, :));
%                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                end
            end
        end
    end
end
%xlim([0 0.7])
linkaxes(axh, 'xy');
%xlim([0 0.5]);
xlim([0 2])
sgtitle(titlestr);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PULSE DURING MOVEMENT FOR VARIABLE TIME               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this case, I did the pulse during movement for variable time. The
% pulse peaks at 0.15:0.05:0.60s.  
% plot the Velocity 

% I'm interested in: 
% 1. What are the each plot against the original one. 
% 2. How to isolate the perturbation results from the original one?

%% plot them in 3 figures 
sstmp(1) = SessionScan(3793); 
sstmp(2) = SessionScan(3794); 
sstmp(3) = SessionScan(3795);

close all;
sstmp(1).plotTrialfyVelocityh();
xlim([-0.1 1.0]);
sstmp(2).plotTrialfyVelocityh();
xlim([-0.1 1.0])
sstmp(3).plotTrialfyVelocityh();
xlim([-0.1 1.0])

%% plot them and tidy them up in format 
ss_num = [   3793        3794        3795];
            % 3793        3793        3793];
            %3785        3784        3779];     % 18N,20N cases
pert_f = 12; % only use 6N perturb
pertT_num = 1 + 10;     % 1 without pert, and 10 perturbation time
data = cell(size(ss_num,1), size(ss_num,2), 5, pertT_num);
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        ss_tmp = SessionScan(ss_num(frc_i, dist_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(frc_i,dist_i,1:5,1:pertT_num) = celltmp1(1:5,1:pertT_num);
    end
end
save('data/processedData/ss3793_3795.mat', 'data'); % 12N perturbation, various time
 
 %% plot the curve in a 2-D form and draw the video  
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3793_3795.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 10;     % 1 without pert, and 10 perturbation time
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/changing_Pert.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);
for fce_i = 1:size(data,1)
    fh(fce_i) = figure();
    for pi = 1:pertT_num
       % 1. plot the perturbed force in the first panel 
       clf; hold on;
       axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(fce_i,1,:,:),5,pertT_num);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       % 2. plot the perturbed velocity in the following panels 
        for dist_i = 1:size(data,2) % for each  
            axh(dist_i+1) = subplot(4,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(fce_i,dist_i,:,:),5,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
            end
        end
        
        % plot notes here: 
        linkaxes(axh, 'x'); xlim([-0.1 1]);
        linkaxes(axh(2:end), 'y'); ylim([-0.8 1]); 
        subplot(4,1,1); title('F-pert'); ylabel('N');
        subplot(4,1,2); title('velocity'); ylabel('m/s');
        subplot(4,1,3); title('velocity'); ylabel('m/s');
        subplot(4,1,4); title('velocity'); ylabel('m/s');
        xlabel('time');
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
end
close(v);

 %%
 %%%%%%%%%%%% plot the curve in a 2-D form and draw the video, the velocity subtraction 
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3793_3795.mat', 'data');
% v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/changing_Pert_SpringVel_subtract.mp4', 'MPEG-4');
% pertT_num = 1 + 10;     % 1 without pert, and 10 perturbation time
color_arr = colormap('lines');
close all;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3913_3921.mat', 'data'); pertT_num = 1+7; % 7 perturbs
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/SubjectTimeChangingPert_3-by-3_200ms.mp4', 'MPEG-4');
data = reshape(data(1,1,:,:,:,:),3,3,10,8);
num_trials = 10;
v.FrameRate = 1;
open(v);
for fce_i = 1:size(data,1)
    fh(fce_i) = figure();
    for pi = 1:pertT_num
       % 1. plot the perturbed force in the first panel 
       clf; hold on;
       axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(fce_i,1,:,:),num_trials,pertT_num);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       % 2. plot the perturbed velocity in the following panels 
       
            
        for dist_i = 1:size(data,2) % for each spring 
            axh(dist_i+1) = subplot(4,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(fce_i,dist_i,:,:),num_trials,pertT_num);
            
            % Get the un-perturbed avg velocity 
            v_mean = zeros(1,500);  
            v_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
            end
            v_mean = mean(v_mean_mat);
            t_mean = 0.002:0.002:1; % 500 data points 
            
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                plot(t_mean, celltmp1{ti,1}.v(2,idx_aftrelease) - v_mean,  'color', [0.5 0.5 0.5]);
%                 idx_release = find(celltmp1{ti,1}.ts == 5);
%                 t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
                

            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi}.ts == 5);
%                 t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_aftrelease = find(celltmp1{ti,pi}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                plot(t_mean, celltmp1{ti,pi}.v(2,idx_aftrelease) - v_mean,  'color',color_arr(4+dist_i,:));
            end
        end
        
        % plot notes here: 
        linkaxes(axh, 'x'); xlim([-0.1 1]);
        linkaxes(axh(2:end), 'y'); ylim([-0.1 0.1]); 
        subplot(4,1,1); title('F-pert'); ylabel('N');
        subplot(4,1,2); title('\Delta velocity'); ylabel('m/s');
        subplot(4,1,3); title('\Delta velocity'); ylabel('m/s');
        subplot(4,1,4); title('\Delta velocity'); ylabel('m/s');
        xlabel('time');
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
close(v);

 %%
 %%%%%%%%%%%% plot the curve in a 3-D form and draw the video, the velocity SUBTRACTION 
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3793_3795.mat', 'data');
% pertT_num = 1 + 10;     % 1 without pert, and 10 perturbation time
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);

color_arr = colormap('lines');
close all;

for dist_i = 1:size(data,2) % for each spring 
    fh(dist_i) = figure(); hold on;
    for fce_i = 1:size(data,1)

        for pi = 1:pertT_num
            % 1. get the data tobe plotted
            celltmp1 = reshape(data(fce_i,1,:,:),size(data,3),pertT_num);
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            
            % Get the un-perturbed avg velocity
            v_mean = zeros(1,500);
            v_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
            end
            v_mean = mean(v_mean_mat);
            t_mean = 0.002:0.002:1; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);  
                
                plot3 (-pi*ones(size(t_mean)), t_mean, celltmp1{ti,1}.v(2,idx_aftrelease) - v_mean,  'color', [0.5 0.5 0.5]);
                
                dat_mat_unpert(ti,:) = celltmp1{ti,1}.v(2,idx_aftrelease) - v_mean;
                
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,pi}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                
                dat_mat_unpert(ti,:) = celltmp1{ti,pi}.v(2,idx_aftrelease) - v_mean;
                
                plot3( -pi*ones(size(t_mean)), t_mean, celltmp1{ti,pi}.v(2,idx_aftrelease) - v_mean,  'color',color_arr(dist_i+4,:));
            end
        end
%         view(120, 57);
        view(58.2, 68.4);
         xlabel('pert time');
         ylabel('time');
         zlabel('velocity');
%         frame = getframe(gcf);
%         writeVideo(v,frame);
    end
    
end
% close(v);

%% Velocity subtraction of longer pulse spring data 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);

color_arr = colormap('lines');
close all;

for dist_i = 1:size(data,2) % for each spring 
    fh(dist_i) = figure(); hold on;
    for fce_i = 1:size(data,1)

        for pi = 1:pertT_num
            % 1. get the data tobe plotted
            celltmp1 = reshape(data(fce_i,1,:,:),size(data,3),pertT_num);
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            
            % Get the un-perturbed avg velocity
            sec = 3; 
            freq = 500;
            v_mean = zeros(1,sec*freq);
            v_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:sec*freq);
%                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                v_mean_mat = [v_mean_mat; celltmp1{ti,1}.x(2, idx_aftrelease)];
            end
            v_mean = mean(v_mean_mat);
            t_mean = 0.002:0.002:sec; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:sec*freq);  
                
%                 plot3 (-pi*ones(size(t_mean)), t_mean, celltmp1{ti,1}.v(2,idx_aftrelease) - v_mean,  'color', [0.5 0.5 0.5]);
                plot3 (-pi*ones(size(t_mean)), t_mean, celltmp1{ti,1}.x(2,idx_aftrelease) - v_mean,  'color', [0.5 0.5 0.5]);
                
                dat_mat_unpert(ti,:) = celltmp1{ti,1}.v(2,idx_aftrelease) - v_mean;
                
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,pi}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:sec*freq);
                
                dat_mat_unpert(ti,:) = celltmp1{ti,pi}.v(2,idx_aftrelease) - v_mean;
                
                plot3( -pi*ones(size(t_mean)), t_mean, celltmp1{ti,pi}.v(2,idx_aftrelease) - v_mean,  'color',color_arr(dist_i+4,:));
            end
        end
%         view(120, 57);
%         view(58.2, 68.4);
         view(90, 0);
         xlabel('pert time');
         ylabel('time');
         zlabel('velocity');
%         frame = getframe(gcf);
%         writeVideo(v,frame);
    end
    
end
% close(v);


%% plot the curve out in a 3-D form
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3793_3795.mat', 'data');
color_arr = colormap('lines');
close all;
for fce_i = 1:size(data,1)
    for dist_i = 1:size(data,2) % for each spring 
        figure(); hold on;
        for pi = 1:11%length(pertT_unq)
%        fh(pi) = figure(); hold on;
       
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(fce_i,1,:,:),5,11);
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
        
%             axh(dist_i+1) = subplot(3,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(fce_i,dist_i,:,:),5,11);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)); 
                    % 4, for consistant with previous color 
            end
            %xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        ylim([-0.1 1]);
        xlabel('perturb positions');
        ylabel('release time');
        zlabel('endpoint velocity(m/s)');
        view(120, 57);
    end
    
end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   NOW THE REAL SUBJECT REPLACE THE SPIRING              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cpDatarg2(3796);
% rough look here 


% ss_tmp = SessionScan(3796);
% ss_tmp.plotTrialfyVelocityh();
celltmp = ss_tmp.export_as_formatted(1);


% Try if this part can tidy it up???
celltmp_varT = reshape(celltmp(1,1,:,2), 1, size(celltmp,3));
        clear pertT
        for ti = 1:length(celltmp_varT)
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts >= 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        
        celltmp1(1:trials_num_max,1) = reshape(celltmp(1,1,1:trials_num_max,1), trials_num_max,1);
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp_varT(idx_trialsPertT{pi});
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp_varT(idx_trialsPertT{pi});
            end
        end
%         data(frc_i,dist_i,1:5,1:11) = celltmp1(1:5,1:11);

%% plot out on different delay

color_arr = colormap('lines');
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/changing_Pert_subj.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);
% for fce_i = 1:size(data,1)
    for pi = 1:11%length(pertT_unq)
%        fh(pi) = figure(); hold on;
       clf; hold on;
       axh(1) = subplot(2,1,1); hold on;                     % plot PF
%        celltmp1 = reshape(data(fce_i,1,:,:),5,11);
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
%         for dist_i = 1:size(data,2) % for each spring 
            axh(2) = subplot(2,1,2); hold on;         % plot each response
%             celltmp1 = reshape(data(fce_i,dist_i,:,:),5,11);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
            end
            %xlim([-0.1 1]); ylim([-0.8 1.0]);
%         end
        
        % plot notes here: 
        linkaxes(axh, 'x'); xlim([-0.1 1]);
        linkaxes(axh(2:end), 'y'); ylim([-0.8 1]); 
        subplot(2,1,1); title('F-pert'); ylabel('N');
        subplot(2,1,2); title('velocity'); ylabel('m/s');
        xlabel('time');
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
% end
close(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%  Do a multiple target distance session 
ss_num  = [ 3801   3800   3802];
for dist_i = 1:3
    
ss_tmp = SessionScan(ss_num(dist_i));
% ss_tmp.plotTrialfyVelocityh();
celltmp = ss_tmp.export_as_formatted(1);


% Try if this part can tidy it up???
celltmp_varT = reshape(celltmp(1,1,:,2), 1, size(celltmp,3));
        clear pertT
        for ti = 1:length(celltmp_varT)
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts >= 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        
        celltmp1(1:trials_num_max,1) = reshape(celltmp(1,1,1:trials_num_max,1), trials_num_max,1);
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp_varT(idx_trialsPertT{pi});
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp_varT(idx_trialsPertT{pi});
            end
        end
        
        data(1,dist_i,1:5,1:12) = celltmp1(1:5,1:12);
end
% save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3800_3802.mat', 'data');

% then, save the data into the 6-D format 
data1 = cell(1,1,1,3,5,12); % subj-direction-frc-dist-trials-pert
data1(1,1,1,1:3,1:5,1:12) = data;
clear data
data = data1;
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3800_3802.mat', 'data');

%% plot out on different delay
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3800_3802.mat', 'data');
color_arr = colormap('lines');
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/changing_Pert_subj_1f3dist.mp4', 'MPEG-4');
v.FrameRate = 1;
% open(v);
% for fce_i = 1:size(data,1)

    for pi = 1:12%length(pertT_unq)
%        fh(pi) = figure(); hold on;
       clf; hold on;
       axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(1,1,1,1,:,:),5,12);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
       for dist_i = 1:size(data,4) % for each spring 
            axh(dist_i+1) = subplot(4,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(1,1,1,dist_i,:,:),5,12);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);


            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
            end
%             xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        linkaxes(axh, 'x'); xlim([-0.1 1]);
        linkaxes(axh(2:end), 'y'); ylim([-0.3 0.6]); 
%         ylim([0.47 0.60]);
        subplot(4,1,1); title('F-pert'); ylabel('N');
        subplot(4,1,2); title('velocity'); ylabel('m/s');
        subplot(4,1,3); title('velocity'); ylabel('m/s');
        subplot(4,1,4); title('velocity'); ylabel('m/s');
        
        xlabel('time');
        frame = getframe(gcf);
%          writeVideo(v,frame);
    end
    
% end
% close(v);


%% plot on different delay on 3D
color_arr = colormap('lines');
close all;
for fce_i = 1:size(data,3)
    for dist_i = 1:size(data,4) % for each spring 
        figure(); hold on;
        for pi = 1:12%length(pertT_unq)
%        fh(pi) = figure(); hold on;
       
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(1,1,fce_i,1,:,:),5,12);
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
        
%             axh(dist_i+1) = subplot(3,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(1,1,fce_i,dist_i,:,:),5,12);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.x(2,:), 'color', [0.5 0.5 0.5]);
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)); 
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)); 
                    % 4, for consistant with previous color 
            end
            %xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        ylim([-0.1 1]);
        xlabel('perturb positions');
        ylabel('release time');
        zlabel('endpoint velocity(m/s)');
        view(120, 57);
    end
    
end
%% plot the velocity difference (the current velocity - the unperturbed avg) velocity subtraction
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3800_3802.mat', 'data');
color_arr = colormap('lines');
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/changing_Pert_subj_1f3dist_vsubst.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);
% for fce_i = 1:size(data,1)

    for pi = 1:12%length(pertT_unq)
%        fh(pi) = figure(); hold on;
       clf; hold on;
       axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(1,1,1,1,:,:),5,12);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
       for dist_i = 1:size(data,4) % for each spring 
            axh(dist_i+1) = subplot(4,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(1,1,1,dist_i,:,:),5,12);
            
            % Get the un-perturbed avg velocity 
            v_mean = zeros(1,500);  
            v_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
            end
            v_mean = mean(v_mean_mat);
            t_mean = 0.002:0.002:1; % 500 data points 
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
%                 idx_release = find(celltmp1{ti,1}.ts == 5);
%                 t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                plot(t_mean, celltmp1{ti,1}.v(2,idx_aftrelease)-v_mean, 'color', [0.5 0.5 0.5]);

            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
%                 idx_release = find(celltmp1{ti,pi}.ts == 5);
%                 t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_aftrelease = find(celltmp1{ti,pi}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                plot(t_mean, celltmp1{ti,pi}.v(2,idx_aftrelease)-v_mean, 'color',color_arr(dist_i+4,:));
            end
%             xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        linkaxes(axh, 'x'); xlim([-0.1 1]);
        linkaxes(axh(2:end), 'y'); ylim([-0.1 0.1]); 
%         ylim([0.47 0.60]);
        subplot(4,1,1); title('F-pert'); ylabel('N');
        subplot(4,1,2); title('\Delta velocity'); ylabel('m/s');
        subplot(4,1,3); title('\Delta velocity'); ylabel('m/s');
        subplot(4,1,4); title('\Delta velocity'); ylabel('m/s');
        
        xlabel('time');
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
% end
close(v);

%% % plot the velocity difference (the current velocity - the unperturbed avg) velocity SUBTRACTION
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3800_3802.mat', 'data');
color_arr = colormap('lines');
close all;
% v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/changing_Pert_subj_1f3dist_vsubst.mp4', 'MPEG-4');
% v.FrameRate = 1;
% open(v);
% for fce_i = 1:size(data,1)
for dist_i = 1:size(data,4) % for each spring 
    fh(dist_i) = figure(); hold on;
    for pi = 1:12%length(pertT_unq)
%        fh(pi) = figure(); hold on;
%        clf; hold on;
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(1,1,1,1,:,:),5,12);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
       
%             axh(dist_i+1) = subplot(4,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(1,1,1,dist_i,:,:),5,12);
            
            % Get the un-perturbed avg velocity 
            v_mean = zeros(1,500);  
            v_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
                v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
            end
            v_mean = mean(v_mean_mat);
            t_mean = 0.002:0.002:1; % 500 data points 
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
%                 idx_release = find(celltmp1{ti,1}.ts == 5);
%                 t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
                idx_aftrelease = find(celltmp1{ti,1}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
%                 plot(t_mean, celltmp1{ti,1}.v(2,idx_aftrelease)-v_mean, 'color', [0.5 0.5 0.5]);
%                 plot3(-pi*ones(size(t_mean)),t_mean, celltmp1{ti,1}.v(2,idx_aftrelease)-v_mean, 'color', [0.5 0.5 0.5]); 
                  plot3(-pi*ones(size(t_mean)),t_mean, celltmp1{ti,1}.v(2,idx_aftrelease)-v_mean, 'color', [0.9 0.9 0.9]); 
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
%                 idx_release = find(celltmp1{ti,pi}.ts == 5);
%                 t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_aftrelease = find(celltmp1{ti,pi}.ts >= 5);
                idx_aftrelease = idx_aftrelease(1:500);
%                 plot(t_mean, celltmp1{ti,pi}.v(2,idx_aftrelease)-v_mean, 'color',color_arr(dist_i+1,:));
%                 plot3(-pi*ones(size(t_mean)),t_mean, celltmp1{ti,pi}.v(2,idx_aftrelease)-v_mean, 'color', color_arr(dist_i+1,:)); 
            plot3(-pi*ones(size(t_mean)),t_mean, ...
                smooth(celltmp1{ti,pi}.v(2,idx_aftrelease)-v_mean, 20), 'color', color_arr(dist_i+4,:)); 
            end
%             xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        view(120, 57);
        % plot notes here: 
%         linkaxes(axh, 'x'); xlim([-0.1 1]);
%         linkaxes(axh(2:end), 'y'); ylim([-0.1 0.1]); 
% %         ylim([0.47 0.60]);
%         subplot(4,1,1); title('F-pert'); ylabel('N');
%         subplot(4,1,2); title('\Delta velocity'); ylabel('m/s');
%         subplot(4,1,3); title('\Delta velocity'); ylabel('m/s');
%         subplot(4,1,4); title('\Delta velocity'); ylabel('m/s');
%         
        xlabel('pert time');
        ylabel('time');
        zlabel('velocity')
%         frame = getframe(gcf);
%         writeVideo(v,frame);
        view(58.2, 68.4);
    end
    
% end
% close(v);



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   TIDY UP DATA AND FORWARD TO JAMES                     %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Tidy up for the spring data 
ss_num = {  3814        3813        3803
            3815        3811        [3804 3805]
            3816        [3807 3808 3809 3810]        3806};     % 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 12;     % 1 without pert, and 12 perturbation time
data = cell(1, size(ss_num,1), size(ss_num,2), 5, pertT_num+1); % The stochastic ones are attached at the end 
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
            celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            
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
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1,frc_i,dist_i,1:5,1:pertT_num) = celltmp1(1:5,1:pertT_num);
        data(1,frc_i,dist_i,1:5,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3803_3812.mat', 'data'); % 12N perturbation, various time



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Tidy up for the subject data
ss_num = {  3825        [3826 3827]        3828
            [3819 3820] 3818               3821
            3824        3823               3822};     % 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 12;     % 1 without pert, and 12 perturbation time
data = cell(1, 1, size(ss_num,1), size(ss_num,2), 5, pertT_num); % NO stochastic ones
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
            celltmptmp = ss_tmp.export_as_formatted(1);
                % dim: ifsucess - #targets - #trials #pert
            % check the size of celltmptmp
            cell_avail_num = zeros(1,3);
            for i = 1:size(celltmptmp,3)            % # of trials
                for j = 1:size(celltmptmp,4)        % # of perturbations
                    if ~isempty(celltmptmp{1,1,i,j})
                        cell_avail_num(j) = cell_avail_num(j) + 1;
                    end
                end
            end
            
            for j = 1:size(celltmptmp,4)
            celltmp(cell_idx_from(j)+(1:cell_avail_num(j)),j) = ...
                celltmptmp(1,1,1:cell_avail_num(j),j);
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts >= 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1, 1,frc_i,dist_i,1:5,1:pertT_num) = celltmp1(1:5,1:pertT_num);
        %data(1,frc_i,dist_i,1:5,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3818_3828.mat', 'data'); % 12N perturbation, various time



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Visualize tidied up data to James                     %

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Visualize for the spring data 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3803_3812.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 12 + 1;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/SpringTimeChangingPert_3-by-3.mp4', 'MPEG-4');
v.FrameRate = 1;
% open(v);

fh(1) = figure();
hold on;
for pi = 1:(pertT_num-1)
    clf;
    for fce_i = 1:size(data,2)
       % 1. plot the perturbed force in the first panel 
       
       axh(1,fce_i) = subplot(4,3,fce_i); hold on;                     % plot PF
       celltmp1 = reshape(data(1,fce_i,1,:,:),5,pertT_num);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       % 2. plot the perturbed velocity in the following panels 
        for dist_i = 1:size(data,3) % for each spring 
            axh(dist_i+1,fce_i) = subplot(4,3,fce_i+(dist_i*3)); hold on;         % plot each response
            celltmp1 = reshape(data(1,fce_i,dist_i,:,:),5,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
                plot(t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(dist_i+1,:));
            end
        end
    end
        % plot notes here: 
        linkaxes(axh(:), 'x'); xlim([-0.1 1]);
        linkaxes(axh(2:end,:), 'y'); ylim([-1 1.2]); 
        for fce_i = 1:3
        subplot(4,3,fce_i); title('F-pert'); ylabel('N');
        subplot(4,3,fce_i*3+1); title('velocity'); ylabel('m/s');
        subplot(4,3,fce_i*3+2); title('velocity'); ylabel('m/s');
        subplot(4,3,fce_i*3+3); title('velocity'); ylabel('m/s');
        xlabel('time');
        end
        set(gcf, 'position', [0,0, 1080, 680]);
        frame = getframe(gcf);
%         writeVideo(v,frame);
end
% close(v);  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Visualize for the subject data 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3818_3828.mat', 'data');
pertT_num = 1 + 12;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
num_trials = 5;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/SubjectTimeChangingPert_3-by-3_100ms.mp4', 'MPEG-4');
color_arr = colormap('lines');
close all;
% v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/SubjectTimeChangingPert_3-by-3.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

fh(1) = figure();
hold on;
for pi = 1:(pertT_num)
    clf;
    for fce_i = 1:size(data,3)
       % 1. plot the perturbed force in the first panel 
       
       axh(1,fce_i) = subplot(4,3,fce_i); hold on;                     % plot PF
       celltmp1 = reshape(data(1,1,fce_i,1,:,:),num_trials,pertT_num);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       % 2. plot the perturbed velocity in the following panels 
        for dist_i = 1:size(data,4) % for each target 
            axh(dist_i+1,fce_i) = subplot(4,3,fce_i+(dist_i*3)); hold on;         % plot each response
            celltmp1 = reshape(data(1,1,fce_i,dist_i,:,:),5,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
            end
        end
    end
        % plot notes here: 
        linkaxes(axh(:), 'x'); xlim([-0.1 1]);
        linkaxes(axh(2:end,:), 'y'); ylim([-0.15 0.6]); 
        for fce_i = 1:3
        subplot(4,3,fce_i); title('F-pert'); ylabel('N');
        subplot(4,3,fce_i*3+1); title('velocity'); ylabel('m/s');
        subplot(4,3,fce_i*3+2); title('velocity'); ylabel('m/s');
        subplot(4,3,fce_i*3+3); title('velocity'); ylabel('m/s');
        xlabel('time');
        end
        set(gcf, 'position', [0,0, 1080, 680]);
        frame = getframe(gcf);
        writeVideo(v,frame);
end
close(v);


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The previous data with 100ms (end-to-end) pulse was said to be too small
% to analysis. Hence this time I widen the pulse to 200ms (end-to-end). To
% let the thing happens faster, I get the data only 5 pulse times (100,
% 300, 500, 700, 900ms)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   TIDY UP DATA AND FORWARD TO JAMES                     %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Tidy up for the spring data 
ss_num = {  3873,             3881,         3882       
            [3874,3875,3876], [3879,3880],  3883  
            3877,             3878,         3884       };     % 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 5;     % 1 without pert, and 5 perturbation time
data = cell(1, size(ss_num,1), size(ss_num,2), 5, pertT_num+1); % The stochastic ones are attached at the end 
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
            
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
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1,frc_i,dist_i,1:5,1:pertT_num) = celltmp1(1:5,1:pertT_num);
        data(1,frc_i,dist_i,1:5,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3873_3884.mat', 'data'); % 12N perturbation, various time



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Tidy up for the subject data
ss_num = {  3888        [3885 3886]        3887
            3889        3890               3891
            3893        3894               3892};     % 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 5;     % 1 without pert, and 12 perturbation time
data = cell(1, 1, size(ss_num,1), size(ss_num,2), 5, pertT_num); % NO stochastic ones
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
            celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
                % dim: ifsucess - #targets - #trials #pert
            % check the size of celltmptmp
            cell_avail_num = zeros(1,3);
            for i = 1:size(celltmptmp,1)            % # of trials
                for j = 1:size(celltmptmp,2)        % # of perturbations
                    if ~isempty(celltmptmp{i,j})
                        cell_avail_num(j) = cell_avail_num(j) + 1;
                    end
                end
            end
            
            for j = 1:size(celltmptmp,2)
            celltmp(cell_idx_from(j)+(1:cell_avail_num(j)),j) = ...
                celltmptmp(1:cell_avail_num(j),j);
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts >= 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1, 1,frc_i,dist_i,1:5,1:pertT_num) = celltmp1(1:5,1:pertT_num);
        %data(1,frc_i,dist_i,1:5,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3885_3894.mat', 'data'); % 12N perturbation, various time

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Visualize for the spring data 
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3925_3937.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 5 + 1;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
close all;
% v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/SpringTimeChangingPert_3-by-3_200ms.mp4', 'MPEG-4');
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3925_3937/v_compare.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

fh(1) = figure();
hold on;
for pi = 1:(pertT_num-1)
    clf;
    for fce_i = 1:size(data,2)
       % 1. plot the perturbed force in the first panel 
       trial_num = size(data,4); % 4 for spring data, 5 for human data
       axh(1,fce_i) = subplot(4,3,fce_i); hold on;                     % plot PF
       celltmp1 = reshape(data(1,fce_i,1,:,:),trial_num,pertT_num);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       % 2. plot the perturbed velocity in the following panels 
        for dist_i = 1:size(data,3) % for each spring 
            axh(dist_i+1,fce_i) = subplot(4,3,fce_i+(dist_i*3)); hold on;         % plot each response
            celltmp1 = reshape(data(1,fce_i,dist_i,:,:),trial_num,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot(t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
%                 plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(dist_i+1,:));
            end
        end
    end
        % plot notes here: 
        linkaxes(axh(:), 'x'); xlim([-0.1 1.5]);
        
%         title_str = 'force'; label_str = 'N';
%         linkaxes(axh(2:end,:), 'y'); ylim([-12 28]); %force 
        title_str = 'velocity'; label_str = 'm/s';
        linkaxes(axh(2:end,:), 'y'); ylim([-1 1.2]); % velocity
        
        for fce_i = 1:3
        subplot(4,3,fce_i); title('F-pert'); ylabel('N');
        subplot(4,3,fce_i*3+1); title(title_str); ylabel(label_str);
        subplot(4,3,fce_i*3+2); title(title_str); ylabel(label_str);
        subplot(4,3,fce_i*3+3); title(title_str); ylabel(label_str);
        xlabel('time');
        end
        set(gcf, 'position', [0,0, 1080, 680]);
        frame = getframe(gcf);
%         writeVideo(v,frame);
end
close(v);  

%% plot in 3D curve 
figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
color_arr = colormap('lines');
close all;
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring 
        subplot(3,3, (fce_i-1)*3+dist_i); hold on;
%         figure; hold on;
        for pi = 1:6%1:length(pertT_unq)
%        fh(pi) = figure(); hold on;
       
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
%        celltmp1 = reshape(data(1,fce_i,1,:,1:6),5,6); % 1 no -ert and 5 pert
        celltmp1 = reshape(data(1,fce_i,1,:,1:6),15,6); % 1 no -ert and 5 pert
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
        
%             axh(dist_i+1) = subplot(3,1,dist_i+1); hold on;         % plot each response
%             celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),5,6);
        celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),15,6);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.x(2,:), 'color', [0.5 0.5 0.5]);
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)); 
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)); 
                    % 4, for consistant with previous color 
            end
            %xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        title(['Force ' num2str(F_list(fce_i)) ' K ' num2str(K_list(dist_i))]);
        ylim([-0.1 2.0]);
        zlim([-1.5 1.5]); % velocity
%         zlim([0.47 0.72]);   % position
        xlabel('perturb positions');
        ylabel('release time');
        zlabel('endpoint velocity(m/s)');
%         view(120, 57);
%         set(gca, 'View', [90, 0]); % this view all lines are overlapped.  
%         set(gca, 'View', [57, 65]); % this view all pert are well aligned. 
        set(gca, 'View', [112, 47]); % according to James' figure. 
%         saveas(gcf, ['data/processedData/dataDescriptions/ss3873_3884/velF' num2str(F_list(fce_i)) 'K' num2str(K_list(dist_i)) '.png']);
%         close all;
    end
    
    
end

%% %%%%%%%%%%%%%%%%%% Andy's simplist idea on the dF/dx %%%%%%%%%%%%%%%%%%%

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
dist_i = 1; % the 640N/m springs
K_est_cell = cell(3,1);
K_est_all = [];
for frc_i = 1:3    
    for trial_i = 1:5
        dattmp = data{1,frc_i,dist_i,trial_i,6};
        dattmp.f(2,:) = smooth(dattmp.f(2,:), 20);
        
        valid_id = find(dattmp.ts >= 5 & ...    % at move
            dattmp.Fp(2,:) ~= 0);       % at pert 
        avg_idx = [valid_id(1) + [-100:0], valid_id(end)+[250:500]];
        f_avg = mean(dattmp.f(2,avg_idx));
        x_avg = mean(dattmp.x(2,avg_idx));
        
        clf;
        axh(1) = subplot(2,1,1); hold on;
        plot(dattmp.t, dattmp.f(2,:));
        plot(dattmp.t(avg_idx), dattmp.f(2,avg_idx), 'r.');
        plot(dattmp.t(avg_idx), ones(1,length(avg_idx))*f_avg);
        
        axh(2) = subplot(2,1,2);
        plot(dattmp.t, dattmp.x(2,:)); hold on;
        plot(dattmp.t(avg_idx), dattmp.x(2,avg_idx), 'r.');
        plot(dattmp.t(avg_idx), ones(1,length(avg_idx))*x_avg);
        linkaxes(axh, 'x');
        
        % this part seems not good
        % axh(3) = subplot(3,1,3);  hold on;
        % df = smooth(dattmp.f(2,valid_id)) - f_avg;
        % dx = dattmp.x(2,valid_id) - x_avg;
        % plot(dattmp.t(valid_id), df./dx, '.');
        
        % max over max
        vlid_id_long = valid_id(100) : valid_id(end) + 250;
        [pks, locs] = findpeaks(abs(dattmp.f(2,vlid_id_long) - f_avg), ...
            'MinPeakDistance', 50, ...
            'MinPeakHeight',1);
        f_dat = dattmp.f(2,vlid_id_long) - f_avg;
        df = f_dat(locs(1:2));  
        
        subplot(axh(1)); 
        plot(dattmp.t(vlid_id_long(locs(1:2))), f_dat(locs(1:2)) + f_avg, 'mo', ...
            'MarkerSize', 5);
        title('Force censored');
        
        [pks, locs] = findpeaks(abs(dattmp.x(2,vlid_id_long) - x_avg), ...
            'MinPeakDistance', 50, ...
            'MinPeakHeight',0.001);
        x_dat = dattmp.x(2,vlid_id_long) - x_avg;
        dx = x_dat(locs(1:2));
        
        subplot(axh(2)); 
        plot(dattmp.t(vlid_id_long(locs(1:2))), x_dat(locs(1:2)) + x_avg, 'mo', ...
            'MarkerSize', 5);
        title('position censored');
        
        linkaxes(axh, 'x'); 
        xlim([dattmp.t(vlid_id_long(1)) - 0.02, dattmp.t(vlid_id_long(end)) + 0.02]);
        K_est = df./dx;
        K_est
        K_est_all = [K_est_all, K_est];
        
    end
    K_est_cell{dist_i,1} = K_est_all;
end

%% the code just duplicate from last segment, and do roughly the same thing on different data
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
dist_i = 2; % the 320N/m springs
% K_est_cell = cell(3,3);
K_est_all = [];
for frc_i = 1:3         % I guess only force 1, 15 N sorks for 
    for trial_i = 1:5
        dattmp = data{1,frc_i,dist_i,trial_i,6};
        dattmp.f(2,:) = smooth(dattmp.f(2,:), 20);
        
        valid_id = find(dattmp.ts >= 5 & ...    % at move
            dattmp.Fp(2,:) ~= 0);       % at pert 
        avg_idx = [valid_id(end)+[250:500]];
        f_avg = mean(dattmp.f(2,avg_idx));
        x_avg = mean(dattmp.x(2,avg_idx));
        
        clf;
        axh(1) = subplot(2,1,1); hold on;
        plot(dattmp.t, dattmp.f(2,:));
        plot(dattmp.t(avg_idx), dattmp.f(2,avg_idx), 'r.');
        plot(dattmp.t(avg_idx), ones(1,length(avg_idx))*f_avg);
        
        axh(2) = subplot(2,1,2);
        plot(dattmp.t, dattmp.x(2,:)); hold on;
        plot(dattmp.t(avg_idx), dattmp.x(2,avg_idx), 'r.');
        plot(dattmp.t(avg_idx), ones(1,length(avg_idx))*x_avg);
        linkaxes(axh, 'x');
        
        % this part seems not good
        % axh(3) = subplot(3,1,3);  hold on;
        % df = smooth(dattmp.f(2,valid_id)) - f_avg;
        % dx = dattmp.x(2,valid_id) - x_avg;
        % plot(dattmp.t(valid_id), df./dx, '.');
        
        % max over max
        vlid_id_long = valid_id(end)+50 : valid_id(end) + 250;
        [pks, locs] = findpeaks(abs(dattmp.f(2,vlid_id_long) - f_avg), ...
            'MinPeakDistance', 50, ...
            'MinPeakHeight',1);
        f_dat = dattmp.f(2,vlid_id_long) - f_avg;
        df = f_dat(locs(1));  
        
        subplot(axh(1)); 
        plot(dattmp.t(vlid_id_long(locs(1))), f_dat(locs(1)) + f_avg, 'm.');
        title('Force censored');
        
        [pks, locs] = findpeaks(abs(dattmp.x(2,vlid_id_long) - x_avg), ...
            'MinPeakDistance', 50, ...
            'MinPeakHeight',0.001);
        x_dat = dattmp.x(2,vlid_id_long) - x_avg;
        dx = x_dat(locs(1));
        
        subplot(axh(2)); 
        plot(dattmp.t(vlid_id_long(locs(1))), x_dat(locs(1)) + x_avg, 'm.');
        title('position censored');
        
        linkaxes(axh, 'x'); 
        xlim([dattmp.t(vlid_id_long(1)) - 0.02, dattmp.t(vlid_id_long(end)) + 0.02]);
        K_est = df./dx;
        K_est
        K_est_all = [K_est_all, K_est];
        
    end
     K_est_cell{dist_i,1} = K_est_all;
end

%% simple show the 'sanityCheck stiffness' and the 'real stiffness'  
figure(); hold on;
% have a bar plot here 
K_est_all = -K_est_cell{1,1}; 
scatter(1*ones(size(K_est_all)), K_est_all, 5);
K_est_mean1 = mean(K_est_all); 
K_est_std1 = std(K_est_all);

K_est_all = -K_est_cell{2,1};
scatter(2*ones(size(K_est_all)), K_est_all, 5);
K_est_mean2 = mean(K_est_all); 
K_est_std2 = std(K_est_all);

yline(640);
yline(320);

xlim([0.5, 2.5]);
ylim([0 700]);
xticks([1,2]);
xticklabels({'K640', 'K320'});
title('stiffness of dF/dx');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Visualize for the subject data 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3885_3894.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 5;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/SubjectTimeChangingPert_3-by-3.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

fh(1) = figure();
hold on;
for pi = 1:(pertT_num)
    clf;
    for fce_i = 1:size(data,3)
       % 1. plot the perturbed force in the first panel 
       
       axh(1,fce_i) = subplot(4,3,fce_i); hold on;                     % plot PF
       celltmp1 = reshape(data(1,1,fce_i,1,:,:),5,pertT_num);
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
       plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       % 2. plot the perturbed velocity in the following panels 
        for dist_i = 1:size(data,4) % for each target 
            axh(dist_i+1,fce_i) = subplot(4,3,fce_i+(dist_i*3)); hold on;         % plot each response
            celltmp1 = reshape(data(1,1,fce_i,dist_i,:,:),5,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+1,:));
            end
        end
    end
        % plot notes here: 
        linkaxes(axh(:), 'x'); xlim([-0.1 1.5]);
        linkaxes(axh(2:end,:), 'y'); ylim([-0.15 0.6]); 
        for fce_i = 1:3
        subplot(4,3,fce_i); title('F-pert'); ylabel('N');
        subplot(4,3,fce_i*3+1); title('velocity'); ylabel('m/s');
        subplot(4,3,fce_i*3+2); title('velocity'); ylabel('m/s');
        subplot(4,3,fce_i*3+3); title('velocity'); ylabel('m/s');
        xlabel('time');
        end
        set(gcf, 'position', [0,0, 1080, 680]);
        frame = getframe(gcf);
        writeVideo(v,frame);
end
close(v);

%% plot in 3D curve 
figure(); 
pt_num = 1 + 5; % 5 perts
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3896_3905.mat', 'data'); pt_num = 1+7; % 7 perturbs
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
for fce_i = 1:size(data,3)
    for dist_i = 1:size(data,4) % for each spring 
        subplot(3,3, (fce_i-1)*3+dist_i); hold on;
% %         figure; hold on;
        for pi = 1:pt_num%1:length(pertT_unq)
%        fh(pi) = figure(); hold on;
       
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
%        celltmp1 = reshape(data(1,1,fce_i,1,:,1:6),5,6); % 1 no -ert and 5 pert
        celltmp1 = reshape(data(1,1,fce_i,1,:,1:8),10,8); % 1 no -ert and 7 pert
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
        
%             axh(dist_i+1) = subplot(3,1,dist_i+1); hold on;         % plot each response
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:6),5,6);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.x(2,:), 'color', [0.5 0.5 0.5]);
%                   plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)); 
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)); 
%                   plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)); 
                    % 4, for consistant with previous color 
            end
            %xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        title(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i))]);
%         ylim([-0.1 1.36]);
        ylim([-0.1 0.5]);
        zlim([-0.2 0.55]); % velocity
%         zlim([0.47 0.58]);   % position
        xlabel('perturb positions');
        ylabel('release time');
        zlabel('endpoint velocity(m/s)');
%         view(120, 57);
%         set(gca, 'View', [90, 0]); % this view all lines are overlapped.
%         set(gca, 'View', [57, 65]); % this view all pert are well aligned. 
        set(gca, 'View', [112, 47]); % according to James' figure. 
%         saveas(gcf, ['data/processedData/dataDescriptions/ss3896_3905/velF' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) '.png']);
%         close all;
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part show the the stiffness variate during movement (for subject). 
% Untimate goal is to see the df/dt throughout moving. (even on a single
% case).  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see what is the x and f change during the movement  
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3818_3828.mat', 'data');
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3913_3921.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1:size(data,3)
    for dist_i = 1:size(data,4) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
        psudoK_mat = zeros(5,7); % 
        for pi = 1:8%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:13),5,13); % 1 no -ert and 5 pert
            celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:8),10,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            f_avg = f_avg/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(3,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(3,2,3);
            plot(t_grids, x_avg);
            axh(5) = subplot(3,2,5);
            plot(t_grids, f_avg);
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(3,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(3,2,4); hold on;% subtracted x
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(3,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(3,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(3,2,4); % subtracted x
                plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(3,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(3,2,4); % subtracted x
                scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
                axh(6) = subplot(3,2,6); % subtracted F
                scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
                
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(3,2,4); % subtracted x
                scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
                axh(6) = subplot(3,2,6); % subtracted F
                scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
        close all;
    end
    
end


 
%% plot psudoK_cell
figure('unit', 'inch', 'position', [0,0,12,12]); 
t_ppeak = 50:50:600;
for fce_i = 1:3
    for dist_i = 1:3
        psudoK_mat = psudoK_cell{fce_i,dist_i};
        axh((fce_i-1)*3 + dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i); 
        var_x = repmat(t_ppeak,5,1); var_y = psudoK_mat;
        scatter(var_x(:), var_y(:), 20, 'MarkerFaceColor', 'b');
        title([num2str(F_list(fce_i)) 'N ' num2str(K_list(dist_i)) 'cm ']);
        ylabel('?N/m?');
        xlabel('pulse peak (ms)');
        grid on;
    end
end
sgtitle('psudo stiffness after release');

linkaxes(axh, 'y');
ylim(axh(1), [-100, 10000]);
saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/PsudoStiffnessDuringMovement.png']);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do more pulses during movement as Andy request 
% 2. Tidy up for the subject data
% clear; clc; 
ss_num = {  3898    3896    3897
            3899    3900    3901
            3904    3905    [3902 3903]};     % 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 7;     % 1 without pert, and 7 perturbation time
data = cell(1, 1, size(ss_num,1), size(ss_num,2), 5, pertT_num); % NO stochastic ones
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
            celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
                % dim: ifsucess - #targets - #trials #pert
            % check the size of celltmptmp
            cell_avail_num = zeros(1,3);
            for i = 1:size(celltmptmp,1)            % # of trials
                for j = 1:size(celltmptmp,2)        % # of perturbations
                    if ~isempty(celltmptmp{i,j})
                        cell_avail_num(j) = cell_avail_num(j) + 1;
                    end
                end
            end
            
            for j = 1:size(celltmptmp,2)
            celltmp(cell_idx_from(j)+(1:cell_avail_num(j)),j) = ...
                celltmptmp(1:cell_avail_num(j),j);
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts >= 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/12.5)*12.5;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1, 1,frc_i,dist_i,1:10,1:pertT_num) = celltmp1(1:10,1:pertT_num);
        %data(1,frc_i,dist_i,1:5,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3896_3905.mat', 'data'); % 12N perturbation, various time


%% plot the psudoStiffness result 

load('data/processedData/ss3896_3905.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
%F_list = [15];
%K_list = [5.0];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1:size(data,3)
    for dist_i = 1:size(data,4) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
        psudoK_mat = zeros(10,7); % 
        for pi = 1:8%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
            celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:8),10,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            f_avg = f_avg/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(3,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(3,2,3);
            plot(t_grids, x_avg);
            axh(5) = subplot(3,2,5);
            plot(t_grids, f_avg);
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(3,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(3,2,4); hold on;% subtracted x
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(3,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(3,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(3,2,4); % subtracted x
                plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(3,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(3,2,4); % subtracted x
                scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
                axh(6) = subplot(3,2,6); % subtracted F
                scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
                
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(3,2,4); % subtracted x
                scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
                axh(6) = subplot(3,2,6); % subtracted F
                scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
%         saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
        close all;
    end
    
end

% t_pert  = 50:25:200;
% figure(); plot(t_pert, -psudoK_mat', 'o'); 
% xlim([-100 500]);  
% ylabel('? N/m?');  
% 
% ylim([0 10000])

%
figure('unit', 'inch', 'position', [0,0,12,12]); 
t_ppeak = 50:25:200;
for fce_i = 1:3
    for dist_i = 1:3
        psudoK_mat = psudoK_cell{fce_i,dist_i};
        axh((fce_i-1)*3 + dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i); 
        var_x = repmat(t_ppeak,10,1); var_y = psudoK_mat;
        scatter(var_x(:), var_y(:), 20, 'MarkerFaceColor', 'b');
        title([num2str(F_list(fce_i)) 'N ' num2str(K_list(dist_i)) 'cm ']);
        ylabel('?N/m?');
        xlabel('pulse peak (ms)');
        grid on;
    end
end
sgtitle('psudo stiffness after release');

linkaxes(axh, 'y');
ylim(axh(1), [-100, 10000]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% test if the data works fine
ss_num = [3919    3920    3921
            3918    3917    3916
            3914    3913    3915];
ss_num = ss_num(:);
for i = 1:9
    ss_num(i)
    sstmp(i) = SessionScan(ss_num(i));
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do more pulses during movement?100ms:25ms:250ms
% 2. Tidy up for the subject data
% clear; clc; 
ss_num = {  3919    3920    3921
            3918    3917    3916
            3914    3913    3915};     % 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 7;     % 1 without pert, and 7 perturbation time
data = cell(1, 1, size(ss_num,1), size(ss_num,2), 5, pertT_num); % NO stochastic ones
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
                % dim: ifsucess - #targets - #trials #pert
            % check the size of celltmptmp
            cell_avail_num = zeros(1,3);
            for i = 1:size(celltmptmp,1)            % # of trials
                for j = 1:size(celltmptmp,2)        % # of perturbations
                    if ~isempty(celltmptmp{i,j})
                        cell_avail_num(j) = cell_avail_num(j) + 1;
                    end
                end
            end
            
            for j = 1:size(celltmptmp,2)
            celltmp(cell_idx_from(j)+(1:cell_avail_num(j)),j) = ...
                celltmptmp(1:cell_avail_num(j),j);
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts >= 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/12.5)*12.5;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1, 1,frc_i,dist_i,1:10,1:pertT_num) = celltmp1(1:10,1:pertT_num);
        
        %data(1,frc_i,dist_i,1:5,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3913_3921.mat', 'data'); % 12N perturbation, various time

% % check if all have 10
% for fce_i = 1:3
%     for dist_i = 1:3
%         pert_cell = reshape(data(1,1,fce_i,dist_i,1:10,1:pertT_num), 10, pertT_num);
%         disp(['ss ' num2str(ss_num{fce_i,dist_i})]);
%         pert_cell
%     end
% end

% show as plotted

%% plot in 3D curve 
figure(); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data'); pt_num = size(data,5); % 7 perturbs
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
% close all;
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for target 
for fce_i = 1:size(data,2)
     for dist_i = 1:size(data,3) % for each spring 
        subplot(3,3, (fce_i-1)*3+dist_i); hold on;
%         figure; hold on;
        for pi = 1:pt_num%1:length(pertT_unq)
%        fh(pi) = figure(); hold on;
       
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
%        celltmp1 = reshape(data(1,1,fce_i,1,:,1:6),5,6); % 1 no -ert and 5 pert
%         celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:8),10,8); % 1 no -ert and 7 pert
         celltmp1 = reshape(data(1,fce_i,dist_i,:,:),size(data,[4 5])); % for sprigns
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
        
%             axh(dist_i+1) = subplot(3,1,dist_i+1); hold on;         % plot each response
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:6),5,6);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.x(2,:), 'color', [0.5 0.5 0.5]);
                  plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)); 
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)); 
                  plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)); 
                    % 4, for consistant with previous color 
            end
            %xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        title(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i))]);
        ylim([-0.1 1.36]);
%         ylim([-0.1 0.5]);
%         zlim([-0.2 0.55]); % velocity
%         zlim([0.47 0.58]);   % position
        xlabel('perturb positions');
        ylabel('release time');
        zlabel('endpoint velocity(m/s)');
%         view(120, 57);
         set(gca, 'View', [90, 0]); % this view all lines are overlapped.
%         set(gca, 'View', [57, 65]); % this view all pert are well aligned. 
%         set(gca, 'View', [112, 47]); % according to James' figure. 
%         saveas(gcf, ['data/processedData/dataDescriptions/ss3913_3921/velF' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) '.png']);
%         close all;
    end 
end
% saveas(gcf, ['data/processedData/dataDescriptions/ss3913_3921/velAll_flat.png']);


%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
% Sanity Check: does the right force being generated?
%%%%%%%%%%%%%%%%%%%%
%% cpDatarg2([3906 3907]);
sstmp1 = SessionScan(3908); 
sstmp2 = SessionScan(3907);
%%

t_idx = find(sstmp1.data.Fp(2,:)~=0 & ...
    (sstmp1.data.ts==4 | sstmp1.data.ts==5));
t_shift = sstmp1.data.t(t_idx(1)); 
sstmp1.data.t_shift = sstmp1.data.t - t_shift;

figure; 
axh(1) = subplot(2,1,1);  grid on; hold on;
plot(sstmp1.data.t_shift, -sstmp1.data.Fp, 'Marker','.'); 
% axh(3) = subplot(4,1,2); grid on;
plot(sstmp1.data.t_shift, sstmp1.data.f, 'Marker','.');
ylim([0 20]);
% linkaxes(axh([1,3]), 'x')  
title('pulse of 100ms');
legend('command x', 'command y', 'command z', 'sensor x', 'sensor x', 'sensor z');


t_idx = find(sstmp2.data.Fp(2,:)~=0 & ...
    (sstmp2.data.ts==4 | sstmp2.data.ts==5));
t_shift = sstmp2.data.t(t_idx(1)); 
xlim([t_shift, t_shift + 0.5]);
sstmp2.data.t_shift = sstmp2.data.t - t_shift;


axh(2) = subplot(2,1,2); grid on; hold on;
plot(sstmp2.data.t_shift, -sstmp2.data.Fp, 'Marker','.'); 
% axh(4) = subplot(4,1,4); grid on;
plot(sstmp2.data.t_shift, sstmp2.data.f, 'Marker','.'); 
ylim([0 20]);
% linkaxes(axh([2,4]), 'x')

linkaxes(axh, 'x'); 
xlim([0 0.5]);
linkaxes(axh, 'y');
ylim([0 18]);

title('pulse of 200ms');
legend('command x', 'command y', 'command z', 'sensor x', 'sensor x', 'sensor z');

%% plot by trials 
% sstmp = sstmp1;  
% 3910, 3909, 3912, 3911
ss_list = [3910, 3909, 3912, 3911];
dur_name = {'100ms', '200ms', '300ms', '400ms'};

for session_i = 1:4
    sstmp(session_i) = SessionScan(ss_list(session_i));
end

%%
% The force exertion and censored during hold 
fh(1) = figure();
color_arr = colormap('lines');
for session_i = 1:4
    sstmp1 = sstmp(session_i); 
    
for trial_i = 1:min(20,length(sstmp1.trials))
    t_idx = find(sstmp1.trials(trial_i).data.Fp(2,:)~=0 & ...
        (sstmp1.trials(trial_i).data.ts==4)); % hold 
    if isempty(t_idx) 
        continue;
    end
    if (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 6)
        axh((session_i-1)*2 + 1) = subplot(4,2,(session_i-1)*2 + 1); % left one;
        hold on;
    elseif (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 12)
        axh((session_i-1)*2 + 2) = subplot(4,2,(session_i-1)*2 + 2); % right one;
        hold on;
    end
    
    t_shift = sstmp1.trials(trial_i).data.t(t_idx(1));
    sstmp1.trials(trial_i).data.t_shift = sstmp1.trials(trial_i).data.t - t_shift;
    
%     axh(1) = subplot(2,1,1);  grid on; hold on;
%     plot(sstmp1.trials(trial_i).data.t_shift, -sstmp1.trials(trial_i).data.Fp(2,:), 'Marker','.', 'Color', color_arr(1,:));
    plot(sstmp1.trials(trial_i).data.t_shift, -sstmp1.trials(trial_i).data.Fp(2,:),  'Color', color_arr(2,:), 'LineWidth', 2);
    % axh(3) = subplot(4,1,2); grid on;
    plot(sstmp1.trials(trial_i).data.t_shift, sstmp1.trials(trial_i).data.f(2,:), 'Marker','.', 'Color', color_arr(3,:));
    ylim([0 15]);
    % linkaxes(axh([1,3]), 'x')
%     title('pulse of 100ms');
    if((session_i-1)*2 + 1 == 1)
        legend('command y', 'sensor y');
    end
    if ((session_i-1) == 3)
        xlabel('t (s)' );
    end
    ylabel('force (N)');
end 
    subplot(4,2,(session_i-1)*2 + 1); grid on; title(['6N max' dur_name{session_i}]);
    subplot(4,2,(session_i-1)*2 + 2); grid on; title(['12N max' dur_name{session_i}]);
end
linkaxes(axh, 'x');
xlim([0 0.4]); % perturb holding 
sgtitle('Force during hold');
% sgtitle('during move');

%% The force exertion and censored during hold and release 
clear axh lnh
fh(2) = figure('Position', [0 0 600 800]);
for ts = 4:5
for session_i = 1:4
    sstmp1 = sstmp(session_i); 
    
for trial_i = 1:min(20,length(sstmp1.trials))
    t_idx = find(sstmp1.trials(trial_i).data.Fp(2,:)~=0 & ...
        (sstmp1.trials(trial_i).data.ts==ts)); % release
    if isempty(t_idx) 
        continue;
    end
    if (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 6)
        axh((session_i-1)*2 + 1) = subplot(4,2,(session_i-1)*2 + 1); % left one;
        hold on;
    elseif (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 12)
        axh((session_i-1)*2 + 2) = subplot(4,2,(session_i-1)*2 + 2); % right one;
        hold on;
    end
    
    t_shift = sstmp1.trials(trial_i).data.t(t_idx(1));
    sstmp1.trials(trial_i).data.t_shift = sstmp1.trials(trial_i).data.t - t_shift;
    
%     axh(1) = subplot(2,1,1);  grid on; hold on;
%     plot(sstmp1.trials(trial_i).data.t_shift, -sstmp1.trials(trial_i).data.Fp(2,:), 'Marker','.', 'Color', color_arr(1,:));
    lnh{session_i}(1) = plot(sstmp1.trials(trial_i).data.t_shift, -sstmp1.trials(trial_i).data.Fp(2,:),  'Color', color_arr(2,:), 'LineWidth', 2);
    % axh(3) = subplot(4,1,2); grid on;
    switch ts
        case 4
    lnh{session_i}(2) = plot(sstmp1.trials(trial_i).data.t_shift, sstmp1.trials(trial_i).data.f(2,:), 'Marker','.', 'Color', color_arr(ts-1,:));
        case 5
            lnh{session_i}(3) = plot(sstmp1.trials(trial_i).data.t_shift, sstmp1.trials(trial_i).data.f(2,:), 'Marker','.', 'Color', color_arr(ts-1,:));
    end
    ylim([0 15]);
    % linkaxes(axh([1,3]), 'x')
%     title('pulse of 100ms');
%     if((session_i-1)*2 + 1 == 1)
%         legend('command y', 'sensor y');
%     end
    if ((session_i-1) == 3)
        xlabel('t (s)' );
    end
    ylabel('force (N)');
end 
    subplot(4,2,(session_i-1)*2 + 1); grid on; title(['6N max' dur_name{session_i}]);
    subplot(4,2,(session_i-1)*2 + 2); grid on; title(['12N max' dur_name{session_i}]);
end
linkaxes(axh, 'x');
xlim(axh(1), [0 0.4]); % perturb during release 
end
legend(lnh{1}, {'command Force', 'censored during hold', 'censored during release'});
% sgtitle('during hold');
sgtitle(fh(2), 'Force command and censored in gaussian pulse');

%% %% The position measurement during hold 
fh(4) = figure('Position', [0 0 600 800]);
color_arr = colormap('lines');

for ts = 4:5 % hold and release
for session_i = 1:4
    sstmp1 = sstmp(session_i); 
    
for trial_i = 2:min(20,length(sstmp1.trials))
    t_idx = find(sstmp1.trials(trial_i).data.Fp(2,:)~=0 & ...
        (sstmp1.trials(trial_i).data.ts==ts)); % hold 
    if isempty(t_idx) 
        continue;
    end
    if (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 6)
        axh((session_i-1)*2 + 1) = subplot(4,2,(session_i-1)*2 + 1); % left one;
        hold on;
    elseif (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 12)
        axh((session_i-1)*2 + 2) = subplot(4,2,(session_i-1)*2 + 2); % right one;
        hold on;
    end
    
    t_shift = sstmp1.trials(trial_i).data.t(t_idx(1));
    sstmp1.trials(trial_i).data.t_shift = sstmp1.trials(trial_i).data.t - t_shift;
    
    switch ts
        case 4
            lnh{session_i}(1) = plot(sstmp1.trials(trial_i).data.t_shift, sstmp1.trials(trial_i).data.x(2,:),  'Color', color_arr(ts-1,:), 'LineWidth', 2);
        case 5
            lnh{session_i}(2) = plot(sstmp1.trials(trial_i).data.t_shift, sstmp1.trials(trial_i).data.x(2,:),  'Color', color_arr(ts-1,:), 'LineWidth', 2);
    end
    ylim([0.483, 0.485]);
    if ((session_i-1) == 3)
        xlabel('t (s)' );
    end
    ylabel('position (m)');
end 
    subplot(4,2,(session_i-1)*2 + 1); grid on; title(['6N max' dur_name{session_i}]);
    subplot(4,2,(session_i-1)*2 + 2); grid on; title(['12N max' dur_name{session_i}]);

end
linkaxes(axh, 'x');
xlim([0 0.4]); % perturb holding 
sgtitle('position measurement in clamp experiment');
% sgtitle('during move');
end

    
legend(lnh{1}, {'during hold', 'during release'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   TIDY UP DATA AND FORWARD TO JAMES                     %
ss_num = [3937, 3931, 3925, ...
    3935,3936, 3930, 3928,3929, ...
    3934, 3932,3933, 3926,3927 ];
not_working_list = [];
for si = 1:length(ss_num)
    try
        SessionScan(ss_num(si));
    catch
        disp(['NOT WORK ss' num2str(ss_num(si))]);
        not_working_list = [not_working_list ss_num(si)];
    end
end
    
not_working_list
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Tidy up for the spring data 
        
ss_num = {  3937,           3931,        3925 
            [3935,3936],    3930,        [3928,3929]
            3934,           [3932,3933], [3926,3927]} 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 5;     % 1 without pert, and 5 perturbation time
data = cell(1, size(ss_num,1), size(ss_num,2), 5, pertT_num+1); % The stochastic ones are attached at the end 
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
            
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
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                idx_PF_peak = floor(idx_PF_peak/25)*25;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1,frc_i,dist_i,1:15,1:pertT_num) = celltmp1(1:15,1:pertT_num);
        data(1,frc_i,dist_i,1:15,pertT_num+1) = celltmp(1:15,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3925_3937.mat', 'data'); % 12N perturbation, various time
%% 
figure(); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3925_3937.mat', 'data'); pt_num = size(data,5); % 7 perturbs
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
% close all;
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for target 
for fce_i = 1:size(data,2)
     for dist_i = 1:size(data,3) % for each spring 
        subplot(3,3, (fce_i-1)*3+dist_i); hold on;
%         figure; hold on;
        for pi = 1:pt_num%1:length(pertT_unq)
%        fh(pi) = figure(); hold on;
       
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
%        celltmp1 = reshape(data(1,1,fce_i,1,:,1:6),5,6); % 1 no -ert and 5 pert
%         celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:8),10,8); % 1 no -ert and 7 pert
         celltmp1 = reshape(data(1,fce_i,dist_i,:,:),size(data,[4 5])); % for sprigns
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
        
%             axh(dist_i+1) = subplot(3,1,dist_i+1); hold on;         % plot each response
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:6),5,6);
            % plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.x(2,:), 'color', [0.5 0.5 0.5]);
                  plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
%                   plot3(-pi*ones(size(t)),t, celltmp1{ti,1}.Fp(2,:), 'color', [0.5 0.5 0.5]);
            end
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,pi+1})
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
%                 idx_release = find(celltmp1{ti,pi+1}.ts == 5);
%                 t = celltmp1{ti,pi+1}.t - celltmp1{ti,pi+1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi+1}.v(2,:), 'color', color_arr(dist_i+1,:));
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)); 
%                 plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)); 
%                   plot3(-pi*ones(size(t)),t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)); 
                  plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)); 
                    % 4, for consistant with previous color 
            end
            %xlim([-0.1 1]); ylim([-0.8 1.0]);
        end
        
        % plot notes here: 
        title(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i))]);
%         ylim([-0.1 1.36]);
          xlim([-0.1 1.36]);
%         ylim([-0.1 0.5]);
%         zlim([-0.2 0.55]); % velocity
%         zlim([0.47 0.58]);   % position
%         xlabel('perturb positions');
%         ylabel('release time');
%         zlabel('endpoint velocity(m/s)');
%         view(120, 57);
%          set(gca, 'View', [90, 0]); % this view all lines are overlapped.
%         set(gca, 'View', [57, 65]); % this view all pert are well aligned. 
%         set(gca, 'View', [112, 47]); % according to James' figure. 
%         saveas(gcf, ['data/processedData/dataDescriptions/ss3913_3921/velF' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) '.png']);
%         close all;
    end 
end
% saveas(gcf, ['data/processedData/dataDescriptions/ss3913_3921/velAll_flat.png']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. SPRING DATA COMPARE Visualize for the spring data 

fce_list = [15 20 25];
dist_list = [2.5 5.0 7.5];
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3925_3937.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 5 + 1;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3925_3937/f_compare.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

fh(1) = figure();
hold on;
for pi = 1:(pertT_num-1)
    clf;
    for dist_i = 1:size(data,3) % for each spring
        % 1. plot the perturbed force in the first panel
        trial_num = size(data,4); % 4 for spring data, 5 for human data
        axh(1,dist_i) = subplot(4,3,dist_i); hold on;                     % plot PF
        celltmp1 = reshape(data(1,1,dist_i,:,:),trial_num,pertT_num);
        idx_release = find(celltmp1{1,pi}.ts == 5);
        t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
        ylabel('command force')
        % 2. plot the perturbed velocity in the following panels
        for fce_i = 1:size(data,2)
            axh(dist_i+1,fce_i) = subplot(4,3,dist_i+(fce_i*3)); hold on;         % plot each response
            celltmp1 = reshape(data(1,fce_i,dist_i,:,:),trial_num,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
%                 plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot(t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
                 plot(t, celltmp1{ti,1}.x(2,:), 'color', [0.5 0.5 0.5]);
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+4,:));
%                 plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(dist_i+4,:));
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(dist_i+4,:));
                plot(t, -celltmp1{ti,pi}.Fp(2,:));
            end
        end
    end
        % plot notes here: 
        linkaxes(axh(:), 'x'); xlim([-0.1 1.5]);
        
        sgtitle_str = 'unperturbed and perturbed force'; label_str = 'censored force (N)';
%         linkaxes(axh(2:end,:), 'y'); ylim([-12 28]); %force 
        linkaxes(axh(2:end,:), 'y'); ylim([0.45 0.8]); %position 
%         sgtitle_str = 'unperturbed and perturbed velocity'; label_str = 'velocity (m/s)';
%         linkaxes(axh(2:end,:), 'y'); ylim([-1 1.2]); % velocity
%         
        for fce_i = 1:3
            for dist_i = 1:3
            subplot(4,3,fce_i*3 + dist_i); 
            title_str = (['force' num2str(fce_list(fce_i)) 'N dist' num2str(dist_list(dist_i)) 'cm']); 
            ylabel('N');
            subplot(4,3,fce_i*3 + dist_i); title(title_str); ylabel(label_str);
            xlabel('time');
            end
        end
        set(gcf, 'position', [0,0, 1080, 680]);
        sgtitle(sgtitle_str);
        fname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3925_3937/f_comparep' num2str(pi) '.png'];
%         saveas(gcf, fname);
        frame = getframe(gcf);
%         writeVideo(v,frame);
end
close(v);  

%% Tricks to get tidier data 
% After dumpped, some condition (f1d1p1, f3d1p1, do not have enough 15
% trials), I will copy and paste then...
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3925_3937.mat', 'data');
data(1,1,1,15,1) = data(1,1,1,6,1);
data(1,3,1,15,1) = data(1,3,1,6,1);
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3925_3937.mat', 'data');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subject data: perturbation during movement, randomize trials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpDatarg2(3923);
        
ss_num = {  [3938,3939],    3931,        3925 
            [3935,3936],    3930,        [3928,3929]
            3934,           [3932,3933], [3926,3927]} 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 7;     % 1 without pert, and 7 perturbation time
data = cell(1, size(ss_num,1), size(ss_num,2), 5, pertT_num); % The stochastic ones are not attached at the end 
for frc_i = 1%:size(ss_num,1) % actually force 
    for dist_i = 1%:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
            
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
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5 | celltmp_varT{ti}.ts == 6;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                ifplot = 0;
                if (ifplot)
                    plot(celltmp_varT{ti}.Fp(2,idx_ts5));
                end
                 idx_PF_peak = floor(idx_PF_peak/12.5)*12.5;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1,frc_i,dist_i,1:8,1:pertT_num) = celltmp1(1:8,1:pertT_num);
%         data(1,frc_i,dist_i,1:15,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3938_3939.mat', 'data'); % 12N perturbation, various time

%%                   TIDY UP DATA AND FORWARD TO JAMES                     %
% subject
ss_num = [3938, 3939, 3944, 3945, 3949, ...
    3940,3943, 3947, 3948, ...
    3941, 3942,3946];
not_working_list = [];
for si = 1:length(ss_num)
    try
        SessionScan(ss_num(si));
    catch
        disp(['NOT WORK ss' num2str(ss_num(si))]);
        not_working_list = [not_working_list ss_num(si)];
    end
end
    
not_working_list
%%
for si = 1:length(ss_num)
    sstmp = SessionScan(ss_num(si));
    idx = sstmp.getDelayedTrialIdx;
    disp(['ss' num2str(ss_num(si)) 'delay' num2str(idx)]);
end

    %%
ss_num = {  [3938,3939],    [3944,3945],        3949 
            3940,           3943,        [3947,3948]
            3941,           3942,       3946} 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 7;     % 1 without pert, and 7 perturbation time
data = cell(1, size(ss_num,1), size(ss_num,2), 5, pertT_num); % The stochastic ones are not attached at the end 
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
            
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
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5 | celltmp_varT{ti}.ts == 6;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                ifplot = 0;
                if (ifplot)
                    plot(celltmp_varT{ti}.Fp(2,idx_ts5));
                end
                 idx_PF_peak = floor(idx_PF_peak/12.5)*12.5;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1,frc_i,dist_i,1:7,1:pertT_num) = celltmp1(1:7,1:pertT_num);
%         data(1,frc_i,dist_i,1:15,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3938_3949.mat', 'data'); % 12N perturbation, various time


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The data with optotrak recording, still pulse during movement! 
% %                   TIDY UP DATA AND FORWARD TO JAMES                     %
% subject
ss_num = [3992,3993,    3988,        3987, ... 
            3994,           3989,        3997, ...
          3995,3996,    3990,3991,   3998, 3999];
not_working_list = [];
for si = 1:length(ss_num)
    try
        SessionScan(ss_num(si));
    catch
        disp(['NOT WORK ss' num2str(ss_num(si))]);
        not_working_list = [not_working_list ss_num(si)];
    end
end
    
not_working_list
%%
for si = 1:length(ss_num)
    sstmp = SessionScan(ss_num(si));
    idx = sstmp.getDelayedTrialIdx;
    disp(['ss' num2str(ss_num(si)) 'delay' num2str(idx)]);
end

    %% export data when there are optotrak message recorded. 
clear; clc; close all; 
ss_num = {  [3992,3993],    3988,        3987 
            3994,           3989,        3997
            [3995,3996],    [3990,3991],   [3998,3999]};
pert_f = 12; % only use 12N pert
pertT_num = 1 + 5;     % 1 without pert, and 5 perturbation time
data = cell(1, size(ss_num,1), size(ss_num,2), 10, pertT_num); % The stochastic ones are not attached at the end 
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % stiffness levels
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
            
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
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5 | celltmp_varT{ti}.ts == 6;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                ifplot = 0;
                if (ifplot)
                    plot(celltmp_varT{ti}.Fp(2,idx_ts5));
                end
                 idx_PF_peak = floor(idx_PF_peak/12.5)*12.5;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1,frc_i,dist_i,1:7,1:pertT_num) = celltmp1(1:7,1:pertT_num);
%         data(1,frc_i,dist_i,1:15,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss3987_3999.mat', 'data'); % 12N perturbation, various time

%% get the data where no release exist ( only perturb)
clear; clc; close all; 
ss_num = {  [4000, 4001, 4002]};
pert_f = 12; % only use 12N pert
pertT_num = 1 + 5;     % 1 without pert, and 5 perturbation time
data = cell(1, size(ss_num,1), size(ss_num,2), 10, pertT_num); % The stochastic ones are not attached at the end 
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % stiffness levels
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_from = [0 0 0];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
            
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
            end
            cell_idx_from = cell_idx_from + cell_avail_num;
            % save data in celltmp;
            
        end
        
        
        % detect the pulse time after cell tmp
            % 
        celltmp_varT = celltmp(:,2);
        clear pertT
        for ti = 1:length(celltmp_varT)
            if isempty(celltmp_varT{ti})
                continue;
            end
            pertT(ti) = 0;
            idx_PF_peak = nan;
            ifplot = 1;
            if (ifplot) 
                clf;
                axh(1) = subplot(2,1,1); 
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.ts);
                axh(2) = subplot(2,1,2);
                plot(celltmp_varT{ti}.t, celltmp_varT{ti}.Fp);
            end
            if max(abs(celltmp_varT{ti}.Fp(2,:))) ~= 0
                idx_ts5 = celltmp_varT{ti}.ts == 5 | celltmp_varT{ti}.ts == 6;
                idx_PF_peak = find([abs(celltmp_varT{ti}.Fp(2,idx_ts5)) == max(abs(celltmp_varT{ti}.Fp(2,idx_ts5)))]);
                ifplot = 0;
                if (ifplot)
                    plot(celltmp_varT{ti}.Fp(2,idx_ts5));
                end
                 idx_PF_peak = floor(idx_PF_peak/12.5)*12.5;
            end
            pertT(ti) = idx_PF_peak * 0.002; % 500Hz
        end
        % classify the pulse time into cells that have different time
        
        % get the index of each delay interval 
        [pertT_unq, ia, ic] = unique(pertT);
        pertT_unq
        idx_trialsPertT = cell(1, length(pertT_unq));
        trials_num_max = 0;
        for pi = 1:length(pertT_unq) 
            idx_trialsPertT{pi} = find(pertT == pertT_unq(pi));
            trials_num_max = max(trials_num_max, length(idx_trialsPertT{pi}));
        end
        
        % put the trials in a new cell mat, n_trials * n_pertT
        celltmp1 = cell(trials_num_max, length(pertT_unq) + 1);
        % save the 1st column as un-perturbed 
        celltmp1(1:trials_num_max,1) = celltmp(1:trials_num_max,1);
        % save the other columns as perturbation according to the pert Time
        for pi = 1:length(pertT_unq)
            if trials_num_max == length(idx_trialsPertT{pi})
            celltmp1(1:trials_num_max,pi+1) = ...
                celltmp(idx_trialsPertT{pi},2);
            else 
                celltmp1(1:length(idx_trialsPertT{pi}),pi+1) = ...
                    celltmp(idx_trialsPertT{pi},2);
            end
        end
        data(1,frc_i,dist_i,1:7,1:pertT_num) = celltmp1(1:7,1:pertT_num);
%         data(1,frc_i,dist_i,1:15,pertT_num+1) = celltmp(1:5,3); % The stochastic ones pert
    end
end
save('data/processedData/ss4000_4002.mat', 'data'); % 12N perturbation, various time

%%
fce_list = [15 20 25];
dist_list = [2.5 5.0 7.5];
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 5;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3987_3999/f_compare.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

fh(1) = figure();
hold on;
for pi = 1:(pertT_num)
    clf;
    for dist_i = 1:size(data,3) % for each spring
        % 1. plot the perturbed force in the first panel
        trial_num = size(data,4); % 4 for spring data, 5 for human data
        axh(1,dist_i) = subplot(4,3,dist_i); hold on;                     % plot PF
        celltmp1 = reshape(data(1,1,dist_i,:,:),trial_num,pertT_num);
        idx_release = find(celltmp1{1,pi}.ts == 5);
        t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
        ylabel('command force')
        % 2. plot the perturbed velocity in the following panels
        for fce_i = 1:size(data,2)
            axh(dist_i+1,fce_i) = subplot(4,3,dist_i+(fce_i*3)); hold on;         % plot each response
            celltmp1 = reshape(data(1,fce_i,dist_i,:,:),trial_num,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5])
                plot(t, celltmp1{ti,1}.ov(2,:), 'color', [0.5 0.5 0.5]);
%                plot(t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
%                 plot(t, celltmp1{ti,1}.ox(2,:), 'color', [0.5 0.5 0.5]);
%                  plot(t, celltmp1{ti,1}.x(2,:), 'color', [0.5 0.5 0.5]);
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+4,:));
                plot(t, celltmp1{ti,pi}.ov(2,:), 'color', [1 1 1] - color_arr(dist_i+4,:));
%                 plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(dist_i+4,:));
%                 plot(t, celltmp1{ti,pi}.ox(2,:), 'color', [1 1 1] - color_arr(dist_i+4,:));
%                 plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(dist_i+4,:));
                plot(t, -celltmp1{ti,pi}.Fp(2,:));
            end
        end
    end
        % plot notes here: 
        linkaxes(axh(:), 'x'); xlim([-0.1 1.5]);
        
        sgtitle_str = 'unperturbed and perturbed force'; label_str = 'censored force (N)';
%         linkaxes(axh(2:end,:), 'y'); ylim([-12 28]); %force 
%         linkaxes(axh(2:end,:), 'y'); ylim([0.45 0.8]); %position 
%         sgtitle_str = 'unperturbed and perturbed velocity'; label_str = 'velocity (m/s)';
        linkaxes(axh(2:end,:), 'y'); ylim([-1 1.2]); % velocity
%         
        for fce_i = 1:3
            for dist_i = 1:3
            subplot(4,3,fce_i*3 + dist_i); 
            title_str = (['force' num2str(fce_list(fce_i)) 'N dist' num2str(dist_list(dist_i)) 'cm']); 
            ylabel('N');
            subplot(4,3,fce_i*3 + dist_i); title(title_str); ylabel(label_str);
            xlabel('time');
            end
        end
        set(gcf, 'position', [0,0, 1080, 680]);
        sgtitle(sgtitle_str);
        fname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3987_3999/optv_comparep' num2str(pi) '.png'];
         saveas(gcf, fname);
        frame = getframe(gcf);
%         writeVideo(v,frame);
end
close(v);  

%%%%%%%
%% see what is the x and f change during the movement  
% plot: 4-row * 2 -col 
% 1,1: The original and perturbed command force;  1,2: The command force net effect 
% 2,1: The original and perturbed position;       2,2: The resulted position net effect;
% 3,1: The original and perturbed censored force; 3,2: The resulted censored net force;
% 4,1: None;                                      4.2: The psudo-stiffness: dF/dx
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 3%1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:13),5,13); % 1 no -ert and 5 pert
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;  % this f_avg might not be right as it is 'centralized' after the release.
            
            % A TRICK TO UPDATE F_AVG HERE
            if_favgupdate = 1;
            if (if_favgupdate)
                f_avg_br_pert = 0;  % the force before released.
                t_tmp = 0;          % the number of trials
                for tti = 1:size(celltmp1,1)
                    if isempty(celltmp1{tti,pi}) || pi == 1
                        continue;
                    end
                    idx_release = find(celltmp1{tti,pi}.ts == 5);
                    t = celltmp1{tti,pi}.t - celltmp1{tti,pi}.t(idx_release(1));
                    idx_tmp = find(t>=-0.1 & t<0);
                    f_tmp = mean(celltmp1{tti,pi}.f(2,idx_tmp));
                    t_tmp = t_tmp + 1;
                    f_avg_br_pert = f_avg_br_pert + f_tmp;
                    celltmp1{tti,pi}.tshift = t;
                end
                f_avg_br_pert = f_avg_br_pert/t_tmp;
                f_avg_bfr = mean(f_avg(t_grids > -0.1 & t_grids<0));    % before release
                f_diff = f_avg_br_pert - f_avg_bfr;
                f_avg_upd = f_avg + f_diff;         % the before-release value are same now
                
                if (ifplot)
                    clf;
                    hold on;
                    plot(t_grids, f_avg_upd, 'r.');
                    plot(t_grids, f_avg, 'b.');
                    for tti = 1:size(celltmp1,1)
                        if isempty(celltmp1{tti,pi}) || pi == 1
                            continue;
                        end
                        plot(celltmp1{tti,pi}.tshift, celltmp1{tti,pi}.f(2,:), 'Color', [0.5 0.5 0.5]);
                    end
                    
                end
                if (pi~=1)
                    f_avg = f_avg_upd;
                end
            end
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
            plot(t_grids, x_avg);
%             plot(t_grids, v_avg);
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; % subtracted dF/dx
                plot(t_grids, (f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted x
                plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                ylim([-0.05 0.04]);
%                 plot(t_grids(pert_idx), v_dat(pert_idx) - v_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
                scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
                axh(6) = subplot(4,2,6); % subtracted F
                scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
                
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
                scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
                axh(6) = subplot(4,2,6); % subtracted F
                scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        if (pi~=1)
            ylabel(axh(8), 'K (N/m)'); subplot(axh(8));ylim(-1000+[0 1000]);
        end
%         saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness300ms' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
        close all;
    end
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validating the method to calculate stiffness by dF/dx
%% see what is the x and f change during the movement  
% plot: 4-row * 2 -col 
% 1,1: The original and perturbed command force;  1,2: The command force net effect 
% 2,1: The original and perturbed position;       2,2: The resulted position net effect;
% 3,1: The original and perturbed censored force; 3,2: The resulted censored net force;
% 4,1: psudo-stiffness of dF/dx, release;         4.2: The psudo-stiffness: dF/dx
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
k_stfcoef = 13/20;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1:3%1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:13),5,13); % 1 no -ert and 5 pert
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
%                 x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;  % this f_avg might not be right as it is 'centralized' after the release.
            
            % A TRICK TO UPDATE F_AVG HERE
            if_favgupdate = 1;
            if (if_favgupdate)
                f_avg_br_pert = 0;  % the force before released.
                t_tmp = 0;          % the number of trials
                for tti = 1:size(celltmp1,1)
                    if isempty(celltmp1{tti,pi}) || pi == 1
                        continue;
                    end
                    idx_release = find(celltmp1{tti,pi}.ts == 5);
                    t = celltmp1{tti,pi}.t - celltmp1{tti,pi}.t(idx_release(1));
                    idx_tmp = find(t>=-0.1 & t<0);
                    f_tmp = mean(celltmp1{tti,pi}.f(2,idx_tmp));
                    t_tmp = t_tmp + 1;
                    f_avg_br_pert = f_avg_br_pert + f_tmp;
                    celltmp1{tti,pi}.tshift = t;
                end
                f_avg_br_pert = f_avg_br_pert/t_tmp;
                f_avg_bfr = mean(f_avg(t_grids > -0.1 & t_grids<0));    % before release
                f_diff = f_avg_br_pert - f_avg_bfr;
                f_avg_upd = f_avg + f_diff;         % the before-release value are same now
                
                if (ifplot)
                    clf;
                    hold on;
                    plot(t_grids, f_avg_upd, 'r.');
                    plot(t_grids, f_avg, 'b.');
                    for tti = 1:size(celltmp1,1)
                        if isempty(celltmp1{tti,pi}) || pi == 1
                            continue;
                        end
                        plot(celltmp1{tti,pi}.tshift, celltmp1{tti,pi}.f(2,:), 'Color', [0.5 0.5 0.5]);
                    end
                    
                end
                if (pi~=1)
                    f_avg = f_avg_upd;
                end
            end
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
            plot(t_grids, x_avg);
%             plot(t_grids, v_avg);
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
 
            axh(7) = subplot(4,2,7);
            k_avgest = f_avg ./ (x_avg - x_avg(end));
            plot(t_grids, -k_avgest);
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                     plot(t, celltmp1{ti,pi}.ox(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
%                 x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; % subtracted dF/dx
                plot(t_grids, -(f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted x
                plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                ylim([-0.05 0.04]);
%                 plot(t_grids(pert_idx), v_dat(pert_idx) - v_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
                scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
                axh(6) = subplot(4,2,6); % subtracted F
                scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
                
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
                scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
                axh(6) = subplot(4,2,6); % subtracted F
                scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'K (N/m)'); subplot(axh(7));ylim(+[0 K_list(dist_i)*1.5]);
        yline(K_list(dist_i)*k_stfcoef);
        yline(K_list(dist_i));
        if (pi~=1)
            ylabel(axh(8), 'K (N/m)'); subplot(axh(8));ylim(+[0 K_list(dist_i)*1.5]);
            yline(K_list(dist_i)*k_stfcoef); yline(K_list(dist_i));
        end
%           saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness300ms' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
%         close all;
    end
    
end

%% Compare the displacement difference at different position
% Also, compare the condition where there is no release (left) and with
% release (right)
% plot: 1-row * 2 -col 
% col1: without release: 
%           x-pert time, 
%           y-displacement before perturb, 
%           z-displacement different in release

% The left column, only plot the displacement (raw) 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4000_4002.mat', 'data');
figure(); 
color_arr = colormap('lines');
close all;
fh1 = figure();
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1 %1:size(data,2)
    for dist_i = 1 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));

            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1);
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                axh(1) = subplot(1,2,1); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
                plot3(t_grids - t_offset, pos_offset ,x_dat - pos_offset, 'color', color_arr(4+2,:));
                

                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 

                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  

%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end

        

        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
%         close all;
    end
    
end
subplot(axh(1));
view([0, 80]); 
xlabel('time after perturb');
ylabel('robot position');
zlabel('induced position difference');
title('perturb no release')

% The right column, only plot displacement difference 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
F_list = [15, 20, 25];
K_list = [640, 320, 160];
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 2 %1:size(data,2)
    for dist_i = 2 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1);
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   

                axh(2) = subplot(1,2,2); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
                plot3(t_grids - t_offset, pos_offset ,x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
%                 
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            
            % plot the subtracted position and force, in another
            % figure/panel
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
    end
    
end
subplot(axh(2));
view([0, 80]); 
xlabel('time after perturb');
ylabel('robot position');
zlabel('induced position difference');
title('perturb during release');
saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossPositions_quene.png')

subplot(axh(1));
view([0, 0]);  zlim([-0.04 0.04]);
subplot(axh(2));
view([0, 0]); zlim([-0.04 0.04]);
saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossPositions_overlay.png')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sanity check: why the position difference is so huge in some of the pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the displacement difference at different position
% Also, compare the condition where there is no release (left) and with
% release (right)
% plot: 1-row * 2 -col 
% col1: without release: 
%           x-pert time, 
%           y-displacement before perturb, 
%           z-displacement different in release

% The left column, only plot the displacement (raw) 
% The right column, only plot displacement difference 

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
F_list = [15, 20, 25];
K_list = [640, 320, 160];
color_arr = colormap('lines');
close all;
fh1 = figure();
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 2 %1:size(data,2)
    for dist_i = 2 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1);
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   

%                 axh(2) = subplot(1,2,2); hold on;% subtracted x
                hold on;
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
%                 plot3(t_grids - t_offset, pos_offset ,x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                plot3(t_grids - t_offset, pos_offset ,x_dat, 'color', color_arr(4+dist_i,:));
                plot3(t_grids - t_offset, pos_offset ,x_avg, 'color', [0.5 0.5 0.5]);
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
%                 
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            
            % plot the subtracted position and force, in another
            % figure/panel
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
    end
    
end
% subplot(axh(2));
view([0, 80]); 
xlabel('time after perturb');
ylabel('robot position');
zlabel('induced position difference');
title('perturb during release');
saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossPositions_SanityCheck.png')
% 
% subplot(axh(1));
% view([0, 0]);  zlim([-0.04 0.04]);
% subplot(axh(2));
% view([0, 0]); zlim([-0.04 0.04]);
% saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossPositions_overlay.png')



%% Compare the displacement difference at different Velocity
% Also, compare the condition where there is no release (left) and with
% release (right) 
% Arrange the displacement difference when there is release on the velocity
% dot product with the perturb force 
% plot: 1-row * 2 -col 
% col1: without release: 
%           x-pert time, 
%           y-velocity before perturb, 
%           z-displacement different in release

% The left column, only plot the displacement (raw) 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4000_4002.mat', 'data');
figure(); 
color_arr = colormap('lines');
close all;
fh1 = figure();
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1 %1:size(data,2)
    for dist_i = 1 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));

            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,3}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,3}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,3}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1);
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                v_offset = v_dat(offset_idx);
                
                fp_effectv = fp_dat .* v_dat;
                v_offset = sum(fp_effectv);
                
                ifplot = 1;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                axh(1) = subplot(1,2,1); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
                vel_offset = v_offset * ones(size(v_avg));
                plot3(t_grids - t_offset, vel_offset ,x_dat - pos_offset, 'color', color_arr(4+2,:));
                

                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 

                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  

%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end

        

        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
%         close all;
    end
    
end
subplot(axh(1));
view([0, 80]); ylim([-100 400]);
xlabel('time after perturb');
ylabel('robot velocity');
zlabel('induced position difference');
title('perturb no release')

% The right column, only plot displacement difference 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
F_list = [15, 20, 25];
K_list = [640, 320, 160];
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 2 %1:size(data,2)
    for dist_i = 2 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1) + 150; % 100ms after perturb~=0
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                v_offset = v_dat(offset_idx); % change this into the convolution between velocity and perturb force. 
                
                fp_effectv = fp_dat .* v_dat;
                v_offset = sum(fp_effectv);
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   

                axh(2) = subplot(1,2,2); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
                vel_offset = v_offset * ones(size(v_avg));
                
                plot3(t_grids - t_offset, vel_offset ,x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
%                 
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            
            % plot the subtracted position and force, in another
            % figure/panel
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
    end
    
end
subplot(axh(2));
view([0, 80]); 
ylim([-100 400]);
xlabel('time after perturb');
ylabel('Fp effect on v');   % sum(f_offset.*fp)
zlabel('induced position difference');
title('perturb during release');
% saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossPositions_quene.png')

% subplot(axh(1));
% view([0, 0]);  zlim([-0.04 0.04]);
% subplot(axh(2));
% view([0, 0]); zlim([-0.04 0.04]);
% saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossForceEffectVelocity_overlay.png')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do dF/dx change in the subject data

%% see what is the x and f change during the movement  
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3818_3828.mat', 'data');
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3913_3921.mat', 'data');
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
        psudoK_mat = zeros(5,7); % 
        for pi = 1:8%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:13),5,13); % 1 no -ert and 5 pert
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:8),10,8); % 1 no -ert and 5 pert
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            f_avg = f_avg/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(3,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(3,2,3);
            plot(t_grids, x_avg);
%             plot(t_grids, x_avg - x_avg(1));
            axh(5) = subplot(3,2,5);
            plot(t_grids, f_avg);
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(3,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(3,2,4); hold on;% subtracted x
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(3,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(3,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(3,2,4); % subtracted x
                plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(3,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(3,2,4); % subtracted x
                scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
                axh(6) = subplot(3,2,6); % subtracted F
                scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
                
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(3,2,4); % subtracted x
                scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
                axh(6) = subplot(3,2,6); % subtracted F
                scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
%         saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
        close all;
    end
    
end


%% The following code answers the question:
% Does the mass change across conditions? 
% The mass measureing method is using the point immediate after release to
% get the mass estimation. The conclusion is based on the following 
% assumption:
% 1. The force transducer measures the net tention (force) before relese,
%    and only measures the acceleration of the free part after release.  
% 2. The acceleration of two parts of force transducer should be the same. 
%
% Hence, the before release force 
%   F0- = (m1 + m2)*a           ...(1)
%   F0+ = m2*a                  ...(2)
% F0-: force measured before release. m1, the mass at subject end, m2, the
% mass at robot end, which will be a free end after release. a, the two
% ends acceleration
% F0+: force measured immediately after release. 
% Hence, the mass at the subject end can be calculated as: 
% m1/m2 = (F(0-) - F(0+))/(F0+)
% Which means, the mass at the subject end can be described by the
% perprotion of the force measured change and the force immediately after
% release. 
%
% The code using the aforementioned theory, and get the conclusion that
% mass have increased by the condition
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
f_ratio_cell = cell(3,3);
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); %
        pi = 1; % only see the unperturbed trials
%         fh(pi,1) = figure(); hold on;
        
        celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
        idx_release = find(celltmp1{1,pi}.ts == 5);
        t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
        
        
        % calculate the Unperturbed situation, x and f
        x_avg = zeros(1, length(t_grids));
        f_avg = zeros(1, length(t_grids));
        fp_avg= zeros(1, length(t_grids));
        cti = 0; % count how many trials are added up
        f_ratio = zeros(1,size(celltmp1,1));
        for ti = 1:1:size(celltmp1,1)
            if isempty(celltmp1{ti,1})
                continue;
            end
            cti = cti + 1;
            idx_release = find(celltmp1{ti,1}.ts == 5);
            t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
            idx_t = find(t>=t_interest(1) & t<=t_interest(2));
            length(idx_t)
            % intropolate (x, f, Fp) to t_grids
            
            x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            
            ifplot = 1; % controls whether plot or not
            t_idx = find(t_grids>0 & t_grids < 0.1); 
            [f_pk_neg,locs1] = findpeaks(-f_dat(t_idx), 'NPeaks', 1);
            [f_pk_pos,locs2] = findpeaks(f_dat(t_idx), 'NPeaks', 1);
            if (ifplot)
%                 clf;
                figure('unit', 'inch', 'position', [12 0 5 5] );
                axh(1) = subplot(3,1,1);  hold on;
                plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                plot(t_grids, fp_dat, 'r', 'Marker', '.');
                
                axh(2) = subplot(3,1,2);  hold on;
                plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                plot(t_grids, x_dat, 'r', 'Marker', '.');
                
                axh(3) = subplot(3,1,3); hold on;
                plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                plot(t_grids, f_dat, 'r', 'Marker', '.'); hold on;
                
                plot(t_grids(t_idx(locs1)),-f_pk_neg, 'marker', '.', 'markersize', 20);
                plot(t_grids(t_idx(locs2)), f_pk_pos, 'marker', '.', 'markersize', 20);
                
                linkaxes(axh, 'x');
                xlim(axh(3), [-0.1 0.1]);
            end    
            f0_bef = mean(f_dat(t_grids<0)); 
%             f0_aft = mean(f_dat(t_grids>0.01 & t_grids<0.04));
            f0_aft = 1/2*(f_pk_pos + (-f_pk_neg));
            f_ratio(ti) = f0_bef/f0_aft;
            
           f_ratio_cell{fce_i,dist_i} = f_ratio;
        end 
    end    
end
close all;
%
fh = figure();
clear axh;
for fce_i = 1:3 
    for dist_i = 1:3 
        % for each panel
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        f_ratio = f_ratio_cell{fce_i,dist_i};
        fratio_num = length(f_ratio); 
        scatter(1:fratio_num, f_ratio, 5, 'MarkerFaceColor' , 'b'); 
        title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
    end
end
linkaxes(axh(:), 'xy');
sgtitle('F0-/F0+ on different conditions');  


% if I convert the relationship into the mass...
% Knowing from previous calculation, Mr = 1.27, Ms0 = 0.33
% F0-/F0+ = (Madd + Ms0 + Mr) / Mr
% Madd = F0-/F0+ * Mr - Mr - Ms0
%      = F0-/F0+ * 0.94-1.27
fh = figure();
clear axh;
for fce_i = 1:3 
    for dist_i = 1:3 
        % for each panel
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        f_ratio = f_ratio_cell{fce_i,dist_i};
        fratio_num = length(f_ratio); 
        mass_add = f_ratio * 0.94 - 1.27; % kg
        scatter(1:fratio_num, mass_add, 5, 'MarkerFaceColor' , 'b'); 
        title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
        xlabel('trial');
        ylabel('mass (kg)');
    end
end
linkaxes(axh(:), 'xy');
sgtitle('mass add on different conditions');  

%%%%%%%%%%%%%%%%%%%%%% do a barplot or a errorbar with line
fh = figure('unit', 'inch', 'position', [0 0 3 3]);
clear axh;
mass_mean = zeros(3,3); 
mass_std  = zeros(3,3); 
for fce_i = 1:3 
    for dist_i = 1:3 
        % for each panel
%         axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        f_ratio = f_ratio_cell{fce_i,dist_i};
        fratio_num = length(f_ratio); 
        mass_add = f_ratio * 0.94 - 1.27; % kg
%         scatter(1:fratio_num, mass_add, 5, 'MarkerFaceColor' , 'b'); 
        %.......!!!!!!!!!!!!! TODO!!!!!
        mass_mean(fce_i,dist_i) = mean(mass_add); 
        mass_std(fce_i,dist_i)  = std(mass_add);
    end
end

hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, mass_mean(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, mass_mean(:,dist_i), mass_std(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([0 10]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
legend(lnh, {'2.5cm', '5cm', '7.5cm'});
ylabel('mass (kg)');
% linkaxes(axh(:), 'xy');
title('mass cross conditions (subject)');  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The same measurement using the spring data
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
f_ratio_cell = cell(3,3);
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); %
        pi = 1; % only see the unperturbed trials
%         fh(pi,1) = figure(); hold on;
        
        celltmp1 = reshape(data(1,fce_i,dist_i,1:7,1:6),7,6); % 1 no -ert and 5 pert
        idx_release = find(celltmp1{1,pi}.ts == 5);
        t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
        
        
        % calculate the Unperturbed situation, x and f
        x_avg = zeros(1, length(t_grids));
        f_avg = zeros(1, length(t_grids));
        fp_avg= zeros(1, length(t_grids));
        cti = 0; % count how many trials are added up
        f_ratio = zeros(1,size(celltmp1,1));
        for ti = 1:1:size(celltmp1,1)
            if isempty(celltmp1{ti,1})
                continue;
            end
            cti = cti + 1;
            idx_release = find(celltmp1{ti,1}.ts == 5);
            t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
            idx_t = find(t>=t_interest(1) & t<=t_interest(2));
            length(idx_t);
            % intropolate (x, f, Fp) to t_grids
            
            x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            
            ifplot = 1; % controls whether plot or not
            t_idx = find(t_grids>0 & t_grids < 0.1); 
            [f_pk_neg,locs1] = findpeaks(-f_dat(t_idx), 'NPeaks', 1);
            [f_pk_pos,locs2] = findpeaks(f_dat(t_idx), 'NPeaks', 1);
            if (ifplot)
%                 clf;
                figure('unit', 'inch', 'position', [12 0 5 5] );
                axh(1) = subplot(3,1,1);  hold on;
                plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                plot(t_grids, fp_dat, 'r', 'Marker', '.');
                
                axh(2) = subplot(3,1,2);  hold on;
                plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                plot(t_grids, x_dat, 'r', 'Marker', '.');
                
                axh(3) = subplot(3,1,3); hold on;
                plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                plot(t_grids, f_dat, 'r', 'Marker', '.'); hold on;
                
                plot(t_grids(t_idx(locs1)),-f_pk_neg, 'marker', '.', 'markersize', 20);
                plot(t_grids(t_idx(locs2)), f_pk_pos, 'marker', '.', 'markersize', 20);
                
                linkaxes(axh, 'x');
                xlim(axh(3), [-0.1 0.1]);
            end    
            f0_bef = mean(f_dat(t_grids<0)); 
%             f0_aft = mean(f_dat(t_grids>0.01 & t_grids<0.04));
            f0_aft = 1/2*(f_pk_pos + (-f_pk_neg));
            f_ratio(ti) = f0_bef/f0_aft;
            
           f_ratio_cell{fce_i,dist_i} = f_ratio;
        end 
    end    
end
close all;
%
fh = figure();
clear axh;
for fce_i = 1:3 
    for dist_i = 1:3 
        % for each panel
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        f_ratio = f_ratio_cell{fce_i,dist_i};
        fratio_num = length(f_ratio); 
        scatter(1:fratio_num, f_ratio, 5, 'MarkerFaceColor' , 'b'); 
        title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
    end
end
linkaxes(axh(:), 'xy');
sgtitle('F0-/F0+ on different conditions');  


% % if I convert the relationship into the mass...
% % Knowing from previous calculation, Mr = 1.27, Ms0 = 0.33
% % F0-/F0+ = (Madd + Ms0 + Mr) / Mr
% % Madd = F0-/F0+ * Mr - Mr - Ms0
% %      = F0-/F0+ * 0.94-1.27
% fh = figure();
% clear axh;
% for fce_i = 1:3 
%     for dist_i = 1:3 
%         % for each panel
%         axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
%         f_ratio = f_ratio_cell{fce_i,dist_i};
%         fratio_num = length(f_ratio); 
%         mass_add = f_ratio * 0.94 - 1.27; % kg
%         scatter(1:fratio_num, mass_add, 5, 'MarkerFaceColor' , 'b'); 
%         title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
%         xlabel('trial');
%         ylabel('mass (kg)');
%     end
% end
% ylim([0 0.5]);
% linkaxes(axh(:), 'xy');
% sgtitle('mass add on different conditions');  

%%%%%%%%%%%%%%%%%%%%%% do a barplot or a errorbar with line
fh = figure('unit', 'inch', 'position', [0 0 3 3]);
clear axh;
mass_mean = zeros(3,3); 
mass_std  = zeros(3,3); 
for fce_i = 1:3 
    for dist_i = 1:3 
        % for each panel
%         axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        f_ratio = f_ratio_cell{fce_i,dist_i};
        fratio_num = length(f_ratio); 
        mass_add = f_ratio * 0.94 - 1.27; % kg
%         scatter(1:fratio_num, mass_add, 5, 'MarkerFaceColor' , 'b'); 
        %.......!!!!!!!!!!!!! TODO!!!!!
        mass_mean(fce_i,dist_i) = mean(mass_add); 
        mass_std(fce_i,dist_i)  = std(mass_add);
    end
end

hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, mass_mean(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, mass_mean(:,dist_i), mass_std(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([0 10]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
legend(lnh, {'2.5cm', '5cm', '7.5cm'});
ylabel('mass (kg)');
% linkaxes(axh(:), 'xy');
title('mass cross conditions (spring)');  



%% use the x and f, do the dF/dx to show how the stiffness works here
% Stiffness measurement using definition dF/dx, have 72 figures * 8 pannels
% assume hand is only 3kg. 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
for fce_i = 1:size(data,2)
% for fce_i = 3
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
        psudoK_mat = zeros(7,7); % 
        for pi = 1:8%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            clear axh;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            f_avg = f_avg/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
            plot(t_grids, x_avg);
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
%             plot(t_grids, x_avg - x_avg(1));
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
            axh(7) = subplot(4,2,7);
            plot(t_grids, -f_avg ./ (x_avg - x_settled));
            
            % find the reference value 
            psudo_stiffness = -f_avg ./ (x_avg - x_settled); 
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            yline(psudo_stiffness1, 'linewidth', 2);
%             ylim([0 2000]);
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; 
                plot(t_grids, -(f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_avg; 
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = 0;
                [~, x_net_idx] = min(x_net_tmp);
                k_est = -(f_dat - f_avg)./(x_dat - x_avg); 
                k_est_pt = k_est(x_net_idx); 
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                plot(t_grids(x_net_idx), k_est_pt, 'marker', 'o', 'markersize', 10); 
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted x
                plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
% %                 % Find and plot on the original data point 
% %                 pert0idx  = find(abs(fp_dat) > 0.01);
% %                 pert0idx = pert0idx(1);  
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
% %                 
% %                 x0 = x_dat(pert0idx) - x_avg(pert0idx);
% %                 f0 = f_dat(pert0idx) - f_avg(pert0idx);
% %                 
% %                 % Find and plot on the peak data point 
% %                 [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
% %                 
% %                 x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
% %                 f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
% %                 
% %                 psudoK = (f1-f0)/(x1-x0);
% %                 psudoK_mat(ti,pi-1) = psudoK;
            end
            try % if has axh8, plot, ifnot, noplot 
                linkaxes(axh(7:8), 'y'); 
                yline(axh(8),psudo_stiffness1, 'linewidth', 2);
            catch
            end
            ylim(axh(7), [0 1000]);
            linkaxes(axh, 'x');
            % plot notes here: 
            xlim(axh(1), [-0.1 1.36]);
        
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'df/dx');
        try
        ylabel(axh(2), 'dFp'); 
        ylabel(axh(4), 'dx'); 
        ylabel(axh(6), 'dF'); 
        ylabel(axh(8), 'dF/dx'); 
        catch
        end
        
%          saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = psudoK_mat;
%         close all;
    end
    
end

%%%%%%%%%%%%%%%%%
figure(); 
clear axh; 
dat_mean_cc = zeros(3,3); 
dat_std_cc  = zeros(3,3); 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        psudoK_Mat = psudoK_cell{fce_i,dist_i};
        % need to detect outlaiers before get mean and std
        dat_arr = psudoK_Mat(:); 
        dat_arr_outlairidx = isoutlier(dat_arr);
        dat_mean_cc(fce_i,dist_i) = mean(dat_arr(~dat_arr_outlairidx));
        dat_std_cc(fce_i,dist_i) = std(dat_arr(~dat_arr_outlairidx));
        ifplot = 1; 
        if (ifplot)
            fh_tmp = figure(); hold on;
            plot(dat_arr, '.'); 
            plot(find(dat_arr_outlairidx), dat_arr(dat_arr_outlairidx), 'marker', 'o', 'markersize', 10, 'color', 'r');
            title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
        end
        try 
            close(fh_tmp)
        catch
        end
        plot(psudoK_Mat); 
        xlabel('pert time'); 
        ylabel('damping estimation');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('mass estimation using dF/da');
linkaxes(axh);
ylim([-200 1000]);

% better figure; 
peak_time = [0.1:0.025:0.25]; % s... Need to change here! 
fh = figure('unit', 'inch', 'position', [0 0 5 5]); 
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoK_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        xlabel('pert time'); 
        ylabel('K (N/m)');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('Subject stiffness after release, dF/dx');
linkaxes(axh);
ylim([-200 1000]);

% corss condition plot 
fh = figure('unit', 'inch', 'position', [0 0 3 3]); 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc(:,dist_i), dat_std_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([-200 1000]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
% legend(lnh, {'640N/m', '320N/m', '160N/m'});
legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
ylabel('stiffness (N/m)');
% linkaxes(axh(:), 'xy');
title('stiffness cross conditions (subject)');  


%% use the x and f, do the dF/dx to show how the stiffness works here
% Stiffness measurement using definition dF/dx, have 72 figures * 8 pannels
% IN THIS CODE BLOCK, TRY ANDY'S IDEA ON AVERAGE THE DATA THEN DO THE df/dx
% assume hand is only 3kg. 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
for fce_i = 1:size(data,2)
% for fce_i = 3
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
        psudoK_mat = zeros(7,7); % 
        for pi = 2:8%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            clear axh;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            
            x_avg_pert = zeros(1, length(t_grids));
            f_avg_pert = zeros(1, length(t_grids));
            fp_avg_pert= zeros(1, length(t_grids));
            
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,1})
                if isempty(celltmp1{ti,pi})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                idx_tp = find(t>=t_interest(1) & t<=t_interest(2));
                x_dat_pert = interp1(t(idx_tp), celltmp1{ti,pi}.x(2,idx_tp), t_grids, 'linear', 'extrap'); % check...
                f_dat_pert = interp1(t(idx_tp), celltmp1{ti,pi}.f(2,idx_tp), t_grids, 'linear', 'extrap'); % check...
                fp_dat_pert= interp1(t(idx_tp), celltmp1{ti,pi}.Fp(2,idx_tp), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    plot(t_grids, fp_dat_pert, 'g', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    plot(t_grids, x_dat_pert, 'g', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                    plot(t_grids, f_dat_pert, 'g', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
                
                x_avg_pert = x_avg_pert + x_dat_pert;
                f_avg_pert = f_avg_pert + f_dat_pert;
                fp_avg_pert = fp_avg_pert + fp_dat_pert;
                
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            f_avg = f_avg/cti;
            
            x_avg_pert = x_avg_pert/cti;
            f_avg_pert = f_avg_pert/cti;
            fp_avg_pert = fp_avg_pert/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); hold on;
            plot(t_grids, fp_avg);
            plot(t_grids, fp_avg_pert); 
            axh(3) = subplot(4,2,3); hold on;
            plot(t_grids, x_avg);
            plot(t_grids, x_avg_pert);
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
%             plot(t_grids, x_avg - x_avg(1));
            axh(5) = subplot(4,2,5); hold on;
            plot(t_grids, f_avg);
            plot(t_grids, f_avg_pert); 
            axh(7) = subplot(4,2,7);
            plot(t_grids, -f_avg ./ (x_avg - x_settled));
            
            axh(2) = subplot(4,2,2); 
            plot(t_grids, fp_avg_pert - fp_avg);
            axh(4) = subplot(4,2,4); 
            plot(t_grids, x_avg_pert - x_avg);
            axh(6) = subplot(4,2,6); 
            plot(t_grids, f_avg_pert - f_avg);
            axh(8) = subplot(4,2,8); 
            plot(t_grids, (f_avg_pert - f_avg)./(x_avg_pert - x_avg));
            
            % find the reference value 
            psudo_stiffness = -f_avg ./ (x_avg - x_settled); 
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            yline(psudo_stiffness1, 'linewidth', 2);
%             ylim([0 2000]);

            % find the stiffness prediction here 
            x_net = x_avg_pert - x_avg;
            f_net = f_avg_pert - f_avg;
            k_net = -(f_avg_pert - f_avg)./(x_avg_pert - x_avg);
            [~,pert_peak_idx] = max(abs(fp_avg_pert));
            x_net_tmp = x_net; x_net_tmp(1:pert_peak_idx) = 0; 
            [~, pos_peak_idx] = min(x_net_tmp); 
            subplot(axh(4)); hold on; 
            plot(t_grids(pos_peak_idx), x_net(pos_peak_idx), '*');
            subplot(axh(6)); hold on;
            % option one, just use the exact time value as the f
%             plot(t_grids(pos_peak_idx), f_net(pos_peak_idx), '*');
            % option two, use the +-50ms mininum as the f value, +-50ms is
            % +- 25 datapoints
            f_net_tmp = f_net; 
            t_range = 25; % points
            f_net_tmp([1:pos_peak_idx-t_range, pos_peak_idx+t_range:end]) = nan;
            [~, fce_peak_idx] = nanmin(f_net_tmp); 
            plot(t_grids(fce_peak_idx), f_net(fce_peak_idx), '*');
            psudo_K_tmp =-f_net(fce_peak_idx)/x_net(pos_peak_idx);
            
            subplot(axh(8)); hold on;
            plot(t_grids(pos_peak_idx), k_net(pos_peak_idx), '*');
            
% %             psudoK_mat(pi-1,:) = k_net(pos_peak_idx); % ... some averaged value, exact point 
            psudoK_mat(pi-1,:) = psudo_K_tmp; % ... some averaged value, shift minimum
            % plot the perturbed one, -perturbed
%             for ti = 1%:size(celltmp1,1) ... looks like only averaged viable trial 
%                 if isempty(celltmp1{ti,pi}) || pi == 1
%                     continue;
%                 end
%                 
%                 idx_release = find(celltmp1{ti,pi}.ts == 5);
%                 t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 subplot(axh(1)); hold on;
%                 plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
%                 
%                 subplot(axh(3)); hold on;
%                 x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
%                 plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
% %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
%                 
%                 subplot(axh(5)); hold on;
%                 plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
%                 
%                 idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t)
% 
%                 fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
%                 x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
%                 f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
% %                 linkaxes(axh(1:3:5), 'x');
%                 
%                 ifplot = 0;
%                 if (ifplot)
%                     clf;
%                     subplot(2,1,1);  hold on;
%                     plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
%                     plot(t_grids, x_dat, 'r', 'Marker', '.');
%                     
%                     subplot(2,1,2); hold on;
%                     plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
%                     plot(t_grids, f_dat, 'r', 'Marker', '.');
%                 end
%                    
%                 % plot the subtraction in other panels 
%                 axh(2) = subplot(4,2,2); hold on; % subtracted Fp
%                 plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
%                 axh(4) = subplot(4,2,4); hold on;% subtracted x
%                 plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
%                 axh(6) = subplot(4,2,6); hold on;% subtracted F
%                 plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
%                 axh(8) = subplot(4,2,8); hold on; 
%                 plot(t_grids, -(f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
%                 
%                 [~,fp_max_idx] = max(abs(fp_dat));
%                 x_net = x_dat - x_avg; 
%                 x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = 0;
%                 [~, x_net_idx] = min(x_net_tmp);
%                 k_est = -(f_dat - f_avg)./(x_dat - x_avg); 
%                 k_est_pt = k_est(x_net_idx); 
%                 psudoK_mat(pi-1,ti) = k_est_pt;
%                 
%                 plot(t_grids(x_net_idx), k_est_pt, 'marker', 'o', 'markersize', 10); 
%                 
%                 % There will be a force and position peak at 0~0.4s after
%                 % the start of the perturbation 
%                 pert_idx = find(abs(fp_dat) > 0.5); 
%                 pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
%                 % plot out the perturbed position 
%                 axh(2) = subplot(4,2,2);  % subtracted Fp
%                 plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
%                 axh(4) = subplot(4,2,4); % subtracted x
%                 plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
%                 axh(6) = subplot(4,2,6); % subtracted F
%                 plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
%                 
% % %                 % Find and plot on the original data point 
% % %                 pert0idx  = find(abs(fp_dat) > 0.01);
% % %                 pert0idx = pert0idx(1);  
% % %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% % %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% % %                 axh(4) = subplot(3,2,4); % subtracted x
% % %                 scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
% % %                 axh(6) = subplot(3,2,6); % subtracted F
% % %                 scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
% % %                 
% % %                 x0 = x_dat(pert0idx) - x_avg(pert0idx);
% % %                 f0 = f_dat(pert0idx) - f_avg(pert0idx);
% % %                 
% % %                 % Find and plot on the peak data point 
% % %                 [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
% % %                 [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
% % %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% % %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% % %                 axh(4) = subplot(3,2,4); % subtracted x
% % %                 scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
% % %                 axh(6) = subplot(3,2,6); % subtracted F
% % %                 scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
% % %                 
% % %                 x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
% % %                 f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
% % %                 
% % %                 psudoK = (f1-f0)/(x1-x0);
% % %                 psudoK_mat(ti,pi-1) = psudoK;
%             end
            try % if has axh8, plot, ifnot, noplot 
                linkaxes(axh(7:8), 'y'); 
                yline(axh(8),psudo_stiffness1, 'linewidth', 2);
            catch
            end
            ylim(axh(7), [0 1000]);
            linkaxes(axh, 'x');
            % plot notes here: 
            xlim(axh(1), [-0.1 1.36]);
        
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'df/dx');
        try
        ylabel(axh(2), 'dFp'); 
        ylabel(axh(4), 'dx'); 
        ylabel(axh(6), 'dF'); 
        ylabel(axh(8), 'dF/dx'); 
        catch
        end
        
%          saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = psudoK_mat;
%         close all;
    end
    
end

%%%%%%%%%%%%%%%%%
figure(); 
clear axh; 
dat_mean_cc = zeros(3,3); 
dat_std_cc  = zeros(3,3); 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        psudoK_Mat = psudoK_cell{fce_i,dist_i};
        % need to detect outlaiers before get mean and std
        dat_arr = psudoK_Mat(:); 
        dat_arr_outlairidx = isoutlier(dat_arr);
        dat_mean_cc(fce_i,dist_i) = mean(dat_arr(~dat_arr_outlairidx));
        dat_std_cc(fce_i,dist_i) = std(dat_arr(~dat_arr_outlairidx));
        ifplot = 1; 
        if (ifplot)
            fh_tmp = figure(); hold on;
            plot(dat_arr, '.'); 
            plot(find(dat_arr_outlairidx), dat_arr(dat_arr_outlairidx), 'marker', 'o', 'markersize', 10, 'color', 'r');
            title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
        end
        try 
            close(fh_tmp)
        catch
        end
        plot(psudoK_Mat); 
        xlabel('pert time'); 
        ylabel('damping estimation');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('mass estimation using dF/da');
linkaxes(axh);
ylim([-200 1000]);

% better figure; 
peak_time = [0.1:0.025:0.25]; % s... Need to change here! 
fh = figure('unit', 'inch', 'position', [0 0 5 5]); 
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoK_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        xlabel('pert time'); 
        ylabel('K (N/m)');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('Subject stiffness after release, dF/dx');
linkaxes(axh);
ylim([-200 1000]);

% corss condition plot 
fh = figure('unit', 'inch', 'position', [0 0 3 3]); 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc(:,dist_i), dat_std_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([-200 1000]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
% legend(lnh, {'640N/m', '320N/m', '160N/m'});
legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
ylabel('stiffness (N/m)');
% linkaxes(axh(:), 'xy');
title('stiffness cross conditions (subject)');  
%% a question needs to be answer (adapt code here)... The code answer 2nd of these questions
% damping (viscosity) measurement here, use the definition dF/dv
% 1. use the x and f, do the dF/dx to show how the stiffness works here 
% 2. use a, v and f, do the dF/dv show how the damping works here 

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoD_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
        psudoD_mat = zeros(7,7); % 
        for pi = 1:8%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); 
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
            plot(t_grids, v_avg);
            v_settled = nanmean(v_avg(t_grids>1.0 & t_grids<1.1));
%             plot(t_grids, x_avg - x_avg(1));
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
            axh(7) = subplot(4,2,7);
            plot(t_grids, f_avg ./ (v_avg - v_settled));
%             ylim([0 2000]);
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
%                 x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
%                 plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
%                 plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; 
                plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                
                v_net = v_dat - v_avg; 
                [~, pert_peak_idx] = max(abs(fp_dat));
                v_net_tmp = v_net; v_net_tmp(1:pert_peak_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part 
                damping_est = (f_dat - f_avg)./(v_dat - v_avg);
                damping_est_pt = damping_est(v_peak_idx); 
                psudoD_mat(pi-1,ti) = damping_est_pt;
                
                plot(t_grids(v_peak_idx), damping_est_pt, 'marker', 'o', 'markersize', 10); 
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted x
%                 plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
                plot(t_grids(pert_idx), v_dat(pert_idx) - v_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
% %                 % Find and plot on the original data point 
% %                 pert0idx  = find(abs(fp_dat) > 0.01);
% %                 pert0idx = pert0idx(1);  
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
% %                 
% %                 x0 = x_dat(pert0idx) - x_avg(pert0idx);
% %                 f0 = f_dat(pert0idx) - f_avg(pert0idx);
% %                 
% %                 % Find and plot on the peak data point 
% %                 [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
% %                 
% %                 x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
% %                 f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
% %                 
% %                 psudoK = (f1-f0)/(x1-x0);
% %                 psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            try 
                linkaxes(axh(7:8), 'y'); 
            catch 
            end
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
%         yline(axh(7), 39.5, 'linewidth', 2);
        try 
%             yline(axh(8), 39.5, 'linewidth', 2); % no need to have
%             reference cause I don't know the right value. 
        catch
        end
        
        xlim(axh(1), [-0.1 1.36]);
        ylim(axh(7), [0 200]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(7), 'origin');
        try
            title(axh(8), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'df/dx');
%         saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoD_cell{fce_i,dist_i} = psudoD_mat;
%         close all;
    end
    
end

%%%%%%%%%%%%%%%%%
figure(); 
clear axh; 
dat_mean_cc = zeros(3,3);
dat_std_cc  = zeros(3,3);
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        psudoD_Mat = psudoD_cell{fce_i,dist_i};
        plot(psudoD_Mat); 
        dat_mean_cc(fce_i,dist_i) = mean(psudoD_Mat(:));
        dat_std_cc(fce_i,dist_i)  = std(psudoD_Mat(:));
        xlabel('pert time'); 
        ylabel('damping estimation');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('damping estimation using dF/dv');
linkaxes(axh);

% better figure; 
peak_time = [0.1:0.025:0.25]; % s... Need to change here! 
fh = figure('unit', 'inch', 'position', [0 0 5 5]); 
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoD_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        xlabel('pert time'); 
        ylabel('D (Ns/m)');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('damping estimation using dF/dv, Subject');
linkaxes(axh);
ylim([-10 100]);

%%corss condition plot 
fh = figure('unit', 'inch', 'position', [0 0 3 3]); 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc(:,dist_i), dat_std_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([-10 100]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
% legend(lnh, {'640N/m', '320N/m', '160N/m'});
legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
ylabel('damping (Ns/m)');
% linkaxes(axh(:), 'xy');
title('damping cross conditions (subject)');  

%% a question needs to be answer (adapt code here)... The code answer 2nd of these questions
% 1. use the x and f, do the dF/dx to show how the stiffness works here 
% 2. use a, v and f, do the dF/dv show how the damping works here
% (acceleration in this code block)
% inertia measurement (mass) here, use dF/da 

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
% psudoK_cell = cell(3,3);
psudoM_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
%         psudoK_mat = zeros(5,7); % 
        psudoM_mat = zeros(7,7);
        for pi = 1:8%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); 
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                
                
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, a_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
%             plot(t_grids, v_avg);
%             v_settled = nanmean(v_avg(t_grids>1.0 & t_grids<1.1));
            plot(t_grids, a_avg);
            a_settled = nanmean(a_avg(t_grids>1.0 & t_grids<1.1));
%             plot(t_grids, x_avg - x_avg(1));
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
            axh(7) = subplot(4,2,7);
            plot(t_grids, f_avg ./ (a_avg - a_settled));
%             ylim([0 2000]);
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
%                 x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
%                 plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                 v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                a_shift = mean(diff(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)))./diff(celltmp1{ti,pi}.t(find(t>-0.1 & t<0))));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                plot(t(2:end), diff(celltmp1{ti,pi}.v(2,:))./diff(celltmp1{ti,pi}.t), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
%                 plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                plot(t_grids, a_dat - a_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; 
%                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                f_net = f_dat - f_avg;
                m_est = -(f_dat - f_avg)./(a_dat - a_avg);
                m_est_idx = find(f_net == max(f_net));
                psudoM_Mat(pi-1,ti) = m_est(m_est_idx);
                plot(t_grids, m_est, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted x
%                 plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
%                 plot(t_grids(pert_idx), v_dat(pert_idx) - v_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
                plot(t_grids(pert_idx), a_dat(pert_idx) - a_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(8) = subplot(4,2,8); % mass estimation
                plot(t_grids(m_est_idx), m_est(m_est_idx), 'marker', 'o', 'markerSize', 5);
                
% %                 % Find and plot on the original data point 
% %                 pert0idx  = find(abs(fp_dat) > 0.01);
% %                 pert0idx = pert0idx(1);  
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
% %                 
% %                 x0 = x_dat(pert0idx) - x_avg(pert0idx);
% %                 f0 = f_dat(pert0idx) - f_avg(pert0idx);
% %                 
% %                 % Find and plot on the peak data point 
% %                 [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
% %                 
% %                 x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
% %                 f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
% %                 
% %                 psudoK = (f1-f0)/(x1-x0);
% %                 psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            try
                yline(axh(8), 3.15, 'linewidth', 2);
                linkaxes(axh(7:8), 'y');
            catch
            end
            ylim(axh(7), [-20 20]);
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'df/dx');
%         saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
%         psudoK_cell{fce_i,dist_i} = -psudoK_mat;
%         close all;
        psudoM_cell{fce_i,dist_i} = psudoM_Mat;
    end
    
end

%% %%%%%%%%%%%%%%%
figure(); 
clear axh; 
dat_mean_cc = zeros(3,3);  % cross condition
dat_std_cc = zeros(3,3); 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        psudoM_Mat = psudoM_cell{fce_i,dist_i};
        dat_mean_cc(fce_i, dist_i) = mean(psudoM_Mat(:)); 
        dat_std_cc(fce_i, dist_i) = std(psudoM_Mat(:)); 
        plot(psudoM_Mat); 
        xlabel('pert time'); 
        ylabel('mass estimation');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('mass estimation using dF/da');
linkaxes(axh);


peak_time = [0.1:0.025:0.25]; % s
fh = figure('unit', 'inch', 'position', [0 0 5 5]); 
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoM_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        xlabel('pert time'); 
        ylabel(' I (kg)');
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
    end
end
sgtitle('Subject mass after release, dF/da');
linkaxes(axh);
ylim([0 20]);

fh = figure('unit', 'inch', 'position', [0 0 3 3]); 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc(:,dist_i), dat_std_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([0 10]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
legend(lnh, {'2.5cm', '5cm', '7.5cm'});
ylabel('mass (kg)');
% linkaxes(axh(:), 'xy');
title('mass cross conditions (subject)');  



%% For the spring measurement data   
% plot: 4-row * 2 -col 
% stiffness measurement here, use dF/dx
clear; clc;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
k_stfcoef = 13/20;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1:3%1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:13),5,13); % 1 no -ert and 5 pert
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
%                 x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;  % this f_avg might not be right as it is 'centralized' after the release.
            
            % A TRICK TO UPDATE F_AVG HERE
            if_favgupdate = 1;
            if (if_favgupdate)
                f_avg_br_pert = 0;  % the force before released.
                t_tmp = 0;          % the number of trials
                for tti = 1:size(celltmp1,1)
                    if isempty(celltmp1{tti,pi}) || pi == 1
                        continue;
                    end
                    idx_release = find(celltmp1{tti,pi}.ts == 5);
                    t = celltmp1{tti,pi}.t - celltmp1{tti,pi}.t(idx_release(1));
                    idx_tmp = find(t>=-0.1 & t<0);
                    f_tmp = mean(celltmp1{tti,pi}.f(2,idx_tmp));
                    t_tmp = t_tmp + 1;
                    f_avg_br_pert = f_avg_br_pert + f_tmp;
                    celltmp1{tti,pi}.tshift = t;
                end
                f_avg_br_pert = f_avg_br_pert/t_tmp;
                f_avg_bfr = mean(f_avg(t_grids > -0.1 & t_grids<0));    % before release
                f_diff = f_avg_br_pert - f_avg_bfr;
                f_avg_upd = f_avg + f_diff;         % the before-release value are same now
                
                if (ifplot)
                    clf;
                    hold on;
                    plot(t_grids, f_avg_upd, 'r.');
                    plot(t_grids, f_avg, 'b.');
                    for tti = 1:size(celltmp1,1)
                        if isempty(celltmp1{tti,pi}) || pi == 1
                            continue;
                        end
                        plot(celltmp1{tti,pi}.tshift, celltmp1{tti,pi}.f(2,:), 'Color', [0.5 0.5 0.5]);
                    end
                    
                end
                if (pi~=1)
                    f_avg = f_avg_upd;
                end
            end
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
            plot(t_grids, x_avg);
%             plot(t_grids, v_avg);
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
 
            axh(7) = subplot(4,2,7);
            k_avgest = f_avg ./ (x_avg - x_avg(end));
            plot(t_grids, -k_avgest);
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                     plot(t, celltmp1{ti,pi}.ox(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
%                 x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; % subtracted dF/dx
                plot(t_grids, -(f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                
                k_est = -(f_dat - f_avg)./(x_dat - x_avg); 
                x_net = x_dat - x_avg; 
                [~, x_net_idx] = min(x_net); 
                k_est_pt = k_est(x_net_idx); 
                plot(t_grids(x_net_idx), k_est_pt, 'marker', 'o', 'markersize', 10); 
                psudoK_mat(pi-1, ti) = k_est_pt;
                

                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted x
                plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                ylim([-0.05 0.04]);
%                 plot(t_grids(pert_idx), v_dat(pert_idx) - v_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
                scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
                axh(6) = subplot(4,2,6); % subtracted F
                scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
                
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
                scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
                axh(6) = subplot(4,2,6); % subtracted F
                scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
%                 psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'K (N/m)'); subplot(axh(7));ylim(+[0 K_list(dist_i)*1.5]);
        yline(K_list(dist_i)*k_stfcoef);
        yline(K_list(dist_i));
        if (pi~=1)
            ylabel(axh(8), 'K (N/m)'); subplot(axh(8));ylim(+[0 K_list(dist_i)*1.5]);
            yline(K_list(dist_i)*k_stfcoef); yline(K_list(dist_i));
        end
%           saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness300ms' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
%         close all;
    end
    
end

%%%%%%%%%%%%%%%%%
figure(); 
clear axh; 
dat_mean_cc = zeros(3,3);
dat_std_cc  = zeros(3,3); 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        psudoK_Mat = -psudoK_cell{fce_i,dist_i};
        dat_mean_cc(fce_i, dist_i) = mean(psudoK_Mat(:));
        dat_std_cc(fce_i ,dist_i)  = std(psudoK_Mat(:));
        plot(psudoK_Mat); 
        xlabel('pert time'); 
        ylabel('damping estimation');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('stiffness estimation using dF/dx, spring');
linkaxes(axh);
ylim([-200 1000]);


%%%%%%%%%%%%%%%%%%%
peak_time = [0.15:0.1:0.55]; % s
fh = figure('unit', 'inch', 'position', [0 0 5 5]); 
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = -psudoK_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        xlabel('pert time'); 
        ylabel(' K (N/m)');
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
    end
end
sgtitle('Spring stiffness after release, dF/dx');
linkaxes(axh);
ylim([-200 1000]);


% corss condition plot 
fh = figure('unit', 'inch', 'position', [0 0 3 3]); 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc(:,dist_i), dat_std_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([-200 1000]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
legend(lnh, {'640N/m', '320N/m', '160N/m'});
% legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
ylabel('stiffness (N/m)');
% linkaxes(axh(:), 'xy');
title('stiffness cross conditions (spring)');  


%% For the spring measurement data   
% plot: 4-row * 2 -col 
% damping measurement (viscosity) here, use dF/dv
clear; clc;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
k_stfcoef = 13/20;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoD_cell = cell(3,3);
for fce_i = 1:3%1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
        psudoD_mat = zeros(5,7); % 
        for pi = 1:6%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
%             celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:13),5,13); % 1 no -ert and 5 pert
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
%                 x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;  % this f_avg might not be right as it is 'centralized' after the release.
            
            % A TRICK TO UPDATE F_AVG HERE
            if_favgupdate = 1;
            if (if_favgupdate)
                f_avg_br_pert = 0;  % the force before released.
                t_tmp = 0;          % the number of trials
                for tti = 1:size(celltmp1,1)
                    if isempty(celltmp1{tti,pi}) || pi == 1
                        continue;
                    end
                    idx_release = find(celltmp1{tti,pi}.ts == 5);
                    t = celltmp1{tti,pi}.t - celltmp1{tti,pi}.t(idx_release(1));
                    idx_tmp = find(t>=-0.1 & t<0);
                    f_tmp = mean(celltmp1{tti,pi}.f(2,idx_tmp));
                    t_tmp = t_tmp + 1;
                    f_avg_br_pert = f_avg_br_pert + f_tmp;
                    celltmp1{tti,pi}.tshift = t;
                end
                f_avg_br_pert = f_avg_br_pert/t_tmp;
                f_avg_bfr = mean(f_avg(t_grids > -0.1 & t_grids<0));    % before release
                f_diff = f_avg_br_pert - f_avg_bfr;
                f_avg_upd = f_avg + f_diff;         % the before-release value are same now
                
                if (ifplot)
                    clf;
                    hold on;
                    plot(t_grids, f_avg_upd, 'r.');
                    plot(t_grids, f_avg, 'b.');
                    for tti = 1:size(celltmp1,1)
                        if isempty(celltmp1{tti,pi}) || pi == 1
                            continue;
                        end
                        plot(celltmp1{tti,pi}.tshift, celltmp1{tti,pi}.f(2,:), 'Color', [0.5 0.5 0.5]);
                    end
                    
                end
                if (pi~=1)
                    f_avg = f_avg_upd;
                end
            end
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
%             plot(t_grids, x_avg);
            plot(t_grids, v_avg);
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
 
            axh(7) = subplot(4,2,7);
%             k_avgest = f_avg ./ (x_avg - x_avg(end));
            d_avgest = f_avg ./ (v_avg - v_avg(end));
            plot(t_grids, -d_avgest);
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
%                 plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                     plot(t, celltmp1{ti,pi}.ox(2,:), 'color', color_arr(4+dist_i,:));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
%                 x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
%                 plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; % subtracted dF/dx
                plot(t_grids, -(f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                
                d_est = (f_dat - f_avg)./(v_dat - v_avg); 
                v_net = v_dat - v_avg; 
                [~, v_net_idx] = max(v_net); 
                d_est_pt = d_est(v_net_idx); 
                plot(t_grids(v_net_idx), d_est_pt, 'marker', 'o', 'markersize', 10); 
                psudoD_mat(pi-1, ti) = d_est_pt;
                

                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted v
%                 plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
%                 ylim([-0.05 0.04]);
                plot(t_grids(pert_idx), v_dat(pert_idx) - v_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
%                 scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
                axh(6) = subplot(4,2,6); % subtracted F
%                 scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
                
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                %axh(2) = subplot(3,2,2);  % subtracted Fp
                %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
                axh(4) = subplot(4,2,4); % subtracted x
%                 scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
                axh(6) = subplot(4,2,6); % subtracted F
%                 scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
%                 psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            
            % plot the subtracted position and force, in another
            % figure/panel
        
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'K (N/m)'); subplot(axh(7));ylim(+[0 K_list(dist_i)*1.5]);
%         yline(K_list(dist_i)*k_stfcoef);
%         yline(K_list(dist_i));
        if (pi~=1)
            ylabel(axh(8), 'K (N/m)'); subplot(axh(8));ylim(+[0 K_list(dist_i)*1.5]);
            yline(K_list(dist_i)*k_stfcoef); yline(K_list(dist_i));
        end
%           saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness300ms' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoD_cell{fce_i,dist_i} = -psudoD_mat;
%         close all;
    end
    
end

%%%%%%%%%%%%%%%%%
figure(); 
clear axh; 
dat_mean_cc = zeros(3,3);
dat_std_cc = zeros(3,3);
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        psudoD_Mat = -psudoD_cell{fce_i,dist_i};
        plot(psudoD_Mat); 
        dat_mean_cc(fce_i,dist_i) = mean(psudoD_Mat(:));
        dat_std_cc(fce_i,dist_i)  = std(psudoD_Mat(:));
        xlabel('pert time'); 
        ylabel('damping estimation');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('damping estimation using dF/dv, spring');
linkaxes(axh);
ylim([-10 10]);

% better figure; 
peak_time = [0.15:0.1:0.55]; % s
fh = figure('unit', 'inch', 'position', [0 0 5 5]); 
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = -psudoD_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        xlabel('pert time'); 
        ylabel('damping estimation');
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
    end
end
sgtitle('damping estimation using dF/dv, spring');
linkaxes(axh);
ylim([-10 100]);


% corss condition plot 
fh = figure('unit', 'inch', 'position', [0 0 3 3]); 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc(:,dist_i), dat_std_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([-10 100]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
legend(lnh, {'640N/m', '320N/m', '160N/m'});
% legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
ylabel('damping (Ns/m)');
% linkaxes(axh(:), 'xy');
title('damping cross conditions (spring)');  

%% a question needs to be answer (adapt code here)... The code answer 2nd of these questions
% 1. use the x and f, do the dF/dx to show how the stiffness works here 
% 2. use a, v and f, do the dF/dv show how the damping works here
% (acceleration in this code block)
% inertia measurement (mass) here, use dF/da 
clear; clc; close all; 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
% psudoK_cell = cell(3,3);
psudoM_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        %         figure; hold on;
%         psudoK_mat = zeros(5,12); % 
%         psudoK_mat = zeros(5,7); % 
        psudoM_mat = zeros(7,7);
        for pi = 1:6%13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); 
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                
                
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, a_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
%                 % also, plot out the origin
%                 axh(1) = subplot(3,2,1); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(3) = subplot(3,2,3); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
%                 axh(5) = subplot(3,2,5); hold on;
%                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            % plot out the avg
%             fh(pi,2) = figure();
            axh(1) = subplot(4,2,1); 
            plot(t_grids, fp_avg);
            axh(3) = subplot(4,2,3);
%             plot(t_grids, v_avg);
%             v_settled = nanmean(v_avg(t_grids>1.0 & t_grids<1.1));
            plot(t_grids, a_avg);
            a_settled = nanmean(a_avg(t_grids>1.0 & t_grids<1.1));
%             plot(t_grids, x_avg - x_avg(1));
            axh(5) = subplot(4,2,5);
            plot(t_grids, f_avg);
            axh(7) = subplot(4,2,7);
            plot(t_grids, f_avg ./ (a_avg - a_settled));
%             ylim([0 2000]);
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
%                 x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
%                 plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                 v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                a_shift = mean(diff(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)))./diff(celltmp1{ti,pi}.t(find(t>-0.1 & t<0))));
%                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                plot(t(2:end), diff(celltmp1{ti,pi}.v(2,:))./diff(celltmp1{ti,pi}.t), 'color', color_arr(4+dist_i,:));
%                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                % plot the subtraction in other panels 
                axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                axh(4) = subplot(4,2,4); hold on;% subtracted x
%                 plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                plot(t_grids, a_dat - a_avg, 'color', color_arr(4+dist_i,:));
                axh(6) = subplot(4,2,6); hold on;% subtracted F
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:));
                axh(8) = subplot(4,2,8); hold on; 
%                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                f_net = f_dat - f_avg;
                a_net = a_dat - a_avg; 
                [~, fp_peak_idx] = max(abs(fp_dat)); 
                a_net_tmp = a_net; a_net_tmp(1:fp_peak_idx) = 0; 
                % pick the first peak from the a_net_tmp 
                [PKS,LOCS]= findpeaks(abs(a_net_tmp)); 
                a_net_peakidx = LOCS(1); 
                m_est = -(f_dat - f_avg)./(a_dat - a_avg);
%                 m_est_idx = find(f_net == max(f_net));
                m_est_idx = a_net_peakidx;
                psudoM_Mat(pi-1,ti) = m_est(m_est_idx);
                plot(t_grids, m_est, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
                % plot out the perturbed position 
                axh(2) = subplot(4,2,2);  % subtracted Fp
                plot(t_grids(pert_idx), fp_dat(pert_idx) - fp_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(4) = subplot(4,2,4); % subtracted x
%                 plot(t_grids(pert_idx), x_dat(pert_idx) - x_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
%                 plot(t_grids(pert_idx), v_dat(pert_idx) - v_avg(pert_idx),...
%                     'color', color_arr(4+dist_i,:), ...
%                     'LineWidth', 2);
                plot(t_grids(pert_idx), a_dat(pert_idx) - a_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(6) = subplot(4,2,6); % subtracted F
                plot(t_grids(pert_idx), f_dat(pert_idx) - f_avg(pert_idx),...
                    'color', color_arr(4+dist_i,:), ...
                    'LineWidth', 2);
                axh(8) = subplot(4,2,8); % mass estimation
                plot(t_grids(m_est_idx), m_est(m_est_idx), 'marker', 'o', 'markerSize', 5);
                
% %                 % Find and plot on the original data point 
% %                 pert0idx  = find(abs(fp_dat) > 0.01);
% %                 pert0idx = pert0idx(1);  
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert0idx), x_dat(pert0idx) - x_avg(pert0idx), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert0idx), f_dat(pert0idx) - f_avg(pert0idx), 10);
% %                 
% %                 x0 = x_dat(pert0idx) - x_avg(pert0idx);
% %                 f0 = f_dat(pert0idx) - f_avg(pert0idx);
% %                 
% %                 % Find and plot on the peak data point 
% %                 [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
% %                 %axh(2) = subplot(3,2,2);  % subtracted Fp
% %                 %scatter(t_grids(pert0idx), fp_dat(pert0idx) - fp_avg(pert0idx), 10);
% %                 axh(4) = subplot(3,2,4); % subtracted x
% %                 scatter(t_grids(pert_idx(locsx(1))), x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1))), 10);
% %                 axh(6) = subplot(3,2,6); % subtracted F
% %                 scatter(t_grids(pert_idx(locsf(1))), f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1))), 10);
% %                 
% %                 x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
% %                 f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
% %                 
% %                 psudoK = (f1-f0)/(x1-x0);
% %                 psudoK_mat(ti,pi-1) = psudoK;
            end
            linkaxes(axh, 'x');
            try
                yline(axh(8), 3.15, 'linewidth', 2);
                linkaxes(axh(7:8), 'y');
            catch
            end
            ylim(axh(7), [-20 20]);
        
        % plot notes here: 
        xlim(axh(1), [-0.1 1.36]);
        sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]); 
        title(axh(1), 'origin');
        try
            title(axh(2), 'subtracted avg');
        catch 
        end
        ylabel(axh(1), 'Fp');
        ylabel(axh(3), 'x');
        ylabel(axh(5), 'f');
        ylabel(axh(7), 'df/dx');
%         saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
%         psudoK_cell{fce_i,dist_i} = -psudoK_mat;
%         close all;
        psudoM_cell{fce_i,dist_i} = psudoM_Mat;
    end
    
end

%%%%%%%%%%%%%%%%%
figure(); 
clear axh;
dat_mean_cc = zeros(3,3);
dat_std_cc = zeros(3,3); 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        psudoM_Mat = psudoM_cell{fce_i,dist_i};
        plot(psudoM_Mat); 
        dat_mean_cc(fce_i, dist_i) = mean(psudoM_Mat(:)); 
        dat_std_cc(fce_i, dist_i)  = std(psudoM_Mat(:)); 
        xlabel('pert time'); 
        ylabel('mass estimation');
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('mass estimation using dF/da');
linkaxes(axh);


peak_time = [0.15:0.1:0.55]; % s
fh = figure('unit', 'inch', 'position', [0 0 5 5]); 
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoM_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        xlabel('pert time'); 
        ylabel(' I (kg)');
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
    end
end
sgtitle('Spring mass after release, dF/da');
linkaxes(axh, 'xy');
ylim([0 10]);


fh = figure('unit', 'inch', 'position', [0 0 3 3]); 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc(:,dist_i), dat_std_cc(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
ylim([0 10]);
xlim([13 27])

% title(['fce' num2str(fce_i) 'dist' num2str(dist_i)]);
xlabel('force threshold');
legend(lnh, {'640N/m', '320N/m', '160N/m'});
ylabel('mass (kg)');
% linkaxes(axh(:), 'xy');
title('mass cross conditions (spring)');  
