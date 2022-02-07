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
        for dist_i = 1:size(data,2) % for each spring 
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
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3793_3795.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 10;     % 1 without pert, and 10 perturbation time
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/changing_Pert_SpringVel_subtract.mp4', 'MPEG-4');
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
       
            
        for dist_i = 1:size(data,2) % for each spring 
            axh(dist_i+1) = subplot(4,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(fce_i,dist_i,:,:),5,pertT_num);
            
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
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3793_3795.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 10;     % 1 without pert, and 10 perturbation time
close all;

for dist_i = 1:size(data,2) % for each spring 
    fh(dist_i) = figure(); hold on;
    for fce_i = 1:size(data,1)

        for pi = 1:pertT_num
            % 1. get the data tobe plotted
            celltmp1 = reshape(data(fce_i,1,:,:),5,pertT_num);
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            celltmp1 = reshape(data(fce_i,dist_i,:,:),5,pertT_num);
            
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
open(v);

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
        writeVideo(v,frame);
end
close(v);  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Visualize for the subject data 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3818_3828.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 12;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
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
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 5 + 1;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/SpringTimeChangingPert_3-by-3_200ms.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

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
        writeVideo(v,frame);
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
%         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        figure; hold on;
        for pi = 1:6%1:length(pertT_unq)
%        fh(pi) = figure(); hold on;
       
%        axh(1) = subplot(4,1,1); hold on;                     % plot PF
       celltmp1 = reshape(data(1,fce_i,1,:,1:6),5,6); % 1 no -ert and 5 pert
%        idx_release = find(celltmp1{1,pi+1}.ts == 5);
%        t = celltmp1{1,pi+1}.t - celltmp1{1,pi+1}.t(idx_release(1));
%        plot(t, celltmp1{1,pi+1}.Fp(2,:), 'color', color_arr(1,:));
       idx_release = find(celltmp1{1,pi}.ts == 5);
       t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
%        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
       %xlim([-0.1 1]); 
        
%             axh(dist_i+1) = subplot(3,1,dist_i+1); hold on;         % plot each response
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),5,6);
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
        saveas(gcf, ['data/processedData/dataDescriptions/ss3873_3884/velF' num2str(F_list(fce_i)) 'K' num2str(K_list(dist_i)) '.png']);
        close all;
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

% the code just duplicate from last segment, and do roughly the same thing on different data
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
%         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        figure; hold on;
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
        saveas(gcf, ['data/processedData/dataDescriptions/ss3896_3905/velF' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) '.png']);
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
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3818_3828.mat', 'data');
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
        psudoK_mat = zeros(5,12); % 
        for pi = 1:13%1:length(pertT_unq)
            fh(pi,1) = figure(); hold on;
            
            axh(1) = subplot(3,1,1); hold on;                     % plot PF
            celltmp1 = reshape(data(1,1,fce_i,dist_i,:,1:13),5,13); % 1 no -ert and 5 pert
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
        saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
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
pt_num = 1 + 7; % 5 perts
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3913_3921.mat', 'data'); pt_num = 1+7; % 7 perturbs
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
for fce_i = 1:size(data,3)
    for dist_i = 1:size(data,4) % for each spring 
%         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        figure; hold on;
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
        saveas(gcf, ['data/processedData/dataDescriptions/ss3913_3921/velF' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) '.png']);
        close all;
    end
%     saveas(gcf, ['data/processedData/dataDescriptions/ss3913_3921/velAll.png']);
end


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
figure()
for session_i = 1:4
    sstmp1 = sstmp(session_i); 
    
for trial_i = 1:min(20,length(sstmp1.trials))
    t_idx = find(sstmp1.trials(trial_i).data.Fp(2,:)~=0 & ...
        (sstmp1.trials(trial_i).data.ts==5)); % hold 
    if isempty(t_idx) 
        continue;
    end
    if (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 6)
        axh((session_i-1)*2 + 1) = subplot(4,2,(session_i-1)*2 + 1); % left one;
        hold on;
    elseif (round(max(abs(sstmp1.trials(trial_i).data.Fp(2,:)))) == 12)
        axh((session_i-1)*2 + 2) = subplot(4,2,(session_i-1)*2 + 2); % left one;
        hold on;
    end
    
    t_shift = sstmp1.trials(trial_i).data.t(t_idx(1));
    sstmp1.trials(trial_i).data.t_shift = sstmp1.trials(trial_i).data.t - t_shift;
    
%     axh(1) = subplot(2,1,1);  grid on; hold on;
    plot(sstmp1.trials(trial_i).data.t_shift, -sstmp1.trials(trial_i).data.Fp(2,:), 'Marker','.', 'Color', 'b');
    % axh(3) = subplot(4,1,2); grid on;
    plot(sstmp1.trials(trial_i).data.t_shift, sstmp1.trials(trial_i).data.f(2,:), 'Marker','.', 'Color', 'r');
    ylim([0 15]);
    % linkaxes(axh([1,3]), 'x')
%     title('pulse of 100ms');
    legend('command y', 'sensor y');
end 
    subplot(4,2,(session_i-1)*2 + 1); grid on; title(['6N max' dur_name{session_i}]);
    subplot(4,2,(session_i-1)*2 + 2); grid on; title(['12N max' dur_name{session_i}]);
end
linkaxes(axh, 'x');
xlim([0 0.4]); % perturb holding 
% sgtitle('during hold');
sgtitle('during move');
