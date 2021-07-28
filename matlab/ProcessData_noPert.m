% four example sesions in the ballistic-release testing. 
clc; clear;
%%
% 2153: subject +y, robot +y; (robot-global)
% 2154: subject -y, robot -y;
% 2155: subject +x, robot -y;
% 2156: subject -x, robot +y;
ss2153 = SessionScan(2153); 
ss2154 = SessionScan(2154); 
ss2155 = SessionScan(2155); 
ss2156 = SessionScan(2156); 
ss_list = [2153, 2154, 2155, 2156];
% see the trial amount used in each experiment 
for ss_i = 1:4
    %trial_tot = ss2153.trials_num;
    ss_tmp = eval(['ss' num2str(ss_list(ss_i))]);
    trial_tot = ss_tmp.trials_num;
    trial_fin = sum([ss_tmp.trials.outcome]==1);
    display(['ss' num2str(ss_list(ss_i)) ', trials: ' num2str(trial_fin) '/' num2str(trial_tot)]);
end

% see the trial amount used in each task condition
for ss_i = 1:4
    ss_tmp = eval(['ss' num2str(ss_list(ss_i))]);
    [sT, tT, sR] = ss_tmp.getConditionalSucessTrials();
    sR_2d = reshape(sR(1,:,:), size(sR, 2), size(sR, 3));
    sR_table = [[ss_tmp.fThs]', sR_2d'];
    % display rate using table
    display(['For session' num2str(ss_list(ss_i))]);
    VarNames = {'Force (N)', 'tar 2.5 (cm)', 'tar 5.0 (cm)', 'tar 7.5 (cm)', 'tar 10.0 (cm)'};
    T = table(sR_table(:,1), sR_table(:,2), sR_table(:,3), sR_table(:,4), sR_table(:,5), 'VariableNames', VarNames)
end

% see the time amount used in each task condition
% ss2155 21N 10cm is significantly longer. (robot insteability).
for ss_i = 1:4
    ss_tmp = eval(['ss' num2str(ss_list(ss_i))]);
    ss_tmp.getConditionaltime();
end

%% See the velocity differences when robot configuration is the same
axhv1 = ss2153.plotMeantrialVel_sameCond_overlap(1, 1, 1);
axhv1 = ss2156.plotMeantrialVel_sameCond_overlap(1, axhv1, 2);
suptitle('front (red) and left (green)');
axhv2 = ss2154.plotMeantrialVel_sameCond_overlap(-1, 1, 3);
axhv2 = ss2155.plotMeantrialVel_sameCond_overlap(-1, axhv2, 4);
suptitle('back (blue) and right (cyne)');


%% Test non-linearity of robot using spring test 
% 1. No-calibrated. 
ss2333 = SessionScan(2333); 
ss2336 = SessionScan(2336);
axhv1 = ss2333.plotMeantrialVel_sameCond_overlap(1, 1, 1);
axhv1 = ss2336.plotMeantrialVel_sameCond_overlap(-1, axhv1, 2);
ylim([-0.7 0.7]);
title('movement difference before calibration');
legend('forward', 'backward (invert)');
% 2. After-calibration
ss2338 = SessionScan(2338);
ss2340 = SessionScan(2340);
axhv2 = ss2340.plotMeantrialVel_sameCond_overlap(1, 1, 2);
axhv2 = ss2338.plotMeantrialVel_sameCond_overlap(-1, axhv2, 1);
ylim([-0.7 0.7]);
title('movement difference after calibration');
legend('forward', 'backward (invert)');
%% 3. New sessions
ssTestList = [2341, 2342];
ssFront = SessionScan(ssTestList(1));
ssBack = SessionScan(ssTestList(2));
axhv2 = ssFront.plotMeantrialVel_sameCond_overlap(1, 1, 1);
axhv2 = ssBack.plotMeantrialVel_sameCond_overlap(-1, axhv2, 2);
ylim([-0.7 0.7]);
title('movement difference after calibration');
legend('forward', 'backward (invert)');

%% Conduct the perturbation checking 
% front:

%--------------------------------------------------------------------------
% | KingKong.02256.mat	| KingKongFT02256.csv	| KingKongWAM02256.csv	| randomize step pert, 3N 5cm release. |  
% | KingKong.02268.mat	| KingKongFT02268.csv	| KingKongWAM02268.csv	| randomize step pert, 3N 7.5cm release. | 
% | KingKong.02262.mat	| KingKongFT02262.csv	| KingKongWAM02262.csv	| randomize step pert, 3N 10cm release. |   
% | KingKong.02257.mat	| KingKongFT02257.csv	| KingKongWAM02257.csv	| randomize step pert, 9N 5cm release. |  
% | KingKong.02259.mat	| KingKongFT02259.csv	| KingKongWAM02259.csv	| randomize step pert, 9N 7.5cm release. |  
% | KingKong.02266.mat	| KingKongFT02266.csv	| KingKongWAM02266.csv	| randomize step pert, 9N 10cm release. |  
% | KingKong.02263.mat	| KingKongFT02263.csv	| KingKongWAM02263.csv	| randomize step pert, 15N 7.5cm release. | ??? no pert data! 
% | KingKong.02258.mat	| KingKongFT02258.csv	| KingKongWAM02258.csv	| randomize step pert, 15N 5cm release. |  
% | KingKong.02267.mat	| KingKongFT02267.csv	| KingKongWAM02267.csv	| randomize step pert, 15N 10cm release. |  
% | KingKong.02264.mat	| KingKongFT02264.csv	| KingKongWAM02264.csv	| randomize step pert, 21N 5cm release. |  
% | KingKong.02265.mat	| KingKongFT02265.csv	| KingKongWAM02265.csv	| randomize step pert, 21N 7.5cm release. |  
% | KingKong.02261.mat	| KingKongFT02261.csv	| KingKongWAM02261.csv	| randomize step pert, 21N 10cm release. |  
%--------------------------------------------------------------------------
%   5cm: 2256 2257 2258 2264
% 7.5cm: 2268 2259 2263 2265
%  10cm: 2262 2266 2267 2261

% back:
% | KingKong.2355.mat	| KingKongFT2355.csv	| KingKongWAM2355.csv	| backward ballistic release with pert,   5cm   3N |
% | KingKong.2358.mat	| KingKongFT2358.csv	| KingKongWAM2358.csv	| backward ballistic release with pert, 7.5cm   9N |
% | KingKong.2359.mat	| KingKongFT2359.csv	| KingKongWAM2359.csv	| backward ballistic release with pert,  10cm   15N |
% | KingKong.2360.mat	| KingKongFT2360.csv	| KingKongWAM2360.csv	| backward ballistic release with pert,  10cm   21N |
% | KingKong.2361.mat	| KingKongFT2361.csv	| KingKongWAM2361.csv	| backward ballistic release with pert, 7.5cm   21N |
% | KingKong.2362.mat	| KingKongFT2362.csv	| KingKongWAM2362.csv	| backward ballistic release with pert, 7.5cm   15N |
% | KingKong.2363.mat	| KingKongFT2363.csv	| KingKongWAM2363.csv	| backward ballistic release with pert,  10cm   9N |
% | KingKong.2365.mat	| KingKongFT2365.csv	| KingKongWAM2365.csv	| backward ballistic release with pert, 7.5cm   3N |
% | KingKong.2366.mat	| KingKongFT2366.csv	| KingKongWAM2366.csv	| backward ballistic release with pert,   5cm   9N |
% | KingKong.2367.mat	| KingKongFT2367.csv	| KingKongWAM2367.csv	| backward ballistic release with pert,   5cm   15N |
% | KingKong.2368.mat	| KingKongFT2368.csv	| KingKongWAM2368.csv	| backward ballistic release with pert,   5cm   21N |
% | KingKong.2369.mat	| KingKongFT2369.csv	| KingKongWAM2369.csv	| backward ballistic release with pert,  10cm   3N |
%   5cm: 2355 2366 2367 2368
% 7.5cm: 2365 2358 2362 2361
%  10cm: 2369 2363 2359 2360
sessions_mat_f = [2256 2257 2258 2264; 2268 2259 2263 2265; 2262 2266 2267 2261]; % front 
sessions_mat_b = [2355 2366 2367 2368; 2365 2358 2362 2361; 2369 2363 2359 2360]; % back
sessions_all = sessions_mat_b(:);
%sessions_all = [2256 2257 2258 2264 2268 2259 2263 2265 2262 2266 2267 2261];
col_array = colormap(lines);
for session_i = 1:length(sessions_all)
    sessions_idx = sessions_all(session_i);
    eval(['ss' num2str(sessions_idx) ' = SessionScan(' num2str(sessions_idx) ');']);
end

%%  plot stacked figure;
axh1 = figure();
force_arr = [3 9 15 21];
for dist = 1:4
sessions_allf = sessions_mat_f(:,dist)'; % 3, 9, 15, 21N
sessions_allb = sessions_mat_b(:,dist)'; % 3, 9, 15, 21N
%sessions_allf = sessions_mat_f(3,:); % 5 ,7.5, 10cm
%sessions_allb = sessions_mat_b(3,:); % 5 ,7.5, 10cm
%axh = figure();
axh = subplot(4, 2, (dist-1)*2+1);
for session_i = 1:length(sessions_allf)
    sessions_idx = sessions_allf(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
ylabel(['force:' num2str(force_arr(dist)) 'N']);
if dist == 1
    title('front')
elseif dist == 4
    xlabel('time (s');
end
%axh = figure();
axh = subplot(4, 2, (dist-1)*2+2);
for session_i = 1:length(sessions_allb)
    sessions_idx = sessions_allb(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
if dist == 1
    title('back')
elseif dist == 4
    xlabel('time (s');
end
end

axh2 = figure();
tar_arr = [5 7.5 15];
for target = 1:3
sessions_allf = sessions_mat_f(target,:); % 5 ,7.5, 10cm
sessions_allb = sessions_mat_b(target,:); % 5 ,7.5, 10cm
axh = subplot(3, 2, (target-1)*2+1);
for session_i = 1:length(sessions_allf)
    sessions_idx = sessions_allf(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
ylabel(['target:' num2str(tar_arr(target)) 'cm']);
if target == 1
    title('front')
elseif target == 3
    xlabel('time (s');
end
%axh = figure();
axh = subplot(3, 2, (target-1)*2+2);
for session_i = 1:length(sessions_allb)
    sessions_idx = sessions_allb(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
if target == 1
    title('back')
elseif target == 3
    xlabel('time (s');
end
end
