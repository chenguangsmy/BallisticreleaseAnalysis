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
sessions_mat_b = [2457 2458 2459; 2460 2461 2462; 2455 2463 2456];
col_array = colormap(lines);
for session_i = 1:length(sessions_all)
    sessions_idx = sessions_all(session_i);
    eval(['ss' num2str(sessions_idx) ' = SessionScan(' num2str(sessions_idx) ');']);
end

%%  plot stacked figure;
axh1 = figure();
%force_arr = [3 9 15 21];
force_arr = [3 9 15];
for dist = 1:4
%for fce = 1:3
%sessions_allf = sessions_mat_f(:,dist)'; % 3, 9, 15, 21N
%sessions_allb = sessions_mat_b(:,dist)'; % 3, 9, 15, 21N
%sessions_allf = sessions_mat_f(fce,:); % 5 ,7.5, 10cm
sessions_allb = sessions_mat_b(fce,:); % 5 ,7.5, 10cm
%axh = figure();
axh = subplot(1, 3, fce);
for session_i = 1:length(sessions_allb)
    sessions_idx = sessions_allb(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
ylabel(['force:' num2str(force_arr(fce)) 'N']);
if dist == 1
    title('back')
elseif dist == 4
    xlabel('time (s)');
end
%axh = figure();
%axh = subplot(4, 2, (dist-1)*2+2);
% for session_i = 1:length(sessions_allb)
%     sessions_idx = sessions_allb(session_i);
%     eval(['sstmp = ss' num2str(sessions_idx) ';']);
%     
%     plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
% end
% if dist == 1
%     title('back')
% elseif dist == 4
%     xlabel('time (s');
% end
end

axh2 = figure();
tar_arr = [5 7.5 10];
for target = 1:3
sessions_allf = sessions_mat_f(target,:); % 5 ,7.5, 10cm
sessions_allb = sessions_mat_b(target,:); % 5 ,7.5, 10cm
%axh = subplot(3, 2, (target-1)*2+1);
axh = subplot(1, 3, target);
for session_i = 1:length(sessions_allf)
    sessions_idx = sessions_allf(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
%ylabel(['target:' num2str(tar_arr(target)) 'cm']);
title(['target:' num2str(tar_arr(target)) 'cm']);
% if target == 1
%     title('front')
% elseif target == 3
%     xlabel('time (s');
% end
%axh = figure();
% axh = subplot(3, 2, (target-1)*2+2);
% for session_i = 1:length(sessions_allb)
%     sessions_idx = sessions_allb(session_i);
%     eval(['sstmp = ss' num2str(sessions_idx) ';']);
%     
%     plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
% end
% if target == 1
%     title('back')
% elseif target == 3
%     xlabel('time (s');
% end
end

%% look the difference between low-impedance robot setting

% see 3N series
ss2468 = SessionScan(2468);
ss2457 = SessionScan(2457);
ss2458 = SessionScan(2458);
ss2459 = SessionScan(2459);
sessions_all = [2457, 2458, 2459];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);

[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2457, axh2, col_array(1,:));
[axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2458, axh, col_array(2,:));
[axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2459, axh, col_array(3,:));
legend([lnh1, lnh2, lnh3], {'3N 5cm', '3N 7.5cm', '3N 10cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');

% 9N series
ss2469 = SessionScan(2469);
ss2460 = SessionScan(2460);
ss2461 = SessionScan(2461);
ss2462 = SessionScan(2462);
sessions_all = [2460, 2461, 2462];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);
[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2460, axh2, col_array(1,:));
[axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2461, axh, col_array(2,:));
[axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2462, axh, col_array(3,:));
legend([lnh1, lnh2, lnh3], {'9N 5cm', '9N 7.5cm', '9N 10 cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');

% 15N series
% ss2470 = SessionScan(2470);
% ss2456 = SessionScan(2456);
% ss2463 = SessionScan(2463);
% ss2455 = SessionScan(2455);
ss2479 = SessionScan(2479);
ss2480 = SessionScan(2480);
sessions_all = [2456, 2463, 2455];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);
% [axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2456, axh2, col_array(1,:));
% [axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2463, axh, col_array(2,:));
% [axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2455, axh, col_array(3,:));
[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2479, axh2, col_array(1,:));
[axh, val3, lnh2] = plotStepPertResponse_raw_subavg(ss2480, axh, col_array(3,:));
%legend([lnh1, lnh2, lnh3], {'15N 5cm', '15N 7.5cm', '15N 10 cm'});
legend([lnh1, lnh2], {'15N 5cm', '15N 10 cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');
pos_mean = mean([val2{1} val1{1}]); 
neg_mean = mean([val2{2} val1{2}]);
pos_1 = mean([val1{1}]); 
neg_1 = mean([val1{2}]);
val12_var = std(val1{2});
display('');
mean(val1{2}) - mean(val2{2})

% 21N series
ss2471 = SessionScan(2471);
ss2464 = SessionScan(2464);
ss2466 = SessionScan(2466);
ss2467 = SessionScan(2467);
sessions_all = [2464, 2466, 2467];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);
[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2464, axh2, col_array(1,:));
[axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2466, axh, col_array(2,:));
[axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2467, axh, col_array(3,:));
legend([lnh1, lnh2, lnh3], {'21N 5cm', '21N 7.5cm', '21N 10 cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');

%% figure for position against time categorized by target as subplots 
%sessions_all = [[2468, 2457, 2458, 2459]; [2469, 2460, 2461, 2462]; [2470, 2456, 2463, 2455]; [2471, 2464, 2466, 2467]];
sessions_all = [[2481, 2482, 2483]; [2484, 2486, 2487]; [2479, 2488, 2480]; [2489, 2490, 2492]];
sessions_all_arr = sessions_all(:);
for session_i = 1:length(sessions_all_arr)
    sstmp = sessions_all_arr(session_i);
    %['ss' num2str(sstmp) ' = SessionScan(' num2str(sstmp) ');']
    eval(['ss' num2str(sstmp) ' = SessionScan(' num2str(sstmp) ');']);
end
    
force_all = 3:6:21; 
figure();
for fce_i = 1:size(sessions_all, 1)
    axh = subplot(1, size(sessions_all, 1), fce_i); hold on; 
    for dist_i = 1:size(sessions_all, 2)
        sstmp = eval(['ss' num2str(sessions_all(fce_i, dist_i))]);
        axh = plotMeantrialPosPert(sstmp, axh, dist_i);
    end
    %ylim([-0.045, 0.005]);
    ylim([-0.005, 0.045]);
    xlabel('time'); 
    if (fce_i == 1)
        ylabel('perturbation response position');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(force_all(fce_i)) 'N']);
end

% figure for position against time categorized by target as subplots 
%sessions_all = [[2468, 2457, 2458, 2459]; [2469, 2460, 2461, 2462]; [2470, 2456, 2463, 2455]; [2471, 2464, 2466, 2467]]';
sessions_all = [[2481, 2482, 2483]; [2484, 2486, 2487]; [2479, 2488, 2480]; [2489, 2490, 2492]]';
%dist_all = [0.025 0.05, 0.075, 0.10]; 
dist_all = [0.05, 0.075, 0.10];
figure();
for dist_i = 1:size(sessions_all, 1)
    axh = subplot(1, size(sessions_all, 1), dist_i); hold on; 
    for fce_i = 1:size(sessions_all, 2)
        sstmp = eval(['ss' num2str(sessions_all(dist_i, fce_i))]);
        axh = plotMeantrialPosPert(sstmp, axh, fce_i);
    end
    %ylim([-0.045, 0.005]);
    ylim([-0.005, 0.045]);
    xlabel('time'); 
    if (dist_i == 1)
        ylabel('perturbation response position');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(dist_all(dist_i)) 'm']);
end

%%
figure();
for dist_i = 1:size(sessions_all, 1)
    axh = subplot(1, size(sessions_all, 1), dist_i); hold on; 
    for fce_i = 1:size(sessions_all, 2)
        sstmp = eval(['ss' num2str(sessions_all(dist_i, fce_i))]);
        axh = plotMeantrialForcePert(sstmp, axh, fce_i);
    end
    %ylim([-0.045, 0.005]);
    ylim([-8, 3]);
    xlabel('time'); 
    if (dist_i == 1)
        ylabel('perturbation response position');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(dist_all(dist_i)) 'm']);
end
figure();
for fce_i = 1:size(sessions_all, 2)
    axh = subplot(1, size(sessions_all, 2), fce_i); hold on; 
    for dist_i = 1:size(sessions_all, 1)
        sstmp = eval(['ss' num2str(sessions_all(dist_i, fce_i))]);
        axh = plotMeantrialForcePert(sstmp, axh, dist_i);
    end
    %ylim([-0.045, 0.005]);
    ylim([-8, 3]);
    xlabel('time'); 
    if (dist_i == 1)
        ylabel('perturbation response position');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(force_all(fce_i)) 'N']);
end

%% %% Look through the trajectory of Burdet's experiment (using WAM) 
ss2455 = SessionScan(2455);
ss2456 = SessionScan(2456);
%% plot all the trajectory out
axh = figure();
lnh_lgd = zeros(2,1); 
arr_lgd = {'150/m','300/m'};
for trial_i = 1:length(ss2455.trials)
    [axh, lnh] = ss2455.trials(trial_i).plotRobotEndpointTraj(axh, 'b');
    if trial_i == 1
        lnh_lgd(1) = lnh;
    end
end
for trial_i = 1:length(ss2456.trials)
    [axh, lnh] = ss2456.trials(trial_i).plotRobotEndpointTraj(axh, 'r');
    if trial_i == 1
        lnh_lgd(2) = lnh;
    end
end
legend(lnh_lgd, arr_lgd);
title('trajectories of "insteability dyamics"');

%% plot spring test with perturbations 
% 2524, 2527, 2529, 2530
% color 
color_arr = colormap('lines');
% read data
ss2524 = SessionScan(2524); ss2527 = SessionScan(2527); 
ss2529 = SessionScan(2529); ss2530 = SessionScan(2530);
% plot out
figure();
ss2524.plotStepPertResponse_raw(subplot(1,2,1), color_arr(1,:));
ss2524.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(1,:));
figure();
ss2527.plotStepPertResponse_raw(subplot(1,2,1), color_arr(2,:));
ss2527.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(2,:));
figure();
ss2529.plotStepPertResponse_raw(subplot(1,2,1), color_arr(3,:));
ss2529.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(3,:));
figure();
ss2530.plotStepPertResponse_raw(subplot(1,2,1), color_arr(4,:));
ss2530.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(4,:));
% mean subtract overlap
% ??? to be continued...

% ss2536; ss2537; ss2538; ss2539; ss2541; ss2540
ss2536 = SessionScan(2536); ss2537 = SessionScan(2537); ss2538 = SessionScan(2538); 
ss2539 = SessionScan(2539); ss2541 = SessionScan(2541); ss2540 = SessionScan(2540); 
figure();
ss2536.plotStepPertResponse_raw(subplot(1,2,1), color_arr(1,:));
ss2536.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(1,:));
figure();
ss2537.plotStepPertResponse_raw(subplot(1,2,1), color_arr(2,:));
ss2537.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(2,:));
figure();
ss2538.plotStepPertResponse_raw(subplot(1,2,1), color_arr(3,:));
ss2538.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(3,:));
figure();
ss2539.plotStepPertResponse_raw(subplot(1,2,1), color_arr(4,:));
ss2539.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(4,:));
figure();
ss2541.plotStepPertResponse_raw(subplot(1,2,1), color_arr(5,:));
ss2541.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(5,:));
figure();
ss2540.plotStepPertResponse_raw(subplot(1,2,1), color_arr(6,:));
ss2540.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(6,:));