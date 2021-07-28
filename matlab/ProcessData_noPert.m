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
axhv2 = ssFront.plotMeantrialVel_sameCond_overlap(1, 1, 2);
axhv2 = ssBack.plotMeantrialVel_sameCond_overlap(-1, axhv2, 1);
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
% | KingKong.02263.mat	| KingKongFT02263.csv	| KingKongWAM02263.csv	| randomize step pert, 15N 7.5cm release. |  
% | KingKong.02258.mat	| KingKongFT02258.csv	| KingKongWAM02258.csv	| randomize step pert, 15N 5cm release. |  
% | KingKong.02267.mat	| KingKongFT02267.csv	| KingKongWAM02267.csv	| randomize step pert, 15N 10cm release. |  
% | KingKong.02264.mat	| KingKongFT02264.csv	| KingKongWAM02264.csv	| randomize step pert, 21N 5cm release. |  
% | KingKong.02265.mat	| KingKongFT02265.csv	| KingKongWAM02265.csv	| randomize step pert, 21N 7.5cm release. |  
% | KingKong.02261.mat	| KingKongFT02261.csv	| KingKongWAM02261.csv	| randomize step pert, 21N 10cm release. |  
%--------------------------------------------------------------------------
%   5cm: 2256 2257 2258 2264
% 7.5cm: 2268 2259 2263 2265
%  10cm: 2262 2266 2267 2261
sessions_all = [2256 2257 2258 2264 2268 2259 2263 2265 2262 2266 2267 2261];
for session_i = 1:length(sessions_all)
    sessions_idx = sessions_all(session_i);
    eval(['ss' num2str(sessions_idx) ' = SessionScan(' num2str(sessions_idx) ');']);
end

