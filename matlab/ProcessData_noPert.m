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
axhv2 = ss2338.plotMeantrialVel_sameCond_overlap(-1, 1, 1);
axhv2 = ss2340.plotMeantrialVel_sameCond_overlap(1, axhv2, 2);
ylim([-0.7 0.7]);
title('movement difference after calibration');
legend('forward', 'backward (invert)');