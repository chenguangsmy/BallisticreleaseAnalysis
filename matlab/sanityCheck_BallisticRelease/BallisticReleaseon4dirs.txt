% This file is no-longer being supported because of the old data format
% cannot fit in the new parser. (time synchrony problem)

% four example sesions in the ballistic-release testing. 
clc; clear;

%%
% 2153: subject +y, robot +y; (robot-global)
% 2154: subject -y, robot -y;
% 2155: subject +x, robot -y;
% 2156: subject -x, robot +y;
ss_list = [2153 2154;    % front and back
          2156 2155];   % left and right

ss_list = ss_list(:);%[2153, 2154, 2155, 2156];
for ss_i = 1:length(ss_list)
    eval(['ss' num2str(ss_list(ss_i)) '=SessionScan(' num2str(ss_list(ss_i)) ');']);
end

% see the trial amount used in each experiment 
disp('---------------------------------------------------------------');
disp('show trial SUCESS RATE as trials ');
for ss_i = 1:length(ss_list)
    %trial_tot = ss2153.trials_num;
    ss_tmp = eval(['ss' num2str(ss_list(ss_i))]);
    trial_tot = ss_tmp.trials_num;
    trial_fin = sum([ss_tmp.trials.outcome]==1);
    display(['ss' num2str(ss_list(ss_i)) ', trials: ' num2str(trial_fin) '/' num2str(trial_tot)]);
end

disp('---------------------------------------------------------------');
disp('show trial SUCESS RATE as forms ');
% see the trial amount used in each task condition
for ss_i = 1:length(ss_list)
    ss_tmp = eval(['ss' num2str(ss_list(ss_i))]);
    [sT, tT, sR] = ss_tmp.getConditionalSucessTrials();
    sR_2d = reshape(sR(1,:,:), size(sR, 2), size(sR, 3));
    sR_table = [[ss_tmp.fThs]', sR_2d'];
    % display rate using table
    display(['For session' num2str(ss_list(ss_i))]);
    VarNames = {'Force (N)', 'tar 2.5 (cm)', 'tar 5.0 (cm)', 'tar 7.5 (cm)', 'tar 10.0 (cm)'};
    T = table(sR_table(:,1), sR_table(:,2), sR_table(:,3), sR_table(:,4), sR_table(:,5), 'VariableNames', VarNames);
end

% see the time amount used in each task condition
% ss2155 21N 10cm is significantly longer. (robot insteability).
for ss_i = 1:4
    ss_tmp = eval(['ss' num2str(ss_list(ss_i))]);
    ss_tmp.getConditionaltime();
end

% See the velocity differences when robot configuration is the same
% plot front
axhv1 = eval(['ss' num2str(ss_list(1,1)) '.plotMeantrialVel_sameCond_overlap(1, 1, 1);']);
axhv1 = eval(['ss' num2str(ss_list(2,1)) '.plotMeantrialVel_sameCond_overlap(1, axhv1, 2);']);
suptitle('front (red) and left (green)');
% plot back
axhv1 = eval(['ss' num2str(ss_list(1,2)) '.plotMeantrialVel_sameCond_overlap(-1, 1, 1);']);
axhv1 = eval(['ss' num2str(ss_list(2,2)) '.plotMeantrialVel_sameCond_overlap(-1, axhv1, 2);']);
suptitle('back (blue) and right (cyne)');

