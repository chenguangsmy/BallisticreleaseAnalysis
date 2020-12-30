% superimpose force data on a trial force data part
% author: cleave
% data: 2020-03-11
%%% This is a trial code, inorder to detect data collection quality

trial_i = 2;
force_i = 1;
file_dir = '../rawdata';
fname = 'KingKong.01348';
force_arr = {'Fx','Fy','Fz','Dx','Tx','Ty','Tz'};

% load force file
load([file_dir '/' fname '.mat']);
Table = readtable([file_dir '/' fname '.force.dat']);

% extract force_t, and force
force_idx_bgn = min(trial(trial_i).force_FtSeq); % in trial data
force_idx_end = max(trial(trial_i).force_FtSeq);
hfforce_idx_bgn = find(Table.FT == force_idx_bgn); % in force file
hfforce_idx_end = find(Table.FT == force_idx_end);

force_t = Table.FT(hfforce_idx_bgn:hfforce_idx_end);
Force_var = eval(['Table.' (force_arr{force_i})]);
force = Force_var(hfforce_idx_bgn:hfforce_idx_end);

% plot 
figure;
hold on;
plot(trial(trial_i).force_FtSeq, trial(trial_i).force(force_i,:));
plot(force_t, force);
legend('RTMA\_force', 'hf\_force');
title('comparation force: RTMA and high frequency');