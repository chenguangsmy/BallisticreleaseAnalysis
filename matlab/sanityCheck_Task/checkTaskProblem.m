% Meet a problem of subject BH conducting experiemnt. Several of the
% condition was not conducted in certain blocks. These blocks include: 
% -[x] 4402
%    -[ ] 5.0cm 15N
%    -[ ] 5.0cm 20N
% -[x] 4404
%    -[ ] 7.5cm 15N
%    -[ ] 5.0cm 20N 
% -[x] 4405
%    -[ ] 7.5cm 20N
%    -[ ] 5.0cm 25N
%    -[ ] 7.5cm 25N
%    -[ ] 5.0cm 15N
%    -[ ] 5.0cm 20N

% Now I'm investigating this problem

clear; close all; clc;

ss_list = [4402 4404 4405];
for ss_i = 1:3
load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate/KingKong.0' num2str(ss_list(ss_i)) '.mat'])
load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/KingKongTSync.0' num2str(ss_list(ss_i)) '.mat']);

t_trial = Data.QL.Headers.TRIAL_CONFIG.recv_time; 
dat_trial = Data.QL.Data.TRIAL_CONFIG.trial_no;

t_state = Data.QL.Headers.TASK_STATE_CONFIG.recv_time; 
dat_state = Data.QL.Data.TASK_STATE_CONFIG.id;

t_target = Data.QL.Headers.TASK_STATE_CONFIG.recv_time; 
dat_targetF = Data.QL.Data.TASK_STATE_CONFIG.target(4,:);
dat_targetL = Data.QL.Data.TASK_STATE_CONFIG.target(6,:);

t_ditital = Data.QL.Headers.DIGITAL_EVENT.recv_time;
dat_ditital = Data.QL.Data.DIGITAL_EVENT.data;
% figure(); 
% axh(1) = subplot(5,1,1); 
% plot(t_trial, dat_trial, '*')
% title('trials'); ylabel('trial#');
% 
% axh(2) = subplot(5,1,2);
% plot(t_state, dat_state, '*')
% title('states'); ylabel('states#');
% 
% axh(3) = subplot(5,1,3); 
% plot(t_state, dat_targetF, '*')
% title('targetF'); ylabel('Force(N)');
% 
% axh(4) = subplot(5,1,4); 
% plot(t_state, dat_targetL, '*')
% title('targetL'); ylabel('Distance(m)');
% 
% axh(5) = subplot(5,1,5); 
% plot(t_ditital, dat_ditital, '*')
% title('EVENT'); ylabel('PRESS==4');
% 
% linkaxes(axh, 'x');
% xlabel('time');

% Question: did the pulse being receieved?  
% I hope not cause if pulse recieved but not go advance trial, shows my
% program error. 
% Or else it could because the pulse receieved but not go because of in
% certain state

figure('unit', 'inch', 'position', [0 0 3.5 5]); 
clear axh
axh(1) = subplot(2,1,1);
plot(data.eventsT, data.eventsTrials, '*');
title('trials'); ylabel('trial#');
axh(2) = subplot(2,1,2); 
plot(data.eventsT, data.eventsL, '*');
title('events'); ylabel('data of the 3rd ch');
linkaxes(axh, 'x');
xlabel('t'); 
sgtitle(['ss' num2str(ss_list(ss_i))]);
end