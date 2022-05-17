%% message Accuracy,
% Check if the accuracy of the certain message is accurate enough for the
% judging. 

% for example, plot the raw data and the message data out. If the raw data
% is different from the raw data (in terms of less digits), then I should
% change the experiment setting somehow. 


%% plot the difference in the raw data and formatted data
% use the sessionScan function 
%ssnum = 3629;
ssnum = 3682;
sstmp = SessionScan(ssnum);
sstmp.plotTaskEndpointForce();
xlim([5372 5377]);

%% plot the difference in the raw data nad the intermediate data
% still use the force as an example

fname_itm = sprintf('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate/KingKong.%05d.mat', ssnum);
load(fname_itm)
dat = Data.QL.Data.FORCE_SENSOR_DATA.data;
t = Data.QL.Headers.FORCE_SENSOR_DATA.recv_time;
figure(); 
subplot(3,1,1); 
plot(t, dat(2,:), '-.');
title('Force Sensor Data');
dat = Data.QL.Data.RAW_FORCE_SENSOR_DATA.data;
t = Data.QL.Headers.RAW_FORCE_SENSOR_DATA.recv_time;
subplot(3,1,2); 
plot(t, dat(2,:), '-.');
title('Force Sensor Data');

% see the wam signal 
dat = Data.QL.Data.BURT_STATUS.pos_y(:);
t = Data.QL.Headers.BURT_STATUS.recv_time;
subplot(3,1,3); 
plot(t, dat, '.');
title('Burt x');