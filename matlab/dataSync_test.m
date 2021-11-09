%tmp2
% Look for the time skewed problem in the code;
% Also looking for its solutions

clc; clf; 
clear;

ss2798 = SessionScan(2798); 

obj = ss2798;
axh(1) = subplot(2,1,1);
plot(obj.wam_t, obj.wamp_h(:,2));
title('after time aligned (fake)');

axh(2) = subplot(2,1,2);
plot(obj.wam.time, obj.wam.tp(:,2));
title('before time aligned (raw)');

linkaxes(axh, 'x');

%% look for the intermediate file
clear; clf; clc;

ssnum = 2894;
fname = sprintf('KingKong.%05d.mat', ssnum);
fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate';
load([fdir '/' fname]);

wam_sdt = Data.QL.Headers.BURT_STATUS.send_time;
wam_cvt = Data.QL.Headers.BURT_STATUS.recv_time;

t_offset(1) = wam_sdt(1); 
t_offset(2) = wam_cvt(1); 

% only look at the time
%axh(1) = subplot(2,1,1);
plot(wam_sdt - t_offset(1), 'b*'); hold on;
title('send time');
%axh(2) = subplot(2,1,2);
plot(wam_cvt - t_offset(2), 'r.');
title('recive time');
legend('send time', 'receieve time');

%linkaxes(axh, 'x');

% see the interval between send and receieve
clf; 
timediff = [wam_cvt - t_offset(2)] - [wam_sdt - t_offset(1)];
subplot(1,2,1); 
plot(timediff, '.');
grid on;
title('time difference (s)');
subplot(1,2,2); 
histogram(timediff, 500)

% see the interval of sending or receieving
figure(); 
subplot(1,2,1); 
plot(1./diff(wam_cvt), '.');
title('receieve time');
subplot(1,2,2);
plot(1./diff(wam_sdt), '.');
title('send time');

%% check the corresponding sample generator
samp_time = Data.QL.Headers.SAMPLE_GENERATED.send_time;
clf; 
plot(1./diff(samp_time), '.');

%% check the force send and receive
ft_sdt = Data.QL.Headers.FORCE_SENSOR_DATA.send_time;
ft_cvt = Data.QL.Headers.FORCE_SENSOR_DATA.recv_time;
t_offset(1) = ft_sdt(1); 
t_offset(2) = ft_cvt(1); 

% see the interval between send and receieve
clf; 
timediff = [ft_cvt - t_offset(2)] - [ft_sdt - t_offset(1)];
subplot(1,2,1); 
plot(timediff, '.');
grid on;
title('time difference (s)');
subplot(1,2,2); 
histogram(timediff, 500)