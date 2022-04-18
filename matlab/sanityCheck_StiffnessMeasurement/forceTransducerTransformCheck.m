% Check the transformation from the force transducer readout to the WAM
% actual position. 

% In this test, I tried a new robot configuration. A new joint- and
% cartesian- position was used. Therefore, a new force-position
% repationship are needed. 

% jp = [0,  0,      0,      1.5708]; %m
% cp = [0,  0.482,  0.516]; %(m)
clear; close all; clc;
ss_tmp = SessionScan(4025)

%% pair the position and the force 
fh = figure();
t_range = [1666, 1686];
idx = [ss_tmp.data.t > t_range(1) & ss_tmp.data.t < t_range(2)];
t = ss_tmp.data.t(:,idx);
x = ss_tmp.data.x(:,idx);
f = ss_tmp.data.f(:,idx);
axh(1) = subplot(2,1,1); % position
plot(t,x - x(:,1));
legend('x', 'y', 'z');
axh(2) = subplot(2,1,2); % force 
plot(t,f); 
legend('x', 'y', 'z');
linkaxes(axh, 'x');