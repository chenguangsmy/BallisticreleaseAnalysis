% The force is calibrate along each of the trials. 
% This might be concerning as I cannot make sure whether subject is put his
% hand on the handle or not, during the drifting. 
% 
% Check this by an example datapoint 
% ss4422, tr125 

%% 
clear; clc; close all;
fcetmp = SessionScanFT(4422); 
figure(); 
axh(1) = subplot(4,1,1); hold on; % x
plot(fcetmp.elapse, fcetmp.force_origin(1,:)); 
plot(fcetmp.elapse, fcetmp.force_net(1,:)); 
xlabel('t(s)'); ylabel('x (N)');

axh(2) = subplot(4,1,2); hold on; % y
plot(fcetmp.elapse, fcetmp.force_origin(2,:)); 
plot(fcetmp.elapse, fcetmp.force_net(2,:)); 
xlabel('t(s)'); ylabel('y (N)');

axh(3) = subplot(4,1,3); hold on; % z
plot(fcetmp.elapse, fcetmp.force_origin(3,:)); 
plot(fcetmp.elapse, fcetmp.force_net(3,:)); 
xlabel('t(s)'); ylabel('z (N)');

axh(4) = subplot(4,1,4);  % offset
plot(fcetmp.elapse, fcetmp.force_net - fcetmp.force_origin); 
xlabel('t(s)'); ylabel('force offset (N)');
linkaxes(axh, 'x')
sgtitle('raw force and net force');

% 
sstmp = SessionScan(4422); 
figure();
axh2(1) = subplot(4,1,1); 
plot(sstmp.data.t, sstmp.data.f(1,:)); 
grid on; 
xlabel('t(s)'); ylabel('Force (N)');

axh2(2) = subplot(4,1,2); 
plot(sstmp.data.t, sstmp.data.ox(1,:,1));
grid on; 
xlabel('t(s)'); ylabel('OPT-x (m)');

axh2(3) = subplot(4,1,3); 
plot(sstmp.data.t, sstmp.data.x(1,:));
grid on; 
xlabel('t(s)'); ylabel('robot-x (m)');

axh2(4) = subplot(4,1,4); 
plot(sstmp.data.t, sstmp.data.tq(:,:));s
legend('j1', 'j2', 'j3', 'j4');

sstmp.plotAddTrialMark(axh(2))

grid on; 
xlabel('t(s)'); ylabel('torque (N/m)');

linkaxes(axh2, 'x');
sgtitle('net force and displacment');

% 12346.1894143135 - sstmp.data.t(1) + fcetmp.elapse(1)

%% See if the force still change in the proper force direction.  
% ... todo, find the jacobian matrix of the WAM according to each joint
% position 
% convert the force into the world-axis force. 

figure(); 
plot(sstmp.data.t, sstmp.data.jp)

theta0 = 1.5708;
theta = sstmp.data.jp(4,:) - theta0;
fce_endpt = sstmp.data.f(1:2,:);
theta_mat = zeros(2,2,length(theta));

theta_mat(1,1,:) = cos(theta);
theta_mat(1,2,:) =-sin(theta);
theta_mat(2,1,:) = sin(theta); 
theta_mat(2,2,:) = cos(theta);

fce_globl = zeros(size(fce_endpt));
for i = 1:length(sstmp.data.t)
    fce_globl(:,i) = theta_mat(:,:,i) * fce_endpt(:,i);
end

clear axh;
figure(); hold on; 
axh(1) = subplot(3,1,1); 
plot(sstmp.data.t, theta);
title('rotation angle'); xlabel('t (s)'); ylabel('\Delta\theta (rad)');
axh(2) = subplot(3,1,2); hold on; 
plot(sstmp.data.t, fce_endpt(1,:), '--');  % x
plot(sstmp.data.t, fce_globl(1,:), '.');   % x
title('F(x)'); xlabel('t (s)'); ylabel('Force (N)');
legend('endpoint frame', 'global frame');
axh(3) = subplot(3,1,3); hold on; 
plot(sstmp.data.t, fce_endpt(2,:), '--');  % y
plot(sstmp.data.t, fce_globl(2,:), '.');   % y
title('F(y)'); xlabel('t (s)'); ylabel('Force (N)');
legend('endpoint frame', 'global frame');
linkaxes(axh, 'x');
sgtitle('compare force in endpoint axis and global axis')

% Conclusion: little difference when the force being released. 

%% The problematic trial: 4427
clear; clc; close all;
fcetmp = SessionScanFT(4427); 
figure(); 
axh(1) = subplot(4,1,1); hold on; % x
plot(fcetmp.elapse, fcetmp.force_origin(1,:)); 
plot(fcetmp.elapse, fcetmp.force_net(1,:)); 
xlabel('t(s)'); ylabel('x (N)');

axh(2) = subplot(4,1,2); hold on; % y
plot(fcetmp.elapse, fcetmp.force_origin(2,:)); 
plot(fcetmp.elapse, fcetmp.force_net(2,:)); 
xlabel('t(s)'); ylabel('y (N)');

axh(3) = subplot(4,1,3); hold on; % z
plot(fcetmp.elapse, fcetmp.force_origin(3,:)); 
plot(fcetmp.elapse, fcetmp.force_net(3,:)); 
xlabel('t(s)'); ylabel('z (N)');

axh(4) = subplot(4,1,4);  % offset
plot(fcetmp.elapse, fcetmp.force_net - fcetmp.force_origin); 
xlabel('t(s)'); ylabel('force offset (N)');
linkaxes(axh, 'x')
sgtitle('raw force and net force');

% 
sstmp = SessionScan(4427); 
figure();
axh2(1) = subplot(4,1,1); 
plot(sstmp.data.t, sstmp.data.f(1,:)); 
grid on; 
xlabel('t(s)'); ylabel('Force (N)');

axh2(2) = subplot(4,1,2); 
plot(sstmp.data.t, sstmp.data.ox(1,:,1));
grid on; 
xlabel('t(s)'); ylabel('OPT-x (m)');

axh2(3) = subplot(4,1,3); 
plot(sstmp.data.t, sstmp.data.x(1,:));
grid on; 
xlabel('t(s)'); ylabel('robot-x (m)');

axh2(4) = subplot(4,1,4); 
plot(sstmp.data.t, sstmp.data.tq(:,:));s
legend('j1', 'j2', 'j3', 'j4');

sstmp.plotAddTrialMark(axh(2))

grid on; 
xlabel('t(s)'); ylabel('torque (N/m)');

linkaxes(axh2, 'x');
sgtitle('net force and displacment');

% 12346.1894143135 - sstmp.data.t(1) + fcetmp.elapse(1)

%% See if the force still change in the proper force direction.  
% ... todo, find the jacobian matrix of the WAM according to each joint
% position 
% convert the force into the world-axis force. 

figure(); 
plot(sstmp.data.t, sstmp.data.jp)

theta0 = 1.5708;
theta = sstmp.data.jp(4,:) - theta0;
fce_endpt = sstmp.data.f(1:2,:);
theta_mat = zeros(2,2,length(theta));

theta_mat(1,1,:) = cos(theta);
theta_mat(1,2,:) =-sin(theta);
theta_mat(2,1,:) = sin(theta); 
theta_mat(2,2,:) = cos(theta);

fce_globl = zeros(size(fce_endpt));
for i = 1:length(sstmp.data.t)
    fce_globl(:,i) = theta_mat(:,:,i) * fce_endpt(:,i);
end

clear axh;
figure(); hold on; 
axh(1) = subplot(3,1,1); 
plot(sstmp.data.t, theta);
title('rotation angle'); xlabel('t (s)'); ylabel('\Delta\theta (rad)');
axh(2) = subplot(3,1,2); hold on; 
plot(sstmp.data.t, fce_endpt(1,:), '--');  % x
plot(sstmp.data.t, fce_globl(1,:), '.');   % x
title('F(x)'); xlabel('t (s)'); ylabel('Force (N)');
legend('endpoint frame', 'global frame');
axh(3) = subplot(3,1,3); hold on; 
plot(sstmp.data.t, fce_endpt(2,:), '--');  % y
plot(sstmp.data.t, fce_globl(2,:), '.');   % y
title('F(y)'); xlabel('t (s)'); ylabel('Force (N)');
legend('endpoint frame', 'global frame');
linkaxes(axh, 'x');
sgtitle('compare force in endpoint axis and global axis')

% Conclusion: little difference when the force being released. 