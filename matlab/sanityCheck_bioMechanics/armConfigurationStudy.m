% armConfiguration study

% Study the configuration of the arm depend on the optotrak capture. 
% From the arm configuration to investigate the mass distribution of the
% entire movement, and discuss the possibility of the stiffness changing
% over movement.  

clear; close all; clc;

% 1. Study the configuration on the optotrak capture 

% load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/' ...
%     'sanityCheck_StiffnessMeasurement/pulseParameterSelection/' ... 
%     'ss4229_4231_pulseParameterCompare_.mat'])  % not good! No mark2 & 3

load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/' ...
    'data/processedData/ss4216_4239.mat'])  % not good! No mark2 & 3


%% position in moving rightward
trialtmp = data{1,1,1,3,1,1};
% Figure 1, plot the xyz according to the time. 
fh(1) = figure('unit', 'inch', 'position', [0 0 3 6]); hold on;
t_range = [-0.5, 2];
t_idx = trialtmp.t > t_range(1) & trialtmp.t < t_range(2);
xyz_arr = 'xyz';
for xyz_i = 1:3
    axh{1}(xyz_i) = subplot(3,1,xyz_i); hold on;
    plot(trialtmp.t(t_idx), trialtmp.ox(xyz_i,t_idx,1), 'LineStyle', '-') % marker 1
    plot(trialtmp.t(t_idx), trialtmp.ox(xyz_i,t_idx,2), 'LineStyle', '--') % marker 2
    plot(trialtmp.t(t_idx), trialtmp.ox(xyz_i,t_idx,3), 'LineStyle', '-.') % marker 3
    title(['opt - ' xyz_arr(xyz_i)]);
    xlabel('time (s)'); ylabel('position (m)');
end
linkaxes(axh{1}, 'x');
legend('marker1', 'marker2', 'marker3');

% Figure 2, plot the position in the xyz axis 
t_range = [-0.5, 2];
t_idx = trialtmp.t > t_range(1) & trialtmp.t < t_range(2);
fh(2) = figure('unit', 'inch', 'position', [0 0 3 3]); hold on;
plot3(trialtmp.ox(1,t_idx,1), trialtmp.ox(2,t_idx,1), trialtmp.ox(3,t_idx,1), 'LineStyle', '-');
plot3(trialtmp.ox(1,t_idx,2), trialtmp.ox(2,t_idx,2), trialtmp.ox(3,t_idx,2), 'LineStyle', '--');
plot3(trialtmp.ox(1,t_idx,3), trialtmp.ox(2,t_idx,3), trialtmp.ox(3,t_idx,3), 'LineStyle', '-.');

% Figure 3, plot the limb length 
trialtmp.ll(1,:) = vecnorm(trialtmp.ox(:,:,1) - trialtmp.ox(:,:,2));
trialtmp.ll(2,:) = vecnorm(trialtmp.ox(:,:,2) - trialtmp.ox(:,:,3));
t_range = [-0.5, 2];
t_idx = trialtmp.t > t_range(1) & trialtmp.t < t_range(2);
fh(3) = figure('unit', 'inch', 'position', [0 0 3 3]); hold on;
axh{3}(1) = subplot(2,1,1); 
plot(trialtmp.t(t_idx), trialtmp.ll(1,t_idx), 'Marker', 'o');
legend('link12');
axh{3}(2) = subplot(2,1,2); 
plot(trialtmp.t(t_idx), trialtmp.ll(2,t_idx), 'Marker', '+');
legend('link23');
linkaxes(axh{3}, 'x');

%% position in moving foreward
trialtmp = data{1,2,1,3,1,1};
% Figure 1, plot the xyz according to the time. 
fh(1) = figure('unit', 'inch', 'position', [0 0 3 6]); hold on;
t_range = [-0.5, 2];
t_idx = trialtmp.t > t_range(1) & trialtmp.t < t_range(2);
xyz_arr = 'xyz';
for xyz_i = 1:3
    axh{1}(xyz_i) = subplot(3,1,xyz_i); hold on;
    plot(trialtmp.t(t_idx), trialtmp.ox(xyz_i,t_idx,1), 'LineStyle', '-') % marker 1
    plot(trialtmp.t(t_idx), trialtmp.ox(xyz_i,t_idx,2), 'LineStyle', '--') % marker 2
    plot(trialtmp.t(t_idx), trialtmp.ox(xyz_i,t_idx,3), 'LineStyle', '-.') % marker 3
    title(['opt - ' xyz_arr(xyz_i)]);
    xlabel('time (s)'); ylabel('position (m)');
end
linkaxes(axh{1}, 'x');
legend('marker1', 'marker2', 'marker3');

% Figure 2, plot the position in the xyz axis 
t_range = [-0.5, 2];
t_idx = trialtmp.t > t_range(1) & trialtmp.t < t_range(2);
fh(2) = figure('unit', 'inch', 'position', [0 0 3 3]); hold on;
plot3(trialtmp.ox(1,t_idx,1), trialtmp.ox(2,t_idx,1), trialtmp.ox(3,t_idx,1), 'LineStyle', '-');
plot3(trialtmp.ox(1,t_idx,2), trialtmp.ox(2,t_idx,2), trialtmp.ox(3,t_idx,2), 'LineStyle', '--');
plot3(trialtmp.ox(1,t_idx,3), trialtmp.ox(2,t_idx,3), trialtmp.ox(3,t_idx,3), 'LineStyle', '-.');
% plot3(trialtmp.ox(1,t_idx,1), trialtmp.ox(2,t_idx,1), trialtmp.ox(3,t_idx,1), 'Marker', '*');
% plot3(trialtmp.ox(1,t_idx,2), trialtmp.ox(2,t_idx,2), trialtmp.ox(3,t_idx,2), 'Marker', '*');
% plot3(trialtmp.ox(1,t_idx,3), trialtmp.ox(2,t_idx,3), trialtmp.ox(3,t_idx,3), 'Marker', '*');
xlim([-0.65 -0.35]); ylim([-0.9 -0.4])

% Figure 3, plot the limb length 
trialtmp.ll(1,:) = vecnorm(trialtmp.ox(:,:,1) - trialtmp.ox(:,:,2));
trialtmp.ll(2,:) = vecnorm(trialtmp.ox(:,:,2) - trialtmp.ox(:,:,3));
t_range = [-0.5, 2];
t_idx = trialtmp.t > t_range(1) & trialtmp.t < t_range(2);
fh(3) = figure('unit', 'inch', 'position', [0 0 3 3]); hold on;
axh{3}(1) = subplot(2,1,1); 
plot(trialtmp.t(t_idx), trialtmp.ll(1,t_idx), 'Marker', 'o');
legend('link12');
axh{3}(2) = subplot(2,1,2); 
plot(trialtmp.t(t_idx), trialtmp.ll(2,t_idx), 'Marker', '+');
legend('link23');
linkaxes(axh{3}, 'x');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
... Pending in mass distribution

