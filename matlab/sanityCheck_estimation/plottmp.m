% this plot is only for showing at the brainbag 20221016

% using data ss4310_4356.mat
% halt at fede_sys_id ln 69

figure(); 
i = 1; 
axh(1) = subplot(2,1,1); 
plot(time_new, force_interp_t(i,:), 'linewidth', 2); 
xlabel('t (s)'); ylabel('force (N)');
grid on; 

axh(2) = subplot(2,1,2); 
plot(time_new, disp_interp_t(i,:), 'linewidth', 2); 
xlabel('t (s)'); ylabel('displacement (m)');
grid on;

sgtitle('Example of force and displacmeent');

%% 
% plot the stiffness, damping, and inertia across directions 
clc, clear, close all
load('data/processedData/ss4310_4356.mat');
figure('name', ['subj' num2str(subj), ', dir' num2str(dir)]),
for subj = 1%:4
    for dir = 1:4
        results = Results(subj,dir);
% tpause = 0.05; %Speed of plotting
% wind_start = 20; %From which time window to start showing

ff = [15 20 25];
xx = [25 50 75];
% ff = {'F1' 'F2' 'F3'};
% xx = categorical({'900' '600' '300'});

widx =23;

%Stiffness 

subplot(1,4,dir)
set(gcf,'color','w');
b = bar(xx,results.K_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.K_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.K_up_avg',results.K_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title(['direction' num2str(dir)])
    end
end


figure('name', ['subj' num2str(subj), ', dir' num2str(dir)]),
for subj = 1%:4
    for dir = 1:4
        results = Results(subj,dir);
%Damping

subplot(1,4,dir)
set(gcf,'color','w');
b = bar(xx,results.B_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.B_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.B_up_avg',results.B_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title(['direction' num2str(dir)])
    end
end


figure('name', ['subj' num2str(subj), ', dir' num2str(dir)]),
for subj = 1%:4
    for dir = 1:4
        results = Results(subj,dir);
%Inertia 

subplot(1,4,dir)
set(gcf,'color','w');
b = bar(xx,results.M_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.M_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.M_up_avg',results.M_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title(['direction' num2str(dir)])

sgtitle(['subj' num2str(subj), ', dir' num2str(dir)]);
    end
end