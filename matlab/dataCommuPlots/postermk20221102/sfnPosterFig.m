% Figure 3 plot an example of data 
clear; 

close all; clc;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4379_4438.mat')

subj = 5; 

triali = 2;
% plot single trial 
trtmp = data{subj,1,2,2,triali,1};
figure('unit', 'inch', 'position', [0 0 8 3]);
axh(1) = subplot(2,1,1); 
plot(trtmp.t, trtmp.f(1,:), 'linewidth', 2);
xlabel('t (s)'); ylabel('F (N)');
title('force'); 

axh(2) = subplot(2,1,2); 
plot(trtmp.t, trtmp.ox(1,:,1), 'linewidth', 2);
xlabel('t (s)'); ylabel('x (m)');
title('displacement');

linkaxes(axh, 'x');
xlim([-2.5 1])
sgtitle('data of ballistic release')

%% plot the raw data across trials 

%% Figure 5.
% plot an exampel data for a bunch of trials, emphasize the raw data and
% the prediction
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4379_4438.mat'); % new data with 5 subject 
figure('unit', 'centimeters', 'position', [0,0,16,14])
f_minn = -5;
f_maxx = 25;
d_minn = -0.01;
d_maxx = 0.05;
t_offset = 0.5; 
t_minn = -0.25; 
t_maxx = 0.75;
% for subj_i = 5%1:r
subj_i = 5;%1:r
%     for dir_i = 1 % 1:c
dir_i = 1; % 1:c
results = Results(subj_i,dir_i);
% Unperturbed
%         figure(),
set(gcf,'color','w');
axh(1) = subplot(2,1,1); plot(results.FD_UP{2, 1}{1, 1}- t_offset,results.FD_UP{2, 1}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1.5), hold on,
plot(results.avg_FD_UP{2,1}(1,:) - t_offset,results.avg_FD_UP{2,1}(2,:),'r','LineWidth',3), grid on
ylim([f_minn f_maxx])
ylabel('Force (N)');
xticks([-0.2 0 0.5]);
lgd = legend({'individual trials','','','','','','','','','trial average'},'Location','best');
lgd.FontSize = 20;
lgd.Box = 'on';
% sgtitle('Force [N] in Unperturbed Ballistic Release')


% figure(),
% set(gcf,'color','w');
axh(2) = subplot(2,1,2); plot(results.FD_UP{2, 1}{1, 1}-t_offset,results.FD_UP{2, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1.5), hold on,
plot(results.avg_FD_UP{2,1}(1,:) - t_offset,results.avg_FD_UP{2,1}(3,:),'r','LineWidth',3), hold on
plot(results.tt_mod{2, 1}-t_offset,results.disp_mod{2, 1},'b--','LineWidth',3),
%Unperturbed 7 Trials
%lgd = legend('','','','','','','','',strcat('FIT: ',num2str(round(results.FIT_up_avg(2,1))),'%'),'Location','best');
%Pulses 20 Trials
lgd = legend({'individual trials','','','','','','','','','trial average', 'model prediction'},'Location','best');
lgd.FontSize = 20;
lgd.Box = 'on';
grid on
ylim([d_minn d_maxx])
ylabel('displacement (m)');
xticks([-0.2 0 0.5]);
linkaxes(axh, 'x');
xlim([t_minn t_maxx]);
sgtitle('model prediction')
%     end
% end

%% figure 6
% The mean and std for different condiitons 
clear dat
%  for subj_i  = 1:4
obj.data = data
cond.subj = 1:size(data,1);
cond.dir  = 1:size(data,2);
cond.fce  = 1:size(data,3);
cond.dist = 1:size(data,4);
cond.trial= 1:size(data,5);
cond.pert = 1:size(data,6);

obj.cond = cond;
subj_i = 5
fh = figure('name', 'EMG', 'unit', 'centimeters', 'position', [0 0 56 20]);
Fs = 500;
col_type_f = [230 137 67;
              221 87 52; 
              139 51 37]/255;
col_type_x = [11 110 64;
              51 165 87;
              120 187 117]/255;
t_range = [-0.2 0.8];
% 
% emg_pair = [1 2 3 4 5 6 7 8];
% emg_pair_label = {'wrist flexor', 'wrist extensor', ...
%     'elbow flexor', 'elbow extensor', ...
%     'anterior deltoid', 'posterior deltiod', ...
%     'pectoralis', 'Trapezius'};
cols = 3;
rows = 2; % only plot position and the muscles
axh = zeros(rows, cols);

dir_i   = 1;
pert_i  = 1;
for dist_i = 1:3
    for fce_i = 1:3
        trials_list = obj.cond.trial;
        % make enough space for data
        dat.row = length(trials_list);
        dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
            obj.data{subj_i,1,fce_i,dist_i,1,1}.t < t_range(2));
        %                     dat.t = t_range(1):1/Fs:t_range(2);
        dat.t = linspace(t_range(1), t_range(2), dat.col);
        dat.pos = nan(dat.row, dat.col);
        dat.fce = nan(dat.row, dat.col);
        dat.emg = nan(8, dat.row, dat.col);

        trial_idx = 0;
        for trial_i = trials_list
            % get the index
            trial_idx = trial_idx + 1;
            % stack the data into matrices
            index_t = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t > t_range(1) & ...
                obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t < t_range(2);
            %                         dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,index_t),dat.t,'spline');
            dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,index_t),dat.t,'spline');
            dat.fce(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.f(1,index_t),dat.t,'spline');
            %                         dat.emg(:,trial_idx,1:sum(index_t))=obj.data{subj_i,dir_i,fce_i,tar_i,trial_i,pert_i}.emg(:,index_t);
            for ch_i = 1:8
                try
                    dat.emg(ch_i,trial_idx,:)=interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,index_t)',dat.t,'spline')';
                catch
                    disp(['no EMG this condition! fce' num2str(fce_i) ' dist' num2str(dist_i)]);
                end
            end
        end

        % plot the mean and average of it

        %
        axh(1,dist_i) = subplot(rows,cols,dist_i); hold on;       % force
        plot(dat.t, mean(dat.fce, 'omitnan'), ...
            'Color', col_type_f(fce_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
%         tmp1 = mean(dat.fce - dat.fce(:,1), 'omitnan')+std(dat.fce, 'omitnan');
%         tmp2 = mean(dat.fce - dat.fce(:,1), 'omitnan')-std(dat.fce, 'omitnan');
        tmp1 = mean(dat.fce, 'omitnan')+std(dat.fce, 'omitnan');
        tmp2 = mean(dat.fce, 'omitnan')-std(dat.fce, 'omitnan');
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type_f(fce_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        if (dist_i == 1)
            ylabel('Force (N)'); 
        end
        %title('force');
        grid on;
        %

        axh(2,dist_i) = subplot(rows,cols,dist_i+3); hold on;       % position
        plot(dat.t, mean(dat.pos - dat.pos(:,1), 'omitnan'), ...
            'Color', col_type_f(fce_i,:), ... col_type_x(dist_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
        tmp1 = mean(dat.pos - dat.pos(:,1), 'omitnan')+std(dat.pos, 'omitnan');
        tmp2 = mean(dat.pos - dat.pos(:,1), 'omitnan')-std(dat.pos, 'omitnan');
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type_f(fce_i,:), ...col_type_x(dist_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        if (dist_i == 1)
            ylabel('Displacement (m)'); 
        end

%         if (dist_i == 2)
%             xlabel('time (s)'); 
%         end
        %title('position');
        grid on;

    end
    linkaxes(axh, 'x');
end

linkaxes(axh(1,:), 'y');
linkaxes(axh(2,:), 'y');
sgtitle('Raw data');
%             end

%% figure 7 
% estimations on single direction

clc, close all

subj = 5%1:6
dir = 1%:4
results = Results(subj,dir);
% tpause = 0.05; %Speed of plotting
% wind_start = 20; %From which time window to start showing

ff = [15 20 25];
xx = [25 50 75];
% ff = {'F1' 'F2' 'F3'};
% xx = categorical({'900' '600' '300'});

widx =23;

%Stiffness 

% figure('name', ['subj' num2str(subj), ', dir' num2str(dir)]),
figure('unit', 'centimeters', 'position', [0 0 50 20]),
subplot(1,3,1)
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
errorbar(x',results.K_up_avg',results.K_up_std','k','linestyle','none', 'linewidth', 3);
grid on
% xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
% title('Ballistic Release')



%Damping

subplot(1,3,2)
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
errorbar(x',results.B_up_avg',results.B_up_std','k','linestyle','none', 'linewidth', 3);
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
% title('Ballistic Release')


%Inertia 

subplot(1,3,3)
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
errorbar(x',results.M_up_avg',results.M_up_std','k','linestyle','none', 'linewidth', 3);
grid on
% xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
% title('Ballistic Release')

sgtitle(['subj' num2str(subj), ', dir' num2str(dir)]);
%     end
% end

%% figure 8
% plot the different direction here 
subj = 5%1:6
% for dir = [1 3 2 4]1:4
% tpause = 0.05; %Speed of plotting
% wind_start = 20; %From which time window to start showing

ff = [15 20 25];
xx = [25 50 75];
% ff = {'F1' 'F2' 'F3'};
% xx = categorical({'900' '600' '300'});

widx =23;
ttitle_str = {'+x', '+y', '-x', '-y'};
%Stiffness 
figure('unit', 'centimeters', 'position', [0 0 60 45]),
dir_dispidx = 0;
for dir = [1 3 2 4]
    results = Results(subj,dir);
dir_dispidx = dir_dispidx + 1;
results = Results(subj,dir);
% figure('name', ['subj' num2str(subj), ', dir' num2str(dir)]),
subplot(2,4,dir_dispidx)
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
errorbar(x',results.K_up_avg',results.K_up_std','k','linestyle','none', 'linewidth', 3);
grid on
% xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 800])
if (dir == 1)
ylabel('Stiffness [N/m]')
end 

if (dir == 4)
    lgd = legend('15 N','20 N','25 N')
end
% legend('S','M','L')
title(ttitle_str(dir));

end 

% %Damping
% 
% subplot(1,3,2)
% set(gcf,'color','w');
% b = bar(xx,results.B_up_avg'); hold on
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(results.B_up_avg');
% % Get the x coordinate of the bars
% x = nan(nbars, ngroups);
% for i = 1:nbars
%     x(i,:) = b(i).XEndPoints;
% end
% % Plot the errorbars
% errorbar(x',results.B_up_avg',results.B_up_std','k','linestyle','none', 'linewidth', 3);
% grid on
% xlabel('Displacement [mm]')
% % xlabel('Target Stiffness [N/m]')
% ylim([0 40])
% ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
% % legend('S','M','L')
% % title('Ballistic Release')
% 

%Inertia 
dir_dispidx = 0;
for dir = [1 3 2 4]
    results = Results(subj,dir);
    dir_dispidx = dir_dispidx + 1;
subplot(2,4,4+dir_dispidx);
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
errorbar(x',results.M_up_avg',results.M_up_std','k','linestyle','none', 'linewidth', 3);
grid on

% xlabel('Target Stiffness [N/m]')
ylim([0 5])


% legend('S','M','L')
% title('Ballistic Release')
if (dir == 1)
ylabel('Mass [kg]')
end 
if (dir == 2)
xlabel('Displacement [mm]')
end 
if (dir == 4)
legend('15 N','20 N','25 N')
end

%     end
% end
end

sgtitle(['Impedance across directions']);