% This file check the ballistic-release on 4 directional movements. 
% It takes the both OPT and WAM positions 

% sanity check contents: 
% 1. Compare the y directional OPTOTRAK position and the WAM position 
% 2. Compare the y directional force vs x directional force 







% |sessions:| 2.5cm    | 5cm        | 7.5cm     | 
% | ------- | -------- | ---------  | --------- | 
% | 15N     | 4062     | 4061       | 4060      | 
% | 20N     | 4063     | 4065       | 4064      | 
% | 25N     | 4067     | 4066       | 4068      | 

% four example sesions in the ballistic-release testing. 

%%
% 1. Compare the y directional OPTOTRAK position and the WAM position 
% (ss4060)
sstmp = SessionScan(4060);
celltmp = sstmp.export_as_formatted(1);
celltmp1 = reshape(celltmp(1,1,:,2), size(celltmp,3),1);

%% see the position at the time of release 
figure('unit', 'inch', 'position', [0 0 3 8], ...
        'name', 'release_displacement'); 
hold on;
for trial_i = 1:13 %-valid trials
    axh(1) = subplot(3,1,1);
    hold on;
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.x(1,:), 'b');
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.ox(1,:), 'r');
    ylabel('x (m)');
    
    axh(2) = subplot(3,1,2);
    hold on;
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.x(2,:), 'b');
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.ox(2,:), 'r');
    ylabel('y (m)');
    
    axh(3) = subplot(3,1,3);
    hold on;
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.x(3,:), 'b');
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.ox(3,:), 'r');
    ylabel('z (m)');
end
linkaxes(axh(1:3), 'x');
xlabel('t (s)');
xlim([-0.1 0.5])
sgtitle('WAM(x) vs OPT (x)')

% see the position at the time of perturb
figure('unit', 'inch', 'position', [0 0 3 8], ...
        'name', 'pert_displacement'); 
hold on;
for trial_i = 1:13 %-valid trials
    pert_idx = find(celltmp1{trial_i}.Fp(2,:) ~= 0);
    t_shift = celltmp1{trial_i}.t - celltmp1{trial_i}.t(pert_idx(1));
    axh(1) = subplot(3,1,1);
    hold on;
    plot(t_shift, celltmp1{trial_i}.x(1,:), 'b');
    plot(t_shift, celltmp1{trial_i}.ox(1,:), 'r');
    ylabel('x (m)');
    
    axh(2) = subplot(3,1,2);
    hold on;
    plot(t_shift, celltmp1{trial_i}.x(2,:), 'b');
    plot(t_shift, celltmp1{trial_i}.ox(2,:), 'r');
    ylabel('y (m)');
    
    axh(3) = subplot(3,1,3);
    hold on;
    plot(t_shift, celltmp1{trial_i}.x(3,:), 'b');
    plot(t_shift, celltmp1{trial_i}.ox(3,:), 'r');
    ylabel('z (m)');
end
linkaxes(axh(1:3), 'x');
xlabel('t (s)');
xlim([-0.1 0.5])
sgtitle('WAM(x) vs OPT (x)')

% see the position overlap at the time of perturb
figure('unit', 'inch', 'position', [0 0 3 8], ...
        'name', 'pert_displacement'); 
hold on;
for trial_i = 1:13 %-valid trials
    pert_idx = find(celltmp1{trial_i}.Fp(2,:) ~= 0);
    t_shift = celltmp1{trial_i}.t - celltmp1{trial_i}.t(pert_idx(1));
    x_shift = celltmp1{trial_i}.x(:,:) - celltmp1{trial_i}.x(:,pert_idx(1));
    ox_shift = celltmp1{trial_i}.ox(:,:) - celltmp1{trial_i}.ox(:,pert_idx(1));
    axh(1) = subplot(3,1,1);
    hold on;
    plot(t_shift, x_shift(1,:), 'b');
    plot(t_shift, ox_shift(1,:), 'r');
    ylabel('x (m)');
    
    axh(2) = subplot(3,1,2);
    hold on;
    plot(t_shift, x_shift(2,:), 'b');
    plot(t_shift, ox_shift(2,:), 'r');
    ylabel('y (m)');
    
    axh(3) = subplot(3,1,3);
    hold on;
    plot(t_shift, x_shift(3,:), 'b');
    plot(t_shift, ox_shift(3,:), 'r');
    ylabel('z (m)');
end
linkaxes(axh(1:3), 'x');
xlabel('t (s)');
xlim([-0.1 0.5])
sgtitle('WAM(x) vs OPT (x), subtracted 0')

%%
% 2. Compare the y directional force vs x directional force 
% (ss4065 v.s. ss4069)

%% 4-directional only release data, check the files 
dist_list = [2.5 5.0 7.5];
fce_list = [15 20 25];

% load('data/processedData/ss4216_4220.mat', 'data');       % 3K values - Himanshu
% load('data/processedData/ss4222_4226.mat', 'data');       % only release -Himanshu
load('data/processedData/ss4216_4226.mat', 'data');       % chenguang & Himanshu
Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 
colors = colormap('lines');
close(fh);
r = size(Data, 1); % subj
d = size(Data, 2); % direction
f = size(Data, 3); % force
l = size(Data, 4); % distance
t = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;

epoc_type = 2;  % 1 perturb
                % 2 release
plot_type = 9;  % 1 displacement
                % 2 force 
                % 3 feedforward force
                % 4 velocity
                % 5 torque
                % 6 norm displacement
                % 7 norm force
                % 8 opt displacement
                % 9 opt velocity
                % 10 opt 1
                % 11 opt 2
                % 12 opt 3
pert_type = 1; % choose option [2 3 4]
axh = zeros(d,r);
xyi = 1;        % x, y

% fh = figure();
% pert_type = 2; colors = colors(4:end,:)
for ri = 1:r % subj
    for di = 1:d % direction
        fh = figure('Name', ['direction' num2str(di)]);
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        %         axh(ri, di) = subplot(d,r,r*(di-1) + ri);grid on;hold on;
        %         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
        for fi = 1:f % target force
            for li = 1:3 % target distance
                axh(fi, li) = subplot(f,l,f*(li-1) + fi);grid on;hold on; %?
                for pi = pert_type%1:p % perturbation
                    trial_num = length(Data(ri,di,fi,li,:,pi));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,di,fi,li,ti,pi}))
                            continue;
                        end
                        
                        switch epoc_type
                            case 1
%                                 idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0 & Data{ri,di,fi,li,ti,pi}.ts==4);  % pert at y
                                idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0);
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if pi == 1
                                    disp('ERROR: should use pi == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,di,fi,li,ti,pi}.ts==5 | Data{ri,di,fi,li,ti,pi}.ts==6);
                                idx = (idx(1)-100):idx(end);
                                %idx = (idx(1)):idx(end);
                        end
                        %plot(Data{ri,ci,di,ti,li}.Fp(xyi,:));
                        %idx = find(Data{ri,ci,di,ti,li}.Fp(xyi,:)~=0);
                        %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                        %idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                        %idx = (idx-9):idx(end);
                        %idx = (idx(1)):(idx(end)+100);
                        time = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,di,fi,li,ti,pi}.x(xyi,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,di,fi,li,ti,pi}.f(xyi,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,di,fi,li,ti,pi}.Fp(xyi,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,di,fi,li,ti,pi}.v(xyi,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,di,fi,li,ti,pi}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,di,fi,li,ti,pi}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,di,fi,li,ti,pi}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm force';
                            case 8 % optotrak x
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ox))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                end
                                titlestr = 'opto-position';
                            case 9 % optotrak v
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ov))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx,1);
                                end
                                titlestr = 'opto-velocity';
                            case 10 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                titlestr = 'opto1-position';
                            case 11
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,2);
                                titlestr = 'opto2-position';
                            case 12 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,3);
                                titlestr = 'opto3-position';
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4+li, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
                title(['fce' num2str(fce_list(fi)) 'dist' num2str(dist_list(li))]);

            end
        end
        linkaxes(axh, 'xy');
        xlim([0 1.0])
    end
end
linkaxes(axh, 'xy');
xlim([0 1.0])
% xlim([0 0.5]);
% xlim([0 2])
sgtitle(titlestr);



%% New data running Federico's analysis on the 4-directional data
%%%%%%%%%%%%%%%%%
% Change ln20, to choose plot things: force, displacement, velocity, etc.
clc,  close all, %clear,

% load('ss3938_3949.mat', 'data');  % 6N perturbation on x, 
% load('ss4129_4137.mat', 'data');  % 6N and 12N perturbation on x, only chenguang,
% load('ss4129_4179.mat', 'data');  % 6N and 12N perturbation on x,
% himanshu and chenguang
% load('ss4181_4202.mat', 'data');  % 6N and 12N perturbation on x, himanshu and chenguang
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4216_4220.mat', 'data');  % just release, Chenguang
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4222_4226.mat', 'data');  % just release, Himanshu
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4216_4226.mat', 'data');

Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 

colors = colormap('lines');
r = size(Data, 1); % subj
c = size(Data, 2); % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;
epoc_type = 2;
plot_type = 1;          % 1 displacement
                        % 2 force
                        % 3 force command
                        % 4 velocity
                        % 5 torque of 3rd joint
                        % 6 displacement (vector length in cartesian space)
                        % 7 force (vector length in cartesian space)
axh = zeros(f,r);



for ri = 1:r % subj
    for ci = 1:c
        for fi = 1:f % fce
            for di = 1:d % target distance
                axh(ri, fi) = subplot(d,f,d*(fi-1) + di);grid on;hold on;
                %for di = 3 % target distance
                for li = 1:p % perturbation
                    trial_num = length(Data(ri,ci,fi,di,:,li));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,ci,fi,di,ti,li}))
                            continue;
                        end

                        switch epoc_type
                            case 1
                                idx = find(Data{ri,ci,fi,di,ti,li}.Fp(2,:)~=0 & ...
                                    Data{ri,ci,fi,di,ti,li}.ts==4);  % pert at y
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if li == 1
                                    disp('ERROR: should use li == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
                                idx = (idx(1)-250):idx(end);
                                %idx = (idx(1)):idx(end);
                        end

                        time = t_step*(idx-idx(1));
                        idx_t{ri,ci,fi,di,ti,li} = idx;
                        time_t{ri,ci,fi,di,ti,li} = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,ci,fi,di,ti,li}.ox(1,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,ci,fi,di,ti,li}.f(1,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,ci,fi,di,ti,li}.Fp(1,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,ci,fi,di,ti,li}.v(1,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,ci,fi,di,ti,li}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,ci,fi,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,ci,fi,di,ti,li}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm force';

                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
            end
        end
    end
end
%xlim([0 0.7])
linkaxes(axh, 'xy');
%xlim([0 0.5]);
xlim([0 2])
sgtitle(titlestr);

%% %
subj_i = 2;
Results = cell(4,1);
for dir_i = 1:4
    Data = data(subj_i,dir_i,:,:,:,:);
    idx_t_tmp = idx_t(subj_i,dir_i,:,:,:,:);
    time_t_tmp = time_t(subj_i,dir_i,:,:,:,:);

    % shift the direction to negative when the direction is 3 or 4
    switch dir_i
        case [3,4]
        for fce_i = 1:3
            for dist_i = 1:3
                for trial_i = 1:15
                    trialtmp = Data{1,1,fce_i,dist_i,trial_i,1};
                    % as only these 2 used here, I only rotate these 2.
                    trialtmp.ox = -trialtmp.ox;
                    trialtmp.f = -trialtmp.f; 
                    Data{1,1,fce_i,dist_i,trial_i,1} = trialtmp;
                end
            end
        end
    end

    %Function Computing system Identification
    result = sys_id_pitt(Data,idx_t_tmp,time_t_tmp,1);
    Results{dir_i} = result;
end

%%% plot out on a 4-by-1 plot 

% stiffness
xx = [25 50 75];
ff = [15 20 25];
k_mat = ff'*((1./(xx/1000)));
fh = figure('name', 'stiffness across directions');
direction_array = {'right', 'front', 'left', 'back'};
for axh_i = 1:4
    axh(axh_i) = subplot(1,4,axh_i);
    
    set(gcf,'color','w');
    results = Results{axh_i};
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
    plot(x', k_mat', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5,...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor','none');

    grid on
    xlabel('Displacement [mm]')
    ylim([0 1100])
    ylabel('Stiffness [N/m]')
    title(direction_array{axh_i});
end
legend('15 N','20 N','25 N')

% damping
fh = figure('name', 'damping across directions');
direction_array = {'right', 'front', 'left', 'back'};
for axh_i = 1:4
    axh(axh_i) = subplot(1,4,axh_i);
    
    set(gcf,'color','w');
    results = Results{axh_i};
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
    ylim([0 50])
    ylabel('Damping [N/m*s]')
    title(direction_array{axh_i});

end
legend('15 N','20 N','25 N')

% mass
fh = figure('name', 'mass across directions');
direction_array = {'right', 'front', 'left', 'back'};
for axh_i = 1:4
    axh(axh_i) = subplot(1,4,axh_i);
    
    set(gcf,'color','w');
    results = Results{axh_i};
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
    ylim([0 6])
    ylabel('mass [kg]')
    title(direction_array{axh_i});

end
legend('15 N','20 N','25 N')

%% %
%%%% new way to interpret the K, B, and M
% use TF, but only use the "magnifier", "natural frequency", "damping coef"