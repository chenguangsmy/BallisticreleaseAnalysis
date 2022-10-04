% investigate the EMG values with the stiffness values. 
% Just plot them as the regression of either displacement (x), force (f) or
% stiffness (k)

% function fh = plotEMG_release_dir1_acrosssubject(obj)
%%%%%%%%%%%%%%%%%%% do something the same with the class function 
clear; clc; close all;
%
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
obj.data_name = 'ss4310_4356';

ifplot = 1; % the debug plot switch 
%         data
%         data_idx_ss
%         data_idx_tr
%         cond


            load([obj.data_dir '/' obj.data_name], 'data*', 'Results');
            obj.data = data;
            obj.data_idx_ss = data_index_ss;
            obj.data_idx_tr = data_index_tr;
            cond.subj = 1:size(data,1);
            cond.dir  = 1:size(data,2);
            cond.fce  = 1:size(data,3);
            cond.dist = 1:size(data,4);
            cond.trial= 1:size(data,5);
            cond.pert = 1:size(data,6);

            obj.cond = cond;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subtraction_mat = [ -1 1 -1 1 -1 1 -1 1;
                    -1 1 -1 1 -1 1 -1 1;
                    1 -1 1 -1 1 -1 1 -1;
                    1 -1 1 -1 1 -1 1 -1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 500;
fhtmp = figure();
col_type = colormap('lines');
close(fhtmp);

t_range = [-0.5 1]; % why???

emg_pair = [1 2 3 4 5 6 7 8];
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
    'elbow flexor', 'elbow extensor', ...
    'anterior deltoid', 'posterior deltiod', ...
    'pectoralis', 'Trapezius'};
cols = 3;
rows = 2 + length(emg_pair); % only plot position and the muscles
axh = zeros(rows, cols);
Dat_f_cell = cell(4,4);
Dat_x_cell = cell(4,4);
Dat_k_cell = cell(4,4);
Dat_b_cell = cell(4,4);
Dat_EMG_cell = cell(4,4);

fh = zeros(4,4,3);
for subj_i  = 1:4
    Dat_f = [];
    Dat_x = [];
    Dat_k = [];
    Dat_b = [];
    Dat_EMG = [];

    time_window_fhold = [-0.5 0]; % before release 0.5s take average
    time_window_rel = [-0.5 0.1]; % before hold and immediately after release

    % subj_i = 1;
    % fh(subj_i) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
    for dir_i   = 1:4
        pert_i  = 1;
        % across trials, get the EMG values as well as the x, f, k values
        % temperarily use k values as f/x, later substitute it with the estimation
        % values.

        % dat_avg % the data of averaged values; (for force and displacement, use
        % the absolute values)
        %

        % force, use actual force of 0.5s average;
        % displacement, use actual displacement from start point to 0.5s after
        % release
        % stiffness, use required force and required displacement.

        if(ifplot)
            fh(subj_i, dir_i, 1) = figure('name', ['pos subj' num2str(subj_i) 'dir' num2str(dir_i)]); % for position 
            fh(subj_i, dir_i, 2) = figure('name', ['fce subj' num2str(subj_i) 'dir' num2str(dir_i)]); % for displacement
            fh(subj_i, dir_i, 3) = figure('name', ['EMG subj' num2str(subj_i) 'dir' num2str(dir_i)]); % for EMG
            axh = cell(3,1); 
            axh{1} = zeros(3,3);
            axh{2} = zeros(3,3);
            axh{3} = zeros(3,3);
        end
        ifplotEMG_scatter = 1;
        if(ifplotEMG_scatter)
            markers = '*so';
            for chi = 1:8
                figure(chi); hold on;
                title(['scatter of each channel EMG - F/x/k Ch' num2str(chi)]);
            end
        end
        for fce_i = 1:3
            for dist_i = 1:3
                
                if (ifplot)
                    set(0,'currentfigure',fh(subj_i,dir_i,1)); 
                    axh{1}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                    set(0,'currentfigure',fh(subj_i,dir_i,2)); 
                    axh{2}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                    set(0,'currentfigure',fh(subj_i,dir_i,3)); 
                    axh{3}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                end
                trials_list = obj.cond.trial;
                % make enough space for data
                dat.row = length(trials_list);
                dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
                    obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
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

                dat_idx_fhold = dat.t>time_window_fhold(1) & dat.t<time_window_fhold(2); % the time at force hold
                [~, dat_idx_phold] = min(abs(dat.t - 0.5));  % the time at position hold
                                                            % Assuming at this time, the handle was stopped
                dat_avg.x = abs(dat.pos(:,dat_idx_phold) - mean(dat.pos(:,dat_idx_fhold),2));
                dat_avg.f = abs(mean(dat.fce(:,dat_idx_fhold), 2));
                dat_avg.emg_fhd = mean(dat.emg(:,:,dat_idx_fhold), 3);

                %         hold on;
                %         plot(dat_avg.f, dat_avg.emg_fhd(1,:), '.');
                %         plot(dat_avg.x, dat_avg.emg_fhd(1,:), '.');
                %         plot(dat_avg.f./dat_avg.x, dat_avg.emg_fhd(1,:), '.');

                Dat_f = [Dat_f; dat_avg.f];
                Dat_x = [Dat_x; dat_avg.x];
%                 Dat_k = [Dat_k; dat_avg.f./dat_avg.x];
                Dat_k = [Dat_k reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_b = [Dat_b reshape(Results(subj_i,dir_i).B_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_EMG = [Dat_EMG dat_avg.emg_fhd];
                Dat_EMG_tmp = dat_avg.emg_fhd;
                Dat_k_tmp = reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9]);
                Dat_f_tmp = dat_avg.f;
                Dat_x_tmp = dat_avg.x;
                
                if(ifplotEMG_scatter)
                    for ch_i = 1:8
                        figure(ch_i);
                        plot(Dat_k_tmp, Dat_EMG_tmp(ch_i,:), markers(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
                        xlabel('K (m/s)'); ylabel('relative EMG'); 
%                         plot(Dat_x_tmp, Dat_EMG_tmp(ch_i,:), markers(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
%                         plot(Dat_f_tmp, Dat_EMG_tmp(ch_i,:), markers(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
                    end
                end
                
                if(ifplot)
                    % for position plot, plot original data line, the phold
                    % and fhold 
                    
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t, dat.pos);
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t(dat_idx_phold), dat.pos(:,dat_idx_phold), 'g.'); % position hold 
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t(dat_idx_fhold), dat.pos(:,dat_idx_fhold), 'g.'); % force hold
                    yline(subplot(axh{1}(fce_i,dist_i)),mean(dat.pos(:,dat_idx_fhold),2));

                    % for force plot, plot original data line, the average 
                    % at fhold
                    
                    plot(subplot(axh{2}(fce_i,dist_i)),dat.t, dat.fce);
                    plot(subplot(axh{2}(fce_i,dist_i)),dat.t(dat_idx_fhold), dat.fce(:,dat_idx_fhold), 'g.'); % force hold
                    yline(subplot(axh{2}(fce_i,dist_i)),mean(dat.fce(:,dat_idx_fhold),2));

                    % for emg plot, plot original data line, and the
                    % average at fhold
                    %% option 1, plot the all channels with all trials 
%                     for ch_i = 1:8
%                         plot(subplot(axh{3}(fce_i,dist_i)),dat.t, reshape(dat.emg(ch_i,:,:),size(dat.emg,[2,3])), 'color', col_type(ch_i,:));
%                         plot(subplot(axh{3}(fce_i,dist_i)),dat.t(dat_idx_fhold), reshape(dat.emg(ch_i,:,dat_idx_fhold),[size(dat.emg,2), sum(dat_idx_fhold)]), '.', 'MarkerEdgeColor', col_type(ch_i,:)); % force hold
%                         yline(subplot(axh{3}(fce_i,dist_i)),reshape(mean(dat.emg(ch_i,:,dat_idx_fhold),3),[1,size(dat.emg,2)]), 'color', col_type(ch_i,:));
%                     end
                    %% option 2, plot 8th channel, different trial different line
                    legend_arr = cell(9,1); 
                    for trial_i = 1:9
                        plot(subplot(axh{3}(fce_i,dist_i)),dat.t, reshape(dat.emg(8,trial_i,:),1, size(dat.emg,3)), 'linewidth', 2);
                        legend_arr{trial_i} = ['ss' num2str(data_index_ss(subj_i,dir_i,fce_i,dist_i,trial_i)) 'tr' num2str(data_index_tr(subj_i,dir_i,fce_i,dist_i,trial_i))];
                    end
                    legend(legend_arr);
                end
            end
        end
        Dat_f_cell{subj_i,dir_i} = Dat_f;
        Dat_x_cell{subj_i,dir_i} = Dat_x;
        Dat_k_cell{subj_i,dir_i} = Dat_k;
        Dat_b_cell{subj_i,dir_i} = Dat_b;
        Dat_EMG_cell{subj_i,dir_i} = Dat_EMG;
    end
end

%% The figure of EMG tuning
% In this section, do a figure of the EMG tuning through F, x and K. 
% a figure of 8 (channel-num) -by - 3 (x,F,K) matrix. 
% Each plot x is the Force/distance value; and y is the EMG activities;
% In each plot, the 3-force and 3-distance levels using different
% marker-color conbinations;
% In each plot, do a regression and show the R^2 indicating the goodness of
% fitting
% Compare the F and K plot, I can know whether the before-release EMG tune
% to stiffness levels 
%%%%%%%%%%%%%%%%%%% do something the same with the class function 
clear; clc; close all;
%
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
obj.data_name = 'ss4310_4356';
% obj.data_name = 'ss4310_4314';

ifplot = 0; % the debug plot switch 
%         data
%         data_idx_ss
%         data_idx_tr
%         cond


load([obj.data_dir '/' obj.data_name], 'data*', 'Results');
obj.data = data;
obj.data_idx_ss = data_index_ss;
obj.data_idx_tr = data_index_tr;
cond.subj = 1:size(data,1);
cond.dir  = 1:size(data,2);
cond.fce  = 1:size(data,3);
cond.dist = 1:size(data,4);
cond.trial= 1:size(data,5);
cond.pert = 1:size(data,6);

obj.cond = cond;
t_range = [-0.5 1]; % why???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subtraction_mat = [ -1 1 -1 1 -1 1 -1 1;
                    -1 1 -1 1 -1 1 -1 1;
                    1 -1 1 -1 1 -1 1 -1;
                    1 -1 1 -1 1 -1 1 -1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitting: a structure of multiple slope, intercept, R2... 
Dat_f = nan(4,4,3,3,9);
Dat_x = nan(4,4,3,3,9);
Dat_k = nan(4,4,3,3,9);
Dat_b = nan(4,4,3,3,9);
Dat_EMG_fhd = nan(4,4,3,3,9,8);
fh = zeros(4,4,3);
mks_type = '*so';
Fs = 500;
fhtmp = figure();
col_type = colormap('lines');
close(fhtmp);

time_window_fhold = [-0.5 0]; % before release 0.5s take average
time_window_atfly = [0 0.1];  % after release but before 0.1s 
time_window_phold = [0.5 0.8];% the EMG activities at position hold

for subj_i  = 1:4    
    for dir_i   = 1:4
        pert_i  = 1;
        % across trials, get the EMG values as well as the x, f, k values
        % temperarily use k values as f/x, later substitute it with the estimation
        % values.

        % dat_avg % the data of averaged values; (for force and displacement, use
        % the absolute values)
        %
        % force, use actual force of 0.5s average;
        % displacement, use actual displacement from start point to 0.5s after
        % release
        % stiffness, use required force and required displacement.

        ifplotEMG_scatter = 0;

        Dat_f_pool = [];
        Dat_x_pool = [];
        Dat_k_pool = [];
        Dat_b_pool = [];
        Dat_EMG_pool = [];

        if(ifplotEMG_scatter) % plot it latter
           
            for chi = 1:8
                figure(chi); hold on;
                title(['scatter of each channel EMG - F/x/k Ch' num2str(chi)]);
            end
        end

        if(ifplot)
            fh(subj_i,dir_i,1) = figure();
            fh(subj_i,dir_i,2) = figure();
            fh(subj_i,dir_i,3) = figure();
        end

        for fce_i = 1:3
            for dist_i = 1:3
                
                if (ifplot) % debug plot 
                    set(0,'currentfigure',fh(subj_i,dir_i,1)); 
                    axh{1}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                    set(0,'currentfigure',fh(subj_i,dir_i,2)); 
                    axh{2}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                    set(0,'currentfigure',fh(subj_i,dir_i,3)); 
                    axh{3}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                end
                trials_list = obj.cond.trial;
                % make enough space for data
                dat.row = length(trials_list);
                dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
                    obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
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

                dat_idx_fhold = dat.t>time_window_fhold(1) & dat.t<time_window_fhold(2); % the time at force hold
                [~, dat_idx_phold] = min(abs(dat.t - 0.5));  % the time at position hold
                                                            % Assuming at this time, the handle was stopped
                dat_reg.x = abs(dat.pos(:,dat_idx_phold) - mean(dat.pos(:,dat_idx_fhold),2));
                dat_reg.f = abs(mean(dat.fce(:,dat_idx_fhold), 2));
                dat_reg.emg_fhd = mean(dat.emg(:,:,dat_idx_fhold), 3);

                %         hold on;
                %         plot(dat_avg.f, dat_avg.emg_fhd(1,:), '.');
                %         plot(dat_avg.x, dat_avg.emg_fhd(1,:), '.');
                %         plot(dat_avg.f./dat_avg.x, dat_avg.emg_fhd(1,:), '.');

                Dat_f_pool = [Dat_f_pool; dat_reg.f];
                Dat_x_pool = [Dat_x_pool; dat_reg.x];
                Dat_k_pool = [Dat_k_pool reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_b_pool = [Dat_b_pool reshape(Results(subj_i,dir_i).B_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_EMG_pool = [Dat_EMG_pool dat_reg.emg_fhd];

                Dat_f(subj_i,dir_i,fce_i,dist_i,:) = dat_reg.f;
                Dat_x(subj_i,dir_i,fce_i,dist_i,:) = dat_reg.x;
                Dat_k(subj_i,dir_i,fce_i,dist_i,:) = reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9]);
                Dat_b(subj_i,dir_i,fce_i,dist_i,:) = reshape(Results(subj_i,dir_i).B_up_tr(fce_i, dist_i, :), [1,9]);
                fce_i;
                dist_i;
                dat_reg.emg_fhd;
                Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,:) = dat_reg.emg_fhd';

                Dat_EMG_tmp = Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,:);
                Dat_k_tmp = reshape(Dat_k(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                Dat_f_tmp = reshape(Dat_f(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                Dat_x_tmp = reshape(Dat_x(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                
                if(ifplotEMG_scatter)
                    for ch_i = 1:8
                        figure(ch_i);
                        plt_x_tmp = Dat_k_tmp;
                        plt_y_tmp = reshape(Dat_EMG_tmp(subj_i,dir_i,fce_i,dist_i,:,ch_i),1,size(Dat_EMG_tmp,5));
                        plot(plt_x_tmp, plt_y_tmp, mks_type(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
                        xlabel('K (m/s)'); ylabel('relative EMG'); 
%                         plot(Dat_x_tmp, Dat_EMG_tmp(ch_i,:), mks_type(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
%                         plot(Dat_f_tmp, Dat_EMG_tmp(ch_i,:), mks_type(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
                    end
                end
                
                if(ifplot)
                    % for position plot, plot original data line, the phold
                    % and fhold 
                    
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t, dat.pos);
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t(dat_idx_phold), dat.pos(:,dat_idx_phold), 'g.'); % position hold 
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t(dat_idx_fhold), dat.pos(:,dat_idx_fhold), 'g.'); % force hold
                    yline(subplot(axh{1}(fce_i,dist_i)),mean(dat.pos(:,dat_idx_fhold),2));

                    % for force plot, plot original data line, the average 
                    % at fhold
                    
                    plot(subplot(axh{2}(fce_i,dist_i)),dat.t, dat.fce);
                    plot(subplot(axh{2}(fce_i,dist_i)),dat.t(dat_idx_fhold), dat.fce(:,dat_idx_fhold), 'g.'); % force hold
                    yline(subplot(axh{2}(fce_i,dist_i)),mean(dat.fce(:,dat_idx_fhold),2));

                    % for emg plot, plot original data line, and the
                    % average at fhold
                    % option 1, plot the all channels with all trials 
                    for ch_i = 1:8
                        plot(subplot(axh{3}(fce_i,dist_i)),dat.t, reshape(dat.emg(ch_i,:,:),size(dat.emg,[2,3])), 'color', col_type(ch_i,:));
                        plot(subplot(axh{3}(fce_i,dist_i)),dat.t(dat_idx_fhold), reshape(dat.emg(ch_i,:,dat_idx_fhold),[size(dat.emg,2), sum(dat_idx_fhold)]), '.', 'MarkerEdgeColor', col_type(ch_i,:)); % force hold
                        yline(subplot(axh{3}(fce_i,dist_i)),reshape(mean(dat.emg(ch_i,:,dat_idx_fhold),3),[1,size(dat.emg,2)]), 'color', col_type(ch_i,:));
                    end
                    % option 2, plot 8th channel, different trial different line
%                     legend_arr = cell(9,1); 
%                     for trial_i = 1:9
%                         plot(subplot(axh{3}(fce_i,dist_i)),dat.t, reshape(dat.emg(8,trial_i,:),1, size(dat.emg,3)), 'linewidth', 2);
%                         legend_arr{trial_i} = ['ss' num2str(data_index_ss(subj_i,dir_i,fce_i,dist_i,trial_i)) 'tr' num2str(data_index_tr(subj_i,dir_i,fce_i,dist_i,trial_i))];
%                     end
%                     legend(legend_arr);
                end
            end
        end
        Dat_f_cell{subj_i,dir_i} = Dat_f_pool;
        Dat_x_cell{subj_i,dir_i} = Dat_x_pool;
        Dat_k_cell{subj_i,dir_i} = Dat_k_pool;
        Dat_b_cell{subj_i,dir_i} = Dat_b_pool;
        Dat_EMG_cell{subj_i,dir_i} = Dat_EMG_pool;
    end
end
%% % plot out the figure
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
    'elbow flexor', 'elbow extensor', ...
    'anterior deltoid', 'posterior deltiod', ...
    'pectoralis', 'Trapezius'};
clear axh
for subj_i = 1:4
    for dir_i = 1:4
        fh = figure('name', ['EMG tuning subj' num2str(subj_i) 'dirc' num2str(dir_i)], ...
            'position', [0 0 600, 1440]);
        for ch_i = 1:8 % EMG channels
            for indpv = 1:3 % 1. Force; 2. x; 3. K
                axh(ch_i, indpv) = subplot(8, 3, (ch_i-1)*3+indpv);
                % for each force and target distance, plot the regressed
                % EMG activity with each independent variable
                hold on;
                for fce_i = 1:3
                    for dist_i = 1:3
                        y_pts = reshape(Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i),1,9);
                        switch indpv
                            case 1 % Force
                                x_pts = reshape(Dat_f(subj_i,dir_i,fce_i,dist_i,:),1,9);
                                xlabel_str = 'Force (N)';
                            case 2 % x
                                x_pts = reshape(Dat_x(subj_i,dir_i,fce_i,dist_i,:),1,9);
                                xlabel_str = 'displacement (m)';
                            case 3 % K
                                x_pts = reshape(Dat_k(subj_i,dir_i,fce_i,dist_i,:),1,9);
                                xlabel_str = 'Stiffness (N/m)';
                        end
                        plot(axh(ch_i,indpv), x_pts, y_pts, ...
                            'Marker', mks_type(fce_i), ...
                            'Linestyle', 'none', ...
                            'Color', col_type(dist_i+4,:));
                        % display to see the outlairs
                        y_pts;
                        reshape([data_index_ss(subj_i,dir_i,fce_i,dist_i,:,1)],1,9);
                        reshape([data_index_tr(subj_i,dir_i,fce_i,dist_i,:,1)],1,9);
                    end
                end

                switch indpv 
                    case 1
%                         x = Dat_f_cell{subj_i,dir_i};
                        x = reshape(Dat_f(subj_i,dir_i,:,:,:),1,81);
                        ylabel_str = {emg_pair_label{ch_i}; 'normalized EMG'}; 
                    case 2
%                         x = Dat_x_cell{subj_i,dir_i};
                        x = reshape(Dat_x(subj_i,dir_i,:,:,:),1,81);
                        ylabel_str = 'normalized EMG'; 

                    case 3
%                         x = Dat_k_cell{subj_i,dir_i};
                        x = reshape(Dat_k(subj_i,dir_i,:,:,:),1,81);
                        ylabel_str = 'normalized EMG'; 
                end

                
                xlabel(xlabel_str);
                ylabel(ylabel_str);

%                 y  = Dat_EMG_cell{subj_i,dir_i}(ch_i,:);  
                y = reshape(Dat_EMG_fhd(subj_i,dir_i,:,:,:,ch_i),1,81);

                y_validx = ~isnan(y); 
                % do regression part
                x = x(y_validx);
                y = y(y_validx); 
%                 plot(x, y, '.', 'color', [0.2 0.2 0.2]); % sanity check making sure the same dtpts

                [P, S] = polyfit(x,y,1); 
                R2 = 1-(S.normr/norm(y-mean(y)))^2;
                R2_qualified = R2>0.2;
                [b, bint, r, rint, stats] = regress(y', [ones(size(x)); x]' );
                p_val = stats(3);
                R2_qualified = p_val<128; %(bufforoni correction)
                yfit1 = polyval(P,x); 

                plot(x, yfit1, 'linewidth', 2, 'color', [0.7 0.7 0.7]); 
                if (R2_qualified)
                    title(['R^2:' num2str(R2)], 'color', [0 0 0]);
                else
                    title(['R^2:' num2str(R2)], 'color', [0.5 0.5 0.5]);
                end

            end
        end
        sgtitle(['subj' num2str(subj_i) 'direc' num2str(dir_i) 'EMG fhold']);
        if (subj_i == 1 && dir_i == 1)
        exportgraphics(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_regress.pdf');
        else
        exportgraphics(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_regress.pdf', 'Append',true);
        end
    end
end
%% Use the ANOVA to see if the EMG before release is correlated with any indp vars
Dat_EMG_fhd; 
Dat_f;
Dat_x;
Dat_k;
Dat_subj = zeros(4,4,3,3,9); 
Dat_dir = zeros(4,4,3,3,9);
for subj_i = 1:4 
    for dir_i = 1:4
        Dat_subj(subj_i,dir_i,:,:,:) = subj_i;
        Dat_dir(subj_i,dir_i,:,:,:) = dir_i;
    end
end

% subject 
% direction 
%EMG ~ 1+ Subj + Dir + Force + displacement %... option1
%EMG ~ 1+ Subj + Dir + Force + k            %... option2
blr.f = Dat_f(:); 
blr.x = Dat_x(:);
blr.k = Dat_k(:);
blr.subj = Dat_subj(:);
blr.dir = Dat_dir(:);

Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,1);
blr.emg1 = Dat_EMG_tmp(:); % fhd
Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,2);
blr.emg2 = Dat_EMG_tmp(:); % fhd
Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,3);
blr.emg3 = Dat_EMG_tmp(:); % fhd
Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,4);
blr.emg4 = Dat_EMG_tmp(:); % fhd
Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,5);
blr.emg5 = Dat_EMG_tmp(:); % fhd
Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,6);
blr.emg6 = Dat_EMG_tmp(:); % fhd
Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,7);
blr.emg7 = Dat_EMG_tmp(:); % fhd
Dat_EMG_tmp = Dat_EMG_fhd(:,:,:,:,:,8);
blr.emg8 = Dat_EMG_tmp(:); % fhd

tbl = table(blr.subj, blr.dir, blr.f, blr.x, blr.k, ...
    blr.emg1, blr.emg2, blr.emg3, blr.emg4, blr.emg5, blr.emg6, blr.emg7, blr.emg8, ...
    'VariableNames', {'subj', 'dir', 'fce', 'disp', 'stiffness', ...
    'emg1', 'emg2', 'emg3', 'emg4', 'emg5', 'emg6', 'emg7', 'emg8'});

mdl1 = fitlm(tbl, 'emg1 ~ subj + dir + fce + stiffness')

%% plot F, x, EMG on each of the muscle
% if (1)
%     for subj_i = 1:4
%         Dat_f = Dat_f_cell{subj_i,dir_i};
%         Dat_x = Dat_x_cell{subj_i,dir_i};
%         Dat_k = Dat_k_cell{subj_i,dir_i};
%         Dat_b = Dat_b_cell{subj_i,dir_i};
%         Dat_EMG = Dat_EMG_cell{subj_i,dir_i};
%         figure('name', 'single EMG with F, x and k', 'position', [0 0 1200 300]);
% 
% %         subplot(4,1,1); hold on;
%         for ch_i = 1:8
%             figure(); hold on;
%             dat_EMG_tmp = Dat_EMG(ch_i,:);
%             dat_EMG_tmp_mean = mean(Dat_EMG(ch_i,:), 'omitnan');
%             dat_EMG_tmp_std = std(Dat_EMG(ch_i,:), 'omitnan');
%             dat_EMG_tmp_validx = dat_EMG_tmp > (dat_EMG_tmp_mean-3*dat_EMG_tmp_std) & ...
%                 dat_EMG_tmp < (dat_EMG_tmp_mean+3*dat_EMG_tmp_std);
% %             dat_f_tmp = Dat_f;
%             dat_f_tmp = Dat_k;
%             plot(dat_EMG_tmp, dat_f_tmp, '.', 'Color', col_type(ch_i,:));
% %             refline;
%             x = dat_EMG_tmp(dat_EMG_tmp_validx);
%             y = dat_f_tmp(dat_EMG_tmp_validx);
%             plot(x, y, 'o', 'Color', col_type(ch_i,:));
%             refline;
% 
% %             plot(Dat_EMG(ch_i,:),Dat_f, '.', 'Color', col_type(ch_i,:));  
% 
% %             p = polyfit(Dat_EMG(ch_i,:),Dat_f, 1);
% %             tmp_x = 10:0.1:30;
% %             tmp_y = p(1) * tmp_x + p(2); 
% %             lnh(ch_i) = plot(tmp_x, tmp_y, 'Color', col_type(ch_i,:), 'LineWidth', 1.5);
%         end
%         lnh = refline;
%         legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
%         xlabel('Force exertion (N)');
%         ylabel('EMG values perportion to MVF EMG');
%         title(['EMG corresponding to F, subj' num2str(subj_i)]);
% 
% 
%         %%%
%         subplot(4,1,2); hold on;
%         for ch_i = 1:8
%             plot(Dat_EMG(ch_i,:), Dat_x, '.', 'Color', col_type(ch_i,:));
%         end
%         lnh = refline;
%         legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
%         xlabel('actual displacement (m)');
%         ylabel('EMG values perportion to MVF EMG');
%         title(['EMG corresponding to x, subj' num2str(subj_i)]);
% 
%         %%%
%         subplot(4,1,3); hold on;
%         for ch_i = 1:8
%             plot(Dat_EMG(ch_i,:), Dat_k, '.', 'Color', col_type(ch_i,:));
%         end
%         lnh = refline;
%         legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
%         xlabel('K_{est} (N/m)');
%         ylabel('EMG values perportion to MVF EMG');
%         title(['EMG corresponding to K, subj' num2str(subj_i)]);
% 
%         subplot(4,1,4); hold on;
%         for ch_i = 1:8
%             plot(Dat_EMG(ch_i,:), Dat_b, '.', 'Color', col_type(ch_i,:));
%         end
%         lnh = refline;
%         legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
%         xlabel('B_{est} (N/m)');
%         ylabel('EMG values perportion to MVF EMG');
%         title(['EMG corresponding to B, subj' num2str(subj_i)]);
%     end
% end
% Think: how to represent the strength of tuning in the muscles??? 

%% plot EMG summation and subtraction to see 'stiffness' and 'force' values
if(1)
    for subj_i = 1:4
        Dat_f = Dat_f_cell{subj_i,dir_i};
        Dat_x = Dat_x_cell{subj_i,dir_i};
        Dat_k = Dat_k_cell{subj_i,dir_i};
        Dat_b = Dat_b_cell{subj_i,dir_i};
        Dat_EMG_fhd = Dat_EMG_cell{subj_i,dir_i};
        figure('name', ['EMG summations dir' num2str(dir_i)], 'position', [0 0 1200 300]);

        subplot(1,3,1); hold on;
        for emg_pr_i = 1:4
            %     EMG_m = Dat_EMG((emg_pr_i-1)*2+1,:) - Dat_EMG((emg_pr_i-1)*2+2,:);
            EMG_p = Dat_EMG_fhd((emg_pr_i-1)*2+1,:) + Dat_EMG_fhd((emg_pr_i-1)*2+2,:);
            plot(Dat_f, EMG_p, '.');
        end
        lnh = refline;
        legend(lnh, {'pair1','pair2','pair3','pair4'});
        xlabel('Force exertion (N)');
        ylabel('EMG pair summation');
        title(['EMG corresponding to F, subj' num2str(subj_i)]);


        %%%
        subplot(1,3,2); hold on;
        for emg_pr_i = 1:4
            EMG_p = Dat_EMG_fhd((emg_pr_i-1)*2+1,:) + Dat_EMG_fhd((emg_pr_i-1)*2+2,:);
            plot(Dat_x, EMG_p, '.');
        end
        lnh = refline;
        legend(lnh, {'pair1','pair2','pair3','pair4'});
        xlabel('actual displacement (m)');
        ylabel('EMG pair summation');
        title(['EMG corresponding to x, subj' num2str(subj_i)]);

        %%%
        subplot(1,3,3); hold on;
        for emg_pr_i = 1:4
            EMG_p = Dat_EMG_fhd((emg_pr_i-1)*2+1,:) + Dat_EMG_fhd((emg_pr_i-1)*2+2,:);
            plot(Dat_k, EMG_p, '.');
        end
        lnh = refline;
        legend(lnh, {'pair1','pair2','pair3','pair4'});
        xlabel('dF/dx (N/m)');
        ylabel('EMG pair summation');
        title(['EMG corresponding to K, subj' num2str(subj_i)]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for subj_i = 1:4
    if(1)
        Dat_f = Dat_f_cell{subj_i,dir_i};
        Dat_x = Dat_x_cell{subj_i,dir_i};
        Dat_k = Dat_k_cell{subj_i,dir_i};
        Dat_EMG_fhd = Dat_EMG_cell{subj_i,dir_i};
        figure('name', ['EMG subtractions dir' num2str(dir_i)], 'position', [0 0 1200 300]);

        subplot(1,3,1); hold on;
        for emg_pr_i = 1:4
            idx1 = (emg_pr_i-1)*2+1;
            idx2 = (emg_pr_i-1)*2+2;
            EMG_m = Dat_EMG_fhd(idx1,:)*subtraction_mat(dir_i,idx1) + Dat_EMG_fhd(idx2,:)*subtraction_mat(dir_i,idx2);
            plot(Dat_f, EMG_m, '.');
        end
        lnh = refline;
        legend(lnh, {'pair1','pair2','pair3','pair4'});
        xlabel('Force exertion (N)');
        ylabel('EMG pair subtraction');
        title(['EMG corresponding to F, subj' num2str(subj_i)]);


        %%%
        subplot(1,3,2); hold on;
        for emg_pr_i = 1:4
            switch emg_pr_i
                case {1 2}
                    EMG_m = Dat_EMG_fhd((emg_pr_i-1)*2+1,:) - Dat_EMG_fhd((emg_pr_i-1)*2+2,:);
                case {3 4}
                    EMG_m = Dat_EMG_fhd((emg_pr_i-1)*2+2,:) - Dat_EMG_fhd((emg_pr_i-1)*2+1,:);
            end
            plot(Dat_x, EMG_m, '.');
        end
        lnh = refline;
        legend(lnh, {'pair1','pair2','pair3','pair4'});
        xlabel('actual displacement (m)');
        ylabel('EMG pair subtraction');
        title(['EMG corresponding to x, subj' num2str(subj_i)]);

        %%%
        subplot(1,3,3); hold on;
        for emg_pr_i = 1:4
            switch emg_pr_i
                case {1 2}
                    EMG_m = Dat_EMG_fhd((emg_pr_i-1)*2+1,:) - Dat_EMG_fhd((emg_pr_i-1)*2+2,:);
                case {3 4}
                    EMG_m = Dat_EMG_fhd((emg_pr_i-1)*2+2,:) - Dat_EMG_fhd((emg_pr_i-1)*2+1,:);
            end
            plot(Dat_k, EMG_m, '.');
        end
        lnh = refline;
        legend(lnh, {'pair1','pair2','pair3','pair4'});
        xlabel('dF/dx (N/m)');
        ylabel('EMG pair subtraction');
        title(['EMG corresponding to K, subj' num2str(subj_i)]);


    end
end

%% Plot the EMG tuning at release time with the F, x and k

clear; clc; close all;
%
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
obj.data_name = 'ss4310_4356';

ifplot = 0; % the debug plot switch 
%         data
%         data_idx_ss
%         data_idx_tr
%         cond


load([obj.data_dir '/' obj.data_name], 'data*', 'Results');
obj.data = data;
obj.data_idx_ss = data_index_ss;
obj.data_idx_tr = data_index_tr;
cond.subj = 1:size(data,1);
cond.dir  = 1:size(data,2);
cond.fce  = 1:size(data,3);
cond.dist = 1:size(data,4);
cond.trial= 1:size(data,5);
cond.pert = 1:size(data,6);

obj.cond = cond;
t_range = [-0.5 1]; % why???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subtraction_mat = [ -1 1 -1 1 -1 1 -1 1;
                    -1 1 -1 1 -1 1 -1 1;
                    1 -1 1 -1 1 -1 1 -1;
                    1 -1 1 -1 1 -1 1 -1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitting: a structure of multiple slope, intercept, R2... 
Dat_f = nan(4,4,3,3,9);
Dat_x = nan(4,4,3,3,9);
Dat_k = nan(4,4,3,3,9);
Dat_b = nan(4,4,3,3,9);
Dat_EMG_fhd = nan(4,4,3,3,9,8);
Dat_EMG_afl = nan(4,4,3,3,9,8);
Dat_EMG_phd = nan(4,4,3,3,9,8);

fh = zeros(4,4,3); % figure handles
mks_type = '*so';
Fs = 500;
fhtmp = figure();
col_type = colormap('lines');
close(fhtmp);

time_window_fhold = [-0.5 0]; % before release 0.5s take average
time_window_atfly = [0.03 0.2];  % after release but before 0.1s 
time_window_phold = [0.5 0.8];% the EMG activities at position hold

EMG_POL_DEF_SINGLE = reshape([ 1 -1 1 -1 1 -1 1 -1; ...
                               1 -1 1 -1 1 -1 1 -1; ...
                               -1 1 -1 1 -1 1 -1 1; ...
                               -1 1 -1 1 -1 1 -1 1 ], [1,4,8]);
% EMG_POL_DEF = repmat(EMG_POL_DEF_SINGLE, [4,1,1]); % subj_i, dir_i, ch_i
EMG_POL_DEF(1,:,:) = [1 -1 1 -1 1 -1 1 -1;
                      1 -1 1 -1 -1 -1 1 1;
                      -1 1 -1 1 -1 1 -1 1;
                      -1 1 1 1 1 1 -1 -1];
EMG_POL_DEF(2,:,:) = [-1 -1 1 -1 1 -1 1 1;
                      -1 1 1 -1 1 1 -1 -1;
                       1 1 -1 1 1 1 -1 1;
                      -1 1 1 1 1 1 -1 1];
EMG_POL_DEF(3,:,:) = [1 1 1 -1 -1 -1 1 -1;
                      1 -1 1 -1 1 -1 1 1;
                     -1 1 -1 1 1 1 1 1;
                     -1 1 -1 1 1 1 1 1];
EMG_POL_DEF(4,:,:) = [-1 1 1 -1 1 -1 1 -1;
                      1 -1 1 -1 -1 1 -1 1;
                      1 1 -1 1 -1 1 -1 1;
                      1 1 -1 1 1 1 1 -1];
% problem matrix!!! Think it latter on the sequence!!! 


for subj_i  = 1:4    
    for dir_i   = 1:4
        pert_i  = 1;
        % across trials, get the EMG values as well as the x, f, k values
        % temperarily use k values as f/x, later substitute it with the estimation
        % values.

        % dat_avg % the data of averaged values; (for force and displacement, use
        % the absolute values)
        %
        % force, use actual force of 0.5s average;
        % displacement, use actual displacement from start point to 0.5s after
        % release
        % stiffness, use required force and required displacement.

        ifplotEMG_scatter = 0;

        if(ifplotEMG_scatter) % plot it latter
           
            for chi = 1:8
                figure(chi); hold on;
                title(['scatter of each channel EMG - F/x/k Ch' num2str(chi)]);
            end
        end

        if(ifplot)
            fh(subj_i,dir_i,1) = figure();
            fh(subj_i,dir_i,2) = figure();
            for ch_i = 1:8
            fh(subj_i,dir_i,2+ch_i) = figure();
            end
        end

        for fce_i = 1:3
            for dist_i = 1:3
                
                if (ifplot) % debug plot 
                    set(0,'currentfigure',fh(subj_i,dir_i,1)); 
                    axh{1}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                    set(0,'currentfigure',fh(subj_i,dir_i,2)); 
                    axh{2}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                    hold on;
                    for ch_i = 1:8
                        set(0,'currentfigure',fh(subj_i,dir_i,2+ch_i)); 
                        axh{2+ch_i}(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i);
                        hold on;
                    end
                end
                trials_list = obj.cond.trial;
                % make enough space for data
                dat.row = length(trials_list);
                dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
                    obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
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

                dat_idx_fhold = dat.t>time_window_fhold(1) & dat.t<time_window_fhold(2); % the time at force hold
                dat_idx_atfly = dat.t>time_window_atfly(1) & dat.t<time_window_atfly(2);
                dat_idx_phold = dat.t>time_window_phold(1) & dat.t<time_window_phold(2);

                % variables in figure plot
                dat_reg.x = abs(mean(dat.pos(:,dat_idx_phold),2) - mean(dat.pos(:,dat_idx_fhold),2));
                dat_reg.f = abs(mean(dat.fce(:,dat_idx_fhold), 2));
                dat_reg.emg_fhd = mean(dat.emg(:,:,dat_idx_fhold), 3);
                dat_reg.emg_phd = mean(dat.emg(:,:,dat_idx_phold), 3);
                    % whether get max or min is defined on the top 
                dat_reg.emg_afl = zeros(size(dat_reg.emg_fhd)); % 8-by-9
                dat_reg.emg_aflidx = zeros(size(dat_reg.emg_fhd)); % 8-by-9

                for ch_i = 1:8
                    pol = EMG_POL_DEF(subj_i,dir_i,ch_i);
                    for tr_i = 1:9
                        switch pol
                            case -1
                                [val_tmp, idx_tmp] = min(dat.emg(ch_i,tr_i,dat_idx_atfly),[],3);
                            case 1
                                [val_tmp, idx_tmp] = max(dat.emg(ch_i,tr_i,dat_idx_atfly),[],3);
                        end
                        dat_idx_atfly_num = find(dat_idx_atfly);
                        dat_reg.emg_afl(ch_i,tr_i) = val_tmp;
                        dat_reg.emg_aflidx(ch_i,tr_i) = dat_idx_atfly_num(idx_tmp);
                    end
                end

                %         hold on;
                %         plot(dat_avg.f, dat_avg.emg_fhd(1,:), '.');
                %         plot(dat_avg.x, dat_avg.emg_fhd(1,:), '.');
                %         plot(dat_avg.f./dat_avg.x, dat_avg.emg_fhd(1,:), '.');


                Dat_f(subj_i,dir_i,fce_i,dist_i,:) = dat_reg.f;
                Dat_x(subj_i,dir_i,fce_i,dist_i,:) = dat_reg.x;
                Dat_k(subj_i,dir_i,fce_i,dist_i,:) = reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9]);
                Dat_b(subj_i,dir_i,fce_i,dist_i,:) = reshape(Results(subj_i,dir_i).B_up_tr(fce_i, dist_i, :), [1,9]);
                fce_i;
                dist_i;
                dat_reg.emg_fhd;
                Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,:) = dat_reg.emg_fhd';
                Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,:) = dat_reg.emg_afl';
                Dat_EMG_aflidx(subj_i,dir_i,fce_i,dist_i,:,:) = dat_reg.emg_aflidx';
                Dat_EMG_phd(subj_i,dir_i,fce_i,dist_i,:,:) = dat_reg.emg_phd';

                Dat_EMG_tmp = Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,:);
                Dat_k_tmp = reshape(Dat_k(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                Dat_f_tmp = reshape(Dat_f(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                Dat_x_tmp = reshape(Dat_x(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                
                if(ifplotEMG_scatter)
                    for ch_i = 1:8
                        figure(ch_i);
                        plt_x_tmp = Dat_k_tmp;
                        plt_y_tmp = reshape(Dat_EMG_tmp(subj_i,dir_i,fce_i,dist_i,:,ch_i),1,size(Dat_EMG_tmp,5));
                        plot(plt_x_tmp, plt_y_tmp, mks_type(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
                        xlabel('K (m/s)'); ylabel('relative EMG'); 
%                         plot(Dat_x_tmp, Dat_EMG_tmp(ch_i,:), mks_type(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
%                         plot(Dat_f_tmp, Dat_EMG_tmp(ch_i,:), mks_type(fce_i), 'Color', col_type(dist_i,:)) % channel 1  
                    end
                end
                
                if(ifplot)
                    % for position plot, plot original data line, the phold
                    % and fhold 
                    
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t, dat.pos);
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t(dat_idx_phold), dat.pos(:,dat_idx_phold), 'g.'); % position hold 
                    plot(subplot(axh{1}(fce_i,dist_i)),dat.t(dat_idx_fhold), dat.pos(:,dat_idx_fhold), 'g.'); % force hold
                    yline(subplot(axh{1}(fce_i,dist_i)),mean(dat.pos(:,dat_idx_fhold),2));

                    % for force plot, plot original data line, the average 
                    % at fhold
                    
                    plot(subplot(axh{2}(fce_i,dist_i)),dat.t, dat.fce);
                    plot(subplot(axh{2}(fce_i,dist_i)),dat.t(dat_idx_fhold), dat.fce(:,dat_idx_fhold), 'g.'); % force hold
                    yline(subplot(axh{2}(fce_i,dist_i)),mean(dat.fce(:,dat_idx_fhold),2));

                    % for emg plot, plot original data line, and the
                    % average at fhold
                    % option 1, plot the all channels with all trials 
                    for ch_i = 1:8
                        % also plot the maximum (minimum) at the release
                        % time 
                        for trial_i = 1:9 
%                             plot(subplot(axh{2+ch_i}(fce_i,dist_i)),dat.t, reshape(dat.emg(ch_i,trial_i,:),[1, size(dat.emg,3)]), 'color', col_type(ch_i,:));
%                             plot(subplot(axh{2+ch_i}(fce_i,dist_i)),dat.t(dat_idx_fhold), reshape(dat.emg(ch_i,trial_i,dat_idx_fhold),[1, sum(dat_idx_fhold)]), '.', 'MarkerEdgeColor', col_type(ch_i,:)); % force hold
%                             plot(subplot(axh{2+ch_i}(fce_i,dist_i)),dat.t(dat_idx_phold), reshape(dat.emg(ch_i,trial_i,dat_idx_phold),[1, sum(dat_idx_phold)]), '.', 'MarkerEdgeColor', col_type(ch_i,:)); % position hold
%                             yline(subplot(axh{2+ch_i}(fce_i,dist_i)),mean(dat.emg(ch_i,trial_i,dat_idx_fhold)), 'color', col_type(ch_i,:));
%                             yline(subplot(axh{2+ch_i}(fce_i,dist_i)),mean(dat.emg(ch_i,trial_i,dat_idx_fhold)), 'color', col_type(ch_i,:));
                            lnh(1) = plot(subplot(axh{2+ch_i}(fce_i,dist_i)),dat.t, reshape(dat.emg(ch_i,trial_i,:),[1, size(dat.emg,3)]), 'color', col_type(ch_i,:));
                            plot(subplot(axh{2+ch_i}(fce_i,dist_i)),dat.t(dat_idx_fhold), reshape(dat.emg(ch_i,trial_i,dat_idx_fhold),[1, sum(dat_idx_fhold)]), '.', 'MarkerEdgeColor', col_type(1,:)); % force hold
                            plot(subplot(axh{2+ch_i}(fce_i,dist_i)),dat.t(dat_idx_phold), reshape(dat.emg(ch_i,trial_i,dat_idx_phold),[1, sum(dat_idx_phold)]), '.', 'MarkerEdgeColor', col_type(2,:)); % position hold
                            lnh(2) = yline(subplot(axh{2+ch_i}(fce_i,dist_i)),mean(dat.emg(ch_i,trial_i,dat_idx_fhold)), 'color', col_type(1,:));
                            lnh(3) = yline(subplot(axh{2+ch_i}(fce_i,dist_i)),mean(dat.emg(ch_i,trial_i,dat_idx_phold)), 'color', col_type(2,:));
                            idx_tmp = Dat_EMG_aflidx(subj_i,dir_i,fce_i,dist_i,trial_i,ch_i);
                            val_tmp = Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,trial_i,ch_i);
                            t_tmp = dat.t(idx_tmp);
%                             plot(subplot(axh{2+ch_i}(fce_i,dist_i)),t_tmp, val_tmp, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', col_type(ch_i,:), 'MarkerFaceColor', [0.5 0.5 0.5]); % force hold
                            lnh(4) = plot(subplot(axh{2+ch_i}(fce_i,dist_i)),t_tmp, val_tmp, 'o', 'MarkerSize', 5, 'MarkerEdgeColor', col_type(ch_i,:), 'MarkerFaceColor', [0.5 0.5 0.5]); % force hold
                            legend(lnh, {'raw data', 'force hold average', 'position hold average', 'peak after release'});
                            xlabel('time');
                            ylabel('EMG value');
                        end
                    end
                    % option 2, plot 8th channel, different trial different line
%                     legend_arr = cell(9,1); 
%                     for trial_i = 1:9
%                         plot(subplot(axh{3}(fce_i,dist_i)),dat.t, reshape(dat.emg(8,trial_i,:),1, size(dat.emg,3)), 'linewidth', 2);
%                         legend_arr{trial_i} = ['ss' num2str(data_index_ss(subj_i,dir_i,fce_i,dist_i,trial_i)) 'tr' num2str(data_index_tr(subj_i,dir_i,fce_i,dist_i,trial_i))];
%                     end
%                     legend(legend_arr);
                end
            end
        end
    end
end
%% plot out the figure
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
    'elbow flexor', 'elbow extensor', ...
    'anterior deltoid', 'posterior deltiod', ...
    'pectoralis', 'Trapezius'};
clear axh
for subj_i = 1:4
    R2_list{subj_i} = [];
    for dir_i = 1:4
        fh = figure('name', ['EMG tuning subj' num2str(subj_i) 'dirc' num2str(dir_i)], ...
            'position', [0 0 600, 1440]);
        for ch_i = 1:8 % EMG channels
            for indpv = 1:3 % 1. Force; 2. x; 3. K
                axh(ch_i, indpv) = subplot(8, 3, (ch_i-1)*3+indpv);
                % for each force and target distance, plot the regressed
                % EMG activity with each independent variable
                hold on;
                for fce_i = 1:3
                    for dist_i = 1:3
                        emg_afl = reshape(Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i),1,9);
                        emg_fhd = reshape(Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i),1,9);
                        y_pts = emg_afl - emg_fhd;
                        switch indpv
                            case 1 % Force
                                x_pts = reshape(Dat_f(subj_i,dir_i,fce_i,dist_i,:),1,9);
                                xlabel_str = 'Force (N)';
                            case 2 % x
                                x_pts = reshape(Dat_x(subj_i,dir_i,fce_i,dist_i,:),1,9);
                                xlabel_str = 'displacement (m)';
                            case 3 % K
                                x_pts = reshape(Dat_k(subj_i,dir_i,fce_i,dist_i,:),1,9);
                                xlabel_str = 'Stiffness (N/m)';
                        end
                        plot(axh(ch_i,indpv), x_pts, y_pts, ...
                            'Marker', mks_type(fce_i), ...
                            'Linestyle', 'none', ...
                            'Color', col_type(dist_i+4,:));
                        % display to see the outlairs
                        y_pts;
                        reshape([data_index_ss(subj_i,dir_i,fce_i,dist_i,:,1)],1,9);
                        reshape([data_index_tr(subj_i,dir_i,fce_i,dist_i,:,1)],1,9);
                    end
                end


                switch indpv 
                    case 1
%                         x = Dat_f_cell{subj_i,dir_i};
                        x = reshape(Dat_f(subj_i,dir_i,:,:,:),1,81);
                        ylabel_str = {emg_pair_label{ch_i}; 'Delta EMG'}; 
                    case 2
%                         x = Dat_x_cell{subj_i,dir_i};
                        x = reshape(Dat_x(subj_i,dir_i,:,:,:),1,81);
                        ylabel_str = 'Delta EMG'; 

                    case 3
%                         x = Dat_k_cell{subj_i,dir_i};
                        x = reshape(Dat_k(subj_i,dir_i,:,:,:),1,81);
                        ylabel_str = 'Delta EMG'; 
                end

                xlabel(xlabel_str);
                ylabel(ylabel_str);

%                 y  = Dat_EMG_cell{subj_i,dir_i}(ch_i,:);  
                y = reshape(Dat_EMG_afl(subj_i,dir_i,:,:,:,ch_i)-Dat_EMG_fhd(subj_i,dir_i,:,:,:,ch_i),1,81);

                y_validx = ~isnan(y); 
                % do regression part
                x = x(y_validx);
                y = y(y_validx); 
%                 plot(x, y, '.', 'color', [0.2 0.2 0.2]); % sanity check making sure the same dtpts

                [P, S] = polyfit(x,y,1);
                R2 = 1-(S.normr/norm(y-mean(y)))^2;
                R2_list{subj_i} = [R2_list{subj_i}, R2];
                R2_qualified = R2>0.2;
                [b, bint, r, rint, stats] = regress(y', [ones(size(x)); x]' );
                p_val = stats(3);
                R2_qualified = p_val<128; %(bufforoni correction)
                yfit1 = polyval(P,x); 

                plot(x, yfit1, 'linewidth', 2, 'color', [0.7 0.7 0.7]); 
                if (R2_qualified)
                    title(['R^2:' num2str(R2)], 'color', [0 0 0]);
                else
                    title(['R^2:' num2str(R2)], 'color', [0.5 0.5 0.5]);
                end

            end
        end
        sgtitle(['subj' num2str(subj_i) 'direc' num2str(dir_i), 'EMG release peak']);
        exportgraphics(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_regress.pdf', 'Append', true);
    end
end

%% 
Dat_EMG_afl; 
Dat_f;
Dat_x;
Dat_k;
Dat_subj = zeros(4,4,3,3,9); 
Dat_dir = zeros(4,4,3,3,9);
Dat_fce_cmd = zeros(4,4,3,3,9); 
Dat_dis_cmd = zeros(4,4,3,3,9); 
for subj_i = 1:4 
    for dir_i = 1:4
        Dat_subj(subj_i,dir_i,:,:,:) = subj_i;
        Dat_dir(subj_i,dir_i,:,:,:) = dir_i;
        fce_list = [15 20 25];
        dist_list = [0.025 0.05 0.075];
        for fce_i = 1:3
            for dist_i = 1:3
                Dat_fce_cmd(subj_i,dir_i,fce_i,dist_i,:) = fce_list(fce_i);
                Dat_dis_cmd(subj_i,dir_i,fce_i,dist_i,:) = dist_list(dist_i);
            end
        end
    end
end

% subject 
% direction 
%EMG ~ 1+ Subj + Dir + Force + displacement %... option1
%EMG ~ 1+ Subj + Dir + Force + k            %... option2
blr.f = Dat_f(:); 
blr.x = Dat_x(:);
blr.k = Dat_k(:);
blr.subj = Dat_subj(:);
blr.dir = Dat_dir(:);
blr.fce_cmd = Dat_fce_cmd(:);
blr.dist_cmd = Dat_dis_cmd(:);

%%%%%%%%%%%%%%% use function anovan (n-way anova
p = anovan(blr.emg1, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});
p = anovan(blr.emg2, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});
p = anovan(blr.emg3, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});
p = anovan(blr.emg4, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});
p = anovan(blr.emg5, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});
p = anovan(blr.emg6, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});
p = anovan(blr.emg7, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});
p = anovan(blr.emg8, [blr.subj,blr.dir,blr.fce_cmd,blr.dist_cmd], 'varnames',{'subj','direction', 'force', 'distance'});

%%%%%%%%%%%%%%% consider using the 2-way anova 
p_mat_fce = zeros(4,4,8); % subj_num*dir_num*muscle_num;
p_mat_dist = zeros(4,4,8); % subj_num*dir_num*muscle_num;
emgtmp_cell = cell(1,8);
for subj_i = 1:4
     for dir_i = 1:4
        fcetmp = reshape(Dat_fce_cmd(subj_i,dir_i,:,:,:),3,3,9);
        disttmp = reshape(Dat_dis_cmd(subj_i,dir_i,:,:,:),3,3,9);

        for chi = 1:8
%         emgtmp_cell{chi} = reshape(Dat_EMG_fhd(subj_i,dir_i,:,:,:,chi),3,3,9);
        emgtmp_cell{chi} = reshape(Dat_EMG_afl(subj_i,dir_i,:,:,:,chi),3,3,9);
          emg1 = [reshape(dattmp(1,:,:), 3, 9)'; reshape(dattmp(2,:,:), 3, 9)'; reshape(dattmp(3,:,:), 3, 9)'];
%         [p, tbl, stats] = anova2(emg1, 3); % anova2 not work due to have nan values 
        [p, tbl, stats] = anovan(emgtmp_cell{chi}(:), [fcetmp(:), disttmp(:)], 'varnames', {'force','distance'});
        
%         [p, tbl, stats] = anovan(emgtmp{chi}(:), [fcetmp(:), disttmp(:), fcetmp(:).*disttmp(:)], 'varnames', {'force','distance', 'force*distance'});
        p_mat_fce(subj_i,dir_i,chi) = tbl{2,7};
        p_mat_dist(subj_i,dir_i,chi) = tbl{3,7};
        end
     end
end
% for the force p-vals 
p_mat_fce_2d = [reshape(p_mat_fce(1,:,:), size(p_mat_fce,[2,3])); 
                reshape(p_mat_fce(2,:,:), size(p_mat_fce,[2,3]));
                reshape(p_mat_fce(3,:,:), size(p_mat_fce,[2,3]));
                reshape(p_mat_fce(4,:,:), size(p_mat_fce,[2,3]));]
% for the dsitance p-vals
p_mat_dist_2d = [reshape(p_mat_dist(1,:,:), size(p_mat_dist,[2,3])); 
                reshape(p_mat_dist(2,:,:), size(p_mat_dist,[2,3]));
                reshape(p_mat_dist(3,:,:), size(p_mat_dist,[2,3]));
                reshape(p_mat_dist(4,:,:), size(p_mat_dist,[2,3]));]

% % %%%%%%%%%%%%%%% use function anova
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,1);
% % blr.emg1 = Dat_EMG_tmp(:); % fhd
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,2);
% % blr.emg2 = Dat_EMG_tmp(:); % fhd
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,3);
% % blr.emg3 = Dat_EMG_tmp(:); % fhd
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,4);
% % blr.emg4 = Dat_EMG_tmp(:); % fhd
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,5);
% % blr.emg5 = Dat_EMG_tmp(:); % fhd
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,6);
% % blr.emg6 = Dat_EMG_tmp(:); % fhd
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,7);
% % blr.emg7 = Dat_EMG_tmp(:); % fhd
% % Dat_EMG_tmp = Dat_EMG_afl(:,:,:,:,:,8);
% % blr.emg8 = Dat_EMG_tmp(:); % fhd
% % 
% % tbl = table(blr.subj, blr.dir, blr.f, blr.x, blr.k, blr.fce_cmd, blr.dist_cmd, ...
% %     blr.emg1, blr.emg2, blr.emg3, blr.emg4, blr.emg5, blr.emg6, blr.emg7, blr.emg8, ...
% %     'VariableNames', {'subj', 'dir', 'fce', 'disp', 'stiffness', 'fce_cmd', 'dist_cmd',...
% %     'emg1', 'emg2', 'emg3', 'emg4', 'emg5', 'emg6', 'emg7', 'emg8'});
% % tbl.subj = categorical(tbl.subj);
% % tbl.dir = categorical(tbl.dir);
% % tbl.fce_cmd = categorical(tbl.fce_cmd);
% % tbl.dist_cmd = categorical(tbl.dist_cmd);
% % 
% % % mdl1 = fitlm(tbl, 'emg1 ~ subj + dir + fce + stiffness');
% % % mdl2 = fitlm(tbl, 'emg2 ~ subj + dir + fce + stiffness');
% % % mdl3 = fitlm(tbl, 'emg3 ~ subj + dir + fce + stiffness');
% % % mdl4 = fitlm(tbl, 'emg4 ~ subj + dir + fce + stiffness');
% % % mdl5 = fitlm(tbl, 'emg5 ~ subj + dir + fce + stiffness');
% % % mdl6 = fitlm(tbl, 'emg6 ~ subj + dir + fce + stiffness');
% % % mdl7 = fitlm(tbl, 'emg7 ~ subj + dir + fce + stiffness');
% % % mdl8 = fitlm(tbl, 'emg8 ~ subj + dir + fce + stiffness');
% % 
% % mdl1 = fitlm(tbl, 'emg1 ~ subj + dir + fce_cmd + dist_cmd');
% % mdl2 = fitlm(tbl, 'emg2 ~ subj + dir + fce_cmd + dist_cmd');
% % mdl3 = fitlm(tbl, 'emg3 ~ subj + dir + fce_cmd + dist_cmd');
% % mdl4 = fitlm(tbl, 'emg4 ~ subj + dir + fce_cmd + dist_cmd');
% % mdl5 = fitlm(tbl, 'emg5 ~ subj + dir + fce_cmd + dist_cmd');
% % mdl6 = fitlm(tbl, 'emg6 ~ subj + dir + fce_cmd + dist_cmd');
% % mdl7 = fitlm(tbl, 'emg7 ~ subj + dir + fce_cmd + dist_cmd');
% % mdl8 = fitlm(tbl, 'emg8 ~ subj + dir + fce_cmd + dist_cmd');
% % 
% % anova(mdl1, 'summary')
% % anova(mdl2, 'summary')
% % anova(mdl3, 'summary')
% % anova(mdl4, 'summary')
% % anova(mdl5, 'summary')
% % anova(mdl6, 'summary')
% % anova(mdl7, 'summary')
% % anova(mdl8, 'summary')
% % 
% % mdl1 = fitlm(tbl, 'emg1 ~ subj + dir + fce + disp');
% % mdl2 = fitlm(tbl, 'emg2 ~ subj + dir + fce + disp');
% % mdl3 = fitlm(tbl, 'emg3 ~ subj + dir + fce + disp');
% % mdl4 = fitlm(tbl, 'emg4 ~ subj + dir + fce + disp');
% % mdl5 = fitlm(tbl, 'emg5 ~ subj + dir + fce + disp');
% % mdl6 = fitlm(tbl, 'emg6 ~ subj + dir + fce + disp');
% % mdl7 = fitlm(tbl, 'emg7 ~ subj + dir + fce + disp');
% % mdl8 = fitlm(tbl, 'emg8 ~ subj + dir + fce + disp');
% % 
% % anova(mdl1, 'summary')
% % anova(mdl2, 'summary')
% % anova(mdl3, 'summary')
% % anova(mdl4, 'summary')
% % anova(mdl5, 'summary')
% % anova(mdl6, 'summary')
% % anova(mdl7, 'summary')
% % anova(mdl8, 'summary')
% % 
% % mdl1 = fitlm(tbl, 'emg1 ~ subj + dir + fce');
% % mdl2 = fitlm(tbl, 'emg2 ~ subj + dir + fce');
% % mdl3 = fitlm(tbl, 'emg3 ~ subj + dir + fce');
% % mdl4 = fitlm(tbl, 'emg4 ~ subj + dir + fce');
% % mdl5 = fitlm(tbl, 'emg5 ~ subj + dir + fce');
% % mdl6 = fitlm(tbl, 'emg6 ~ subj + dir + fce');
% % mdl7 = fitlm(tbl, 'emg7 ~ subj + dir + fce');
% % mdl8 = fitlm(tbl, 'emg8 ~ subj + dir + fce');
% % 
% % anova(mdl1, 'summary')
% % anova(mdl2, 'summary')
% % anova(mdl3, 'summary')
% % anova(mdl4, 'summary')
% % anova(mdl5, 'summary')
% % anova(mdl6, 'summary')
% % anova(mdl7, 'summary')
% % anova(mdl8, 'summary')
% % 
% % mdl1 = fitlm(tbl, 'emg1 ~ subj + dir + stiffness');
% % mdl2 = fitlm(tbl, 'emg2 ~ subj + dir + stiffness');
% % mdl3 = fitlm(tbl, 'emg3 ~ subj + dir + stiffness');
% % mdl4 = fitlm(tbl, 'emg4 ~ subj + dir + stiffness');
% % mdl5 = fitlm(tbl, 'emg5 ~ subj + dir + stiffness');
% % mdl6 = fitlm(tbl, 'emg6 ~ subj + dir + stiffness');
% % mdl7 = fitlm(tbl, 'emg7 ~ subj + dir + stiffness');
% % mdl8 = fitlm(tbl, 'emg8 ~ subj + dir + stiffness');
% % 
% % anova(mdl1, 'summary')
% % anova(mdl2, 'summary')
% % anova(mdl3, 'summary')
% % anova(mdl4, 'summary')
% % anova(mdl5, 'summary')
% % anova(mdl6, 'summary')
% % anova(mdl7, 'summary')
% % anova(mdl8, 'summary')
% % 
% % mdl1 = fitlm(tbl, 'emg1 ~ stiffness');
% % mdl2 = fitlm(tbl, 'emg2 ~ stiffness');
% % mdl3 = fitlm(tbl, 'emg3 ~ stiffness');
% % mdl4 = fitlm(tbl, 'emg4 ~ stiffness');
% % mdl5 = fitlm(tbl, 'emg5 ~ stiffness');
% % mdl6 = fitlm(tbl, 'emg6 ~ stiffness');
% % mdl7 = fitlm(tbl, 'emg7 ~ stiffness');
% % mdl8 = fitlm(tbl, 'emg8 ~ stiffness');
% % 
% % anova(mdl1, 'summary')
% % anova(mdl2, 'summary')
% % anova(mdl3, 'summary')
% % anova(mdl4, 'summary')
% % anova(mdl5, 'summary')
% % anova(mdl6, 'summary')
% % anova(mdl7, 'summary')
% % anova(mdl8, 'summary')
% % 
% % 
% % mdl1 = fitlm(tbl, 'emg1 ~ subj');
% % mdl2 = fitlm(tbl, 'emg2 ~ subj');
% % mdl3 = fitlm(tbl, 'emg3 ~ subj');
% % mdl4 = fitlm(tbl, 'emg4 ~ subj');
% % mdl5 = fitlm(tbl, 'emg5 ~ subj');
% % mdl6 = fitlm(tbl, 'emg6 ~ subj');
% % mdl7 = fitlm(tbl, 'emg7 ~ subj');
% % mdl8 = fitlm(tbl, 'emg8 ~ subj');
% % 
% % anova(mdl1, 'summary')
% % anova(mdl2, 'summary')
% % anova(mdl3, 'summary')
% % anova(mdl4, 'summary')
% % anova(mdl5, 'summary')
% % anova(mdl6, 'summary')
% % anova(mdl7, 'summary')
% % anova(mdl8, 'summary')

%% plot out the EMG using bar plot 
xx = [25 50 75];
x = [20 25 30; 45 50 55; 70 75 80];
ff = [15 20 25];
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
    'elbow flexor', 'elbow extensor', ...
    'anterior deltoid', 'posterior deltiod', ...
    'pectoralis', 'Trapezius'};
for subj_i = 1:4 
    clear axh
    for dir_i = 1:4
        fh(subj_i,dir_i) = figure('position', [0 0 1200 200]);
        
        for ch_i = 1:8
            axh(dir_i, ch_i) = subplot(1,8,(ch_i));
            emg_tmp_means = zeros(3,3);
            emg_tmp_stds = zeros(3,3);
            for fce_i = 1:3
                for dist_i = 1:3
%                     emg_tmp_means(fce_i,dist_i) = mean(Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
%                     emg_tmp_stds(fce_i,dist_i) = std(Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                     emg_tmp_means(fce_i,dist_i) = mean(Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                     emg_tmp_stds(fce_i,dist_i) = std(Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                end
            end
            lgd_fh = bar(xx, emg_tmp_means'); hold on;
%             errorbar([xx; xx; xx]', emg_tmp_means, emg_tmp_stds); % is it right?
            errorbar(x, emg_tmp_means', emg_tmp_stds', 'LineWidth',2, 'LineStyle', 'none', 'Color', 'k'); % is it right?
            
            title(emg_pair_label{ch_i});
            xlabel('distance (mm)'); ylabel('EMG (normalized)'); grid on;
        end

        if dir_i == 4
            for ch_i = 1:8
                linkaxes(axh(:,ch_i), 'y');
                if (ch_i == 8)
                    legend(lgd_fh, {'15N', '20N', '25N'});
                end
            end
        end
    sgtitle(['EMG activities of subject ' num2str(subj_i) ' direction' num2str(dir_i)]);
    end
end


%% Do a t-teset for the f-hold and p-hold values on each of the pairs 
% for each muscle activity 
clear axh
figure(); 
% for ch_i = 1:8
%     fhd_val = Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i)
%     phd_val = Dat_EMG_phd(subj_i,dir_i,fce_i,dist_i,:,ch_i)
%     fhd_val = Dat_EMG_fhd(:,:,:,:,:,ch_i);
%     phd_val = Dat_EMG_phd(:,:,:,:,:,ch_i);
    fhd_val = Dat_EMG_fhd(:,:,:,:,:,:);
    phd_val = Dat_EMG_phd(:,:,:,:,:,:);
    chg_val = phd_val - fhd_val;
    axh(ch_i) = subplot(8,1,ch_i);
    
    [h,p] = ttest(fhd_val(:), phd_val(:));

%     figure(); hold on;
%     histogram(fhd_val(:)); 
%     histogram(phd_val(:));
%     legend('force hold', 'poistion hold');
    figure();
    histogram(chg_val(:));
    xlim([-0.5 0.5]);
    grid on;
    xlabel('EMG_{phold} - EMG_{fhold}'); 
    ylabel('times');
    title('EMG change after release');
% end
linkaxes(axh, 'x');
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if response tunes with stiffness levels 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if EMG levels change before and after release 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the fitting of F and EMG

dat_EMG_tmp = Dat_EMG_fhd(ch_i,:);
dat_EMG_tmp_mean = mean(Dat_EMG_fhd(ch_i,:));
dat_EMG_tmp_std = std(Dat_EMG_fhd(ch_i,:));
dat_EMG_tmp_validx = dat_EMG_tmp > dat_EMG_tmp_mean-3*dat_EMG_tmp_std & ...
                    dat_EMG_tmp < dat_EMG_tmp_mean+3*dat_EMG_tmp_std;
dat_f_tmp = Dat_f; 

x = dat_EMG_tmp(dat_EMG_tmp_validx); 
y = dat_f_tmp(dat_EMG_tmp_validx); 
scatter(x, y); 
refline;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chenugnag: I saw the EMG data looks filled in clusters. I'm wondering
% whether they are the same condition or different conditions 
