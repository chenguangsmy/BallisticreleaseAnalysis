% investigate the EMG values with the stiffness values. 
% Just plot them as the regression of either displacement (x), force (f) or
% stiffness (k)

% function fh = plotEMG_release_dir1_acrosssubject(obj)
%%%%%%%%%%%%%%%%%%% do something the same with the class function 
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
        ifplotEMG_scatter = 0;
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
                dat_avg.emg = mean(dat.emg(:,:,dat_idx_fhold), 3);

                %         hold on;
                %         plot(dat_avg.f, dat_avg.emg(1,:), '.');
                %         plot(dat_avg.x, dat_avg.emg(1,:), '.');
                %         plot(dat_avg.f./dat_avg.x, dat_avg.emg(1,:), '.');

                Dat_f = [Dat_f; dat_avg.f];
                Dat_x = [Dat_x; dat_avg.x];
%                 Dat_k = [Dat_k; dat_avg.f./dat_avg.x];
                Dat_k = [Dat_k reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_b = [Dat_b reshape(Results(subj_i,dir_i).B_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_EMG = [Dat_EMG dat_avg.emg];
                Dat_EMG_tmp = dat_avg.emg;
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

% fitting: a structure of multiple slope, intercept, R2... 
Dat_f = nan(4,4,3,3,9);
Dat_x = nan(4,4,3,3,9);
Dat_k = nan(4,4,3,3,9);
Dat_b = nan(4,4,3,3,9);
Dat_EMG = nan(4,4,3,3,9,8);
mks_type = '*so';
for subj_i  = 1:4    
    for dir_i   = 1:4
        pert_i  = 1;
        % across trials, get the EMG values as well as the x, f, k values
        % temperarily use k values as f/x, later substitute it with the estimation
        % values.

        % dat_avg % the data of averaged values; (for force and displacement, use
        % the absolute values)
        %
        time_window_fhold = [-0.5 0]; % before release 0.5s take average
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
                dat_avg.x = abs(dat.pos(:,dat_idx_phold) - mean(dat.pos(:,dat_idx_fhold),2));
                dat_avg.f = abs(mean(dat.fce(:,dat_idx_fhold), 2));
                dat_avg.emg = mean(dat.emg(:,:,dat_idx_fhold), 3);

                %         hold on;
                %         plot(dat_avg.f, dat_avg.emg(1,:), '.');
                %         plot(dat_avg.x, dat_avg.emg(1,:), '.');
                %         plot(dat_avg.f./dat_avg.x, dat_avg.emg(1,:), '.');

                Dat_f_pool = [Dat_f_pool; dat_avg.f];
                Dat_x_pool = [Dat_x_pool; dat_avg.x];
                Dat_k_pool = [Dat_k_pool reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_b_pool = [Dat_b_pool reshape(Results(subj_i,dir_i).B_up_tr(fce_i, dist_i, :), [1,9])];
                Dat_EMG_pool = [Dat_EMG_pool dat_avg.emg];

                Dat_f(subj_i,dir_i,fce_i,dist_i,:) = dat_avg.f;
                Dat_x(subj_i,dir_i,fce_i,dist_i,:) = dat_avg.x;
                Dat_k(subj_i,dir_i,fce_i,dist_i,:) = reshape(Results(subj_i,dir_i).K_up_tr(fce_i, dist_i, :), [1,9]);
                Dat_b(subj_i,dir_i,fce_i,dist_i,:) = reshape(Results(subj_i,dir_i).B_up_tr(fce_i, dist_i, :), [1,9]);
                fce_i
                dist_i
                dat_avg.emg
                Dat_EMG(subj_i,dir_i,fce_i,dist_i,:,:) = dat_avg.emg';

                Dat_EMG_tmp = Dat_EMG(subj_i,dir_i,fce_i,dist_i,:,:);
                Dat_k_tmp = reshape(Dat_k(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                Dat_f_tmp = reshape(Dat_f(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                Dat_x_tmp = reshape(Dat_x(subj_i,dir_i,fce_i,dist_i,:), [1,9]);
                
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
        Dat_f_cell{subj_i,dir_i} = Dat_f_pool;
        Dat_x_cell{subj_i,dir_i} = Dat_x_pool;
        Dat_k_cell{subj_i,dir_i} = Dat_k_pool;
        Dat_b_cell{subj_i,dir_i} = Dat_b_pool;
        Dat_EMG_cell{subj_i,dir_i} = Dat_EMG_pool;
    end
end
%% % plot out the figure
for subj_i = 1:4
    for dir_i = 1:4
        figure('name', ['EMG tuning subj' num2str(subj_i) 'dirc' num2str(dir_i)], ...
            'position', [0 0 600, 1440]);
        for ch_i = 1:8 % EMG channels
            for indpv = 1:3 % 1. Force; 2. x; 3. K
                axh(ch_i, indpv) = subplot(8, 3, (ch_i-1)*3+indpv);
                % for each force and target distance, plot the regressed
                % EMG activity with each independent variable
                hold on;
                for fce_i = 1:3
                    for dist_i = 1:3
                        y_pts = reshape(Dat_EMG(subj_i,dir_i,fce_i,dist_i,:,ch_i),1,9);
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
                    case 2
%                         x = Dat_x_cell{subj_i,dir_i};
                        x = reshape(Dat_x(subj_i,dir_i,:,:,:),1,81);
                    case 3
%                         x = Dat_k_cell{subj_i,dir_i};
                        x = reshape(Dat_k(subj_i,dir_i,:,:,:),1,81);
                end

                ylabel_str = 'normalized EMG'; 
                xlabel(xlabel_str);
                ylabel(ylabel_str);

%                 y  = Dat_EMG_cell{subj_i,dir_i}(ch_i,:);  
                y = reshape(Dat_EMG(subj_i,dir_i,:,:,:,ch_i),1,81);

                y_validx = ~isnan(y); 
                % do regression part
                x = x(y_validx);
                y = y(y_validx); 
%                 plot(x, y, '.', 'color', [0.2 0.2 0.2]); % sanity check making sure the same dtpts

                [P, S] = polyfit(x,y,1);
                R2 = 1-(S.normr/norm(y-mean(y)))^2;
                R2_qualified = R2>0.2;
                yfit1 = polyval(P,x); 

                plot(x, yfit1, 'linewidth', 2, 'color', [0.7 0.7 0.7]); 
                if (R2_qualified)
                    title(['R^2:' num2str(R2)], 'color', [0 0 0]);
                else
                    title(['R^2:' num2str(R2)], 'color', [0.5 0.5 0.5]);
                end

            end
        end
        sgtitle(['subj' num2str(subj_i) 'direc' num2str(dir_i)]);
    end
end

%% plot F, x, EMG on each of the muscle
if (1)
    for subj_i = 1:4
        Dat_f = Dat_f_cell{subj_i,dir_i};
        Dat_x = Dat_x_cell{subj_i,dir_i};
        Dat_k = Dat_k_cell{subj_i,dir_i};
        Dat_b = Dat_b_cell{subj_i,dir_i};
        Dat_EMG = Dat_EMG_cell{subj_i,dir_i};
        figure('name', 'single EMG with F, x and k', 'position', [0 0 1200 300]);

%         subplot(4,1,1); hold on;
        for ch_i = 1:8
            figure(); hold on;
            dat_EMG_tmp = Dat_EMG(ch_i,:);
            dat_EMG_tmp_mean = mean(Dat_EMG(ch_i,:), 'omitnan');
            dat_EMG_tmp_std = std(Dat_EMG(ch_i,:), 'omitnan');
            dat_EMG_tmp_validx = dat_EMG_tmp > (dat_EMG_tmp_mean-3*dat_EMG_tmp_std) & ...
                dat_EMG_tmp < (dat_EMG_tmp_mean+3*dat_EMG_tmp_std);
%             dat_f_tmp = Dat_f;
            dat_f_tmp = Dat_k;
            plot(dat_EMG_tmp, dat_f_tmp, '.', 'Color', col_type(ch_i,:));
%             refline;
            x = dat_EMG_tmp(dat_EMG_tmp_validx);
            y = dat_f_tmp(dat_EMG_tmp_validx);
            plot(x, y, 'o', 'Color', col_type(ch_i,:));
            refline;

%             plot(Dat_EMG(ch_i,:),Dat_f, '.', 'Color', col_type(ch_i,:));  

%             p = polyfit(Dat_EMG(ch_i,:),Dat_f, 1);
%             tmp_x = 10:0.1:30;
%             tmp_y = p(1) * tmp_x + p(2); 
%             lnh(ch_i) = plot(tmp_x, tmp_y, 'Color', col_type(ch_i,:), 'LineWidth', 1.5);
        end
        lnh = refline;
        legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
        xlabel('Force exertion (N)');
        ylabel('EMG values perportion to MVF EMG');
        title(['EMG corresponding to F, subj' num2str(subj_i)]);


        %%%
        subplot(4,1,2); hold on;
        for ch_i = 1:8
            plot(Dat_EMG(ch_i,:), Dat_x, '.', 'Color', col_type(ch_i,:));
        end
        lnh = refline;
        legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
        xlabel('actual displacement (m)');
        ylabel('EMG values perportion to MVF EMG');
        title(['EMG corresponding to x, subj' num2str(subj_i)]);

        %%%
        subplot(4,1,3); hold on;
        for ch_i = 1:8
            plot(Dat_EMG(ch_i,:), Dat_k, '.', 'Color', col_type(ch_i,:));
        end
        lnh = refline;
        legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
        xlabel('K_{est} (N/m)');
        ylabel('EMG values perportion to MVF EMG');
        title(['EMG corresponding to K, subj' num2str(subj_i)]);

        subplot(4,1,4); hold on;
        for ch_i = 1:8
            plot(Dat_EMG(ch_i,:), Dat_b, '.', 'Color', col_type(ch_i,:));
        end
        lnh = refline;
        legend(lnh, {'ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8'});
        xlabel('B_{est} (N/m)');
        ylabel('EMG values perportion to MVF EMG');
        title(['EMG corresponding to B, subj' num2str(subj_i)]);
    end
end
% Think: how to represent the strength of tuning in the muscles??? 

%% plot EMG summation and subtraction to see 'stiffness' and 'force' values
if(0)
    for subj_i = 1:4
        Dat_f = Dat_f_cell{subj_i,dir_i};
        Dat_x = Dat_x_cell{subj_i,dir_i};
        Dat_k = Dat_k_cell{subj_i,dir_i};
        Dat_b = Dat_b_cell{subj_i,dir_i};
        Dat_EMG = Dat_EMG_cell{subj_i,dir_i};
        figure('name', ['EMG summations dir' num2str(dir_i)], 'position', [0 0 1200 300]);

        subplot(1,3,1); hold on;
        for emg_pr_i = 1:4
            %     EMG_m = Dat_EMG((emg_pr_i-1)*2+1,:) - Dat_EMG((emg_pr_i-1)*2+2,:);
            EMG_p = Dat_EMG((emg_pr_i-1)*2+1,:) + Dat_EMG((emg_pr_i-1)*2+2,:);
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
            EMG_p = Dat_EMG((emg_pr_i-1)*2+1,:) + Dat_EMG((emg_pr_i-1)*2+2,:);
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
            EMG_p = Dat_EMG((emg_pr_i-1)*2+1,:) + Dat_EMG((emg_pr_i-1)*2+2,:);
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
        Dat_EMG = Dat_EMG_cell{subj_i,dir_i};
        figure('name', ['EMG subtractions dir' num2str(dir_i)], 'position', [0 0 1200 300]);

        subplot(1,3,1); hold on;
        for emg_pr_i = 1:4
            idx1 = (emg_pr_i-1)*2+1;
            idx2 = (emg_pr_i-1)*2+2;
            EMG_m = Dat_EMG(idx1,:)*subtraction_mat(dir_i,idx1) + Dat_EMG(idx2,:)*subtraction_mat(dir_i,idx2);
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
                    EMG_m = Dat_EMG((emg_pr_i-1)*2+1,:) - Dat_EMG((emg_pr_i-1)*2+2,:);
                case {3 4}
                    EMG_m = Dat_EMG((emg_pr_i-1)*2+2,:) - Dat_EMG((emg_pr_i-1)*2+1,:);
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
                    EMG_m = Dat_EMG((emg_pr_i-1)*2+1,:) - Dat_EMG((emg_pr_i-1)*2+2,:);
                case {3 4}
                    EMG_m = Dat_EMG((emg_pr_i-1)*2+2,:) - Dat_EMG((emg_pr_i-1)*2+1,:);
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

dat_EMG_tmp = Dat_EMG(ch_i,:);
dat_EMG_tmp_mean = mean(Dat_EMG(ch_i,:));
dat_EMG_tmp_std = std(Dat_EMG(ch_i,:));
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
