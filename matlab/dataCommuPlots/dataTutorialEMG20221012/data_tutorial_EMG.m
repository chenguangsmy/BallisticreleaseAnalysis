% This code shows the tidied up cells of EMG 

clear; close all; clc; 

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4310_4356_EMG.mat'); 

% variable comments
% data: 6-D cell, including data of each trial 
% data_index_ss: 6-D double, data come from which session; 
% data_index_tr: 6-D double, data come from which trial; 

% ... There are reapeating 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% HERE COMES THE PLOTTING ON THE EMG ACTIVITIES
%% 
% Plot1, plot single task condition EMG with displacement/velocity
obj.data = data;
obj.cond.subj = 1:size(data,1);
obj.cond.dir  = 1:size(data,2);
obj.cond.fce  = 1:size(data,3);
obj.cond.dist = 1:size(data,4);
obj.cond.trial= 1:size(data,5);
obj.cond.pert = 1:size(data,6);

clear dat

Fs = 500;
fh1 = figure();
col_type = colormap('lines');
close(fh1);
t_range = [-0.5 1];

emg_pair = [1 2 3 4 5 6 7 8];
emg_pair_label = {'wrist', 'wrist', ...
    'elbow', 'elbow', ...
    'deltoid', 'deltoid', ...
    'shoulder', 'shoulder'};
cols = 1;
rows = 2 + length(emg_pair)/2; % only plot position and the muscles
axh = zeros(rows, cols);

pert_i  = 1;
for subj_i  = 1
    for dir_i   = 2

        for fce_i = 2 % 1:3
            for dist_i = 2 %1:3
                fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
                trials_list = obj.cond.trial;
                % make enough space for data
                dat.row = length(trials_list);
                dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
                    obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
                dat.t = linspace(t_range(1), t_range(2), dat.col);
                dat.pos = nan(dat.row, dat.col);
                dat.fce = nan(dat.row, dat.col);
                dat.vel = nan(dat.row, dat.col);
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
                    dat.vel(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ov(1,index_t),dat.t,'spline');
                    for ch_i = 1:8
                        try
                            dat.emg(ch_i,trial_idx,:)=interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,index_t)',dat.t,'spline')';
%                             dat.emg(ch_i,trial_idx,:)=interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emgrtf(ch_i,index_t)',dat.t,'spline')';
                        catch
                            disp(['no EMG this condition! fce' num2str(fce_i) ' dist' num2str(dist_i)]);
                        end
                    end
                end

                % plot the mean and average of it
                %
                axh(1, 1) = subplot(rows,cols,1); hold on;       % position
                plot(dat.t, mean(dat.pos, 'omitnan'), ...
                    'Color', col_type(4+dist_i,:), ...
                    'LineWidth', 2);
                patch_x = [dat.t, dat.t(end:-1:1)];
                tmp1 = mean(dat.pos, 'omitnan')+std(dat.pos, 'omitnan');
                tmp2 = mean(dat.pos, 'omitnan')-std(dat.pos, 'omitnan');
                patch_y = [tmp1, tmp2(end:-1:1)];
                patch(patch_x, patch_y, ...
                    col_type(4+dist_i,:), ...
                    'FaceAlpha', 0.3, ...
                    'EdgeColor', 'none');
                ylabel('x (m)'); title('position');
                xlabel('t (s)');
                grid on;

                %
                axh(2, 1) = subplot(rows,cols,2); hold on;       % force
                plot(dat.t, mean(dat.vel, 'omitnan'), ...
                    'Color', col_type(4+dist_i,:), ...
                    'LineWidth', 2);
                patch_x = [dat.t, dat.t(end:-1:1)];
                tmp1 = mean(dat.vel, 'omitnan')+std(dat.vel,'omitnan');
                tmp2 = mean(dat.vel, 'omitnan')-std(dat.vel,'omitnan');
                patch_y = [tmp1, tmp2(end:-1:1)];
                patch(patch_x, patch_y, ...
                    col_type(4+dist_i,:), ...
                    'FaceAlpha', 0.3, ...
                    'EdgeColor', 'none');
                ylabel('v (m/s)'); title('velocity');
                xlabel('t (s)');
                grid on;

                %
                for muscle_i = 1:length(emg_pair)
                    if mod(muscle_i,2) == 1
                        axh(2+muscle_i, 1) = subplot(rows,cols,2 + ceil(muscle_i/2)); hold on;       % ch 3: elbow flexor
                        ag_ant_i = 0;
                    else
                        axh(2+muscle_i, 1) = subplot(rows,cols,2 + ceil(muscle_i/2)); hold on;       % ch 3: elbow flexor
                        ag_ant_i = 1;
                    end
                    lnh(ag_ant_i+1) = plot(dat.t, mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                        'Color', col_type(4+dist_i+ag_ant_i,:), ...
                        'LineWidth', 2);
                    %             ylabel('EMG (mV)');
                    ylabel('EMG (portion)');
                    title(emg_pair_label{emg_pair(muscle_i)});

                    patch_x = [dat.t, dat.t(end:-1:1)];
                    tmp1 = mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
                        std(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
                    tmp2 = mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
                        std(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
                    patch_y = [tmp1, tmp2(end:-1:1)];
                    patch(patch_x, patch_y, ...
                        col_type(4+dist_i+ag_ant_i,:), ...
                        'FaceAlpha', 0.3, ...
                        'EdgeColor', 'none');

                    xline(0.05, 'LineStyle','--');
                    xlabel('t (s)');
                    grid on;


                end
                legend(lnh, 'flexor', 'extensor');

                linkaxes(axh(:), 'x');
                xlim([-0.5 0.5]); % sec
                sgtitle(['EMG demo subj' num2str(subj_i), 'dir' num2str(dir_i) 'f' num2str(fce_i) 'd' num2str(dist_i)]);
            end

        end

    end
end


%% 
% plot2, single condition with agonist/antagonist muscle pair (each trace)

obj.data = data;
obj.cond.subj = 1:size(data,1);
obj.cond.dir  = 1:size(data,2);
obj.cond.fce  = 1:size(data,3);
obj.cond.dist = 1:size(data,4);
obj.cond.trial= 1:size(data,5);
obj.cond.pert = 1:size(data,6);

clear dat

for emg_muscle_pairs = 1:4
fh(emg_muscle_pairs) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
Fs = 500;
col_type = colormap('lines');
t_range = [-0.5 2];

emg_pair = (emg_muscle_pairs-1)*2 + [1 2];

emg_pair_label = {'wrist flexor', 'wrist exensor', ...
    'elbow flexor', 'elbow extensor', ...
    'deltoid flexor', 'deltoid extensor', ...
    'shoulder flexor', 'shoulder extensor'};
cols = 1;
rows = 2 + length(emg_pair); % only plot position and the muscles
axh = zeros(rows, cols);
dir_i   = 1;
subj_i = 1;
pert_i  = 1;
for fce_i = 2
    for dist_i = 2
        trials_list = obj.cond.trial;
        % make enough space for data
        dat.row = length(trials_list);
        dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
            obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
        dat.t = linspace(t_range(1), t_range(2), dat.col);
        dat.pos = nan(dat.row, dat.col);
        dat.fce = nan(dat.row, dat.col);
        dat.vel = nan(dat.row, dat.col);
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
            dat.vel(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ov(1,index_t),dat.t,'spline');
            for ch_i = 1:8
                try
                    dat.emg(ch_i,trial_idx,:)=interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,index_t)',dat.t,'spline')';
%                      dat.emg(ch_i,trial_idx,:)=interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emgrtf(ch_i,index_t)',dat.t,'spline')';
                catch
                    disp(['no EMG this condition! fce' num2str(fce_i) ' dist' num2str(dist_i)]);
                end
            end
        end

        % plot the mean and average of it


        axh(1, 1) = subplot('position', [0.1 0.7 0.8 0.16]); hold on;       % position
        plot(dat.t, dat.pos, ...
            'Color', col_type(4+dist_i,:), ...
            'LineWidth', 2);
        ylabel('x (m)'); title('position');

        %
        axh(2, 1) = subplot('position',  [0.1 0.5 0.8 0.16]); hold on;       % position
        plot(dat.t, dat.vel, ...
            'Color', col_type(4+dist_i,:), ...
            'LineWidth', 2);
        ylabel('v (m/s)'); title('velocity');

        %
        for muscle_i = 1:length(emg_pair)

            if mod(muscle_i,2) == 1
                axh(2+muscle_i, 1) = subplot('position', [0.1 0.3 0.8 0.16]); hold on;       % ch 3: elbow flexor
                ag_ant_i = 0;
            else
                axh(2+muscle_i, 1) = subplot('position', [0.1 0.1 0.8 0.16]); hold on;       % ch 3: elbow flexor
                ag_ant_i = 1;
            end

            muscle_idx = emg_pair(muscle_i);
            emg_tmp = reshape(dat.emg(muscle_idx,:,:), size(dat.emg,2), size(dat.emg,3))';
            plot(dat.t, emg_tmp, ... % muscle-trial-time
                 'Color', col_type(4+dist_i+ag_ant_i,:)); %;%, ...
            
            ylabel('EMG (mV)'); title(emg_pair_label{emg_pair(muscle_i)});

            lnh(1) = xline(0.05, 'LineStyle','--');
            lnh(2) = xline(0);
            grid on;

            
        end
            legend(lnh, '50ms', 'release');

        linkaxes(axh(:), 'x');
        % xlim([-3 2]); % sec
        set(axh(1, 1), 'YLim', [-0.50 -0.40 ])
        set(axh(2, 1), 'YLim', [-0.1 0.6]);
        xlim([-0.2 0.8]); % sec
        ylim_emg_list((emg_muscle_pairs-1)*2+1,:)=get(axh(3,1), 'YLim');
        ylim_emg_list((emg_muscle_pairs-1)*2+2,:)=get(axh(4,1), 'YLim');
    end
end
sgtitle(['subject' num2str(subj_i) ' EMG']);
end

%% 
% plot3, all the condition in 1 plot 
fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
Fs = 500;
col_type = colormap('lines');
t_range = [-0.5 1];

emg_pair = [1 2 3 4 5 6 7 8];
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
    'elbow flexor', 'elbow extensor', ...
    'anterior deltoid', 'posterior deltiod', ...
    'pectoralis', 'Trapezius'};
cols = 3;
rows = 2 + length(emg_pair); % only plot position and the muscles
axh = zeros(rows, cols);
%             subj_i  = 4;
subj_i  = 1;
dir_i   = 1;
pert_i  = 1;
for fce_i = 1:3
    for dist_i = 1:3
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
            dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,index_t),dat.t,'spline');
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
        axh(1, fce_i) = subplot(rows,cols,fce_i); hold on;       % position
        plot(dat.t, mean(dat.pos - dat.pos(:,1), 'omitnan'), ...
            'Color', col_type(4+dist_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
        tmp1 = mean(dat.pos - dat.pos(:,1), 'omitnan')+std(dat.pos, 'omitnan');
        tmp2 = mean(dat.pos - dat.pos(:,1), 'omitnan')-std(dat.pos, 'omitnan');
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type(4+dist_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        ylabel('m'); title('position');

        %
        axh(2, fce_i) = subplot(rows,cols,fce_i+3); hold on;       % force
        plot(dat.t, mean(dat.fce, 'omitnan'), ...
            'Color', col_type(4+dist_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
        tmp1 = mean(dat.fce, 'omitnan')+std(dat.fce,'omitnan');
        tmp2 = mean(dat.fce, 'omitnan')-std(dat.fce,'omitnan');
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type(4+dist_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        ylabel('N'); title('force');

        %
        for muscle_i = 1:length(emg_pair)
            axh(2+muscle_i, fce_i) = subplot(rows,cols,fce_i + (muscle_i+2-1)*3); hold on;       % ch 3: elbow flexor
            plot(dat.t, mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i,:), ...
                'LineWidth', 2);
            ylabel('EMG activity'); title(emg_pair_label{emg_pair(muscle_i)});

            patch_x = [dat.t, dat.t(end:-1:1)];
            tmp1 = mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
                std(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
            tmp2 = mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
                std(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
            patch_y = [tmp1, tmp2(end:-1:1)];
            patch(patch_x, patch_y, ...
                col_type(4+dist_i,:), ...
                'FaceAlpha', 0.3, ...
                'EdgeColor', 'none');

            xline(0.07);
            grid on;
        end
        %     end


        % xlim([-3 2]); % sec
        %set(axh(1, fce_i), 'YLim', [-0.50 -0.40 ])

        % ylim([-1 5]);

        %                     set(axh(3, fce_i), 'YLim', [0 0.1]);

    end
end
linkaxes(axh(:), 'x');
xlim([-0.5 0.8]); % sec
for muscle_i = 1:length(emg_pair)
    linkaxes(axh(2+muscle_i,:), 'y');
end
sgtitle(['Subject' num2str(subj_i) ' Direction' num2str(dir_i)]);

%% 
% plot4, grab data from behavioral landmark & do the barplot
%% Get the EMG averaged value based on the behavior landmarks 
% In this section, get the EMG averaged balue based on behavior landmarks.
% (velocity)
% EMG_fhd: Force_hold EMG that averaged 500ms before movement (which is 5%
% of the maximum velocity); 
% EMG_afl: release (at_fly) EMG that averaged from velocity cross 50% of
% its peak velocity to hit 50% of its peak velocity again. 
% EMG_phd: Position_hold EMG that averaged 200ms after the velocity
% below the 5% of its peak
%
%%%%%%%%%%%%%%%%%%% do something the same with the class function 

    %% generate data part
ifplot = 0; % the debug plot switch 
load([obj.data_dir '/' obj.data_name], 'data*');
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

dir_sig = [1 1 -1 -1];         % direction 3 and 4 have - pos and vel;
t_range = [-0.5 0.8];          % all the behavior landmark should come from this time. If not, it is task-unrelated movement.
time_window_fhold = [-0.5 0];  % When time aligned at beginging of movement, average EMG as the EMG_fhd
time_window_atfly= [nan nan];  % To find the peak velocity here, devined later
                               % When time aligned at the peak of velocity, average EMG as the EMG_afl
                               % not nessesary tobe the same across trials, as each velocity is different
time_window_phold = [0 0.1];   % When time aigned at the end of movement, defined by velocity <0.05v_max


Dat_EMG_fhd = nan(4,4,3,3,9,8);
Dat_EMG_afl = nan(4,4,3,3,9,8);
Dat_EMG_phd = nan(4,4,3,3,9,8);
fh = zeros(4,4,3);
mks_type = '*so';
Fs = 500;
fhtmp = figure();
col_type = colormap('lines');
close(fhtmp);


for subj_i  = 1:4    
    for dir_i   = 1:4
        pert_i  = 1;

        ifplotEMG_scatter = 0;

        Dat_EMG_pool = [];

        for fce_i = 1:3
            for dist_i = 1:3
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
                    % deal with exception from data collection problem
                    if (subj_i == 3 && dir_i == 2)
                        obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.v = ...
                            reshape(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ov(:,:,1), ...
                            size(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ov(:,:,1),[1 2]));
                    end
                    % get the index
                    trial_idx = trial_idx + 1;
                    % stack the data into matrices
                    index_t = find(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t > t_range(1) & ...
                        obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t < t_range(2));
                    v = smooth(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.v(1,index_t));
                    
                    t_ftrial = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t;
                    v_ftrial = smooth(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.v(1,:));
                    t = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t);
                    [v_max, idx_v_max] = max(dir_sig(dir_i) * v);
                    
                    % force-hold;
                    v_threshold1 = 0.05*v_max;
                    index_t_fhd = find(t > 0 & t < 0.1); % only could be after release
                    [v_closest1, idx_v_closestafr] = min(abs(v(index_t_fhd)-v_threshold1));
                    idx_v_closest1 = index_t_fhd(idx_v_closestafr);
                    t_fhd = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t - ...
                        obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t(idx_v_closest1));
                    fhd_idx = t_fhd > time_window_fhold(1) & t_fhd < time_window_fhold(2);

                    % at-fly
                    v_threshold2 = 0.5*v_max;
                    v_threshold3 = 0.5*v_max;
                    [v_closest2, idx_v_closest2] = min(abs(dir_sig(dir_i)*v(1:idx_v_max)-v_threshold2));
                    [v_closest3, idx_v_closest3] = min(abs(dir_sig(dir_i)*v(idx_v_max:end)-v_threshold2));
                    idx_v_closest3 = idx_v_max + idx_v_closest3-1;
                    t_afl = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t - ...
                        obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t(idx_v_max));
                    time_window_atfly = [t_afl(index_t(idx_v_closest2)), t_afl(index_t(idx_v_closest3))];
                    afl_idx = t_afl > time_window_atfly(1) & t_afl < time_window_atfly(2);
                        

                    % position-hold;
                    v_threshold4 = 0.05*v_max;
                    [v_closest4, idx_v_closest4] = min(abs(v(idx_v_max:end)-v_threshold4));
                    idx_v_closest4 = idx_v_max + idx_v_closest4-1;
                    t_phd = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t - ...
                        obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t(idx_v_closest4));
                    phd_idx = t_phd > time_window_phold(1) & t_phd < time_window_phold(2);
                
%                     ifplot = 0;
                    if (ifplot & 0) % check the velocity and its corresponding marks
                        clf; 

                        axh(1) = subplot(4,1,1); title('original time');
                        hold on; 
                        plot(t, v, 'linewidth', 2); 
                        yline(v_threshold1);
                        yline(v_threshold2);
                        plot(t(idx_v_closest1), v(idx_v_closest1), 'ro');
                        plot(t(idx_v_closest2), v(idx_v_closest2), 'go');
                        plot(t(idx_v_closest3), v(idx_v_closest3), 'bo');
                        plot(t(idx_v_closest4), v(idx_v_closest4), 'co');
                        xlabel('time (s)'); 
                        ylabel('velocity (m/s)');

                        axh(2) = subplot(4,1,2); title('force hold time');
                        hold on; 
                        plot(t_fhd, v_ftrial, 'linewidth', 2); 
                        yline(v_threshold1);
                        yline(v_threshold2);
                        xline(time_window_fhold);
                        xlabel('time (s)'); 
                        ylabel('velocity (m/s)');

                        axh(3) = subplot(4,1,3); title('peak velocity time');
                        hold on; 
                        plot(t_afl, v_ftrial, 'linewidth', 2); 
                        yline(v_threshold1);
                        yline(v_threshold2);
                        xline(time_window_atfly);
                        xlabel('time (s)'); 
                        ylabel('velocity (m/s)');

                        axh(4) = subplot(4,1,4); title('position hold time');
                        hold on; 
                        plot(t_phd, v_ftrial, 'linewidth', 2); 
                        yline(v_threshold1);
                        yline(v_threshold2);
                        xline(time_window_phold);
                        xlabel('time (s)'); 
                        ylabel('velocity (m/s)');
                        sgtitle('test if the velocity landmark alright');

                        linkaxes(axh, 'x'); xlim([-0.6 0.5]);
                    end

                    for ch_i = 1:8
                            emg_fhd = mean(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,fhd_idx), 'omitnan');
                            emg_afl = mean(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,afl_idx), 'omitnan');
                            emg_phd = mean(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,phd_idx), 'omitnan');
    
                            Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,trial_i,ch_i) = emg_fhd;
                            Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,trial_i,ch_i) = emg_afl;
                            Dat_EMG_phd(subj_i,dir_i,fce_i,dist_i,trial_i,ch_i) = emg_phd;
    
%                             ifplot = 0;
                            if (ifplot)
                                clf; clear axh;
                                axh(1) = subplot(2,1,1); 
                                hold on;
                                title('velocity'); 
                                plot(t_ftrial,v_ftrial); 
                                plot(t_ftrial(fhd_idx), v_ftrial(fhd_idx), 'r.');
                                plot(t_ftrial(afl_idx), v_ftrial(afl_idx), 'g.');
                                plot(t_ftrial(phd_idx), v_ftrial(phd_idx), 'b.');
                                xlabel('time (s)'); ylabel('v (m/s)');
                                legend('all', 'fhold', 'release', 'phold');
                                grid on;
    
                                axh(2) = subplot(2,1,2); 
                                hold on;
                                title(['EMG channel' num2str(ch_i)]); 
                                plot(t_ftrial,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,:)); 
                                plot(t_ftrial(fhd_idx), obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,fhd_idx), 'r.');
                                plot(t_ftrial(afl_idx), obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,afl_idx), 'g.');
                                plot(t_ftrial(phd_idx), obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,phd_idx), 'b.');
                                xlabel('time (s)'); ylabel('EMG portion');
                                grid on;
                        
                                linkaxes(axh, 'x');

                                xlim([-0.6 0.8]);
                                sgtitle('example of velocity and EMG');

                                figure('unit', 'inch', 'position', [0 0 2 3]); hold on;
                                val = [mean(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,fhd_idx)), ...
                                    mean(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,afl_idx)), ...
                                    mean(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,phd_idx))];
                                bar(1, val(1), 'r');
                                bar(2, val(2), 'g');
                                bar(3, val(3), 'b');
                                xticks([1:3]); xticklabels({'F', 'T', 'P'}); xlabel('epochs');
                                ylabel('EMG');
                                title('EMG during release');
                        end
                    end
                    
                end
            end
        end
    end
end


    %% barplot part
xx = [25 50 75];
x = [20 25 30; 45 50 55; 70 75 80];
ff = [15 20 25];
emg_pair_label = {'FCR', 'ECU', ...
    'BI', 'TRI', ...
    'AD', 'PD', ...
    'PEC', 'TPZ'};
epoch_key_idx = 6; 
epoch_key_names = {'force hold', 'transient', 'position hold', 'transient-forcehold',...
    'transient-forcehold\/forcehold', 'positionhold-forcehold'};
%               1 force_hold
%               2 transient (at-fly)
%               3 position_hold
%               4 transient - force_hold
%               5 (transient - force_hold)/force_hold
%               6 position_hold - force_hold
for subj_i = 1%:4%1:4 
    clear axh
    fh(subj_i) = figure('position', [0 0 600 1000]);
    for ch_i = 1:8
        for dir_i = 1:4
            axh(dir_i, ch_i) = subplot(8,4,(ch_i-1)*4+dir_i);
            emg_tmp_means = zeros(3,3);
            emg_tmp_stds = zeros(3,3);
            for fce_i = 1:3
                for dist_i = 1:3
                    switch epoch_key_idx
                        case 1                      % EMG_fhd
                            emg_tmp_means(fce_i,dist_i) = mean(Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                            emg_tmp_stds(fce_i,dist_i) = std(Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                        case 2                      % EMG_afl
                            emg_tmp_means(fce_i,dist_i) = mean(Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                            emg_tmp_stds(fce_i,dist_i) = std(Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                        case 3                      % EMG_phd
                            emg_tmp_means(fce_i,dist_i) = mean(Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                            emg_tmp_stds(fce_i,dist_i) = std(Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                        case 4                      % EMG_afl - EMG_fhd
                            emg_tmp_means(fce_i,dist_i) = mean((Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i) - Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i)), 'omitnan');
                            emg_tmp_stds(fce_i,dist_i)  =  std((Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i) - Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i)), 'omitnan');

                        case 5                     % (EMG_afl - EMG_fhd)/EMG_fhd
                            emg_tmp_means(fce_i,dist_i) = mean((Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i) - Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i))./Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                            emg_tmp_stds(fce_i,dist_i) = std((Dat_EMG_afl(subj_i,dir_i,fce_i,dist_i,:,ch_i) - Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i))./Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i), 'omitnan');
                        case 6                      % EMG_phd - EMG_fhd
                            emg_tmp_means(fce_i,dist_i) = mean((Dat_EMG_phd(subj_i,dir_i,fce_i,dist_i,:,ch_i) - Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i)), 'omitnan');
                            emg_tmp_stds(fce_i,dist_i)  =  std((Dat_EMG_phd(subj_i,dir_i,fce_i,dist_i,:,ch_i) - Dat_EMG_fhd(subj_i,dir_i,fce_i,dist_i,:,ch_i)), 'omitnan');
                    end
                end
            end
            lgd_fh = bar(xx, emg_tmp_means'); hold on;
            errorbar(x, emg_tmp_means', emg_tmp_stds', 'LineWidth',1, 'LineStyle', 'none', 'Color', 'k'); % is it right?
            
            grid on;


                switch ch_i 
                    case 1
                        title(['direction' num2str(dir_i)]);
                    case 8
                        xlabel('distance (mm)');
                end 

            if dir_i == 1
              ylabel({emg_pair_label{ch_i};'EMG' });
                if (ch_i == 1)
                    legend(lgd_fh, {'15N', '20N', '25N'});
                end
            end 
        end
    end
    for ch_i = 1:8
            linkaxes(axh(:,ch_i), 'y');
    end
%     sgtitle(['EMG activities of subject ' num2str(subj_i)]);
    sgtitle(fh(subj_i), 'String', {['EMG activities of subject ' num2str(subj_i)]; epoch_key_names{epoch_key_idx}});
end
