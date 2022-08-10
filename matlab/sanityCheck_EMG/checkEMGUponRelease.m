%% 1. plot data according to endpoint position/velocity 
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
% obj.data_name = 'ss4253_4263';
% obj.data_name = 'ss4253_4274';
% obj.data_name = 'ss4265_4274';
obj.data_name = 'ss4310_4341';

load([obj.data_dir '/' obj.data_name], 'data');
obj.data = data;
obj.cond.subj = 1:size(data,1);
obj.cond.dir  = 1:size(data,2);
obj.cond.fce  = 1:size(data,3);
obj.cond.dist = 1:size(data,4);
obj.cond.trial= 1:size(data,5);
obj.cond.pert = 1:size(data,6);

clear dat

% fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
Fs = 500;
fh1 = figure()
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
for subj_i  = 1:3%2;%2;
for dir_i   = 1:4
pert_i  = 1;
for fce_i = 2%1:3
    for dist_i = 2%1:3
        fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
        trials_list = obj.cond.trial;
        % make enough space for data
        dat.row = length(trials_list);
        dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
            obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
        %                     dat.t = t_range(1):1/Fs:t_range(2);
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
            ylabel('EMG (mV)'); title(emg_pair_label{emg_pair(muscle_i)});

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
            grid on;

            
        end
        legend(lnh, 'flexor', 'extensor');
        %     end

        linkaxes(axh(:), 'x');
        % xlim([-3 2]); % sec
%         set(axh(1, 1), 'YLim', [-0.50 -0.40 ])
        set(axh(1, 1), 'YLim', [-0.60 -0.40 ])
%         set(axh(2, 1), 'YLim', [-0.1 0.6]);
        set(axh(2, 1), 'YLim', [-0.6 0.6]);
        xlim([-0.5 0.5]); % sec
%         xlim([-0.5 2]); % sec
        % ylim([-1 5]);
        %                     linkaxes(axh(3:end,:), 'y');
        %                     set(axh(3, fce_i), 'YLim', [0 0.1]);
        sgtitle(['EMG demo subj' num2str(subj_i), 'dir' (dir_i)]);
    end
    
end

end
end
%% 1.2. plot data according to endpoint position/velocity + plot single raw EMG data
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
% obj.data_name = 'ss4253_4263';
% obj.data_name = 'ss4253_4274';
obj.data_name = 'ss4310_4341';
% 

load([obj.data_dir '/' obj.data_name], 'data');
obj.data = data;
obj.cond.subj = 1:size(data,1);
obj.cond.dir  = 1:size(data,2);
obj.cond.fce  = 1:size(data,3);
obj.cond.dist = 1:size(data,4);
obj.cond.trial= 1:size(data,5);
obj.cond.pert = 1:size(data,6);

clear dat
% subj_i  = 1; emgtmp = SessionScanEMG(4257);
% subj_i  = 2; emgtmp = SessionScanEMG(4269);
subj_i  = 3; emgtmp = SessionScanEMG(4336);

for emg_muscle_pairs = 1:4
fh(emg_muscle_pairs) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
Fs = 500;
col_type = colormap('lines');
t_range = [-0.5 2];

emg_pair = (emg_muscle_pairs-1)*2 + [1 2];
% emg_pair = [1 2 3 4 5 6 7 8];
% emg_pair = [1 2];
% emg_pair = [3 4];
% emg_pair = [5 6];
% emg_pair = [7 8];
emg_pair_label = {'wrist flexor', 'wrist exensor', ...
    'elbow flexor', 'elbow extensor', ...
    'deltoid flexor', 'deltoid extensor', ...
    'shoulder flexor', 'shoulder extensor'};
cols = 1;
rows = 2 + length(emg_pair); % only plot position and the muscles
axh = zeros(rows, cols);
dir_i   = 1;
pert_i  = 1;
for fce_i = 3%2
    for dist_i = 1%2
        trials_list = obj.cond.trial;
        % make enough space for data
        dat.row = length(trials_list);
        dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
            obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
        %                     dat.t = t_range(1):1/Fs:t_range(2);
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
%             dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,index_t),dat.t,'spline');
%             dat.vel(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.v(1,index_t),dat.t,'spline');
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
%         axh(1, 1) = subplot(rows,cols,1); hold on;       % position
        axh(1, 1) = subplot('position', [0.1 0.7 0.8 0.16]); hold on;       % position
%         plot(dat.t, dat.fce, ...
        plot(dat.t, dat.pos, ...
            'Color', col_type(4+dist_i,:), ...
            'LineWidth', 2);
%         patch_x = [dat.t, dat.t(end:-1:1)];
%         tmp1 = mean(dat.pos, 'omitnan')+std(dat.pos, 'omitnan');
%         tmp2 = mean(dat.pos, 'omitnan')-std(dat.pos, 'omitnan');
%         patch_y = [tmp1, tmp2(end:-1:1)];
%         patch(patch_x, patch_y, ...
%             col_type(4+dist_i,:), ...
%             'FaceAlpha', 0.3, ...
%             'EdgeColor', 'none');
        ylabel('x (m)'); title('position');

        %
%         axh(2, 1) = subplot(rows,cols,2); hold on;       % force
        axh(2, 1) = subplot('position',  [0.1 0.5 0.8 0.16]); hold on;       % position
        plot(dat.t, dat.vel, ...
            'Color', col_type(4+dist_i,:), ...
            'LineWidth', 2);
%         patch_x = [dat.t, dat.t(end:-1:1)];
%         tmp1 = mean(dat.vel, 'omitnan')+std(dat.vel,'omitnan');
%         tmp2 = mean(dat.vel, 'omitnan')-std(dat.vel,'omitnan');
%         patch_y = [tmp1, tmp2(end:-1:1)];
%         patch(patch_x, patch_y, ...
%             col_type(4+dist_i,:), ...
%             'FaceAlpha', 0.3, ...
%             'EdgeColor', 'none');
        ylabel('v (m/s)'); title('velocity');

        %
        for muscle_i = 1:length(emg_pair)

            if mod(muscle_i,2) == 1
%                 axh(2+muscle_i, 1) = subplot(rows,cols,2 + ceil(muscle_i)); hold on;       % ch 3: elbow flexor
                axh(2+muscle_i, 1) = subplot('position', [0.1 0.3 0.8 0.16]); hold on;       % ch 3: elbow flexor
                ag_ant_i = 0;
            else
%                 axh(2+muscle_i, 1) = subplot(rows,cols,2 + ceil(muscle_i)); hold on;       % ch 3: elbow flexor
                axh(2+muscle_i, 1) = subplot('position', [0.1 0.1 0.8 0.16]); hold on;       % ch 3: elbow flexor
                ag_ant_i = 1;
            end
%             lnh(ag_ant_i+1) = plot(dat.t, mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
%                 'Color', col_type(4+dist_i+ag_ant_i,:), ...
%                 'LineWidth', 2);
%             lnh(ag_ant_i+1) = 

            muscle_idx = emg_pair(muscle_i);
            emg_tmp = reshape(dat.emg(muscle_idx,:,:), size(dat.emg,2), size(dat.emg,3))';
            plot(dat.t, emg_tmp, ... % muscle-trial-time
                 'Color', col_type(4+dist_i+ag_ant_i,:)); %;%, ...
%                 'LineWidth', 2);
            
%             dat.emg_env = envelope(emg_tmp, 50, 'peak'); 
            ylabel('EMG (mV)'); title(emg_pair_label{emg_pair(muscle_i)});

%             patch_x = [dat.t, dat.t(end:-1:1)];
%             tmp1 = mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
%                 std(reshape(dat.emg(muscle_i,:,:), s 
% ize(dat.emg,2), size(dat.emg,3)),'omitnan');
%             tmp2 = mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
%                 std(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
%             patch_y = [tmp1, tmp2(end:-1:1)];
%             patch(patch_x, patch_y, ...
%                 col_type(4+dist_i+ag_ant_i,:), ...
%                 'FaceAlpha', 0.3, ...
%                 'EdgeColor', 'none');

            lnh(1) = xline(0.05, 'LineStyle','--');
            lnh(2) = xline(0);
            grid on;

            
        end
            legend(lnh, '50ms', 'release');
%         legend(lnh, 'flexor', 'extensor');
        %     end

        linkaxes(axh(:), 'x');
        % xlim([-3 2]); % sec
        set(axh(1, 1), 'YLim', [-0.50 -0.40 ])
        set(axh(2, 1), 'YLim', [-0.1 0.6]);
        xlim([-0.2 0.8]); % sec
        ylim_emg_list((emg_muscle_pairs-1)*2+1,:)=get(axh(3,1), 'YLim');
        ylim_emg_list((emg_muscle_pairs-1)*2+2,:)=get(axh(4,1), 'YLim');
%         xlim([-0.5 0.8]); % sec
        % ylim([-1 5]);
        %                     linkaxes(axh(3:end,:), 'y');
        %                     set(axh(3, fce_i), 'YLim', [0 0.1]);

    end
end
sgtitle(['subject' num2str(subj_i) ' EMG']);
end

% plot out the minimum activities as referencing 
% subject 1
figure(); 
clear t_range
t_range{1} = [1.0514 1.0515
           1.0514 1.0515
           1.0587, 1.0588
           1.0587, 1.0588
           1.05632, 1.05642
           1.05632, 1.05642
           1.0563, 1.0564
           1.0563, 1.0564]*1e4;
t_range{2} = [  1.4417 1.4418
                1.4417 1.4418
                1.4422 1.4423
                1.4531 1.4532
                1.4552 1.4553
                1.4737 1.4738
                1.44015 1.44025
                1.4737 1.4738]*1e4;
for i = 1:8 % muscles 
    if (mod(i,2)==1)
        figure(fh(ceil(i/2)));
        axes('position', [0.7 0.38 0.2 0.08]);
    else
        axes('position', [0.7, 0.18 0.2 0.08])
    end
    
    t_idx = find(emgtmp.data.t>t_range{subj_i}(i,1) & emgtmp.data.t<t_range{subj_i}(i,2));
    sum(t_idx)
    plot(emgtmp.data.t(t_idx) - emgtmp.data.t(t_idx(1)), ...
        emgtmp.data.emg(i,t_idx));
    grid on;
%     ylim([0 0.1])
    ylim(ylim_emg_list(i,:));
end

%% Plot the behavior against joint position/angular vilocity 
% just use the mean activities
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
obj.data_name = 'ss4253_4263';
% obj.data_name = 'ss4265_4274';
% obj.data_name = 'ss4253_4274';

load([obj.data_dir '/' obj.data_name], 'data');
obj.data = data;
obj.cond.subj = 1:size(data,1);
obj.cond.dir  = 1:size(data,2);
obj.cond.fce  = 1:size(data,3);
obj.cond.dist = 1:size(data,4);
obj.cond.trial= 1:size(data,5);
obj.cond.pert = 1:size(data,6);

fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
Fs = 500;
col_type = colormap('lines');
t_range = [-0.5 1];

emg_pair = [1 2 3 4 5 6 7 8];
emg_pair_label = {'wrist', 'wrist', ...
    'elbow', 'elbow', ...
    'anterior', 'posterior', ...
    'shoulder', 'shoulder'};
cols = 1;
rows = 2 + length(emg_pair)/2; % only plot position and the muscles
axh = zeros(rows, cols);
subj_i  = 1;%2;
dir_i   = 1;
pert_i  = 1;
for fce_i = 2
    for dist_i = 2
        trials_list = obj.cond.trial;
        % make enough space for data
        dat.row = length(trials_list);
        dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
            obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t < t_range(2));
        %                     dat.t = t_range(1):1/Fs:t_range(2);
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
            dat.pos_m1(:,trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(:,index_t,1)',dat.t,'spline')';
            dat.pos_m2(:,trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(:,index_t,2)',dat.t,'spline')';
            dat.pos_m3(:,trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(:,index_t,3)',dat.t,'spline')';
            dat.fce(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.f(1,index_t),dat.t,'spline');
            dat.vel(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ov(1,index_t),dat.t,'spline');
            %                         dat.emg(:,trial_idx,1:sum(index_t))=obj.data{subj_i,dir_i,fce_i,tar_i,trial_i,pert_i}.emg(:,index_t);
            for ch_i = 1:8
                try
                    dat.emg(ch_i,trial_idx,:)=interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,index_t)',dat.t,'spline')';
                catch
                    disp(['no EMG this condition! fce' num2str(fce_i) ' dist' num2str(dist_i)]);
                end
            end
        end
    end
end

%% take the average for all of them 
dat.t_ = dat.t;      % keep the same 
for xyzi = 1:3
dat.pos_m1_(xyzi,:) = mean(dat.pos_m1(xyzi,:,:), 'omitnan');
dat.pos_m2_(xyzi,:) = mean(dat.pos_m2(xyzi,:,:), 'omitnan');
dat.pos_m3_(xyzi,:) = mean(dat.pos_m3(xyzi,:,:), 'omitnan');
end
dat.fce_ = mean(dat.fce, 'omitnan');
dat.vel_ = mean(dat.vel, 'omitnan');
for ch_i = 1:8
    dat.emg_(ch_i,:) = mean(dat.emg(ch_i,:,:),'omitnan');
end
% calculate the two link vector 
index_t = dat.t>-0.5 & dat.t<0.5;
dat.vec_hum = dat.pos_m2_ - dat.pos_m3_;
dat.vec_rad = dat.pos_m1_ - dat.pos_m2_; 
[dat.vec_hum_theta,dat.vec_hum_len] = cart2pol(dat.vec_hum(1,:), dat.vec_hum(2,:));
[dat.vec_rad_theta, dat.vec_rad_len] = cart2pol(dat.vec_rad(1,:), dat.vec_rad(2,:));
dat.vec_rad_theta = dat.vec_rad_theta - dat.vec_hum_theta;
% calculate the radius velocity 
dat.jw1 = diff([0 dat.vec_hum_theta])./diff([dat.t(1) dat.t]);
dat.jw2 = diff([0 dat.vec_rad_theta])./diff([dat.t(1) dat.t]);

% sanity check on whether the length is the same 
ifplot = 1;
if (ifplot)
    fh = figure('name', 'Sanity Check: len & theta of humerus & radius');
    axh(1) = subplot(4,1,1); hold on;
    plot(dat.t(index_t), dat.vec_hum(1,index_t), '.');
    plot(dat.t(index_t), dat.vec_hum(2,index_t), '.');
    plot(dat.t(index_t), dat.vec_rad(1,index_t), '.');
    plot(dat.t(index_t), dat.vec_rad(2,index_t), '.');
    legend('h-x', 'h-y', 'r-x', 'r-y');
    axh(2) = subplot(4,1,2); hold on;
    plot(dat.t(index_t), dat.vec_hum_len(index_t));
    plot(dat.t(index_t), dat.vec_rad_len(index_t));
    legend('humerus', 'radius');
    ylabel('length (m)');
    axh(3) = subplot(4,1,3); hold on;
    plot(dat.t(index_t), dat.vec_hum_theta(index_t));
    plot(dat.t(index_t), dat.vec_rad_theta(index_t));
    legend('humerus', 'radius');
    ylabel('theta (rad)');
    axh(4) = subplot(4,1,4); hold on;
    plot(dat.t(index_t), dat.jw1(index_t));
    plot(dat.t(index_t), dat.jw2(index_t));
    legend('humerus', 'radius');
    ylabel('angular velocity (rad/s)');
end
linkaxes(axh, 'x');
xlim(axh(3), [-0.5 0.5]);


%% plot the joint angular position with the EMG activities

% 1. for the elbow joint, the bicep and tricep 
figure('name', 'elbow joint', 'unit', 'inch', 'position', [0 0 3 4]);
axh(1) = subplot(2,1,1); 
plot(dat.t(index_t), dat.vec_rad_theta(index_t), 'LineWidth', 3);
xlabel('time (s)');
ylabel('elbow joint angle (rad)');
title('joint angle');
axh(2) = subplot(2,1,2); hold on;
lnh(0+1) = plot(dat.t, mean(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+0,:), ...
                'LineWidth', 2);
lnh(1+1) = plot(dat.t, mean(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+1,:), ...
                'LineWidth', 2);
ylabel('EMG activity (mV)'); title('EMG for elbow joint');

patch_x = [dat.t, dat.t(end:-1:1)];
tmp1 = mean(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+0,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

tmp1 = mean(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+1,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

xlim([-0.5 0.5])
xline(0.07, 'LineStyle','--');
grid on;
legend(lnh, 'flexor', 'extensor');


% 2. for the shoulder joint velocity with the EMG activities 
figure('name', 'shoulder joint', 'unit', 'inch', 'position', [0 0 3 4]);
axh(1) = subplot(3,1,1); 
plot(dat.t(index_t), dat.vec_hum_theta(index_t), 'LineWidth', 3);
xlabel('time (s)');
ylabel('shoulder joint angle (rad)');
title('joint angle');
    % muscles: ant & post deltaloid
axh(2) = subplot(3,1,2); hold on;
lnh(0+1) = plot(dat.t, mean(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+0,:), ...
                'LineWidth', 2);
lnh(1+1) = plot(dat.t, mean(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+1,:), ...
                'LineWidth', 2);
ylabel('EMG activity (mV)'); title('EMG deltaloid ant/post');
patch_x = [dat.t, dat.t(end:-1:1)];
tmp1 = mean(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+0,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

tmp1 = mean(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+1,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');
xlim([-0.5 0.5])
xline(0.07);
grid on;
legend(lnh, 'flexor', 'extensor');
    % muscle: pec & trip
axh(2) = subplot(3,1,3); hold on;
lnh(0+1) = plot(dat.t, mean(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+0,:), ...
                'LineWidth', 2);
lnh(1+1) = plot(dat.t, mean(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+1,:), ...
                'LineWidth', 2);
ylabel('EMG activity (mV)'); title('EMG pectoralis/trapedius');
patch_x = [dat.t, dat.t(end:-1:1)];
tmp1 = mean(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+0,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

tmp1 = mean(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+1,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');


xlim([-0.5 0.5])
xline(0.07);
grid on;
legend(lnh, 'flexor', 'extensor');




%% plot the joint angular velocity with the EMG activities

% 1. for the elbow joint, the bicep and tricep 
figure('name', 'elbow joint', 'unit', 'inch', 'position', [0 0 3 4]);
axh(1) = subplot(2,1,1); 
plot(dat.t(index_t), dat.jw2(index_t), 'LineWidth', 3);
xlabel('time (s)');
ylabel('elbow joint angle velocity (rad/s)');
title('joint velocity');
axh(2) = subplot(2,1,2); hold on;
lnh(0+1) = plot(dat.t, mean(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+0,:), ...
                'LineWidth', 2);
lnh(1+1) = plot(dat.t, mean(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+1,:), ...
                'LineWidth', 2);
ylabel('EMG activity (mV)'); title('EMG for elbow joint');

patch_x = [dat.t, dat.t(end:-1:1)];
tmp1 = mean(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(3,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+0,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

tmp1 = mean(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(4,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+1,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

xlim([-0.5 0.5])
xline(0.07);
grid on;
legend(lnh, 'flexor', 'extensor');


% 2. for the shoulder joint velocity with the EMG activities 
figure('name', 'shoulder joint', 'unit', 'inch', 'position', [0 0 3 4]);
axh(1) = subplot(3,1,1); 
plot(dat.t(index_t), dat.jw1(index_t), 'LineWidth', 3);
xlabel('time (s)');
ylabel('shoulder joint angle velocity (rad/s)');
title('joint velocity');
    % muscles: ant & post deltaloid
axh(2) = subplot(3,1,2); hold on;
lnh(0+1) = plot(dat.t, mean(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+0,:), ...
                'LineWidth', 2);
lnh(1+1) = plot(dat.t, mean(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+1,:), ...
                'LineWidth', 2);
ylabel('EMG activity (mV)'); title('EMG deltaloid ant/post');
patch_x = [dat.t, dat.t(end:-1:1)];
tmp1 = mean(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(5,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+0,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

tmp1 = mean(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(6,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+1,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');
xlim([-0.5 0.5])
xline(0.07);
grid on;
legend(lnh, 'flexor', 'extensor');
    % muscle: pec & trip
axh(2) = subplot(3,1,3); hold on;
lnh(0+1) = plot(dat.t, mean(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+0,:), ...
                'LineWidth', 2);
lnh(1+1) = plot(dat.t, mean(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                'Color', col_type(4+dist_i+1,:), ...
                'LineWidth', 2);
ylabel('EMG activity (mV)'); title('EMG pectoralis/trapedius');
patch_x = [dat.t, dat.t(end:-1:1)];
tmp1 = mean(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(7,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+0,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

tmp1 = mean(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')+ ...
    std(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
tmp2 = mean(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan')- ...
    std(reshape(dat.emg(8,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan');
patch_y = [tmp1, tmp2(end:-1:1)];
patch(patch_x, patch_y, ...
    col_type(4+dist_i+1,:), ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');


xlim([-0.5 0.5])
xline(0.07);
grid on;
legend(lnh, 'flexor', 'extensor');




%%
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
            ylabel('EMG activity'); title(emg_pair_label{emg_pair(muscle_i)});

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

            xline(0.07);
            grid on;

            
        end
