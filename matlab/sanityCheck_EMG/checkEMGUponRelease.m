%% 1. plot data according to endpoint position/velocity 
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
obj.data_name = 'ss4253_4263';
% obj.data_name = 'ss4265_4274';

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
subj_i  = 1;
dir_i   = 1;
pert_i  = 1;
for fce_i = 2
    for dist_i = 2
        trials_list = obj.cond.trial;
        % make enough space for data
        dat.row = length(trials_list);
        dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
            obj.data{1,1,fce_i,dist_i,1,1}.t < t_range(2));
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
        ylabel('m'); title('position');

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
        ylabel('N'); title('force');

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
        legend(lnh, 'flexor', 'extensor');
        %     end

        linkaxes(axh(:), 'x');
        % xlim([-3 2]); % sec
        set(axh(1, 1), 'YLim', [-0.50 -0.40 ])
        set(axh(2, 1), 'YLim', [-0.1 0.6]);
        xlim([-0.5 0.5]); % sec
        % ylim([-1 5]);
        %                     linkaxes(axh(3:end,:), 'y');
        %                     set(axh(3, fce_i), 'YLim', [0 0.1]);

    end
end
sgtitle('variate force');


%% Plot the behavior against joint position/angular vilocity 
% just use the mean activities
obj.data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
% obj.data_name = 'ss4253_4263';
% obj.data_name = 'ss4265_4274';
obj.data_name = 'ss4253_4274';

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
subj_i  = 1;
dir_i   = 1;
pert_i  = 1;
for fce_i = 2
    for dist_i = 2
        trials_list = obj.cond.trial;
        % make enough space for data
        dat.row = length(trials_list);
        dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
            obj.data{1,1,fce_i,dist_i,1,1}.t < t_range(2));
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
ylabel('EMG activity (mV)'); title('EMG pectorilod/traledeous');
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
