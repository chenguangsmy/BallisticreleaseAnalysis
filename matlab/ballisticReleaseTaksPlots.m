classdef ballisticReleaseTaksPlots
    %BALLISTICRELEASETAKSPLOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
%         data_name = 'ss4253_4263'
        data_name = 'ss4253_4274'
        data
        cond
    end
    
    methods
        function obj = ballisticReleaseTaksPlots()
            %BALLISTICRELEASETAKSPLOTS Construct an instance of this class
            %   read data to construct the plot 
            load([obj.data_dir '/' obj.data_name], 'data');
            obj.data = data;
            cond.subj = 1:size(data,1);
            cond.dir  = 1:size(data,2);
            cond.fce  = 1:size(data,3);
            cond.dist = 1:size(data,4);
            cond.trial= 1:size(data,5);
            cond.pert = 1:size(data,6);

            obj.cond = cond;
        end
        
        function fh = plotEMG_release(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of release
            % use the mean and std
            %   Detailed explanation goes here
            fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
            Fs = 500;
            col_type = colormap('lines');
            t_range = [-0.5 1];

            emg_pair = [1 2 3 4 5 6 7 8];
            emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                'elbow flexor', 'elbow extensor', ...
                'anterior deltaloid', 'posterior deltaloid', ...
                'pectoralis', 'Trapezius'};
            cols = 3;
            rows = 2 + length(emg_pair); % only plot position and the muscles 
            axh = zeros(rows, cols);
            subj_i  = 1;
            dir_i   = 1;
            pert_i  = 1;
            for fce_i = 1:3
                for dist_i = 1:3
                    trials_list = obj.cond.trial;
                    % make enough space for data
                    dat.row = length(trials_list);
                    dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
                        obj.data{1,1,fce_i,dist_i,1,1}.t < t_range(2));
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

                    linkaxes(axh(:), 'x');
                    % xlim([-3 2]); % sec
                    set(axh(1, fce_i), 'YLim', [-0.50 -0.40 ])
                    xlim([-0.5 0.8]); % sec
                    % ylim([-1 5]);
%                     linkaxes(axh(3:end,:), 'y');
%                     set(axh(3, fce_i), 'YLim', [0 0.1]);

                end
            end
            sgtitle('variate force');

        end
    
        function fh = plotEMG_release_1(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of release
            % use the mean and std
            %   Detailed explanation goes here
            fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
            Fs = 500;
            col_type = colormap('lines');
            t_range = [-0.5 1];

            emg_pair = [1 2 3 4 5 6 7 8];
            emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                'elbow flexor', 'elbow extensor', ...
                'anterior deltaloid', 'posterior deltaloid', ...
                'pectoralis', 'Trapezius'};
            cols = 3;
            rows = 2 + length(emg_pair); % only plot position and the muscles 
            axh = zeros(rows, cols);
            subj_i  = 2;
            dir_i   = 1;
            pert_i  = 1;
            for fce_i = 1:3
                for dist_i = 1:3
                    trials_list = obj.cond.trial;
                    % make enough space for data
                    dat.row = length(trials_list);
                    dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
                        obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t  < t_range(2));
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
                    axh(1, fce_i) = subplot(rows,cols,fce_i); hold on;       % position
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
                    axh(2, fce_i) = subplot(rows,cols,fce_i+3); hold on;       % force
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

                    linkaxes(axh(:), 'x');
                    % xlim([-3 2]); % sec
                    set(axh(1, fce_i), 'YLim', [-0.50 -0.40 ])
                    set(axh(2, fce_i), 'YLim', [-0.1 0.6]);
                    xlim([-0.5 0.8]); % sec
                    % ylim([-1 5]);
%                     linkaxes(axh(3:end,:), 'y');
%                     set(axh(3, fce_i), 'YLim', [0 0.1]);

                end
            end
            sgtitle('variate force');

        end
    
    end
end

