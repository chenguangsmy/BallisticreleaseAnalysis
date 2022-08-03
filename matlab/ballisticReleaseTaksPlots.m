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
        
        function fh = plotx_samef_release(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of release
            % use the mean and std
            %   Detailed explanation goes here
            fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
            Fs = 500;
            col_type = colormap('lines');
            t_range = [-0.1 0.5];

            emg_pair = [1 2 3 4 5 6 7 8];
            emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                'elbow flexor', 'elbow extensor', ...
                'anterior deltoid', 'posterior deltiod', ...
                'pectoralis', 'Trapezius'};
            cols = 1;
            rows = 3; % only plot position and the muscles 
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
                    ylabel('x (m)'); title('position');
                    grid on;
                    %  
                    
                end
                linkaxes(axh, 'x');
            end
            sgtitle('variate force');
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
                'anterior deltoid', 'posterior deltiod', ...
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

                    linkaxes(axh(:), 'x');
                    % xlim([-3 2]); % sec
                    %set(axh(1, fce_i), 'YLim', [-0.50 -0.40 ])
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
            fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
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
                    ylabel('x (m)'); title('position');
                    
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
                    ylabel('v (m/s)'); title('velocity');

                    % 
                    for muscle_i = 1:length(emg_pair)
                        axh(2+muscle_i, fce_i) = subplot(rows,cols,fce_i + (muscle_i+2-1)*3); hold on;       % ch 3: elbow flexor
                        plot(dat.t, mean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)),'omitnan'), ...
                            'Color', col_type(4+dist_i,:), ...
                            'LineWidth', 2);
                        ylabel('EMG (mV)'); title(emg_pair_label{emg_pair(muscle_i)});

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
                    set(axh(1, fce_i), 'YLim', [-0.50 -0.40 ])
                    set(axh(2, fce_i), 'YLim', [-0.1 0.6]);
                    xlim([-0.5 0.8]); % sec
                    % ylim([-1 5]);
%                     linkaxes(axh(3:end,:), 'y');
%                     set(axh(3, fce_i), 'YLim', [0 0.1]);

                end
            end

            for muscle_i = 1:length(emg_pair)
                 linkaxes(axh(2+muscle_i,:), 'xy');
            end
            sgtitle('variate force');

        end
    
        function fh = plotEMG_comparepulseEffect(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of perturb
            % use the mean and std
            %   Detailed explanation goes here
%             fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
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
            t_step = 1/500;
            colors = colormap('lines');
               
            Data = obj.data;

            r = size(Data, 1); % subj
            d = size(Data, 2); % direction
            f = size(Data, 3); % force
            l = size(Data, 4); % distance
            t = size(Data, 5); % trials
            p = size(Data, 6); % perturbation type
            idx_last = 200;
            if_subtract = 0;
            fce_list = [15 20 25];
            dist_list = [2.5 5.0 7.5]; % cm

            epoc_type = 1;  % 1 perturb
            pert_type = 2:3; % choose [2 3]
            axh = zeros(d,r);
            xyi = 1;        % x, y

            % fh = figure();
            % pert_type = 2; colors = colors(4:end,:)
            for ri = 1:r % subj
                fh = figure('Name', ['subject' num2str(ri)]);
                for di = 1%:d % direction
                    %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
                    %         axh(ri, di) = subplot(d,r,r*(di-1) + ri);grid on;hold on;
                    %         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
                    % first, plot the force pulse
                    for fi = 1:f % target force
                        for li = 1:3 % target distance
                            for pi = pert_type%1:p % perturbation
                                switch epoc_type
                                    case 1
                                        %                                 idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0 & Data{ri,di,fi,li,ti,pi}.ts==4);  % pert at y
                                        idx = find(Data{ri,di,fi,li,1,pi}.Fp(xyi,:)~=0);
                                        idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                        if pi == 1
                                            disp('ERROR: should use pi == 2!!!');
                                        end
                                    case 2
                                        idx = find(Data{ri,di,fi,li,1,pi}.ts==5 | Data{ri,di,fi,li,ti,pi}.ts==6);
                                        idx = (idx(1)-100):idx(end);
                                        %                                 idx = (idx(1)-500):(idx(end)+100);
                                        %idx = (idx(1)):idx(end);
                                end
                                



%                                 dat = Data{ri,1,fi,li,1,pi}.Fp(1,idx);
                                dat = Data{ri,1,fi,li,1,pi}.f(1,idx);
                                time = t_step*(idx-idx(1));
                                axh(1,fi) = subplot(8+1,f,fi); grid on; hold on;
                                plot(time, dat, 'Color', [0.5 0.5 0.5]);
                                title(['fce' num2str(fce_list(fi))]);


                                for mi = 1:8 % muscles

                                    axh(1+mi, fi) = subplot(8+1,f,f*(mi) + fi);grid on;hold on; %?



                                    trial_num = length(Data(ri,di,fi,li,:,pi));
                                    clear Dat;
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
                                                %                                 idx = (idx(1)-500):(idx(end)+100);
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

                                        channeli = mi;
                                        dat = Data{ri,di,fi,li,ti,pi}.emg(channeli,idx);
                                        Dat(ti,:) = dat;
%                                         titlestr = ['EMG' num2str(channeli)];

                                        if (if_subtract)
                                            dat = dat - mean(dat(1:50));
                                        end
                                        %                                     plot(time, dat, 'Color', colors(4+li, :));

                                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                                    end
                                    
                                    plot(time, mean(Dat, 'omitnan'), 'Color', colors(4+li, :));

                                    patch_x = [time, time(end:-1:1)];
                                    tmp1 = mean(Dat,'omitnan')+ std(Dat,'omitnan');
                                    tmp2 = mean(Dat,'omitnan')- std(Dat,'omitnan');
                                    patch_y = [tmp1, tmp2(end:-1:1)];
                                    patch(patch_x, patch_y, ...
                                        col_type(4+li,:), ...
                                        'FaceAlpha', 0.3, ...
                                        'EdgeColor', 'none');

                                    
                                    

                                    
                                end
                                
                                
                                xticks([0.2 0.45 0.7 1.0]);
                                xticklabels({'0' '0.25' '0.5' '0.8'});
                            end
                        end
                        
                        %         xlim([0.2 1.2])
                        %         ylim([0 0.1]);
                        %         ylim([0 0.1]);
                    end
%                     linkaxes(axh(2:9,:), 'xy');
%                     ylim([0 0.1])
%                     xlim([0 1.0])

                    sgtitle('EMG after pulse');
                end
            % xlim([-2 2])
            % xlim([0 0.5]);
            % xlim([0 2])
            
            end




        end
    
% ss
    end
end

