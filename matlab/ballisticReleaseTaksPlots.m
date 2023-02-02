classdef ballisticReleaseTaksPlots
    %BALLISTICRELEASETAKSPLOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
%         data_name = 'ss4253_4263'
%         data_name = 'ss4253_4274'
%         data_name = 'ss4310_4341'
%         data_name = 'ss4310_4356'
%         data_name = 'ss4310_4314'
%         data_name = 'ss4351_4356';
%         data_name = 'ss4310_4356_emg'
%         data_name = 'ss4446_4467'       % other 3 subjects
%         data_name = 'ss4379_4414'       % chenguang test with optotrak markers on
%         data_name = 'ss4379_4438'       % 6-subject data
%         data_name = 'ss4472_4495'       % subject 10, 11, and 12
%         data_name = 'ss4491_4495'       % subject 12
%         data_name = 'ss4500_4524'       % subject 13, 14, 15 
%         data_name = 'ss4512_4524'       % subject 14, and 15
%         data_name = 'ss4530_4563'       % subject 16, 17, 18
          data_name = 'ss4573_4587'       % subject 19, 20
%         data_name = 'ss4385_4388'       % subjeH data
%         data_name = 'ss4379_4383'       % subjeC data
%         data_name = 'ss4401_4427'       % subjeB data
%         data_name = 'ss4408_4414'       % subjeM data
%         data_name = 'ss4418_4422'       % subjeM data
%         data_name = 'ss4379_4422'       % subjeM data
        data
        data_idx_ss
        data_idx_tr
        cond
    end
    
    methods
        function obj = ballisticReleaseTaksPlots()
            %BALLISTICRELEASETAKSPLOTS Construct an instance of this class
            %   read data to construct the plot
            load([obj.data_dir '/' obj.data_name], 'data*');
            obj.data = data;
            cond.subj = 1:size(data,1);
            cond.dir  = 1:size(data,2);
            cond.fce  = 1:size(data,3);
            cond.dist = 1:size(data,4);
            cond.trial= 1:size(data,5);
            cond.pert = 1:size(data,6);

            obj.cond = cond;

            obj.data_idx_ss = data_index_ss;
            obj.data_idx_tr = data_index_tr;
        end
        
        function fh = plotx_samef_release(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of release
            % use the mean and std
            %   Detailed explanation goes here
            for subj_i  = 1:4
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
            
            dir_i   = 1;
            pert_i  = 1;
            for fce_i = 1:3
                for dist_i = 1:3
                    trials_list = obj.cond.trial;
                    % make enough space for data
                    dat.row = length(trials_list);
                    dat.col = sum(obj.data{subj_i,dir_i,fce_i,dist_i,1,pert_i}.t > t_range(1) & ...
                        obj.data{subj_i,1,fce_i,dist_i,1,1}.t < t_range(2));
                    % dat.t = t_range(1):1/Fs:t_range(2);
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
                        % dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,index_t),dat.t,'spline');
                        dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,index_t),dat.t,'spline');
                        dat.fce(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.f(1,index_t),dat.t,'spline');
                        % dat.emg(:,trial_idx,1:sum(index_t)) = obj.data{subj_i,dir_i,fce_i,tar_i,trial_i,pert_i}.emg(:,index_t);
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
        end

        function fh = plotx_raw_subj_dir(obj)

            Data = obj.data(:,:,:,:,:,:);
            sbj_idx = 4;%[1:6];%[1 2];%[1:7];  %%[1 4 5 6];
%             subj_names = subj_names_all(sbj_idx);
%             Data = data(sbj_idx,:,:,:,:,:);
            fce_levels = [15 20 25];
            dist_levels = [2.5 5.0 7.5];
            clear axh
            Freq = 2000;
            t_step = 1/Freq;
            fh = figure();
            colors = colormap('lines');
            close(fh);
            r = size(Data, 1); %r = 2;% subj
            c = size(Data, 2); % direction
            f = size(Data, 3); % force
            d = size(Data, 4); % target
            l = size(Data, 5); % trials
            p = size(Data, 6); % perturbation type
            idx_lead = Freq/5;
            idx_last = Freq/2;
            if_subtract = 0;
            epoc_type = 2;
            plot_type = 8;%4;  
            % 1 displacement
            % 2 force
            % 3 Fp
            % 4 vlocity
            ifavg = false;
            %
            if (~ifavg)
                for si = 1:r % subj
                    for ci = 1:c % direction
                        figure(); 
                        for fi = 1:f%1:2 % target force
                            % axh(ri, ci) = subplot(f,c,f*(ci-1) + ci);
                            % axh(1, ci) = subplot(1,f,ci); grid on; hold on; % plot on columns
                            % axh(ci, 1) = subplot(f,1,ci); grid on; hold on; % plot on rows
                            % axh(ci, si) = subplot(f,r,(r*(ci-1)+si)); grid on; hold on; % plot on rows
                            
                            if (ci==1)
%                                 title(subj_names{si});
%                                 title(['subj#' num2str(si)]);
                            end
                            for di = 1:d % target distance
                                axh(di, fi) = subplot(d,f,f*(di-1) + fi); grid on; hold on;
                                %                   axh(di) = subplot(d,1,di); grid on; hold on;
                                %                 axh(di) = subplot(1,1,1); grid on; hold on;
                                for li = 1%:p % perturbation
                                    trial_num = length(Data(si,ci,fi,di,:,li));
                                    for ti = 1:trial_num % each trial
                                        if (isempty(Data{si,ci,fi,di,ti,li}))
                                            continue;
                                        end
                                        switch epoc_type
                                            case 1
                                                idx = find(Data{si,ci,fi,di,ti,li}.Fp(2,:)~=0);
                                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                            case 2
                                                idx = find(Data{si,ci,fi,di,ti,li}.ts==5 | Data{si,ci,fi,di,ti,li}.ts==6);
                                                idx = (idx-idx_lead):idx(end);
                                                %                                 idx = (idx(1)):idx(end);
                                        end
                                        time = t_step*(idx-idx(idx_lead));
%                                             time = t_step*(idx);
                                        %                                         switch plot_type
                                        %                                             case 2
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.f(2,idx);
                                        %                                                 titlestr = 'force';
                                        %                                             case 1
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.x(2,idx);
                                        %                                                 dat = dat - dat(1);
                                        %                                                 titlestr = 'displacement';
                                        %                                             case 3
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.Fp(2,idx);
                                        %                                                 titlestr = 'Fp';
                                        %                                             case 4
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.v(2,idx);
                                        %                                                 titlestr = 'velocity';
                                        %                                             case 5
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.tq(4,idx);
                                        %                                                 titlestr = 'torque4';
                                        %                                             case 6
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.x(:,idx);
                                        %                                                 dat_submean = dat - mean(dat(:,1:50),2);
                                        %                                                 dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                        %                                                 dat = dat_norm .* sign(dat_submean(2,:));
                                        %                                                 titlestr = 'norm displacement';
                                        %                                             case 7 % the force mode
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.f(:,idx);
                                        %                                                 dat_submean = dat;% - mean(dat(:,1:50),2);
                                        %                                                 dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                        %                                                 dat = dat_norm .* sign(dat_submean(2,:));
                                        %                                                 titlestr = 'norm force';
                                        %
                                        %                                         end

                                        switch plot_type
                                            case 1
                                                dat = Data{si,ci,fi,di,ti,li}.x(1,idx);
                                                titlestr = 'displacement';
                                            case 2
                                                dat = Data{si,ci,fi,di,ti,li}.f(xyi,idx);
                                                titlestr = 'force';
                                            case 3
                                                dat = Data{si,ci,fi,di,ti,li}.Fp(xyi,idx);
                                                titlestr = 'Fp';
                                            case 4
                                                dat = Data{si,ci,fi,di,ti,li}.v(xyi,idx);
                                                titlestr = 'velocity';
                                            case 5
                                                dat = Data{si,ci,fi,di,ti,li}.tq(3,idx);
                                                titlestr = 'torque3';
                                            case 6
                                                dat = Data{si,ci,fi,di,ti,li}.x(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm displacement';
                                            case 7 % the force mode
                                                dat = Data{si,ci,fi,di,ti,li}.f(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm force';
                                            case 8 % optotrak x
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ox))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(1,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(1,idx,1);    % marker 1, x
%                                                             dat = Data{si,ci,fi,di,ti,li}.ox(3,idx,1);    % marker 1, z
%                                                         dat = Data{si,ci,fi,di,ti,li}.ox(1,idx,2); % marker 2,x
%                                                         dat = Data{si,ci,fi,di,ti,li}.ox(1,idx,3); % marker 3,x
                                                end
                                                titlestr = 'opto-position';
                                            case 9 % optotrak v
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ov))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx,1);
                                                end
                                                titlestr = 'opto-velocity';
                                            case 10
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,1);
                                                titlestr = 'opto1-position';
                                            case 11
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,2);
                                                titlestr = 'opto2-position';
                                            case 12
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,3);
                                                titlestr = 'opto3-position';
                                        end




                                        if (if_subtract)
                                            dat = dat - mean(dat(1:50));
                                        end

%                                         plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                                        plot(time, dat, 'Color', colors(4*(2-1)+di, :));
                                        %                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                                        xlabel('t (s)');
                                        ylabel('x (m)');
                                    end
                                end
                                title(['fce' num2str(fce_levels(fi)) 'dist' num2str(dist_levels(di))], ...
                                    'FontSize', 10);
                            end
                        end
                        sgtitle(['subj#' num2str(si) ' direction#' num2str(ci)]);
                        linkaxes(axh,'xy');
                        xlim([-0.2 0.8]);
%                         xlim([0.0 1.0]);
                        switch ci
                            case 1
                            case 2
%                                 ylim([-0.5 -0.38]);
                            case 3
                            case 4
%                                 ylim([-0.57 -0.45]);
                        end
                    end
                end

            elseif (ifavg)

                for si = 1:r % subj
                    for ci = 1:c % direction
                        for fi = 1:f%1:2 % target force
                            %axh(ri, ci) = subplot(f,c,f*(ci-1) + ci);
                            %             axh(1, ci) = subplot(1,f,ci); grid on; hold on; % plot on columns
                            %            axh(ci, 1) = subplot(f,1,ci); grid on; hold on; % plot on rows
                            axh(ci, si) = subplot(f,r,(r*(ci-1)+si)); grid on; hold on; % plot on rows
                            if (ci==1)
                                title(subj_names{si});
                            end
                            for di = 1:d % target distance
                                %                   axh(di) = subplot(d,1,di); grid on; hold on;
                                %                 axh(di) = subplot(1,1,1); grid on; hold on;
                                for li = 2%1:p % perturbation
                                    trial_num = length(Data(si,ci,ci,di,:,li));
                                    dat_mat = [];
                                    dat_avg = [];
                                    dat_std = [];
                                    for ti = 1:trial_num % each trial
                                        if (isempty(Data{si,ci,fi,di,ti,li}))
                                            continue;
                                        end
                                        switch epoc_type
                                            case 1
                                                idx = find(Data{si,ci,fi,di,ti,li}.Fp(2,:)~=0);
                                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                            case 2
                                                idx = find(Data{si,ci,fi,di,ti,li}.ts==5 | Data{si,ci,fi,di,ti,li}.ts==6);
                                                if (length(idx)<300)
                                                    continue
                                                else
                                                    idx = (idx-50):idx(300);
                                                end
                                                % limit idx length
                                        end
                                        time = t_step*(idx-idx(1));
                                        switch plot_type
                                            case 1
                                                dat = Data{si,ci,fi,di,ti,li}.x(xyi,idx);
                                                titlestr = 'displacement';
                                            case 2
                                                dat = Data{si,ci,fi,di,ti,li}.f(xyi,idx);
                                                titlestr = 'force';
                                            case 3
                                                dat = Data{si,ci,fi,di,ti,li}.Fp(xyi,idx);
                                                titlestr = 'Fp';
                                            case 4
                                                dat = Data{si,ci,fi,di,ti,li}.v(xyi,idx);
                                                titlestr = 'velocity';
                                            case 5
                                                dat = Data{si,ci,fi,di,ti,li}.tq(3,idx);
                                                titlestr = 'torque3';
                                            case 6
                                                dat = Data{si,ci,fi,di,ti,li}.x(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm displacement';
                                            case 7 % the force mode
                                                dat = Data{si,ci,fi,di,ti,li}.f(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm force';
                                            case 8 % optotrak x
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ox))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,1);
                                                end
                                                titlestr = 'opto-position';
                                            case 9 % optotrak v
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ov))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx,1);
                                                end
                                                titlestr = 'opto-velocity';
                                            case 10
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,1);
                                                titlestr = 'opto1-position';
                                            case 11
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,2);
                                                titlestr = 'opto2-position';
                                            case 12
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,3);
                                                titlestr = 'opto3-position';
                                        end

                                        if (if_subtract)
                                            dat = dat - mean(dat(1:50));
                                        end

                                        dat_mat = [dat_mat; dat];
                                        %plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                                        %                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                                        %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
                                    end
                                    dat_avg = nanmean(dat_mat);
                                    dat_std = nanstd(dat_mat);
                                    % plot the center
                                    plot(time, dat_avg, 'Color', colors(4*(li-1)+di, :), 'LineWidth', 3);
                                    patch_x = [time time(end:-1:1)];
                                    patch_y = [dat_avg+dat_std dat_avg(end:-1:1)-dat_std(end:-1:1)];
                                    patch('XData', patch_x, 'YData', patch_y, ...
                                        'FaceColor', colors(4*(li-1)+di, :), ...
                                        'FaceAlpha', 0.3);

                                    %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
                                    %                    %plot(figure(11), time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :), '--');
                                end
                            end
                        end
                    end
                end

            end

        end
        function fh = plotf_raw_subj_dir(obj)

            Data = obj.data(:,:,:,:,:,:);
            sbj_idx = 1;%[1:6];%[1 2];%[1:7];  %%[1 4 5 6];
%             subj_names = subj_names_all(sbj_idx);
%             Data = data(sbj_idx,:,:,:,:,:);
            fce_levels = [15 20 25];
            dist_levels = [2.5 5.0 7.5];
            clear axh
            Freq = 2000;
            t_step = 1/Freq;
            fh = figure();
            colors = colormap('lines');
            close(fh);
            r = size(Data, 1); %r = 2;% subj
            c = size(Data, 2); % direction
            f = size(Data, 3); % force
            d = size(Data, 4); % target
            l = size(Data, 5); % trials
            p = size(Data, 6); % perturbation type
            idx_last = Freq/2;
            idx_lead = Freq/5;
            if_subtract = 0;
            epoc_type = 2;
            plot_type = 2;%4;  % 1 displacement
            % 2 force
            % 3 Fp
            % 4 vlocity
            ifavg = false;
            %
            if (~ifavg)
                for si = 1:r % subj
                    for ci = 1:c % direction
                        figure(); 
                        for fi = 1:f%1:2 % target force
                            % axh(ri, ci) = subplot(f,c,f*(ci-1) + ci);
                            % axh(1, ci) = subplot(1,f,ci); grid on; hold on; % plot on columns
                            % axh(ci, 1) = subplot(f,1,ci); grid on; hold on; % plot on rows
                            % axh(ci, si) = subplot(f,r,(r*(ci-1)+si)); grid on; hold on; % plot on rows
                            
                            if (ci==1)
%                                 title(subj_names{si});
%                                title(['subj#' num2str(si)]);
                            end
                            for di = 1:d % target distance
                                axh(di, fi) = subplot(d,f,f*(di-1) + fi); grid on; hold on;
                                %                   axh(di) = subplot(d,1,di); grid on; hold on;
                                %                 axh(di) = subplot(1,1,1); grid on; hold on;
                                for li = 1%:p % perturbation
                                    trial_num = length(Data(si,ci,fi,di,:,li));
                                    for ti = 1:trial_num % each trial
                                        if (isempty(Data{si,ci,fi,di,ti,li}))
                                            continue;
                                        end
                                        switch epoc_type
                                            case 1
                                                idx = find(Data{si,ci,fi,di,ti,li}.Fp(2,:)~=0);
                                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                            case 2
                                                idx = find(Data{si,ci,fi,di,ti,li}.ts==5 | Data{si,ci,fi,di,ti,li}.ts==6);
                                                idx = (idx-idx_lead):idx(end);
                                                %                                 idx = (idx(1)):idx(end);
                                        end
                                        time = t_step*(idx-idx(idx_lead));
%                                             time = t_step*(idx);
                                        %                                         switch plot_type
                                        %                                             case 2
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.f(2,idx);
                                        %                                                 titlestr = 'force';
                                        %                                             case 1
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.x(2,idx);
                                        %                                                 dat = dat - dat(1);
                                        %                                                 titlestr = 'displacement';
                                        %                                             case 3
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.Fp(2,idx);
                                        %                                                 titlestr = 'Fp';
                                        %                                             case 4
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.v(2,idx);
                                        %                                                 titlestr = 'velocity';
                                        %                                             case 5
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.tq(4,idx);
                                        %                                                 titlestr = 'torque4';
                                        %                                             case 6
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.x(:,idx);
                                        %                                                 dat_submean = dat - mean(dat(:,1:50),2);
                                        %                                                 dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                        %                                                 dat = dat_norm .* sign(dat_submean(2,:));
                                        %                                                 titlestr = 'norm displacement';
                                        %                                             case 7 % the force mode
                                        %                                                 dat = Data{si,ri,ci,di,ti,li}.f(:,idx);
                                        %                                                 dat_submean = dat;% - mean(dat(:,1:50),2);
                                        %                                                 dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                        %                                                 dat = dat_norm .* sign(dat_submean(2,:));
                                        %                                                 titlestr = 'norm force';
                                        %
                                        %                                         end

                                        switch plot_type
                                            case 1
                                                dat = Data{si,ci,fi,di,ti,li}.x(xyi,idx);
                                                titlestr = 'displacement';
                                            case 2
%                                                 dat = Data{si,ci,fi,di,ti,li}.f(xyi,idx);
                                                dat = Data{si,ci,fi,di,ti,li}.f(1,idx);
                                                titlestr = 'force';
                                            case 3
                                                dat = Data{si,ci,fi,di,ti,li}.Fp(xyi,idx);
                                                titlestr = 'Fp';
                                            case 4
                                                dat = Data{si,ci,fi,di,ti,li}.v(xyi,idx);
                                                titlestr = 'velocity';
                                            case 5
                                                dat = Data{si,ci,fi,di,ti,li}.tq(3,idx);
                                                titlestr = 'torque3';
                                            case 6
                                                dat = Data{si,ci,fi,di,ti,li}.x(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm displacement';
                                            case 7 % the force mode
                                                dat = Data{si,ci,fi,di,ti,li}.f(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm force';
                                            case 8 % optotrak x
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ox))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(1,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(1,idx,1);
                                                end
                                                titlestr = 'opto-position';
                                            case 9 % optotrak v
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ov))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx,1);
                                                end
                                                titlestr = 'opto-velocity';
                                            case 10
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,1);
                                                titlestr = 'opto1-position';
                                            case 11
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,2);
                                                titlestr = 'opto2-position';
                                            case 12
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,3);
                                                titlestr = 'opto3-position';
                                        end




                                        if (if_subtract)
                                            dat = dat - mean(dat(1:50));
                                        end

%                                         plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                                        plot(time, dat, 'Color', colors(4*(2-1)+di, :));
                                        %                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                                    end
                                end
                                title(['fce' num2str(fce_levels(fi)) 'dist' num2str(dist_levels(di))], ...
                                    'FontSize', 10);
                                xlabel('t (s)');
                                ylabel('F (N)');
                            end
                        end
                        sgtitle(['subj#' num2str(si) ' direction#' num2str(ci)]);
                        linkaxes(axh,'xy');
                        xlim([-0.2 0.8]);
%                         xlim([0.0 1.0]);
                        switch ci
                            case 1
                            case 2
%                                 ylim([-0.5 -0.38]);
                            case 3
                            case 4
%                                 ylim([-0.57 -0.45]);
                        end
                    end
                end

            elseif (ifavg)

                for si = 1:r % subj
                    for ci = 1:c % direction
                        for fi = 1:f%1:2 % target force
                            %axh(ri, ci) = subplot(f,c,f*(ci-1) + ci);
                            %             axh(1, ci) = subplot(1,f,ci); grid on; hold on; % plot on columns
                            %            axh(ci, 1) = subplot(f,1,ci); grid on; hold on; % plot on rows
                            axh(ci, si) = subplot(f,r,(r*(ci-1)+si)); grid on; hold on; % plot on rows
                            if (ci==1)
                                title(subj_names{si});
                            end
                            for di = 1:d % target distance
                                %                   axh(di) = subplot(d,1,di); grid on; hold on;
                                %                 axh(di) = subplot(1,1,1); grid on; hold on;
                                for li = 2%1:p % perturbation
                                    trial_num = length(Data(si,ci,ci,di,:,li));
                                    dat_mat = [];
                                    dat_avg = [];
                                    dat_std = [];
                                    for ti = 1:trial_num % each trial
                                        if (isempty(Data{si,ci,fi,di,ti,li}))
                                            continue;
                                        end
                                        switch epoc_type
                                            case 1
                                                idx = find(Data{si,ci,fi,di,ti,li}.Fp(2,:)~=0);
                                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                            case 2
                                                idx = find(Data{si,ci,fi,di,ti,li}.ts==5 | Data{si,ci,fi,di,ti,li}.ts==6);
                                                if (length(idx)<300)
                                                    continue
                                                else
                                                    idx = (idx-50):idx(300);
                                                end
                                                % limit idx length
                                        end
                                        time = t_step*(idx-idx(1));
                                        switch plot_type
                                            case 1
                                                dat = Data{si,ci,fi,di,ti,li}.x(xyi,idx);
                                                titlestr = 'displacement';
                                            case 2
                                                dat = Data{si,ci,fi,di,ti,li}.f(xyi,idx);
                                                titlestr = 'force';
                                            case 3
                                                dat = Data{si,ci,fi,di,ti,li}.Fp(xyi,idx);
                                                titlestr = 'Fp';
                                            case 4
                                                dat = Data{si,ci,fi,di,ti,li}.v(xyi,idx);
                                                titlestr = 'velocity';
                                            case 5
                                                dat = Data{si,ci,fi,di,ti,li}.tq(3,idx);
                                                titlestr = 'torque3';
                                            case 6
                                                dat = Data{si,ci,fi,di,ti,li}.x(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm displacement';
                                            case 7 % the force mode
                                                dat = Data{si,ci,fi,di,ti,li}.f(:,idx);
                                                dat_submean = dat - mean(dat(:,1:50),2);
                                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                                titlestr = 'norm force';
                                            case 8 % optotrak x
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ox))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,1);
                                                end
                                                titlestr = 'opto-position';
                                            case 9 % optotrak v
                                                switch length(size(Data{si,ci,fi,di,ti,li}.ov))
                                                    case 2
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx);
                                                    case 3
                                                        dat = Data{si,ci,fi,di,ti,li}.ov(xyi,idx,1);
                                                end
                                                titlestr = 'opto-velocity';
                                            case 10
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,1);
                                                titlestr = 'opto1-position';
                                            case 11
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,2);
                                                titlestr = 'opto2-position';
                                            case 12
                                                dat = Data{si,ci,fi,di,ti,li}.ox(xyi,idx,3);
                                                titlestr = 'opto3-position';
                                        end

                                        if (if_subtract)
                                            dat = dat - mean(dat(1:50));
                                        end

                                        dat_mat = [dat_mat; dat];
                                        %plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                                        %                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                                        %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
                                    end
                                    dat_avg = nanmean(dat_mat);
                                    dat_std = nanstd(dat_mat);
                                    % plot the center
                                    plot(time, dat_avg, 'Color', colors(4*(li-1)+di, :), 'LineWidth', 3);
                                    
                                    patch_x = [time time(end:-1:1)];
                                    patch_y = [dat_avg+dat_std dat_avg(end:-1:1)-dat_std(end:-1:1)];
                                    patch('XData', patch_x, 'YData', patch_y, ...
                                        'FaceColor', colors(4*(li-1)+di, :), ...
                                        'FaceAlpha', 0.3);

                                    %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
                                    %                    %plot(figure(11), time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :), '--');
                                end
                            end
                        end
                    end
                end

            end

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
%             subj_i  = 4;
            subj_i  = 4;
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

                    % set(axh(3, fce_i), 'YLim', [0 0.1]);

                end
            end
            linkaxes(axh(:), 'x');
            xlim([-0.5 0.8]); % sec
            for muscle_i = 1:length(emg_pair)
                linkaxes(axh(2+muscle_i,:), 'y');
            end
            sgtitle(['Subject' num2str(subj_i) ' Direction' num2str(dir_i)]);

        end
    
        function fh = plotEMG_release_dir1_acrosssubject(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of release
            % use the mean and std
            %   Detailed explanation goes here
            
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
            for subj_i  = 1:4
                fh(subj_i) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
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
            end
        end
    
        function fh = plotEMG_release_acrossdir(obj, subj_i)
            %PLOTEMG_RELEASE plot out the emg data at the time of release
            % use the mean and std
            %   Detailed explanation goes here

            if (~exist('subj_i', 'var'))
                subj_i = 1;
            end
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
%             subj_i  = 1:4
            for dir_i = 1:4
                fh(subj_i) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 6 12]);
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
                            try
                                dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,index_t),dat.t,'spline');
                            catch
                                display(['trial', num2str(trial_idx)]);
                                dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,index_t),dat.t,'spline');
                            end
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
            end
        end

        function fh = plotEMG_release_1(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of release
            % use the mean and std
            %   Detailed explanation goes here
            fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
            Fs = 500;
            col_type = colormap('lines');
            close(fh); clear fh;
            t_range = [-0.5 1];

            emg_pair = [1 2 3 4 5 6 7 8];
            emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                'elbow flexor', 'elbow extensor', ...
                'anterior deltaloid', 'posterior deltaloid', ...
                'pectoralis', 'Trapezius'};
            cols = 3;
            rows = 2 + length(emg_pair); % only plot position and the muscles 
            axh = zeros(rows, cols);
            for subj_i  = 1:4 %;
                fh(subj_i) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
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
            sgtitle(['variate force subj'  num2str(subj_i) 'dir' num2str(dir_i)]);
            end
        end
    
        function fh = plotEMG_release_specificCond(obj, subj_i, dir_i, fce_i, tar_i, fh)

            if ~exist('subj_i', 'var')
                subj_i = 1;
            end
            if ~exist('dir_i', 'var')
                dir_i = 1;
            end
            if ~exist('fce_i', 'var')
                fce_i = 1;
            end
            if ~exist('tar_i', 'var')
                tar_i = 1;
            end

            if (~exist('fh', 'var'))
                for emg_muscle_pairs = 1:4
                    fh(emg_muscle_pairs) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
                end
            else
                for emg_muscle_pairs = 1:4
                    set(fh(emg_muscle_pairs), 'name', 'EMG');
                    set(fh(emg_muscle_pairs), 'unit', 'inch');
                    set(fh(emg_muscle_pairs), 'position', [0 0 8.5 11]);
                end
            end

            
            

            for emg_muscle_pairs = 1:4
                Fs = 500;
                col_type = colormap('lines');
%                 t_range = [-0.5 2];
                t_range = [-5 3];
%                 figure(fh(emg_muscle_pairs)); 
                set(0, 'CurrentFigure', fh(emg_muscle_pairs));
                emg_pair = (emg_muscle_pairs-1)*2 + [1 2];
                emg_pair_label = {'wrist flexor', 'wrist exensor', ...
                    'elbow flexor', 'elbow extensor', ...
                    'deltoid flexor', 'deltoid extensor', ...
                    'shoulder flexor', 'shoulder extensor'};
                cols = 1;
                rows = 2 + length(emg_pair); % only plot position and the muscles
                axh = zeros(rows, cols);
                pert_i  = 1;
                    dist_i = tar_i;
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
                        legend_str = cell(1,length(trials_list));
                        for trial_i = trials_list
                            % get the index
                            trial_idx = trial_idx + 1;
                            legend_str{trial_i} = ['ss' num2str(obj.data_idx_ss(subj_i,dir_i,fce_i,dist_i,trial_i,pert_i)) ...
                                'tr' num2str(obj.data_idx_tr(subj_i,dir_i,fce_i,dist_i,trial_i,pert_i))];
                            % stack the data into matrices
                            index_t = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t > t_range(1) & ...
                                obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t < t_range(2);
                            dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,index_t),dat.t,'spline');
                            dat.fce(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.f(1,index_t),dat.t,'spline');
                            dat.vel(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ov(1,index_t),dat.t,'spline');
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
                        
                        axh(1, 1) = subplot('position', [0.1 0.7 0.8 0.16]); hold on;       % position
                        for trial_i = trials_list
%                             col_tmp = 1-(1-col_type(4+dist_i,:))/length(trials_list)*trial_i;
                            col_tmp = col_type(trial_i,:);
                            plot(dat.t, dat.pos(trial_i,:), ...
                                'Color', col_tmp, ...
                                'LineWidth', 2);
                        end
%                         plot(dat.t, dat.pos, ...
%                                 'Color', col_type(4+dist_i,:)/length(trials_list)*trial_i, ...
%                                 'LineWidth', 2);
                        ylabel('x (m)'); title('position');


                        axh(2, 1) = subplot('position',  [0.1 0.5 0.8 0.16]); hold on;       % velocity
                        for trial_i = trials_list
%                             col_tmp = 1-(1-col_type(4+dist_i,:))/length(trials_list)*trial_i;
                            col_tmp = col_type(trial_i,:);
                            plot(dat.t, dat.vel(trial_i,:), ...
                                'Color', col_tmp, ...
                                'LineWidth', 2);
                        end
%                         plot(dat.t, dat.vel, ...
%                                 'Color', col_type(4+dist_i,:)/length(trials_list)*trial_i, ...
%                                 'LineWidth', 2);
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
                            
                              for trial_i = trials_list
%                                   col_tmp = 1-(1-col_type(4+dist_i+ag_ant_i,:))/length(trials_list)*trial_i;
                                  col_tmp = col_type(trial_i,:);
                                  lnh_emg(trial_i) = plot(dat.t, emg_tmp(:,trial_i), ... % muscle-trial-time
                                      'Color', col_tmp, ...
                                      'LineWidth', 1); %;%, ...
                              end
                          

%                             plot(dat.t, emg_tmp, ... % muscle-trial-time
%                                 'Color', col_type(4+dist_i+ag_ant_i,:), ...
%                                 'LineWidth', 1); %;%, ...
                            

                            ylabel('EMG (% of MVF)'); 
                            title(emg_pair_label{emg_pair(muscle_i)});

                            lnh(1) = xline(0.05, 'LineStyle','--');
                            lnh(2) = xline(0);
                            grid on;

                            legend([lnh lnh_emg] , {'50ms', 'release', legend_str{:}});
                        end
%                         legend(lnh , '50ms', 'release');
                          
                        
                        linkaxes(axh(:), 'x');
                        switch (dir_i)
                            case 1 
                            case 2
                                set(axh(1, 1), 'YLim', [-0.50 -0.36])
                                set(axh(2, 1), 'YLim', [-0.1 0.6]);
                            case 3
                                set(axh(1, 1), 'YLim', [-0.60 -0.46])
                                set(axh(2, 1), 'YLim', [-0.6 0.1]);
                            case 4

                        end
%                         xlim([-0.2 0.8]); % sec
%                         xlim([-0.5 1.0]); % sec
                        xlim([-5.0 3.0]); % sec
                        ylim_emg_list((emg_muscle_pairs-1)*2+1,:)=get(axh(3,1), 'YLim');
                        ylim_emg_list((emg_muscle_pairs-1)*2+2,:)=get(axh(4,1), 'YLim');
                sgtitle(['EMG of subject' num2str(subj_i) 'dir' num2str(dir_i) ...
                    'fce' num2str(fce_i) 'dist' num2str(dist_i)]);
            end
        end

        function fh = plotEMG_release_specificCond_fce(obj, subj_i, dir_i, fce_i, tar_i, fh)

            if ~exist('subj_i', 'var')
                subj_i = 1;
            end
            if ~exist('dir_i', 'var')
                dir_i = 1;
            end
            if ~exist('fce_i', 'var')
                fce_i = 1;
            end
            if ~exist('tar_i', 'var')
                tar_i = 1;
            end

            if (~exist('fh', 'var'))
                for emg_muscle_pairs = 1:4
                    fh(emg_muscle_pairs) = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
                end
            else
                for emg_muscle_pairs = 1:4
                    set(fh(emg_muscle_pairs), 'name', 'EMG');
                    set(fh(emg_muscle_pairs), 'unit', 'inch');
                    set(fh(emg_muscle_pairs), 'position', [0 0 8.5 11]);
                end
            end

            
            

            for emg_muscle_pairs = 1:4
                Fs = 500;
                col_type = colormap('lines');
%                 t_range = [-0.5 2];
                t_range = [-5 3];
%                 figure(fh(emg_muscle_pairs)); 
                set(0, 'CurrentFigure', fh(emg_muscle_pairs));
                emg_pair = (emg_muscle_pairs-1)*2 + [1 2];
                emg_pair_label = {'wrist flexor', 'wrist exensor', ...
                    'elbow flexor', 'elbow extensor', ...
                    'deltoid flexor', 'deltoid extensor', ...
                    'shoulder flexor', 'shoulder extensor'};
                cols = 1;
                rows = 2 + length(emg_pair); % only plot position and the muscles
                axh = zeros(rows, cols);
                pert_i  = 1;
                    dist_i = tar_i;
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
                        legend_str = cell(1,length(trials_list));
                        for trial_i = trials_list
                            % get the index
                            trial_idx = trial_idx + 1;
                            legend_str{trial_i} = ['ss' num2str(obj.data_idx_ss(subj_i,dir_i,fce_i,dist_i,trial_i,pert_i)) ...
                                'tr' num2str(obj.data_idx_tr(subj_i,dir_i,fce_i,dist_i,trial_i,pert_i))];
                            % stack the data into matrices
                            index_t = obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t > t_range(1) & ...
                                obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t < t_range(2);
                            dat.pos(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,index_t),dat.t,'spline');
                            dat.fce(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.f(1,index_t),dat.t,'spline');
                            dat.vel(trial_idx,:) = interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ov(1,index_t),dat.t,'spline');
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
                        
                        axh(1, 1) = subplot('position', [0.1 0.7 0.8 0.16]); hold on;       % position
                        for trial_i = trials_list
%                             col_tmp = 1-(1-col_type(4+dist_i,:))/length(trials_list)*trial_i;
                            col_tmp = col_type(trial_i,:);
                            plot(dat.t, dat.pos(trial_i,:), ...
                                'Color', col_tmp, ...
                                'LineWidth', 2);
                        end
%                         plot(dat.t, dat.pos, ...
%                                 'Color', col_type(4+dist_i,:)/length(trials_list)*trial_i, ...
%                                 'LineWidth', 2);
                        ylabel('x (m)'); title('position');


                        axh(2, 1) = subplot('position',  [0.1 0.5 0.8 0.16]); hold on;       % velocity
                        for trial_i = trials_list
%                             col_tmp = 1-(1-col_type(4+dist_i,:))/length(trials_list)*trial_i;
                            col_tmp = col_type(trial_i,:);
                            plot(dat.t, dat.fce(trial_i,:), ...
                                'Color', col_tmp, ...
                                'LineWidth', 2);
                        end
%                         plot(dat.t, dat.vel, ...
%                                 'Color', col_type(4+dist_i,:)/length(trials_list)*trial_i, ...
%                                 'LineWidth', 2);
                        ylabel('F (N)'); title('Force');

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
                            
                              for trial_i = trials_list
%                                   col_tmp = 1-(1-col_type(4+dist_i+ag_ant_i,:))/length(trials_list)*trial_i;
                                  col_tmp = col_type(trial_i,:);
                                  lnh_emg(trial_i) = plot(dat.t, emg_tmp(:,trial_i), ... % muscle-trial-time
                                      'Color', col_tmp, ...
                                      'LineWidth', 1); %;%, ...
                              end
                          

%                             plot(dat.t, emg_tmp, ... % muscle-trial-time
%                                 'Color', col_type(4+dist_i+ag_ant_i,:), ...
%                                 'LineWidth', 1); %;%, ...
                            

                            ylabel('EMG (portion of MVF)'); 
                            title(emg_pair_label{emg_pair(muscle_i)});

                            lnh(1) = xline(0.05, 'LineStyle','--');
                            lnh(2) = xline(0);
                            grid on;

                            legend([lnh lnh_emg] , {'50ms', 'release', legend_str{:}});
                        end
%                         legend(lnh , '50ms', 'release');
                          
                        
                        linkaxes(axh(:), 'x');
                        switch (dir_i)
                            case 1 
                            case 2
                                set(axh(1, 1), 'YLim', [-0.50 -0.36])
                                set(axh(2, 1), 'YLim', [-0.1 0.6]);
                            case 3
                                set(axh(1, 1), 'YLim', [-0.60 -0.46])
                                set(axh(2, 1), 'YLim', [-0.6 0.1]);
                            case 4

                        end
%                         xlim([-0.2 0.8]); % sec
%                         xlim([-0.5 1.0]); % sec
                        xlim([-1.5 2.0]); % sec
                        ylim_emg_list((emg_muscle_pairs-1)*2+1,:)=get(axh(3,1), 'YLim');
                        ylim_emg_list((emg_muscle_pairs-1)*2+2,:)=get(axh(4,1), 'YLim');
                sgtitle(['EMG of subject' num2str(subj_i) 'dir' num2str(dir_i) ...
                    'fce' num2str(fce_i) 'dist' num2str(dist_i)]);
            end
        end


        function fh = plotEMG_release_specificCond_Shade(obj, subj_i, dir_i, fce_i, tar_i, fh)

            if ~exist('subj_i', 'var')
                subj_i = 1;
            end
            if ~exist('dir_i', 'var')
                dir_i = 1;
            end
            if ~exist('fce_i', 'var')
                fce_i = 1;
            end
            if ~exist('tar_i', 'var')
                tar_i = 1;
            end


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


            pert_i = 1;
            dist_i = tar_i;

            if (~exist('fh', 'var'))
                fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
            else
                fh = figure(fh); 
                set(fh, 'name', 'EMG');
                set(fh, 'unit', 'inch');
                set(fh, 'position', [0 0 8.5 11]);
            end
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
                %                         dat.emg(:,trial_idx,1:sum(index_t))=obj.data{subj_i,dir_i,fce_i,tar_i,trial_i,pert_i}.emg(:,index_t);
                for ch_i = 1:8
                    try
                        dat.emg(ch_i,trial_idx,:)=interp1(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t(index_t),obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,index_t)',dat.t,'spline')';
                    catch
                        disp(['no EMG this condition! fce' num2str(fce_i) ' dist' num2str(dist_i)]);
                    end
                end
            end


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
                ylabel('EMG %');
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
                grid on;


            end
            legend(lnh, 'flexor', 'extensor');
            %     end

            linkaxes(axh(:), 'x');
            set(axh(1, 1), 'YLim', [-0.60 -0.40 ])
            set(axh(2, 1), 'YLim', [-0.6 0.6]);
            xlim([-0.5 0.5]); % sec
            sgtitle(['EMG demo subj' num2str(subj_i), 'dir' num2str(dir_i) 'f' num2str(fce_i) 'd' num2str(dist_i)]);
        end


        function fh = plotEMG_comparepulseEffect(obj)
            %PLOTEMG_RELEASE plot out the emg data at the time of perturb
            % use the mean and std
            %   Detailed explanation goes here
%             fh = figure('name', 'EMG', 'unit', 'inch', 'position', [0 0 7 12]);
            Freq = 2000;
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
            t_step = 1/Freq;
            colors = colormap('lines');
               
            Data = obj.data;

            r = size(Data, 1); % subj
            d = size(Data, 2); % direction
            f = size(Data, 3); % force
            l = size(Data, 4); % distance
            t = size(Data, 5); % trials
            p = size(Data, 6); % perturbation type
            idx_last = Freq/2;
            idx_lead = Freq/5;
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
                                        idx = (idx(1)-idx_lead):idx(end);
                                        %                                 idx = (idx(1)-500):(idx(end)+idx_lead);
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
                                                idx = (idx(1)-idx_lead):idx(end);
                                                %                                 idx = (idx(1)-500):(idx(end)+idx_lead);
                                                %idx = (idx(1)):idx(end);
                                        end
                                        %plot(Data{ri,ci,di,ti,li}.Fp(xyi,:));
                                        %idx = find(Data{ri,ci,di,ti,li}.Fp(xyi,:)~=0);
                                        %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                        %idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                                        %idx = (idx-9):idx(end);
                                        %idx = (idx(1)):(idx(end)+idx_lead);
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
                        
                        %         xlim([0.2 1.2]);
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

