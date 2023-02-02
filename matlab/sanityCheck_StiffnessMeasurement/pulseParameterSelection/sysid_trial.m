% sysid_trial
% Chenguang write system identification code dealing with the data he get. 
% Trying Federico's approach to see if I can get some close-related result.


% 1. Get data 
fce_list = [15 20 25];
dist_list = [0.025 0.050 0.075];
load('pulseParameterCompare', 'data');  % just release, Chenguang

K_est_mat = zeros(3,3,3,7); % subj-fce-disp_subj
for subj_i  = 1:3     % thinnest pulse
    dir_i   = 1;
    for fce_i   = 1:3
        for disp_i  = 1:3
            for trial_i = 1:7

                origin_freq     = 500;

                resample_freq   = 500;

                resample_time   = [-0.1 0.5];
                resample_tgrid  = 0:1/resample_freq:range(resample_time);


                % find pulse
                fce_cmd = abs(data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.Fp(1,:));
                [~, idx_Fpeak] = max(fce_cmd);
                idx_sel = origin_freq*min(resample_time) : ...
                    1 : ...
                    origin_freq*max(resample_time);
                idx_sel = idx_sel+idx_Fpeak;

                % resample
                t_preitp = data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.t(idx_sel) - ...
                    data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.t(idx_sel(1));
                f_preitp = data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.f(1,idx_sel);
                fp_preitp = data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.Fp(1,idx_sel);
                x_preitp = data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.ox(1,idx_sel);

                %%%%%%%%%%%%de-trend the x_preitp
                detrend.x = 2;
                if (detrend.x)
                    switch detrend.x
                        case 1 % linear
                            x_preitp_trend = (x_preitp(end) - x_preitp(1))/length(x_preitp) ... % slope
                                * [1:length(x_preitp)];
                        case 2 % piece-wise
                            % Get the velocity threshold
                            vel_threshold = 0.01;
                            v_preitp = diff(x_preitp)./diff(t_preitp);
                            v_preitp_idx = find(abs(v_preitp)>vel_threshold);
                            v_preitp_idx_length = v_preitp_idx(end) - v_preitp_idx(1) + 1;
                            x_preitp_trend = x_preitp(1)*ones(size(x_preitp));
                            % piece 1, before velocity increase
                            x_preitp_trend(1:(v_preitp_idx(1)-1)) = x_preitp(v_preitp_idx(1));
                            % piece 2, at the middle
                            x_preitp_trend(v_preitp_idx(1):v_preitp_idx(end)) = ...
                                (x_preitp(v_preitp_idx(end)) - x_preitp(v_preitp_idx(1))) / (v_preitp_idx_length) ... % slope
                                * (1:v_preitp_idx_length) ...
                                + x_preitp(1);
                            % piece 3, at the end
                            x_preitp_trend(v_preitp_idx(end):end) = x_preitp(v_preitp_idx(end));
                            if (0)
                                clf; hold on;
                                plot(x_preitp, 'r.');
                                plot(x_preitp_trend, 'b.');
                            end
                    end
                    x_preitp_detrend = x_preitp - x_preitp_trend;
                    x_preitp_raw = x_preitp;
                    x_preitp = x_preitp_detrend; % use the detrend x
                end
                %%%%%%%%%%%%

                %%%%%%%%%%%%de-trend the f_preitp
                detrend.f = 1;
                if (detrend.f)
                    f_sf_threshold = 1; %0.5; %StaticFriction
                    f_command_validx = abs(fp_preitp) > f_sf_threshold;
                    fp_preitp_filter = f_preitp;
                    fp_preitp_filter(f_command_validx) = nan;
                    fp_preitp_offset = f_preitp - mean(fp_preitp_filter,'omitnan'); % remove offset
                    fp_preitp_offset(abs(fp_preitp_offset)<f_sf_threshold) = 0;
                    fp_preitp_filted = fp_preitp_offset;
                    fp_preitp_filted(~f_command_validx & fp_preitp_offset>=f_sf_threshold) = ...
                        fp_preitp_offset(~f_command_validx & fp_preitp_offset>=f_sf_threshold) - f_sf_threshold;
                    fp_preitp_filted(~f_command_validx & fp_preitp_offset<=-f_sf_threshold) = ...
                        fp_preitp_offset(~f_command_validx & fp_preitp_offset<=-f_sf_threshold) + f_sf_threshold;

                    f_preitp_raw = f_preitp;
                    f_preitp = fp_preitp_filted;

                    ifplot = 0;
                    if(ifplot)
                        clf;
                        hold on;
                        plot(f_preitp, 'b.');
                        plot(fp_preitp_offset, 'g.');
                        plot(fp_preitp_filted, 'r.');
                        legend('origin', 'offset', 'filted');
                    end
                end
                %%%%%%%%%%%%

                f_preitp_raw = -(f_preitp_raw - f_preitp_raw(1));
                f_preitp = -(f_preitp - f_preitp(1));
                x_preitp = x_preitp - x_preitp(1);
                x_preitp_raw = x_preitp_raw - x_preitp_raw(1);
                t_aftitp = resample_tgrid;
                f_aftitp = interp1(t_preitp,f_preitp,resample_tgrid,'spline');
                x_aftitp = interp1(t_preitp,x_preitp,resample_tgrid,'spline');
                % remove-out the offset




                % plot the original data;
                ifplot = 0;
                if (ifplot)
                    figure();
                    clf;
                    subplot(1,2,1);
                    plot(data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.t(idx_sel), ...
                        data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.f(1,idx_sel));
                    subplot(1,2,2);
                    plot(data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.t(idx_sel), ...
                        data{subj_i,dir_i,fce_i,disp_i,trial_i,2}.ox(1,idx_sel));
                end
                % plot the time-resampled data;

                ifplot = 1;
                clear lnh
                if (ifplot)
                    fh = figure();
                    clf;
                    subplot(1,2,1); hold on;
                    plot(t_preitp, f_preitp_raw, 'b-.');
                    plot(t_preitp, f_preitp, 'k-.');
                    plot(t_aftitp, f_aftitp, 'r.')
                    subplot(1,2,2); hold on;
                    plot(t_preitp, x_preitp_raw, 'b-.');
                    lnh(1) = plot(t_preitp, x_preitp, 'k-.');
                    lnh(2) = plot(t_aftitp, x_aftitp, 'r.');
                end


                % 2. do system id

                data_est = iddata(x_aftitp',f_aftitp',1/resample_freq);

                sys_pulse = tfest(data_est,2,0);
                [NUM_P,DEN_P] = tfdata(sys_pulse);
                K_est = DEN_P{1}(3)/NUM_P{1}(3);
                B_est = DEN_P{1}(2)/NUM_P{1}(3);
                M_est = DEN_P{1}(1)/NUM_P{1}(3);
                FIT = sys_pulse.Report.Fit.FitPercent;

                est.G = tf(NUM_P, DEN_P);
                est.U = f_aftitp';
                est.t = t_aftitp;%time_up{i}(idxmin:idxmax);
                est.y = lsim(est.G, est.U, est.t);

                if (ifplot)
                    figure(fh);
                    subplot(1,2,1); title('force');
                    subplot(1,2,2); title('displacement');
                    lnh(3) = plot(est.t, est.y, '.g');
                    sgtitle({['f' num2str(fce_i) 'd' num2str(disp_i) 't' num2str(trial_i)];...
                        ['Kest' num2str(K_est)]});

                    legend(lnh, 'raw', 'pre-processed', 'fitting');
                end

                K_est_mat(subj_i, fce_i, disp_i, trial_i) = K_est;
                B_est_mat(subj_i, fce_i, disp_i, trial_i) = B_est;
                M_est_mat(subj_i, fce_i, disp_i, trial_i) = M_est;
                F_est_mat(subj_i, fce_i, disp_i, trial_i) = FIT;
            end
        end
    end

    % plot out the parameters
    x_fce = repmat(repmat(fce_list',1,3),1,1,7);
    x_dist= repmat(repmat(dist_list, 3,1),1,1,7);
    K_mat = reshape(K_est_mat(subj_i,:,:,:), 1, 3, 3, 7);
    K_arr = K_mat(:);
    B_mat = reshape(B_est_mat(subj_i,:,:,:), 1, 3, 3, 7);
    B_arr = B_mat(:);
    M_mat = reshape(M_est_mat(subj_i,:,:,:), 1, 3, 3, 7);
    M_arr = M_mat(:);
    F_mat = reshape(F_est_mat(subj_i,:,:,:), 1, 3, 3, 7);
    F_arr = F_mat(:);

    % if remove values
    ifremove = 0; 
    if (ifremove)
    K_arr(F_arr<80) = nan;
    B_arr(F_arr<80) = nan;
    M_arr(F_arr<80) = nan;
    F_arr(F_arr<80) = nan;
    end

    ifplot = 0;
    if (ifplot)
        figure();
        plot3(x_dist(:), x_fce(:), K_arr(:), 'x'); title('K');
        figure();
        plot3(x_dist(:), x_fce(:), B_arr(:), 'x'); title('B');
        figure();
        plot3(x_dist(:), x_fce(:), M_arr(:), 'x'); title('M');

    end
    
        K_arr_cell{subj_i} = K_arr;
        B_arr_cell{subj_i} = B_arr;
        M_arr_cell{subj_i} = M_arr;
        F_arr_cell{subj_i} = F_arr;
        
        x_fce_arr = x_fce(:);
        x_dist_arr= x_dist(:);

end
save('pulseParameterCompare.mat', 'K_arr_cell', 'B_arr_cell',...
    'M_arr_cell', 'F_arr_cell', ...
    'x_fce_arr', 'x_dist_arr', '-append');


%% compare the different settings 
pulse_width_list = [300 150 75];
distributionPlot(F_arr_cell);
title('valid-ness of fitting');

%% barplot the parameters 
% stiffness
figure('unit', 'inch', 'position', [0 0 6 3]);
for subj_i = 1:3
    axh(subj_i) = subplot(1,3,subj_i); hold on;
    data_mean = zeros(3,3);
    dat_std = zeros(3,3);
    clear lnh
    for fce_i = 1:3
        for dist_i = 1:3
            dat_idx = x_fce == fce_list(fce_i) ...
                & x_dist == dist_list(dist_i);
            dat_mean(fce_i,dist_i) = mean(K_arr_cell{subj_i}(dat_idx), 'omitnan');
            dat_std(fce_i,dist_i) = std(K_arr_cell{subj_i}(dat_idx), 'omitnan');
        end
%         lnh(fce_i) = plot(1:3, dat_mean(fce_i,:),  'linewidth', 3);
        lnh(fce_i) = errorbar(1:3, dat_mean(fce_i,:), dat_std(fce_i,:), 'linewidth', 3);
    end
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels({'2.5', '5.0', '7.5'});
    xlabel('target dist (cm)');
    legend(lnh, '15N', '20N', '25N');
    title(['pulse width ' num2str(pulse_width_list(subj_i))]);
    if (subj_i==1)
        ylabel('stiffness (N/m)');
    end
end
sgtitle('stiffness');
linkaxes(axh, 'xy');

% damping 
figure('unit', 'inch', 'position', [0 0 6 3]);
for subj_i = 1:3
    axh(subj_i) = subplot(1,3,subj_i); hold on;
    data_mean = zeros(3,3);
    dat_std = zeros(3,3);
    for fce_i = 1:3
        for dist_i = 1:3
            dat_idx = x_fce == fce_list(fce_i) ...
                & x_dist == dist_list(dist_i);
            dat_mean(fce_i,dist_i) = mean(B_arr_cell{subj_i}(dat_idx), 'omitnan');
            dat_std(fce_i,dist_i) = std(B_arr_cell{subj_i}(dat_idx), 'omitnan');
        end
%         plot(1:3, dat_mean(fce_i,:),  'linewidth', 3);
        lnh(fce_i) = errorbar(1:3, dat_mean(fce_i,:), dat_std(fce_i,:), 'linewidth', 3);
    end
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels({'2.5', '5.0', '7.5'});
    xlabel('target dist (cm)');
    legend(lnh, '15N', '20N', '25N');
    title(['pulse width ' num2str(pulse_width_list(subj_i))]);
    if (subj_i==1)
        ylabel('damping (N/m*s)');
    end
end
sgtitle('damping');
linkaxes(axh, 'xy');

% mass

figure('unit', 'inch', 'position', [0 0 6 3]);
for subj_i = 1:3
    axh(subj_i) = subplot(1,3,subj_i); hold on;
    data_mean = zeros(3,3);
    dat_std = zeros(3,3);
    for fce_i = 1:3
        for dist_i = 1:3
            dat_idx = x_fce == fce_list(fce_i) ...
                & x_dist == dist_list(dist_i);
            dat_mean(fce_i,dist_i) = mean(M_arr_cell{subj_i}(dat_idx), 'omitnan');
            dat_std(fce_i,dist_i) = std(M_arr_cell{subj_i}(dat_idx), 'omitnan');
        end
%         plot(dist_list, dat_mean(fce_i,:),  'linewidth', 3);
%         xlabel('target distance');
        lnh(fce_i) = errorbar(1:3, dat_mean(fce_i,:), dat_std(fce_i,:), 'linewidth', 3);
    end
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels({'2.5', '5.0', '7.5'});
    xlabel('target dist (cm)');
    legend(lnh, '15N', '20N', '25N');
    title(['pulse width ' num2str(pulse_width_list(subj_i))]);
    if (subj_i==1)
        ylabel('damping (N/m*s)');
    end
end
sgtitle('mass');
linkaxes(axh, 'xy');

% fitting

figure('unit', 'inch', 'position', [0 0 6 3]);
for subj_i = 1:3
    axh(subj_i) = subplot(1,3,subj_i); hold on;
    data_mean = zeros(3,3);
    dat_std = zeros(3,3);
    for fce_i = 1:3
        for dist_i = 1:3
            dat_idx = x_fce == fce_list(fce_i) ...
                & x_dist == dist_list(dist_i);
            dat_mean(fce_i,dist_i) = mean(F_arr_cell{subj_i}(dat_idx), 'omitnan');
            dat_std(fce_i,dist_i) = std(F_arr_cell{subj_i}(dat_idx), 'omitnan');
        end
%         plot(dist_list, dat_mean(fce_i,:),  'linewidth', 3);
        lnh(fce_i) = errorbar(1:3, dat_mean(fce_i,:), dat_std(fce_i,:), 'linewidth', 3);
    end
    xlim([0.5 3.5]);
    xticks([1 2 3]);
    xticklabels({'2.5', '5.0', '7.5'});
    xlabel('target dist (cm)');
    legend(lnh, '15N', '20N', '25N');
    title(['pulse width ' num2str(pulse_width_list(subj_i))]);
    if (subj_i==1)
        ylabel('fitting %');
    end
end
sgtitle('Fitting percentage');
linkaxes(axh, 'xy');

%% 
% See how these estimations overlap with the release measurements.  
% color_arr = colormap('lines'); close all;
load('pulseParameterCompare.mat', 'results1');
K_mean = results1.K_up_avg;
K_std  = results1.K_up_std;
B_mean = results1.B_up_avg; 
B_std  = results1.B_up_std;
M_mean = results1.M_up_avg; 
M_std  = results1.M_up_std; 


% plot the means and stds using errorbars... 
fh1 = figure();  hold on;
for fce_i = 1:3
%     for dist_i = 1:3
%     plot(1:3,K_mean(fce_i,:),'linewidth', 3, 'color', color_arr(fce_i,:));
    lnh(fce_i) = errorbar(1:3, K_mean(fce_i,:), K_std(fce_i,:), 'linewidth', 3, 'color', color_arr(fce_i,:));
%     end
end


% overlap the pulse-measurement using different marks...

figure(fh1);
for fce_i = 1:3
    for dist_i = 1:3
        dat_idx = x_fce_arr == fce_list(fce_i) & ...
            x_dist_arr == dist_list(dist_i);
        for subj_i = 1:3 % different pulses
            dat_val = K_arr_cell{subj_i}(dat_idx);
        
        % plot them 
        markers_arr = 'xos';
        lnh_tmp = plot(dist_i + subj_i/10 - 2/10, dat_val, 'Marker',markers_arr(subj_i), ...
            'color', color_arr(fce_i,:));
        lnh(3+subj_i) = lnh_tmp(1);
        end
    end
end

title('stiffness estimation');
xticks([1 2 3]);
xticklabels({'2.5' '5.0', '7.5'});
xlabel('target displacement (cm)');
ylabel('stiffness estimation');
legend(lnh, 'release 15N', 'release 20N', 'release 25N','300ms pulse', '150ms pulse', '75ms pulse');
