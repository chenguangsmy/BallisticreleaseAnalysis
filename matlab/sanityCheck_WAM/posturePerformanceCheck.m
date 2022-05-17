% Re-mount WAM to a posture of "joystick", that mount horizontally and the
% small arm and handle is pointing up. 
% contents: 
% figure 1, check the force exertion on different axis;
% figure 2, check the displacement readout (hysteresis) on different axis;
% figure 3, check the release response on different axis;

% %% 
% % Tidy up the data for the plot
% 
% % tidy up data in 6D (direction * force)
% % Direction: 0, 2, 4, 6~ 0, 90, 180, 270 degrees
% %   corresponding pulse direction: -x, -y, +x, +y
% ss_list = [4048 4052 4053];
% celltmp_combine = {};
% for ss_i = 1:length(ss_list)
% sstmp = SessionScan(ss_list(ss_i)); 
% celltmp = sstmp.export_as_formatted(1);
% celltmp1 = reshape(celltmp(1,1,:,2), size(celltmp,3),1);
% celltmp_combine = {celltmp_combine{:} celltmp1{:}};
% end
% %
% pertfce_list = [ 6  12  18  24  30];
% trials_list = cell(4,5);
% xyi = [1 2 1 2];
% xy_npi = [-1 -1 1 1];
% 
% % get the index from the pilled-up data 
% for diri = 1:4
%     for fcei = 1:5
%         for triali = 1:length(celltmp_combine)
% %             max(abs(celltmp_combine{trial_i}.Fp(xyi(dir_i),:))) * xy_npi(dir_i)
%             [~,fp_idx] = max(abs(celltmp_combine{triali}.Fp(xyi(diri),:)));
%             fp = celltmp_combine{triali}.Fp(xyi(diri),fp_idx);
%             pertidx = find(pertfce_list(fcei)*xy_npi(diri) == round(fp));
%             if isempty(pertidx)
%                 continue;
%             end
%             trials_list{diri,fcei} = [trials_list{diri,fcei}, triali];
%         end
%     end
% end
% % get the data tidied up 
% data = cell(1,4,5,1,10,3);
% % 1 subject, 4 directions, 5 force levels, 1 distance, 10 trials, 3 pert
% for diri = 1:4
%     for fcei = 1:5
%         trials_num = length(trials_list{diri, fcei});
%         data(1,diri,fcei,1,1:trials_num,2) = ...
%             celltmp_combine(trials_list{diri,fcei});
%     end
% end
% save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4048_4053.mat', 'data')

% % The data for the slow perturbation 
% ss_list = [4049 4051];
% celltmp_combine = {};
% for ss_i = 1:length(ss_list)
%     sstmp = SessionScan(ss_list(ss_i));
%     celltmp = sstmp.export_as_formatted(1);
%     celltmp1 = reshape(celltmp(1,1,:,2), size(celltmp,3),1);
%     celltmp_combine = {celltmp_combine{:} celltmp1{:}};
% end
% %
% pertfce_list = [ 6  12  18  24  30];
% trials_list = cell(4,5);
% xyi = [1 2 1 2];
% xy_npi = [-1 -1 1 1];
% %
% % get the index from the pilled-up data 
% for diri = 1:4
%     for fcei = 1:5
%         for triali = 1:length(celltmp_combine)
%             Fp_nonzeroidx = find(sum(celltmp_combine{triali}.Fp(:,:))~=0);
% %             max(abs(celltmp_combine{trial_i}.Fp(xyi(dir_i),:))) * xy_npi(dir_i)
%             [~,fp_idx] = max(abs(celltmp_combine{triali}.Fp(xyi(diri),Fp_nonzeroidx(1)+[0:500])));
% %             clf; plot(abs(celltmp_combine{triali}.Fp(:,Fp_nonzeroidx(1)+[0:500]))');
%             fp = celltmp_combine{triali}.Fp(xyi(diri),Fp_nonzeroidx(fp_idx));
%             pertidx = find(pertfce_list(fcei)*xy_npi(diri) == round(fp));
%             if isempty(pertidx)
%                 continue;
%             end
%             trials_list{diri,fcei} = [trials_list{diri,fcei}, triali];
%         end
%     end
% end
% % get the data tidied up 
% data = cell(1,4,5,1,10,3);
% % 1 subject, 4 directions, 5 force levels, 1 distance, 10 trials, 3 pert
% for diri = 1:4
%     for fcei = 1:5
%         trials_num = length(trials_list{diri, fcei});
%         data(1,diri,fcei,1,1:trials_num,2) = ...
%             celltmp_combine(trials_list{diri,fcei});
%     end
% end
% save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4049_4051.mat', 'data')


%% Plot1, check the pulse force exertion on different axis 
clear; 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4048_4053.mat', 'data');

fh = figure(); 
color_arr = colormap('lines'); 
close(fh); 
clear axh lnh
fh1 = figure('unit', 'inch', 'position', [0 0 4 8]);
xyi = [1 2 1 2];
fce_peaks = cell(4,5);
for fcei = 1:5
for diri = 1:4
data1 = reshape(data(1,diri,fcei,1,1:10,2), 10, 1);
for triali = 1:length(data1)
    trialtmp = data1{triali};
    if isempty(trialtmp) % not a trial
        continue;
    end
    t_idx = find(trialtmp.Fp(xyi(diri),:)~=0); % perturb
    axh(fcei,xyi(diri)) = subplot(5,2,(fcei-1)*2 + xyi(diri)); 
    hold on;
    
    t_shift = trialtmp.t(t_idx(1));
    trialtmp.t_shift = trialtmp.t - t_shift;
    t_avg_range = [-0.1 0]; % before perturbation happens
    f_avg_idx = trialtmp.t_shift > t_avg_range(1) & trialtmp.t_shift < t_avg_range(2);
    f_avg = mean(trialtmp.f(:,f_avg_idx),2);
    
    lnh{fcei}(1) = plot(trialtmp.t_shift, -trialtmp.Fp(xyi(diri),:),  'Color', color_arr(2,:), 'LineWidth', 2);
    lnh{fcei}(2) = plot(trialtmp.t_shift, trialtmp.f(xyi(diri),:)-f_avg(xyi(diri)), 'Marker','.', 'Color', color_arr(1,:));
    
    % get fce_peaks here
    t_peak_range = [0, 0.4];
    f_peak_idx = trialtmp.t_shift > t_peak_range(1) & trialtmp.t_shift < t_peak_range(2);
    f_subavg =  trialtmp.f(xyi(diri),f_peak_idx)-f_avg(xyi(diri));
    [~, f_peakval_idx] = max(abs(f_subavg));
    fce_peaks{diri,fcei} = [fce_peaks{diri,fcei} f_subavg(f_peakval_idx)];  
end 
    if ((fcei-1) == 4)
        xlabel('t (s)' );
    end
    ylabel('force (N)');
    subplot(5,2,(fcei-1)*2 + 1); grid on; title(['x+- pert' ]);
    subplot(5,2,(fcei-1)*2 + 2); grid on; title(['y+- pert' ]);
end
linkaxes(axh(fcei,:), 'xy');
xlim(axh(fcei,:), [0 0.4]); % perturb during release 
end
linkaxes(axh(:), 'x');
legend(lnh{1}, {'command Force', 'censored force'});
% sgtitle('during hold');
sgtitle(fh1, 'Force generation in gaussian pulse');

% save the figure to a place

% Plot2, check the pulse force exertion on different axis, summary
fh = figure('unit','inch','position',[0 0 3 3]); 
hold on;
fce_sign = [1 1 -1 -1];
fce_val = [6:6:30];
% concatinate the x-direction peaks 
Xx = []; Xy = [];
for diri = [1 3]
    for fcei = 1:5
        Xy = [Xy fce_peaks{diri,fcei}];
        Xx = [Xx ones(size(fce_peaks{diri,fcei}))* fce_val(fcei) * fce_sign(diri)];
    end
end
plot(Xx, Xy, '.', 'markersize', 10)
% concatinate the y-direction peaks
Yx = []; Yy = [];
for diri = [2 4]
    for fcei = 1:5
        Yy = [Yy fce_peaks{diri,fcei}];
        Yx = [Yx ones(size(fce_peaks{diri,fcei}))* fce_val(fcei) * fce_sign(diri)];
    end
end
plot(Yx, Yy, '.', 'markersize', 10)
grid on;
line([-30 30], [-30 30], 'color', [0.5 0.5 0.5]);
legend('-x', '-y', 'reference');
xlabel('F command (N)');
ylabel('F censored(N)');
title('force peak in gaussian pulse');

% save the figure to a place

%% Figure 3, the position hysterisis in two directions
clear; 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4049_4051.mat', 'data');
fh = figure(); 
color_arr = colormap('lines'); 
close(fh); 
clear axh lnh
fh1 = figure('unit', 'inch', 'position', [0 0 4 8]);
xyi = [1 2 1 2];
fce_peaks = cell(4,5);
for fcei = 1:5
for diri = 1:4
data1 = reshape(data(1,diri,fcei,1,1:10,2), 10, 1);
for triali = 1:length(data1)
    trialtmp = data1{triali};
    if isempty(trialtmp) % not a trial
        continue;
    end
    t_idx = find(trialtmp.Fp(xyi(diri),:)~=0); % perturb
    axh(fcei,xyi(diri)) = subplot(5,2,(fcei-1)*2 + xyi(diri)); 
    hold on;
    
    t_shift = trialtmp.t(t_idx(1));
    trialtmp.t_shift = trialtmp.t - t_shift;
    t_avg_range = [-0.1 0]; % before perturbation happens
    x_avg_idx = trialtmp.t_shift > t_avg_range(1) & trialtmp.t_shift < t_avg_range(2);
    x_avg = mean(trialtmp.x(:,x_avg_idx),2);
%     x_avg = x_avg*0;
    
    lnh{fcei}(1) = plot(trialtmp.t_shift, zeros(size(trialtmp.t_shift)),  'Color', color_arr(2,:), 'LineWidth', 2);
    lnh{fcei}(2) = plot(trialtmp.t_shift, trialtmp.x(xyi(diri),:)-x_avg(xyi(diri)), 'Marker','.', 'Color', color_arr(1,:));
    
    % get fce_peaks here
%     t_peak_range = [0, 0.4];
%     f_peak_idx = trialtmp.t_shift > t_peak_range(1) & trialtmp.t_shift < t_peak_range(2);
%     f_subavg =  trialtmp.f(xyi(diri),f_peak_idx)-f_avg(xyi(diri));
%     [~, f_peakval_idx] = max(abs(f_subavg));
%     fce_peaks{diri,fcei} = [fce_peaks{diri,fcei} f_subavg(f_peakval_idx)];  
end 
    if ((fcei-1) == 4)
        xlabel('t (s)' );
    end
    ylabel('force (N)');
    subplot(5,2,(fcei-1)*2 + 1); grid on; title(['x+- pert' ]);
    subplot(5,2,(fcei-1)*2 + 2); grid on; title(['y+- pert' ]);
end
linkaxes(axh(fcei,:), 'xy');
xlim(axh(fcei,:), [0 4]); % perturb during release 
end
linkaxes(axh(:), 'x');
legend(lnh{1}, {'theory position', 'read position'});
% sgtitle('during hold');
sgtitle(fh1, 'position reading in gaussian pulse');

%% %
% now plot the force against the displacement
clear axh lnh
fh1 = figure('unit', 'inch', 'position', [0 0 6 3]); 
for diri = 1:2
    for fcei = 1:5
data1 = reshape(data(1,diri,fcei,1,1:10,2), 10, 1);
for triali = 1%:length(data1)
    trialtmp = data1{triali};
    if isempty(trialtmp) % not a trial
        continue;
    end
    t_idx = find(trialtmp.Fp(xyi(diri),:)~=0); % perturb
%     axh(fcei,xyi(diri)) = subplot(5,2,(fcei-1)*2 + xyi(diri)); 
    axh(diri) = subplot(1,2,diri);  grid on;
    hold on;
    
    t_shift = trialtmp.t(t_idx(1));
    trialtmp.t_shift = trialtmp.t - t_shift;
    t_avg_range = [-0.1 0]; % before perturbation happens
    x_avg_idx = trialtmp.t_shift > t_avg_range(1) & trialtmp.t_shift < t_avg_range(2);
    x_avg = mean(trialtmp.x(:,x_avg_idx),2);
    f_avg = mean(trialtmp.f(:,x_avg_idx),2);
%     x_avg = x_avg*0;
    
%     lnh{fcei}(1) = plot(trialtmp.t_shift, zeros(size(trialtmp.t_shift)),  'Color', color_arr(2,:), 'LineWidth', 2);
%     lnh{fcei}(2) = plot(trialtmp.t_shift, trialtmp.x(xyi(diri),:)-x_avg(xyi(diri)), 'Marker','.', 'Color', color_arr(1,:));
    lnh{1}(fcei) = plot(-(trialtmp.f(xyi(diri),:)-f_avg(xyi(diri))),...
            trialtmp.x(xyi(diri),:)-x_avg(xyi(diri)),...
            'Marker','.', 'Color', color_arr(fcei,:));
    
    % get fce_peaks here
%     t_peak_range = [0, 0.4];
%     f_peak_idx = trialtmp.t_shift > t_peak_range(1) & trialtmp.t_shift < t_peak_range(2);
%     f_subavg =  trialtmp.f(xyi(diri),f_peak_idx)-f_avg(xyi(diri));
%     [~, f_peakval_idx] = max(abs(f_subavg));
%     fce_peaks{diri,fcei} = [fce_peaks{diri,fcei} f_subavg(f_peakval_idx)];  
end 
    if ((fcei-1) == 4)
        xlabel('t (s)' );
    end
    xlabel('censored force (N)');
    ylabel('x_{err} (m)');
%     subplot(5,2,(fcei-1)*2 + 1); grid on; title(['x+- pert' ]);
%     subplot(5,2,(fcei-1)*2 + 2); grid on; title(['y+- pert' ]);
end
% linkaxes(axh(fcei,:), 'xy');
% xlim(axh(fcei,:), [0 4]); % perturb during release 
if (diri == 1)
    title('direction x')
else
    title('direction y')
end
end
linkaxes(axh(:), 'xy');

legend(lnh{1}, {'6N', '12N', '18N', '24N', '30N'});
% sgtitle('during hold');
sgtitle(fh1, 'position error when force exerted');

%% 
% The ballistic-release task when WAM is mounting horizontally posture. 
% x- major joint4;
% y- major joint3, (also could be joint 1, need some data) 

% plot out the data of the displacement, velocity, force and trajectory for
% each of the directions 

% For considering of time, I did 4 directions with only one force -
% displacement combinations (20N * 5cm). 
% Subject was sitting on the -x direction of the robot (face the control
% panel of the robot). People was face the wall. The below mentioned (FLBR)
% correspond to (+x, +y, -x, -y) specifically. 

% 4027 front; 4028: left; 4029: back; 4030: right;

% load data and put them into the format 


% tidy up data in 6D (direction * force)
% Direction: 0, 2, 4, 6~ 0, 90, 180, 270 degrees
ss_list = [4035 4036 4037 4038];
data = cell(1,4,5,1,30,3);
% celltmp_combine = {};
for ss_i = 1:length(ss_list)
    sstmp = SessionScan(ss_list(ss_i));
    celltmp = sstmp.export_as_formatted(1);
    celltmp1 = reshape(celltmp(1,1,:,1), size(celltmp,3),1);
    % celltmp_combine = {celltmp_combine{:} celltmp1{:}};
    trialnum = size(celltmp,3);
    if trialnum > 30
        data(1,ss_i,1,1,:,1) = celltmp1(1:30);
    else
        data(1,ss_i,1,1,1:trialnum,1) = celltmp1(:);
    end
end
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4035_4038.mat', 'data');
%% plot the trajectory, velocity and force 
% 
% plot the trajectory
clear; 
fh = figure(); 
color_arr = colormap('lines');
close(fh);  

fh = figure('unit', 'inch', 'position', [0 0 3 3]); hold on;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4035_4038.mat', 'data');
for dir_i = 1:4 
    data1 = reshape(data(1,dir_i,1,1,:,1), size(data,5), 1);
    for trial_i = 1:length(data1)
        if (isempty(data1{trial_i})) 
            continue;
        end
        % find the release data 
        t_range = [-0.1 1];
        t_idx = data1{trial_i}.t>t_range(1) &  data1{trial_i}.t<t_range(2);
%         sum(t_idx) % check the length of tiem , should be ~600
        posx = data1{trial_i}.x(1,t_idx);
        posy = data1{trial_i}.x(2,t_idx);
        posz = data1{trial_i}.x(3,t_idx);  % when consider z
%         plot(posx, posy, 'Color', color_arr(dir_i, :));
        plot3(posx, posy, posz, 'Color', color_arr(dir_i, :));
    end
end
xlim([-0.56 -0.40]);
ylim([-0.08 0.08]); 

xlabel('x (m)'); 
ylabel('y (m)');
try 
    zlabel('z (m)');
    zlim([0.43 0.59]); 
catch 
    
end
title('trajectory in the horizontal WAM mounting')

% set(gca, 'view', [90 0])

%% plot the velocity and force 


fh(1) = figure('unit', 'inch', 'position', [0 0 4 4]); hold on; % every trial
fh(2) = figure('unit', 'inch', 'position', [0 0 4 4]); hold on; % mean plot
fh_test = figure(); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4035_4038.mat', 'data');
dir_list = [1 2 1 2];
sign_list = [1 1 -1 -1];
lnh = zeros(4,1);
freq = 500;

        
t_range = [-0.2 1.5];
t_grid = t_range(1):1/freq:t_range(2);
trials_num = 30;
t_grid_num = length(t_grid);


for dir_i = 1:4 
    data1 = reshape(data(1,dir_i,1,1,:,1), size(data,5), 1);
    
    x_mat = nan(trials_num,t_grid_num); 
    v_mat = nan(trials_num,t_grid_num); 
    f_mat = nan(trials_num,t_grid_num); 
    
    for trial_i = 1:length(data1)
        if isempty(data1{trial_i})
            continue;
        end
        % find the release data 
        t_idx = data1{trial_i}.t>t_range(1) &  data1{trial_i}.t<t_range(2);
%         sum(t_idx) % check the length of time , should be ~600
        t = data1{trial_i}.t(t_idx);
        x = data1{trial_i}.x(dir_list(dir_i),t_idx);
        v = data1{trial_i}.v(dir_list(dir_i),t_idx) * sign_list(dir_i);
        f = data1{trial_i}.f(dir_list(dir_i),t_idx) * sign_list(dir_i);  
        
        x_itp = interp1(data1{trial_i}.t, data1{trial_i}.x(dir_list(dir_i),:), t_grid, 'linear','extrap');
        v_itp = interp1(data1{trial_i}.t, data1{trial_i}.v(dir_list(dir_i),:), t_grid, 'linear','extrap');
        f_itp = interp1(data1{trial_i}.t, data1{trial_i}.f(dir_list(dir_i),:), t_grid, 'linear','extrap');
        x_mat(trial_i,:) = x_itp;
        v_mat(trial_i,:) = v_itp;
        f_mat(trial_i,:) = f_itp;
        
        ifplot = 0;
        if (ifplot)
            figure(fh_test); 
            clf; 
            subplot(3,1,1); hold on; 
            plot(t,x); 
            plot(t_grid, x_itp, '.');
            subplot(3,1,2); hold on; 
            plot(t,v * sign_list(dir_i)); 
            plot(t_grid, v_itp, '.');
            subplot(3,1,3); hold on; 
            plot(t,f * sign_list(dir_i)); 
            plot(t_grid, f_itp, '.');
        end
        
        figure(fh(1)); % each trial
        subplot(2,1,1); hold on;
        lnh(dir_i) = plot(t, v, 'Color', color_arr(dir_i, :));
        subplot(2,1,2); hold on;
        plot(t, f, 'Color', color_arr(dir_i, :));
        
    end
    figure(fh(2));
    f_mean = nanmean(f_mat);
    f_std = nanstd(f_mat);
    %         x_mean = nanmean(x_mat);
    %         x_std = nanstd(x_mat);
    v_mean = nanmean(v_mat);
    v_std = nanstd(v_mat);
    subplot(2,1,1); hold on;
    lnh(dir_i) = plot(t_grid, v_mean * sign_list(dir_i), 'Color', color_arr(dir_i, :), 'LineWidth', 3);
    patch([t_grid t_grid(end:-1:1)], [v_mean-v_std v_mean(end:-1:1)+v_std(end:-1:1)] * sign_list(dir_i), color_arr(dir_i, :), 'FaceAlpha',.3);
    subplot(2,1,2); hold on;
    lnh(dir_i) = plot(t_grid, f_mean * sign_list(dir_i), 'Color', color_arr(dir_i, :), 'LineWidth', 3);
    patch([t_grid t_grid(end:-1:1)], [f_mean-f_std f_mean(end:-1:1)+f_std(end:-1:1)] * sign_list(dir_i), color_arr(dir_i, :), 'FaceAlpha',.3);
end

for fhi = 1:2
    figure(fh(fhi));
    subplot(2,1,1);
    ylim([-0.2 0.6]); 
    xlabel('t (s)'); ylabel('v (m/s)');
    legend(lnh, '+x', '+y', '-x', '-y');

    subplot(2,1,2); 
    xlim([-0.1 0.6]);
    ylim([-10 20]); 
    xlabel('t (s)'); ylabel('f (N)');
    
    if (fhi == 1)
        sgtitle('velocity and force on 4 directions');
    else
        sgtitle('avrage v and f on 4 directions');
    end
end