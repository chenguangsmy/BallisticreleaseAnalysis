%% Check the EMG data before release, see if they change through conditions...

% data format: 
% 20 (subj)*4 (dir)*3 (force)*3 (displacement)*9 (trials)*8 (muscles)

%% 1. format the into required format 
clear; clc; close all;

fdir = '../data/processedData/';
fname = {   'ss4379_4438_OPTupdate.mat';
            'ss4446_4467.mat';
            'ss4472_4495.mat';
            'ss4500_4524.mat';
            'ss4530_4563.mat';
            'ss4573_4587.mat'};
fsubj = {   1:6;
            7:9;
            10:12;
            13:15;
            16:18;
            19:20;};
dat_emg = cell(20,4,3,3,9);

for file_i = 1:length(fname)
    clearvars data
    load([fdir, fname{file_i}]);
    for s_i = 1:length(fsubj{file_i})
        s_i2 = fsubj{file_i}(s_i)
        for d_i = 1:4
            for f_i = 1:3
                for x_i = 1:3
                    for t_i = 1:9
                        dat_emg{s_i2,d_i,f_i,x_i,t_i}.t = data{s_i,d_i,f_i,x_i,t_i,1}.t;
                        dat_emg{s_i2,d_i,f_i,x_i,t_i}.emg = data{s_i,d_i,f_i,x_i,t_i,1}.emg;
                    end
                end
            end
        end
    end
end

dat_emg_pt = zeros(20,4,3,3,9,8);
t_range = [-0.5 0]; % right before release
for s_i = 1:20
    for d_i = 1:4
        for f_i = 1:3
            for x_i = 1:3
                for t_i = 1:9
                    t_sel = dat_emg{s_i,d_i,f_i,x_i,t_i}.t > t_range(1) & ...
                        dat_emg{s_i,d_i,f_i,x_i,t_i}.t < t_range(2);
                    dat_emg_pt(s_i,d_i,f_i,x_i,t_i,:) = mean(dat_emg{s_i,d_i,f_i,x_i,t_i}.emg(:,t_sel), 2);
                end
            end
        end
    end
end

save('data/emg_beforeRelease.mat', 'dat_emg', 'dat_emg_pt', '-v7.3');

%% 2. comparing plots...
clear; close all; clc
load('data/emg_beforeRelease.mat', 'dat_emg_pt');

% figure1, errorbar plot on each participant, each muscles (20*8)
for s_i = 1:20
    for m_i = 5:8
        for d_i = 1:1
            dat2plot = reshape(dat_emg_pt(s_i,d_i,:,:,:,m_i), 3, 3, 9);
            dat2plot_mean = mean(dat2plot,3);
            dat2plot_std = std(dat2plot,[],3);

            figure();
            colormap('lines');
            hold on;
            errorbar([25 50 75], dat2plot_mean(1,:), dat2plot_std(1,:), 'linewidth',2 );
            errorbar([25 50 75], dat2plot_mean(2,:), dat2plot_std(1,:), 'linewidth',2 );
            errorbar([25 50 75], dat2plot_mean(3,:), dat2plot_std(1,:), 'linewidth',2 );
            legend('15N', '20N', '25N');
            xlim([10 90]);
            xlabel('mm');
            title(['subj', num2str(s_i), ' muscle', num2str(m_i), ' direction', num2str(d_i)])
        end
    end
end


%% 3. (possible) statistics...