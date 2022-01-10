% checkEMGwithx
% plot the EMG data (ch3 and ch4) with the motion data parallel, see the
% result

sstmp = SessionScan(3717);

%%
%t_range = [2.60 2.63]*1e4;
%t_range = [3.465 3.501]*1e4;
%t_range = [3.5260 3.5624]*1e4;
t_range = [min(sstmp.data.t) + 60+120, max(sstmp.data.t) - 30];
idx = sstmp.data.t>t_range(1) & sstmp.data.t<t_range(2);
t = sstmp.data.t(idx);
x = sstmp.data.x(2,idx);
f = sstmp.data.f(2,idx);
ts =sstmp.data.ts(idx);
emg=sstmp.data.emg(:,idx);
emg_processed = zeros(size(emg));

%idx_tsinvalid = ~(ts>3 & ts<7);
%emg(:,idx_tsinvalid) = 0;

%% try to do a bunch of pre-processing... 
chi = 6;
% 1. do the 100Hz low pass filter
fs = 500; % data frequency
lpf1 = 100;
%emg1 = lowpass(emg(chi,:), lpf1, fs);
emg_stp1 = highpass(emg(chi,:), lpf1, fs);
% 2. mean-centered and scaled by standard deviation... 
emg_stp2 = (emg_stp1 - mean(emg_stp1)) / std(emg_stp1);
% 3. squared, and do low-pass filter of 30Hz
lpf2 = 30;
emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
% 3. squre root transform and *2
emg_stp4 = sqrt(emg_stp3)*2;

axh(1) = subplot(2,1,1); 
plot(t,emg_stp1);
axh(2) = subplot(2,1,2);
plot(t,emg_stp4);
linkaxes(axh, 'x');

%% try to do bunch of pre-processing, iterate with channels
for chi = 1:8
    % 1. do the 100Hz low pass filter
    fs = 500; % data frequency
    lpf1 = 100;
    % emg1 = bandstop(emg(ch_i,:),[118.5 121.5],fs);
    emg_stp1 = lowpass(emg(chi,:), lpf1, fs);
    % 2. mean-centered and scaled by standard deviation...
    emg_stp2 = (emg_stp1 - mean(emg_stp1)) / std(emg_stp1);
    % 3. squared, and do low-pass filter of 30Hz
    lpf2 = 30;
    emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
    % 3. squre root transform and *2
    emg_stp4 = sqrt(emg_stp3)*2;                                                   % bad name, need chagne
    
    emg_processed(chi,:) = emg_stp4;
    
    figure();
    axh(1) = subplot(2,1,1);
    plot(t,emg_stp1);
    title(['channel' num2str(chi)]);
    axh(2) = subplot(2,1,2);
    plot(t,emg_stp4);
    linkaxes(axh, 'x');
end


%% % a simple plot of rms emg, to see if it co-variate with the data
fh = figure();
axh(1) = subplot(3,1,1);
plot(t,x); title('position'); ylabel('m');
axh(2) = subplot(3,1,2);
plot(t,f); title('force'); ylabel('N');
axh(3) = subplot(3,1,3);
ylim([-10 20])
%plot(t,abs(emg)); title('EMG'); ylabel('?');
%dat = smooth(abs(emg(1,:)));
rms = sqrt((emg(1,:) - mean(emg(1,:))).^2);
dat = smooth(rms);
plot(t,dat); title('EMG'); ylabel('?');
ylim([-5000 5000])

linkaxes(axh, 'x');
%xlim([2.6 2.61]*10e6)

%% Also plot the ts out to see the task states 
fh = figure();
axh(1) = subplot(4,1,1);
plot(t,x); title('position'); ylabel('m');
axh(2) = subplot(4,1,2);
plot(t,f); title('force'); ylabel('N');
axh(3) = subplot(4,1,3);
ylim([-10 20])
%plot(t,abs(emg)); title('EMG'); ylabel('?');
%dat = smooth(abs(emg(1,:)));
rms = sqrt((emg(2,:) - mean(emg(2,:))).^2);
%dat = smooth(rms);
%dat = smooth(rms, 'lowess');
dat = smooth(rms, 'moving', 500);
%plot(t,dat); title('EMG'); ylabel('?');
plot(t, emg_stp4); title('EMG'); ylabel('EMG');
%ylim([-5000 5000])

%xlim([2.6 2.61]*10e6)
axh(4) = subplot(4,1,4);
plot(t,ts);
ylim([0 8])
title('task state');
linkaxes(axh, 'x');

%% plot 4 pairs of muscles activities against the force level 
emg_means = zeros(9,8);
emg_stds = zeros(9,8);
ss_list = [3710 3709 3711, 3715, 3714, 3713, 3716, 3718, 3717];
for ssi = 1:9
fh = figure();
ssnum = ss_list(ssi);
sstmp = SessionScan(ssnum);
clf;

% data
t_range = [min(sstmp.data.t) + 60+120, max(sstmp.data.t) - 50];
idx = sstmp.data.t>t_range(1) & sstmp.data.t<t_range(2);
t = sstmp.data.t(idx);
x = sstmp.data.x(2,idx);
f = sstmp.data.f(2,idx);
ts =sstmp.data.ts(idx);
emg=sstmp.data.emg(:,idx);
emg_processed = zeros(size(emg));


% precess emg
for chi = 1:8 % iterate through channels
    % 1. do the 100Hz low pass filter
    fs = 500; % data frequency
    lpf1 = 100;
    % emg1 = bandstop(emg(ch_i,:),[118.5 121.5],fs);
    emg_stp1 = lowpass(emg(chi,:), lpf1, fs);
    % 2. mean-centered and scaled by standard deviation...
    emg_means(ssi,chi) = mean(emg_stp1);
    emg_stds(ssi,chi) = std(emg_stp1);
    emg_stp2 = (emg_stp1 - mean(emg_stp1)) / std(emg_stp1);
    % 3. squared, and do low-pass filter of 30Hz
    lpf2 = 30;
    emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
    % 3. squre root transform and *2
    emg_stp4 = sqrt(emg_stp3)*2;                                                   % bad name, need chagne
    
    emg_processed(chi,:) = emg_stp4;
    
%     figure();
%     axh(1) = subplot(2,1,1);
%     plot(t,emg_stp1);
%     title(['channel' num2str(chi)]);
%     axh(2) = subplot(2,1,2);
%     plot(t,emg_stp4);
%     linkaxes(axh, 'x');
end

% plot
axh(1) = subplot(5,1,1);
plot(t,f); title('force'); ylabel('N');
axh(2) = subplot(5,1,2);
plot(t,emg_processed(1:2,:)); %title('EMG12'); %ylabel('N');
axh(3) = subplot(5,1,3);
plot(t,emg_processed(3:4,:)); %title('EMG34'); %ylabel('N');
axh(4) = subplot(5,1,4);
plot(t,emg_processed(5:6,:)); %title('EMG56'); %ylabel('N');
axh(5) = subplot(5,1,5);
plot(t,emg_processed(7:8,:)); %title('EMG78'); %ylabel('N');
linkaxes(axh, 'x');
linkaxes(axh(2:end), 'y');
ylim(axh(5), [0, 10]);
sgtitle(num2str(ssnum));
end

%% append the data for several sessions
ss_list = [3710 3709 3711, 3715, 3714, 3713, 3716, 3718, 3717];
emg_all = [];
for ssi = 1:9
    ssnum = ss_list(ssi);
    sstmp = SessionScan(ssnum);
    % start from the seccond sucessful trial
    trialsSuccList = find([sstmp.trials.outcome]==1);
    t_range = [sstmp.trials(trialsSuccList(2)).time_orn(1), ...
               sstmp.trials(trialsSuccList(end-2)).time_orn(end)];
    %t_range = [min(sstmp.data.t) + 60+120, max(sstmp.data.t) - 50];
    idx = sstmp.data.t>t_range(1) & sstmp.data.t<t_range(2);
    emg=sstmp.data.emg(:,idx);
    emg_all = [emg_all, emg];
end

% get the mean and std coefficient
emg_means = mean(emg_all,2);
emg_stds = std(emg_all,0,2);

% now process the data with the coefficient: 

%% 
for ssi = 1:9
fh = figure();
ssnum = ss_list(ssi);
sstmp = SessionScan(ssnum);
clf;

% data
t = sstmp.data.t(:,idx);
f = sstmp.data.f(:,idx);
emg=sstmp.data.emg(:,idx);
emg_processed0= zeros(size(emg));
emg_processed = zeros(size(emg));

% precess emg
for chi = 1:8 % iterate through channels
    % 1. do the 100Hz low pass filter
    fs = 500; % data frequency
    lpf1 = 100;
    % emg1 = bandstop(emg(ch_i,:),[118.5 121.5],fs);
    emg_stp1 = highpass(emg(chi,:), lpf1, fs);
    % 2. mean-centered and scaled by standard deviation...
    emg_stp2 = (emg_stp1 - emg_means(chi)) / emg_stds(chi);
    % 3. squared, and do low-pass filter of 30Hz
    lpf2 = 30;
    emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
    % 3. squre root transform and *2
    emg_stp4 = sqrt(emg_stp3)*2;                                                   % bad name, need chagne
    
    [emg_stp5, emg_stp5n] = envelope(real(emg_stp4), 5, 'peak');
    emg_processed0(chi,:)=emg_stp4;
    emg_processed(chi,:) = emg_stp5;%emg_stp4;
    
end

% plot
axh(1) = subplot(5,1,1);
plot(t,f); title('force'); ylabel('N');
ylim([0 30]);
axh(2) = subplot(5,1,2); hold on;
plot(t,emg_processed(1:2,:)); %title('EMG12'); %ylabel('N');
plot(t,emg_processed0(1:2,:), '.'); %title('EMG12'); %ylabel('N');
axh(3) = subplot(5,1,3); hold on;
plot(t,emg_processed(3:4,:)); %title('EMG34'); %ylabel('N');
plot(t,emg_processed0(3:4,:), '.'); %title('EMG12'); %ylabel('N');
% hold on; plot(t,emg_processed(3,:)); 
% plot(t,emg_processed0(3,:)); %title('EMG34'); %ylabel('N');

axh(4) = subplot(5,1,4); hold on;
plot(t,emg_processed(5:6,:)); %title('EMG56'); %ylabel('N');
plot(t,emg_processed0(5:6,:), '.'); %title('EMG12'); %ylabel('N');
axh(5) = subplot(5,1,5); hold on;
plot(t,emg_processed(7:8,:)); %title('EMG78'); %ylabel('N');
plot(t,emg_processed0(7:8,:), '.'); %title('EMG12'); %ylabel('N');
linkaxes(axh, 'x');
linkaxes(axh(2:end), 'y');
ylim(axh(5), [0, 10]);
sgtitle(num2str(ssnum));
end

%% plot the emg activity according to different task conditions 
ss_list = [...
    3710    3709    3711;   %   15N 
    3715    3714    3713;   %   20N
    3716    3718    3717];  %   25N
for fce_i = 1:3
    for tar_i = 1:3
        sstmp((fce_i-1)*3+tar_i) = SessionScan(ss_list(fce_i,tar_i));
    end
end
%% ... continues
fh = figure();
% fh1 = figure();
col_type = colormap('lines');
for fce_i = 1:3 
    for tar_i = 1:3
        sstmp1 = sstmp((fce_i-1)*3+tar_i);
        % just plot the raw data overlay eath other
%         figure(fh1);
        %sstmp = SessionScan(ss_list(fce_i, tar_i));
        figure(fh);
        for chi = 1:8
        axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
        trials_list = find([sstmp1.trials.outcome] == 1);

        for trial_i = 1:trials_list
            %chi = 3;    % only plot one muscle here
            plot(sstmp1.trials(trial_i).data.t_shift, ...
                 sstmp1.trials(trial_i).data.emg(chi,:), ...
                 'Color', col_type(tar_i,:));
        end
        end
    end
    xlim([-3 2]); % sec
    ylim([-1 20]);
end