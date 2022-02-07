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
% ss_list = [...
%     3710    3709    3711;   %   15N 
%     3715    3714    3713;   %   20N
%     3716    3718    3717];  %   25N
ss_list = [...
    3866    3865    3864;   %   15N 
    3867    3869    3868;   %   20N
    3871    3870    3872];  %   25N
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
        ylabel(['channel ' num2str(chi)]);
        if (chi == 1) 
            title(['force #' num2str(fce_i)]);
        end
        
        ylim([0 0.5]);
        end
        
    end
    xlim([-3 2]); % sec
    ylim([-1 20]);
end

%% plot several sessions % ss3716, 3718, ss4

%% ... continues
% sstmp1 = SessionScan(3710);
close all;
fh = figure();
% fh1 = figure();
col_type = colormap('lines');
% just plot the raw data overlay eath other
%         figure(fh1);
%sstmp = SessionScan(ss_list(fce_i, tar_i));
for fce_i = 3:-1:1
    sstmp1 = sstmp((fce_i-1)*3+1);
%     sstmp1 = sstmp(fce_i);
figure(fh);
for chi = 8
    %axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
    trials_list = find([sstmp1.trials.outcome] == 1);
    
    for trial_i = trials_list
        axh(1) = subplot(2,1,1); hold on;
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.x(2,:), ...
            'Color', col_type(fce_i,:));
        %chi = 3;    % only plot one muscle here
        axh(2) = subplot(2,1,2); hold on;
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.emg(chi,:), ...
            'Color', col_type(fce_i,:));
    end
end
subplot(axh(1));
title('position');
subplot(axh(2));
title(['EMG-Signal, ch' num2str(chi)]);
linkaxes(axh, 'x');
xlim([-3 2]); % sec
% ylim([-1 5]);
ylim([-1 5]/10);
end

%% 
%%%%%%%% try the EMG data with only ch3 (bicep) and ch4 (tricep) when doing
%%%%%%%% the data
%%%%%%%% 25N, 5cm, 2.5cm, 7.5cm
%%%%%%%% 5cm, 15N, 20N, 25N
% sstmp1 = SessionScan(3856);
% close all;
col_type = colormap('lines');

ss_num = [3857 3856 3858];
for ssi = 1:length(ss_num) 
    sstmp(ssi) = SessionScan(ss_num(ssi));
end
%
% just plot the raw data overlay eath other
fh1 = figure('Unit', 'inch', 'Position', [0 0 3 6]);
figure(fh1);
%sstmp = SessionScan(ss_list(fce_i, tar_i));
for fce_i = 1:3
     sstmp1 = sstmp(fce_i);
% fce_i = 1
% figure(fh);
    %axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
    trials_list = find([sstmp1.trials.outcome] == 1);
    
    for trial_i = trials_list
        axh(1) = subplot(4,1,1); hold on;       % posiitoin
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.x(2,:), ...
            'Color', col_type(fce_i,:));
        ylabel('m'); title('position');
        %chi = 3;    % only plot one muscle here
        axh(2) = subplot(4,1,2); hold on;       % force 
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.f(2,:), ...
            'Color', col_type(fce_i,:));
        ylabel('N'); title('force');
        axh(3) = subplot(4,1,3); hold on;       % ch 3: elbow flexor
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.emg(3,:), ...
            'Color', col_type(fce_i,:));    
        ylabel('EMG activity'); title('Elbow Flesor');
        axh(4) = subplot(4,1,4); hold on;       % % ch 4: elbow flexor
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.emg(4,:), ...
            'Color', col_type(fce_i,:));
        ylabel('EMG activity'); title('Elbow Extensor')
    end
    
linkaxes(axh, 'x');
% xlim([-3 2]); % secxlim([-1 0.5]); % sec
xlim([-1 0.5]); % sec
% ylim([-1 5]);
set(axh(3), 'YLim', [0 40]);
set(axh(4), 'YLim', [0 20]);
end
sgtitle('variate target');

%% variate target position
ss_num = [3860 3859 3856];
for ssi = 1:length(ss_num) 
%     sstmp2(ssi) = SessionScan(ss_num(ssi));
end 

% just plot the raw data overlay eath other
fh1 = figure('Unit', 'inch', 'Position', [0 0 3 6]);
figure(fh1);
%sstmp = SessionScan(ss_list(fce_i, tar_i));
for fce_i = 3:-1:1
     sstmp1 = sstmp2(fce_i);
% fce_i = 1
% figure(fh);
    %axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
    trials_list = find([sstmp1.trials.outcome] == 1);
    
    for trial_i = trials_list
        axh(1) = subplot(4,1,1); hold on;       % posiitoin
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.x(2,:), ...
            'Color', col_type(fce_i+3,:));
        ylabel('m'); title('position');
        %chi = 3;    % only plot one muscle here
        axh(2) = subplot(4,1,2); hold on;       % force 
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.f(2,:), ...
            'Color', col_type(fce_i+3,:));
        ylabel('N'); title('force');
        axh(3) = subplot(4,1,3); hold on;       % ch 3: elbow flexor
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.emg(3,:), ...
            'Color', col_type(fce_i+3,:));    
        ylabel('EMG activity'); title('Elbow Flexor');
        axh(4) = subplot(4,1,4); hold on;       % % ch 4: elbow flexor
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.emg(4,:), ...
            'Color', col_type(fce_i+3,:));
        ylabel('EMG activity'); title('Elbow Extensor')
    end
    
linkaxes(axh, 'x');
% xlim([-3 2]); % sec
set(axh(1), 'YLim', [0.45 0.58])
xlim([-1 0.5]); % sec
% ylim([-1 5]);
set(axh(3), 'YLim', [0 40]);
set(axh(4), 'YLim', [0 20]);
end
sgtitle('variate force');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EMG analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Try to get the potent-null space EMG activities %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the first test, just use 1 pair of muscle on one session data
%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Get the tidied up data (already trialfied) 
% % % sstmp = SessionScan(3859);
% % % celltmp = sstmp.export_as_formatted(1);
% 1. Get the data from the force ramp 
%%
ts_fcrmp = 3;  
fh(1) = figure(); 
axh(1) = subplot(4,1,1); hold on;% position 
axh(2) = subplot(4,1,2); hold on; % force 
axh(3) = subplot(4,1,3); hold on; % emg-3
axh(4) = subplot(4,1,4); hold on; % emg-4
binrange = [0, 3];
y_fce = [];
X_emg = [];

for triali = 1:15 
    
    % find the first continuously ts_fcrmp in the data 
    trialtmp = celltmp{1,1,triali,2}; 
    tsidx = trialtmp.ts == find(ts_fcrmp); 
    trialtmp.t_shift = trialtmp.t - trialtmp.t(ts_fcrmp(1) + 200);
    plot(axh(1),trialtmp.t_shift,trialtmp.ts);
    % plot them (ifplot)
    plot(axh(2), trialtmp.t_shift, trialtmp.f(2,:));
    
    plot(axh(3), trialtmp.t_shift, smooth(trialtmp.emg(3,:), 200));
    
    plot(axh(4), trialtmp.t_shift, smooth(trialtmp.emg(4,:), 200)); 
    
    dat_idx = trialtmp.t_shift >= binrange(1) & trialtmp.t_shift <= binrange(2);
    y_fce = [y_fce, trialtmp.f(2,dat_idx)];
    X_emg = [X_emg, trialtmp.emg(3:4,dat_idx)];
    
end
linkaxes(axh, 'x');
xlim([0, 3]);

% plot the force and emg
fh(2) = figure();
axh(1) = subplot(2,1,1); 
scatter(X_emg(1,:), y_fce, 3); 
title('bicep with force'); 
ylabel('N');
axh(2) = subplot(2,1,2); 
scatter(X_emg(2,:), y_fce, 3); 
title('tricep with force'); 
ylabel('N');
xlabel('activity value');



% 2. regress the EMG data during the force ramp to the actual force 
% do regression
X = [ones(1, size(X_emg,2)); X_emg]';
y = y_fce';
mdl = fitlm(X_emg', y_fce');

% 3. Plot the predicted force and the actual force                         
% (Scott thesis, figure 3-13) 

y_fce_pred = (mdl.predict(X_emg'))';
figure(); hold on;
scatter(y_fce_pred, y_fce, 3);
plot(y_fce_pred, y_fce_pred, 'LineWidth', 3);
xlabel('pred fce'); 
ylabel('real fce'); 


% 4. plot the potent and null activities against the stiffness predictions
% (Scott thesis, figure 3-7)

% F_t = \beta_0 + W*EMG_t 
W = mdl.Coefficients.Estimate;
W = W(2:end)';
[U,S,V] = svd(W);

V_pot = V(1,:);
V_null= V(2:end,:);
EMG_pot = V_pot * X_emg;
EMG_null = V_null * X_emg;

subplot(1,2,1); plot(EMG_pot); 
subplot(1,2,2); plot(EMG_pot); 


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the second test, just use 1 pair of muscle on multiple sessions data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Get the tidied up data (already trialfied) 
data = cell(1,1,3,1,15,3);
ss_num = [3860 3859 3856];
for ssi = 1:length(ss_num) 
    sstmp2(ssi) = SessionScan(ss_num(ssi));
    celltmp = sstmp2(ssi).export_as_formatted();
    data(1,1,ssi,1,:,:) = reshape(celltmp(1,1,1:15,1:3),1,15,3);
end 

% 1. Get the data from the force ramp 

ts_fcrmp = 3;  
fh(1) = figure(); 
col_type = colormap('lines');
axh(1) = subplot(4,1,1); hold on;% position 
axh(2) = subplot(4,1,2); hold on; % force 
axh(3) = subplot(4,1,3); hold on; % emg-3
axh(4) = subplot(4,1,4); hold on; % emg-4
binrange = [0, 3];
y_fce = [];
X_emg = [];

for fce_i = 1:3
    celltmp = reshape(data(1,1,fce_i,1,:,:), 1,1,15,3);
for triali = 1:15 
    
    % find the first continuously ts_fcrmp in the data 
    trialtmp = celltmp{1,1,triali,2}; 
    tsidx = trialtmp.ts == find(ts_fcrmp); 
    trialtmp.t_shift = trialtmp.t - trialtmp.t(ts_fcrmp(1) + 200);
    plot(axh(1),trialtmp.t_shift,trialtmp.ts, 'Color', col_type(fce_i,:));
    % plot them (ifplot)
    plot(axh(2), trialtmp.t_shift, trialtmp.f(2,:), 'Color', col_type(fce_i,:));
    
    plot(axh(3), trialtmp.t_shift, smooth(trialtmp.emg(3,:), 200), 'Color', col_type(fce_i,:));
    
    plot(axh(4), trialtmp.t_shift, smooth(trialtmp.emg(4,:), 200), 'Color', col_type(fce_i,:)); 
    
    dat_idx = trialtmp.t_shift >= binrange(1) & trialtmp.t_shift <= binrange(2);
    y_fce = [y_fce, trialtmp.f(2,dat_idx)];
    X_emg = [X_emg, trialtmp.emg(3:4,dat_idx)];
    
end
end
linkaxes(axh, 'x');
xlim([0, 3]);

% plot the force and emg
fh(2) = figure();
axh(1) = subplot(2,1,1); 
scatter(X_emg(1,:), y_fce, 3); 
title('bicep with force'); 
ylabel('N');
axh(2) = subplot(2,1,2); 
scatter(X_emg(2,:), y_fce, 3); 
title('tricep with force'); 
ylabel('N');
xlabel('activity value');



% 2. regress the EMG data during the force ramp to the actual force 
% do regression
X = [ones(1, size(X_emg,2)); X_emg]';
y = y_fce';
mdl = fitlm(X_emg', y_fce');

% 3. Plot the predicted force and the actual force                         
% (Scott thesis, figure 3-13) 

y_fce_pred = (mdl.predict(X_emg'))';
figure(); hold on;
scatter(y_fce_pred, y_fce, 3);
plot(y_fce_pred, y_fce_pred, 'LineWidth', 3);
xlabel('pred fce'); 
ylabel('real fce'); 


% 4. plot the potent and null activities against the stiffness predictions
% (Scott thesis, figure 3-7)

% F_t = \beta_0 + W*EMG_t 
W = mdl.Coefficients.Estimate;
W = W(2:end)';
[U,S,V] = svd(W);

V_pot = V(1,:);
V_null= V(2:end,:);
EMG_pot = V_pot * X_emg;
EMG_null = V_null * X_emg;

subplot(1,2,1); plot(EMG_pot); 
subplot(1,2,2); plot(EMG_pot);  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  After collected a new batch of data of 8 muscles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variate target position
ss_num = [...
    3866    3865    3864;   %   15N 
    3867    3869    3868;   %   20N
    3871    3870    3872];  %   25N
ss_num = ss_num(:);
for ssi = 1:size(ss_num,1) 
%     for tar_i = 1:size(ss_num,2)
%         sstmp2(ssi) = SessionScan(ss_num(fce_i,tar_i));
    sstmp2(ssi) = SessionScan(ss_num(ssi));
%     end
end 

%% 
% export and save 
data = cell(1,1,3,1,15,3);
ss_num = [...
    3866    3865    3864;   %   15N 
    3867    3869    3868;   %   20N
    3871    3870    3872];  %   25N
for fce_i = 1:size(ss_num,1) 
    for tar_i = 1:size(ss_num,2)
        sstmp = SessionScan(ss_num(fce_i,tar_i));
        celltmp = sstmp.export_as_formatted();
        data(1,1,fce_i,tar_i,:,:) = reshape(celltmp(1,1,1:15,1:3),1,15,3);
    end 
end
save('data/processedData/ss3864_3872.mat', 'data');%% just plot the raw data overlay eath other

%%
fh1 = figure('Unit', 'inch', 'Position', [0 0 3 6]);
figure(fh1);
%sstmp = SessionScan(ss_list(fce_i, tar_i));
emg_pair = [5 6];
fce_i = 2
for tar_i = 1:3
% for fce_i = 1:3
     sstmp1 = sstmp2((tar_i-1)*3+fce_i);
% fce_i = 1
% figure(fh);
    %axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
    trials_list = find([sstmp1.trials.outcome] == 1);
    
    for trial_i = trials_list
        axh(1) = subplot(4,1,1); hold on;       % posiitoin
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.x(2,:), ...
            'Color', col_type(fce_i+3+tar_i,:));
        ylabel('m'); title('position');
        %chi = 3;    % only plot one muscle here
        axh(2) = subplot(4,1,2); hold on;       % force 
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.f(2,:), ...
            'Color', col_type(fce_i+3+tar_i,:));
        ylabel('N'); title('force');
        axh(3) = subplot(4,1,3); hold on;       % ch 3: elbow flexor
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.emg(emg_pair(1),:), ...
            'Color', col_type(fce_i+3+tar_i,:));    
        ylabel('EMG activity'); title('Elbow Flexor');
        axh(4) = subplot(4,1,4); hold on;       % % ch 4: elbow flexor
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.emg(emg_pair(2),:), ...
            'Color', col_type(fce_i+3+tar_i,:));
        ylabel('EMG activity'); title('Elbow Extensor')
    end
    
linkaxes(axh, 'x');
% xlim([-3 2]); % sec
set(axh(1), 'YLim', [0.45 0.58])
xlim([-1 0.5]); % sec
% ylim([-1 5]);
set(axh(3), 'YLim', [0 0.2]);
set(axh(4), 'YLim', [0 0.2]);
end
sgtitle('variate force');


%% just plot the raw data overlay eath other, 8 muscles and 3*3 conditiosn
fh1 = figure('Unit', 'inch', 'Position', [0 0 3 6]);
figure(fh1);
%sstmp = SessionScan(ss_list(fce_i, tar_i));
% % emg_pair = [1 2 3 4 5 6 7 8];
% emg_pair = [1 2 3 4 5 6];
% emg_pair = [3 4];
emg_pair = [1 2];
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                  'elbow flexor', 'elbow extensor', ...
                  'anterior deltaloid', 'posterior deltaloid', ...
                  'pectoralis', 'Trapezius'}; 
cols = 3;
rows = 2 + length(emg_pair);
for fce_i = 1:3
for tar_i = 1:3
% for fce_i = 1:3
     sstmp1 = sstmp2((tar_i-1)*3+fce_i);
% fce_i = 1
% figure(fh);
    %axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
    trials_list = find([sstmp1.trials.outcome] == 1);
    
    for trial_i = trials_list
        axh(1, fce_i) = subplot(rows,cols,fce_i); hold on;       % posiitoin
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.x(2,:), ...
            'Color', col_type(tar_i,:));
        ylabel('m'); title('position');
        %chi = 3;    % only plot one muscle here
        axh(2, fce_i) = subplot(rows,cols,fce_i+3); hold on;       % force 
        plot(sstmp1.trials(trial_i).data.t_shift, ...
            sstmp1.trials(trial_i).data.f(2,:), ...
            'Color', col_type(tar_i,:));
        ylabel('N'); title('force');
        for muscle_i = 1:length(emg_pair)
            axh(2+muscle_i, fce_i) = subplot(rows,cols,fce_i + (muscle_i+2-1)*3); hold on; grid on       % ch 3: elbow flexor
            plot(sstmp1.trials(trial_i).data.t_shift, ...
                sstmp1.trials(trial_i).data.emg(emg_pair(muscle_i),:), ...
                'Color', col_type(tar_i,:));    
            ylabel('EMG activity'); title(emg_pair_label{emg_pair(muscle_i)});
        end
    end
    
linkaxes(axh(:), 'x');
% xlim([-3 2]); % sec
set(axh(1, fce_i), 'YLim', [0.45 0.58])
xlim([-1 0.5]); % sec
% ylim([-1 5]);
linkaxes(axh(3:end,fce_i), 'y');
set(axh(3, fce_i), 'YLim', [0 0.2]);

end
end
sgtitle('variate force');


%% just plot the raw data overlay eath other, use tidied up format
fh1 = figure('Unit', 'inch', 'Position', [0 0 3 6]);
figure(fh1);
%sstmp = SessionScan(ss_list(fce_i, tar_i));
% % emg_pair = [1 2 3 4 5 6 7 8];
% emg_pair = [1 2 3 4 5 6];
% emg_pair = [3 4];
emg_pair = [1 2];
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                  'elbow flexor', 'elbow extensor', ...
                  'anterior deltaloid', 'posterior deltaloid', ...
                  'pectoralis', 'Trapezius'}; 
cols = 3;
rows = 2 + length(emg_pair);
for fce_i = 1:3
for tar_i = 1:3
% for fce_i = 1:3
     datatmp = reshape(data(1,1,fce_i,tar_i,:,1:2),1,30);
     
% fce_i = 1
% figure(fh);
    %axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
    trials_list = 1:30;
    
    for trial_i = trials_list
        axh(1, fce_i) = subplot(rows,cols,fce_i); hold on;       % posiitoin
        if isempty(datatmp{trial_i})
            continue;
        end
        idx_tshift = find((datatmp{trial_i}.ts == 5));
        t_shift = datatmp{trial_i}.t - datatmp{trial_i}.t(idx_tshift(1));
        plot(t_shift, ...
            datatmp{trial_i}.x(2,:), ...
            'Color', col_type(tar_i,:));
        ylabel('m'); title('position');
        %chi = 3;    % only plot one muscle here
        axh(2, fce_i) = subplot(rows,cols,fce_i+3); hold on;       % force 
        plot(t_shift, ...
            datatmp{trial_i}.f(2,:), ...
            'Color', col_type(tar_i,:));
        ylabel('N'); title('force');
        for muscle_i = 1:length(emg_pair)
            axh(2+muscle_i, fce_i) = subplot(rows,cols,fce_i + (muscle_i+2-1)*3); hold on; grid on       % ch 3: elbow flexor
            plot(t_shift, ...
                datatmp{trial_i}.emg(emg_pair(muscle_i),:), ...
                'Color', col_type(tar_i,:));    
            ylabel('EMG activity'); title(emg_pair_label{emg_pair(muscle_i)});
        end
    end
    
linkaxes(axh(:), 'x');
% xlim([-3 2]); % sec
set(axh(1, fce_i), 'YLim', [0.45 0.58])
xlim([-1 0.5]); % sec
% ylim([-1 5]);
linkaxes(axh(3:end,fce_i), 'y');
set(axh(3, fce_i), 'YLim', [0 0.2]);

end
end
sgtitle('variate force');

%% calculate the average and std across trials, plot them 

fh1 = figure('Unit', 'inch', 'Position', [0 0 3 6]);
figure(fh1);
t_range = [-3 1];
%sstmp = SessionScan(ss_list(fce_i, tar_i));
emg_pair = [1 2 3 4 5 6 7 8];
% emg_pair = [1 2 3 4 5 6];
% emg_pair = [3 4];
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                  'elbow flexor', 'elbow extensor', ...
                  'anterior deltaloid', 'posterior deltaloid', ...
                  'pectoralis', 'Trapezius'}; 
cols = 3;
rows = 2 + length(emg_pair);
for fce_i = 1:3
for tar_i = 1:3
% for fce_i = 1:3
     sstmp1 = sstmp2((tar_i-1)*3+fce_i);
% fce_i = 1
% figure(fh);
    %axh(fce_i) = subplot(8,3,(chi-1)*3+fce_i); hold on; grid on;
    trials_list = find([sstmp1.trials.outcome] == 1);
    
    % make enough space for data 
    dat.row = length(trials_list);
    dat.col = sum(sstmp1.trials(trials_list(1)).data.t_shift > t_range(1) & ...
                    sstmp1.trials(trials_list(1)).data.t_shift < t_range(2));
%     dat.t   = sstmp1.trials(trials_list(1)).data.t_shift(...
%                     sstmp1.trials(trials_list(1)).data.t_shift > t_range(1) & ... % the index of time
%                     sstmp1.trials(trials_list(1)).data.t_shift < t_range(2));
    dat.t = t_range(1):1/500:t_range(2);
    dat.t = dat.t(2:end-1);
    dat.pos = nan(dat.row, dat.col);
    dat.fce = nan(dat.row, dat.col);
    dat.emg = nan(8, dat.row, dat.col);
    
    trial_idx = 0;
    for trial_i = trials_list 
        % get the index
        trial_idx = trial_idx + 1;
        % stack the data into matrices
        index_t = sstmp1.trials(trial_i).data.t_shift > t_range(1) & ...
                  sstmp1.trials(trial_i).data.t_shift < t_range(2);
        dat.pos(trial_idx,1:sum(index_t)) = sstmp1.trials(trial_i).data.x(2,index_t);
        dat.fce(trial_idx,1:sum(index_t)) = sstmp1.trials(trial_i).data.f(2,index_t);
        dat.emg(:,trial_idx,1:sum(index_t))=sstmp1.trials(trial_i).data.emg(:,index_t);
    end
    
    
    % plot the mean and average of it 
    
%     for trial_i = trials_list
        axh(1, fce_i) = subplot(rows,cols,fce_i); hold on;       % posiitoin
        plot(dat.t, nanmean(dat.pos), ...
            'Color', col_type(fce_i+3+tar_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
        tmp1 = nanmean(dat.pos)+nanstd(dat.pos); tmp2 =  nanmean(dat.pos)-nanstd(dat.pos);
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type(fce_i+3+tar_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        ylabel('m'); title('position');
        %chi = 3;    % only plot one muscle here
        axh(2, fce_i) = subplot(rows,cols,fce_i+3); hold on;       % force 
        plot(dat.t, nanmean(dat.fce), ...
            'Color', col_type(fce_i+3+tar_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
        tmp1 = nanmean(dat.fce)+nanstd(dat.fce); tmp2 =  nanmean(dat.fce)-nanstd(dat.fce);
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type(fce_i+3+tar_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        ylabel('N'); title('force');
        for muscle_i = 1:length(emg_pair)
            axh(2+muscle_i, fce_i) = subplot(rows,cols,fce_i + (muscle_i+2-1)*3); hold on;       % ch 3: elbow flexor
            plot(dat.t, nanmean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3))), ...
                'Color', col_type(fce_i+3+tar_i,:), ...
                'LineWidth', 2);    
            ylabel('EMG activity'); title(emg_pair_label{emg_pair(muscle_i)});  
            
            patch_x = [dat.t, dat.t(end:-1:1)];
            tmp1 = nanmean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)))+ ...
                nanstd(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)));
            tmp2 = nanmean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)))- ...
                nanstd(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)));
            patch_y = [tmp1, tmp2(end:-1:1)];
            patch(patch_x, patch_y, ...
                col_type(fce_i+3+tar_i,:), ...
                'FaceAlpha', 0.3, ...
                'EdgeColor', 'none');
            
            xline(0.07); 
            grid on;
        end
%     end
    
linkaxes(axh(:), 'x');
% xlim([-3 2]); % sec
set(axh(1, fce_i), 'YLim', [0.45 0.58])
xlim([-1 0.5]); % sec
% ylim([-1 5]);
linkaxes(axh(3:end,:), 'y');
set(axh(3, fce_i), 'YLim', [0 0.1]);

end
end
sgtitle('variate force');

%% Do the same plot, but used the tidied up data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A good plot here showing how muscles react during the movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh1 = figure();
col_type = colormap('lines');
figure(fh1);
t_range = [-3 1];
%sstmp = SessionScan(ss_list(fce_i, tar_i));
emg_pair = [1 2 3 4 5 6 7 8];
% emg_pair = [1 2 3 4 5 6];
% emg_pair = [3 4];
emg_pair_label = {'wrist flexor', 'wrist extensor', ...
                  'elbow flexor', 'elbow extensor', ...
                  'anterior deltaloid', 'posterior deltaloid', ...
                  'pectoralis', 'Trapezius'}; 
emg_pair_Label = {'FCR', 'EXD', 'BIC', 'TIC', 'DTA', 'DTP', 'PEC', 'TPZ'};
cols = 3;
y_mag = [0.55, 30, .15, .15, .1, .1, .2, .05, 0.03, .05];
rows = 2 + length(emg_pair);
load('data/processedData/ss3864_3872.mat', 'data');%
force_list = [15 20 25];
for fce_i = 1:3
for tar_i = 1:3

    datass = reshape(data(1,1,fce_i, tar_i,:,1:2), 1, 30);

    trials_list = 1:30;
    % get t_shift
    idx_tshift = find((datass{1}.ts == 5));
    t_shift = datass{1}.t - datass{1}.t(idx_tshift(1));
    % make enough space for data 
    dat.row = length(trials_list);
    dat.col = sum(t_shift > t_range(1) & ...
                    t_shift < t_range(2));
    dat.t = t_range(1):1/500:t_range(2);
    dat.t = dat.t(2:end-1);
    dat.pos = nan(dat.row, dat.col);
    dat.fce = nan(dat.row, dat.col);
    dat.emg = nan(8, dat.row, dat.col);
    
    trial_idx = 0;
    for trial_i = trials_list 
        % get t-shift
        idx_tshift = find((datass{trial_i}.ts == 5));
        t_shift = datass{trial_i}.t - datass{trial_i}.t(idx_tshift(1));
        % get the index
        trial_idx = trial_idx + 1;
        % stack the data into matrices
        index_t = t_shift > t_range(1) & t_shift < t_range(2);
        dat.pos(trial_idx,1:sum(index_t)) = datass{trial_i}.x(2,index_t);
        dat.fce(trial_idx,1:sum(index_t)) = datass{trial_i}.f(2,index_t);
        dat.emg(:,trial_idx,1:sum(index_t))=datass{trial_i}.emg(:,index_t);
    end
    
    
    % plot the mean and average of it 
    
%     for trial_i = trials_list
        axh(1, fce_i) = subplot(rows,cols,fce_i); hold on;       % posiitoin
        plot(dat.t, nanmean(dat.pos), ...
            'Color', col_type(tar_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
        tmp1 = nanmean(dat.pos)+nanstd(dat.pos); tmp2 =  nanmean(dat.pos)-nanstd(dat.pos);
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type(tar_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        if (fce_i == 1)
            ylabel('m'); 
        else 
            set(gca, 'yticklabels', {});
        end
        set(gca, 'xticklabels', {});
        title(['force' num2str(force_list(fce_i)) 'N']);
        %chi = 3;    % only plot one muscle here
        axh(2, fce_i) = subplot(rows,cols,fce_i+3); hold on;       % force 
        plot(dat.t, nanmean(dat.fce), ...
            'Color', col_type(tar_i,:), ...
            'LineWidth', 2);
        patch_x = [dat.t, dat.t(end:-1:1)];
        tmp1 = nanmean(dat.fce)+nanstd(dat.fce); tmp2 =  nanmean(dat.fce)-nanstd(dat.fce);
        patch_y = [tmp1, tmp2(end:-1:1)];
        patch(patch_x, patch_y, ...
            col_type(tar_i,:), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none');
        if (fce_i == 1)
            ylabel('N');     
        else 
            set(gca, 'yticklabels', {});
        end
        set(gca, 'xticklabels', {});
        %title('force');
        for muscle_i = 1:length(emg_pair)
            axh(2+muscle_i, fce_i) = subplot(rows,cols,fce_i + (muscle_i+2-1)*3); hold on;       % ch 3: elbow flexor
            plot(dat.t, nanmean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3))), ...
                'Color', col_type(tar_i,:), ...
                'LineWidth', 2);    
            if (fce_i == 1)
%                 ylabel('EMG'); %title(emg_pair_label{emg_pair(muscle_i)}); 
                ylabel(emg_pair_Label{muscle_i});
            else 
                set(gca, 'yticklabels', {});
            end
            if (muscle_i == length(emg_pair))
                xlabel('time (s)');
            else
                set(gca, 'xticklabels', {});
            end
            
            patch_x = [dat.t, dat.t(end:-1:1)];
            tmp1 = nanmean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)))+ ...
                nanstd(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)));
            tmp2 = nanmean(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)))- ...
                nanstd(reshape(dat.emg(muscle_i,:,:), size(dat.emg,2), size(dat.emg,3)));
            patch_y = [tmp1, tmp2(end:-1:1)];
            patch(patch_x, patch_y, ...
                col_type(tar_i,:), ...
                'FaceAlpha', 0.3, ...
                'EdgeColor', 'none');
            
            xline(0, '-');
            xline(0.07, '--'); 
            grid on;
        end
%     end
    


end
end
legend(axh(1,1), '2.5cm', '5cm', '7.5cm');
linkaxes(axh(:), 'x');
% xlim([-3 2]); % sec
set(axh(1, fce_i), 'YLim', [0.45 0.58])
xlim([-0.5 0.5]); % sec
% ylim([-1 5]);
% linkaxes(axh(3:end,:), 'y');
% set(axh(3, fce_i), 'YLim', [0 0.1]);
% reset the panel positions 
for pci = 1:3
    for pri = 1:10
        ri = pri;
        ci = pci;
        set(axh(pri, pci), 'position', [0.1+0.283*(ci-1), 0.05 + 0.09*((10-ri+1)-1), 0.28, 0.08]);
        if (pri >= 3)
            set(axh(pri,pci), 'YLim', [0, y_mag(pri)]);
        end
        linkaxes(axh(pri,:), 'y');
    end
end

sgtitle('EMG data ');


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Using SVD to find potent and null EMG activities 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Get the data from the force ramp 

ts_fcrmp = 3; 
fh(1) = figure(); 
col_type = colormap('lines');
axh(1) = subplot(4,1,1); hold on;% position 
axh(2) = subplot(4,1,2); hold on; % force 
axh(3) = subplot(4,1,3); hold on; % emg-3
axh(4) = subplot(4,1,4); hold on; % emg-4
binrange = [0, 3];
y_fce = [];
X_emg = [];

for fce_i = 1:3
    celltmp = reshape(data(1,1,fce_i,1,:,:), 1,1,15,3);
for triali = 1:15 
    
    % find the first continuously ts_fcrmp in the data 
    trialtmp = celltmp{1,1,triali,2}; 
    tsidx = trialtmp.ts == find(ts_fcrmp); 
    trialtmp.t_shift = trialtmp.t - trialtmp.t(ts_fcrmp(1) + 200);
    plot(axh(1),trialtmp.t_shift,trialtmp.ts, 'Color', col_type(fce_i,:));
    % plot them (ifplot)
    plot(axh(2), trialtmp.t_shift, trialtmp.f(2,:), 'Color', col_type(fce_i,:));
    
    plot(axh(3), trialtmp.t_shift, smooth(trialtmp.emg(3,:), 200), 'Color', col_type(fce_i,:));
    
    plot(axh(4), trialtmp.t_shift, smooth(trialtmp.emg(4,:), 200), 'Color', col_type(fce_i,:)); 
    
    dat_idx = trialtmp.t_shift >= binrange(1) & trialtmp.t_shift <= binrange(2);
    y_fce = [y_fce, trialtmp.f(2,dat_idx)];
    X_emg = [X_emg, trialtmp.emg(1:8,dat_idx)];
    
end
end
X_emg = X_emg(:,1:4.4e4);
y_fce = y_fce(1:4.4e4);
linkaxes(axh, 'x');
xlim([0, 3]);

% plot the force and emg
fh(2) = figure();
axh(1) = subplot(4,2,1); 
scatter(X_emg(1,:), y_fce, 3); 
title('FCR & f'); 
ylabel('N');
axh(2) = subplot(4,2,2); 
scatter(X_emg(2,:), y_fce, 3); 
title('EXD & F'); 
ylabel('N');
xlabel('activity value');
axh(1) = subplot(4,2,3); 
scatter(X_emg(3,:), y_fce, 3); 
title('bicep with force'); 
ylabel('N');
axh(2) = subplot(4,2,4); 
scatter(X_emg(4,:), y_fce, 3); 
title('tricep with force'); 
ylabel('N');
xlabel('activity value');
axh(1) = subplot(4,2,5); 
scatter(X_emg(5,:), y_fce, 3); 
title('DTA & f'); 
ylabel('N');
axh(2) = subplot(4,2,6); 
scatter(X_emg(6,:), y_fce, 3); 
title('DTP & f'); 
ylabel('N');
xlabel('activity value');
axh(1) = subplot(4,2,7); 
scatter(X_emg(7,:), y_fce, 3); 
title('PEC & f'); 
ylabel('N');
axh(2) = subplot(4,2,8); 
scatter(X_emg(8,:), y_fce, 3); 
title('tricep with force'); 
ylabel('N');
xlabel('TPZ & f');



% 2. regress the EMG data during the force ramp to the actual force 
% do regression
X = [ones(1, size(X_emg,2)); X_emg]';
y = y_fce';
mdl = fitlm(X_emg', y_fce');

% 3. Plot the predicted force and the actual force                         
% (Scott thesis, figure 3-13) 

%%% figure of predicted force and real force
y_fce_pred = (mdl.predict(X_emg'))';
figure(); hold on;
scatter(y_fce_pred, y_fce, 3);
plot(y_fce_pred, y_fce_pred, 'LineWidth', 3);
legend('real', 'predicted');
title('force regression result');
xlabel('pred fce'); 
ylabel('real fce'); 

%%% figure of predicted force and real force, across time
figure(); hold on;
plot(1:length(y_fce_pred), y_fce_pred, 'r');
plot(1:length(y_fce), y_fce, 'b'); 
legend('predicted', 'actual');
title('predicted/real vlue across time');


% 4. plot the potent and null activities against the stiffness predictions
% (Scott thesis, figure 3-7)

% F_t = \beta_0 + W*EMG_t 
W = mdl.Coefficients.Estimate;
mdl_coef = mdl.Coefficients.Estimate;
W = W(2:end)';
[U,S,V] = svd(W);

V_pot = V(:,1);
V_null= V(:,2:8);
EMG_pot = S(1)*V_pot' * X_emg;
EMG_null = S(1)*V_null' * X_emg;    % keep it the same scale with EMG_pot

[coeff, score, latent] = pca((V_null'*X_emg)');


subplot(1,2,1); plot(EMG_pot'); 
subplot(1,2,2); plot(EMG_null'); 

figure(); % show how the PCA results give us
scores_EMG = (V_null' * X_emg)' *coeff;  
scores_EMG = scores_EMG';
plot(1:size(scores_EMG,2), scores_EMG);
legend('1', '2', '3', '4', '5', '6', '7');
title('PC scores of the output-null values');

% only plot the principle projection of V_null
subplot(2,1,1);
hold on; 
plot(1:length(EMG_pot), EMG_pot + mdl_coef(1));
plot(1:length(EMG_pot), y_fce); 
title('net force with potent EMG');

subplot(2,1,2); 
hold on; plot(1:length(EMG_pot), S(1)*scores_EMG(1,:));
plot(1:length(EMG_pot), y_fce); 
title('net force with null Principal EMG');

% get the value average of the EMG_null priciple 
EMG_np = zeros(3,3,15); 
EMG_pot= zeros(3,3,15);
for fce_i = 1:3
    for tar_i = 1:3
        for trial_i = 1:15
            idx = (data{1,1,fce_i,tar_i,trial_i,2}.ts == 4); % hold 
            X_emg = data{1,1,fce_i,tar_i,trial_i,2}.emg(:,idx);
            scores_EMG = (V_null' * X_emg)' *coeff;  
            scores_EMG = scores_EMG';
            emg_np = S(1)*scores_EMG(1,:); 
            emg_pot= S(1)*V_pot' * X_emg;
            EMG_np(fce_i,tar_i,trial_i) = mean(emg_np);
            EMG_pot(fce_i,tar_i,trial_i) = mean(emg_pot);
        end
    end
end


%% compare the null-principle-EMG with the estimated stiffness
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3864_3872.mat', 'data');
data_human = reshape(data(1,1,:,:,:,:),1,3,3,15,3);
dexSubject = 1; % [first subject]
dexForce = 1:3; % [15N, 20N, 25N]
dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]

depMeasures_TestHuman = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
depMeasures_TestHuman.k_nonNan


%% connect it with EMG_np
emg_np = reshape(EMG_np, 1,3*3*15);
emg_pot= reshape(EMG_pot,1,3*3*15);
khat_pulse = reshape(depMeasures_TestHuman.k_hat_pulse, 1, 3*3*15);
fh = figure('unit', 'inch', 'position', [0 0 3 3]); hold on; 
scatter(khat_pulse, emg_np, 4); 
scatter(khat_pulse, emg_pot, 4, 'r');
xlabel('pulse K hat (N/m)');
ylabel('emg null estimation');
title('pulse stiffness estimation vs emg null');  
% null
lm1 = fitlm(khat_pulse,emg_np);
lm1.Coefficients
coeff1 = lm1.Coefficients.Estimate;
% potent
lm2 = fitlm(khat_pulse,emg_pot);
lm2.Coefficients
coeff2 = lm2.Coefficients.Estimate;

x = [1,1; 200, 1600];
y = x' * coeff1; 
axh(1) = plot(x(2,:)', y, 'LineWidth', 2);
y = x' * coeff2;
axh(2) = plot(x(2,:)', y, 'LineWidth', 2);

legend([axh(1), axh(2)], 'output-null', 'output-potent');

% saveas('')