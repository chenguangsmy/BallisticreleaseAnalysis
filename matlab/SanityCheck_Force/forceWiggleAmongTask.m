% this script is trying to see if the force wiggles among different task
% conditoins.
% The force wiggle might come from the subject, rather than the actual
% force transducer noise

load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/' ...
    'data/processedData/ss4310_4356.mat']);

t_range = [-0.5 0];
fce_list = [15 20 25];
dist_list = [2.5 5.0 7.5];


%% plot the force
fh = figure(); 
col_type = colormap('lines');
close(fh);
for subj_i = 1%1:4
    for dir_i = 1:4
        fh = figure('name', ['subj ' num2str(subj_i) 'direction ' num2str(dir_i)]);
        clear axh;
        for fce_i = 1:3
            for dist_i = 1:3
                axh(fce_i, dist_i) = subplot(3,3,(fce_i-1)*3+dist_i); hold on;
                for trial_i = 1:9
                    dat_trial = data{subj_i,dir_i,fce_i,dist_i,trial_i,1};
                    t_idx = dat_trial.t>t_range(1) & dat_trial.t<t_range(2);
                    f_offset = mean(vecnorm(dat_trial.f(1:3,t_idx)));
                    plot(dat_trial.t, vecnorm(dat_trial.f(1:3,:)) - f_offset, ...
                        'linewidth', 2, 'color', col_type(trial_i,:));

                    title([num2str(fce_list(fce_i)) 'N ' num2str(dist_list(dist_i)) ...
                        'cm']);
                    grid on;
                    xlabel('time upon release (s)');
                    ylabel('force (N)');
                    
                end
                xlim([-1 -0.001]);
            end
        end
        linkaxes(axh, 'xy');
        ylim([-3 3]);
        sgtitle(fh, ['subj ' num2str(subj_i) ...
            'direction ' num2str(dir_i)]);
    end
end

% 
% figure(1); figure(2); figure(3); figure(4);         % subject 1
% figure(5); figure(6); figure(7); figure(8);         % subject 2
% figure(9); figure(10); figure(11); figure(12);      % subject 3
% figure(13); figure(14); figure(15); figure(16);     % subject 4

% % %
% p_mat = []; 
% for subj_i = 1:4
%     for dir_i = 1:4
%         for fce_i = 1:3
%             for dist_i = 1:3
%                 for trial_i = 1:9
%                     dat_trial = data{subj_i,dir_i,fce_i,dist_i,trial_i,1};
%                     t_idx = dat_trial.t>t_range(1) & dat_trial.t<t_range(2);
%                     f_offset = mean(vecnorm(dat_trial.f(1:3,t_idx)));
%                     dat_tmp = vecnorm(dat_trial.f(1:3,:)) - f_offset;
%                     [p,f] = pspectrum(dat_tmp, 500);
%                     p_mat = [p_mat, p];
%                 end
%             end
%         end
%     end
% end
% p_avg = min(p_mat,[],2);

%%
for subj_i = 1%1:4
    for dir_i = 1:4
        fh = figure('name', ['subj ' num2str(subj_i) 'direction ' num2str(dir_i)]);
        clear axh;
        for fce_i = 1:3
            for dist_i = 1:3
                axh(fce_i, dist_i) = subplot(3,3,(fce_i-1)*3+dist_i); hold on;
                for trial_i = 1:9
                    dat_trial = data{subj_i,dir_i,fce_i,dist_i,trial_i,1};
                    t_idx = dat_trial.t>-1 & dat_trial.t<0;
                    f_offset = mean(vecnorm(dat_trial.f(1:3,t_idx)));
%                     plot(dat_trial.t, vecnorm(dat_trial.f(1:3,:)) - f_offset);
                    dat_tmp = vecnorm(dat_trial.f(1:3,t_idx)) - f_offset;
%                     pwelch(dat_tmp, 500)
%                     [p,f] = pspectrum(dat_tmp, 500);
%                     [p,f] = pwelch(dat_tmp, 500);
%                     [p,f] = pwelch(dat_tmp, 100, 90, 1:0.2:500, 500);
% plot(dat_tmp)
                    pwelch(dat_tmp(1:200), hamming(100), 90, 1:0.2:250, 500)
                    % pdb = pow2db(p);%./flg;
%                     flg = log10(f);
%                     pdb = pow2db(p)./f;

%                     plot(flg, pdb, 'linewidth', 2,'color', col_type(trial_i,:)); 
%                     plot(log10(f), 10*log10(p), 'linewidth', 2,'color', col_type(trial_i,:)); 
%                     plot(log10(f), 10*log10(p./f), 'linewidth', 2,'color', col_type(trial_i,:)); 
                    title([num2str(fce_list(fce_i)) 'N ' num2str(dist_list(dist_i)) ...
                        'cm']);
                    grid on; 
%                     ylim([0 0.2]);
%                     xlabel('log(Frequency) (log(Hz))');
%                     xlabel('Frequency (Hz))');
%                     xlabel('Frequency (Hz))');
%                     ylabel('intensity (?)');
%                     ylabel('power');
                    xlabel('log(Hz)');
                    ylabel('dB/Hz');
%                     ylabel('dB');
                end
%                 xlim([-1 0.01]);
%                 xlim(log([0.2 4])); 
%                 xlim([1 15]);
            end
        end
        linkaxes(axh, 'xy');
%         ylim([-3 3]);
%         ylim([-1 5]);
%         ylim([0 0.2]);
        sgtitle(fh, ['subj ' num2str(subj_i) ...
            'direction ' num2str(dir_i)]);
    end
end

% figure(1); figure(2); figure(3); figure(4);         % subject 1
% figure(5); figure(6); figure(7); figure(8);         % subject 2
% figure(9); figure(10); figure(11); figure(12);      % subject 3
% figure(13); figure(14); figure(15); figure(16);     % subject 4

%% does the shaking contains the movement in x and ox? 
fh = figure(); 
col_type = colormap('lines');
t_range = [-2 0]; 
close(fh);
for subj_i = 1%1:4
    for dir_i = 2
        fh = figure('name', ['subj ' num2str(subj_i) 'direction ' num2str(dir_i)]);
        clear axh;
        for fce_i = 3
            for dist_i = 1
                
                for trial_i = 8%1:9
                    figure(); 
                    dat_trial = data{subj_i,dir_i,fce_i,dist_i,trial_i,1};
                    t_idx = dat_trial.t>t_range(1) & dat_trial.t<t_range(2);
                    f_offset = mean(vecnorm(dat_trial.f(1:3,t_idx)));
                    x_offset = mean(vecnorm(dat_trial.x(1:3,t_idx)));
                    ox_offset = mean(vecnorm(dat_trial.ox(1:3,t_idx,1)));
                    axh(1) = subplot(3,1,1); 
                    plot(dat_trial.t, vecnorm(dat_trial.f(1:3,:)) - f_offset, ...
                        'linewidth', 2, 'color', col_type(trial_i,:));
                    title('Force');

%                     title([num2str(fce_list(fce_i)) 'N ' num2str(dist_list(dist_i)) ...
%                         'cm']);
                    grid on;
                    xlabel('time upon release (s)');
                    ylabel('force (N)');

                    axh(2) = subplot(3,1,2); 
                    plot(dat_trial.t, vecnorm(dat_trial.x(1:3,:)) - x_offset, ...
                        'linewidth', 2, 'color', col_type(trial_i,:));
                    grid on;
                    xlabel('time upon release (s)');
                    ylabel('displacement (m)');
                    title('displacement (WAM, m)');

                    axh(3) = subplot(3,1,3); 
                    plot(dat_trial.t, vecnorm(dat_trial.ox(1:3,:,1)) - ox_offset, ...
                        'linewidth', 2, 'color', col_type(trial_i,:));
                    grid on;
                    xlabel('time upon release (s)');
                    ylabel('displacement (m)');
                    title('displacement (OPT, m)');

                    
                    
                end
                xlim([-1 -0.001]);
            end
        end
%         linkaxes(axh, 'xy');
%         ylim([-3 3]);
        linkaxes(axh, 'x');
        sgtitle(fh, ['subj ' num2str(subj_i) ...
            'direction ' num2str(dir_i)]);
    end
end
%% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do some force task among different conditions::
% ss4360: force ramp with seeing 2.5cm target. +-1N 2 trials with simulated arm shaking at 15N (big)
% ss4361: force ramp without seeing target. +-1N.
% ss4362: force ramp without seeing target. +-2N.
% ss4363: force ramp with seeing 2.5cm target. +-2N
% ss4364: force ramp with seeing 2.5cm target. +-2N, longer holding time...
clear; close all; 
data = cell(2,1,15,5);

sstmp(1) = SessionScan(4360);
sstmp(2) = SessionScan(4361);
sstmp(3) = SessionScan(4362);
sstmp(4) = SessionScan(4363);
sstmp(5) = SessionScan(4364);

fh(1) = sstmp(1).plotTrialfyForceh();
fh(2) = sstmp(2).plotTrialfyForceh();
fh(3) = sstmp(3).plotTrialfyForceh();
fh(4) = sstmp(4).plotTrialfyForceh();
fh(5) = sstmp(5).plotTrialfyForceh();
% set(fh(1), 'Title', {'2.5cm target ±1N'});
% title(fh(2), 'No target ±1N')
% title(fh(3), 'No target ±2N')
% title(fh(4), '2.5cm target ±2N')
% title(fh(5), '2.5cm target ±2N, long duration');

% plot out the spectrum of each trial/condition

dat1 = sstmp(1).export_as_formatted(1);
dat2 = sstmp(2).export_as_formatted(1);
dat3 = sstmp(3).export_as_formatted(1);
dat4 = sstmp(4).export_as_formatted(1);
dat5 = sstmp(5).export_as_formatted(1);
% data = cell(1,1,15,5);
data(1,1,:,1) = dat1(1,1,1:15,1);
data(1,1,:,2) = dat2(1,1,1:15,1);
data(1,1,:,3) = dat3(1,1,1:15,1);
data(1,1,:,4) = dat4(1,1,1:15,1);
data(1,1,:,5) = dat5(1,1,1:15,1);
dat_ref = dat1{1,1,17,1};

%%%%%%%% let subject HA do the same thing
% ss4366: force ramp with seeing 2.5cm target. +-1N 2 trials with simulated arm shaking at 15N (big)
% ss4367: force ramp without seeing target. +-1N.
% ss4368: force ramp without seeing target. +-2N.
% ss4369: force ramp with seeing 2.5cm target. +-2N
% ss4371: force ramp with seeing 2.5cm target. +-2N, longer holding time...

sstmp(1) = SessionScan(4366);
sstmp(2) = SessionScan(4367);
sstmp(3) = SessionScan(4368);
sstmp(4) = SessionScan(4369);
sstmp(5) = SessionScan(4371);

fh(1) = sstmp(1).plotTrialfyForceh();
fh(2) = sstmp(2).plotTrialfyForceh();
fh(3) = sstmp(3).plotTrialfyForceh();
fh(4) = sstmp(4).plotTrialfyForceh();
fh(5) = sstmp(5).plotTrialfyForceh();
% set(fh(1), 'Title', {'2.5cm target ±1N'});
% title(fh(2), 'No target ±1N')
% title(fh(3), 'No target ±2N')
% title(fh(4), '2.5cm target ±2N')
% title(fh(5), '2.5cm target ±2N, long duration');

% plot out the spectrum of each trial/condition

dat1 = sstmp(1).export_as_formatted(1);
dat2 = sstmp(2).export_as_formatted(1);
dat3 = sstmp(3).export_as_formatted(1);
dat4 = sstmp(4).export_as_formatted(1);
dat5 = sstmp(5).export_as_formatted(1);
% data = cell(2,1,15,5);
data(2,1,:,1) = dat1(1,1,1:15,1);
data(2,1,:,2) = dat2(1,1,1:15,1);
data(2,1,:,3) = dat3(1,1,1:15,1);
data(2,1,:,4) = dat4(1,1,1:15,1);
data(2,1,:,5) = dat5(1,1,1:15,1);


ref2_t_range = [1784 1800];
ref2_t_idx = sstmp(1).data.t>ref2_t_range(1) & sstmp(1).data.t<ref2_t_range(2); 
ref2_t_shift_idx = 212644;

dat_ref2.t = sstmp(1).data.t(ref2_t_idx) - sstmp(1).data.t(ref2_t_shift_idx);
dat_ref2.x = sstmp(1).data.x(:,ref2_t_idx);
dat_ref2.f = sstmp(1).data.f(:,ref2_t_idx);
dat_ref2.v = sstmp(1).data.v(:,ref2_t_idx);
dat_ref2.ftq= sstmp(1).data.ftq(:,ref2_t_idx);
dat_ref2.Fp = sstmp(1).data.Fp(:,ref2_t_idx);
dat_ref2.ts = sstmp(1).data.ts(:,ref2_t_idx);
dat_ref2.tq = sstmp(1).data.tq(:,ref2_t_idx);
dat_ref2.emg = sstmp(1).data.emg(:,ref2_t_idx);
% dat_ref2.ox = sstmp(1).data.ox(:,ref2_t_idx);
% dat_ref2.ov = sstmp(1).data.ov(:,ref2_t_idx);
% dat_ref2.mvst = sstmp(1).data.mvst(:,ref2_t_idx);

save('20220916cgHATesting.mat','data');

save('20220916cgHATesting.mat','dat_ref','dat_ref2','-append');
%% plot a specific trial data 
clear; 
load('20220916cgHATesting.mat','dat*');
% dat_ref = dat1{1,1,17,1}; % the 'shaking' condition
% dat_ref = dat1{1,1,15,1}; % 1N 2.5 target
dat_ref = dat_ref2;
% dat_ref = data{2,1,11,1};
% dat_ref = data{1,2,3,9};
% dat_ref = data{1,1,12,5};
% dat_ref = data{1,2,3,1,9,1}; 
% t_idx = dat_ref.ts == 4;
t_idx = dat_ref.t>-1 & dat_ref.t<0; 
figure('unit', 'inch', 'position', [0 0 4 4]); 
subplot(2,1,1);
grid on; 
plot(dat_ref.t(t_idx), dat_ref.f(:,t_idx) - mean(dat_ref.f(:,t_idx),2), 'linewidth', 1);
xlabel('t (s)'); ylabel('average subtracted Force (N) ' );
legend('x', 'y', 'z');
title('force in a deliberate shaking');
grid on; 
subplot(2,1,2);
[p,f] = pspectrum(dat_ref.f(:,t_idx)', 500);
plot(f, p, 'linewidth', 2); 
xlim([0 20]);
ylim([0 0.02])
xlabel('Freq (Hz)'); ylabel('power' );
legend('PS-x', 'PS-y', 'PS-z');
title('PowerSpectrum in a deliberate shaking');
grid on; 
sgtitle('An example of a deliberate shaking');

% % plot the magnitude and direction of the shaking...
% t = dat_ref.t(t_idx);
% f_tmp = dat_ref.f(:,t_idx);
% f_tmp_high = highpass(f_tmp', 3, 500)';
% % f_rmvoffset = f_tmp - mean(f_tmp);
% f_rmvoffset = f_tmp_high;
% [f_rmvoffset_d, f_rmvoffset_l] = cart2pol(f_rmvoffset(1,:), f_rmvoffset(2,:));
% subplot(2,2,1); 
% plot(t, f_rmvoffset_l); 
% xlabel('t (s)'); ylabel('\Delta F (N)'); 
% title('shaking force net value'); 
% subplot(2,2,3); 
% plot(t, f_rmvoffset_d); 
% xlabel('t (s)'); ylabel('Degree F (rad)'); 
% title('shaking force direction'); 
% subplot(2,2,2)
% plot(f_tmp_high(1,:), f_tmp_high(2,:)); 
% grid on; 
% xlabel('Fx (N)'); ylabel('Fy (N)');
% title('tremer force plot');

% plot the magnitude and direction of the shaking...
figure('unit', 'inch', 'position', [0 0 5 2]);
t = dat_ref.t(t_idx);
f_tmp = dat_ref.f(:,t_idx);
f_tmp_high = highpass(f_tmp', 3, 500)';
% f_rmvoffset = f_tmp - mean(f_tmp,2);
subplot('position', [0.1 0.2 0.5 0.6]); grid on;
plot(t, f_tmp_high, 'linewidth', 1 ); 
legend('x', 'y', 'z');
xlabel('t (s)'); ylabel('\Delta F (N)'); 
title('shaking force net value'); 
% subplot(1,2,2); 
subplot('position', [0.70 0.2 0.25 0.6]); 
plot(f_tmp_high(1,:), f_tmp_high(2,:), '-.'); 
grid on; 
xlabel('Fx (N)'); ylabel('Fy (N)');
title('tremer force plot');
sgtitle('a deliberately shake');
%% 
% check with specific trials... 
for trial_i = 1:9
    dat_ref = data{4,2,3,1,trial_i,1};
    % t_idx = dat_ref.ts == 4;
    t_idx = dat_ref.t>-1 & dat_ref.t<0; 
    
    
    figure('unit', 'inch', 'position', [0 0 5 2]);
    t = dat_ref.t(t_idx);
    f_tmp = dat_ref.f(:,t_idx);
    f_tmp_high = highpass(f_tmp', 1, 500)';
    % f_rmvoffset = f_tmp - mean(f_tmp,2);
    subplot('position', [0.1 0.2 0.5 0.6]); grid on;
    plot(t, f_tmp_high, 'linewidth', 1 ); 
    legend('x', 'y', 'z');
    xlabel('t (s)'); ylabel('\Delta F (N)'); 
    title('shaking force net value'); 
    % subplot(1,2,2); 
    subplot('position', [0.70 0.2 0.25 0.6]); 
    plot(f_tmp_high(1,:), f_tmp_high(2,:), '-.'); 
    grid on; 
    xlabel('Fx (N)'); ylabel('Fy (N)');
    title('tremer force plot');
    sgtitle(['subj1 dir2 fce3 tar1 trial' num2str(trial_i)]);
end

%% deal with matrix data
fh = figure(); 
clear lnh
col_type = colormap('lines'); 
close(fh); 
load('20220916cgHATesting.mat','data');
data = data(2,:,:,:); % Himanshu
sgtitile_arrays = {'2.5cm target. +-1N'...
                    'without seeing target. +-1N' ...
                    'without seeing target. +-2N' ...
                    '2.5cm target. +-2N'...
                    '2.5cm target. +-2N, hold longer'...
    };


for ss_i = 1:5

    figure('unit', 'inch', 'position', [0 0 3.5 3]); 
    axh(1) = subplot(2,1,1); hold on; grid on;
    xlabel('time (s)'); ylabel('net force (N)');
    title('average-subtracted force during holding');

    axh(2) = subplot(2,1,2); hold on; grid on;
%     xlabel('frequency (Hz)');ylabel('power');
    xlabel('log(Hz)');ylabel('dB');
    title('power spectrum of 1s before release');


    %1:5 
    for trial_i = 15:-1:1
        col_i = floor((trial_i-1)/5) + 1;
%         dat_idx = data{1,1,trial_i,ss_i}.ts == 4 |  data{1,1,trial_i,ss_i}.ts == 3; % when hold 
        dat_trial = data{1,1,trial_i,ss_i};
%         t_idx = data{1,1,trial_i,ss_i}.ts == 4;
        t_idx = dat_trial.t>-1 & dat_trial.t<0; 
%         plot(axh(1), data{1,1,trial_i,ss_i}.t(t_idx), ...
%             data{1,1,trial_i,ss_i}.f(1,t_idx) , ...
%             'color', col_type(col_i,:));
        plot(axh(1), data{1,1,trial_i,ss_i}.t(t_idx), ...
            data{1,1,trial_i,ss_i}.f(1,t_idx)-...
            mean(data{1,1,trial_i,ss_i}.f(1,t_idx)), ...
            'color', col_type(col_i,:), 'linewidth', 1);
        
        
        subplot(axh(2));
%         pspectrum(data{1,1,trial_i,ss_i}.f(1,dat_idx), 500);

        f_offset = mean(vecnorm(dat_trial.f(1:3,t_idx)));
%                     plot(dat_trial.t, vecnorm(dat_trial.f(1:3,:)) - f_offset);
        dat_tmp = vecnorm(dat_trial.f(1:3,t_idx)) - f_offset;
        [p,f] = pspectrum(dat_tmp, 500);
        ydb = pow2db(p); 
        logf = log(f);

        if (mod(trial_i,5) == 1)
            lnh(ceil(trial_i/5)) = plot(logf, ydb, 'linewidth', 2, 'color', col_type(col_i,:)); 
        else 
            plot(logf, ydb, 'linewidth', 2, 'color', col_type(col_i,:)); 
        end
%         xlim([0 20]);
%         ylim([0 0.02]);
%         ylim([0 0.05]);
    end
    sgtitle(sgtitile_arrays{ss_i});
end
% legend(lnh, {'1', '2', '3', '4', '5'});
legend(lnh, {'15N', '20N', '25N'});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check some specific trial with EMG raw data. 
% clear; close all; clc;
% load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/' ...
%     'data/processedData/ss4310_4356.mat']);

% subj_i = 1;
% dir_i = 2;
% fce_i = 3;
% dist_i = 1; 
% trial_i = 1;
% pert_i = 1; 
% 
% data_trial = data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i};
% figure('name', 'see the force vs EMG');
% subplot(3,1,1); 
% plot(data_trial.t, data_trial.f(1,:));
% subplot(3,1,2); 
% plot(data_trial.t, data_trial.emg(3,:));
% subplot(3,1,3); 
% plot(data_trial.t, data_trial.emg(4,:));

% cannot see from the tidied-up data. I'll look through the raw data 
% sstmp = SessionScan(4313);
fce_sel = 25;
tar_sel = 0.025; 
trial_sel = find([sstmp.trials.tarF] == fce_sel & [sstmp.trials.tarL] == tar_sel); 
trial_i = 1; 
%  trial_i = 1; 

trial_idx = trial_sel(trial_i);
% figure('name', 'see the force vs EMG');
% axh(1) = subplot(3,1,1); 
% plot(sstmp.trials(trial_idx).data.t - sstmp.trials(trial_idx).data.t(1),...
%     sstmp.trials(trial_idx).data.f(1,:));
% axh(2) = subplot(3,1,2); 
% plot(sstmp.trials(trial_idx).data.t - sstmp.trials(trial_idx).data.t(1),...
%     sstmp.trials(trial_idx).data.emg(3,:));
% axh(3) = subplot(3,1,3); 
% plot(sstmp.trials(trial_idx).data.t - sstmp.trials(trial_idx).data.t(1),...
%     sstmp.trials(trial_idx).data.emg(4,:));
% linkaxes(axh, 'x');

% figure('name', 'see the force vs EMG');
% axh(1) = subplot(3,1,1); 
% plot(sstmp.trials(trial_idx).data.t,...
%     sstmp.trials(trial_idx).data.f(1,:));
% axh(2) = subplot(3,1,2); 
% plot(sstmp.trials(trial_idx).data.t,...
%     sstmp.trials(trial_idx).data.emg(3,:));
% axh(3) = subplot(3,1,3); 
% plot(sstmp.trials(trial_idx).data.t,...
%     sstmp.trials(trial_idx).data.emg(4,:));
% linkaxes(axh, 'x');


% t_sel = [15741.9755664963 15744.0972554099]; % trial 1
% t_sel = [15766.3364559344 15770.4879385487]; % trial 3
ts_val_idx = find(sstmp.trials(trial_idx).data.ts == 4);
t_sel = [sstmp.trials(trial_idx).data.t(min(ts_val_idx)) ...
        sstmp.trials(trial_idx).data.t(max(ts_val_idx))];
t_idx = sstmp.trials(trial_idx).data.t>t_sel(1) & ...
        sstmp.trials(trial_idx).data.t<t_sel(2);
dattmp.t = sstmp.trials(trial_idx).data.t(t_idx);
dattmp.f = sstmp.trials(trial_idx).data.f(1,t_idx);
dattmp.f = dattmp.f - mean(dattmp.f);
dattmp.emg1 = sstmp.trials(trial_idx).data.emg(1,t_idx);
dattmp.emg2 = sstmp.trials(trial_idx).data.emg(2,t_idx);
dattmp.emg3 = sstmp.trials(trial_idx).data.emg(3,t_idx);
dattmp.emg4 = sstmp.trials(trial_idx).data.emg(4,t_idx);

clear('axh'); 
figure('name', 'see the force vs EMG, powerspectrum');
axh(1,1) = subplot(3,2,1); 
plot(dattmp.t, dattmp.f);
axh(1,2) = subplot(3,2,2); 
% power spectrum here
[P,F] = pspectrum(dattmp.f, 500);
plot(F, P);

axh(2,1) = subplot(3,2,3); 
plot(dattmp.t, dattmp.emg3);
axh(2,2) = subplot(3,2,4); 
[P,F] = pspectrum(dattmp.emg3, 500);
plot(F, P);

axh(3,1) = subplot(3,2,5); 
plot(dattmp.t, dattmp.emg4);
xlabel('time'); 
axh(3,2) = subplot(3,2,6); 
[P,F] = pspectrum(dattmp.emg4, 500);
plot(F, P);

linkaxes(axh(:,1), 'x');

%% 
% mscohere plot

x = dattmp.f;
x = x-mean(x);
% y = dattmp.emg2; 
% y = dattmp.emg3; 
y = dattmp.emg4; 
% y = dattmp.emg2; 

[cxy,f] = mscohere(x,y,hamming(200),180,[1:250],500);
% [cxy,f] = mscohere(x,y,hamming(500),450,[1:250],500);

figure(); 
subplot(3,1,1); 
plot(x); 
subplot(3,1,2); 
plot(y); 
subplot(3,1,3); 
plot(f,cxy);

% use un-rectified emg, and do remove DC component. 

% do timelag
figure();
Fs = 500;
[Pxy,F] = cpsd(x,y,hamming(500),450,1:250,500);

Pxy(cxy < 0.4) = 0;

plot(F,angle(Pxy)/pi)
title('Cross Spectrum Phase')
xlabel('Frequency (Hz)')
ylabel('Lag (\times\pi rad)')
grid

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code for find shaking, see 2 muscles..  
% code for checking specific trials 

% clear; close all; clc;
close all;
% cannot see from the tidied-up data. I'll look through the raw data 
% sstmp = SessionScan(4313); % subject 1
% sstmp = SessionScan(4339); % subject 3
% for trial_i = setdiff([1:9], 7)
for trial_i = setdiff([1:9], 5)
% for trial_i = [3 4] % which I cannot find bi-cep and tri-cep activities 
% for trial_i = setdiff([1:9], [])
% for trial_i = 1
fce_sel = 25;
tar_sel = 0.025; 02
trial_sel = find([sstmp.trials.tarF] == fce_sel & [sstmp.trials.tarL] == tar_sel); 
% trial_i = 3; 

trial_idx = trial_sel(trial_i);

% ts_val_idx = find(sstmp.trials(trial_idx).data.ts == 4);
% t_sel = [sstmp.trials(trial_idx).data.t(min(ts_val_idx)) ...
%         sstmp.trials(trial_idx).data.t(max(ts_val_idx))];
ts_val_idx = find(sstmp.trials(trial_idx).states_arr == 4);
t_sel = [sstmp.trials(trial_idx).time_orn(min(ts_val_idx)) ...
    sstmp.trials(trial_idx).time_orn(max(ts_val_idx))];
t_idx = sstmp.trials(trial_idx).data.t>t_sel(1) & ...
        sstmp.trials(trial_idx).data.t<t_sel(2);

% dattmp.t = sstmp.trials(trial_idx).data.t(t_idx);
% dattmp.f = sstmp.trials(trial_idx).data.f(1,t_idx);
% dattmp.f = dattmp.f - mean(dattmp.f);
% dattmp.emg1 = sstmp.trials(trial_idx).data.emg(1,t_idx);
% dattmp.emg2 = sstmp.trials(trial_idx).data.emg(2,t_idx);
% dattmp.emg3 = sstmp.trials(trial_idx).data.emg(3,t_idx);
% dattmp.emg4 = sstmp.trials(trial_idx).data.emg(4,t_idx);
dattmp.ft_idx = sstmp.ft.brtime>t_sel(1) & sstmp.ft.brtime<t_sel(2);
dattmp.ft = sstmp.ft.brtime(dattmp.ft_idx);
dattmp.f  = sstmp.ft.force(1,dattmp.ft_idx);
dattmp.emgt_idx = sstmp.emg.brtime>t_sel(1) & sstmp.emg.brtime<t_sel(2);
dattmp.emgt = sstmp.emg.brtime(dattmp.emgt_idx);
dattmp.emg = sstmp.emg.data.emgraw(:,dattmp.emgt_idx);

% high sample to 2000Hz
datfmt.t = t_sel(1):1/2000:t_sel(2); 
datfmt.f = interp1(dattmp.ft,dattmp.f,datfmt.t,"spline");
datfmt.emg = interp1(dattmp.emgt,dattmp.emg',datfmt.t,"spline")';

% subtract the average of F 
datfmt.f = datfmt.f - mean(datfmt.f);
datfmt.emg = datfmt.emg - mean(datfmt.emg,2);

% maybe it's good to filter out the <2Hz frequency?
datfmt.f_low = highpass(datfmt.f, 5, 2000);
ch_interest = [3 4]; % used to be [3 4]

ifplot = 1;
if (ifplot)
    clear('axh');
    figure('name', 'see the force vs EMG, powerspectrum');
    axh(1,1) = subplot(3,2,1); hold on;
    plot(datfmt.t - datfmt.t(1), datfmt.f); grid on;
%     plot(datfmt.t - datfmt.t(1), datfmt.f_low); 
    xlabel('t (s)');
    ylabel('force (N)');
    title('Force_x');
    axh(1,2) = subplot(3,2,2); hold on; 
    % power spectrum here
    [P,F] = pspectrum(datfmt.f, 2000);  grid on;
    plot(F, P);
%     [P,F] = pspectrum(datfmt.f_low, 2000);  grid on;
%     plot(F, P);
    xlabel('Frequency (Hz)');
    ylabel('power');
    title('power spectrum force');

    axh(2,1) = subplot(3,2,3); 
    plot(datfmt.t - datfmt.t(1), datfmt.emg(ch_interest(1),:)); grid on;
    xlabel('t (s)');
    ylabel('EMG');
    title('EMG bicep');

    axh(2,2) = subplot(3,2,4); 
%     [P,F] = pspectrum(datfmt.emg(3,:), 2000); grid on;
    [P,F] = pspectrum(datfmt.emg(ch_interest(1),:), 2000); grid on;
    plot(F, P);
    xlabel('Frequency (Hz)');
    ylabel('power');
    title('power spectrum');

    axh(3,1) = subplot(3,2,5); 
    plot(datfmt.t - datfmt.t(1), datfmt.emg(ch_interest(2),:)); grid on;
    xlabel('t (s)');
    ylabel('EMG');
    title('EMG tricep');

    axh(3,2) = subplot(3,2,6); 
%     [P,F] = pspectrum(datfmt.emg(4,:), 2000); grid on;
    [P,F] = pspectrum(datfmt.emg(ch_interest(2),:), 2000); grid on;
    plot(F, P);
    xlabel('Frequency (Hz)');
    ylabel('power');
    title('power spectrum');

    linkaxes(axh(:,1), 'x');
    linkaxes(axh(:,2), 'x');
    xlim([5 50]);

    sgtitle(['trial' num2str(trial_i)]);
end

% mscohere plot

x = datfmt.f;
% x = datfmt.f_low;
x = x-mean(x);
y = datfmt.emg; 


% [cxy,f] = mscohere(x',y',hamming(200),180,[1:250],500);
% [cxy,f] = mscohere(x',y',hamming(800),750,[1:250],500);
% [cxy,f] = mscohere(x',y',hamming(1000),950,[1:250],2000);
interest_band = 5:0.2:20; %% Here I changed the interest band...
[cxy,f] = mscohere(x',y',hamming(1500),1450,interest_band,2000);


% use un-rectified emg, and do remove DC component. 

% do timelag
% figure();
Fs = 2000;


PXY = zeros(size(cxy));
for i = 1:8
    [Pxy,F] = cpsd(x,y(i,:),hamming(1500),450,interest_band,2000);
%     Pxy(cxy(i,:) < 0.2) = 0;
    Pxy(cxy(i,:) < 0.1) = 0;
    PXY(i,:) = Pxy;
end


% plot(F,angle(PXY(3:4,:))/pi)
% title('Cross Spectrum Phase')
% xlabel('Frequency (Hz)')
% ylabel('Lag (\times\pi rad)')
% grid

clear('axh'); 
figure(); 
axh(1) = subplot(4,1,1); 
plot(datfmt.t - datfmt.t(1), x); 
title('force'); 
xlabel('time'); ylabel('F (N)');
axh(2) = subplot(4,1,2); 
plot(datfmt.t - datfmt.t(1), y(ch_interest,:)'); 
title('EMG'); 
xlabel('time'); ylabel('F (N)');
axh(3) = subplot(4,1,3); 
plot(f,cxy(ch_interest,:));
title('Cross Coherence'); 
xlabel('Frequency'); ylabel('?');
axh(4) = subplot(4,1,4); hold on; 
plot(F,angle(PXY(ch_interest,:))/pi);
% highlight the plot 
F_idx = F>=8 & F<=12; 
plot(F(F_idx), angle(PXY(ch_interest(2),F_idx))/pi, '.');
title('Phase of coherece');
xlabel('Frequency'); ylabel('Lag (\times\pi rad)')
legend('bicep', 'tricep');

linkaxes(axh(1:2), 'x');
linkaxes(axh(3:4), 'x');

sgtitle(['trial' num2str(trial_i)]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code for find the power spectrum, see 8 muscles... 
% find the peak of Force around 8-15Hz, and then plot the EMG power
% spectrum peak
% also, plot the cross-coherence according to corresponding peaks...

% clear; close all; clc;
close all;
% cannot see from the tidied-up data. I'll look through the raw data 
% sstmp = SessionScan(4310); trials_all = setdiff([1:9], [2 5 6 7 8 9]);
sstmp = SessionScan(4313);  trials_all = setdiff([1:9], [7]);% subject 1
% sstmp = SessionScan(4311);  trials_all = setdiff([1:9], [1 3 4 5 6 9]);% subject 1
% sstmp = SessionScan(4314);  trials_all = setdiff([1:9], [5 6 7 9]);% subject 1
% sstmp = SessionScan(4339); trials_all = setdiff([1:9], 5); % subject 3, dir 2
% sstmp = SessionScan(4336); trials_all = setdiff([1:9], [9]); % subject 3, dir 2
% sstmp = SessionScan(4341); trials_all = setdiff([1:9], 7);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4337); trials_all = setdiff([1:9], [1 4]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4356); trials_all = setdiff([1:9], [3 7]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4355); trials_all = setdiff([1:9], [8]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4351); trials_all = setdiff([1:9], [1:8]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4353); trials_all = setdiff([1:9], [1 9]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4329); trials_all = setdiff([1:9], [7 8]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4326); trials_all = setdiff([1:9], [1 2 3 4 5 7]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4325); trials_all = setdiff([1:9], [1 2 3 4 8 9 6]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4328); trials_all = setdiff([1:9], [9]);% subject 3, dir 4 trial 8 out 
fh = figure(); 
% for trial_i = setdiff([1:9], 7)
% for trial_i = setdiff([1:9], 5)
% for trial_i = setdiff([1:9], 6, 8)
% for trial_i = [3 4] % which I cannot find bi-cep and tri-cep activities 
% for trial_i = setdiff([1:9], [])
% for trial_i = 1
for trial_i = trials_all
fce_sel = 25;
tar_sel = 0.025; 
trial_sel = find([sstmp.trials.tarF] == fce_sel & [sstmp.trials.tarL] == tar_sel); 
% trial_i = 3; 

trial_idx = trial_sel(trial_i);

% ts_val_idx = find(sstmp.trials(trial_idx).data.ts == 4);
% t_sel = [sstmp.trials(trial_idx).data.t(min(ts_val_idx)) ...
%         sstmp.trials(trial_idx).data.t(max(ts_val_idx))];
ts_val_idx = find(sstmp.trials(trial_idx).states_arr == 4);
t_sel = [sstmp.trials(trial_idx).time_orn(min(ts_val_idx)) ...
    sstmp.trials(trial_idx).time_orn(max(ts_val_idx))];
t_idx = sstmp.trials(trial_idx).data.t>t_sel(1) & ...
        sstmp.trials(trial_idx).data.t<t_sel(2);

% dattmp.t = sstmp.trials(trial_idx).data.t(t_idx);
% dattmp.f = sstmp.trials(trial_idx).data.f(1,t_idx);
% dattmp.f = dattmp.f - mean(dattmp.f);
% dattmp.emg1 = sstmp.trials(trial_idx).data.emg(1,t_idx);
% dattmp.emg2 = sstmp.trials(trial_idx).data.emg(2,t_idx);
% dattmp.emg3 = sstmp.trials(trial_idx).data.emg(3,t_idx);
% dattmp.emg4 = sstmp.trials(trial_idx).data.emg(4,t_idx);
dattmp.ft_idx = sstmp.ft.brtime>t_sel(1) & sstmp.ft.brtime<t_sel(2);
dattmp.ft = sstmp.ft.brtime(dattmp.ft_idx);
dattmp.f  = sstmp.ft.force(1,dattmp.ft_idx);
dattmp.emgt_idx = sstmp.emg.brtime>t_sel(1) & sstmp.emg.brtime<t_sel(2);
dattmp.emgt = sstmp.emg.brtime(dattmp.emgt_idx);
dattmp.emg = sstmp.emg.data.emgraw(:,dattmp.emgt_idx);



% high sample to 2000Hz
datfmt.t = t_sel(1):1/2000:t_sel(2); 
datfmt.f = interp1(dattmp.ft,dattmp.f,datfmt.t,"spline");
datfmt.emg = interp1(dattmp.emgt,dattmp.emg',datfmt.t,"spline")';
% datfmt.x = interp(dattmp.wam.t,datatmp.wam.x, datfmt.t, "spline");

% subtract the average of F 
datfmt.f = datfmt.f - mean(datfmt.f);
datfmt.emg = datfmt.emg - mean(datfmt.emg,2);

% maybe it's good to filter out the <2Hz frequency?
datfmt.f_low = highpass(datfmt.f, 5, 2000);
% ch_interest = [3 4]; % used to be [3 4]
ch_interest = [1:8]; % used to be [3 4]

[P,F] = pspectrum(datfmt.f, 2000);  grid on;
F_interest_idx = find(F > 8 & F<16); 
[p_max, m_idx] = max(P(F_interest_idx));
f_max = F(F_interest_idx(m_idx));

muscles_arr = {'FCR','ECU','BI','TRI','AD','PD','PEC','TPZ'};
ifplot = 1;
if (ifplot)
    clear('axh');
    figure('name', 'see the force vs EMG, powerspectrum', ...
        'position', [0 0 600 1000]);
    axh(1,1) = subplot(9,2,1); hold on;
    plot(datfmt.t - datfmt.t(1), datfmt.f); grid on;
    xlabel('t (s)');
    ylabel('force (N)');
    title('Force_x');
    axh(1,2) = subplot(9,2,2); hold on; 
    % power spectrum here
    [P,F] = pspectrum(datfmt.f, 2000);  grid on;
    plot(F, P);
    plot(f_max, p_max, 's');  % the marker of maximum
    xlabel('Frequency (Hz)');
    ylabel('power');
    title('power spectrum force');

%     subplot(9,2,3); hold on;
%     plot()
    for ii = 1:8
        axh(ii+1,1) = subplot(9,2,2+(ii-1)*2+1); 
        plot(datfmt.t - datfmt.t(1), datfmt.emg(ch_interest(ii),:)); grid on;
        xlabel('t (s)');
        ylabel('EMG');
        title(['EMG ' muscles_arr{ii}]);
        axh(ii+1,2) = subplot(9,2,2+(ii-1)*2+2);  
        [P,F] = pspectrum(datfmt.emg(ch_interest(ii),:), 2000); grid on;
        plot(F, P);
        xline(f_max); 
        xlabel('Frequency (Hz)');
        ylabel('power');
        title('power spectrum');

    end

    linkaxes(axh(:,1), 'x');
    linkaxes(axh(:,2), 'x');
    xlim([5 50]);

    sgtitle(['trial' num2str(trial_i)]);
end


% mscohere plot
x = datfmt.f;
% x = datfmt.f_low;
x = x-mean(x);
y = datfmt.emg; 


% [cxy,f] = mscohere(x',y',hamming(200),180,[1:250],500);
% [cxy,f] = mscohere(x',y',hamming(800),750,[1:250],500);
% [cxy,f] = mscohere(x',y',hamming(1000),950,[1:250],2000);
interest_band = 5:0.2:20; %% Here I changed the interest band...
% interest_band = 5:1:20; %% Here I changed the interest band...
% interest_band = 1:1:250; %% Here I changed the interest band...
[cxy,f] = mscohere(x',y',hamming(1500),1450,interest_band,2000);


% use un-rectified emg, and do remove DC component. 

% do timelag
% figure();
Fs = 2000;


PXY = zeros(size(cxy));
for i = 1:8
    [Pxy,F] = cpsd(x,y(i,:),hamming(1500),1450,interest_band,2000);
%     Pxy(cxy(i,:) < 0.2) = 0;
    Pxy(cxy(i,:) < 0.5) = 0;
    PXY(i,:) = Pxy;
end


% plot(F,angle(PXY(3:4,:))/pi)
% title('Cross Spectrum Phase')
% xlabel('Frequency (Hz)')
% ylabel('Lag (\times\pi rad)')
% grid
% 
clear('axh'); 
figure('position', [600 0 600 1000]); 
% figure(fh); 
for ii = 1:8
axh(1) = subplot(8,2,2*(ii-1)+1); hold on;
cxy_tmp = cxy(ii,:); 
cxy_tmp(cxy_tmp<0.5) = 0;
% plot(f, cxy(ii,:)); 
% plot(f, cxy_tmp); 
plot(f, cxy_tmp, '.'); 
ylim([0 1]);
xline(f_max); grid on;
title(['magnitude ' muscles_arr{ii}]); 
xlabel('Freq (Hz)'); ylabel('power');

axh(2) = subplot(8,2,2*(ii-1)+2); hold on;
% plot(datfmt.t - datfmt.t(1), y(ch_interest,:)'); 
plot(F,angle(PXY(ii,:))/pi);
xline(f_max); grid on; 
title(['Phase ' muscles_arr{ii}]); 
xlabel('Freq (Hz)'); ylabel('rad/\pi');

sgtitle(['trial' num2str(trial_i)]);
end
% Phase<0, force after EMG; Phase>0, force ahead of EMG; 

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code for find the power spectrum, see 8 muscles... 
% find the peak of Force around 8-15Hz, and then plot the EMG power
% spectrum peak
% also, plot the cross-coherence according to corresponding peaks...
% construct a matrix: muscles-by-trials. 
% If a trial is considered to be having shake, check all muscles' EMG
% cross-coherence with in this trial. Take the cross coherence bigger than
% 0.5 as a valid one. 
% For each muscle, if one condition the trial number is greater than half
% of the trial, it will be regarded "responsable to" the force. And the
% number in corresponding matrix will be set to 1. 


% clear; close all; clc;
close all;

x_coh_greater = nan(8,9); 
% sstmp = SessionScan(4313); % subject 1, direction 2
% sstmp = SessionScan(4339); % subject 3, direction 2
fh = figure(); 

% 
% for trial_i = setdiff([1:9], 7)
% for trial_i = setdiff([1:9], 5)
% for trial_i = setdiff([1:9], [8 9])
for trial_i = trials_all
        % for trial_i = [3 4] % which I cannot find bi-cep and tri-cep activities 
        % for trial_i = setdiff([1:9], [])
        % for trial_i = 1
fce_sel = 25;
tar_sel = 0.025; 
trial_sel = find([sstmp.trials.tarF] == fce_sel & [sstmp.trials.tarL] == tar_sel); 

trial_idx = trial_sel(trial_i);

% ts_val_idx = find(sstmp.trials(trial_idx).data.ts == 4);
% t_sel = [sstmp.trials(trial_idx).data.t(min(ts_val_idx)) ...
%         sstmp.trials(trial_idx).data.t(max(ts_val_idx))];
ts_val_idx = find(sstmp.trials(trial_idx).states_arr == 4);
t_sel = [sstmp.trials(trial_idx).time_orn(min(ts_val_idx)) ...
    sstmp.trials(trial_idx).time_orn(max(ts_val_idx))];
t_idx = sstmp.trials(trial_idx).data.t>t_sel(1) & ...
        sstmp.trials(trial_idx).data.t<t_sel(2);

dattmp.ft_idx = sstmp.ft.brtime>t_sel(1) & sstmp.ft.brtime<t_sel(2);
dattmp.ft = sstmp.ft.brtime(dattmp.ft_idx);
dattmp.f  = sstmp.ft.force(1,dattmp.ft_idx);
dattmp.emgt_idx = sstmp.emg.brtime>t_sel(1) & sstmp.emg.brtime<t_sel(2);
dattmp.emgt = sstmp.emg.brtime(dattmp.emgt_idx);
dattmp.emg = sstmp.emg.data.emgraw(:,dattmp.emgt_idx);

% high sample to 2000Hz
datfmt.t = t_sel(1):1/2000:t_sel(2); 
datfmt.f = interp1(dattmp.ft,dattmp.f,datfmt.t,"spline");
datfmt.emg = interp1(dattmp.emgt,dattmp.emg',datfmt.t,"spline")';

% subtract the average of F 
datfmt.f = datfmt.f - mean(datfmt.f);
datfmt.emg = datfmt.emg - mean(datfmt.emg,2);

% maybe it's good to filter out the <2Hz frequency?
datfmt.f_low = highpass(datfmt.f, 5, 2000);
% ch_interest = [3 4]; % used to be [3 4]
ch_interest = [1:8]; % used to be [3 4]

[P,F] = pspectrum(datfmt.f, 2000);  grid on;
F_interest_idx = find(F > 8 & F<16); 
[p_max, m_idx] = max(P(F_interest_idx));
f_max = F(F_interest_idx(m_idx));

muscles_arr = {'FCR','ECU','BI','TRI','AD','PD','PEC','TPZ'};
ifplot = 0;
if (ifplot)
    clear('axh');
    figure('name', 'see the force vs EMG, powerspectrum', ...
        'position', [0 0 600 1000]);
    axh(1,1) = subplot(9,2,1); hold on;
    plot(datfmt.t - datfmt.t(1), datfmt.f); grid on;
    xlabel('t (s)');
    ylabel('force (N)');
    title('Force_x');
    axh(1,2) = subplot(9,2,2); hold on; 
    % power spectrum here
    [P,F] = pspectrum(datfmt.f, 2000);  grid on;
    plot(F, P);
    plot(f_max, p_max, 's');  % the marker of maximum
    xlabel('Frequency (Hz)');
    ylabel('power');
    title('power spectrum force');

    for ii = 1:8
        axh(ii+1,1) = subplot(9,2,2+(ii-1)*2+1); 
        plot(datfmt.t - datfmt.t(1), datfmt.emg(ch_interest(ii),:)); grid on;
        xlabel('t (s)');
        ylabel('EMG');
        title(['EMG ' muscles_arr{ii}]);
        axh(ii+1,2) = subplot(9,2,2+(ii-1)*2+2);  
        [P,F] = pspectrum(datfmt.emg(ch_interest(ii),:), 2000); grid on;
        plot(F, P);
        xline(f_max); 
        xlabel('Frequency (Hz)');
        ylabel('power');
        title('power spectrum');

    end

    linkaxes(axh(:,1), 'x');
    linkaxes(axh(:,2), 'x');
    xlim([5 50]);

    sgtitle(['trial' num2str(trial_i)]);
end


% mscohere plot
x = datfmt.f;
% x = datfmt.f_low;
x = x-mean(x);
y = datfmt.emg; 

interest_band = 5:0.2:20; %% Here I changed the interest band...
[cxy,f] = mscohere(x',y',hamming(1500),1450,interest_band,2000);
[~, f_max_idx] = min(abs(f - f_max));
xcoh_threshold = 0.5; 

for i = 1:8 % do the judgement of the cross coherence 
    x_coh_greater(i, trial_i) = cxy(i,f_max_idx)>xcoh_threshold;
end

% use un-rectified emg, and do remove DC component. 

% do timelag
% figure();
Fs = 2000;


PXY = zeros(size(cxy));
for i = 1:8
    [Pxy,F] = cpsd(x,y(i,:),hamming(1500),1450,interest_band,2000);
%     Pxy(cxy(i,:) < 0.2) = 0;
    Pxy(cxy(i,:) < xcoh_threshold) = 0;
    PXY(i,:) = Pxy;
end


% plot(F,angle(PXY(3:4,:))/pi)
% title('Cross Spectrum Phase')
% xlabel('Frequency (Hz)')
% ylabel('Lag (\times\pi rad)')
% grid
% 
clear('axh'); 
% figure(fh, 'position', [600 0 600 1000]); 
figure(fh); 
for ii = 1:8
axh(1) = subplot(8,2,2*(ii-1)+1); hold on;
cxy_tmp = cxy(ii,:); 
cxy_tmp(cxy_tmp<0.5) = 0;
% plot(f, cxy(ii,:)); 
% plot(f, cxy_tmp); 
plot(f, cxy_tmp, '.', 'MarkerFaceColor', col_type(trial_i,:)); 
ylim([0 1]);
xline(f_max, 'Color', col_type(trial_i,:)); grid on;
title(['magnitude ' muscles_arr{ii}]); 
xlabel('Freq (Hz)'); ylabel('power');

axh(2) = subplot(8,2,2*(ii-1)+2); hold on;
% plot(datfmt.t - datfmt.t(1), y(ch_interest,:)'); 
plot(F,angle(PXY(ii,:))/pi);
xline(f_max); grid on; 
title(['Phase ' muscles_arr{ii}]); 
xlabel('Freq (Hz)'); ylabel('rad/\pi');

sgtitle(['trial' num2str(trial_i)]);
end
% Phase<0, force after EMG; Phase>0, force ahead of EMG; 
end

valid_trials = sum(~isnan(sum(x_coh_greater)));
sum(x_coh_greater, 2, 'omitnan')/valid_trials % if greater than 0.5
                                                 % over half of the trials
                                                 % have the
                                                 % "crosscoherence"


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Himanshu's idea, first seeing if there is coherence among different 
% muslces. If there is coherence, then check the force if there is any
% "shaking" performance 


% clear; close all; clc;
close all;
% cannot see from the tidied-up data. I'll look through the raw data 

% sstmp = SessionScan(4310); trials_all = setdiff([1:9], [7 8 ]);
% sstmp = SessionScan(4311);  trials_all = setdiff([1:9], [1 3 4 5 6 9]);% subject 1
% sstmp = SessionScan(4313);  trials_all = setdiff([1:9], [7]);% subject 1
% sstmp = SessionScan(4314);  trials_all = setdiff([1:9], [5 6 7 9]);% subject 1

% sstmp = SessionScan(4325); trials_all = setdiff([1:9], [1 2 3 4 8 9 6]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4326); trials_all = setdiff([1:9], [1 2 3 4 5 7]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4328); trials_all = setdiff([1:9], [9]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4329); trials_all = setdiff([1:9], [7 8]);% subject 3, dir 4 trial 8 out 

sstmp = SessionScan(4336); trials_all = setdiff([1:9], [9]); % subject 3, dir 2
% sstmp = SessionScan(4337); trials_all = setdiff([1:9], [1 4]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4339); trials_all = setdiff([1:9], 5); % subject 3, dir 2
% sstmp = SessionScan(4341); trials_all = setdiff([1:9], 7);% subject 3, dir 4 trial 8 out 

% sstmp = SessionScan(4351); trials_all = setdiff([1:9], [1:8]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4353); trials_all = setdiff([1:9], [1 9]);% subject 3, dir 4 trial 8 ou
% sstmp = SessionScan(4355); trials_all = setdiff([1:9], [8]);% subject 3, dir 4 trial 8 out 
% sstmp = SessionScan(4356); trials_all = setdiff([1:9], [3 7]);% subject 3, dir 4 trial 8 out 


fh = figure(); 
% for trial_i = setdiff([1:9], 7)
% for trial_i = setdiff([1:9], 5)
% for trial_i = setdiff([1:9], 6, 8)
% for trial_i = [3 4] % which I cannot find bi-cep and tri-cep activities 
% for trial_i = setdiff([1:9], [])
% for trial_i = 1
for trial_i = trials_all
fce_sel = 25;
tar_sel = 0.025; 
trial_sel = find([sstmp.trials.tarF] == fce_sel & [sstmp.trials.tarL] == tar_sel); 
% trial_i = 3; 

trial_idx = trial_sel(trial_i);

% ts_val_idx = find(sstmp.trials(trial_idx).data.ts == 4);
% t_sel = [sstmp.trials(trial_idx).data.t(min(ts_val_idx)) ...
%         sstmp.trials(trial_idx).data.t(max(ts_val_idx))];
ts_val_idx = find(sstmp.trials(trial_idx).states_arr == 4);
t_sel = [sstmp.trials(trial_idx).time_orn(min(ts_val_idx)) ...
    sstmp.trials(trial_idx).time_orn(max(ts_val_idx))];
t_idx = sstmp.trials(trial_idx).data.t>t_sel(1) & ...
        sstmp.trials(trial_idx).data.t<t_sel(2);

% dattmp.t = sstmp.trials(trial_idx).data.t(t_idx);
% dattmp.f = sstmp.trials(trial_idx).data.f(1,t_idx);
% dattmp.f = dattmp.f - mean(dattmp.f);
% dattmp.emg1 = sstmp.trials(trial_idx).data.emg(1,t_idx);
% dattmp.emg2 = sstmp.trials(trial_idx).data.emg(2,t_idx);
% dattmp.emg3 = sstmp.trials(trial_idx).data.emg(3,t_idx);
% dattmp.emg4 = sstmp.trials(trial_idx).data.emg(4,t_idx);
dattmp.ft_idx = sstmp.ft.brtime>t_sel(1) & sstmp.ft.brtime<t_sel(2);
dattmp.ft = sstmp.ft.brtime(dattmp.ft_idx);
dattmp.f  = sstmp.ft.force(1,dattmp.ft_idx);
dattmp.emgt_idx = sstmp.emg.brtime>t_sel(1) & sstmp.emg.brtime<t_sel(2);
dattmp.emgt = sstmp.emg.brtime(dattmp.emgt_idx);
dattmp.emg = sstmp.emg.data.emgraw(:,dattmp.emgt_idx);

% high sample to 2000Hz
datfmt.t = t_sel(1):1/2000:t_sel(2); 
datfmt.f = interp1(dattmp.ft,dattmp.f,datfmt.t,"spline");
datfmt.emg = interp1(dattmp.emgt,dattmp.emg',datfmt.t,"spline")';

% subtract the average of F 
datfmt.f = datfmt.f - mean(datfmt.f);
datfmt.emg = datfmt.emg - mean(datfmt.emg,2);

% maybe it's good to filter out the <2Hz frequency?
datfmt.f_low = highpass(datfmt.f, 5, 2000);
% ch_interest = [3 4]; % used to be [3 4]
ch_interest = [1:8]; % used to be [3 4]

[P,F] = pspectrum(datfmt.f, 2000);  grid on;
F_interest_idx = find(F > 4 & F<16); 
[p_max, m_idx] = max(P(F_interest_idx));
f_max = F(F_interest_idx(m_idx));

muscles_arr = {'FCR','ECU','BI','TRI','AD','PD','PEC','TPZ'};
ifplot = 1;
if (ifplot)
    clear('axh');
    figure('name', 'see the force vs EMG, powerspectrum', ...
        'position', [0 0 600 1000]);
    axh(1,1) = subplot(9,2,1); hold on;
    plot(datfmt.t - datfmt.t(1), datfmt.f); grid on;
    xlabel('t (s)');
    ylabel('force (N)');
    title('Force_x');
    axh(1,2) = subplot(9,2,2); hold on; 
    % power spectrum here
    [P,F] = pspectrum(datfmt.f, 2000);  grid on;
    plot(F, P);
    plot(f_max, p_max, 's');  % the marker of maximum
    xlabel('Frequency (Hz)');
    ylabel('power');
    title('power spectrum force');

    for ii = 1:8
        axh(ii+1,1) = subplot(9,2,2+(ii-1)*2+1); 
        plot(datfmt.t - datfmt.t(1), datfmt.emg(ch_interest(ii),:)); grid on;
        xlabel('t (s)');
        ylabel('EMG');
        title(['EMG ' muscles_arr{ii}]);
        axh(ii+1,2) = subplot(9,2,2+(ii-1)*2+2);  
        [P,F] = pspectrum(datfmt.emg(ch_interest(ii),:), 2000); grid on;
        plot(F, P);
        xline(f_max); 
        xlabel('Frequency (Hz)');
        ylabel('power');
        title('power spectrum');

    end

    linkaxes(axh(:,1), 'x');
    linkaxes(axh(:,2), 'x');
    xlim([5 50]);

    sgtitle(['trial' num2str(trial_i)]);
end


% mscohere plot
x = datfmt.f;
% x = datfmt.f_low;
x = x-mean(x);
y = datfmt.emg; 


% [cxy,f] = mscohere(x',y',hamming(200),180,[1:250],500);
% [cxy,f] = mscohere(x',y',hamming(800),750,[1:250],500);
% [cxy,f] = mscohere(x',y',hamming(1000),950,[1:250],2000);
interest_band = 1:0.2:50; %% Here I changed the interest band...
% interest_band = 5:1:20; %% Here I changed the interest band...
% interest_band = 1:1:250; %% Here I changed the interest band...
[cxy,f] = mscohere(x',y',hamming(1500),1450,interest_band,2000);
if (length(x) < 1500*2)
    disp(['trial' num2str(trial_i) 'err! NOT enough datapoints!']);
end

% use un-rectified emg, and do remove DC component. 

% do timelag
% figure();
Fs = 2000;


PXY = zeros(size(cxy));
for i = 1:8
    [Pxy,F] = cpsd(x,y(i,:),hamming(1500),1450,interest_band,2000);
%     Pxy(cxy(i,:) < 0.2) = 0;
    Pxy(cxy(i,:) < 0.5) = 0;
    PXY(i,:) = Pxy;
end


% plot(F,angle(PXY(3:4,:))/pi)
% title('Cross Spectrum Phase')
% xlabel('Frequency (Hz)')
% ylabel('Lag (\times\pi rad)')
% grid
% 
clear('axh'); 
figure('position', [600 0 600 1000]); 
% figure(fh); 
for ii = 1:8
axh(1) = subplot(8,2,2*(ii-1)+1); hold on;
cxy_tmp = cxy(ii,:); 
cxy_tmp(cxy_tmp<0.5) = 0;
% plot(f, cxy(ii,:)); 
% plot(f, cxy_tmp); 
plot(f, cxy_tmp, '.'); 
ylim([0 1]);
xline(f_max); grid on;
title(['magnitude ' muscles_arr{ii}]); 
xlabel('Freq (Hz)'); ylabel('power');

axh(2) = subplot(8,2,2*(ii-1)+2); hold on;
% plot(datfmt.t - datfmt.t(1), y(ch_interest,:)'); 
plot(F,angle(PXY(ii,:))/pi);
xline(f_max); grid on; 
title(['Phase ' muscles_arr{ii}]); 
xlabel('Freq (Hz)'); ylabel('rad/\pi');

sgtitle(['trial' num2str(trial_i)]);
end
% Phase<0, force after EMG; Phase>0, force ahead of EMG; 

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% see if the force power spectrum is different 
clear; close all; clc;
load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/' ...
    'data/processedData/ss4310_4356.mat']);

% t_range = [-0.5 0];
t_range = [-0.5 0];
fce_list = [15 20 25];
dist_list = [2.5 5.0 7.5];

%%
%% plot force and power spectrum on one condition 
for subj_i = 3
    for dir_i = 2        
        figure(); 
        for fce_i = 1:3
            for dist_i = 1:3
%                 figure(); hold on;
                axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i); hold on;
                
                P_set = [];
                 for trial_i = 1:9
                     datatmp = data{subj_i,dir_i,fce_i,dist_i,trial_i,1};
                     fce_idx = datatmp.t> t_range(1) & datatmp.t<t_range(2);                      
%                     subplot(1,2,1); hold on;
                        plot(datatmp.t(fce_idx), datatmp.f(1,fce_idx));
%                     subplot(1,2,2); hold on;
%                        [P,F] = pspectrum(datatmp.f(fce_idx),500, 'FrequencyLimits',[100 290],'TimeResolution',1);
                        [P,F] = pspectrum(datatmp.f(fce_idx),500);
%                         P_set = [P_set P];
%                         F
%                         plot(F, '.')
%                        plot(log10(F), pow2db(P));
%                         semilogx(F, pow2db(P)); hold on;
                 end
%                  semilogx(F, pow2db(P_set)); grid on;
                 title(['fce' num2str(fce_list(fce_i)) 'dist' num2str(dist_list(dist_i))]);
%                  xlabel('Hz');                 ylabel('power');
                xlabel('log10(Hz)'); ylabel('dB'); grid on;
            end
        end
        linkaxes(axh, 'xy');
%         xlim([0 1.5]); 
        sgtitle(['Force power spectrum ' 'subject' num2str(subj_i) 'direction' num2str(dir_i)]);
    end
end
%% plot force and power spectrum on one condition 
for subj_i = 3
    for dir_i = 2        
        figure(); 
        for fce_i = 1:3
            for dist_i = 1:3
%                 figure(); hold on;
                axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3+dist_i); hold on;
                
                P_set = [];
                 for trial_i = 1:9
                     datatmp = data{subj_i,dir_i,fce_i,dist_i,trial_i,1};
                     fce_idx = datatmp.t> t_range(1) & datatmp.t<t_range(2);                      
%                     subplot(1,2,1); hold on;
%                         plot(datatmp.t(fce_idx), datatmp.f(1,fce_idx));
%                     subplot(1,2,2); hold on;
%                        [P,F] = pspectrum(datatmp.f(fce_idx),500, 'FrequencyLimits',[100 290],'TimeResolution',1);
                        [P,F] = pspectrum(datatmp.f(fce_idx),500);
%                         P_set = [P_set P];
%                         F
%                         plot(F, '.')
                       plot(log10(F), pow2db(P));
%                         semilogx(F, pow2db(P)); hold on;
                 end
%                  semilogx(F, pow2db(P_set)); grid on;
                 title(['fce' num2str(fce_list(fce_i)) 'dist' num2str(dist_list(dist_i))]);
%                  xlabel('Hz');                 ylabel('power');
                xlabel('log10(Hz)'); ylabel('dB'); grid on;
            end
        end
        linkaxes(axh, 'xy');
        xlim([0 1.5]); 
        sgtitle(['Force power spectrum ' 'subject' num2str(subj_i) 'direction' num2str(dir_i)]);
    end
end