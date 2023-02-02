% Check maximum voluntary 

% Check data ss4284 for maximum voluntary force 

sstmp = SessionScanFT(4284);
sstmp.plotForceOrigin(); 

% time for maximum voluntary force 
idx_range{1} = [2.0e4 2.5e4;
                3.3e4 3.8e4];
idx_range{2} = [4.6e4 5.1e4
                6.0e4 6.5e4];
idx_range{3} = [1.24e5 1.29e5
                1.4e5  1.45e5];
idx_range{4} = [1.59e5 1.64e5
                1.79e5 1.84e5];

idx_range{5} = [0.8e5 1.0e5];
% calculate the average force of each direction 
idx_all = 1:size(sstmp.force_origin,2); 
fce_avg = cell(4,1);
for dir_i = 1:4
    for rep_i = 1:2
        idx_dat = idx_all > idx_range{dir_i}(rep_i,1) & ...
                    idx_all < idx_range{dir_i}(rep_i,2);
        idx_mean = idx_all > idx_range{5}(1,1) & ...
                    idx_all < idx_range{5}(1,2);
        fce_avg{dir_i}(rep_i) = vecnorm(...
                                mean(sstmp.force_origin(:,idx_dat),2) - ...
                                mean(sstmp.force_origin(:,idx_mean),2));
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check with Himanshu's data: ss4315 & ss4312

sstmp = SessionScan(4315); 
figure(); 
axh = zeros(5,1);
axh(1) = subplot(5,1,1);
plot(sstmp.data.t, sstmp.data.f); 
for pair_i = 1:4
    axh(pair_i+1) = subplot(5,1,pair_i+1); hold on;
    plot(sstmp.data.t, sstmp.data.emg((pair_i-1)*2+1,:));
    plot(sstmp.data.t,-sstmp.data.emg((pair_i-1)*2+2,:));
end
    linkaxes(axh, 'x');

    %%


sstmp = SessionScan(4312); 
figure(); 
axh = zeros(5,1);
axh(1) = subplot(5,1,1);
plot(sstmp.data.t, sstmp.data.f); 
for pair_i = 1:4
    axh(pair_i+1) = subplot(5,1,pair_i+1); hold on;
    plot(sstmp.data.t, sstmp.data.emg((pair_i-1)*2+1,:));
    plot(sstmp.data.t,-sstmp.data.emg((pair_i-1)*2+2,:));
end
    linkaxes(axh, 'x');

%% from MVF get the devision value of EMG  
ss_num = {[4315, 4312]...
          [4324, 4327]...
          [4335, 4338]...
          [4349, 4354]};

t_range = cell(4,2); 
t_range{1,1} = [18212.3 18221.8; ...
                18231.0 18241.6;...
                18264.5 18273.3; ...
                18280.8 18289.2];
t_range{1,2} = [14835.2 14844.1;...
                14865.9 14874.1;...
                14891.3 14898.1;...
                14910.2 14919.5];

t_range{2,1} = [9643.5198119585 9653.67652953104;
                9664.26119442529 9674.18394907395;
                9687.31029045808 9695.88521726528;
                9703.44427294688 9712.17699027539];
t_range{2,2} = [13101.2422090639 13110.6810426836;
                13120.5498060513 13130.4625489696;
                13139.6373848838 13149.1541786017;
                13156.4312549322 13165.2501463469];

t_range{3,1} = [13180.1671323532 13190.1758726644;
                13204.6460491165 13213.4658144877;
            	13223.5416716495 13231.8066242708;
            	13240.8254880111 13248.3806143766];
t_range{3,2} = [17498.3372776238 17508.5899818657;
                17518.65471288	 17527.4336049202;
                17534.7066835018 17543.6095581598;
                17555.0381165225 17563.2550896691];

t_range{4,1} = [7750.41391900039 7760.41667279335;
                7769.79948913016 7779.25430325858;
                7789.88297454617 7798.9378278812;
            	7811.87219566194 7820.43094156109];
t_range{4,2} = [12167.5679156448 12176.2048347109;
                12184.1898223807 12192.5287732852;
                12202.5875041247 12210.8584637797;
                12220.5252459192 12228.4103541671];

mvf_range_val = cell(4,1);
mvf_range_val{1} = zeros(8,8);
mvf_range_val{2} = zeros(8,8);
mvf_range_val{3} = zeros(8,8);
mvf_range_val{4} = zeros(8,8);

% take an exple 
for subj_i = 1:4
    for ss_i = 1:2 % 2 sessions (sitting positions)
        % 1. check the time
%         SessionScanEMG(ss_num{subj_i}(ss_i), 1)
        sstmp = SessionScan(ss_num{subj_i}(ss_i));
        close all;
        fh = sstmp.plotForceEMG();
        axh = fh.Children;
        subplot(axh(1));
        plot(sstmp.data.t,sstmp.data.emgevl(7,:), 'color', 'k', 'linewidth', 2);
        plot(sstmp.data.t,-sstmp.data.emgevl(8,:), 'color', 'k', 'linewidth', 2);
        xline(t_range{subj_i,ss_i}(:));
        xlabel('time (s)');
        ylabel('EMG'); % or mV
        
%         ylim([-200 200]);
        legend('Flexor', 'Extensor');
        subplot(axh(2));
        xline(t_range{subj_i,ss_i}(:));
        plot(sstmp.data.t,sstmp.data.emgevl(5,:), 'color', 'k', 'linewidth', 2);
        plot(sstmp.data.t,-sstmp.data.emgevl(6,:), 'color', 'k', 'linewidth', 2);
        xlabel('time (s)');
        ylabel('EMG');
        
%         ylim([-200 200]);
        subplot(axh(3));
        xline(t_range{subj_i,ss_i}(:));
        plot(sstmp.data.t,sstmp.data.emgevl(3,:), 'color', 'k', 'linewidth', 2);
        plot(sstmp.data.t,-sstmp.data.emgevl(4,:), 'color', 'k', 'linewidth', 2);
        xlabel('time (s)');
        ylabel('EMG');
%         ylim([-200 200]);
        subplot(axh(4));
        xline(t_range{subj_i,ss_i}(:));
        plot(sstmp.data.t,sstmp.data.emgevl(1,:), 'color', 'k', 'linewidth', 2);
        plot(sstmp.data.t,-sstmp.data.emgevl(2,:), 'color', 'k', 'linewidth', 2);
        xlabel('time (s)');
        ylabel('EMG');
%         ylim([-200 200]);
        subplot(axh(5));
        xline(t_range{subj_i,ss_i}(:));
        xlabel('time (s)');
        ylabel('Force (N)');
        
        % 2. get the value put into mvf_range_val
        for ch_i=1:8
            for rep_i = 1:4
            t_idx = sstmp.data.t >= t_range{subj_i,ss_i}(rep_i,1) & ...
                sstmp.data.t < t_range{subj_i,ss_i}(rep_i,2) ;
            dat_tmp = sstmp.data.emgevl(ch_i,t_idx); 
            mvf_range_val{subj_i}(ch_i,rep_i+(ss_i-1)*4) = mean(dat_tmp);
            end
        end

    end
end

% display every subject MVF_range 
t_mvf = [];
for subj_i = 1:4
    t_mvf_tmp = max(mvf_range_val{subj_i}, [], 2);
    t_mvf = [t_mvf; t_mvf_tmp']; % each row new subject; each column new muscle
end

format long 
t_mvf

%% an exmple to show MVF

sstmp = SessionScan(4315);

fh = sstmp.plotForceEMGEVL();
axh = fh.Children();

% xlim([1.82 1.83] * 1e4);
xlim([440 530]);

subplot(axh(1)); title('Shoulder'); 
legend('flexor', 'extensor');

subplot(axh(2)); title('Deltoid'); 

subplot(axh(3)); title('Elbow'); 

subplot(axh(4)); title('Wrist'); 

sgtitle('MVF measurement and EMG');

% figure 1, the MVF vs the filtered EMG
% figure(1); 
% subplot(2,1,1); % force 
% 
% subplot(2,1,2); % EMG 
% figure 2, the filtered EMG and the rtf EMG

% figure 3, the force of each subject 


