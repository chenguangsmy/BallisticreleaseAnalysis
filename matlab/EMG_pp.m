clear, clc, close all

% LOAD THE DATA
% load('ss4310_4356_EMG.mat'); %4Subject Data
load('ss4379_4438_OPTupdate.mat'); %6Subject Data (only first 4 have EMG)

%3Subject Multiple Direction - Thanksgiving 2022 - Continiue of Human
%Experimentation Data (Subj 7-9)
% load('ss4446_4467.mat');

%3Subject Multiple Direction - Early December 2022 - Continiue of Human
%Experimentation Data (Subj 10-12)
% load('ss4472_4495.mat');

%3Subject Multiple Direction - Mid December 2022 - Continiue of Human
%Experimentation Data (Subj 13-15)

% load('ss4500_4524.mat');

%3Subject Multiple Direction - Early January 2023 - Continiue of Human
%Experimentation Data (Subj 16-18)
% load('ss4530_4563.mat');

%2Subject Multiple Direction - Late January 2023 - Continiue of Human
%Experimentation Data (Subj 19-20)
% load('ss4573_4587.mat');

%%
clearvars -except data data_index_ss data_index_tr
%%
Data = data;
Freq = 2e3;
t_step = 1/Freq;

r = size(Data, 1); % subj
c = size(Data, 2); % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
% 1 displacement
% 2 force
% 3 force command
% 4 velocity
% 5 torque of 3rd joint
% 6 displacement (vector length in cartesian space)
% 7 force (vector length in cartesian space)
% 8 EMG

% In this for loop it is defined a variable called idx_t that trims the
% data only over the interested piece (the ballistic release)
% There are two options for it now, one for pure unperturbed ballistic
% release and another for unpert+pulses.
% REMEMBER to switch between the two to correctly analyze data!!!!

for ri = 1:4 % subj
    for ci = 1:c
        for fi = 1:f % fce
            for di = 1:d % target distance
                for li = 1:1 % perturbation
                    trial_num = length(Data(ri,ci,fi,di,:,li));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,ci,fi,di,ti,li}))
                            continue;
                        end

                        if (isnan(prod(Data{ri,ci,fi,di,ti,li}.ts)) == 1)
                            F = Data{ri,ci,fi,di,ti,li}.f(1,:);
                            Fd = diff(F);
                            idx_i = find(Fd == min(Fd));
                            idx = (idx_i):(idx_i+3e3);
                        else
                            idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
                        end
                        
                        %For Pulses t_in = release - 1.2s
                        %                                 idx = (idx(1)-600):idx(end);

                        %For Multiple t_in = release - 0.5s
                        idx = (idx(1)-1000):idx(end);

                        %For Analysis of Force Hold t_in = release -1s
                        %                                 idx = (idx(1)-500):idx(end);
                        %idx = (idx(1)):idx(end);

                        time = t_step*(idx-idx(1));
                        idx_t{ri,ci,fi,di,ti,li} = idx;
                        time_t{ri,ci,fi,di,ti,li} = t_step*(idx-idx(1));

                    end
                end
            end
        end
    end
end


%% Plot Data
clc, close all

set(0, 'DefaultLineLineWidth', 2);
set(groot,'defaultAxesFontSize',16);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background
set(groot,'defaultAxesBox','on'); % box on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');


%Select Subject and Direction
subj = 1;
dir = 1;
norm = 1; % 1 --- Absolute Magnitude
          % 2 --- Normalization wrt to Maximum Magnitude

fig = figure(); %Force-Displacement Data
set(gcf,'color','w');
sgtitle('Force-Displacement Data')
c_force = [0 0.4470 0.7410];
c_disp = [0.8500 0.3250 0.0980];
idx = 1;
for f_tg = 1:3
    for d_tg = 1:3
        subplot(3,3,idx),
        for tr = 1:9
            yyaxis left
            plot(time_t{subj,dir,f_tg,d_tg,tr,1},Data{subj,dir,f_tg,d_tg,tr,1}.f(1,idx_t{subj,dir,f_tg,d_tg,tr,1}),'-','Color',c_force), hold on
            if (dir == 1) || (dir == 2)
                ylim([0 30])
            else
                ylim([-30 0])
            end
            yyaxis right
            plot(time_t{subj,dir,f_tg,d_tg,tr,1},Data{subj,dir,f_tg,d_tg,tr,1}.ox(1,idx_t{subj,dir,f_tg,d_tg,tr,1})-mean(Data{subj,dir,f_tg,d_tg,tr,1}.ox(1,idx_t{subj,dir,f_tg,d_tg,tr,1}(1:950))),'-','Color',c_disp), hold on
            if (dir == 1) || (dir == 2)
                ylim([0 0.1])
            else
                ylim([-0.1 0])
            end
        end
        idx = idx + 1;
    end
end
% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Force [N] - Displacement [m]');
xlabel(han,'Time [s]');


%Filter Infinite Impulse Response
% EMG_hp = designfilt('highpassfir', ...       % Response type
%     'StopbandFrequency',2, ...    % Frequency constraints [Hz]
%     'PassbandFrequency',6, ...
%     'StopbandAttenuation',65, ...   % Magnitude constraints [dB]
%     'PassbandRipple',0.3, ...
%     'SampleRate',2000);               % Sample rate [Hz]
% fvtool(EMG_hp)

% EMG_lp = designfilt('lowpassfir', ...       % Response type
%     'StopbandFrequency',250, ...    % Frequency constraints [Hz]
%     'PassbandFrequency',200, ...
%     'StopbandAttenuation',40, ...   % Magnitude constraints [dB]
%     'PassbandRipple',0.1, ...
%     'SampleRate',2000);               % Sample rate [Hz]
% % fvtool(EMG_lp)


% EMG PostProcessing
for subj = 1:4
    for dir = 1:4
        for d_tg = 1:3
            for f_tg = 1:3
                for tr = 1:9
                    for emg_ch = 1:8
                        emg_raw = Data{subj,dir,f_tg,d_tg,tr,1}.emg(emg_ch,idx_t{subj,dir,f_tg,d_tg,tr,1});
                        time_raw = time_t{subj,dir,f_tg,d_tg,tr,1};

                        emg_mean = mean(emg_raw);
                        if (isnan(emg_mean))
                            emg_rect = 0;
                            EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch} = zeros(1,length(time_raw));
                            EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch} = zeros(1,length(time_raw));
                            EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch} = zeros(1,length(time_raw));
                            nan_test{subj,dir,f_tg,d_tg,tr,emg_ch} = 1;
                        else
                            emg_rect = abs(emg_raw-mean(emg_raw));
                            if norm == 1
                                emg_rect_c = emg_rect;
                            elseif norm == 2
                                emg_rect_c = emg_rect/max(emg_rect);
                            end
                            EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch} = emg_rect_c;
                            EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch} = movmean(emg_rect_c,250); %250 = 125 ms for ENVELOPE %20 = 10 ms for CUSUM
                            EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch} = cumsum(EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch}-mean(EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch}));
                            nan_test{subj,dir,f_tg,d_tg,tr,emg_ch} = 0;
                        end

                    end
                end

            end
        end
    end
end

for subj = 1:4
    for dir = 1:4
        for d_tg = 1:3
            for f_tg = 1:3
                for emg_ch = 1:8
                    for tr = 1:9
                        L(tr) = length(EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch});
                    end
                    Lmin = min(L);
                    for tr = 1:9
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,tr,emg_ch} = EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch}(1:Lmin);
                    end
                    EMG_mv_cr_lu{subj,dir,f_tg,d_tg,emg_ch} = [...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,1,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,2,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,3,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,4,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,5,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,6,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,7,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,8,emg_ch};...
                        EMG_mv_cropped{subj,dir,f_tg,d_tg,9,emg_ch}];
                    EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch} = mean(EMG_mv_cr_lu{subj,dir,f_tg,d_tg,emg_ch});
                    EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch} = cumsum(EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch}-mean(EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch}));
                end
            end
        end
    end
end

%% PLOTS
tr_max = 9;
clc, close all

% %Select Subject
% subj = 1;
% %Select Direction
% dir = 4;


for subj = 1:1
    for dir = 1:4
%         clc, close all
        if dir == 1
            dirz = '+x';
        elseif dir == 2
            dirz = '+y';
        elseif dir == 3
            dirz = '-x';
        elseif dir == 4
            dirz = '-y';
        end

        switch norm

            case 1 %Absolute Magnitude

                % Figure - Option 1 - Absolute Magnitude
                triplet = [0.3010 0.7450 0.9330];
                triplet_f = [0 0.4470 0.7410];

                triplet2 = [0.8500 0.3250 0.0980];
                triplet2_f = [0.6350 0.0780 0.1840];

                for emg_ch = 1:2:8
                    fig2 = figure(); %EMG Post-Processing
                    set(gcf,'color','w');
                    if emg_ch == 1
                        sgtitle('Flexor vs Extensor Carpii Radialis')
                    elseif emg_ch == 3
                        sgtitle('Biceps Short vs Triceps Long Head')
                    elseif emg_ch == 5
                        sgtitle('Deltoid Anterior vs Posterior')
                    elseif emg_ch == 7
                        sgtitle('Pectoralis vs Trapezius')
                    end
                    idx = 1;
                    for f_tg = 1:3
                        for d_tg = 1:3
                            subplot(3,6,idx),
                            for tr = 1:tr_max
                                plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet,'LineWidth',0.1), hold on
                            end
                            for tr = 1:tr_max
                                plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet_f,'LineWidth',2), hold on
                                plot([0.5 0.5],[0 max(EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch})],'r','LineWidth',0.2), hold on
                            end
                            plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch})),EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch},'r','LineWidth',2), hold on
                            ylim([0 0.2])

                            subplot(3,6,idx+1),
                            for tr = 1:tr_max
                                plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch+1},'Color',triplet2,'LineWidth',0.1), hold on
                            end
                            for tr = 1:tr_max
                                plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch+1},'Color',triplet2_f,'LineWidth',2), hold on
                                plot([0.5 0.5],[0 max(EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch+1})],'b','LineWidth',0.2), hold on
                            end
                            plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch+1})),EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch+1},'b','LineWidth',2), hold on
                            ylim([0 0.2])

                            idx = idx + 2;
                        end
                    end
                    % Give common xlabel, ylabel and title to your figure
                    han2=axes(fig2,'visible','off');
                    han2.Title.Visible='on';
                    han2.XLabel.Visible='on';
                    han2.YLabel.Visible='on';
                    ylabel(han2,'EMG Signal [mV]');
                    xlabel(han2,'Time [s]');
                end


                % Figure - Cumulative EMG - Absolute Magnitude
                triplet = [0.3010 0.7450 0.9330];
                triplet_f = [0 0.4470 0.7410];

                triplet2 = [0.8500 0.3250 0.0980];
                triplet2_f = [0.6350 0.0780 0.1840];


%                 subjdir = figure(); axis off, hold on, title(strcat('Subject ',num2str(subj+18)),'FontSize',26), subtitle(strcat('Direction ',dirz),'FontSize',24);
% %                 exportgraphics(subjdir,'EMG_pp.pdf','Resolution',600,'Append',true);
%                 for emg_ch = 1:2:8
%                     fig2 = figure(); %EMG Post-Processing
%                     set(gcf,'color','w');
%                     if emg_ch == 1
%                         sgtitle('Flexor vs Extensor Carpii Radialis')
%                     elseif emg_ch == 3
%                         sgtitle('Biceps Short vs Triceps Long Head')
%                     elseif emg_ch == 5
%                         sgtitle('Deltoid Anterior vs Posterior')
%                     elseif emg_ch == 7
%                         sgtitle('Pectoralis vs Trapezius')
%                     end
%                     idx = 1;
%                     for f_tg = 1:3
%                         for d_tg = 1:3
%                             subplot(3,6,idx),
%                             for tr = 1:tr_max
%                                 plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet_f,'LineWidth',1), hold on
%                                 plot([0.5 0.5],[min(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch}) max(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch})],'k--','LineWidth',1), hold on
%                             end
%                             %                     plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch})),EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch},'r','LineWidth',2), hold on
% 
%                             %                     ylim([0 20])
% 
%                             subplot(3,6,idx+1),
%                             for tr = 1:tr_max
%                                 plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch+1},'Color',triplet2_f,'LineWidth',1), hold on
%                                 plot([0.5 0.5],[min(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch+1}) max(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch+1})],'k--','LineWidth',1), hold on
%                             end
%                             %                     plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch+1})),EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch+1},'b','LineWidth',2), hold on
%                             %                     ylim([0 20])
% 
%                             idx = idx + 2;
%                         end
%                     end
%                     % Give common xlabel, ylabel and title to your figure
%                     han2=axes(fig2,'visible','off');
%                     han2.Title.Visible='on';
%                     han2.XLabel.Visible='on';
%                     han2.YLabel.Visible='on';
%                     ylabel(han2,'Cum. EMG Signal [mV$\cdot$s]');
%                     xlabel(han2,'Time [s]');
% 
%                     emg_plot = gcf;
% %                     exportgraphics(emg_plot,'EMG_pp.pdf','Resolution',600,'Append',true);
%                 end



            case 2 %Normalized wrt Average Rectified

                % Figure - EMG - Normalized Magnitude with respect to Rectified Mean
                triplet = [0.3010 0.7450 0.9330];
                triplet_f = [0 0.4470 0.7410];

                triplet2 = [0.8500 0.3250 0.0980];
                triplet2_f = [0.6350 0.0780 0.1840];

                %         for emg_ch = 1:2:8
                %             fig2 = figure(); %EMG Post-Processing
                %             set(gcf,'color','w');
                %             if emg_ch == 1
                %                 sgtitle('Flexor vs Extensor Carpii Radialis')
                %             elseif emg_ch == 3
                %                 sgtitle('Biceps Short vs Triceps Long Head')
                %             elseif emg_ch == 5
                %                 sgtitle('Deltoid Anterior vs Posterior')
                %             elseif emg_ch == 7
                %                 sgtitle('Pectoralis vs Trapezius')
                %             end
                %             idx = 1;
                %             for f_tg = 1:3
                %                 for d_tg = 1:3
                %                     subplot(3,6,idx),
                %                     for tr = 1:tr_max
                %                         plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet,'LineWidth',0.1), hold on
                %                     end
                %                     for tr = 1:tr_max
                %                         plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet_f,'LineWidth',2), hold on
                %                         plot([0.5 0.5],[0 max(EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch})],'r','LineWidth',0.2), hold on
                %                     end
                %                     plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch})),EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch},'r','LineWidth',2), hold on
                %                     ylim([0 1])
                %
                %                     subplot(3,6,idx+1),
                %                     for tr = 1:tr_max
                %                         plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch+1},'Color',triplet2,'LineWidth',0.1), hold on
                %                     end
                %                     for tr = 1:tr_max
                %                         plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGfilt_mv{subj,dir,f_tg,d_tg,tr,emg_ch+1},'Color',triplet2_f,'LineWidth',2), hold on
                %                         plot([0.5 0.5],[0 max(EMGpp{subj,dir,f_tg,d_tg,tr,emg_ch+1})],'b','LineWidth',0.2), hold on
                %                     end
                %                     plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch+1})),EMG_mv_avg{subj,dir,f_tg,d_tg,emg_ch+1},'b','LineWidth',2), hold on
                %                     ylim([0 1])
                %
                %                     idx = idx + 2;
                %                 end
                %             end
                %             % Give common xlabel, ylabel and title to your figure
                %             han2=axes(fig2,'visible','off');
                %             han2.Title.Visible='on';
                %             han2.XLabel.Visible='on';
                %             han2.YLabel.Visible='on';
                %             ylabel(han2,'Norm. EMG Signal [-]');
                %             xlabel(han2,'Time [s]');
                %         end

                % Figure - Cumulative EMG - Normalized Magnitude with respect to Rectified Mean
                triplet = [0.3010 0.7450 0.9330];
                triplet_f = [0 0.4470 0.7410];

                triplet2 = [0.8500 0.3250 0.0980];
                triplet2_f = [0.6350 0.0780 0.1840];

                for emg_ch = 1:2:8
                    fig2 = figure(); %EMG Post-Processing
                    set(gcf,'color','w');
                    if emg_ch == 1
                        sgtitle('Flexor vs Extensor Carpii Radialis')
                    elseif emg_ch == 3
                        sgtitle('Biceps Short vs Triceps Long Head')
                    elseif emg_ch == 5
                        sgtitle('Deltoid Anterior vs Posterior')
                    elseif emg_ch == 7
                        sgtitle('Pectoralis vs Trapezius')
                    end
                    idx = 1;
                    for f_tg = 1:3
                        for d_tg = 1:3
                            subplot(3,6,idx),
                            for tr = 1:tr_max
                                plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet_f,'LineWidth',1), hold on
                                plot([0.5 0.5],[min(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch}) max(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch})],'k--','LineWidth',1), hold on
                            end
                            plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch})),EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch},'r','LineWidth',1), hold on
                            ylim([-60 100])

                            subplot(3,6,idx+1),
                            for tr = 1:tr_max
                                plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch+1},'Color',triplet2_f,'LineWidth',1), hold on
                                plot([0.5 0.5],[min(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch+1}) max(EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch+1})],'k--','LineWidth',1), hold on
                            end
                            plot(time_t{subj,dir,f_tg,d_tg,tr,1}(1:length(EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch+1})),EMG_cs_avg{subj,dir,f_tg,d_tg,emg_ch+1},'b','LineWidth',1), hold on
                            ylim([-60 100])

                            idx = idx + 2;
                        end
                    end
                    % Give common xlabel, ylabel and title to your figure
                    han2=axes(fig2,'visible','off');
                    han2.Title.Visible='on';
                    han2.XLabel.Visible='on';
                    han2.YLabel.Visible='on';
                    ylabel(han2,'Norm. Cum. EMG Signal [s]');
                    xlabel(han2,'Time [s]');
                end
        end

    end
end

%% REFLEX PICKER

clc, close all
idx = 1;
triplet_f = [0 0.4470 0.7410];

%Training Data Size
subj_max = 4;
dir_max = 4;
d_tg_max = 3;
f_tg_max = 3;
emg_ch_max = 7;
tr_max = 3;

Ntr = subj_max*dir_max*d_tg_max*f_tg_max*emg_ch_max*tr_max;

%%

for subj = 1:subj_max%4
    for dir = 1:dir_max%4
        for d_tg = 1:d_tg_max%3
            for f_tg = 1:f_tg_max%3
                for emg_ch = 1:emg_ch_max%7 %(7 = no trapezius)
                    for tr = 1:tr_max %used 4-9 as validation data-set
                        %Iteration Counter Percentage
                        percentage = 100*(idx/Ntr)
                        time_matrix(idx,:) = time_t{subj,dir,f_tg,d_tg,tr,1}(1:1500);
                        EMGcusum_matrix(idx,:) = EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch}(1:1500);
                        figure(),
                        subplot(2,1,1),
                        plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet_f,'LineWidth',2), hold on
                        plot(0.5,0,'rx','MarkerSize',12,'LineWidth',2), grid on
                        ylabel('CuSum EMG [mV$\cdot$s]');
                        subplot(2,1,2),
                        plot(time_t{subj,dir,f_tg,d_tg,tr,1},EMGcs{subj,dir,f_tg,d_tg,tr,emg_ch},'Color',triplet_f,'LineWidth',2), grid on
%                         plot(0.5,0,'r','MarkerSize',12,'LineWidth',2), grid on
                        xlabel('Time [s]')
                        ylabel('CuSum EMG [mV$\cdot$s]');
                        xlim([0.5 0.7])
                        [x,y,button] = ginput(1);
                        Subject(idx,1) = subj;
                        Direction(idx,1) = dir;
                        Force(idx,1) = f_tg;
                        Displacement(idx,1) = d_tg;
                        Trial(idx,1) = tr;
                        if button == 1 %Press Left-Mouse Button
                            ReflexCheck(idx,1) = "Evident";
                            ReflexOnset(idx,1) = x;
                        elseif button == 2 %Press Central-Mouse Button
                            ReflexCheck(idx,1) = "Unsure";
                            ReflexOnset(idx,1) = x;
                        elseif button == 3 %Press Right-Mouse Button
                            ReflexCheck(idx,1) = "Non-Evident";
                            ReflexOnset(idx,1) = NaN;
                        else
                            ReflexCheck(idx,1) = "Non-Sampled";
                            ReflexOnset(idx,1) = NaN;
                        end

                        if emg_ch == 1
                            EMGch(idx,1) = "Flexor Radialis";
                        elseif emg_ch == 2
                            EMGch(idx,1) = "Extensor Radialis";
                        elseif emg_ch == 3
                            EMGch(idx,1) = "Biceps";
                        elseif emg_ch == 4
                            EMGch(idx,1) = "Triceps";
                        elseif emg_ch == 5
                            EMGch(idx,1) = "Deltoid Anterior";
                        elseif emg_ch == 6
                            EMGch(idx,1) = "Deltoid Posterior";
                        elseif emg_ch == 7
                            EMGch(idx,1) = "Pectoralis";
                        elseif emg_ch == 8
                            EMGch(idx,1) = "Trapezius";
                        end
                        idx = idx+1;
                        close all

                        if button == 48 %It is 0 on the keypad
                            break;
                        end
                    end
                    if button == 48
                        break;
                    end
                end
                if button == 48
                    break;
                end
            end
            if button == 48
                break;
            end
        end
        if button == 48
            break;
        end
    end
    if button == 48
        break;
    end
end

reflex = table(Subject,Direction,Force,Displacement,EMGch,Trial,ReflexCheck,ReflexOnset)
%%
idx = 1;
for subj = 1:subj_max%4
    for dir = 1:dir_max%4
        for d_tg = 1:d_tg_max%3
            for f_tg = 1:f_tg_max%3
                for emg_ch = 1:emg_ch_max%7 %(7 = no trapezius)
                    for tr = 1:tr_max %used 4-9 as validation data-set

                        EMGdata_cs{idx,1} = EMGcusum_matrix(idx,:);
                        idx = idx+1;
                        close all

                    end
                end
            end
        end
    end
end
%% Results Table
reflex_results = table(Subject,Direction,Force,Displacement,EMGch,Trial,ReflexOnset,ReflexCheck,EMGcusum_matrix);
% SHUFFLED MATRIX
reflex_res_rand = reflex_results(randperm(size(reflex_results ,1)), :);