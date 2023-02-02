clc, clear, close all
%Optimal Figure Setting to use in any script
set(0,'defaultaxesfontname','Times New Roman');
set(0,'defaulttextfontname','Times New Roman');
set(0,'defaultaxesfontsize',8); % 8 for paper images 12 for normal images
set(0,'defaulttextfontsize',8); % 8 for paper images 12 for normal images

set(gcf,'color','w');

%Rika Figure Options
set(0, 'DefaultLineLineWidth', 1);
set(groot,'defaultAxesFontSize',16);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background
set(groot,'defaultAxesBox','on'); % box on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

% LOAD THE DATA
%4Subject Multiple Direction - Summer 2022 - Not Human Experimentation Data
% load('ss4310_4356.mat');

%6Subject Multiple Direction - Fall 2022 - Beginning of Human
%Experimentation Data (Subj 1-6)
% load('ss4379_4438_OPTupdate.mat');

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
load('ss4573_4587.mat');

%Check on Pulse with Ballistic Models
% load('ss4253_4274.mat')
%%

clearvars -except data data_index_ss data_index_tr subjprop

%%
clc, close all
Data = data;
Freq = 2000;
% Freq = 200; %old data
t_step = 1/Freq;

r = size(Data, 1); % subj
c = size(Data, 2); % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
% if_subtract = 0;
% epoc_type = 2;
% plot_type = 2;          % 1 displacement
%                         % 2 force
%                         % 3 force command
%                         % 4 velocity
%                         % 5 torque of 3rd joint
%                         % 6 displacement (vector length in cartesian space)
%                         % 7 force (vector length in cartesian space)
% axh = zeros(f,r);

for ri = 1:r % subj
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
%                         idx = (idx(1)-600):idx(end);

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


%% OLD 
% % In this for loop it is defined a variable called idx_t that trims the
% % data only over the interested piece (the ballistic release)
% % There are two options for it now, one for pure unperturbed ballistic
% % release and another for unpert+pulses.
% % REMEMBER to switch between the two to correctly analyze data!!!!
% 
% for ri = 1:r % subj
%     for ci = 1:c
%         for fi = 1:f % fce
%             for di = 1:d % target distance
%                 axh(ri, fi) = subplot(d,f,d*(fi-1) + di);grid on;hold on;
%                 %for di = 3 % target distance
%                 for li = 1:p % perturbation
%                     trial_num = length(Data(ri,ci,fi,di,:,li));
%                     for ti = 1:trial_num % each trial
%                         if (isempty(Data{ri,ci,fi,di,ti,li}))
%                             continue;
%                         end
% 
%                         switch epoc_type
%                             case 1
%                                 idx = find(Data{ri,ci,fi,di,ti,li}.Fp(2,:)~=0 & ...
%                                     Data{ri,ci,fi,di,ti,li}.ts==4);  % pert at y
%                                 idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
%                                 if li == 1
%                                     disp('ERROR: should use li == 2!!!');
%                                 end
%                             case 2
%                                 idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
%                               
%                                 %For Multiple t_in = release - 0.5s
%                                 idx = (idx(1)-1000):idx(end); %1000 @2 kHz // 250 @500 Hz
% 
%                                 %For Analysis of Force Hold t_in = release -1s 
% %                                 idx = (idx(1)-500):idx(end);
%                                 %idx = (idx(1)):idx(end);
%                         end
% 
%                         time = t_step*(idx-idx(1));
%                         idx_t{ri,ci,fi,di,ti,li} = idx;
%                         time_t{ri,ci,fi,di,ti,li} = t_step*(idx-idx(1));
%                         %time = idx-idx(1);
%                         switch plot_type
%                             case 1
%                                 dat = Data{ri,ci,fi,di,ti,li}.ox(1,idx);
%                                 %dat = dat - dat(1);
%                                 titlestr = 'displacement';
%                             case 2
%                                 dat = Data{ri,ci,fi,di,ti,li}.f(1,idx);
%                                 titlestr = 'force';
%                             case 3
%                                 dat = Data{ri,ci,fi,di,ti,li}.Fp(1,idx);
%                                 titlestr = 'Fp';
%                             case 4
%                                 dat = Data{ri,ci,fi,di,ti,li}.v(1,idx);
%                                 titlestr = 'velocity';
%                             case 5
%                                 dat = Data{ri,ci,fi,di,ti,li}.tq(3,idx);
%                                 titlestr = 'torque3';
%                             case 6
%                                 dat = Data{ri,ci,fi,di,ti,li}.x(:,idx);
%                                 dat_submean = dat - mean(dat(:,1:50),2);
%                                 dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
%                                 dat = dat_norm .* sign(dat_submean(2,:));
%                                 titlestr = 'norm displacement';
%                             case 7 % the force mode
%                                 dat = Data{ri,ci,fi,di,ti,li}.f(:,idx);
%                                 dat_submean = dat - mean(dat(:,1:50),2);
%                                 dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
%                                 dat = dat_norm .* sign(dat_submean(2,:));
%                                 titlestr = 'norm force';
% 
%                         end
%                         if (if_subtract)
%                             dat = dat - mean(dat(1:50));
%                         end
%                         plot(time, dat, 'Color', colors(4*(li-1)+di, :));
%                         %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
%                     end
%                 end
%             end
%         end
%     end
% end
% %xlim([0 0.7])
% linkaxes(axh, 'xy');
% %xlim([0 0.5]);
% % xlim([0 2])
% sgtitle(titlestr);

%% System Identification

clc, close all
%Select Subject
% subj = 1;
% Select Direction (1=+x,2=+y,3=-x,4=-y)
% dir = 1;

%Here you have the functions to perfrom system identification depending 
%on the kind of test you performed: Unpeturbed+Pulses or UnperturbedOnly

%------Function Computing system Identification Unpert+Pulses (remember to
%change idx_t above)
% results = sys_id_pp(Data,idx_t,time_t,subj);

%------Function Computing system Identification Unpert Only (remember to
%change idx_t above)
% results = sys_id_ubr(Data,idx_t,time_t,subj,dir);

%The following function is to compute the kinematic analysis based on the Optotrak
%markers on hand, elbow and shoulder
%(As far as I have seen now, it works well only for direction "+x" e
%subject Chenguamng)

%------Function Computing Arm Kinematic Analysis (Joint Angles, Jacobian, Mass Matrix)
% results = sys_kin(Data,idx_t,time_t,subj,dir);

%------Function Computing PSD of Force Hold
% results = sys_force_hold_analysis(Data,idx_t,time_t,subj,dir);

%------Imp. Ident and Kin. Analysis of All Subject and Directions
for subj = 1:2
    for dir = 1:4
        BMI{subj} = subjprop.weight(subj)/((0.01*subjprop.height(subj))^2);
        results{subj+18,dir} = sys_id_ubr(Data,idx_t,time_t,subj,dir,subjprop);
        resultsK{subj,dir} = sys_kin(Data,idx_t,time_t,subj,dir,subjprop);
    end
    subj+18
end
%% Specific Data Removals for Anomalies

%Subject 10, All Directions, Force 3, Displacement 1 --> not following
%experimental protocol, remove estimates

subj = 10;
for dir = 1:4
    results{subj,dir}.K_up_avg(3,1) = NaN;
    results{subj,dir}.K_up_std(3,1) = NaN;
    results{subj,dir}.B_up_avg(3,1) = NaN;
    results{subj,dir}.B_up_std(3,1) = NaN;
    results{subj,dir}.M_up_avg(3,1) = NaN;
    results{subj,dir}.M_up_std(3,1) = NaN;
    results{subj,dir}.FIT_up_avg(3,1) = NaN;
    results{subj,dir}.FIT_up_std(3,1) = NaN;
    results{subj,dir}.up_omegan_avg(3,1) = NaN;
    results{subj,dir}.up_omegan_std(3,1) = NaN;
    results{subj,dir}.up_dampr_avg(3,1) = NaN;
    results{subj,dir}.up_dampr_std(3,1) = NaN;

    results{subj,dir}.K_fr_avg(3,1) = NaN;
    results{subj,dir}.B_fr_avg(3,1) = NaN;
    results{subj,dir}.M_fr_avg(3,1) = NaN;
    results{subj,dir}.Fc_fr_avg(3,1) = NaN;
    results{subj,dir}.K_fr_std(3,1) = NaN;
    results{subj,dir}.B_fr_std(3,1) = NaN;
    results{subj,dir}.M_fr_std(3,1) = NaN;
    results{subj,dir}.Fc_fr_std(3,1) = NaN;
end

%Subject 14, Direction +/- y, All Forces and displacement are bi-phasic (complete wrong mass estimates) --> not following
%experimental protocol, remove estimates

subj = 14;
for dir = [2,4]
    results{subj,dir}.K_up_avg(:,:) = NaN;
    results{subj,dir}.K_up_std(:,:) = NaN;
    results{subj,dir}.B_up_avg(:,:) = NaN;
    results{subj,dir}.B_up_std(:,:) = NaN;
    results{subj,dir}.M_up_avg(:,:) = NaN;
    results{subj,dir}.M_up_std(:,:) = NaN;
    results{subj,dir}.FIT_up_avg(:,:) = NaN;
    results{subj,dir}.FIT_up_std(:,:) = NaN;
    results{subj,dir}.up_omegan_avg(:,:) = NaN;
    results{subj,dir}.up_omegan_std(:,:) = NaN;
    results{subj,dir}.up_dampr_avg(:,:) = NaN;
    results{subj,dir}.up_dampr_std(:,:) = NaN;

    results{subj,dir}.K_fr_avg(:,:) = NaN;
    results{subj,dir}.B_fr_avg(:,:) = NaN;
    results{subj,dir}.M_fr_avg(:,:) = NaN;
    results{subj,dir}.Fc_fr_avg(:,:) = NaN;
    results{subj,dir}.K_fr_std(:,:) = NaN;
    results{subj,dir}.B_fr_std(:,:) = NaN;
    results{subj,dir}.M_fr_std(:,:) = NaN;
    results{subj,dir}.Fc_fr_std(:,:) = NaN;
end
%% SUBJEECT GROUPING
clc, close all

kk = [];
bb = [];
mm = [];
oo = [];
dd = [];

for dir = 1:4
    for f_sel = 1:3
        for d_sel = 1:3
            for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]
                if subj == 10 && f_sel == 3 && d_sel == 1
                    kk = [kk,NaN(1,length(results{subj,dir}.K_up_tr{f_sel,d_sel}))];
                    bb = [bb,NaN(1,length(results{subj,dir}.B_up_tr{f_sel,d_sel}))];
                    mm = [mm,NaN(1,length(100*results{subj,dir}.M_up_tr{f_sel,d_sel}/results{subj,dir}.subj_mass))];
                    oo = [oo,NaN(1,length(results{subj,dir}.up_omegan_tr{f_sel,d_sel}))];
                    dd = [dd,NaN(1,length(results{subj,dir}.up_dampr_tr{f_sel,d_sel}))];
                elseif subj == 14 && (dir == 2 || dir == 4)
                    kk = [kk,NaN(1,length(results{subj,dir}.K_up_tr{f_sel,d_sel}))];
                    bb = [bb,NaN(1,length(results{subj,dir}.B_up_tr{f_sel,d_sel}))];
                    mm = [mm,NaN(1,length(100*results{subj,dir}.M_up_tr{f_sel,d_sel}/results{subj,dir}.subj_mass))];
                    oo = [oo,NaN(1,length(results{subj,dir}.up_omegan_tr{f_sel,d_sel}))];
                    dd = [dd,NaN(1,length(results{subj,dir}.up_dampr_tr{f_sel,d_sel}))];
                else 
                    kk = [kk,results{subj,dir}.K_up_tr{f_sel,d_sel}];
                    bb = [bb,results{subj,dir}.B_up_tr{f_sel,d_sel}];
                    mm = [mm,100*results{subj,dir}.M_up_tr{f_sel,d_sel}/results{subj,dir}.subj_mass];
                    oo = [oo,results{subj,dir}.up_omegan_tr{f_sel,d_sel}];
                    dd = [dd,results{subj,dir}.up_dampr_tr{f_sel,d_sel}];
                end
            end
            imp_tot{dir}.K_tot{f_sel,d_sel} = kk;
            imp_tot{dir}.B_tot{f_sel,d_sel} = bb;
            imp_tot{dir}.M_tot{f_sel,d_sel} = mm;
            imp_tot{dir}.omega_tot{f_sel,d_sel} = oo;
            imp_tot{dir}.damp_tot{f_sel,d_sel} = dd;

            imp_tot{dir}.K_tot_avg(f_sel,d_sel) = nanmean(kk);
            imp_tot{dir}.B_tot_avg(f_sel,d_sel) = nanmean(bb);
            imp_tot{dir}.M_tot_avg(f_sel,d_sel) = nanmean(mm);
            imp_tot{dir}.omega_tot_avg(f_sel,d_sel) = nanmean(oo);
            imp_tot{dir}.damp_tot_avg(f_sel,d_sel) = nanmean(dd);


            imp_tot{dir}.K_tot_std(f_sel,d_sel) = nanstd(kk);
            imp_tot{dir}.B_tot_std(f_sel,d_sel) = nanstd(bb);
            imp_tot{dir}.M_tot_std(f_sel,d_sel) = nanstd(mm);
            imp_tot{dir}.omega_tot_std(f_sel,d_sel) = nanstd(oo);
            imp_tot{dir}.damp_tot_std(f_sel,d_sel) = nanstd(dd);
            kk = [];
            bb = [];
            mm = [];
            oo = [];
            dd = [];
        end
    end
end

%% PLOTTING (the rest of the code concerns only results plotting)

% IN the following subsections you can find the different plots for the
% analyzed data, run the specific subsection depending on what you are
% interested in:
% - Time plots;
% - Frequency plots;
% - Bar plots (for Unperturbed only)
% - Kinematic plots
% After the plotting you have also a little piece of code to perform ANOVA
% on the unperturbed ballistic release data

%% OVERALL FIT PERFOMANCE TABLE
clc
for dir = 1:4
    idx = 1;
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        FIT_avg_s(:,:,idx) = results{subj,dir}.FIT_up_avg;
        idx = idx+1;
    end

    FIT_avg(:,:,dir) = nanmean(FIT_avg_s,3);
    FIT_std(:,:,dir) = nanstd(FIT_avg_s,[],3);
end


%% TIME PLOTS

clc, close all

for subj = 20:20
    for dir = 1:4
        %Uncomment for pulses data
        % wind_start = 20; %From which time window to start showing

        % Axes Limits (POSITIVE directions)
        if (dir == 1) || (dir == 2)
            f_maxx = 30;
            f_minn = -5;
            d_minn = -0.01;
            d_maxx = +0.1;
        elseif (dir == 3) || (dir == 4)
            % Axes Limits (NEGATIVE directions)
            f_maxx = +5;
            f_minn = -30;
            d_minn = -0.1;
            d_maxx = +0.01;
        end

        fig = figure(); %Force-Displacement Data
        set(gcf,'color','w');
        sgtitle('Force-Displacement Data')
        c_force = [0 0.4470 0.7410];
        c_disp = [0.8500 0.3250 0.0980];
        idx = 1;
        for f_sel = 1:3
            for d_sel = 1:3
                subplot(3,3,idx),
                yyaxis left
                plot(results{subj,dir}.FD_UP{f_sel,d_sel}{1,1},results{subj,dir}.FD_UP{f_sel, d_sel}{2, 1},'-','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results{subj,dir}.avg_FD_UP{f_sel,d_sel}(1,:),results{subj,dir}.avg_FD_UP{f_sel,d_sel}(2,:),'b-','LineWidth',2), grid on
                ylim([f_minn f_maxx])

                yyaxis right
                plot(results{subj,dir}.FD_UP{f_sel, d_sel}{1, 1},results{subj,dir}.FD_UP{f_sel, d_sel}{3, 1},'-','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results{subj,dir}.avg_FD_UP{f_sel, d_sel}(1,:),results{subj,dir}.avg_FD_UP{f_sel, d_sel}(3,:),'r-','LineWidth',2), hold on
                plot(results{subj,dir}.tt_mod{f_sel, d_sel},results{subj,dir}.disp_mod{f_sel, d_sel},'g--','LineWidth',2),
                %Unperturbed 9Trials
                lgd = legend({'','','','','','','','','','','','','','','','','','','','',strcat('FIT: ',num2str(round(results{subj,dir}.FIT_up_avg(f_sel,d_sel))),'\%')},'Location','best');
                %Pulses 20 Trials
                %lgd = legend({'','','','','','','','','','','','','','','','','','','','','',strcat('FIT: ',num2str(round(results{subj,dir}.FIT_up_avg(1,1))),'%')},'Location','best');
                lgd.FontSize = 7;
                lgd.Box = 'off';
                grid on
                ylim([d_minn d_maxx])
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

    end
end
%% Force-vs-Displacement

clc, close all

for subj = 1:1
    for dir = 1:4
        %Uncomment for pulses data
        % wind_start = 20; %From which time window to start showing

        % Axes Limits (POSITIVE directions)
        if (dir == 1) || (dir == 2)
            f_maxx = 30;
            f_minn = -5;
            d_minn = -0.01;
            d_maxx = +0.1;
        elseif (dir == 3) || (dir == 4)
            % Axes Limits (NEGATIVE directions)
            f_maxx = +5;
            f_minn = -30;
            d_minn = -0.1;
            d_maxx = +0.01;
        end

        fig = figure(); %Force-Displacement Data
        set(gcf,'color','w');
        sgtitle('Force-Displacement Data')
        c_force = [0 0.4470 0.7410];
        c_disp = [0.8500 0.3250 0.0980];
        idx = 1;
        for f_sel = 1:3
            for d_sel = 1:3
                subplot(3,3,idx),
                plot(results{subj,dir}.FD_UP{f_sel, d_sel}{3, 1},results{subj,dir}.FD_UP{f_sel, d_sel}{2, 1},'-','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results{subj,dir}.avg_FD_UP{f_sel, d_sel}(3,:),results{subj,dir}.avg_FD_UP{f_sel,d_sel}(2,:),'b-','LineWidth',2), grid on
                ylim([f_minn f_maxx])

                %Unperturbed 9Trials
%                 lgd = legend({'','','','','','','','','','','','','','','','','','','','',strcat('FIT: ',num2str(round(results{subj,dir}.FIT_up_avg(f_sel,d_sel))),'\%')},'Location','best');
%                 lgd.FontSize = 7;
%                 lgd.Box = 'off';
                grid on
                xlim([d_minn d_maxx])
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

    end
end


%% FREQUENCY PLOTS
%Bode Plots
clc, close all

options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
% options.MagVisible = 'off';
options.YLim = {[-80 -40]};

subj = 20;
for dir = 1:4
    figure(),
    set(gcf,'color','w');
    sgtitle('Bode Plot - Magnitude - Identified Models')
    idx = 1;

    for f_sel = 1:3
        for d_sel = 1:3
            subplot(3,3,idx),
            bodemag(results{subj,dir}.TF.up{f_sel, d_sel},{0.1*2*pi,10*2*pi},'b',options), hold on,
            grid on
            idx = idx+1;
        end
    end
end

%%
% %PSD Plots (Need to Run sys_force_hold_analysis)
% clc, close all
% figure(),
% set(gcf,'color','w');
% sgtitle('Welch Method PSD')
% idx = 1;
% for f_sel = 1:3
%     for d_sel = 1:3
%         subplot(3,3,idx),
%         sizeV = size(results{subj,dir}.freq_psd{f_sel,d_sel});
%         for i = 1:sizeV(1)
%             semilogx(results{subj,dir}.freq_psd{f_sel,d_sel}(i,:),10*log10(results{subj,dir}.force_psd{f_sel,d_sel}(i,:)),'k','LineWidth',1), hold on
%         end
%         xlabel('f (Hz)')
%         ylabel('PSD [dB/Hz]')
%         grid on       
%         
%         idx = idx+1;
%     end
% end
% 
% %FFT Plots
% figure(),
% set(gcf,'color','w');
% sgtitle('Fast Fourier Transform')
% idx = 1;
% for f_sel = 1:3
%     for d_sel = 1:3
%         subplot(3,3,idx),
%         sizeV = size(results{subj,dir}.freq_fft{f_sel,d_sel});
%         for i = 1:sizeV(1)
%             plot(results{subj,dir}.freq_fft{f_sel,d_sel}(i,:),results{subj,dir}.force_fft{f_sel,d_sel}(i,:),'k','LineWidth',1), hold on
%         end
%         xlabel("f (Hz)")
%         ylabel("|P1(f)|")
%         grid on       
%         
%         idx = idx+1;
%     end
% end
% 
% figure(),
% set(gcf,'color','w');
% sgtitle('Fast Fourier Transform')
% idx = 1;
% for f_sel = 1:3
%     for d_sel = 1:3
%         subplot(3,3,idx),
%         sizeV = size(results{subj,dir}.freq_fft{f_sel,d_sel});
%         for i = 1:sizeV(1)
%             plot(results{subj,dir}.freq_fft{f_sel,d_sel}(i,:),results{subj,dir}.force_fft{f_sel,d_sel}(i,:),'k','LineWidth',1), hold on
%         end
%         xlabel("f (Hz)")
%         ylabel("|P1(f)|")
%         grid on       
%         xlim([0 20])
%         idx = idx+1;
%     end
% end
% 
% 
% 
% %Correlated Time Plot
% figure(),
% set(gcf,'color','w');
% sgtitle('Force Hold - Time Plots')
% idx = 1;
% for f_sel = 1:3
%     for d_sel = 1:3
%         subplot(3,3,idx),
%         sizeV = size(results{subj,dir}.force_hold{f_sel,d_sel});
%         for i = 1:sizeV(1)
%             plot(results{subj,dir}.time,results{subj,dir}.force_hold{f_sel,d_sel}(i,:),'k','LineWidth',2), hold on
%         end
%         xlabel('Time (s)')
%         ylabel('Force [N]')
%         grid on       
%         
%         idx = idx+1;
%     end
% end

%% 2-D Bar Plots Unperturbed 
clc, close all

subj = 20;
dir = 4;

ff = [15 20 25];
xx = [25 50 75];
% ff = {'F1' 'F2' 'F3'};
% xx = categorical({'900' '600' '300'});



%Stiffness 

figure(),
subplot(1,3,1)
set(gcf,'color','w');
b = bar(xx,results{subj,dir}.K_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results{subj,dir}.K_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results{subj,dir}.K_up_avg',results{subj,dir}.K_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 1000])
ylabel('Stiffness [N/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')



%Damping

subplot(1,3,2)
set(gcf,'color','w');
b = bar(xx,results{subj,dir}.B_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results{subj,dir}.B_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results{subj,dir}.B_up_avg',results{subj,dir}.B_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')


%Inertia 

subplot(1,3,3)
set(gcf,'color','w');
b = bar(xx,results{subj,dir}.M_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results{subj,dir}.M_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results{subj,dir}.M_up_avg',results{subj,dir}.M_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')

% Natural Frequency and Damping Ratio

figure(),
set(gcf,'color','w');

% Natural Frequency

subplot(1,2,1)
set(gcf,'color','w');
b = bar(xx,results{subj,dir}.up_omegan_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results{subj,dir}.up_omegan_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results{subj,dir}.up_omegan_avg',results{subj,dir}.up_omegan_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
% ylim([0 20])
ylabel('Natural Frequency [rad/s]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')

% Damping Ratio

subplot(1,2,2),
set(gcf,'color','w');
b = bar(xx,results{subj,dir}.up_dampr_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results{subj,dir}.up_dampr_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results{subj,dir}.up_dampr_avg',results{subj,dir}.up_dampr_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 1])
ylabel('Damping Ratio [-]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')
%% 2-D Bar Plots Unperturbed showing multiple subjects GROUPED (No Friction)
clc, close all

ff = [15 20 25];
xx = [25 50 75];

figure(),
for dir = 1:4
    subplot(2,2,dir)
    %Stiffness
    set(gcf,'color','w');
    errorbar(xx,imp_tot{dir}.K_tot_avg(1,:),imp_tot{dir}.K_tot_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.K_tot_avg(2,:),imp_tot{dir}.K_tot_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.K_tot_avg(3,:),imp_tot{dir}.K_tot_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

    grid on
    xlabel('Displacement [mm]')
    xlim([20 80])
    ylim([0 1000])
    ylabel('Stiffness [N/m]')
    legend('15 N','20 N','25 N')
    if dir == 1
        title('+x')
    elseif dir == 2
        title('+y')
    elseif dir == 3
        title('-x')
    elseif dir == 4
        title('-y')
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)

    %Damping
    set(gcf,'color','w');
    errorbar(xx,imp_tot{dir}.B_tot_avg(1,:),imp_tot{dir}.B_tot_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.B_tot_avg(2,:),imp_tot{dir}.B_tot_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.B_tot_avg(3,:),imp_tot{dir}.B_tot_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

    grid on
    xlabel('Displacement [mm]')
    xlim([20 80])
    ylim([0 40])
    ylabel('Damping [Ns/m]')
    legend('15 N','20 N','25 N')
    if dir == 1
        title('+x')
    elseif dir == 2
        title('+y')
    elseif dir == 3
        title('-x')
    elseif dir == 4
        title('-y')
    end

end

figure(),
for dir = 1:4
    subplot(2,2,dir)

    %Mass Percentage
    set(gcf,'color','w');
    errorbar(xx,imp_tot{dir}.M_tot_avg(1,:),imp_tot{dir}.M_tot_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.M_tot_avg(2,:),imp_tot{dir}.M_tot_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.M_tot_avg(3,:),imp_tot{dir}.M_tot_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

    grid on
    xlabel('Displacement [mm]')
    xlim([20 80])
    ylim([0 7])
    ylabel('Perc. Body Mass [\%]')
    legend('15 N','20 N','25 N')
    if dir == 1
        title('+x')
    elseif dir == 2
        title('+y')
    elseif dir == 3
        title('-x')
    elseif dir == 4
        title('-y')
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    %Natural Frequency
    set(gcf,'color','w');
    errorbar(xx,imp_tot{dir}.omega_tot_avg(1,:),imp_tot{dir}.omega_tot_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.omega_tot_avg(2,:),imp_tot{dir}.omega_tot_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.omega_tot_avg(3,:),imp_tot{dir}.omega_tot_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

    grid on
    xlabel('Displacement [mm]')
    xlim([20 80])
    ylim([0 35])
    ylabel('Natural Frequency [rad/s]')
    legend('15 N','20 N','25 N')
    if dir == 1
        title('+x')
    elseif dir == 2
        title('+y')
    elseif dir == 3
        title('-x')
    elseif dir == 4
        title('-y')
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    %Damping Ratio
    set(gcf,'color','w');
    errorbar(xx,imp_tot{dir}.damp_tot_avg(1,:),imp_tot{dir}.damp_tot_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.damp_tot_avg(2,:),imp_tot{dir}.damp_tot_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
    errorbar(xx,imp_tot{dir}.damp_tot_avg(3,:),imp_tot{dir}.damp_tot_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

    grid on
    xlabel('Displacement [mm]')
    xlim([20 80])
    ylim([0 1.5])
    ylabel('Damping Ratio [-]')
    legend('15 N','20 N','25 N')
    if dir == 1
        title('+x')
    elseif dir == 2
        title('+y')
    elseif dir == 3
        title('-x')
    elseif dir == 4
        title('-y')
    end
end


%% 2-D Bar Plots Unperturbed showing multiple subjects together (No Friction)
clc, close all

ff = [15 20 25];
xx = [25 50 75];

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Stiffness        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.K_up_avg(1,:),results{subj,dir}.K_up_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.K_up_avg(2,:),results{subj,dir}.K_up_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.K_up_avg(3,:),results{subj,dir}.K_up_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 1000])
        ylabel('Stiffness [N/m]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Damping        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.B_up_avg(1,:),results{subj,dir}.B_up_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.B_up_avg(2,:),results{subj,dir}.B_up_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.B_up_avg(3,:),results{subj,dir}.B_up_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 40])
        ylabel('Damping [Ns/m]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

% figure(),
% for dir = 1:4
%     subplot(2,2,dir)
%     for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18]%1:9
%         %Mass        
%         set(gcf,'color','w');
%         errorbar(xx,results{subj,dir}.M_up_avg(1,:),results{subj,dir}.M_up_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
%         errorbar(xx,results{subj,dir}.M_up_avg(2,:),results{subj,dir}.M_up_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
%         errorbar(xx,results{subj,dir}.M_up_avg(3,:),results{subj,dir}.M_up_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on
% 
%         grid on
%         xlabel('Displacement [mm]')
%         xlim([20 80])
%         ylim([0 5])
%         ylabel('Mass [kg]')
%         legend('15 N','20 N','25 N')
%         if dir == 1
%             title('+x')
%         elseif dir == 2
%             title('+y')
%         elseif dir == 3
%             title('-x')
%         elseif dir == 4
%             title('-y')
%         end
%     end
% end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Mass Percentage        
        set(gcf,'color','w');
        errorbar(xx,100*results{subj,dir}.M_up_avg(1,:)/results{subj,dir}.subj_mass,100*results{subj,dir}.M_up_std(1,:)/results{subj,dir}.subj_mass,'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,100*results{subj,dir}.M_up_avg(2,:)/results{subj,dir}.subj_mass,100*results{subj,dir}.M_up_std(2,:)/results{subj,dir}.subj_mass,':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,100*results{subj,dir}.M_up_avg(3,:)/results{subj,dir}.subj_mass,100*results{subj,dir}.M_up_std(3,:)/results{subj,dir}.subj_mass,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 7])
        ylabel('Perc. Body Mass [\%]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
%         subj
%         pause;
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Natural Frequency        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.up_omegan_avg(1,:),results{subj,dir}.up_omegan_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.up_omegan_avg(2,:),results{subj,dir}.up_omegan_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.up_omegan_avg(3,:),results{subj,dir}.up_omegan_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 35])
        ylabel('Natural Frequency [rad/s]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Damping Ratio        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.up_dampr_avg(1,:),results{subj,dir}.up_dampr_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.up_dampr_avg(2,:),results{subj,dir}.up_dampr_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.up_dampr_avg(3,:),results{subj,dir}.up_dampr_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 1.5])
        ylabel('Damping Ratio [-]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

%% 2-D Bar Plots Unperturbed showing multiple subjects together (Friction)
clc, close all

ff = [15 20 25];
xx = [25 50 75];

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Stiffness        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.K_fr_avg(1,:),results{subj,dir}.K_fr_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.K_fr_avg(2,:),results{subj,dir}.K_fr_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.K_fr_avg(3,:),results{subj,dir}.K_fr_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 1000])
        ylabel('Stiffness [N/m]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Damping        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.B_fr_avg(1,:),results{subj,dir}.B_fr_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.B_fr_avg(2,:),results{subj,dir}.B_fr_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.B_fr_avg(3,:),results{subj,dir}.B_fr_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 40])
        ylabel('Damping [Ns/m]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Mass        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.M_fr_avg(1,:),results{subj,dir}.M_fr_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.M_fr_avg(2,:),results{subj,dir}.M_fr_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.M_fr_avg(3,:),results{subj,dir}.M_fr_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 5])
        ylabel('Mass [kg]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Mass Percentage        
        set(gcf,'color','w');
        errorbar(xx,100*results{subj,dir}.M_fr_avg(1,:)/results{subj,dir}.subj_mass,100*results{subj,dir}.M_fr_std(1,:)/results{subj,dir}.subj_mass,'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,100*results{subj,dir}.M_fr_avg(2,:)/results{subj,dir}.subj_mass,100*results{subj,dir}.M_fr_std(2,:)/results{subj,dir}.subj_mass,':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,100*results{subj,dir}.M_fr_avg(3,:)/results{subj,dir}.subj_mass,100*results{subj,dir}.M_fr_std(3,:)/results{subj,dir}.subj_mass,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 7])
        ylabel('Perc. Body Mass [\%]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Natural Frequency        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.omegan_fr_avg(1,:),results{subj,dir}.omegan_fr_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.omegan_fr_avg(2,:),results{subj,dir}.omegan_fr_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.omegan_fr_avg(3,:),results{subj,dir}.omegan_fr_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 35])
        ylabel('Natural Frequency [rad/s]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end

figure(),
for dir = 1:4
    subplot(2,2,dir)
    for subj = [1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20]%1:9
        %Damping Ratio        
        set(gcf,'color','w');
        errorbar(xx,results{subj,dir}.dampr_fr_avg(1,:),results{subj,dir}.dampr_fr_std(1,:),'-','Color',[0 0.4470 0.7410],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.dampr_fr_avg(2,:),results{subj,dir}.dampr_fr_std(2,:),':','Color',[0.8500 0.3250 0.0980],'LineWidth',2), hold on
        errorbar(xx,results{subj,dir}.dampr_fr_avg(3,:),results{subj,dir}.dampr_fr_std(3,:),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2), hold on

        grid on
        xlabel('Displacement [mm]')
        xlim([20 80])
        ylim([0 1.5])
        ylabel('Damping Ratio [-]')
        legend('15 N','20 N','25 N')
        if dir == 1
            title('+x')
        elseif dir == 2
            title('+y')
        elseif dir == 3
            title('-x')
        elseif dir == 4
            title('-y')
        end
    end
end
%% 2-D Bar Plots Unperturbed - Friction vs Non-Friction
clc, close all

ff = [15 20 25];
xx = [25 50 75];
% ff = {'F1' 'F2' 'F3'};
% xx = categorical({'900' '600' '300'});

for subj = 20:20
    for dir = 1:4

        %Stiffness

        figure(),
        subplot(1,3,1)
        set(gcf,'color','w');
        matrix = [results{subj,dir}.K_up_avg',results{subj,dir}.K_fr_avg'];
        new_idx = [1 4 2 5 3 6];
        matrix_s = matrix(:,new_idx);
        b = bar(xx,matrix_s); hold on
        b(1).FaceColor = [0, 0.4470, 0.7410];
        b(2).FaceColor = [0, 0.4470, 0.7410];
        b(3).FaceColor = [0.8500, 0.3250, 0.0980];
        b(4).FaceColor = [0.8500, 0.3250, 0.0980];
        b(5).FaceColor = [0.9290, 0.6940, 0.1250];
        b(6).FaceColor = [0.9290, 0.6940, 0.1250];
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(results{subj,dir}.K_up_avg');
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        x(1,:) = b(1).XEndPoints;
        x(2,:) = b(3).XEndPoints;
        x(3,:) = b(5).XEndPoints;
        % Plot the errorbars
        errorbar(x',results{subj,dir}.K_up_avg',results{subj,dir}.K_up_std','k','linestyle','none');
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(results{subj,dir}.K_fr_avg');
        % Get the x coordinate of the bars
        xf = nan(nbars, ngroups);
        xf(1,:) = b(2).XEndPoints;
        xf(2,:) = b(4).XEndPoints;
        xf(3,:) = b(6).XEndPoints;
        % Plot the errorbars
        errorbar(xf',results{subj,dir}.K_fr_avg',results{subj,dir}.K_fr_std','k','linestyle','none');
        grid on
        xlabel('Displacement [mm]')
        % xlabel('Target Stiffness [N/m]')
        ylim([0 1000])
        ylabel('Stiffness [N/m]')
        legend('15 N','','20 N','','25 N')
        % legend('S','M','L')
        title('Ballistic Release')


        %Damping

        subplot(1,3,2)
        set(gcf,'color','w');
        matrix = [results{subj,dir}.B_up_avg',results{subj,dir}.B_fr_avg'];
        new_idx = [1 4 2 5 3 6];
        matrix_s = matrix(:,new_idx);
        b = bar(xx,matrix_s); hold on
        b(1).FaceColor = [0, 0.4470, 0.7410];
        b(2).FaceColor = [0, 0.4470, 0.7410];
        b(3).FaceColor = [0.8500, 0.3250, 0.0980];
        b(4).FaceColor = [0.8500, 0.3250, 0.0980];
        b(5).FaceColor = [0.9290, 0.6940, 0.1250];
        b(6).FaceColor = [0.9290, 0.6940, 0.1250];
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(results{subj,dir}.B_up_avg');
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        x(1,:) = b(1).XEndPoints;
        x(2,:) = b(3).XEndPoints;
        x(3,:) = b(5).XEndPoints;
        % Plot the errorbars
        errorbar(x',results{subj,dir}.B_up_avg',results{subj,dir}.B_up_std','k','linestyle','none');
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(results{subj,dir}.B_fr_avg');
        % Get the x coordinate of the bars
        xf = nan(nbars, ngroups);
        xf(1,:) = b(2).XEndPoints;
        xf(2,:) = b(4).XEndPoints;
        xf(3,:) = b(6).XEndPoints;
        % Plot the errorbars
        errorbar(xf',results{subj,dir}.B_fr_avg',results{subj,dir}.B_fr_std','k','linestyle','none');
        grid on
        xlabel('Displacement [mm]')
        % xlabel('Target Stiffness [N/m]')
        ylim([0 40])
        ylabel('Damping [Ns/m]')
        legend('15 N','','20 N','','25 N')
        % legend('S','M','L')
        title('Ballistic Release')

        %Inertia

        subplot(1,3,3)
        set(gcf,'color','w');
        matrix = [results{subj,dir}.M_up_avg',results{subj,dir}.M_fr_avg'];
        new_idx = [1 4 2 5 3 6];
        matrix_s = matrix(:,new_idx);
        b = bar(xx,matrix_s); hold on
        b(1).FaceColor = [0, 0.4470, 0.7410];
        b(2).FaceColor = [0, 0.4470, 0.7410];
        b(3).FaceColor = [0.8500, 0.3250, 0.0980];
        b(4).FaceColor = [0.8500, 0.3250, 0.0980];
        b(5).FaceColor = [0.9290, 0.6940, 0.1250];
        b(6).FaceColor = [0.9290, 0.6940, 0.1250];
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(results{subj,dir}.M_up_avg');
        % Get the x coordinate of the bars
        x = nan(nbars, ngroups);
        x(1,:) = b(1).XEndPoints;
        x(2,:) = b(3).XEndPoints;
        x(3,:) = b(5).XEndPoints;
        % Plot the errorbars
        errorbar(x',results{subj,dir}.M_up_avg',results{subj,dir}.M_up_std','k','linestyle','none');
        % Calculate the number of groups and number of bars in each group
        [ngroups,nbars] = size(results{subj,dir}.M_fr_avg');
        % Get the x coordinate of the bars
        xf = nan(nbars, ngroups);
        xf(1,:) = b(2).XEndPoints;
        xf(2,:) = b(4).XEndPoints;
        xf(3,:) = b(6).XEndPoints;
        % Plot the errorbars
        errorbar(xf',results{subj,dir}.M_fr_avg',results{subj,dir}.M_fr_std','k','linestyle','none');
        grid on
        xlabel('Displacement [mm]')
        % xlabel('Target Stiffness [N/m]')
        ylim([0 5])
        ylabel('Mass [kg]')
        legend('15 N','','20 N','','25 N')
        % legend('S','M','L')
        title('Ballistic Release')

    end
end

%% Friction Analysis - Check
clc, close all

for subj = 1:1
    for dir = 1:1
%         results{subj,dir}.Fc_fr_avg;
%         results{subj,dir}.Fc_fr_std;
% 
%         maximum(subj,dir) = max(max(results{subj,dir}.Fc_fr_avg));
%         [x,y]=find(results{subj,dir}.Fc_fr_avg==maximum(subj,dir));
% 
%         F_max(subj,dir) = x;
%         d_max(subj,dir) = y;
% 
%         minimum(subj,dir) = min(min(results{subj,dir}.Fc_fr_avg));
%         [x,y]=find(results{subj,dir}.Fc_fr_avg==minimum(subj,dir));
% 
%         F_min(subj,dir) = x;
%         d_min(subj,dir) = y;

          results{subj,dir}.e.K
          results{subj,dir}.e.B
          results{subj,dir}.e.M
    end
end

%% STATISTICAL ANALYSIS
clc, close all

subj = 20;
dir = 4;

%ANOVA-test on unperturbed
KKK=[];
BBB=[];
MMM=[];
ooo=[];
ddd=[];
ForceK=[];
ForceB=[];
ForceM=[];
Forceo=[];
Forced=[];
dispK=[];
dispB=[];
dispM=[];
dispo=[];
dispd=[];

for f_sel = 1:3
    for d_sel = 1:3
            KKK = [KKK;results{subj,dir}.K_up_tr{f_sel,d_sel,1}'];
            BBB = [BBB;results{subj,dir}.B_up_tr{f_sel,d_sel,1}'];
            MMM = [MMM;results{subj,dir}.M_up_tr{f_sel,d_sel,1}'];
            ooo = [ooo;results{subj,dir}.up_omegan_tr{f_sel,d_sel,1}'];
            ddd = [ddd;results{subj,dir}.up_dampr_tr{f_sel,d_sel,1}'];

            if f_sel == 1
                ForceK = [ForceK;15*ones(length(results{subj,dir}.K_up_tr{f_sel,d_sel,1}),1)];
                ForceB = [ForceB;15*ones(length(results{subj,dir}.B_up_tr{f_sel,d_sel,1}),1)];
                ForceM = [ForceM;15*ones(length(results{subj,dir}.M_up_tr{f_sel,d_sel,1}),1)];
                Forceo = [Forceo;15*ones(length(results{subj,dir}.up_omegan_tr{f_sel,d_sel,1}),1)];
                Forced = [Forced;15*ones(length(results{subj,dir}.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif f_sel == 2
                ForceK = [ForceK;20*ones(length(results{subj,dir}.K_up_tr{f_sel,d_sel,1}),1)];
                ForceB = [ForceB;20*ones(length(results{subj,dir}.B_up_tr{f_sel,d_sel,1}),1)];
                ForceM = [ForceM;20*ones(length(results{subj,dir}.M_up_tr{f_sel,d_sel,1}),1)];
                Forceo = [Forceo;20*ones(length(results{subj,dir}.up_omegan_tr{f_sel,d_sel,1}),1)];
                Forced = [Forced;20*ones(length(results{subj,dir}.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif f_sel == 3
                ForceK = [ForceK;25*ones(length(results{subj,dir}.K_up_tr{f_sel,d_sel,1}),1)];
                ForceB = [ForceB;25*ones(length(results{subj,dir}.B_up_tr{f_sel,d_sel,1}),1)];
                ForceM = [ForceM;25*ones(length(results{subj,dir}.M_up_tr{f_sel,d_sel,1}),1)];
                Forceo = [Forceo;25*ones(length(results{subj,dir}.up_omegan_tr{f_sel,d_sel,1}),1)];
                Forced = [Forced;25*ones(length(results{subj,dir}.up_dampr_tr{f_sel,d_sel,1}),1)];
            end

            if d_sel == 1
                dispK = [dispK;25*ones(length(results{subj,dir}.K_up_tr{f_sel,d_sel,1}),1)];
                dispB = [dispB;25*ones(length(results{subj,dir}.B_up_tr{f_sel,d_sel,1}),1)];
                dispM = [dispM;25*ones(length(results{subj,dir}.M_up_tr{f_sel,d_sel,1}),1)];
                dispo = [dispo;25*ones(length(results{subj,dir}.up_omegan_tr{f_sel,d_sel,1}),1)];
                dispd = [dispd;25*ones(length(results{subj,dir}.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif d_sel == 2
                dispK = [dispK;50*ones(length(results{subj,dir}.K_up_tr{f_sel,d_sel,1}),1)];
                dispB = [dispB;50*ones(length(results{subj,dir}.B_up_tr{f_sel,d_sel,1}),1)];
                dispM = [dispM;50*ones(length(results{subj,dir}.M_up_tr{f_sel,d_sel,1}),1)];
                dispo = [dispo;50*ones(length(results{subj,dir}.up_omegan_tr{f_sel,d_sel,1}),1)];
                dispd = [dispd;50*ones(length(results{subj,dir}.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif d_sel == 3
                dispK = [dispK;75*ones(length(results{subj,dir}.K_up_tr{f_sel,d_sel,1}),1)];
                dispB = [dispB;75*ones(length(results{subj,dir}.B_up_tr{f_sel,d_sel,1}),1)];
                dispM = [dispM;75*ones(length(results{subj,dir}.M_up_tr{f_sel,d_sel,1}),1)];
                dispo = [dispo;75*ones(length(results{subj,dir}.up_omegan_tr{f_sel,d_sel,1}),1)];
                dispd = [dispd;75*ones(length(results{subj,dir}.up_dampr_tr{f_sel,d_sel,1}),1)];
            end

    end
end
      
p_stiff = anovan(KKK,{ForceK dispK},'model','interaction','varnames',{'Force','Displacement'})
p_damp = anovan(BBB,{ForceB dispB},'model','interaction','varnames',{'Force','Displacement'})
p_mass = anovan(MMM,{ForceM dispM},'model','interaction','varnames',{'Force','Displacement'})
p_omega = anovan(ooo,{Forceo dispo},'model','interaction','varnames',{'Force','Displacement'})
p_damp = anovan(ddd,{Forced dispd},'model','interaction','varnames',{'Force','Displacement'})

%%
clearvars KKK MMM BBB ooo ddd ForceK Forced Forceo ForceM ForceB dispM dispB dispK dispo dispd 
clearvars p_stiff p_damp p_mass p_omega p_damp
clearvars results

%% KINEMATIC PLOTS - run this part to plot the results of the kinematic analysis
clc, close all
trial_l = 9;
subj = 2;
dir = 4;

%Marker Trajectories
figure(),
set(gcf,'color','w');
sgtitle('Trajectories of the Hand-Elbow-Shoulder Joints')
idx = 1;

for f_sel = 1:3
    for d_sel = 1:3
        for i = 1:trial_l
            subplot(3,3,idx),
            plot3(resultsK{subj,dir}.hand{f_sel,d_sel,i}.x,resultsK{subj,dir}.hand{f_sel,d_sel,i}.y,resultsK{subj,dir}.hand{f_sel,d_sel,i}.z,'.k'), hold on
            plot3(resultsK{subj,dir}.elbow{f_sel,d_sel,i}.x,resultsK{subj,dir}.elbow{f_sel,d_sel,i}.y,resultsK{subj,dir}.elbow{f_sel,d_sel,i}.z,'.r'), hold on
            plot3(resultsK{subj,dir}.shoulder{f_sel,d_sel,i}.x,resultsK{subj,dir}.shoulder{f_sel,d_sel,i}.y,resultsK{subj,dir}.shoulder{f_sel,d_sel,i}.z,'.b'), hold on, grid on
            xlabel('X [m]')
            ylabel('Y [m]')
            zlabel('Z [m]')
            xlim([-0.5 0])
            ylim([-0.9 -0.5])
            view(2)
            axis equal
        end
        idx = idx+1;
    end
end

%Mass Ellipsoid
en_f = 8; %Enalrgement factor of Arm Lengths
figure(),
set(gcf,'color','w');
sgtitle('End-Effector Mass Ellipsoid')
idx = 1;
for f_sel = 1:3
    for d_sel = 1:3
        for i = 1:trial_l
            
            subplot(3,3,idx),
            plot(resultsK{subj,dir}.Fxx{f_sel,d_sel,i,1},resultsK{subj,dir}.Fyy{f_sel,d_sel,i,1},'.k'), hold on
            plot(resultsK{subj,dir}.Fxx{f_sel,d_sel,i,end},resultsK{subj,dir}.Fyy{f_sel,d_sel,i,end},'.r'), hold on
            plot(en_f*[0,resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))],en_f*[0,-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1)),resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))], ...
                en_f*[-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1)),-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[0,resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))],en_f*[0,-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            plot(en_f*[resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end)),resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))], ...
                en_f*[-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end)),-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            grid on
            xlabel('Fx [N]')
            ylabel('Fy [N]')
            [t,s] = title(strcat('$M_x(t_0)=$ ',num2str(resultsK{subj,dir}.MEE_avg_t0{f_sel,d_sel}(1,1),'%.1f'),' kg , $M_x(t_f)=$ ',num2str(resultsK{subj,dir}.MEE_avg_tf{f_sel,d_sel}(1,1),'%.1f'),' kg'),...
                strcat('$M_y(t_0)=$ ',num2str(resultsK{subj,dir}.MEE_avg_t0{f_sel,d_sel}(2,2),'%.1f'),' kg , $M_y(t_f)=$ ',num2str(resultsK{subj,dir}.MEE_avg_tf{f_sel,d_sel}(2,2),'%.1f'),' kg'));
            t.FontSize = 9;
            t.FontWeight = 'normal';
            s.FontSize = 9;
            s.FontWeight = 'normal';
            axis equal
        end
        idx = idx+1;
    end
end

%Stiffness Ellipsoid
en_f = 800; %Enalrgement factor of Arm Lengths
figure(),
set(gcf,'color','w');
sgtitle('End-Effector Stiffness Ellipsoid')
idx = 1;
for f_sel = 1:3
    for d_sel = 1:3
        for i = 1:trial_l
            
            subplot(3,3,idx),
            plot(resultsK{subj,dir}.Fxx_k{f_sel,d_sel,i,1},resultsK{subj,dir}.Fyy_k{f_sel,d_sel,i,1},'.k'), hold on
            plot(resultsK{subj,dir}.Fxx_k{f_sel,d_sel,i,end},resultsK{subj,dir}.Fyy_k{f_sel,d_sel,i,end},'.r'), hold on
            plot(en_f*[0,resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))],en_f*[0,-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1)),resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(1).*cosd(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))], ...
                en_f*[-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1)),-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(1)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(1).*sind(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[0,resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))],en_f*[0,-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            plot(en_f*[resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end)),resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(end).*cosd(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))], ...
                en_f*[-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end)),-resultsK{subj,dir}.forearm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_e{f_sel,d_sel,i}(end)-resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))-resultsK{subj,dir}.arm{f_sel,d_sel,i}(end).*sind(resultsK{subj,dir}.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            grid on
            xlabel('Fx [N]')
            ylabel('Fy [N]')
            [t,s] = title(strcat('$K_x(t_0)=$ ',num2str(resultsK{subj,dir}.KEE_avg_t0{f_sel,d_sel}(1,1),'%.1f'),' N/m , $K_x(t_f)=$ ',num2str(resultsK{subj,dir}.KEE_avg_tf{f_sel,d_sel}(1,1),'%.1f'),' N/m'),...
                strcat('$K_y(t_0)=$ ',num2str(resultsK{subj,dir}.KEE_avg_t0{f_sel,d_sel}(2,2),'%.1f'),' N/m , $K_y(t_f)=$ ',num2str(resultsK{subj,dir}.KEE_avg_tf{f_sel,d_sel}(2,2),'%.1f'),' N/m'));
            t.FontSize = 9;
            t.FontWeight = 'normal';
            s.FontSize = 9;
            s.FontWeight = 'normal';
            axis equal
        end
        idx = idx+1;
    end
end

% %Angular Plots
% figure(),
% set(gcf,'color','w');
% sgtitle('Elbow-Shoulder Angles')
% idx = 1;
% for f_sel = 1:3
%     for d_sel = 1:3
%         for i = 1:trial_l
%             
%             subplot(3,3,idx),
%             plot(resultsK{subj,dir}.th_s{f_sel,d_sel,i},resultsK{subj,dir}.th_e{f_sel,d_sel,i}), hold on
%             grid on
%             xlabel('\theta_s [deg]')
%             ylabel('\theta_e [deg]')
% %             [t,s] = title(strcat('M_x(t_0)= ',num2str(resultsK{subj,dir}.MEE_avg_t0{f_sel,d_sel}(1,1),'%.1f'),' kg , M_x(t_f)= ',num2str(resultsK{subj,dir}.MEE_avg_tf{f_sel,d_sel}(1,1),'%.1f'),' kg'),...
% %                 strcat('M_y(t_0)= ',num2str(resultsK{subj,dir}.MEE_avg_t0{f_sel,d_sel}(2,2),'%.1f'),' kg , M_y(t_f)= ',num2str(resultsK{subj,dir}.MEE_avg_tf{f_sel,d_sel}(2,2),'%.1f'),' kg'));
%             t.FontSize = 9;
%             t.FontWeight = 'normal';
%             s.FontSize = 9;
%             s.FontWeight = 'normal';
%             axis equal
%         end
%         idx = idx+1;
%     end
% end

