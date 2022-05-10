% Change ln20, to choose plot things: force, displacement, velocity, etc.
clc, clear, close all

% load('ss3938_3949.mat', 'data');  % 6N perturbation on x, 
load('ss4129_4137.mat', 'data');  % 6N perturbation on x, 

Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 

colors = colormap('lines');
r = size(Data, 1); % subj
c = size(Data, 2); % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;
epoc_type = 2;
plot_type = 1;          % 1 displacement
                        % 2 force
                        % 3 force command
                        % 4 velocity
                        % 5 torque of 3rd joint
                        % 6 displacement (vector length in cartesian space)
                        % 7 force (vector length in cartesian space)
axh = zeros(f,r);



for ri = 1:r % subj
    for ci = 1:c
        for fi = 1:f % fce
            for di = 1:d % target distance
                axh(ri, fi) = subplot(d,f,d*(fi-1) + di);grid on;hold on;
                %for di = 3 % target distance
                for li = 1:p % perturbation
                    trial_num = length(Data(ri,ci,fi,di,:,li));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,ci,fi,di,ti,li}))
                            continue;
                        end

                        switch epoc_type
                            case 1
                                idx = find(Data{ri,ci,fi,di,ti,li}.Fp(2,:)~=0 & ...
                                    Data{ri,ci,fi,di,ti,li}.ts==4);  % pert at y
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if li == 1
                                    disp('ERROR: should use li == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
                                idx = (idx(1)-250):idx(end);
                                %idx = (idx(1)):idx(end);
                        end

                        time = t_step*(idx-idx(1));
                        idx_t{ri,ci,fi,di,ti,li} = idx;
                        time_t{ri,ci,fi,di,ti,li} = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,ci,fi,di,ti,li}.ox(1,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,ci,fi,di,ti,li}.f(1,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,ci,fi,di,ti,li}.Fp(1,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,ci,fi,di,ti,li}.v(1,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,ci,fi,di,ti,li}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,ci,fi,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,ci,fi,di,ti,li}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm force';

                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
            end
        end
    end
end
%xlim([0 0.7])
linkaxes(axh, 'xy');
%xlim([0 0.5]);
xlim([0 2])
sgtitle(titlestr);

%% Sys Ident - THIS SECTION OF CODE TAKES ABOUT 400secs to RUN!!!
clc, close all

%Select Force/Disp combination
f_sel = 1;
d_sel = 1;

tic

for f_sel= 1:3
    for d_sel = 1:3

        f_sel
        d_sel
        switch f_sel
            case 1
                Ftg = 15;
            case 2
                Ftg = 20;
            case 3
                Ftg = 25;
            otherwise
                disp('No Force Value')
        end

        switch d_sel
            case 1
                Dtg = 25;
            case 2
                Dtg = 50;
            case 3
                Dtg = 75;
            otherwise
                disp('No Displacement Value')
        end

        % Select Perturbation (1=unperturbed,2-7=perturbed with different starting
        % times)
        unpert = 1;
        pert = 4;

        trial_l = length(Data(1,1,f_sel,d_sel,:,pert));

        %Indices of Perturbation
        if pert > 1
            for i = 1:trial_l
                if size(Data{1,1,f_sel,d_sel,i,pert}) == [1 1]
                    idx_p(:,i) = find(Data{1,1,f_sel,d_sel,i,pert}.Fp(1,idx_t{1,1,f_sel,d_sel,i,pert}) < -0.25); %.Fp(2,.. for old .Fp(1,.. for new
                end
            end
        end

        %Create Average Profile of Unperturbed
        time_new = linspace(0,1.25,1.25e3);
        for i = 1:trial_l
            time_up{i,:} = time_t{1,1,f_sel,d_sel,i,1};
            force_interp(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,1}.f(1,idx_t{1,1,f_sel,d_sel,i,1}),time_new); %.f(2,.. for old .f(1,.. for new
            disp_interp(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,1}.ox(1,idx_t{1,1,f_sel,d_sel,i,1}),time_new); %.x(2,.. for old .x(1,.. for new

            %Identification on Single Trial
            force_interp_t(i,:) = -(force_interp(i,:)-force_interp(i,1));
            disp_interp_t(i,:) = disp_interp(i,:)-disp_interp(i,1);
            % SYS IDENT
            % Unperturbed
            Ts = (time_new(15)-time_new(14));
            data_est_UPs = iddata(disp_interp_t(i,:)',force_interp_t(i,:)',Ts);
            sysUP_s = tfest(data_est_UPs,2,0);
            [NUM_UPs,DEN_UPs] = tfdata(sysUP_s);
            K_est_up_s(i) = DEN_UPs{1}(3)/NUM_UPs{1}(3);
            B_est_up_s(i) = DEN_UPs{1}(2)/NUM_UPs{1}(3);
            M_est_up_s(i) = DEN_UPs{1}(1)/NUM_UPs{1}(3);
            FIT_up_s(i) = sysUP_s.Report.Fit.FitPercent;

            results.K_up_tr(f_sel,d_sel,i) = K_est_up_s(i);
            results.B_up_tr(f_sel,d_sel,i) = B_est_up_s(i);
            results.M_up_tr(f_sel,d_sel,i) = M_est_up_s(i);
            
        end
        
        %Average and STD of Impedance and FIT over trials
        K_est_up_avg = mean(K_est_up_s);
        K_est_up_std = std(K_est_up_s);
        B_est_up_avg = mean(B_est_up_s);
        B_est_up_std = std(B_est_up_s);
        M_est_up_avg = mean(M_est_up_s);
        M_est_up_std = std(M_est_up_s);
        FIT_est_up_avg = mean(FIT_up_s);
        FIT_est_up_std = std(FIT_up_s);

        force_up_avg = mean(force_interp);
        force_up_avg_t =-(force_up_avg-force_up_avg(1));
        disp_up_avg = mean(disp_interp);
        disp_up_avg_t = disp_up_avg-disp_up_avg(1);

        % SYS IDENT
        % Unperturbed
        Ts = (time_new(15)-time_new(14));
        data_est_UP = iddata(disp_up_avg_t',force_up_avg_t',Ts);
        sysUP = tfest(data_est_UP,2,0);
        [NUM_UP,DEN_UP] = tfdata(sysUP);
        K_est_up = DEN_UP{1}(3)/NUM_UP{1}(3);
        B_est_up = DEN_UP{1}(2)/NUM_UP{1}(3);
        M_est_up = DEN_UP{1}(1)/NUM_UP{1}(3);
        FIT_up = sysUP.Report.Fit.FitPercent;
        % bode(sysUP)
        
        results.FD_UP{f_sel,d_sel} = {time_new;force_interp_t;disp_interp_t};
        results.avg_FD_UP{f_sel,d_sel} = [time_new;force_up_avg_t;disp_up_avg_t];
        results.K_up(f_sel,d_sel) = K_est_up;
        results.B_up(f_sel,d_sel) = B_est_up;
        results.M_up(f_sel,d_sel)= M_est_up;
        results.TF.up{f_sel,d_sel} = sysUP;
        results.FIT_up(f_sel,d_sel) = FIT_up;
        
        %Computing Identified System Response to Ballistic Release
        [yy,tt,xx] = lsim(results.TF.up{f_sel,d_sel},results.avg_FD_UP{f_sel,d_sel}(2,:),results.avg_FD_UP{f_sel,d_sel}(1,:));
        results.tt_mod{f_sel,d_sel} = tt;
        results.disp_mod{f_sel,d_sel} = yy;
               
        results.K_up_avg(f_sel,d_sel) = K_est_up_avg;
        results.K_up_std(f_sel,d_sel) = K_est_up_std;
        results.B_up_avg(f_sel,d_sel) = B_est_up_avg;
        results.B_up_std(f_sel,d_sel) = B_est_up_std;
        results.M_up_avg(f_sel,d_sel) = M_est_up_avg;
        results.M_up_std(f_sel,d_sel) = M_est_up_std;
        results.FIT_up_avg(f_sel,d_sel) = FIT_est_up_avg;
        results.FIT_up_std(f_sel,d_sel) = FIT_est_up_std;
        
        time_new_long = linspace(0,2.5,2.5e3);
        for pert = 2:4%2:length(Data(1,1,1,1,1,:))
            %Create Average Profile of Perturbed
            for i = 1:trial_l
                if pert ~= 4
                    if size(Data{1,1,f_sel,d_sel,i,pert}) == [1 1]
                        time_up{i,:} = time_t{1,1,f_sel,d_sel,i,pert};
                        forcep_interp(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,pert}.f(1,idx_t{1,1,f_sel,d_sel,i,pert}),time_new);
                        dispp_interp(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,pert}.ox(1,idx_t{1,1,f_sel,d_sel,i,pert}),time_new);
                        force_command(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,pert}.Fp(1,idx_t{1,1,f_sel,d_sel,i,pert}),time_new);
                    end
                else
                    if size(Data{1,1,f_sel,d_sel,i,pert}) == [1 1]
                        time_up{i,:} = time_t{1,1,f_sel,d_sel,i,pert};
                        forcep_interp_ph(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,pert}.f(1,idx_t{1,1,f_sel,d_sel,i,pert}),time_new_long);
                        dispp_interp_ph(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,pert}.ox(1,idx_t{1,1,f_sel,d_sel,i,pert}),time_new_long);
                        force_command_ph(i,:) = interp1(time_up{i,:},Data{1,1,f_sel,d_sel,i,pert}.Fp(1,idx_t{1,1,f_sel,d_sel,i,pert}),time_new_long);
                        [max_fc,peak_fc] = (max(-force_command_ph(i,:)));
                        idx_fc(i) = peak_fc;

                        forcep_interp_al(i,:) = forcep_interp_ph(i,idx_fc(i)-150:idx_fc(i)+350);
                        dispp_interp_al(i,:) = dispp_interp_ph(i,idx_fc(i)-150:idx_fc(i)+350);
                        force_command_al(i,:) = force_command_ph(i,idx_fc(i)-150:idx_fc(i)+350);                        
                    end
                end
            end
            
            
            if pert ~= 4
                force_p_avg = mean(forcep_interp);
                disp_p_avg = mean(dispp_interp);
                force_c_avg = -mean(force_command);

                %Force and Displacement Difference
                disp_diff = disp_up_avg-disp_p_avg;
                disp_diff = disp_diff-disp_diff(1);
                force_diff = force_p_avg-force_up_avg;
                force_diff = force_diff - force_diff(1);
            else
                dispp_interp_al = dispp_interp_al-dispp_interp_al(:,1); %To align data
                force_p_avg = mean(forcep_interp_al);
                disp_p_avg = mean(dispp_interp_al);
                force_c_avg = -mean(force_command_al);

                %Force and Displacement Difference
                disp_diff = disp_p_avg(1)-disp_p_avg;
                force_diff = force_p_avg;
                force_diff = force_diff - force_diff(1);
            end

            % Indices of Perturbation Resampled
            idx_p_rs = find(force_c_avg == max(force_c_avg));
            
            % Windowing of Pulse Signals
            wind_v = 30:10:300; %samples after peak // 200 fairly good estimates
            FIT_p = 0;
            idxopt = 1;

            for wind = 30:10:300 %Searching for Best Window
                idxmin = idx_p_rs-150;
                idxmax = idx_p_rs+wind;
                if idxmax >= length(disp_diff)
                    idxmax = length(disp_diff);
                end
                if idx_p_rs <= 150
                    idxmin = 1;
                end

                disp_diff_ponly = disp_diff(idxmin:idxmax); %+200
                force_diff_ponly = force_diff(idxmin:idxmax);
                force_c_ponly = force_c_avg(idxmin:idxmax);

                %SYS ID
                % Perturbed
                Ts = (time_new(15)-time_new(14));
                % data_est_P = iddata(disp_diff',force_diff',Ts);
                data_est_P{idxopt} = iddata(disp_diff_ponly',force_diff_ponly',Ts);
                sysP{idxopt} = tfest(data_est_P{idxopt},2,0);
                [NUM_P,DEN_P] = tfdata(sysP{idxopt});
                K_est_p(idxopt) = DEN_P{1}(3)/NUM_P{1}(3);
                B_est_p(idxopt) = DEN_P{1}(2)/NUM_P{1}(3);
                M_est_p(idxopt) = DEN_P{1}(1)/NUM_P{1}(3);
                FIT_p(idxopt) = sysP{idxopt}.Report.Fit.FitPercent;
                % bode(sysUP)

                %Getting all Windows
                if pert ~= 4
                    results.FD_P{pert-1}.p{f_sel,d_sel} = {forcep_interp;dispp_interp};
                else
                    results.FD_P{pert-1}.p{f_sel,d_sel} = {forcep_interp_al;dispp_interp_al};
                end
                results.avg_FD_P{pert-1}.p{f_sel,d_sel} = [force_p_avg;disp_p_avg];
                results.K_p{pert-1}.p(f_sel,d_sel,idxopt) = K_est_p(idxopt);
                results.B_p{pert-1}.p(f_sel,d_sel,idxopt) = B_est_p(idxopt);
                results.M_p{pert-1}.p(f_sel,d_sel,idxopt)= M_est_p(idxopt);
                results.TF.p{f_sel,d_sel,pert-1,idxopt} = sysP{idxopt};
                results.FIT_p{pert-1}.p(f_sel,d_sel,idxopt) = FIT_p(idxopt);

                [yy,tt,xx] = lsim(results.TF.p{f_sel,d_sel,pert-1,idxopt},results.avg_FD_UP{f_sel,d_sel}(2,:),results.avg_FD_UP{f_sel,d_sel}(1,:));
                results.tt_modP{pert-1}.p{f_sel,d_sel,idxopt} = tt;
                results.disp_modP{pert-1}.p{f_sel,d_sel,idxopt} = yy;
                
                idxopt = idxopt+1;
            end

%             [FIT_p_opt,opt_fit] = max(FIT_p);

            % For ss3938_3949
%             disp_diff_ponly = disp_diff(1:idx_p_rs+wind_v(opt_fit));
%             force_diff_ponly = force_diff(1:idx_p_rs+wind_v(opt_fit));
%             force_c_ponly = force_c_avg(1:idx_p_rs+wind_v(opt_fit));
            
            % For ss4129_4137
%             idxmax_opt = idx_p_rs+wind_v(opt_fit);
%             if idxmax_opt >= length(disp_diff)
%                 idxmax_opt = length(disp_diff);
%             end
%             disp_diff_ponly = disp_diff(idxmin:idxmax_opt);
%             force_diff_ponly = force_diff(idxmin:idxmax_opt);
%             force_c_ponly = force_c_avg(idxmin:idxmax_opt);

%             results.avg_FD_P{pert-1}.p{f_sel,d_sel} = [force_diff_ponly;disp_diff_ponly];
%             results.K_p{pert-1}.p(f_sel,d_sel) = K_est_p(opt_fit);
%             results.B_p{pert-1}.p(f_sel,d_sel) = B_est_p(opt_fit);
%             results.M_p{pert-1}.p(f_sel,d_sel)= M_est_p(opt_fit);
%             results.TF.p{f_sel,d_sel,pert-1} = sysP{opt_fit};
%             results.FIT_p{pert-1}.p(f_sel,d_sel) = FIT_p_opt;


        end
    end
end

toc
% FIGURES
% figure(),
%     set(gcf,'color','w');
%     subplot(3,2,1)
%     plot(time_new,force_diff), hold on
%     plot(time_new,force_c_avg), grid on
%     ylabel('Force Diff. [N]')
%     subplot(3,2,3)
%     plot(time_new,disp_diff), grid on
%     ylabel('Displacemente Diff. [m]')
%     xlabel('Time [s]')
%     subplot(3,2,5)
%     plot(disp_diff,force_diff), grid on
%     xlabel('Displacemente Diff. [m]')
%     ylabel('Force Diff. [N]')
%     
%     subplot(3,2,2)
%     plot(time_new(1:idx_p_rs+200),force_diff(1:idx_p_rs+200)), hold on
%     plot(time_new(1:idx_p_rs+200),force_c_avg(1:idx_p_rs+200)), grid on
%     ylabel('Force Diff. [N]')
%     subplot(3,2,4)
%     plot(time_new(1:idx_p_rs+200),disp_diff(1:idx_p_rs+200)), grid on
%     ylabel('Displacemente Diff. [m]')
%     xlabel('Time [s]')
%     subplot(3,2,6)
%     plot(disp_diff(1:idx_p_rs+200),force_diff(1:idx_p_rs+200)), grid on
%     xlabel('Displacemente Diff. [m]')
%     ylabel('Force Diff. [N]')
% 
% figure(),
%     set(gcf,'color','w');
%     % Plot Release Aligned - Unperturbed
%     subplot(4,2,1),
%     for i = 1:trial_l
%         plot(time_t{1,1,f_sel,d_sel,i,unpert},Data{1,1,f_sel,d_sel,i,unpert}.f(2,idx_t{1,1,f_sel,d_sel,i,unpert})), hold on
%     end
%     plot(time_new,force_up_avg,'r','LineWidth',2),
%     ylabel('Force [N]')
%     grid on
%     
%     subplot(4,2,3),
%     for i = 1:trial_l
%         plot(time_t{1,1,f_sel,d_sel,i,unpert},Data{1,1,f_sel,d_sel,i,unpert}.x(2,idx_t{1,1,f_sel,d_sel,i,unpert})), hold on
%     end
%     plot(time_new,disp_up_avg,'r','LineWidth',2),
%     ylabel('Displacement [m]')
%     xlabel('Time [s]')
%     grid on
%     
%     subplot(4,2,5),
%     plot(time_new,force_up_avg_t), hold on
%     xlabel('Time [s]')
%     ylabel('Force [N]')
%     grid on
%     
%     subplot(4,2,7),
%     plot(time_new,disp_up_avg_t), hold on
%     xlabel('Time [s]')
%     ylabel('Displacement [m]')
%     grid on
%     
%     
%     % Plot Release Aligned - Perturbed (1)
%     subplot(4,2,2),
%     for i = 1:trial_l
%         if size(Data{1,1,f_sel,d_sel,i,pert}) == [1 1]
%             plot(time_t{1,1,f_sel,d_sel,i,pert},Data{1,1,f_sel,d_sel,i,pert}.f(2,idx_t{1,1,f_sel,d_sel,i,pert})), hold on
%         end
%     end
%     plot(time_new,force_p_avg,'r','LineWidth',2),
%     ylabel('Force [N]')
%     grid on
%     
%     subplot(4,2,4),
%     for i = 1:trial_l
%         if size(Data{1,1,f_sel,d_sel,i,pert}) == [1 1]
%             plot(time_t{1,1,f_sel,d_sel,i,pert},Data{1,1,f_sel,d_sel,i,pert}.x(2,idx_t{1,1,f_sel,d_sel,i,pert})), hold on
%         end
%     end
%     plot(time_new,disp_p_avg,'r','LineWidth',2),
%     ylabel('Displacement [m]')
%     xlabel('Time [s]')
%     grid on
% 
% 
% figure(),
%     options = bodeoptions;
%     options.FreqUnits = 'Hz';
%     set(gcf,'color','w');
%     bode(sysUP,options), hold on
%     bode(sysP,options), grid on
%     legend('Unperturbed','Perturbed')


% resume = table(Ftg,Dtg,K_est_up,K_est_p,B_est_up,B_est_p,M_est_up,M_est_p);
% resume.Properties.VariableNames = {'Force [N]','Distance [mm]','Stiff. UP [N/m]','Stiff. P [N/m]','Damp. UP [Ns/m]','Damp. P [Ns/m]','Mass UP [kg]','Mass P [kg]'};
% 
% resume

%% TO RUN ONLY IF YOU WANT TO CANCEL BAD FITS

for f_sel= 1:3
    for d_sel = 1:3

        switch f_sel
            case 1
                Ftg = 15;
            case 2
                Ftg = 20;
            case 3
                Ftg = 25;
            otherwise
                disp('No Force Value')
        end

        switch d_sel
            case 1
                Dtg = 25;
            case 2
                Dtg = 50;
            case 3
                Dtg = 75;
            otherwise
                disp('No Displacement Value')
        end

        for pert = 2:3
            if results.FIT_p{pert-1}.p(f_sel,d_sel) < 85
                results.K_p{pert-1}.p(f_sel,d_sel) = NaN;
                results.B_p{pert-1}.p(f_sel,d_sel) = NaN;
                results.M_p{pert-1}.p(f_sel,d_sel)= NaN;
                results.FIT_p{pert-1}.p(f_sel,d_sel) = NaN;
            end
        end
    end
end

%% STATISTICAL ANALYSIS
clc, close all

%T-test
for i = 1:9
    if i == 1
        f_sel = 1;
        d_sel = 1;
    elseif i == 2
        f_sel = 1;
        d_sel = 2;
    elseif i == 3
        f_sel = 1;
        d_sel = 3;
    elseif i == 4
        f_sel = 2;
        d_sel = 1;
    elseif i == 5
        f_sel = 2;
        d_sel = 2;
    elseif i == 6
        f_sel = 2;
        d_sel = 3;
    elseif i == 7
        f_sel = 3;
        d_sel = 1;
    elseif i == 8
        f_sel = 3;
        d_sel = 2;
    elseif i == 9
        f_sel = 3;
        d_sel = 3;
    end

    for j = 1:9

        if j == 1
            f_selt = 1;
            d_selt = 1;
        elseif j == 2
            f_selt = 1;
            d_selt = 2;
        elseif j == 3
            f_selt = 1;
            d_selt = 3;
        elseif j == 4
            f_selt = 2;
            d_selt = 1;
        elseif j == 5
            f_selt = 2;
            d_selt = 2;
        elseif j == 6
            f_selt = 2;
            d_selt = 3;
        elseif j == 7
            f_selt = 3;
            d_selt = 1;
        elseif j == 8
            f_selt = 3;
            d_selt = 2;
        elseif j == 9
            f_selt = 3;
            d_selt = 3;
        end

        Kxx = results.K_up_tr(f_sel,d_sel,:);
        Kyy = results.K_up_tr(f_selt,d_selt,:);
        Kxx = Kxx(:);
        Kyy = Kyy(:);
        [hK,pK] = ttest2(Kxx,Kyy);

        Bxx = results.B_up_tr(f_sel,d_sel,:);
        Byy = results.B_up_tr(f_selt,d_selt,:);
        Bxx = Bxx(:);
        Byy = Byy(:);
        [hB,pB] = ttest2(Bxx,Byy);

        Mxx = results.M_up_tr(f_sel,d_sel,:);
        Myy = results.M_up_tr(f_selt,d_selt,:);
        Mxx = Mxx(:);
        Myy = Myy(:);
        [hM,pM] = ttest2(Mxx,Myy);

        p_val_K(i,j) = pK;
        p_val_B(i,j) = pB;
        p_val_M(i,j) = pM;

        if p_val_K(i,j) > 0.05
            ver_K(i,j) = 1;
        else
            ver_K(i,j) = 0;
        end

        if p_val_B(i,j) > 0.05
            ver_B(i,j) = 1;
        else
            ver_B(i,j)= 0;
        end

        if p_val_M(i,j) > 0.05
            ver_M(i,j) = 1;
        else
            ver_M(i,j) = 0;
        end
    end

end

%ANOVA-test
idx = 1;
for i = 1:3
    for j = 1:3
        for k = 1:length(results.K_up_tr(i,j,:))
            
            KKK(idx,1) = results.K_up_tr(i,j,k);
            BBB(idx,1) = results.B_up_tr(i,j,k);
            MMM(idx,1) = results.M_up_tr(i,j,k);
            
        
            if i == 1
                Force(idx,1) = 15;
            elseif i == 2
                Force(idx,1) = 20;
            elseif i == 3
                Force(idx,1) = 25;
            end

            if j == 1
                Disp(idx,1) = 25;
            elseif j == 2
                Disp(idx,1) = 50;
            elseif j == 3
                Disp(idx,1) = 75;
            end


            idx = idx+1;
        end
    end
end
      
p_stiff = anovan(KKK,{Force Disp},'model','interaction','varnames',{'Force','Displacement'})
p_damp = anovan(BBB,{Force Disp},'model','interaction','varnames',{'Force','Displacement'})
p_mass = anovan(MMM,{Force Disp},'model','interaction','varnames',{'Force','Displacement'})


% for i = 1:3
%     for j = 1:3
%         if i == 1
%             idx = 1;
%         elseif i == 2
%             idx = 16;
%         elseif i == 3
%             idx = 31;
%         end
%         for k = 1:15
%             KK2(idx,j) = results.K_up_tr(i,j,k);
%             idx = idx+1;
%         end
%     end
% end
% 
% anova2(KK2,15)
%%
clc
close all,

tpause = 0.2; %Speed of plotting
wind_start = 20; %From which time window to start showing

ff = [15 20 25];
xx = [25 50 75];

[XX,FF] = meshgrid(xx,ff);

% Interpolate Impedance Matrix to get Surface

fq = 15:1:25;
xq = 25:1:75;

[Xq,Fq] = meshgrid(xq,fq);

Kup_interp = interp2(XX,FF,results.K_up,Xq,Fq);
Bup_interp = interp2(XX,FF,results.B_up,Xq,Fq);
Mup_interp = interp2(XX,FF,results.M_up,Xq,Fq);
FITup_interp = interp2(XX,FF,results.FIT_up,Xq,Fq);



%--------------------------------------------

% % Initialize video
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
% myvideo.Quality = 100;
% open(myVideo)


figure(),%'Units','normalized','Position',[0 0 1 1])
set(gcf,'color','w');
plot3(FF,XX,results.K_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.K_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.K_up_avg+results.K_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.K_up_avg-results.K_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
s = surf(Fq,Xq,Kup_interp); hold on
s.EdgeColor = 'none';
colorbar
colormap summer
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Stiffness [N/m]')
zlim([0 1000])
for widx = wind_start:length(wind_v)
    for p = 1:3
        pplot = plot3(FF,XX,results.K_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 2
            pplot(1).Marker = "o";
            pplot(2).Marker = "o";
            pplot(3).Marker = "o";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 3
            pplot(1).Marker = "+";
            pplot(2).Marker = "+";
            pplot(3).Marker = "+";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 4
            pplot(1).Marker = "*";
            pplot(2).Marker = "*";
            pplot(3).Marker = "*";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 5
            pplot(1).Marker = "s";
            pplot(2).Marker = "s";
            pplot(3).Marker = "s";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 6
            pplot(1).Marker = "d";
            pplot(2).Marker = "d";
            pplot(3).Marker = "d";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 1000])
    pause(tpause)
% 
%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
end
% close(myVideo)

% 
% % Initialize video
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
% myvideo.Quality = 100;
% open(myVideo)

figure(),
set(gcf,'color','w');
plot3(FF,XX,results.B_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.B_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.B_up_avg+results.B_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.B_up_avg-results.B_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
s = surf(Fq,Xq,Bup_interp); grid on
s.EdgeColor = 'none';
colorbar
colormap cool
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Damping [Ns/m]')
zlim([0 50]) 
for widx = wind_start:length(wind_v)
    for p = 1:3
        pplot = plot3(FF,XX,results.B_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 2
            pplot(1).Marker = "o";
            pplot(2).Marker = "o";
            pplot(3).Marker = "o";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 3
            pplot(1).Marker = "+";
            pplot(2).Marker = "+";
            pplot(3).Marker = "+";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 4
            pplot(1).Marker = "*";
            pplot(2).Marker = "*";
            pplot(3).Marker = "*";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 5
            pplot(1).Marker = "s";
            pplot(2).Marker = "s";
            pplot(3).Marker = "s";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 6
            pplot(1).Marker = "d";
            pplot(2).Marker = "d";
            pplot(3).Marker = "d";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 50])
    pause(tpause)

%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
end
% close(myVideo)

% % Initialize video
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
% myvideo.Quality = 100;
% open(myVideo)

figure(),
set(gcf,'color','w');
plot3(FF,XX,results.M_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.M_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.M_up_avg+results.M_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.M_up_avg-results.M_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
s = surf(Fq,Xq,Mup_interp); grid on
s.EdgeColor = 'none';
colorbar
colormap hot
caxis([1 4])
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Mass [kg]')
zlim([0 4])
for widx = wind_start:length(wind_v)
    for p = 1:3
        pplot = plot3(FF,XX,results.M_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 2
            pplot(1).Marker = "o";
            pplot(2).Marker = "o";
            pplot(3).Marker = "o";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 3
            pplot(1).Marker = "+";
            pplot(2).Marker = "+";
            pplot(3).Marker = "+";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 4
            pplot(1).Marker = "*";
            pplot(2).Marker = "*";
            pplot(3).Marker = "*";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 5
            pplot(1).Marker = "s";
            pplot(2).Marker = "s";
            pplot(3).Marker = "s";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 6
            pplot(1).Marker = "d";
            pplot(2).Marker = "d";
            pplot(3).Marker = "d";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 4])
    pause(tpause)

%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
end
% close(myVideo)

% % Initialize video
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
% myvideo.Quality = 100;
% open(myVideo)

figure()
set(gcf,'color','w');
plot3(FF,XX,results.FIT_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.FIT_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.FIT_up_avg+results.FIT_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.FIT_up_avg-results.FIT_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
s = surf(Fq,Xq,FITup_interp); grid on
s.EdgeColor = 'none';
colorbar
colormap spring
caxis([90 100])
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Model FIT [%]')
zlim([0 100])
for widx = wind_start:length(wind_v)
    for p = 1:3
        pplot = plot3(FF,XX,results.FIT_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 2
            pplot(1).Marker = "o";
            pplot(2).Marker = "o";
            pplot(3).Marker = "o";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 3
            pplot(1).Marker = "+";
            pplot(2).Marker = "+";
            pplot(3).Marker = "+";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 4
            pplot(1).Marker = "*";
            pplot(2).Marker = "*";
            pplot(3).Marker = "*";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 5
            pplot(1).Marker = "s";
            pplot(2).Marker = "s";
            pplot(3).Marker = "s";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        elseif p == 6
            pplot(1).Marker = "d";
            pplot(2).Marker = "d";
            pplot(3).Marker = "d";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 100])
    pause(tpause)
%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);this.
end
% close(myVideo)

%%

% Unperturbed
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_UP{1, 1}{1, 1},results.FD_UP{1, 1}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,1}(1,:),results.avg_FD_UP{1,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_UP{1, 2}{1, 1},results.FD_UP{1, 2}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,2}(1,:),results.avg_FD_UP{1,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_UP{1, 3}{1, 1},results.FD_UP{1, 3}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,3}(1,:),results.avg_FD_UP{1,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_UP{2, 1}{1, 1},results.FD_UP{2, 1}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,1}(1,:),results.avg_FD_UP{2,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_UP{2, 2}{1, 1},results.FD_UP{2, 2}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,2}(1,:),results.avg_FD_UP{2,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_UP{2, 3}{1, 1},results.FD_UP{2, 3}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,3}(1,:),results.avg_FD_UP{2,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_UP{3, 1}{1, 1},results.FD_UP{3, 1}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,1}(1,:),results.avg_FD_UP{3,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_UP{3, 2}{1, 1},results.FD_UP{3, 2}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,2}(1,:),results.avg_FD_UP{3,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_UP{3, 3}{1, 1},results.FD_UP{3, 3}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,3}(1,:),results.avg_FD_UP{3,3}(2,:),'r','LineWidth',2), grid on
sgtitle('Force [N] in Unperturbed Ballistic Release')

figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_UP{1, 1}{1, 1},results.FD_UP{1, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,1}(1,:),results.avg_FD_UP{1,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 1},results.disp_mod{1, 1},'b--','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_UP{1, 2}{1, 1},results.FD_UP{1, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,2}(1,:),results.avg_FD_UP{1,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 2},results.disp_mod{1, 2},'b--','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_UP{1, 3}{1, 1},results.FD_UP{1, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,3}(1,:),results.avg_FD_UP{1,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 3},results.disp_mod{1, 3},'b--','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_UP{2, 1}{1, 1},results.FD_UP{2, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,1}(1,:),results.avg_FD_UP{2,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 1},results.disp_mod{2, 1},'b--','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_UP{2, 2}{1, 1},results.FD_UP{2, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,2}(1,:),results.avg_FD_UP{2,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 2},results.disp_mod{2, 2},'b--','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_UP{2, 3}{1, 1},results.FD_UP{2, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,3}(1,:),results.avg_FD_UP{2,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 3},results.disp_mod{2, 3},'b--','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_UP{3, 1}{1, 1},results.FD_UP{3, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,1}(1,:),results.avg_FD_UP{3,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 1},results.disp_mod{3, 1},'b--','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_UP{3, 2}{1, 1},results.FD_UP{3, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,2}(1,:),results.avg_FD_UP{3,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 2},results.disp_mod{3, 2},'b--','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_UP{3, 3}{1, 1},results.FD_UP{3, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,3}(1,:),results.avg_FD_UP{3,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 3},results.disp_mod{3, 3},'b--','LineWidth',2), grid on
sgtitle('Displacement [mm] in Unperturbed Ballistic Release')

%Ballistic Response of Pulse Identified Models
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_UP{1, 1}{1, 1},results.FD_UP{1, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,1}(1,:),results.avg_FD_UP{1,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 1},results.disp_mod{1, 1},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{1,1,w},results.disp_modP{1}.p{1,1,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{1,1,w},results.disp_modP{2}.p{1,1,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{1,1,w},results.disp_modP{3}.p{1,1,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,2), plot(results.FD_UP{1, 2}{1, 1},results.FD_UP{1, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,2}(1,:),results.avg_FD_UP{1,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 2},results.disp_mod{1, 2},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{1,2,w},results.disp_modP{1}.p{1,2,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{1,2,w},results.disp_modP{2}.p{1,2,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{1,2,w},results.disp_modP{3}.p{1,2,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,3), plot(results.FD_UP{1, 3}{1, 1},results.FD_UP{1, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,3}(1,:),results.avg_FD_UP{1,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 3},results.disp_mod{1, 3},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{1,3,w},results.disp_modP{1}.p{1,3,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{1,3,w},results.disp_modP{2}.p{1,3,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{1,3,w},results.disp_modP{3}.p{1,3,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,4), plot(results.FD_UP{2, 1}{1, 1},results.FD_UP{2, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,1}(1,:),results.avg_FD_UP{2,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 1},results.disp_mod{2, 1},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{2,1,w},results.disp_modP{1}.p{2,1,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{2,1,w},results.disp_modP{2}.p{2,1,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{2,1,w},results.disp_modP{3}.p{2,1,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,5), plot(results.FD_UP{2, 2}{1, 1},results.FD_UP{2, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,2}(1,:),results.avg_FD_UP{2,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 2},results.disp_mod{2, 2},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{2,2,w},results.disp_modP{1}.p{2,2,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{2,2,w},results.disp_modP{2}.p{2,2,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{2,2,w},results.disp_modP{3}.p{2,2,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,6), plot(results.FD_UP{2, 3}{1, 1},results.FD_UP{2, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,3}(1,:),results.avg_FD_UP{2,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 3},results.disp_mod{2, 3},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{2,3,w},results.disp_modP{1}.p{2,3,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{2,3,w},results.disp_modP{2}.p{2,3,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{2,3,w},results.disp_modP{3}.p{2,3,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,7), plot(results.FD_UP{3, 1}{1, 1},results.FD_UP{3, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,1}(1,:),results.avg_FD_UP{3,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 1},results.disp_mod{3, 1},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{3,1,w},results.disp_modP{1}.p{3,1,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{3,1,w},results.disp_modP{2}.p{3,1,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{3,1,w},results.disp_modP{3}.p{3,1,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,8), plot(results.FD_UP{3, 2}{1, 1},results.FD_UP{3, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,2}(1,:),results.avg_FD_UP{3,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 2},results.disp_mod{3, 2},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{3,2,w},results.disp_modP{1}.p{3,2,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{3,2,w},results.disp_modP{2}.p{3,2,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{3,2,w},results.disp_modP{3}.p{3,2,w},'c--','LineWidth',1), hold on
                end
                grid on
subplot(3,3,9), plot(results.FD_UP{3, 3}{1, 1},results.FD_UP{3, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,3}(1,:),results.avg_FD_UP{3,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 3},results.disp_mod{3, 3},'b--','LineWidth',2), hold on
                for w = wind_start:length(wind_v)
                    plot(results.tt_modP{1}.p{3,3,w},results.disp_modP{1}.p{3,3,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{3,3,w},results.disp_modP{2}.p{3,3,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{3,3,w},results.disp_modP{3}.p{3,3,w},'c--','LineWidth',1), hold on
                end
                grid on
sgtitle('Displacement [mm] in Unperturbed Ballistic Release')



% Perturbed on Force Hold
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 1}.p{1, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_P{1, 1}.p{1, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_P{1, 1}.p{1, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_P{1, 1}.p{2, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_P{1, 1}.p{2, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_P{1, 1}.p{2, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_P{1, 1}.p{3, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_P{1, 1}.p{3, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_P{1, 1}.p{3, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,3}(1,:),'r','LineWidth',2), grid on
sgtitle('Force [N] in Pulse @ Force Hold')

figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 1}.p{1, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_P{1, 1}.p{1, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_P{1, 1}.p{1, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_P{1, 1}.p{2, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_P{1, 1}.p{2, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_P{1, 1}.p{2, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_P{1, 1}.p{3, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_P{1, 1}.p{3, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_P{1, 1}.p{3, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,3}(2,:),'r','LineWidth',2), grid on
sgtitle('Displacement [mm] in Pulse @ Force Hold')

% Perturbed on Release
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 2}.p{1, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_P{1, 2}.p{1, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_P{1, 2}.p{1, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_P{1, 2}.p{2, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_P{1, 2}.p{2, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_P{1, 2}.p{2, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_P{1, 2}.p{3, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_P{1, 2}.p{3, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_P{1, 2}.p{3, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,3}(1,:),'r','LineWidth',2), grid on
sgtitle('Force [N] in Pulse @ Relsease')

figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 2}.p{1, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_P{1, 2}.p{1, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_P{1, 2}.p{1, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_P{1, 2}.p{2, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_P{1, 2}.p{2, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_P{1, 2}.p{2, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_P{1, 2}.p{3, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_P{1, 2}.p{3, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_P{1, 2}.p{3, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,3}(2,:),'r','LineWidth',2), grid on
sgtitle('Displacement [mm] in Pulse @ Release')

% Perturbed on Pos Hold
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 3}.p{1, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_P{1, 3}.p{1, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_P{1, 3}.p{1, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_P{1, 3}.p{2, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_P{1, 3}.p{2, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_P{1, 3}.p{2, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_P{1, 3}.p{3, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_P{1, 3}.p{3, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_P{1, 3}.p{3, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,3}(1,:),'r','LineWidth',2), grid on
sgtitle('Force [N] in Pulse @ Position Hold')

figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 3}.p{1, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_P{1, 3}.p{1, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_P{1, 3}.p{1, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_P{1, 3}.p{2, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_P{1, 3}.p{2, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_P{1, 3}.p{2, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,3}(2,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_P{1, 3}.p{3, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,1}(2,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_P{1, 3}.p{3, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,2}(2,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_P{1, 3}.p{3, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,3}(2,:),'r','LineWidth',2), grid on
sgtitle('Displacement [mm] in Pulse @ Position Hold')

%Bode Plots
options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
figure(),
set(gcf,'color','w');
subplot(3,3,1), bodemag(results.TF.up{1, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{1,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{1,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{1,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,2), bodemag(results.TF.up{1, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{1,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{1,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{1,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,3), bodemag(results.TF.up{1, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{1,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{1,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{1,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,4), bodemag(results.TF.up{2, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{2,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{2,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{2,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,5), bodemag(results.TF.up{2, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{2,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{2,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{2,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,6), bodemag(results.TF.up{2, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{2,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{2,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{2,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,7), bodemag(results.TF.up{3, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{3,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{3,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{3,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,8), bodemag(results.TF.up{3, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{3,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{3,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{3,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,9), bodemag(results.TF.up{3, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodemag(results.TF.p{3,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{3,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{3,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
sgtitle('Bode Plot - Magnitude - Identified Models')


options.MagVisible = 'off';
figure(),
set(gcf,'color','w');
subplot(3,3,1), bodeplot(results.TF.up{1, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{1,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{1,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{1,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,2), bodeplot(results.TF.up{1, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{1,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{1,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{1,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,3), bodeplot(results.TF.up{1, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{1,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{1,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{1,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,4), bodeplot(results.TF.up{2, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{2,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{2,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{2,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,5), bodeplot(results.TF.up{2, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{2,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{2,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{2,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,6), bodeplot(results.TF.up{2, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{2,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{2,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{2,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,7), bodeplot(results.TF.up{3, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{3,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{3,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{3,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,8), bodeplot(results.TF.up{3, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{3,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{3,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{3,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,9), bodeplot(results.TF.up{3, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(wind_v)
                    bodeplot(results.TF.p{3,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{3,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{3,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
sgtitle('Bode Plot - Phase - Identified Models')



