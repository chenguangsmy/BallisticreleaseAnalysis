function results = sys_id_pitt(Data,idx_t,time_t,subj)
% Sys Ident - THIS SECTION OF CODE TAKES ABOUT 1000secs to RUN!!!
clc, close all

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

        trial_l = length(Data(subj,1,f_sel,d_sel,:,pert));

        %Indices of Perturbation
        if pert > 1
            for i = 1:trial_l
                if size(Data{subj,1,f_sel,d_sel,i,pert}) == [1 1]
                    idx_p(:,i) = find(Data{subj,1,f_sel,d_sel,i,pert}.Fp(1,idx_t{subj,1,f_sel,d_sel,i,pert}) < -0.25); %.Fp(2,.. for old .Fp(1,.. for new
                end
            end
        end

        %Create Average Profile of Unperturbed
        time_new = linspace(0,1.25,1.25e3);
        for i = 1:trial_l
            if size(Data{subj,1,f_sel,d_sel,i,unpert}) == [1 1]
                time_up{i,:} = time_t{subj,1,f_sel,d_sel,i,1};
                force_interp(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,1}.f(1,idx_t{subj,1,f_sel,d_sel,i,1}),time_new); %.f(2,.. for old .f(1,.. for new
                disp_interp(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,1}.ox(1,idx_t{subj,1,f_sel,d_sel,i,1}),time_new); %.x(2,.. for old .x(1,.. for new

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
                    if size(Data{subj,1,f_sel,d_sel,i,pert}) == [1 1]
                        time_up{i,:} = time_t{subj,1,f_sel,d_sel,i,pert};
                        forcep_interp(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,pert}.f(1,idx_t{subj,1,f_sel,d_sel,i,pert}),time_new);
                        dispp_interp(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,pert}.ox(1,idx_t{subj,1,f_sel,d_sel,i,pert}),time_new);
                        force_command(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,pert}.Fp(1,idx_t{subj,1,f_sel,d_sel,i,pert}),time_new);
                    end
                else
                    if size(Data{subj,1,f_sel,d_sel,i,pert}) == [1 1]
                        time_up{i,:} = time_t{subj,1,f_sel,d_sel,i,pert};
                        forcep_interp_ph(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,pert}.f(1,idx_t{subj,1,f_sel,d_sel,i,pert}),time_new_long);
                        dispp_interp_ph(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,pert}.ox(1,idx_t{subj,1,f_sel,d_sel,i,pert}),time_new_long);
                        force_command_ph(i,:) = interp1(time_up{i,:},Data{subj,1,f_sel,d_sel,i,pert}.Fp(1,idx_t{subj,1,f_sel,d_sel,i,pert}),time_new_long);
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

                %Force and Displacement Average Difference
                disp_diff = disp_up_avg-disp_p_avg;
                disp_diff = disp_diff-disp_diff(1);
                force_diff = force_p_avg-force_up_avg;
                force_diff = force_diff - force_diff(1);

                %Force and Displacement Different over Trials
                disp_diff_t = disp_up_avg - dispp_interp;
                disp_diff_t = disp_diff_t - disp_diff_t(:,1);

                force_diff_t = forcep_interp - force_up_avg;
                force_diff_t = force_diff_t - force_diff_t(:,1);
            else
                dispp_interp_al = dispp_interp_al-dispp_interp_al(:,1); %To align data
                force_p_avg = mean(forcep_interp_al);
                disp_p_avg = mean(dispp_interp_al);
                force_c_avg = -mean(force_command_al);

                %Force and Displacement Difference
                disp_diff = disp_p_avg(1)-disp_p_avg;
                force_diff = force_p_avg;
                force_diff = force_diff - force_diff(1);

                %Force and Displacement Different over Trials
                disp_diff_t = dispp_interp_al(:,1)-dispp_interp_al;

                force_diff_t = forcep_interp_al - forcep_interp_al(:,1);
            end

            % Indices of Perturbation Resampled
            idx_p_rs = find(force_c_avg == max(force_c_avg));
            
            % Windowing of Pulse Signals
            wind_v = 30:10:300; %samples after peak // 200 fairly good estimates
            results.wind_v = wind_v; %Save windowing range for later plotting
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
                
                %Windowing Average Pulse Data
                disp_diff_ponly = disp_diff(idxmin:idxmax); %+200
                force_diff_ponly = force_diff(idxmin:idxmax);
                force_c_ponly = force_c_avg(idxmin:idxmax);

                %Windowing Single Pulse Data
                disp_diff_t_ponly = disp_diff_t(:,idxmin:idxmax);
                force_diff_t_ponly = force_diff_t(:,idxmin:idxmax);

                for i = 1:trial_l
                    if size(Data{subj,1,f_sel,d_sel,i,pert}) == [1 1]
                        %SYS ID
                        % Perturbed
                        Ts = (time_new(15)-time_new(14));
                        data_est_Pt{i,idxopt} = iddata(disp_diff_t_ponly(i,:)',force_diff_t_ponly(i,:)',Ts);
                        sysPt{i,idxopt} = tfest(data_est_Pt{i,idxopt},2,0);
                        [NUM_P,DEN_P] = tfdata(sysPt{i,idxopt});
                        K_est_pt(i,idxopt) = DEN_P{1}(3)/NUM_P{1}(3);
                        B_est_pt(i,idxopt) = DEN_P{1}(2)/NUM_P{1}(3);
                        M_est_pt(i,idxopt) = DEN_P{1}(1)/NUM_P{1}(3);
                        FIT_pt(i,idxopt) = sysPt{i,idxopt}.Report.Fit.FitPercent;
                    end
                end

                %SYS ID
                % Perturbed
                Ts = (time_new(15)-time_new(14));
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

                results.K_pt{pert-1}.p{f_sel,d_sel,idxopt} = K_est_pt(:,idxopt);
                results.K_p_avg{pert-1}.p{f_sel,d_sel,idxopt} = mean(results.K_pt{pert-1}.p{f_sel,d_sel,idxopt});
                results.K_p_std{pert-1}.p{f_sel,d_sel,idxopt} = std(results.K_pt{pert-1}.p{f_sel,d_sel,idxopt});
                results.B_pt{pert-1}.p{f_sel,d_sel,idxopt} = B_est_pt(:,idxopt);
                results.B_p_avg{pert-1}.p{f_sel,d_sel,idxopt} = mean(results.B_pt{pert-1}.p{f_sel,d_sel,idxopt});
                results.B_p_std{pert-1}.p{f_sel,d_sel,idxopt} = std(results.B_pt{pert-1}.p{f_sel,d_sel,idxopt});
                results.M_pt{pert-1}.p{f_sel,d_sel,idxopt}= M_est_pt(:,idxopt);
                results.M_p_avg{pert-1}.p{f_sel,d_sel,idxopt} = mean(results.M_pt{pert-1}.p{f_sel,d_sel,idxopt});
                results.M_p_std{pert-1}.p{f_sel,d_sel,idxopt} = std(results.M_pt{pert-1}.p{f_sel,d_sel,idxopt});
                results.TF.pt{f_sel,d_sel,pert-1,idxopt} = sysPt{:,idxopt};
                results.FIT_pt{pert-1}.p{f_sel,d_sel,idxopt} = FIT_p(:,idxopt);
                results.FIT_p_avg{pert-1}.p{f_sel,d_sel,idxopt} = mean(results.FIT_pt{pert-1}.p{f_sel,d_sel,idxopt});
                results.FIT_p_std{pert-1}.p{f_sel,d_sel,idxopt} = std(results.FIT_pt{pert-1}.p{f_sel,d_sel,idxopt});

                [yy,tt,xx] = lsim(results.TF.p{f_sel,d_sel,pert-1,idxopt},results.avg_FD_UP{f_sel,d_sel}(2,:),results.avg_FD_UP{f_sel,d_sel}(1,:));
                results.tt_modP{pert-1}.p{f_sel,d_sel,idxopt} = tt;
                results.disp_modP{pert-1}.p{f_sel,d_sel,idxopt} = yy;
                
                idxopt = idxopt+1;
            end
        end
        
        results.up_omegan(f_sel,d_sel) = sqrt(results.K_up(f_sel,d_sel)/results.M_up(f_sel,d_sel));
        results.up_dampr(f_sel,d_sel) = results.B_up(f_sel,d_sel)/(2*sqrt(results.K_up(f_sel,d_sel)*results.M_up(f_sel,d_sel)));
    end
end
results.M_avg_tt = mean(results.M_up_avg(:));
results.M_std_tt = std(results.M_up_avg(:));
results.M_se_tt = results.M_std_tt/sqrt(results.M_avg_tt );
toc
end