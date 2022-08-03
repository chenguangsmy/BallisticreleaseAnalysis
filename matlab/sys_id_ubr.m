function results = sys_id_ubr(Data,idx_t,time_t,subj,dir)
% Sys Ident - THIS SECTION OF CODE TAKES a while to run
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

        trial_l = length(Data(subj,dir,f_sel,d_sel,:,unpert));

        %Create Average Profile of Unperturbed
        time_new = linspace(0,1.25,1.25e3);
        for i = 1:trial_l
            if size(Data{subj,dir,f_sel,d_sel,i,unpert}) == [1 1]
                time_up{i,:} = time_t{subj,dir,f_sel,d_sel,i,1};
                force_interp(i,:) = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.f(1,idx_t{subj,dir,f_sel,d_sel,i,1}),time_new); %.f(2,.. for old .f(1,.. for new
                disp_interp(i,:) = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(1,idx_t{subj,dir,f_sel,d_sel,i,1}),time_new); %.x(2,.. for old .x(1,.. for new

                %Identification on Single Trial
                force_interp_t(i,:) = -(force_interp(i,:)-mean(force_interp(i,1:230)));
                disp_interp_t(i,:) = disp_interp(i,:)-mean(disp_interp(i,1:230));
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

        results.up_omegan(f_sel,d_sel) = sqrt(results.K_up(f_sel,d_sel)/results.M_up(f_sel,d_sel));
        results.up_dampr(f_sel,d_sel) = results.B_up(f_sel,d_sel)/(2*sqrt(results.K_up(f_sel,d_sel)*results.M_up(f_sel,d_sel)));
    end
end
results.M_avg_tt = mean(results.M_up_avg(:));
results.M_std_tt = std(results.M_up_avg(:));
results.M_se_tt = results.M_std_tt/sqrt(results.M_avg_tt );
toc
end