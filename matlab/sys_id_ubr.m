function results = sys_id_ubr(Data,idx_t,time_t,subj,dir,subjprop)
% Sys Ident - THIS SECTION OF CODE TAKES a while to run
clc, close all

%Flag for NaN trials
results.idx_nan = 0;

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

        %Reset Force-Displacement Identification Data
        disp_interp_t = [];
        force_interp_t = [];


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

                pausa = 1;
                
                if (sum(isnan(force_interp_t(i,:)))+sum(isnan(disp_interp_t(i,:)))) == 0
                    data_est_UPs = iddata(disp_interp_t(i,:)',force_interp_t(i,:)',Ts);
                    sysUP_s = tfest(data_est_UPs,2,0);
                    [NUM_UPs,DEN_UPs] = tfdata(sysUP_s);
                    K_est_up_s(i) = DEN_UPs{1}(3)/NUM_UPs{1}(3);
                    B_est_up_s(i) = DEN_UPs{1}(2)/NUM_UPs{1}(3);
                    M_est_up_s(i) = DEN_UPs{1}(1)/NUM_UPs{1}(3);
                    FIT_up_s(i) = sysUP_s.Report.Fit.FitPercent;

                    up_omegan_s(i) = sqrt(K_est_up_s(i)/M_est_up_s(i));
                    up_dampr_s(i) = B_est_up_s(i)/(2*sqrt(K_est_up_s(i)*M_est_up_s(i)));

                    %Estimates from Friction Model
                    [param,Fest] = sysID_2nd_Fc(Ts,disp_interp_t(i,:),force_interp_t(i,:));

                    K_fr(i) = param(3);         
                    B_fr(i) = param(2);
                    M_fr(i) = param(1);
                    Fc_fr(i) = param(4);

                    up_omegan_fr(i) = sqrt(K_fr(i)/M_fr(i));
                    up_dampr_fr(i) = B_fr(i)/(2*sqrt(K_fr(i)*M_fr(i)));
                else
                    results.idx_nan = results.idx_nan + 1;
                    K_est_up_s(i) = NaN;
                    B_est_up_s(i) = NaN;
                    M_est_up_s(i) = NaN;
                    FIT_up_s(i) = NaN;

                    up_omegan_s(i) = NaN;
                    up_dampr_s(i) = NaN;

                    K_fr(i) = NaN;
                    B_fr(i) = NaN;
                    M_fr(i) = NaN;
                    Fc_fr(i) = NaN;

                    up_omegan_fr(i) = NaN;
                    up_dampr_fr(i) = NaN;
                end
            end
        end

        %Remove NaNs
        K_est_up_s(isnan(K_est_up_s)) = [];
        B_est_up_s(isnan(B_est_up_s)) = [];
        M_est_up_s(isnan(M_est_up_s)) = [];
        FIT_up_s(isnan(FIT_up_s)) = [];
        up_omegan_s(isnan(up_omegan_s)) = [];
        up_dampr_s(isnan(up_dampr_s)) = [];
        K_fr(isnan(K_fr)) = [];
        B_fr(isnan(B_fr)) = [];
        M_fr(isnan(M_fr)) = [];
        Fc_fr(isnan(Fc_fr)) = [];
        up_omegan_fr(isnan(Fc_fr)) = [];
        up_dampr_fr(isnan(Fc_fr)) = [];
        
        %Remove Outliers 10th-90th percentile
        K_est_up_s = rmoutliers(K_est_up_s,'median');
        B_est_up_s = rmoutliers(B_est_up_s,'median');
        M_est_up_s = rmoutliers(M_est_up_s,'median');
        FIT_up_s = rmoutliers(FIT_up_s,'median');
        up_omegan_s = rmoutliers(up_omegan_s,'median');
        up_dampr_s = rmoutliers(up_dampr_s,'median');
        
        results.K_up_tr{f_sel,d_sel,:} = K_est_up_s;
        results.B_up_tr{f_sel,d_sel,:} = B_est_up_s;
        results.M_up_tr{f_sel,d_sel,:} = M_est_up_s;

        results.up_omegan_tr{f_sel,d_sel,:} = up_omegan_s;
        results.up_dampr_tr{f_sel,d_sel,:} = up_dampr_s;
        
        %Remove Outliers 10th-90th percentile
        K_fr = rmoutliers(K_fr,'median');
        B_fr = rmoutliers(B_fr,'median');
        M_fr = rmoutliers(M_fr,'median');
        Fc_fr = rmoutliers(Fc_fr,'median');
        up_omegan_fr = rmoutliers(up_omegan_fr,'median');
        up_dampr_fr  = rmoutliers(up_dampr_fr,'median');

        results.K_fr{f_sel,d_sel,:} = K_fr;
        results.B_fr{f_sel,d_sel,:} = B_fr;
        results.M_fr{f_sel,d_sel,:} = M_fr;
        results.Fc_fr{f_sel,d_sel,:} = Fc_fr;
        results.up_omegan_fr{f_sel,d_sel,:} = up_omegan_fr;
        results.up_dampr_fr{f_sel,d_sel,:} = up_dampr_fr;

        %Average and STD of Impedance and FIT over trials
        K_est_up_avg = mean(K_est_up_s);
        K_est_up_std = std(K_est_up_s);
        B_est_up_avg = mean(B_est_up_s);
        B_est_up_std = std(B_est_up_s);
        M_est_up_avg = mean(M_est_up_s);
        M_est_up_std = std(M_est_up_s);
        FIT_est_up_avg = mean(FIT_up_s);
        FIT_est_up_std = std(FIT_up_s);
        omegan_avg = mean(up_omegan_s);
        omegan_std = std(up_omegan_s);
        dampr_avg = mean(up_dampr_s);
        dampr_std = std(up_dampr_s);

        K_fr_avg = mean(K_fr);
        K_fr_std = std(K_fr);
        B_fr_avg = mean(B_fr);
        B_fr_std = std(B_fr);
        M_fr_avg = mean(M_fr);
        M_fr_std = std(M_fr);
        Fc_fr_avg = mean(Fc_fr);
        Fc_fr_std = std(Fc_fr);
        omegan_fr_avg = mean(up_omegan_fr);
        omegan_fr_std = std(up_omegan_fr);
        dampr_fr_avg = mean(up_dampr_fr);
        dampr_fr_std = std(up_dampr_fr);

        clearvars K_est_up_s B_est_up_s M_est_up_s up_omegan_s up_dampr_s FIT_up_s K_fr B_fr M_fr Fc_fr

        %Error Computation wrt to Non-Friction Model
        e.K = 100*abs((K_fr_avg-K_est_up_avg)/K_est_up_avg); 
        e.B = 100*abs((B_fr_avg-B_est_up_avg)/B_est_up_avg); 
        e.M = 100*abs((M_fr_avg-M_est_up_avg)/M_est_up_avg); 

        force_up_avg = nanmean(force_interp);
        force_up_avg_t =-(force_up_avg-force_up_avg(1));
        disp_up_avg = nanmean(disp_interp);
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

        results.FD_UP{f_sel,d_sel} = {time_new;force_interp;disp_interp_t}; %Properly Cropped Force-Displacement Data
        results.avg_FD_UP{f_sel,d_sel} = [time_new;force_up_avg;disp_up_avg_t]; %Average Cropped Force-Displacement Data
        FD_UP_avg{f_sel,d_sel} = [time_new;force_up_avg_t;disp_up_avg_t];
        results.K_up(f_sel,d_sel) = K_est_up;
        results.B_up(f_sel,d_sel) = B_est_up;
        results.M_up(f_sel,d_sel)= M_est_up;
        results.TF.up{f_sel,d_sel} = sysUP;
        results.FIT_up(f_sel,d_sel) = FIT_up;

        %Computing Identified System Response to Ballistic Release
        [yy,tt,xx] = lsim(results.TF.up{f_sel,d_sel},FD_UP_avg{f_sel,d_sel}(2,:),FD_UP_avg{f_sel,d_sel}(1,:));
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
        
        results.up_omegan_avg(f_sel,d_sel) = omegan_avg;
        results.up_omegan_std(f_sel,d_sel) = omegan_std;
        results.up_dampr_avg(f_sel,d_sel) = dampr_avg;
        results.up_dampr_std(f_sel,d_sel) = dampr_std;
        
        %Estimates from Friction Model
        results.K_fr_avg(f_sel,d_sel) = K_fr_avg;
        results.B_fr_avg(f_sel,d_sel) = B_fr_avg;
        results.M_fr_avg(f_sel,d_sel) = M_fr_avg;
        results.Fc_fr_avg(f_sel,d_sel) = Fc_fr_avg;

        results.K_fr_std(f_sel,d_sel) = K_fr_std;
        results.B_fr_std(f_sel,d_sel) = B_fr_std;
        results.M_fr_std(f_sel,d_sel) = M_fr_std;
        results.Fc_fr_std(f_sel,d_sel) = Fc_fr_std;

        results.omegan_fr_avg(f_sel,d_sel) = omegan_fr_avg;
        results.omegan_fr_std(f_sel,d_sel) = omegan_fr_std;
        results.dampr_fr_avg(f_sel,d_sel) = dampr_fr_avg;
        results.dampr_fr_std(f_sel,d_sel) = dampr_fr_std;

        %Error Between Friction and Non-Friction Models

        results.e.K{f_sel,d_sel} = e.K;
        results.e.B{f_sel,d_sel} = e.B; 
        results.e.M{f_sel,d_sel} = e.M; 


    end
end
results.M_avg_tt = mean(results.M_up_avg(:));
results.M_std_tt = std(results.M_up_avg(:));
results.M_se_tt = results.M_std_tt/sqrt(results.M_avg_tt );

results.subj_mass = subjprop.weight(subj);
results.subj_height = subjprop.height(subj);
toc
end