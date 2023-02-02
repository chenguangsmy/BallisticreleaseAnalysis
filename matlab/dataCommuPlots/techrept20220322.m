% Techreport for 2022-03-22 
%
% The code here generates the plot serve for the docuemnt. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the figures index: 
% fig1:     The stiffness measurement illustration; 
% fig2:     The damping measurement illustration; 
% fig3:     The inertia measurement illustration; 
% fig4:     The inertia measurement illustration on release time; 
% fig5:     The stiffness measurment cross condition: subj+ spring;
% fig6:     The stiffnes measurement cross time, subj;
% fig7:     The stiffness measurement across time, spring;
% fig8:     The damping measurement cross condition, subj + spring; 
% fig9:     The damping measurement across time, subj;
% fig10:    The damping measurement across time, spring;
% fig11:    The inertia measurement cross condition, subj + spring; 
% fig12:    The inertia measurement across time, subj;
% fig13:    The inertia measurement across time, spring;
% fig14:    The inertia measurement cross condition (using release epoc),
%           subj + spring;
% sup-file1:    The pdf file generate for all subject measurments(K+D+I)
% sup-file2:    The pdf file generate for all spring measurments(K+D+I)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1:     The stiffness measurement illustration; 
% also, generate all the K, D, I for subject in this plot. 
if(0) % use the method of nearest point 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949_tmp.mat', 'data'); % this one nearest subtracted
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data'); %Average subtracted
fh1 = figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
mk_big = 10;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 1;
for fce_i = 1:size(data,2)
    % for fce_i = 3
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 1:8
%             fh(pi,1) = figure(); hold on;
            to_plot.ifplot = (fce_i == to_plot.fce_i && ...
                dist_i == to_plot.dist_i && ...
                pi == to_plot.pi);
%             t0_plot.ifplot = 1;
            if (to_plot.ifplot)
                fh1 = figure('unit', 'inch', 'position', [0 0 6 5]);
            end
            clear axh;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8,1),7,8); % 1 no -ert and 5 pert
            celltmp2 = reshape(data(1,fce_i,dist_i,:,1:8,2),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            fp_avg = fp_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            if (to_plot.ifplot)
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(4,2,1);
                plot(t_grids, fp_avg);
                axh(3) = subplot(4,2,3);
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(4,2,5);
                plot(t_grids, x_avg);
                axh(7) = subplot(4,2,7);
                plot(t_grids, f_avg ./ (x_avg - x_settled));
                yline(psudo_stiffness1, 'linewidth', 2);
                %             ylim([0 2000]);
            end
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                idx_release0= find(celltmp2{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                t0= celltmp2{ti,pi}.t - celltmp2{ti,pi}.t(idx_release0(1));
                
                if (to_plot.ifplot)
                    subplot(axh(1)); hold on;
                    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    
                    subplot(axh(3)); hold on;
                    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(5)); hold on;
                    x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                    plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                end
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                idx_t0= find(t0>=t_interest(1) & t0<=t_interest(2));
%                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % 
                fp_dat0= interp1(t0(idx_t0), celltmp2{ti,pi}.Fp(2,idx_t0), t_grids, 'linear', 'extrap');   %
                x_dat0= interp1(t0(idx_t0), celltmp2{ti,pi}.x(2,idx_t0), t_grids, 'linear', 'extrap');     %
                v_dat0= interp1(t0(idx_t0), celltmp2{ti,pi}.v(2,idx_t0), t_grids, 'linear', 'extrap');     %
                f_dat0=-interp1(t0(idx_t0), celltmp2{ti,pi}.f(2,idx_t0), t_grids, 'linear', 'extrap');     % 

                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                v_filter0= filter(b,a,v_dat0);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                a_dat0= [0 diff(v_filter0)./diff(t_grids)];
                
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_dat0;
                v_net = v_dat - v_dat0;
                a_net = a_dat - a_dat0; 
                fp_net = fp_dat - fp_dat0;
                f_net = f_dat - f_dat0;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_dat0)./(x_dat - x_dat0);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_dat0)./(v_dat - v_dat0);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
%                 a_net_tmp = a_dat0; a_net_tmp(1:fp_max_idx) = 0;
                a_net_tmp = a_dat0; a_net_tmp(fp_dat==0) = 0;
%                 [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                [~, a_peak_idx] = min(a_net);    % only take after perturbation part
%                 m_est = (f_dat - f_dat0)./(a_dat - a_dat0);
                m_est = f_net./a_net;
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                if (to_plot.ifplot)
                    axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                    plot(t_grids, fp_dat - fp_dat0, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(4,2,4); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_dat0), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
%                     plot(t_grids(a_peak_idx),-f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(6) = subplot(4,2,6); hold on;% subtracted F
                    plot(t_grids, x_dat - x_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids, a_net, 'color', color_arr(4+dist_i,:));
                    plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
%                     plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(4,2,8); hold on;
                    plot(t_grids, (f_dat - f_dat0)./(x_dat - x_dat0), 'color', color_arr(4+dist_i,:));
%                     plot(t_grids, f_net./a_net, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big);
                    plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big);
                end
            end
            
            if (to_plot.ifplot)
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(7:8), 'y');
                    yline(axh(8),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                ylim(axh(7), [0 1000]);
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
%                 sgtitle(['Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);
                sgtitle(['Force 15N dist 2.5cm, 1st pulse']);
                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'x');
                ylabel(axh(7), 'df/dx');
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'dx');
                    ylabel(axh(8), 'df/dx');
                catch
                end
            end
            %          saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = psudoK_mat;
        psudoD_cell{fce_i,dist_i} = psudoD_mat;
        psudoI_cell{fce_i,dist_i} = psudoI_mat;
        %         close all;
        
    end
    
end
saveas(fh1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure1.png');
end

%% Get the subtraction from the averaged trials
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data'); %Average subtracted
fh1 = figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
mk_big = 10;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 1;
for fce_i = 1:size(data,2)
    % for fce_i = 3
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 1:8
%             fh(pi,1) = figure(); hold on;
            to_plot.ifplot = (fce_i == to_plot.fce_i && ...
                dist_i == to_plot.dist_i && ...
                pi == to_plot.pi);
%             t0_plot.ifplot = 1;
            if (to_plot.ifplot)
                fh1 = figure('unit', 'inch', 'position', [0 0 6 5]);
            end
            clear axh;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            fp_avg = fp_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            if (to_plot.ifplot)
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(4,2,1);
                plot(t_grids, fp_avg);
                axh(3) = subplot(4,2,3);
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(4,2,5);
                plot(t_grids, x_avg);
                axh(7) = subplot(4,2,7);
                plot(t_grids, f_avg ./ (x_avg - x_settled));
                yline(psudo_stiffness1, 'linewidth', 2);
                %             ylim([0 2000]);
            end
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
%                 idx_release0= find(celltmp2{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 t0= celltmp2{ti,pi}.t - celltmp2{ti,pi}.t(idx_release0(1));
                
                if (to_plot.ifplot)
                    subplot(axh(1)); hold on;
                    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
%                     plot(t0, celltmp2{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    plot(t_grids, fp_avg, 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    
                    subplot(axh(3)); hold on;
                    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t_grids,-f_avg, 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
%                     plot(t0, celltmp2{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(5)); hold on;
                    x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                    plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                     plot(t0, celltmp2{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    plot(t_grids, x_avg, 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                end
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % 
                fp_dat0 = fp_avg;
                x_dat0 = x_avg;
                v_dat0 = v_avg;
                f_dat0 = f_avg;
                
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                v_filter0= filter(b,a,v_dat0);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                a_dat0= [0 diff(v_filter0)./diff(t_grids)];
                
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_dat0;
                v_net = v_dat - v_dat0;
                a_net = a_dat - a_dat0; 
                fp_net = fp_dat - fp_dat0;
                f_net = f_dat - f_dat0;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_dat0)./(x_dat - x_dat0);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_dat0)./(v_dat - v_dat0);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
%                 a_net_tmp = a_dat0; a_net_tmp(1:fp_max_idx) = 0;
                a_net_tmp = a_dat0; a_net_tmp(fp_dat==0) = 0;
%                 [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                [~, a_peak_idx] = min(a_net);    % only take after perturbation part
%                 m_est = (f_dat - f_dat0)./(a_dat - a_dat0);
                m_est = f_net./a_net;
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                if (to_plot.ifplot)
                    axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                    plot(t_grids, fp_dat - fp_dat0, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(4,2,4); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_dat0), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
%                     plot(t_grids(a_peak_idx),-f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(6) = subplot(4,2,6); hold on;% subtracted F
                    plot(t_grids, x_dat - x_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids, a_net, 'color', color_arr(4+dist_i,:));
                    plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
%                     plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(4,2,8); hold on;
                    plot(t_grids, (f_dat - f_dat0)./(x_dat - x_dat0), 'color', color_arr(4+dist_i,:));
%                     plot(t_grids, f_net./a_net, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big);
                    plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big);
                end
            end
            
            if (to_plot.ifplot)
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(7:8), 'y');
                    yline(axh(8),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                ylim(axh(7), [0 1000]);
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
%                 sgtitle(['Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);
                sgtitle(['Force 15N dist 2.5cm, 1st pulse']);
                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'x');
                ylabel(axh(7), 'df/dx');
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'dx');
                    ylabel(axh(8), 'df/dx');
                catch
                end
            end
            %          saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK_cell{fce_i,dist_i} = psudoK_mat;
        psudoD_cell{fce_i,dist_i} = psudoD_mat;
        psudoI_cell{fce_i,dist_i} = psudoI_mat;
        %         close all;
        
    end
    
end
saveas(fh1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure1.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate spring data for the following statistic comparing
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
F_list = [15, 20, 25];
K_list = [640, 320, 160];
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK1_cell = cell(3,3); % for spring 
psudoD1_cell = cell(3,3);
psudoI1_cell = cell(3,3);
to_plot.fce_i = 0;  % not plot here
to_plot.dist_i = 0;
to_plot.pi = 0;
to_plot.ifplot = 0;
for fce_i = 1:size(data,2)
    % for fce_i = 3
    for dist_i = 1:size(data,3) % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        psudoK_mat = zeros(5,7);
        psudoD_mat = zeros(5,7);
        psudoI_mat = zeros(5,7);
        for pi = 1:6
            %             fh(pi,1) = figure(); hold on;
            to_plot.ifplot = (fce_i == to_plot.fce_i && ...
                dist_i == to_plot.dist_i && ...
                pi == to_plot.pi);
            if (to_plot.ifplot)
                fh1 = figure('unit', 'inch', 'position', [0 0 6 5]);
            end
            clear axh;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
%                 x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            if (to_plot.ifplot)
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(4,2,1);
                plot(t_grids, fp_avg);
                axh(3) = subplot(4,2,3);
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(4,2,5);
                plot(t_grids, x_avg);
                axh(7) = subplot(4,2,7);
                plot(t_grids, f_avg ./ (x_avg - x_settled));
                yline(psudo_stiffness1, 'linewidth', 2);
                %             ylim([0 2000]);
            end
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                if (to_plot.ifplot)
                    subplot(axh(1)); hold on;
                    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(3)); hold on;
                    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                    %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(5)); hold on;
                    x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                    plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                end
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
%                 x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_avg;
                v_net = v_dat - v_avg;
                a_net = a_dat - a_avg; 
                fp_net = fp_dat - fp_avg;
                f_net = f_dat - f_avg;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_avg)./(x_dat - x_avg);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_avg)./(v_dat - v_avg);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                m_est = (f_dat - f_avg)./(a_dat - a_avg);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                if (to_plot.ifplot)
                    axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                    plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(4,2,4); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_avg), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                    axh(6) = subplot(4,2,6); hold on;% subtracted F
                    plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                    plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(4,2,8); hold on;
                    plot(t_grids, (f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big);
                end
            end
            
            if (to_plot.ifplot)
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(7:8), 'y');
                    yline(axh(8),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                ylim(axh(7), [0 1000]);
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
%                 sgtitle(['Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);
                sgtitle(['Force 15N dist 2.5cm, 1st pulse']);
                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'x');
                ylabel(axh(7), 'df/dx');
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'dx');
                    ylabel(axh(8), 'df/dx');
                catch
                end
            end
            %          saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
        psudoK1_cell{fce_i,dist_i} = psudoK_mat;
        psudoD1_cell{fce_i,dist_i} = psudoD_mat;
        psudoI1_cell{fce_i,dist_i} = psudoI_mat;
        
%         psudoK1O_cell{fce_i,dist_i} = psudoK_mat; % means for optotrak measurements
        %         close all;
        
    end
    
end


save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat',...
    'psudoK_cell','psudoD_cell','psudoI_cell', ...
    'psudoK1_cell','psudoD1_cell','psudoI1_cell', '-append');
% save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat',...
%     'psudoK1O_cell', '-append');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig2:     The damping measurement illustration; 
% only do plot but not generate data 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949_tmp.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 0;
for fce_i = to_plot.fce_i
    % for fce_i = 3
    for dist_i = to_plot.dist_i % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 1:8
            %             fh(pi,1) = figure(); hold on;
            to_plot.ifplot = (fce_i == to_plot.fce_i && ...
                dist_i == to_plot.dist_i && ...
                pi == to_plot.pi);
            if (to_plot.ifplot)
                fh2 = figure('unit', 'inch', 'position', [0 0 6 5]);
            end
            clear axh;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8,1),7,8); % 1 no -ert and 5 pert
            celltmp2 = reshape(data(1,fce_i,dist_i,:,1:8,2),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));            
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            if (to_plot.ifplot)
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(4,2,1);
%                 plot(t_grids, fp_avg);
                axh(3) = subplot(4,2,3);
%                 plot(t_grids,-f_avg);
                axh(5) = subplot(4,2,5);
%                 plot(t_grids, v_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(7) = subplot(4,2,7);
                plot(t_grids,-f_avg ./ (v_avg - v_avg(end)));
                yline(psudo_stiffness1, 'linewidth', 2);
%                             ylim([0 2000]);
            end
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                idx_release0= find(celltmp2{ti,pi}.ts == 5);
                t0 = celltmp2{ti,pi}.t - celltmp2{ti,pi}.t(idx_release0(1));
                if (to_plot.ifplot)
                    subplot(axh(1)); hold on;
                    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    subplot(axh(3)); hold on;
                    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    
                    subplot(axh(5)); hold on;
                    v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                    plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                end
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                idx_t0= find(t0>=t_interest(1) & t0<=t_interest(2));
%                 length(idx_t);
                

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % 
                fp_dat0= interp1(t0(idx_t0), celltmp2{ti,pi}.Fp(2,idx_t0), t_grids, 'linear', 'extrap');    %
                x_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.x(2,idx_t0), t_grids, 'linear', 'extrap');     %
                v_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.v(2,idx_t0), t_grids, 'linear', 'extrap');     %
                f_dat0 =-interp1(t0(idx_t0), celltmp2{ti,pi}.f(2,idx_t0), t_grids, 'linear', 'extrap');     % 
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                v_filter0 = filter(b,a,v_dat);
                a_dat0 = [0 diff(v_filter)./diff(t_grids)];
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_dat0;
                v_net = v_dat - v_dat0;
                a_net = a_dat - a_dat0; 
                fp_net = fp_dat - fp_dat0;
                f_net = f_dat - f_dat0;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_dat0)./(x_dat - x_dat0);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_dat0)./(v_dat - v_dat0);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                m_est = (f_dat - f_dat0)./(a_dat - a_dat0);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                % plot the subtraction in other panels
                if (to_plot.ifplot)
                    axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                    plot(t_grids, fp_dat - fp_dat0, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(4,2,4); hold on;% subtracted f
                    plot(t_grids, -(f_dat - f_dat0), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(v_peak_idx), -f_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(6) = subplot(4,2,6); hold on;% subtracted v
                    plot(t_grids, v_dat - v_dat0, 'color', color_arr(4+dist_i,:));
                    plot(t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(4,2,8); hold on;
                    plot(t_grids,-(f_dat - f_dat0)./(v_dat - v_dat0), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(v_peak_idx), d_est_pt, 'marker', '.', 'markersize', mk_big);
                end
            end
            
            if (to_plot.ifplot)
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(7:8), 'y');
%                     yline(axh(8),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                ylim(axh(7), [0 200]);
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
%                 sgtitle(['Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);
                sgtitle(['Force 15N dist 2.5cm, 1st pulse']);
                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'v');
                ylabel(axh(7), '-df/dv');
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'dv');
                    ylabel(axh(8), '-df/dv');
                catch
                end
            end
            %          saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
%         psudoK_cell{fce_i,dist_i} = psudoK_mat;
%         psudoD_cell{fce_i,dist_i} = psudoD_mat;
%         psudoI_cell{fce_i,dist_i} = psudoI_mat;
        %         close all;
        
    end
    
end
saveas(fh2, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure2.png');
% save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat', 'psudoK_cell','psudoD_cell','psudoI_cell');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig3:     The inertia measurement illustration; 
% only do plot but not generate data 

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949_tmp.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 0;
for fce_i = to_plot.fce_i
    % for fce_i = 3
    for dist_i = to_plot.dist_i % for each spring
        %         subplot(3,3, (fce_i-1)*3+dist_i); hold on;
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 1:8
            %             fh(pi,1) = figure(); hold on;
            to_plot.ifplot = (fce_i == to_plot.fce_i && ...
                dist_i == to_plot.dist_i && ...
                pi == to_plot.pi);
            if (to_plot.ifplot)
                fh3 = figure('unit', 'inch', 'position', [0 0 6 5]);
            end
            clear axh;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8,1),7,8); % 1 no -ert and 5 pert
            celltmp2 = reshape(data(1,fce_i,dist_i,:,1:8,2),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));            
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            if (to_plot.ifplot)
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(4,2,1);
%                 plot(t_grids, fp_avg);
                axh(3) = subplot(4,2,3);
%                 plot(t_grids,-f_avg);
                axh(5) = subplot(4,2,5);
%                 plot(t_grids, v_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(7) = subplot(4,2,7);
                plot(t_grids,-f_avg ./ (v_avg - v_avg(end)));
                yline(psudo_stiffness1, 'linewidth', 2);
                %             ylim([0 2000]);
            end
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                idx_release0= find(celltmp2{ti,pi}.ts == 5);
                t0 = celltmp2{ti,pi}.t - celltmp2{ti,pi}.t(idx_release0(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                idx_t0= find(t0>=t_interest(1) & t0<=t_interest(2));
%                 length(idx_t);
                

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % 
                fp_dat0= interp1(t0(idx_t0), celltmp2{ti,pi}.Fp(2,idx_t0), t_grids, 'linear', 'extrap');    %
                x_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.x(2,idx_t0), t_grids, 'linear', 'extrap');     %
                v_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.v(2,idx_t0), t_grids, 'linear', 'extrap');     %
                f_dat0 =-interp1(t0(idx_t0), celltmp2{ti,pi}.f(2,idx_t0), t_grids, 'linear', 'extrap');     % 
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                v_filter0 = filter(b,a,v_dat0);
                a_dat0 = [0 diff(v_filter0)./diff(t_grids)];
                %                 linkaxes(axh(1:3:5), 'x');
                
                if (to_plot.ifplot)
                    subplot(axh(1)); hold on;
                    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    subplot(axh(3)); hold on;
                    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                    plot(t0, celltmp2{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                    
                    subplot(axh(5)); hold on;
                    v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                    plot(t_grids, a_dat, 'color', color_arr(4+dist_i,:));
                    plot(t_grids, a_dat0, 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                end
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_dat0;
                v_net = v_dat - v_dat0;
                a_net = a_dat - a_dat0; 
                fp_net = fp_dat - fp_dat0;
                f_net = f_dat - f_dat0;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_dat0)./(x_dat - x_dat0);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_dat0)./(v_dat - v_dat0);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                m_est = (f_dat - f_dat0)./(a_dat - a_dat0);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                % plot the subtraction in other panels
                if (to_plot.ifplot)
                    axh(2) = subplot(4,2,2); hold on; % subtracted Fp
                    plot(t_grids, fp_dat - fp_dat0, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(4,2,4); hold on;% subtracted f
                    plot(t_grids, -(f_dat - f_dat0), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(a_peak_idx), -f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(6) = subplot(4,2,6); hold on;% subtracted v
                    plot(t_grids, a_dat - a_dat0, 'color', color_arr(4+dist_i,:));
                    plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(4,2,8); hold on;
                    plot(t_grids,-(f_dat - f_dat0)./(a_dat - a_dat0), 'color', color_arr(4+dist_i,:));
                    plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big);
                end
            end
            
            if (to_plot.ifplot)
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(7:8), 'y');
%                     yline(axh(8),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                ylim(axh(7), [0 200]);
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
%                 sgtitle(['Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);
                sgtitle(['Force 15N dist 2.5cm, 1st pulse']);
                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'a');
                ylabel(axh(7), '-df/da');
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'da');
                    ylabel(axh(8), 'df/da');
                catch
                end
            end
            %          saveas(gcf, ['sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/psudoStiffness200ms_subj' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i)) 'pert' num2str(pi-1) '.png']);
        end
%         psudoK_cell{fce_i,dist_i} = psudoK_mat;
%         psudoD_cell{fce_i,dist_i} = psudoD_mat;
%         psudoI_cell{fce_i,dist_i} = psudoI_mat;
        %         close all;
        
    end
    
end
saveas(fh3, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure3.png');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig4:     The inertia measurement illustration on release time; 
close all; clear; clc; 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
fh4 = figure('unit', 'inch', 'position', [0 0 3 3]);
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
% for fce_i = 1:size(data,3)
%     for dist_i = 1:size(data,4) % for each spring
f_ratio_cell = cell(3,3);
for fce_i = 1%:size(data,2)
    for dist_i = 1%:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); %
        pi = 1; % only see the unperturbed trials
%         fh(pi,1) = figure(); hold on;
        
        celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
        idx_release = find(celltmp1{1,pi}.ts == 5);
        t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
        
        
        f_ratio = zeros(1,size(celltmp1,1));
        for ti = 1%:1:size(celltmp1,1)
            if isempty(celltmp1{ti,1})
                continue;
            end
            idx_release = find(celltmp1{ti,1}.ts == 5);
            t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
            idx_t = find(t>=t_interest(1) & t<=t_interest(2));
            % intropolate (x, f, Fp) to t_grids
            
            x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
            
            ifplot = 1; % controls whether plot or not
            t_idx = find(t_grids>0 & t_grids < 0.1); 
            [f_pk_neg,locs1] = findpeaks(-f_dat(t_idx), 'NPeaks', 1);
            [f_pk_pos,locs2] = findpeaks(f_dat(t_idx), 'NPeaks', 1);
            if (ifplot)                
%                 axh(3) = subplot(3,1,3); 
                hold on;
                plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                plot(t_grids, f_dat, 'r', 'Marker', '.'); hold on;
                
                plot(t_grids(t_idx(locs1)),-f_pk_neg, 'marker', '.', 'markersize', 10);
                plot(t_grids(t_idx(locs2)), f_pk_pos, 'marker', '.', 'markersize', 10);
                
                xlim([-0.02 0.18]);
                ylim([-2 20]);
            end    
            f0_bef = mean(f_dat(t_grids<0)); 
%             f0_aft = mean(f_dat(t_grids>0.01 & t_grids<0.04));
            f0_aft = 1/2*(f_pk_pos + (-f_pk_neg));
            yline([f0_bef], 'linewidth', 2);
            yline([f0_aft], 'linewidth', 2);
            f_ratio(ti) = f0_bef/f0_aft;
            
           f_ratio_cell{fce_i,dist_i} = f_ratio;
        end 
    end    
end
xlabel('t (s)');
ylabel('f (N)');
title('censored force at release point');
saveas(fh4, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure4.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig5:     The stiffness measurment cross condition: subj+ spring;
% corss condition plot 
% todo ... need to change to remove all the 'outlair values'
clear; clc; 
% close all;
F_list = [15 20 25];
K_list = [2.5 5.0 7.5];
K_list_spr = [640 320 160];
ref_list = F_list' * (1./(K_list/100));
k_stfcoef_subj = 1/3.62; % after relase, the number is corresponding to 3kg hand
k_stfcoef_spr = 13/20;
color_arr = colormap('lines');
% close all;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
dat_mean_cc_subj = zeros(3,3); 
dat_std_cc_subj  = zeros(3,3); 
dat_mean_cc_spr = zeros(3,3); 
dat_std_cc_spr  = zeros(3,3); 
for fce_i = 1:3
     for dist_i = 1:3 
         subj_K = psudoK_cell{fce_i,dist_i}(:); % all perturbation time 
%         subj_K = psudoK_cell{fce_i,dist_i}(1,:);  % the first pert time
         % get rid of outlairer from K 
         dat_arr_outlairidx = isoutlier(subj_K);
         subj_K(dat_arr_outlairidx) = nan; 
         subj_K(subj_K == 0) = nan; 
         dat_arr_nan = reshape(subj_K, size(psudoK_cell{fce_i,dist_i},1), size(psudoK_cell{fce_i,dist_i},2));
         % if I want to only select the first pulse here
         dat_arr_nan = dat_arr_nan(1,:);
%          dat_arr_outlairidx = zeros(1:size(subj_K)); % manuly set to avoid "remove outlaiers:
%          dat_mean_cc_subj(fce_i,dist_i) = mean(subj_K(~dat_arr_outlairidx));
%          dat_std_cc_subj(fce_i,dist_i) = std(subj_K(~dat_arr_outlairidx));
         dat_mean_cc_subj(fce_i,dist_i) = nanmean(dat_arr_nan(:));
         dat_std_cc_subj(fce_i,dist_i) = nanstd(dat_arr_nan(:));

         
         spr_K = psudoK1_cell{fce_i,dist_i}(:); 
         dat_mean_cc_spr(fce_i,dist_i) = mean(spr_K);
         dat_std_cc_spr(fce_i,dist_i) = std(spr_K);
     end
end

fh5 = figure('unit', 'inch', 'position', [0 0 6 3]); 
axh(1) = subplot(1,2,1); % subject 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc_subj(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc_subj(:,dist_i), dat_std_cc_subj(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    
    plot(15:5:25, ref_list(:,dist_i) * k_stfcoef_subj, 'color', [0.8 0.8 0.8], 'linewidth', 2);
end
title('subject');
xlabel('force threshold'); 
ylabel('K (N/m)'); 
legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
axh(2) = subplot(1,2,2); % spring 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc_spr(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc_spr(:,dist_i), dat_std_cc_spr(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    
    plot(15:5:25, K_list_spr(dist_i)* k_stfcoef_spr*[1 1 1], 'color', [0.8 0.8 0.8], 'linewidth', 2);
end
title('springs');
xlabel('force threshold');
legend(lnh, {'640N/m', '320N/m', '160N/m'});

linkaxes(axh, 'xy');
ylim([-100 800]);
xlim([13 27])
sgtitle('cross condition stiffness estimation')
saveas(fh5, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure5_subj1stpulse.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig6:     The stiffnes measurement cross time, subj;
color_arr = colormap('lines');
% close all;
F_list = [15 20 25];
K_list = [2.5 5.0 7.5];
ref_list = F_list' * (1./(K_list/100));
k_stfcoef = 1/3.62; % after relase, the number is corresponding to 3kg hand
fh6 = figure('unit', 'inch', 'position', [0 0 5 5]); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
peak_time = [0.1:0.025:0.25]; % subject
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoK_cell{fce_i,dist_i}(:);
       
        dat_arr_outlairidx = isoutlier(dat);
        dat(dat_arr_outlairidx) = nan; 
        dat(dat == 0) = nan; 
        dat = reshape(dat, size(psudoK_cell{fce_i,dist_i},1), size(psudoK_cell{fce_i,dist_i},2));
        % if I only want the first time 
%         dat(2:end,:) = nan;
        % end if
        dat_mean = nanmean(dat, 2); % each row a perturbation 
        dat_std = nanstd(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        if (fce_i == 3)
            xlabel('pert time'); 
        end
        if (dist_i == 1)
            ylabel('K (N/m)');
        end
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
        yline(k_stfcoef*ref_list(fce_i,dist_i));
    end
end
sgtitle('Subject stiffness after release, dF/dx');
linkaxes(axh);
ylim([-200 400]);
xlim([0.1 0.25]);
saveas(fh6, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure6_subavg.png');
% saveas(fh6, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure6_1stpulse.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig7:     The stiffness measurement across time, spring;
color_arr = colormap('lines');
close all;
F_list = [15 20 25];
K_list = [640 320 160];
fh7 = figure('unit', 'inch', 'position', [0 0 5 5]); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
peak_time = [0.15:0.1:0.55]; % spring
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoK1_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        if (fce_i == 3)
            xlabel('pert time'); 
        end
        if (dist_i == 1)
            ylabel('K (N/m)');
        end
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
    end
end
sgtitle('Spring stiffness after release, dF/dx');
linkaxes(axh);
ylim([-200 1000]);
saveas(fh7, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure7.png');

%% fig7.1:     The stiffness measurement across time, spring;
% using optotrak (dashed) and wam (solid) measurements
color_arr = colormap('lines');
close all;
F_list = [15 20 25];
K_list = [640 320 160];
k_stfcoef = 13/20;
fh7 = figure('unit', 'inch', 'position', [0 0 5 5]); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
peak_time = [0.15:0.1:0.55]; % spring
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        hold on;
        dat = psudoK1_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        
        dato= psudoK1O_cell{fce_i,dist_i};
        dato_mean = mean(dato,2); 
        dato_std = std(dato, [], 2);
        
%         plot(peak_time, dato_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:), 'LineStyle', ':');
        d = errorbar(peak_time, dato_mean', dato_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:), 'LineStyle', ':');
        d.Bar.LineStyle = 'dotted';
        
%         plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        if (fce_i == 3)
            xlabel('pert time'); 
        end
        if (dist_i == 1)
            ylabel('K (N/m)');
        end
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
        yline(K_list(dist_i)*k_stfcoef);
    end
end
sgtitle('Spring stiffness after release, dF/dx');
linkaxes(axh);
ylim([-50 600]);
saveas(fh7, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure7_1.png');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig8:     The damping measurement cross condition, subj + spring; 
% corss condition plot 
clear; close all; clc; 
color_arr = colormap('lines');
close all;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
dat_mean_cc_subj = zeros(3,3); 
dat_std_cc_subj  = zeros(3,3); 
dat_mean_cc_spr = zeros(3,3); 
dat_std_cc_spr  = zeros(3,3); 
for fce_i = 1:3
     for dist_i = 1:3 
%          subj_D = psudoD_cell{fce_i,dist_i}(:); 
        subj_D = psudoD_cell{fce_i,dist_i}(1,:); % the first pulse 
         % get rid of outlairer from K 
%          dat_arr_outlairidx = isoutlier(subj_D);
         dat_mean_cc_subj(fce_i,dist_i) = mean(subj_D);
         dat_std_cc_subj(fce_i,dist_i) = std(subj_D);
         
         spr_D = psudoD1_cell{fce_i,dist_i}(:); 
         dat_mean_cc_spr(fce_i,dist_i) = mean(spr_D);
         dat_std_cc_spr(fce_i,dist_i) = std(spr_D);
     end
end

fh8 = figure('unit', 'inch', 'position', [0 0 6 3]); 
axh(1) = subplot(1,2,1); % subject 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc_subj(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc_subj(:,dist_i), dat_std_cc_subj(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
title('subject');
xlabel('force threshold'); 
ylabel('D (Ns/m)'); 
legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
axh(2) = subplot(1,2,2); % spring 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc_spr(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc_spr(:,dist_i), dat_std_cc_spr(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
title('springs');
xlabel('force threshold');
legend(lnh, {'640N/m', '320N/m', '160N/m'});


linkaxes(axh, 'xy');
ylim([-10 100]);
xlim([13 27])
sgtitle('cross condition damping estimation')
saveas(fh8, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure8.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig9:     The damping measurement across time, subj;
color_arr = colormap('lines');
close all;
F_list = [15 20 25];
K_list = [2.5 5.0 7.5];
fh9 = figure('unit', 'inch', 'position', [0 0 5 5]); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
peak_time = [0.1:0.025:0.25]; % subject
clear axh;
% There is little 'outlaiers' in the damping estimation. 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoD_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        if (fce_i == 3)
            xlabel('pert time'); 
        end
        if (dist_i == 1)
            ylabel('D (Ns/m)');
        end
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('Subject damping after release, dF/dv');
linkaxes(axh);
ylim([-10 100]);
saveas(fh9, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure9.png');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig10:    The damping measurement across time, spring;
color_arr = colormap('lines');
close all;
F_list = [15 20 25];
K_list = [640 320 160];
fh10 = figure('unit', 'inch', 'position', [0 0 5 5]); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
peak_time = [0.15:0.1:0.55]; % spring
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoD1_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        if (fce_i == 3)
            xlabel('pert time'); 
        end
        if (dist_i == 1)
            ylabel('D (Ns/m)');
        end
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
    end
end
sgtitle('Spring damping after release, dF/dv');
linkaxes(axh);
ylim([-10 100]);
saveas(fh10, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure10.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig11:    The inertia measurement cross condition, subj + spring; 
% corss condition plot 
clear; close all; clc; 
color_arr = colormap('lines');
close all;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
dat_mean_cc_subj = zeros(3,3); 
dat_std_cc_subj  = zeros(3,3); 
dat_mean_cc_spr = zeros(3,3); 
dat_std_cc_spr  = zeros(3,3); 
for fce_i = 1:3
     for dist_i = 1:3 
         subj_I = psudoI_cell{fce_i,dist_i}(:); 
         % get rid of outlairer from K 
         dat_arr_outlairidx = isoutlier(subj_I);
         subj_I(dat_arr_outlairidx) = nan;
         subj_I = reshape(subj_I, size(psudoI_cell{fce_i,dist_i},1), size(psudoI_cell{fce_i,dist_i},2));
         subj_I = subj_I(1,:); % first pert
%          subj_I = subj_I(:);
         dat_mean_cc_subj(fce_i,dist_i) = nanmean(subj_I);
         dat_std_cc_subj(fce_i,dist_i) = nanstd(subj_I);
         
         spr_I = psudoI1_cell{fce_i,dist_i}(:); 
         dat_mean_cc_spr(fce_i,dist_i) = mean(spr_I);
         dat_std_cc_spr(fce_i,dist_i) = std(spr_I);
     end
end

fh11 = figure('unit', 'inch', 'position', [0 0 6 3]); 
axh(1) = subplot(1,2,1); % subject 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc_subj(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc_subj(:,dist_i), dat_std_cc_subj(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
title('subject');
xlabel('force threshold'); 
ylabel('I (kg)'); 
legend(lnh, {'2.5cm', '5.0cm', '7.5cm'});
axh(2) = subplot(1,2,2); % spring 
hold on;
clear lnh
for dist_i = 1:3
    lnh(dist_i) = plot(15:5:25, dat_mean_cc_spr(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
    errorbar(15:5:25, dat_mean_cc_spr(:,dist_i), dat_std_cc_spr(:,dist_i), 'linewidth', 2, 'color', color_arr(4+dist_i,:));
end
title('springs');
xlabel('force threshold');
legend(lnh, {'640N/m', '320N/m', '160N/m'});

linkaxes(axh, 'xy');
ylim([-0 10]);
xlim([13 27])
sgtitle('cross condition inertia estimation')
saveas(fh11, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure11.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig12:    The inertia measurement across time, subj;
color_arr = colormap('lines');
close all;
F_list = [15 20 25];
K_list = [2.5 5.0 7.5];
fh12 = figure('unit', 'inch', 'position', [0 0 5 5]); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
peak_time = [0.1:0.025:0.25]; % subject
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoI_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        if (fce_i == 3)
            xlabel('pert time'); 
        end
        if (dist_i == 1)
            ylabel('I (kg)');
        end
        title(['fce' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
    end
end
sgtitle('Subject mass after release, dF/da');
linkaxes(axh);
ylim([-2 10]);
saveas(fh12, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure12.png');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% fig13:    The inertia measurement across time, spring;
color_arr = colormap('lines');
close all;
F_list = [15 20 25];
K_list = [640 320 160];
fh13 = figure('unit', 'inch', 'position', [0 0 5 5]); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/psudoImpPredict.mat')
peak_time = [0.15:0.1:0.55]; % spring
clear axh; 
for fce_i = 1:3
    for dist_i = 1:3
        axh(fce_i,dist_i) = subplot(3,3,(fce_i-1)*3 + dist_i);
        dat = psudoI1_cell{fce_i,dist_i};
        dat_mean = mean(dat,2);
        dat_std = std(dat, [], 2);
        plot(peak_time, dat_mean', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        errorbar(peak_time, dat_mean', dat_std', 'lineWidth', 2, 'color', color_arr(4+dist_i,:));
        if (fce_i == 3)
            xlabel('pert time'); 
        end
        if (dist_i == 1)
            ylabel('I (kg)');
        end
        title(['fce' num2str(F_list(fce_i)) ' K' num2str(K_list(dist_i))]);
    end
end
sgtitle('Spring mass after release, dF/da');
linkaxes(axh);
ylim([-2 10]);
saveas(fh13, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure13.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig14:    The inertia measurement cross condition (using release epoc),
%           subj + spring;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sup-file1:    The pdf file generate for all subject measurments(K+D+I)
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
fh_sp1 = figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
mk_big = 10;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 1;

panels_size = [0.3168 0.1089]*5/6;
% panels_pos_itv = [0.1235 0.0958];
panels_pos_itv = [0.1235 0.0958] * 2/3;
% panels_pos_offset = [0.1461 0.1118]; 
% panels_pos_offset = [0.3461 0.3118]; 
panels_pos_offset = [(1 - (panels_size(1) * 2 + panels_pos_itv(1)*1))/2, ...
                     (1 - (panels_size(2) * 4 + panels_pos_itv(2)*3))/2]; 

for kdi = 1:3 % 1, stiffness, 2, damping, 3, inertia
for fce_i = 3%1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 2:8
            %             fh(pi,1) = figure(); hold on;

            fh_sp1 = figure('unit', 'inch', 'position', [0 0 8.5 11]); 
            clear axh;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot panels here
            for panel_i = 1:8
                panel_ri = ceil((8-panel_i + 1)/2);
                panel_ci = 2- mod(panel_i,2);
                panel_pos = [...
                    panels_pos_offset(1)+ (panel_ci-1)*(panels_pos_itv(1) + panels_size(1)), ...
                    ((panel_ri-1)*(panels_size(2)+panels_pos_itv(2))+panels_pos_offset(2)), ... % y-axis, row
                    panels_size(1), ...
                    panels_size(2)];
                axh(panel_i) = subplot('position',  panel_pos);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));   g
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(axh(1));
                plot(t_grids, fp_avg);
                axh(3) = subplot(axh(3));
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(axh(5));
                switch kdi
                    case 1
                        plot(t_grids, x_avg);
                    case 2
                        plot(t_grids, v_avg);
                    case 3
                        plot(t_grids, a_avg);
                end
                axh(7) = subplot(axh(7));
                switch kdi
                    case 1
                        plot(t_grids, f_avg ./ (x_avg - x_settled));
                        yline(psudo_stiffness1, 'linewidth', 2);
                    case 2
                        plot(t_grids,-f_avg ./ (v_avg));
                    case 3
                        plot(t_grids, f_avg ./ (a_avg));
                end
                %             ylim([0 2000]);
            
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                
                    subplot(axh(1)); hold on;
                    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(3)); hold on;
                    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                    %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(5)); hold on;
                    switch kdi 
                        case 1
                            x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                            plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                        case 2 %...........todo
                            v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                            plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                        case 3
                            a_shift = mean(diff(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)))./diff(celltmp1{ti,pi}.t(find(t>-0.1 & t<0))));
                            plot(t(2:end), diff(celltmp1{ti,pi}.v(2,:))./diff(celltmp1{ti,pi}.t), 'color', color_arr(4+dist_i,:));
                    end
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_avg;
                v_net = v_dat - v_avg;
                a_net = a_dat - a_avg; 
                fp_net = fp_dat - fp_avg;
                f_net = f_dat - f_avg;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_avg)./(x_dat - x_avg);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_avg)./(v_dat - v_avg);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                m_est = (f_dat - f_avg)./(a_dat - a_avg);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                
                    axh(2) = subplot(axh(2)); hold on; % subtracted Fp
                    plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(axh(4)); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_avg), 'color', color_arr(4+dist_i,:));
                    switch kdi
                        case 1
                            plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids(v_peak_idx),-f_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids(a_peak_idx),-f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    end
                    axh(6) = subplot(axh(6)); hold on;% subtracted F
                    switch kdi
                        case 1
                            plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                            plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                            plot(t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids, a_dat - a_avg, 'color', color_arr(4+dist_i,:));
                            plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    end
                    axh(8) = subplot(axh(8)); hold on;
                    switch kdi 
                        case 1
                            plot(t_grids, (f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids,-(f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(v_peak_idx), d_est_pt, 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids, (f_dat - f_avg)./(a_dat - a_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big);
                    end
                
            end
            
            
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(7:8), 'y');
                    yline(axh(8),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                
                switch kdi
                    case 1
                        ylim(axh(7), [0 1000]);
                    case 2
                        ylim(axh(7), [-20 100]);
                    case 3
                        ylim(axh(7), [-2 10]);
                end
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
                kdi_arr = 'KDI'; 
                sgtitle([kdi_arr(kdi) ': Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);

                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                switch kdi 
                    case 1
                        ylabel(axh(5), 'x');
                        ylabel(axh(7), 'df/dx');
                    case 2
                        ylabel(axh(5), 'v');
                        ylabel(axh(7), 'df/dv');
                    case 3
                        ylabel(axh(5), 'a');
                        ylabel(axh(7), 'df/da');
                end
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    switch kdi
                        case 1
                            ylabel(axh(6), 'dx');
                            ylabel(axh(8), 'df/dx');
                        case 2
                            ylabel(axh(6), 'dv');
                            ylabel(axh(8), 'df/dv');
                        case 3
                            ylabel(axh(6), 'da');
                            ylabel(axh(8), 'df/da');
                    end
                catch
                end
            
            if (kdi == 1 && dist_i == 1 && fce_i == 1 && pi == 2)
                %export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI.pdf', '-preserve_size');
            else
                %export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI.pdf', '-preserve_size', '-append');
            end
            clf;
        end
        
    end
    
end

end
% saveas(fh1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure1.png');

%% sup2, spring kdi
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
fh_sp1 = figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
mk_big = 10;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 1;

panels_size = [0.3168 0.1089]*5/6;
% panels_pos_itv = [0.1235 0.0958];
panels_pos_itv = [0.1235 0.0958] * 2/3;
% panels_pos_offset = [0.1461 0.1118]; 
% panels_pos_offset = [0.3461 0.3118]; 
panels_pos_offset = [(1 - (panels_size(1) * 2 + panels_pos_itv(1)*1))/2, ...
                     (1 - (panels_size(2) * 4 + panels_pos_itv(2)*3))/2]; 

for kdi = 1:3%1:3 % 1, stiffness, 2, damping, 3, inertia
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 2:6
            %             fh(pi,1) = figure(); hold on;

            fh_sp1 = figure('unit', 'inch', 'position', [0 0 8.5 11]); 
            clear axh;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot panels here
            for panel_i = 1:8
                panel_ri = ceil((8-panel_i + 1)/2);
                panel_ci = 2- mod(panel_i,2);
                panel_pos = [...
                    panels_pos_offset(1)+ (panel_ci-1)*(panels_pos_itv(1) + panels_size(1)), ...
                    ((panel_ri-1)*(panels_size(2)+panels_pos_itv(2))+panels_pos_offset(2)), ... % y-axis, row
                    panels_size(1), ...
                    panels_size(2)];
                axh(panel_i) = subplot('position',  panel_pos);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(axh(1));
                plot(t_grids, fp_avg);
                axh(3) = subplot(axh(3));
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(axh(5));
                switch kdi
                    case 1
                        plot(t_grids, x_avg);
                    case 2
                        plot(t_grids, v_avg);
                    case 3
                        plot(t_grids, a_avg);
                end
                axh(7) = subplot(axh(7));
                switch kdi
                    case 1
                        plot(t_grids, f_avg ./ (x_avg - x_settled));
                        yline(psudo_stiffness1, 'linewidth', 2);
                    case 2
                        plot(t_grids,-f_avg ./ (v_avg));
                    case 3
                        plot(t_grids, f_avg ./ (a_avg));
                end
                %             ylim([0 2000]);
            
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                
                    subplot(axh(1)); hold on;
                    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(3)); hold on;
                    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                    %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                    
                    subplot(axh(5)); hold on;
                    switch kdi 
                        case 1
                            x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                            plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                        case 2 %...........todo
                            v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                            plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                        case 3
                            a_shift = mean(diff(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)))./diff(celltmp1{ti,pi}.t(find(t>-0.1 & t<0))));
                            plot(t(2:end), diff(celltmp1{ti,pi}.v(2,:))./diff(celltmp1{ti,pi}.t), 'color', color_arr(4+dist_i,:));
                    end
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                x_net = x_dat - x_avg;
                v_net = v_dat - v_avg;
                a_net = a_dat - a_avg; 
                fp_net = fp_dat - fp_avg;
                f_net = f_dat - f_avg;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_avg)./(x_dat - x_avg);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_avg)./(v_dat - v_avg);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                m_est = (f_dat - f_avg)./(a_dat - a_avg);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                
                    axh(2) = subplot(axh(2)); hold on; % subtracted Fp
                    plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(axh(4)); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_avg), 'color', color_arr(4+dist_i,:));
                    switch kdi
                        case 1
                            plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids(v_peak_idx),-f_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids(a_peak_idx),-f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    end
                    axh(6) = subplot(axh(6)); hold on;% subtracted F
                    switch kdi
                        case 1
                            plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                            plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
                            plot(t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids, a_dat - a_avg, 'color', color_arr(4+dist_i,:));
                            plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    end
                    axh(8) = subplot(axh(8)); hold on;
                    switch kdi 
                        case 1
                            plot(t_grids, (f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids,-(f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(v_peak_idx), d_est_pt, 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids, (f_dat - f_avg)./(a_dat - a_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big);
                    end
                
            end
            
            
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(7:8), 'y');
                    yline(axh(8),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                
                switch kdi
                    case 1
                        ylim(axh(7), [0 1000]);
                    case 2
                        ylim(axh(7), [-20 100]);
                    case 3
                        ylim(axh(7), [-2 10]);
                end
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
                kdi_arr = 'KDI'; 
                sgtitle([kdi_arr(kdi) ': Force ' num2str(F_list(fce_i)) 'N K ' num2str(K_list(dist_i)) 'N/m pulse ' num2str(pi-1)]);

                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                switch kdi 
                    case 1
                        ylabel(axh(5), 'x');
                        ylabel(axh(7), 'df/dx');
                    case 2
                        ylabel(axh(5), 'v');
                        ylabel(axh(7), 'df/dv');
                    case 3
                        ylabel(axh(5), 'a');
                        ylabel(axh(7), 'df/da');
                end
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    switch kdi
                        case 1
                            ylabel(axh(6), 'dx');
                            ylabel(axh(8), 'df/dx');
                        case 2
                            ylabel(axh(6), 'dv');
                            ylabel(axh(8), 'df/dv');
                        case 3
                            ylabel(axh(6), 'da');
                            ylabel(axh(8), 'df/da');
                    end
                catch
                end
                
               set(fh_sp1, 'color', 'w');
            
            if (kdi == 1 && dist_i == 1 && fce_i == 1 && pi == 2)
                export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/spr_KDI.pdf', '-preserve_size');
            else
                export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/spr_KDI.pdf', '-preserve_size', '-append');
            end
            close all;
        end
        
    end
    
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sup-file3:    The pdf file generate for all subject measurments(K+D+I)
% have more panels, x, v, a, in all figures 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949_tmp.mat', 'data');
fh_sp1 = figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
mk_big = 10;
lw_ref = 0.1;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 1;

panels_size = [0.3168 0.1089]*5/6;
% panels_pos_itv = [0.1235 0.0958];
% panels_pos_itv = [0.1235 0.0958] * 2/3;
panels_pos_itv = [0.1235 0.0758] * 2/3;
% panels_pos_offset = [0.1461 0.1118]; 
% panels_pos_offset = [0.3461 0.3118]; 
panels_pos_offset = [(1 - (panels_size(1) * 2 + panels_pos_itv(1)*1))/2, ...
                     (1 - (panels_size(2) * 6 + panels_pos_itv(2)*5))/2]; 

for kdi = 1:3 % 1, stiffness, 2, damping, 3, inertia
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 2:8
            %             fh(pi,1) = figure(); hold on;

            fh_sp1 = figure('unit', 'inch', 'position', [0 0 8.5 11]); 
            clear axh;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot panels here
            for panel_i = 1:12
                panel_ri = ceil((12-panel_i + 1)/2);
                panel_ci = 2- mod(panel_i,2);
                panel_pos = [...
                    panels_pos_offset(1)+ (panel_ci-1)*(panels_pos_itv(1) + panels_size(1)), ...
                    ((panel_ri-1)*(panels_size(2)+panels_pos_itv(2))+panels_pos_offset(2)), ... % y-axis, row
                    panels_size(1), ...
                    panels_size(2)];
                axh(panel_i) = subplot('position',  panel_pos);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8,1),7,8); % 1 no -ert and 5 pert
            celltmp2 = reshape(data(1,fce_i,dist_i,:,1:8,2),7,8); % The unperturbed counterpart
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));  
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(axh(1));
                plot(t_grids, fp_avg);
                axh(3) = subplot(axh(3));
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(axh(5));
                plot(t_grids, x_avg);
                axh(7) = subplot(axh(7));
                plot(t_grids, v_avg);
                axh(9) = subplot(axh(9));
                plot(t_grids, a_avg);
                axh(11) = subplot(axh(11));
                switch kdi
                    case 1
                        plot(t_grids, f_avg ./ (x_avg - x_settled));
                        yline(psudo_stiffness1, 'linewidth', 2);
                    case 2
                        plot(t_grids,-f_avg ./ (v_avg));
                    case 3
                        plot(t_grids, f_avg ./ (a_avg));
                end
                %             ylim([0 2000]);
            
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                idx_release0= find(celltmp2{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                t0= celltmp2{ti,pi}.t - celltmp2{ti,pi}.t(idx_release0(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                idx_t0 = find(t0>=t_interest(1) & t0<=t_interest(2));
                %                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
                fp_dat0= interp1(t0(idx_t0), celltmp2{ti,pi}.Fp(2,idx_t0), t_grids, 'linear', 'extrap');    %
                x_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.x(2,idx_t0), t_grids, 'linear', 'extrap');     %
                v_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.v(2,idx_t0), t_grids, 'linear', 'extrap');     %
                f_dat0 =-interp1(t0(idx_t0), celltmp2{ti,pi}.f(2,idx_t0), t_grids, 'linear', 'extrap');     % ...
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                v_filter0 = filter(b,a,v_dat0);
                a_dat0 = [0 diff(v_filter0)./diff(t_grids)];
                
                
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                plot(t0, celltmp2{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                plot(t0, celltmp2{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                plot(t0, celltmp2{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                
                subplot(axh(7)); hold on;
                v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                plot(t0, celltmp2{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                
                subplot(axh(9)); hold on;
                a_shift = mean(diff(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)))./diff(celltmp1{ti,pi}.t(find(t>-0.1 & t<0))));
                a_shift0 = mean(diff(celltmp2{ti,pi}.v(2, find(t0>-0.1 & t0<0)))./diff(celltmp2{ti,pi}.t(find(t0>-0.1 & t0<0))));
                plot(t_grids, a_dat, 'color', color_arr(4+dist_i,:));
                plot(t_grids, a_dat0, 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                
                
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                fp_dat_tmp = fp_dat; fp_dat_tmp(fp_max_idx:end) = 0;
                [~,fp_stt_idx] = min(abs(fp_dat_tmp - fp_dat(fp_max_idx)*0.05));
                x_net = x_dat - x_dat0; %x_avg;
                v_net = v_dat - v_dat0; %v_avg;
                a_net = a_dat - a_dat0; %a_avg; 
                fp_net = fp_dat - fp_dat0; % fp_avg;
                f_net = f_dat - f_dat0; %f_avg;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
%                 k_est = (f_dat - f_avg)./(x_dat - x_avg);
                k_est = (f_dat - f_dat0)./(x_dat - x_dat0);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
%                 d_est =-(f_dat - f_avg)./(v_dat - v_avg);
                d_est =-(f_dat - f_dat0)./(v_dat - v_dat0);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
%                 m_est = (f_dat - f_avg)./(a_dat - a_avg);
                m_est = (f_dat - f_dat0)./(a_dat - a_dat0);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                
                    axh(2) = subplot(axh(2)); hold on; % subtracted Fp
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref); % marks the perturbation start 
                    plot(t_grids, fp_dat - fp_dat0, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(axh(4)); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_dat0), 'color', color_arr(4+dist_i,:));
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    switch kdi
                        case 1
                            plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            for axh_i = 2:12
                                xline(axh(axh_i),t_grids(x_net_idx), 'color', color_arr(ti,:), 'linewidth', lw_ref);
                            end
                        case 2
                            plot(t_grids(v_peak_idx),-f_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            for axh_i = 2:12
                                xline(axh(axh_i),t_grids(v_peak_idx), 'color', color_arr(ti,:), 'linewidth', lw_ref);
                            end
                        case 3
                            plot(t_grids(a_peak_idx),-f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            for axh_i = 2:12
                                xline(axh(axh_i),t_grids(a_peak_idx), 'color', color_arr(ti,:), 'linewidth', lw_ref);
                            end
                    end
                    axh(6) = subplot(axh(6)); hold on;% subtracted x
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, x_dat - x_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(axh(8)); hold on;% subtracted v
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, v_dat - v_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(10) = subplot(axh(10)); hold on;% subtracted a
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, a_dat - a_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    switch kdi 
                        case 1
%                             axh(6) = subplot(axh(6)); hold on;
                            plot(axh(3),t_grids(x_net_idx),-f_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(5),t_grids(x_net_idx), x_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(7),t_grids(x_net_idx), v_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(9),t_grids(x_net_idx), a_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(6),t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(8),t_grids(x_net_idx), v_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(10),t_grids(x_net_idx), a_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 2
%                             axh(8) = subplot(axh(8)); hold on;
                            plot(axh(3),t_grids(v_peak_idx),-f_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(5),t_grids(v_peak_idx), x_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(7),t_grids(v_peak_idx), v_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(9),t_grids(v_peak_idx), a_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(6),t_grids(v_peak_idx), x_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(8),t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(10),t_grids(v_peak_idx), a_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 3
%                             axh(10) = subplot(axh(10)); hold on;
                            plot(axh(3),t_grids(a_peak_idx),-f_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(5),t_grids(a_peak_idx), x_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(7),t_grids(a_peak_idx), v_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(9),t_grids(a_peak_idx), a_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(6),t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(8),t_grids(a_peak_idx), v_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(10),t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                    end
                    
                    axh(12) = subplot(axh(12)); hold on;
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    switch kdi 
                        case 1
                            plot(t_grids, (f_dat - f_dat0)./(x_dat - x_dat0), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 2
                            plot(t_grids,-(f_dat - f_dat0)./(v_dat - v_dat0), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(v_peak_idx), d_est_pt, 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 3
                            plot(t_grids, (f_dat - f_dat0)./(a_dat - a_dat0), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                    end
                
            end
            
            
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(11:12), 'y');
                    yline(axh(12),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                
                switch kdi
                    case 1
                        ylim(axh(11), [0 1000]);
                    case 2
                        ylim(axh(11), [-20 100]);
                    case 3
                        ylim(axh(11), [-2 10]);
                end
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
                kdi_arr = 'KDI'; 
                sgtitle([kdi_arr(kdi) ': Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);

                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'x');
                ylabel(axh(7), 'v');
                ylabel(axh(9), 'a');
                switch kdi 
                    case 1
                        ylabel(axh(11), 'df/dx');
                    case 2
                        ylabel(axh(11), 'df/dv');
                    case 3
                        ylabel(axh(11), 'df/da');
                end
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'dx');
                    ylabel(axh(8), 'dv');
                    ylabel(axh(10), 'da');
                    switch kdi
                        case 1
                            ylabel(axh(12), 'df/dx');
                        case 2
                            ylabel(axh(12), 'df/dv');
                        case 3
                            ylabel(axh(12), 'df/da');
                    end
                catch
                end
            
            yline(axh(7), 0, 'color', [0.5 0.5 0.5], 'linewidth', lw_ref);
            set(fh_sp1, 'Name', [kdi_arr(kdi) 'fce' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i))]);
            set(fh_sp1, 'color', 'w');
            if (kdi == 1 && dist_i == 1 && fce_i == 1 && pi == 2)
                export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_ff.pdf', '-preserve_size', ...
                    '-bookmark');
            else
                if pi == 2
                    export_fig('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_ff.pdf', ...
                         '-preserve_size', '-append', ...
                        '-bookmark');
                else
                    export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_ff.pdf', '-preserve_size', '-append');
                end
%                 export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_f.pdf', '-preserve_size', '-append');
            end
            close all; 
        end
        
    end
    
end

end
% saveas(fh1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure1.png');

%% supfig 3.5 generate pdf from subject, with pulses
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sup-file3:    The pdf file generate for all subject measurments(K+D+I)
% have more panels, x, v, a, in all figures 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
fh_sp1 = figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
mk_big = 10;
lw_ref = 0.1;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 1;

panels_size = [0.3168 0.1089]*5/6;
% panels_pos_itv = [0.1235 0.0958];
% panels_pos_itv = [0.1235 0.0958] * 2/3;
panels_pos_itv = [0.1235 0.0758] * 2/3;
% panels_pos_offset = [0.1461 0.1118]; 
% panels_pos_offset = [0.3461 0.3118]; 
panels_pos_offset = [(1 - (panels_size(1) * 2 + panels_pos_itv(1)*1))/2, ...
                     (1 - (panels_size(2) * 6 + panels_pos_itv(2)*5))/2]; 

for kdi = 1% :3 % 1, stiffness, 2, damping, 3, inertia
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 2:8
            %             fh(pi,1) = figure(); hold on;

            fh_sp1 = figure('unit', 'inch', 'position', [0 0 8.5 11]); 
            clear axh;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot panels here
            for panel_i = 1:12
                panel_ri = ceil((12-panel_i + 1)/2);
                panel_ci = 2- mod(panel_i,2);
                panel_pos = [...
                    panels_pos_offset(1)+ (panel_ci-1)*(panels_pos_itv(1) + panels_size(1)), ...
                    ((panel_ri-1)*(panels_size(2)+panels_pos_itv(2))+panels_pos_offset(2)), ... % y-axis, row
                    panels_size(1), ...
                    panels_size(2)];
                axh(panel_i) = subplot('position',  panel_pos);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:8),7,8); % 1 no -ert and 5 pert
%             celltmp2 = reshape(data(1,fce_i,dist_i,:,1:8,2),7,8); % The unperturbed counterpart
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));  
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(axh(1));
                plot(t_grids, fp_avg);
                axh(3) = subplot(axh(3));
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(axh(5));
                plot(t_grids, x_avg);
                axh(7) = subplot(axh(7));
                plot(t_grids, v_avg);
                axh(9) = subplot(axh(9));
                plot(t_grids, a_avg);
                axh(11) = subplot(axh(11));
                switch kdi
                    case 1
                        plot(t_grids, f_avg ./ (x_avg - x_settled));
                        yline(psudo_stiffness1, 'linewidth', 2);
                    case 2
                        plot(t_grids,-f_avg ./ (v_avg));
                    case 3
                        plot(t_grids, f_avg ./ (a_avg));
                end
                %             ylim([0 2000]);
            
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
%                 idx_release0= find(celltmp2{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
%                 t0= celltmp2{ti,pi}.t - celltmp2{ti,pi}.t(idx_release0(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 idx_t0 = find(t0>=t_interest(1) & t0<=t_interest(2));
                %                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 fp_dat0= interp1(t0(idx_t0), celltmp2{ti,pi}.Fp(2,idx_t0), t_grids, 'linear', 'extrap');    %
%                 x_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.x(2,idx_t0), t_grids, 'linear', 'extrap');     %
%                 v_dat0 = interp1(t0(idx_t0), celltmp2{ti,pi}.v(2,idx_t0), t_grids, 'linear', 'extrap');     %
%                 f_dat0 =-interp1(t0(idx_t0), celltmp2{ti,pi}.f(2,idx_t0), t_grids, 'linear', 'extrap');     % ...
                fp_dat0 = fp_avg;
                x_dat0 = x_avg;
                v_dat0 = v_avg;
                f_dat0 = f_avg;
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                v_filter0 = filter(b,a,v_dat0);
                a_dat0 = [0 diff(v_filter0)./diff(t_grids)];
                
                
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t0, celltmp2{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t0, celltmp2{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t0, celltmp2{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                
                subplot(axh(7)); hold on;
                v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
%                 plot(t0, celltmp2{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                
                subplot(axh(9)); hold on;
                a_shift = mean(diff(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)))./diff(celltmp1{ti,pi}.t(find(t>-0.1 & t<0))));
%                 a_shift0 = mean(diff(celltmp2{ti,pi}.v(2, find(t0>-0.1 & t0<0)))./diff(celltmp2{ti,pi}.t(find(t0>-0.1 & t0<0))));
                plot(t_grids, a_dat, 'color', color_arr(4+dist_i,:));
%                 plot(t_grids, a_dat0, 'color', color_arr(4+dist_i,:)/1.5, 'LineStyle', ':');
                
                
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                fp_dat_tmp = fp_dat; fp_dat_tmp(fp_max_idx:end) = 0;
                [~,fp_stt_idx] = min(abs(fp_dat_tmp - fp_dat(fp_max_idx)*0.05));
                x_net = x_dat - x_dat0; %x_avg;
                v_net = v_dat - v_dat0; %v_avg;
                a_net = a_dat - a_dat0; %a_avg; 
                fp_net = fp_dat - fp_dat0; % fp_avg;
                f_net = f_dat - f_dat0; %f_avg;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
%                 k_est = (f_dat - f_avg)./(x_dat - x_avg);
                k_est = (f_dat - f_dat0)./(x_dat - x_dat0);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
%                 d_est =-(f_dat - f_avg)./(v_dat - v_avg);
                d_est =-(f_dat - f_dat0)./(v_dat - v_dat0);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
%                 m_est = (f_dat - f_avg)./(a_dat - a_avg);
                m_est = (f_dat - f_dat0)./(a_dat - a_dat0);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                
                    axh(2) = subplot(axh(2)); hold on; % subtracted Fp
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref); % marks the perturbation start 
                    plot(t_grids, fp_dat - fp_dat0, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(axh(4)); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_dat0), 'color', color_arr(4+dist_i,:));
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    switch kdi
                        case 1
                            plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            for axh_i = 2:12
                                xline(axh(axh_i),t_grids(x_net_idx), 'color', color_arr(ti,:), 'linewidth', lw_ref);
                            end
                        case 2
                            plot(t_grids(v_peak_idx),-f_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            for axh_i = 2:12
                                xline(axh(axh_i),t_grids(v_peak_idx), 'color', color_arr(ti,:), 'linewidth', lw_ref);
                            end
                        case 3
                            plot(t_grids(a_peak_idx),-f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            for axh_i = 2:12
                                xline(axh(axh_i),t_grids(a_peak_idx), 'color', color_arr(ti,:), 'linewidth', lw_ref);
                            end
                    end
                    axh(6) = subplot(axh(6)); hold on;% subtracted x
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, x_dat - x_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(axh(8)); hold on;% subtracted v
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, v_dat - v_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(10) = subplot(axh(10)); hold on;% subtracted a
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, a_dat - a_dat0, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    switch kdi 
                        case 1
%                             axh(6) = subplot(axh(6)); hold on;
                            plot(axh(3),t_grids(x_net_idx),-f_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(5),t_grids(x_net_idx), x_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(7),t_grids(x_net_idx), v_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(9),t_grids(x_net_idx), a_dat(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(6),t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(8),t_grids(x_net_idx), v_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(10),t_grids(x_net_idx), a_net(x_net_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 2
%                             axh(8) = subplot(axh(8)); hold on;
                            plot(axh(3),t_grids(v_peak_idx),-f_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(5),t_grids(v_peak_idx), x_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(7),t_grids(v_peak_idx), v_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(9),t_grids(v_peak_idx), a_dat(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(6),t_grids(v_peak_idx), x_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(8),t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(10),t_grids(v_peak_idx), a_net(v_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 3
%                             axh(10) = subplot(axh(10)); hold on;
                            plot(axh(3),t_grids(a_peak_idx),-f_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(5),t_grids(a_peak_idx), x_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(7),t_grids(a_peak_idx), v_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(9),t_grids(a_peak_idx), a_dat(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(6),t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(8),t_grids(a_peak_idx), v_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                            plot(axh(10),t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                    end
                    
                    axh(12) = subplot(axh(12)); hold on;
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    switch kdi 
                        case 1
                            plot(t_grids, (f_dat - f_dat0)./(x_dat - x_dat0), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 2
                            plot(t_grids,-(f_dat - f_dat0)./(v_dat - v_dat0), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(v_peak_idx), d_est_pt, 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                        case 3
                            plot(t_grids, (f_dat - f_dat0)./(a_dat - a_dat0), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big, 'color', color_arr(ti,:));
                    end
                
            end
            
            
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(11:12), 'y');
                    yline(axh(12),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                
                switch kdi
                    case 1
                        ylim(axh(11), [0 1000]);
                    case 2
                        ylim(axh(11), [-20 100]);
                    case 3
                        ylim(axh(11), [-2 10]);
                end
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
                kdi_arr = 'KDI'; 
                sgtitle([kdi_arr(kdi) ': Force ' num2str(F_list(fce_i)) 'N dist ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);

                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'x');
                ylabel(axh(7), 'v');
                ylabel(axh(9), 'a');
                switch kdi 
                    case 1
                        ylabel(axh(11), 'df/dx');
                    case 2
                        ylabel(axh(11), 'df/dv');
                    case 3
                        ylabel(axh(11), 'df/da');
                end
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'dx');
                    ylabel(axh(8), 'dv');
                    ylabel(axh(10), 'da');
                    switch kdi
                        case 1
                            ylabel(axh(12), 'df/dx');
                        case 2
                            ylabel(axh(12), 'df/dv');
                        case 3
                            ylabel(axh(12), 'df/da');
                    end
                catch
                end
            
            yline(axh(7), 0, 'color', [0.5 0.5 0.5], 'linewidth', lw_ref);
            set(fh_sp1, 'Name', [kdi_arr(kdi) 'fce' num2str(F_list(fce_i)) 'dist' num2str(K_list(dist_i))]);
            set(fh_sp1, 'color', 'w');
            if (kdi == 1 && dist_i == 1 && fce_i == 1 && pi == 2)
                export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_K.pdf', '-preserve_size', ...
                    '-bookmark');
            else
                if pi == 2
                    export_fig('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_K.pdf', ...
                         '-preserve_size', '-append', ...
                        '-bookmark');
                else
                    export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_K.pdf', '-preserve_size', '-append');
                end
%                 export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_f.pdf', '-preserve_size', '-append');
            end
            close all; 
        end
        
    end
    
end

end
% saveas(fh1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/figure1.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sup-file4:    The pdf file generate for all subject measurments(K+D+I)
% have more panels, x, v, a, in all figures 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
lw_ref = 0.5;
fh_sp1 = figure(); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
mk_big = 10;
color_arr = colormap('lines');
close all;
t_interest = [-0.1 1.3]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
psudoD_cell = cell(3,3);
psudoI_cell = cell(3,3);
to_plot.fce_i = 1; 
to_plot.dist_i = 1;
to_plot.pi = 2;
to_plot.ifplot = 1;

panels_size = [0.3168 0.1089]*5/6;
% panels_pos_itv = [0.1235 0.0958];
% panels_pos_itv = [0.1235 0.0958] * 2/3;
panels_pos_itv = [0.1235 0.0758] * 2/3;
% panels_pos_offset = [0.1461 0.1118]; 
% panels_pos_offset = [0.3461 0.3118]; 
panels_pos_offset = [(1 - (panels_size(1) * 2 + panels_pos_itv(1)*1))/2, ...
                     (1 - (panels_size(2) * 6 + panels_pos_itv(2)*5))/2]; 

for kdi = 1%:3 % 1, stiffness, 2, damping, 3, inertia
for fce_i = 1:size(data,2)
    for dist_i = 1:size(data,3) % for each spring
        psudoK_mat = zeros(7,7);
        psudoD_mat = zeros(7,7);
        psudoI_mat = zeros(7,7);
        for pi = 2:6
            %             fh(pi,1) = figure(); hold on;

            fh_sp1 = figure('unit', 'inch', 'position', [0 0 8.5 11]); 
            clear axh;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot panels here
            for panel_i = 1:12
                panel_ri = ceil((12-panel_i + 1)/2);
                panel_ci = 2- mod(panel_i,2);
                panel_pos = [...
                    panels_pos_offset(1)+ (panel_ci-1)*(panels_pos_itv(1) + panels_size(1)), ...
                    ((panel_ri-1)*(panels_size(2)+panels_pos_itv(2))+panels_pos_offset(2)), ... % y-axis, row
                    panels_size(1), ...
                    panels_size(2)];
                axh(panel_i) = subplot('position',  panel_pos);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation: f, x, v, a
            x_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            a_avg = zeros(1, length(t_grids));  
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat =-interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap');
                % build up a filter to get a cleaner acceleration
                fc = 15;
                fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                % end of the filter
                a_dat = [0 diff(v_filter)./diff(t_grids)]; % as no a_dat here, just get the diffrerentiation of the v
                
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                v_avg = v_avg + v_dat;
                a_avg = a_avg + a_dat;
                f_avg = f_avg + f_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            a_avg = a_avg/cti;
            f_avg = f_avg/cti;
            
            x_settled = nanmean(x_avg(t_grids>1.0 & t_grids<1.1));
            % find the reference value
            psudo_stiffness = -f_avg ./ (x_avg - x_settled);
            psudo_stiffness0= mean(psudo_stiffness(t_grids<0));
            psudo_stiffness1= psudo_stiffness0/3.62; % after relase, the number is corresponding to 3kg hand
            
            
                % plot out the avg
                %             fh(pi,2) = figure();
                axh(1) = subplot(axh(1));
                plot(t_grids, fp_avg);
                axh(3) = subplot(axh(3));
                plot(t_grids,-f_avg);
                %             plot(t_grids, x_avg - x_avg(1));
                axh(5) = subplot(axh(5));
                plot(t_grids, x_avg);
                axh(7) = subplot(axh(7));
                plot(t_grids, v_avg);
                axh(9) = subplot(axh(9));
                plot(t_grids, a_avg);
                axh(11) = subplot(axh(11));
                switch kdi
                    case 1
                        plot(t_grids, f_avg ./ (x_avg - x_settled));
                        yline(psudo_stiffness1, 'linewidth', 2);
                    case 2
                        plot(t_grids,-f_avg ./ (v_avg));
                    case 3
                        plot(t_grids, f_avg ./ (a_avg));
                end
                %             ylim([0 2000]);
            
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
%                 length(idx_t);
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat =-interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:));
                %                 plot(t, celltmp1{ti,pi}.x(2,:) - x_shift, 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                x_shift = mean(celltmp1{ti,pi}.x(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(7)); hold on;
                v_shift = mean(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(9)); hold on;
                plot(t_grids, a_dat, 'color', color_arr(4+dist_i,:));
%                 a_shift = mean(diff(celltmp1{ti,pi}.v(2, find(t>-0.1 & t<0)))./diff(celltmp1{ti,pi}.t(find(t>-0.1 & t<0))));
%                 plot(t(2:end), diff(celltmp1{ti,pi}.v(2,:))./diff(celltmp1{ti,pi}.t), 'color', color_arr(4+dist_i,:));
                
                    
                %                 linkaxes(axh(1:3:5), 'x');
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                [~,fp_max_idx] = max(abs(fp_dat));
                fp_dat_tmp = fp_dat; fp_dat_tmp(fp_max_idx:end) = 0;
                [~,fp_stt_idx] = min(abs(fp_dat_tmp - fp_dat(fp_max_idx)*0.05));
                x_net = x_dat - x_avg;
                v_net = v_dat - v_avg;
                a_net = a_dat - a_avg; 
                fp_net = fp_dat - fp_avg;
                f_net = f_dat - f_avg;
                
                % get the K here
                x_net_tmp = x_net; x_net_tmp(1:fp_max_idx) = nan;
                [~, x_net_idx] = nanmin(x_net_tmp);
                k_est = (f_dat - f_avg)./(x_dat - x_avg);
                k_est_pt = k_est(x_net_idx);
                psudoK_mat(pi-1,ti) = k_est_pt;
                
                % get the D here
                v_net_tmp = v_net; v_net_tmp(1:fp_max_idx) = 0;
                [~, v_peak_idx] = max(v_net_tmp);    % only take after perturbation part
                d_est =-(f_dat - f_avg)./(v_dat - v_avg);
                d_est_pt = d_est(v_peak_idx);
                psudoD_mat(pi-1,ti) = d_est_pt;
                
                % get the I here
                a_net_tmp = a_net; a_net_tmp(1:fp_max_idx) = 0;
                [~, a_peak_idx] = min(a_net_tmp);    % only take after perturbation part
                m_est = (f_dat - f_avg)./(a_dat - a_avg);
                m_est_pt = m_est(a_peak_idx);
                psudoI_mat(pi-1,ti) = m_est_pt;
                
                
                % plot the subtraction in other panels
                
                    axh(2) = subplot(axh(2)); hold on; % subtracted Fp
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref); % marks the perturbation start 
                    plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:));
                    axh(4) = subplot(axh(4)); hold on;% subtracted x
                    plot(t_grids, -(f_dat - f_avg), 'color', color_arr(4+dist_i,:));
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    switch kdi
                        case 1
                            plot(t_grids(x_net_idx),-f_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids(v_peak_idx),-f_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids(a_peak_idx),-f_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    end
                    axh(6) = subplot(axh(6)); hold on;% subtracted x
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                    axh(8) = subplot(axh(8)); hold on;% subtracted v
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                    axh(10) = subplot(axh(10)); hold on;% subtracted a
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    plot(t_grids, a_dat - a_avg, 'color', color_arr(4+dist_i,:));
%                     plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    switch kdi 
                        case 1
                            axh(6) = subplot(axh(6)); hold on;
                            plot(t_grids(x_net_idx), x_net(x_net_idx), 'marker', '.', 'markersize', mk_big);
                        case 2
                            axh(8) = subplot(axh(8)); hold on;
                            plot(t_grids(v_peak_idx), v_net(v_peak_idx), 'marker', '.', 'markersize', mk_big);
                        case 3
                            axh(10) = subplot(axh(10)); hold on;
                            plot(t_grids(a_peak_idx), a_net(a_peak_idx), 'marker', '.', 'markersize', mk_big);
                    end
                    
                    axh(12) = subplot(axh(12)); hold on;
                    xline(t_grids(fp_stt_idx), 'linewidth', lw_ref);
                    switch kdi 
                        case 1
                            plot(t_grids, (f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(x_net_idx), k_est_pt, 'marker', '.', 'markersize', mk_big);
                        case 2
                            plot(t_grids,-(f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(v_peak_idx), d_est_pt, 'marker', '.', 'markersize', mk_big);
                        case 3
                            plot(t_grids, (f_dat - f_avg)./(a_dat - a_avg), 'color', color_arr(4+dist_i,:));
                            plot(t_grids(a_peak_idx), m_est_pt, 'marker', '.', 'markersize', mk_big);
                    end
                
            end
            
            
                try % if has axh8, plot, ifnot, noplot
                    linkaxes(axh(11:12), 'y');
                    yline(axh(12),psudo_stiffness1, 'linewidth', 2);
                catch
                end
                
                switch kdi
                    case 1
                        ylim(axh(11), [0 1000]);
                    case 2
                        ylim(axh(11), [-20 100]);
                    case 3
                        ylim(axh(11), [-2 10]);
                end
                linkaxes(axh, 'x');
                % plot notes here:
                xlim(axh(1), [-0.1 1.36]);
                
                kdi_arr = 'KDI'; 
                sgtitle([kdi_arr(kdi) ': Force ' num2str(F_list(fce_i)) 'N K ' num2str(K_list(dist_i)) 'cm pulse ' num2str(pi-1)]);

                title(axh(1), 'origin');
                try
                    title(axh(2), 'subtracted avg');
                catch
                end
                ylabel(axh(1), 'Fp');
                ylabel(axh(3), '-f');
                ylabel(axh(5), 'x');
                ylabel(axh(7), 'v');
                ylabel(axh(9), 'a');
                switch kdi 
                    case 1
                        ylabel(axh(11), 'df/dx');
                    case 2
                        ylabel(axh(11), 'df/dv');
                    case 3
                        ylabel(axh(11), 'df/da');
                end
                try
                    ylabel(axh(2), 'dFp');
                    ylabel(axh(4), '-df');
                    ylabel(axh(6), 'dx');
                    ylabel(axh(8), 'dv');
                    ylabel(axh(10), 'da');
                    switch kdi
                        case 1
                            ylabel(axh(12), 'df/dx');
                        case 2
                            ylabel(axh(12), 'df/dv');
                        case 3
                            ylabel(axh(12), 'df/da');
                    end
                catch
                end
            
            set(fh_sp1, 'Name', [kdi_arr(kdi) 'fce' num2str(F_list(fce_i)) 'K' num2str(K_list(dist_i))]);
            set(fh_sp1, 'color', 'w');
            if (kdi == 1 && dist_i == 1 && fce_i == 1 && pi == 2)
                export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/spr_KDI_f.pdf', '-preserve_size', ...
                    '-bookmark');
            else
                if pi == 2
                    export_fig('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/spr_KDI_f.pdf', ...
                         '-preserve_size', '-append', ...
                        '-bookmark');
                else
                    export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/spr_KDI_f.pdf', '-preserve_size', '-append');
                end
%                 export_fig(fh_sp1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_f.pdf', '-preserve_size', '-append');
            end
            close all; 
        end
        
    end
    
end

end

%% compare the perturbed after release and the perturbed before release 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3353_3417.mat', 'data');     %chenguang

figure(); 
F_list = [15, 20, 25];
K_list = [2.5, 5.0, 7.5];
color_arr = colormap('lines');
close all;
t_interest = [-0.5 0.55]; % s, calculate average from here 
t_static = [-0.5 0];
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_dat = zeros(3,3,15);
psudoD_dat = zeros(3,3,15);
psudoI_dat = zeros(3,3,15);
kdi_name = 'KDI';
for kdi = 1:3
for fce_i = 1:size(data,3)
    for dist_i = 1:size(data,4) % for each spring
        psudoK_mat = zeros(1,15); % 
%         for pi = 2:8%13%1:length(pertT_unq)
        fh(fce_i,dist_i) = figure('unit', 'inch', 'position', [0 0 5 11]); hold on;
%             fh(pi,1) = figure(); hold on;
            clear axh;
            
            celltmp1 = reshape(data(1,1,fce_i,dist_i,:,2),15,1); % 1 no -ert and 5 pert

            
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
%                 if isempty(celltmp1{ti,1})
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                pert_signal = celltmp1{ti,1}.Fp(2,:);
                pert_signal(abs(pert_signal)<max(abs(pert_signal))*0.05) = 0;
                idx_pert = find(pert_signal);
                
                
                
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_pert(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                idx_static= find(t_grids> t_static(1) & t_grids < t_static(2));
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); 
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...

                fc = 15;    fs = 500;
                [b,a] = butter(2,fc/(fs/2)); % 2nd-order, %? What is cut-off frequency?
                v_filter = filter(b,a,v_dat);
                a_dat = [0 diff(v_filter)./diff(t_grids)];
                
                % get the static error terms here 
                x_avg0 = mean(x_dat(idx_static));
                f_avg0 = mean(f_dat(idx_static)); 
                
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_pert(1));
                idx_tp = find(t>=t_interest(1) & t<=t_interest(2));
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    plot(t_grids, fp_dat_pert, 'g', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    plot(t_grids, x_dat_pert, 'g', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                    plot(t_grids, f_dat_pert, 'g', 'Marker', '.');
                end
                
%                 % also, plot out the origin
                axh(1) = subplot(6,1,1); hold on;
                plot(t_grids, fp_dat, 'color', color_arr(dist_i + 4,:));
                axh(2) = subplot(6,1,2); hold on;
                f_net = -(f_dat - f_avg0);
                plot(t_grids,f_net, 'color', color_arr(dist_i + 4,:));
                axh(3) = subplot(6,1,3); hold on;
                x_net = x_dat - x_avg0;
                plot(t_grids, x_net, 'color', color_arr(dist_i + 4,:));
                axh(4) = subplot(6,1,4); hold on;
                plot(t_grids, v_dat, 'color', color_arr(dist_i + 4,:));
                axh(5) = subplot(6,1,5); hold on;
                plot(t_grids, a_dat, 'color', color_arr(dist_i + 4,:));
                axh(6) = subplot(6,1,6); hold on;
                
                
                % 
                [~,idx_k] = min(x_dat - x_avg0);
                [~,idx_b] = max(v_dat);
                min_acc = min(a_dat(fp_dat < min(fp_dat)*0.05));
                [~,idx_i] = find(a_dat == min_acc);
                
                % get the "psudo-stiffness" values 
                psudo_stiffness =-(f_dat - f_avg0)./(x_dat - x_avg0); 
                psudo_damping = (f_dat - f_avg0)./(v_dat); 
                psudo_inertia = -(f_dat - f_avg0)./(a_dat); 
                
                subplot(axh(2)); 
                plot(t_grids(idx_k), f_net(idx_k), '*', 'color', color_arr(ti,:));
                plot(t_grids(idx_b), f_net(idx_b), 'o', 'color', color_arr(ti,:));
                plot(t_grids(idx_i), f_net(idx_i), 's', 'color', color_arr(ti,:));
                
                subplot(axh(3));
                plot(t_grids(idx_k), x_net(idx_k), '*', 'color', color_arr(ti,:));
                
                subplot(axh(4));
                plot(t_grids(idx_b), v_dat(idx_b), 'o', 'color', color_arr(ti,:));
                
                subplot(axh(5));
                plot(t_grids(idx_i), a_dat(idx_i), 's', 'color', color_arr(ti,:));
                
                subplot(axh(6));
                switch kdi 
                    case 1
                        plot(t_grids, -(f_dat - f_avg0)./(x_dat - x_avg0), 'color',  color_arr(dist_i + 4,:));
                        plot(t_grids(idx_k), psudo_stiffness(idx_k), '*', 'color', color_arr(ti,:));
                        ylabel(axh(6), 'df/dx');
                        ylim(axh(6), [-1000 1000]);
                    case 2
                        psudo_damping_nnan = psudo_damping;
%                         psudo_damping_nnan(abs(psudo_damping_nnan)>100) = nan;
                        plot(t_grids, psudo_damping, 'color',  color_arr(dist_i + 4,:));
                        plot(t_grids(idx_b), psudo_damping(idx_b), 'o', 'color', color_arr(ti,:));
                        ylabel(axh(6), 'df/dv');
                        ylim(axh(6), [-100 100]);
                    case 3
                        plot(t_grids, psudo_inertia, 'color',  color_arr(dist_i + 4,:));
                        plot(t_grids(idx_i), psudo_inertia(idx_i), 's', 'color', color_arr(ti,:));
                        ylabel(axh(6), 'df/da');
                        ylim(axh(6), [-10 10]);
                end
                
                psudoK_dat(fce_i, dist_i, ti) = psudo_stiffness(idx_k);
                psudoD_dat(fce_i, dist_i, ti) = psudo_damping(idx_b);
                psudoI_dat(fce_i, dist_i, ti) = psudo_inertia(idx_i);
                
                sgtitle([kdi_name(kdi) ' force' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
            end
            
            
            ylabel(axh(1), 'Fp');
            ylabel(axh(2), 'f');
            ylabel(axh(3), 'x');
            ylabel(axh(4), 'v');
            ylabel(axh(5), 'a');
            
            linkaxes(axh, 'x');
            xlim([-0.1 1]);
            
            set(gcf, 'name', [kdi_name(kdi) ' force' num2str(F_list(fce_i)) ' dist' num2str(K_list(dist_i))]);
            if (kdi==1 && fce_i == 1 && dist_i == 1)
                export_fig('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_beforeRelease.pdf', ...
                         '-preserve_size','-bookmark');
            else
                export_fig('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220322/subj_KDI_beforeRelease.pdf', ...
                         '-preserve_size', '-append', ...
                        '-bookmark');
            end
    end
    
end
end
