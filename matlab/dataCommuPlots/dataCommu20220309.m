% plots at the data communication at Mar-09-2022

% Topic: The pulse during movement. 
%           1. Can we estimate the spring stiffness during the movement?
%           2. Why the estimation method get different results?
%           3. Can we get consistant results from spring stiffness
%           estimation? If not, why? 


% Figures menu: 
% 1. The dF/dx during release and the dF/dx during hold 
% 2. The dF/dx during release is equal to the dF/dx during perturb
% 3. The t-dx figure along different robot nominal positions 
% 4. The t-dx figure along different energy exerting. 
% 5. (supplimentary), The dF/dx along on all combinations, tidy them in
% ppt.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1 The dF/dx during release and the dF/dx during hold 

% see what is the x and f change during the movement  
% plot: 3-row * 1 -col 
% 1,1: The original and perturbed censored force;
% 2,1: The original and perturbed position;
% 3,1: Stiffness dF/dx;                    
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
fh = figure('unit', 'inch', 'position', [0 0 5 6]); 
k_stfcoef = 13/20;
lw = 1; % pixels
F_list = [15, 20, 25];
K_list = [640, 320, 160];
color_arr = colormap('lines');
close all;
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
% psudoK_cell = cell(3,3);
fce_i = 3;
dist_i = 1;

pi = 1;
fh(pi,1) = figure(); hold on;

axh(1) = subplot(3,1,1); hold on;                     % plot PF
celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); 
idx_release = find(celltmp1{1,pi}.ts == 5);
t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
% plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));

% calculate the Unperturbed situation, x and f
x_avg = zeros(1, length(t_grids));
v_avg = zeros(1, length(t_grids));
f_avg = zeros(1, length(t_grids));
fp_avg= zeros(1, length(t_grids));
cti = 0; % count how many trials are added up
for ti = 1:1:size(celltmp1,1)
    if isempty(celltmp1{ti,1})
        continue;
    end
    cti = cti + 1;
    idx_release = find(celltmp1{ti,1}.ts == 5);
    t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
    idx_t = find(t>=t_interest(1) & t<=t_interest(2));
    length(idx_t)
    % intropolate (x, v, f, Fp) to t_grids
    
%     x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); 
    x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); 
    v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); 
    f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); 
    fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap');
    
    
    x_avg = x_avg + x_dat;
    f_avg = f_avg + f_dat;
    v_avg = v_avg + v_dat;
    fp_avg = fp_avg + fp_dat;
end
x_avg = x_avg/cti;
v_avg = v_avg/cti;
f_avg = f_avg/cti;  % this f_avg might not be right as it is 'centralized' after the release.


axh(1) = subplot(3,1,1);
plot(t_grids, f_avg, 'linewidth', lw);
%             plot(t_grids, v_avg);
axh(2) = subplot(3,1,2);
plot(t_grids, x_avg, 'linewidth', lw);

axh(3) = subplot(3,1,3);
k_avgest = f_avg ./ (x_avg - x_avg(end));
plot(t_grids, -k_avgest, 'linewidth', lw);


yline(K_list(dist_i), 'color', [0.5 0.5 0.5], 'linewidth', lw);
yline(K_list(dist_i)*k_stfcoef, 'color', [0.8 0.3 0.1], 'linewidth', lw);
legend('dF/dx', 'theoretical stiffness', 'portioned sitffness', 'fontSize', 12);
ylim(+[0 K_list(dist_i)*1.5]);

% plot notes here:
linkaxes(axh, 'x');
xlim(axh(1), [-0.1 1.36]);

% sgtitle(['Force ' num2str(F_list(fce_i)) ' dist ' num2str(K_list(dist_i)) 'pulse ' num2str(pi-1)]);
title(axh(1), 'Censored Force', 'fontSize', 15); 
ylabel(axh(1), 'F (N)', 'fontSize', 12);

title(axh(2), 'Handle position', 'fontSize', 15);
ylabel(axh(2), 'x (m)', 'fontSize', 12);

title(axh(3), 'stiffness: dF/{dx}', 'fontSize', 15 )
ylabel(axh(3), 'N/m', 'fontSize', 12);
xlabel(axh(3), 'time (s)', 'fontSize', 12);

ylabel(axh(2), 'ox (m)', 'fontSize', 12);
saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure1_ox.png');
% ylabel(axh(2), 'x (m)', 'fontSize', 12);
% saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure1_x.png');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. The dF/dx during release is equal to the dF/dx during perturb

% see what is the x and f change during the movement  
% plot: 4-row * 2 -col 
% 1,1: The original and perturbed command force;  1,2: The command force net effect 
% 2,1: The original and perturbed position;       2,2: The resulted position net effect;
% 3,1: The original and perturbed censored force; 3,2: The resulted censored net force;
% 4,1: None;                                      4.2: The psudo-stiffness: dF/dx
clear;
color_arr = colormap('lines');
close all; clc;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
fh = figure('unit', 'inch', 'position', [0 0 6 5]); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
k_stfcoef = 13/20;
lw = 1; % pixels
fs_big = 15;
fs_small = 12;
fs_mini = 8;
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
fce_i = 3;
dist_i = 1;
psudoK_mat = zeros(5,7); %
pi = 2;

axh(1) = subplot(3,1,1); hold on;                     % plot PF
celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
idx_release = find(celltmp1{1,pi}.ts == 5);
t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));

% calculate the Unperturbed situation, x and f
x_avg = zeros(1, length(t_grids));
v_avg = zeros(1, length(t_grids));
f_avg = zeros(1, length(t_grids));
fp_avg= zeros(1, length(t_grids));
cti = 0; % count how many trials are added up
for ti = 1:1:size(celltmp1,1)
    if isempty(celltmp1{ti,1})
        continue;
    end
    cti = cti + 1;
    idx_release = find(celltmp1{ti,1}.ts == 5);
    t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
    idx_t = find(t>=t_interest(1) & t<=t_interest(2));
    length(idx_t)
    % intropolate (x, f, Fp) to t_grids
    
    x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
%     x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); % check...
    v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
    f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
    fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
    
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
    f_avg = f_avg + f_dat;
    v_avg = v_avg + v_dat;
    fp_avg = fp_avg + fp_dat;
    %                 % also, plot out the origin
    %                 axh(1) = subplot(3,2,1); hold on;
    %                 plot(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), 'color', [0.5 0.5 0.5]);
    %                 axh(3) = subplot(3,2,3); hold on;
    %                 plot(t(idx_t), celltmp1{ti,1}.x(2,idx_t), 'color', [0.5 0.5 0.5]);
    %                 axh(5) = subplot(3,2,5); hold on;
    %                 plot(t(idx_t), celltmp1{ti,1}.f(2,idx_t), 'color', [0.5 0.5 0.5]);
end
x_avg = x_avg/cti;
v_avg = v_avg/cti;
f_avg = f_avg/cti;  % this f_avg might not be right as it is 'centralized' after the release.

% A TRICK TO UPDATE F_AVG HERE
if_favgupdate = 1;
if (if_favgupdate)
    f_avg_br_pert = 0;  % the force before released.
    t_tmp = 0;          % the number of trials
    for tti = 1:size(celltmp1,1)
        if isempty(celltmp1{tti,pi}) || pi == 1
            continue;
        end
        idx_release = find(celltmp1{tti,pi}.ts == 5);
        t = celltmp1{tti,pi}.t - celltmp1{tti,pi}.t(idx_release(1));
        idx_tmp = find(t>=-0.1 & t<0);
        f_tmp = mean(celltmp1{tti,pi}.f(2,idx_tmp));
        t_tmp = t_tmp + 1;
        f_avg_br_pert = f_avg_br_pert + f_tmp;
        celltmp1{tti,pi}.tshift = t;
    end
    f_avg_br_pert = f_avg_br_pert/t_tmp;
    f_avg_bfr = mean(f_avg(t_grids > -0.1 & t_grids<0));    % before release
    f_diff = f_avg_br_pert - f_avg_bfr;
    f_avg_upd = f_avg + f_diff;         % the before-release value are same now
    
    if (ifplot)
        clf;
        hold on;
        plot(t_grids, f_avg_upd, 'r.');
        plot(t_grids, f_avg, 'b.');
        for tti = 1:size(celltmp1,1)
            if isempty(celltmp1{tti,pi}) || pi == 1
                continue;
            end
            plot(celltmp1{tti,pi}.tshift, celltmp1{tti,pi}.f(2,:), 'Color', [0.5 0.5 0.5]);
        end
        
    end
    if (pi~=1)
        f_avg = f_avg_upd;
    end
end

% plot out the avg
%             fh(pi,2) = figure();
axh(1) = subplot(4,2,1);
plot(t_grids, fp_avg, 'linewidth', lw);
axh(3) = subplot(4,2,3);
plot(t_grids, f_avg, 'linewidth', lw);
%             plot(t_grids, v_avg);
axh(5) = subplot(4,2,5);
plot(t_grids, x_avg, 'linewidth', lw);

axh(7) = subplot(4,2,7);
k_avgest = f_avg ./ (x_avg - x_avg(end));
plot(t_grids, -k_avgest, 'linewidth', lw);
% plot the perturbed one, -perturbed
for ti = 1:size(celltmp1,1)
    if isempty(celltmp1{ti,pi}) || pi == 1
        continue;
    end
    
    idx_release = find(celltmp1{ti,pi}.ts == 5);
    t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
    subplot(axh(1)); hold on;
    plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
     
    subplot(axh(3)); hold on;
    plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
    %                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
    
    subplot(axh(5)); hold on;
    plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
%     plot(t, celltmp1{ti,pi}.ox(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
    idx_t = find(t>=t_interest(1) & t<=t_interest(2));
    length(idx_t)
    
    fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
    x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
%     x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     %
    v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
    f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
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
    
    % plot the subtraction in other panels
    axh(2) = subplot(4,2,2); hold on; % subtracted Fp
    plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:), 'linewidth', lw);
    axh(4) = subplot(4,2,4); hold on;% subtracted x
    plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:), 'linewidth', lw);
    %                 plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
    axh(6) = subplot(4,2,6); hold on;% subtracted F
    plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:), 'linewidth', lw);
    axh(8) = subplot(4,2,8); hold on; % subtracted dF/dx
    plot(t_grids, -(f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
    %                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
    
    
end
linkaxes(axh, 'x');
linkaxes(axh(3:4), 'y');
% linkaxes(axh(5:6), 'y');

% plot notes here:
xlim(axh(1), [-0.1 1.36]);
sgtitle('stiffness measurement comparision');
title(axh(1), 'unperturbed and perturbed', 'fontsize', fs_big);
title(axh(2), 'subtracted', 'fontsize', fs_big);

ylabel(axh(1), 'Fp (N)', 'fontsize', fs_small);
ylabel(axh(3), 'F (N)', 'fontsize', fs_small);
ylabel(axh(5), 'x (m)', 'fontsize', fs_small);
ylabel(axh(7), 'K (N/m)', 'fontsize', fs_small);
xlabel(axh(7), 'time (s)',  'fontsize', fs_small);

ylabel(axh(2), 'dFp (N)', 'fontsize', fs_small);
ylabel(axh(4), 'dF (N)', 'fontsize', fs_small);
ylabel(axh(6), 'dx (m)', 'fontsize', fs_small);
ylabel(axh(8), 'K (N/m)', 'fontsize', fs_small);
xlabel(axh(8), 'time (s)',  'fontsize', fs_small);

subplot(axh(7));
title('dF/dx in release', 'fontsize', fs_small);
ylim(+[0 K_list(dist_i)*1.5]);
yline(K_list(dist_i), 'color', [0.5 0.5 0.5], 'linewidth', lw);
yline(K_list(dist_i)*k_stfcoef, 'color', [0.8 0.3 0.1], 'linewidth', lw);
legend('dF/dx', 'theoretical stiffness', 'portioned sitffness', 'fontSize', fs_mini);

% ylabel(axh(5), 'ox (m)', 'fontsize', fs_small);
% ylabel(axh(6), 'dox (m)', 'fontsize', fs_small);
% legend(axh(7),'dF/dox', 'theoretical stiffness', 'portioned sitffness', 'fontSize', fs_mini);

subplot(axh(8));
title('dF/dx in perturb', 'fontsize', fs_small);
ylim(+[0 K_list(dist_i)*1.5]);
yline(K_list(dist_i), 'color', [0.5 0.5 0.5], 'linewidth', lw);
yline(K_list(dist_i)*k_stfcoef, 'color', [0.8 0.3 0.1], 'linewidth', lw);

saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2_x.png');
% saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2_ox.png');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. The t-dx figure along different robot nominal positions 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sanity check: why the position difference is so huge in some of the pulses

% Compare the displacement difference at different Velocity
% Also, compare the condition where there is no release (left) and with
% release (right) 
% Arrange the displacement difference when there is release on the velocity
% dot product with the perturb force 
% plot: 1-row * 2 -col 
% col1: without release: 
%           x-pert time, 
%           y-velocity before perturb, 
%           z-displacement different in release

% The left column, only plot the displacement (raw) 
clear; 
color_arr = colormap('lines');
close all;
fs_big = 15;
fs_small = 12;
fs_mini = 8;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4000_4002.mat', 'data');    % I don't have optotrak information in these sessions
fh1 = figure('unit', 'inch', 'position', [0 0 6 2.5]); 
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1 %1:size(data,2)
    for dist_i = 1 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t);

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
%                 x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1);
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                v_offset = v_dat(offset_idx);
                
                fp_effectv = fp_dat .* v_dat;
                v_offset = sum(fp_effectv);
                
                ifplot = 0; % check the intropolate 
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                axh(1) = subplot(1,2,1); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_dat)); % get position offset from the raw data
                vel_offset = v_offset * ones(size(v_dat));
                plot3(t_grids - t_offset, pos_offset ,x_dat - pos_offset, 'color', color_arr(4+2,:));

            end
        end

    end
    
end
subplot(axh(1));
view([-60, 60]); 
zlim([-0.04 0.04]);
%ylim([-100 400]);
xlabel('t (s)', 'fontsize', fs_small);
ylabel('x_0 (m)', 'fontsize', fs_small);   % sum(f_offset.*fp)
zlabel('x - x_0 (m)', 'fontsize', fs_small);
title('perturb no release', 'fontsize', fs_big);

% The right column, only plot displacement difference 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
F_list = [15, 20, 25];
K_list = [640, 320, 160];
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 2 %1:size(data,2)
    for dist_i = 2 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t);
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
%                 x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
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
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t);

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
%                 x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1) + 150; % 100ms after perturb~=0
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                v_offset = v_dat(offset_idx); % change this into the convolution between velocity and perturb force. 
                
                fp_effectv = fp_dat .* v_dat;
                v_offset = sum(fp_effectv);
                
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
                   

                axh(2) = subplot(1,2,2); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
                vel_offset = v_offset * ones(size(v_avg));
                
%                 plot3(t_grids - t_offset, vel_offset ,x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                plot3(t_grids - t_offset, pos_offset ,x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
%                 
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            
            % plot the subtracted position and force, in another
            % figure/panel
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
    end
    
end
subplot(axh(2));
view([-60, 60]); 
% ylim([-100 400]);
xlabel('t (s)', 'fontsize', fs_small);
ylabel('x_0 (m)', 'fontsize', fs_small);   % sum(f_offset.*fp)
zlabel('dx (m)', 'fontsize', fs_small);
zlim([-0.04 0.04]);
title('perturb during release', 'fontsize', fs_big);

linkaxes(axh, 'xy');
sgtitle('different robot configuration generate same perturbation')
saveas(fh1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure3.png')


% subplot(axh(1));
% view([0, 0]);  zlim([-0.04 0.04]);
% subplot(axh(2));
% view([0, 0]); zlim([-0.04 0.04]);
% saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossForceEffectVelocity_overlay.png')


%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sanity check: why the position difference is so huge in some of the pulses

% Compare the displacement difference at different Velocity
% Also, compare the condition where there is no release (left) and with
% release (right) 
% Arrange the displacement difference when there is release on the velocity
% dot product with the perturb force 
% plot: 1-row * 2 -col 
% col1: without release: 
%           x-pert time, 
%           y-velocity before perturb, 
%           z-displacement different in release

% The left column, only plot the displacement (raw) fs_big = 15;
clear
color_arr = colormap('lines');
close all;
fs_big = 15;
fs_small = 12;
fs_mini = 8;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4000_4002.mat', 'data');
fh1 = figure('unit', 'inch', 'position', [0 0 6 2.5]); 
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 1 %1:size(data,2)
    for dist_i = 1 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));

            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
                ifplot = 0; % controls whether plot or not
                if (ifplot)
                    clf;
                    subplot(3,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,3}.Fp(2,idx_t), 'b');
                    plot(t_grids, fp_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,2);  hold on;
                    plot(t(idx_t), celltmp1{ti,3}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(3,1,3); hold on;
                    plot(t(idx_t), celltmp1{ti,3}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
 
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1);
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                v_offset = v_dat(offset_idx);
                
                fp_effectv = fp_dat .* v_dat;
                v_offset = sum(fp_effectv);
                
                ifplot = 1;
                if (ifplot)
                    clf;
                    subplot(2,1,1);  hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'b');
                    plot(t_grids, x_dat, 'r', 'Marker', '.');
                    
                    subplot(2,1,2); hold on;
                    plot(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), 'b');
                    plot(t_grids, f_dat, 'r', 'Marker', '.');
                end
                   
                axh(1) = subplot(1,2,1); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
                vel_offset = v_offset * ones(size(v_avg));
                plot3(t_grids - t_offset, vel_offset ,x_dat - pos_offset, 'color', color_arr(4+2,:));
                

                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 

                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  

%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end

        

        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
%         close all;
    end
    
end
subplot(axh(1));
view([-60, 60]); 
zlim([-0.04 0.04]);
ylim([-100 400]);
xlabel('t (s)', 'fontsize', fs_small);
ylabel('\SigmaF*v (Nm/s)', 'fontsize', fs_small);   % sum(f_offset.*fp)
zlabel('x - x_0 (m)', 'fontsize', fs_small);
title('perturb no release', 'fontsize', fs_big);

% The right column, only plot displacement difference 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
F_list = [15, 20, 25];
K_list = [640, 320, 160];
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for fce_i = 2 %1:size(data,2)
    for dist_i = 2 %1:size(data,3) % for each spring
        psudoK_mat = zeros(5,7); % 
        for pi = 1:6 %13%1:length(pertT_unq)
             figure(fh1); hold on;
            
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
               
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
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;

            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;
            
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)

                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    % 
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     % 
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     % 
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
%                 linkaxes(axh(1:3:5), 'x');
                
                % offset index
                offset_idx = find(fp_dat ~= 0);
                offset_idx = offset_idx(1) + 150; % 100ms after perturb~=0
                % offset time 
                t_offset = t_grids(offset_idx);
                % offset position
                x_offset = x_dat(offset_idx);
                v_offset = v_dat(offset_idx); % change this into the convolution between velocity and perturb force. 
                
                fp_effectv = fp_dat .* v_dat;
                v_offset = sum(fp_effectv);
                
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
                   

                axh(2) = subplot(1,2,2); hold on;% subtracted x
                pos_offset = x_offset * ones(size(x_avg)); % get position offset from the raw data
                vel_offset = v_offset * ones(size(v_avg));
                
                plot3(t_grids - t_offset, vel_offset ,x_dat - x_avg, 'color', color_arr(4+dist_i,:));
                
                % There will be a force and position peak at 0~0.4s after
                % the start of the perturbation 
                pert_idx = find(abs(fp_dat) > 0.5); 
                pert_idx = pert_idx(1): min((pert_idx(1) + 0.4*freq), length(t_grids)); 
%                 
                % Find and plot on the original data point 
                pert0idx  = find(abs(fp_dat) > 0.01);
                pert0idx = pert0idx(1);  
%                 
                x0 = x_dat(pert0idx) - x_avg(pert0idx);
                f0 = f_dat(pert0idx) - f_avg(pert0idx);
                
                % Find and plot on the peak data point 
                [pksx,locsx]=findpeaks(-(x_dat(pert_idx) - x_avg(pert_idx)), 'MinPeakDistance', 100);
                [pksf,locsf]=findpeaks((f_dat(pert_idx) - f_avg(pert_idx)), 'MinPeakDistance', 100);
                
                x1 = x_dat(pert_idx(locsx(1))) - x_avg(pert_idx(locsx(1)));
                f1 = f_dat(pert_idx(locsf(1))) - f_avg(pert_idx(locsf(1)));
                
                psudoK = (f1-f0)/(x1-x0);
                psudoK_mat(ti,pi-1) = psudoK;
            end
            
            % plot the subtracted position and force, in another
            % figure/panel
        end
        psudoK_cell{fce_i,dist_i} = -psudoK_mat;
    end
    
end
subplot(axh(2));
view([0, 80]); 
view([-60, 60]); 
zlim([-0.04 0.04]);
ylim([-100 400]);
xlabel('t (s)', 'FontSize', fs_small);
ylabel('\SigmaF*v (Nm/s)', 'FontSize', fs_small);   % sum(f_offset.*fp)
zlabel('dx (m)', 'FontSize', fs_small);
title('perturb during release', 'FontSize', fs_big);
sgtitle('different input energy generate different perturbation');
saveas(fh1, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure4.png')

% subplot(axh(1));
% view([0, 0]);  zlim([-0.04 0.04]);
% subplot(axh(2));
% view([0, 0]); zlim([-0.04 0.04]);
% saveas(fh1, 'sanityCheck_StiffnessMeasurement/pulsePertDuringMovement/sanityCheck_pulseAcrossForceEffectVelocity_overlay.png')

%% save to the extra-saving figures...
% see what is the x and f change during the movement  
% plot: 4-row * 2 -col 
% 1,1: The original and perturbed command force;  1,2: The command force net effect 
% 2,1: The original and perturbed position;       2,2: The resulted position net effect;
% 3,1: The original and perturbed censored force; 3,2: The resulted censored net force;
% 4,1: None;                                      4.2: The psudo-stiffness: dF/dx
clear;
color_arr = colormap('lines');
close all; clc;

fh = figure('unit', 'inch', 'position', [0 0 8.5 11]); 
panels_size = [0.3168 0.1089]*5/6;
% panels_pos_itv = [0.1235 0.0958];
panels_pos_itv = [0.1235 0.0958] * 2/3;
% panels_pos_offset = [0.1461 0.1118]; 
% panels_pos_offset = [0.3461 0.3118]; 
panels_pos_offset = [(1 - (panels_size(1) * 2 + panels_pos_itv(1)*1))/2, ...
                     (1 - (panels_size(2) * 4 + panels_pos_itv(2)*3))/2]; 
% for panel_i = 1:8
% panel_ri = ceil((8-panel_i + 1)/2);
% panel_ci = 2- mod(panel_i,2);
% panel_pos = [...
%     panels_pos_offset(1)+ (panel_ci-1)*(panels_pos_itv(1) + panels_size(1)), ...
%     ((panel_ri-1)*(panels_size(2)+panels_pos_itv(2))+panels_pos_offset(2)), ... % y-axis, row
%     panels_size(1), ...
%     panels_size(2)];
% axh(panel_i) = subplot('position',  panel_pos);
% end



load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
% fh = figure('unit', 'inch', 'position', [0 0 8.5 11]); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
k_stfcoef = 13/20;
lw = 1; % pixels
fs_big = 15;
fs_small = 12;
fs_mini = 8;
t_interest = [-0.1 2]; % s, calculate average from here 
freq = 500; 
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
for dist_i = 1:3
    for fce_i = 1:3
        psudoK_mat = zeros(5,7); %
        for pi = 2:6 % only 6 perturbations
            
%             axh(1) = subplot(3,1,1); hold on;                     % plot PF
            celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
            idx_release = find(celltmp1{1,pi}.ts == 5);
            t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
            
            % calculate the Unperturbed situation, x and f
            x_avg = zeros(1, length(t_grids));
            v_avg = zeros(1, length(t_grids));
            f_avg = zeros(1, length(t_grids));
            fp_avg= zeros(1, length(t_grids));
            cti = 0; % count how many trials are added up
            for ti = 1:1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                cti = cti + 1;
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                % intropolate (x, f, Fp) to t_grids
                
                %     x_dat = interp1(t(idx_t), celltmp1{ti,1}.x(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                x_dat = interp1(t(idx_t), celltmp1{ti,1}.ox(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                v_dat = interp1(t(idx_t), celltmp1{ti,1}.v(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                f_dat = interp1(t(idx_t), celltmp1{ti,1}.f(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                fp_dat= interp1(t(idx_t), celltmp1{ti,1}.Fp(2,idx_t), t_grids, 'linear', 'extrap'); % check...
                
                x_avg = x_avg + x_dat;
                f_avg = f_avg + f_dat;
                v_avg = v_avg + v_dat;
                fp_avg = fp_avg + fp_dat;
            end
            x_avg = x_avg/cti;
            v_avg = v_avg/cti;
            f_avg = f_avg/cti;  % this f_avg might not be right as it is 'centralized' after the release.
            
            % A TRICK TO UPDATE F_AVG HERE
            if_favgupdate = 1;
            if (if_favgupdate)
                f_avg_br_pert = 0;  % the force before released.
                t_tmp = 0;          % the number of trials
                for tti = 1:size(celltmp1,1)
                    if isempty(celltmp1{tti,pi}) || pi == 1
                        continue;
                    end
                    idx_release = find(celltmp1{tti,pi}.ts == 5);
                    t = celltmp1{tti,pi}.t - celltmp1{tti,pi}.t(idx_release(1));
                    idx_tmp = find(t>=-0.1 & t<0);
                    f_tmp = mean(celltmp1{tti,pi}.f(2,idx_tmp));
                    t_tmp = t_tmp + 1;
                    f_avg_br_pert = f_avg_br_pert + f_tmp;
                    celltmp1{tti,pi}.tshift = t;
                end
                f_avg_br_pert = f_avg_br_pert/t_tmp;
                f_avg_bfr = mean(f_avg(t_grids > -0.1 & t_grids<0));    % before release
                f_diff = f_avg_br_pert - f_avg_bfr;
                f_avg_upd = f_avg + f_diff;         % the before-release value are same now
                
                if (0)
                    clf;
                    hold on;
                    plot(t_grids, f_avg_upd, 'r.');
                    plot(t_grids, f_avg, 'b.');
                    for tti = 1:size(celltmp1,1)
                        if isempty(celltmp1{tti,pi}) || pi == 1
                            continue;
                        end
                        plot(celltmp1{tti,pi}.tshift, celltmp1{tti,pi}.f(2,:), 'Color', [0.5 0.5 0.5]);
                    end
                    
                end
                if (pi~=1)
                    f_avg = f_avg_upd;
                end
            end
            
            
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
            
            % plot out the avg
            %             fh(pi,2) = figure();
%             axh(1) = subplot(4,2,1);
            axh(1) = subplot(axh(1));
            plot(t_grids, fp_avg, 'linewidth', lw);
%             axh(3) = subplot(4,2,3);
            axh(3) = subplot(axh(3));
            plot(t_grids, f_avg, 'linewidth', lw);
            %             plot(t_grids, v_avg);
%             axh(5) = subplot(4,2,5);
            axh(5) = subplot(axh(5));
            plot(t_grids, x_avg, 'linewidth', lw);
            
%             axh(7) = subplot(4,2,7);
            axh(7) = subplot(axh(7));
            k_avgest = f_avg ./ (x_avg - x_avg(end));
            plot(t_grids, -k_avgest, 'linewidth', lw);
            % plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                subplot(axh(1)); hold on;
                plot(t, celltmp1{ti,pi}.Fp(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
                
                subplot(axh(3)); hold on;
                plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
                %                 plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(4+dist_i,:));
                
                subplot(axh(5)); hold on;
                %     plot(t, celltmp1{ti,pi}.x(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
                plot(t, celltmp1{ti,pi}.ox(2,:), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
                idx_t = find(t>=t_interest(1) & t<=t_interest(2));
                length(idx_t)
                
                fp_dat= interp1(t(idx_t), celltmp1{ti,pi}.Fp(2,idx_t), t_grids, 'linear', 'extrap');    %
                %     x_dat = interp1(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), t_grids, 'linear', 'extrap');     %
                x_dat = interp1(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), t_grids, 'linear', 'extrap');     %
                v_dat = interp1(t(idx_t), celltmp1{ti,pi}.v(2,idx_t), t_grids, 'linear', 'extrap');     %
                f_dat = interp1(t(idx_t), celltmp1{ti,pi}.f(2,idx_t), t_grids, 'linear', 'extrap');     % ...
                %                 linkaxes(axh(1:3:5), 'x');
                
                total_eng = sum(fp_dat.*v_dat); 
                if total_eng > 0 
                    title_str_appd = '+ energy';
                else
                    title_str_appd = '- energy';
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
                
                % plot the subtraction in other panels
%                 axh(2) = subplot(4,2,2); 
                axh(2) = subplot(axh(2));
                hold on; % subtracted Fp
                plot(t_grids, fp_dat - fp_avg, 'color', color_arr(4+dist_i,:), 'linewidth', lw);
%                 axh(4) = subplot(4,2,4); hold on;% subtracted x
                axh(4) = subplot(axh(4));
                hold on;
                plot(t_grids, f_dat - f_avg, 'color', color_arr(4+dist_i,:), 'linewidth', lw);
                %                 plot(t_grids, v_dat - v_avg, 'color', color_arr(4+dist_i,:));
%                 axh(6) = subplot(4,2,6); hold on;% subtracted F
                axh(6) = subplot(axh(6));
                hold on;
                plot(t_grids, x_dat - x_avg, 'color', color_arr(4+dist_i,:), 'linewidth', lw);
%                 axh(8) = subplot(4,2,8); hold on; % subtracted dF/dx
                axh(8) = subplot(axh(8)); 
                hold on; 
                plot(t_grids, -(f_dat - f_avg)./(x_dat - x_avg), 'color', color_arr(4+dist_i,:), 'linewidth', lw);
                %                 plot(t_grids, (f_dat - f_avg)./(v_dat - v_avg), 'color', color_arr(4+dist_i,:));
                
                
            end
            linkaxes(axh, 'x');
            linkaxes(axh(3:4), 'y');
            % linkaxes(axh(5:6), 'y');
            
            % plot notes here:
            xlim(axh(1), [-0.1 1.36]);
            % sgtitle('stiffness measurement comparision');
            title_str = ['Force ' num2str(F_list(fce_i))  'N stiffness' num2str(K_list(dist_i)) ...
                'N/m pert' num2str(pi-1)];
            sgtitle([title_str ', ' title_str_appd]);
%             title(title_str, 'position', [0.25 0.88], 'fontSize', fs_big + 2);
            title(axh(1), 'unperturbed and perturbed', 'fontsize', fs_big);
            title(axh(2), 'subtracted', 'fontsize', fs_big);
            
            ylabel(axh(1), 'Fp (N)', 'fontsize', fs_small);
            ylabel(axh(3), 'F (N)', 'fontsize', fs_small);
            ylabel(axh(5), 'x (m)', 'fontsize', fs_small);
            ylabel(axh(7), 'K (N/m)', 'fontsize', fs_small);
            xlabel(axh(7), 'time (s)',  'fontsize', fs_small);
            
            ylabel(axh(2), 'dFp (N)', 'fontsize', fs_small);
            ylabel(axh(4), 'dF (N)', 'fontsize', fs_small);
            ylabel(axh(6), 'dx (m)', 'fontsize', fs_small);
            ylabel(axh(8), 'K (N/m)', 'fontsize', fs_small);
            xlabel(axh(8), 'time (s)',  'fontsize', fs_small);
            
            
            subplot(axh(7));
            title('dF/dx in release', 'fontsize', fs_small);
            ylim(+[0 K_list(dist_i)*1.5]);
            yline(K_list(dist_i), 'color', [0.5 0.5 0.5], 'linewidth', lw);
            yline(K_list(dist_i)*k_stfcoef, 'color', [0.8 0.3 0.1], 'linewidth', lw);
            % legend('dF/dx', 'theoretical stiffness', 'portioned sitffness', 'fontSize', fs_mini);
            
            subplot(axh(8));
            title('dF/dx in perturb', 'fontsize', fs_small);
            ylim(+[0 K_list(dist_i)*1.5]);
            yline(K_list(dist_i), 'color', [0.5 0.5 0.5], 'linewidth', lw);
            yline(K_list(dist_i)*k_stfcoef, 'color', [0.8 0.3 0.1], 'linewidth', lw);
            
            % saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.png');
            set(fh, 'color', 'w');
%             set(fh, 'position', [1.25 3 6 5]);
%             set(fh, 'PaperSize', [20 10]);
            if (dist_i == 1 && fce_i == 1 && pi == 2)
%                 export_fig('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.pdf', '-nocrop');
%                 saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.pdf');
%                 print(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.pdf', '-dpdf', '-r0');
                export_fig(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.pdf', '-preserve_size');
            else
%                   print(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.pdf', '-dpdf', '-r0', '-append');
                export_fig(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.pdf', '-preserve_size', '-append');
%                 saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure2.pdf', '-append');
            end
            clf;
        end
    end
end

%% Show the similarity between the optotrak recording and the WAM recording
% As the recording have the difference only when it has the biggest force
% exerted. Here I plan to plot the biggst force, with the highest stiffness
% of the springs. 



fh = figure(); 

clear;
color_arr = colormap('lines');
close all; clc;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3987_3999.mat', 'data');
fh = figure('unit', 'inch', 'position', [0 0 6 3]); 
F_list = [15, 20, 25];
K_list = [640, 320, 160];
k_stfcoef = 13/20;
lw = 1; % pixels
fs_big = 15;
fs_small = 12;
fs_mini = 8;
t_interest = [-0.1 2]; % s, calculate average from here
freq = 500;
t_grids = t_interest(1) : 1/freq : t_interest(2);
psudoK_cell = cell(3,3);
fce_i = 3;
dist_i = 1;
psudoK_mat = zeros(5,7); %
pi = 1;
for panel_i = 1:2
    switch panel_i 
        case 1
            fce_i = 1; 
            dist_i = 3;
        case 2
            fce_i = 3;
            dist_i = 1;
    end
    
axh(panel_i) = subplot(1,2,panel_i);
% axh(1) = subplot(3,1,1); hold on;                     % plot PF
celltmp1 = reshape(data(1,fce_i,dist_i,:,1:6),10,6); % 1 no -ert and 5 pert
idx_release = find(celltmp1{1,pi}.ts == 5);
t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));

% calculate the Unperturbed situation, x and f
x_avg = zeros(1, length(t_grids));
v_avg = zeros(1, length(t_grids));
f_avg = zeros(1, length(t_grids));
fp_avg= zeros(1, length(t_grids));
cti = 0; % count how many trials are added up
% fh = figure(); hold on;
hold on;
for ti = 1:1:size(celltmp1,1)
    if isempty(celltmp1{ti,1})
        continue;
    end
    cti = cti + 1;
    idx_release = find(celltmp1{ti,pi}.ts == 5);
    t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
    idx_t = find(t>=t_interest(1) & t<=t_interest(2));
    length(idx_t);
    
    lnh(ti,1) = plot(t(idx_t), celltmp1{ti,pi}.ox(2,idx_t), 'color', 'b', 'linewidth', 2);
    lnh(ti,2) = plot(t(idx_t), celltmp1{ti,pi}.x(2,idx_t), 'color', 'r', 'linewidth', 2);
    title_str = ['Force' num2str(F_list(fce_i)) 'N Stiffness' num2str(K_list(dist_i)) 'N/m pert' num2str(pi)];
    title(title_str);
end
xlabel('time (s)'); 
ylabel('position (m)')
if (panel_i == 1)
legend(lnh(1,:), {'OPTOTRAK', 'WAM'});
end
grid on;
end
linkaxes(axh, 'y'); 
ylim([0.48 0.65]);

saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/dataCommu20220309/Figure4.png');


%%
% In this scripts, check whether the position readout is accurate. 

% Position histerisis
% I clamped the robot, and do a positive-negative-positive command force,
% I hope the force-position relationship will shows a histerisis curve. 
sstmp = SessionScan(4005);

lw = 1;
t_range = [3767.5, 3771.5];
%plot(sstmp.data.t, sstmp.data.x(2,:));
t_idx = sstmp.data.t>t_range(1) & sstmp.data.t<t_range(2);

fh1 = figure('unit', 'inch', 'position', [0 0 6 3]);  % plot the command force, the censored force and the censored positions
axh(1) = subplot(2,1,1);  hold on;
plot(sstmp.data.t(t_idx), sstmp.data.Fp(2,t_idx), 'linewidth', lw);
plot(sstmp.data.t(t_idx), -sstmp.data.f(2,t_idx), 'linewidth', lw);
grid on;
xlabel('time (s)'); ylabel('Force (N)');
legend('command', 'censored');
title('command and censored force');
axh(2) = subplot(2,1,2); hold on;
plot(sstmp.data.t(t_idx), sstmp.data.x(2,t_idx), 'linewidth', lw);
plot(sstmp.data.t(t_idx), sstmp.data.opty(t_idx), 'linewidth', lw);
grid on;
legend('WAM', 'OPTOTRAK');
title('robot and motion capture position');
xlabel('time (s)'); ylabel('Position (m)');
linkaxes(axh, 'x'); 
xlim(t_range);

saveas(fh1,'/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_WAM/WAMpositionSanityCheck/forcePositionTimeSequence.png');

% fh2 = figure(2); % plot 1) the command force-censored force histerisis and 
%                  %      2) the command force-censored position histerisis
% axh(1) = subplot(1,3,1); hold on; grid on;
% plot(sstmp.data.Fp(2,t_idx), -sstmp.data.f(2,t_idx), 'Marker', '.');
% axh(2) = subplot(1,3,2); hold on; grid on;
% x_offset = sstmp.data.x(2,find(t_idx));
% x_shift = sstmp.data.x(2,t_idx)  - x_offset(1);
% plot(sstmp.data.Fp(2,t_idx), x_shift, 'Marker', '.');
% axh(3) = subplot(1,3,3); hold on; grid on;
% %lot(-sstmp.data.f(2,t_idx), x_shift, 'Marker', '.');

% figure3, plot the censored Force ~ x error historisis 
fh3 = figure('unit', 'inch', 'position', [0 0 4 4]);
axh(1) = subplot('position', [ 0.15 0.1 0.7 0.7]);
fce = sstmp.data.f(2,t_idx); 
x_e = sstmp.data.opty(t_idx) - sstmp.data.x(2,t_idx);
plot(fce, x_e, 'Marker', '.');
xlim([-30 30]); ylim([-4 4]*1e-3);
grid on; 
xlabel('force Censored (N)', 'fontSize', 12)
ylabel('x_{WAM} - x_{OPT} (m)', 'fontSize', 12)
sgtitle('hysteresis property of WAM position readout', 'fontSize', 14);
saveas(fh3,'/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_WAM/WAMpositionSanityCheck/positionHysteresis.png');