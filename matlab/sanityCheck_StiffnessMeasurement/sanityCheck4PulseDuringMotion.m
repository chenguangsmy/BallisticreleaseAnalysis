%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Spring data, checking whetehr the data recording is reliable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spring force pulse = 200ms, 12N

%load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
% data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);

color_arr = colormap('lines');
close all;

% 1. Check the position measurement if repeatable in the release. Plot out
% the each raw data, the average, and the raw-avg
for dist_i = 1:size(data,2) % for each spring 
    for fce_i = 1:size(data,1)
%         fh(dist_i) = figure(); hold on;
        figure(); hold on;
        for pi = 2
            % 1. get the data tobe plotted
            
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            sec = [-1 2];
            freq = 500;
            pair_t = sec(1):1/freq:sec(2);
            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate 
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 1;
                if (ifplot)
                    clf; 
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
%                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
                axh(1) = subplot(2,1,1); hold on;
                plot(t_mean, dat_interest, 'Marker', '.', 'color', [0.8 0.8 0.8]);
                plot(t_mean, dat_mean, 'color', [0 0 0]);
                axh(2) = subplot(2,1,2); hold on;
                plot(t_mean, dat_interest - dat_mean,  'color', [0.5 0.5 0.5]);
                
            end
        end
        linkaxes(axh, 'x');
        title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        title(axh(2), 'difference with unperturbed average');
    end
    
end

%% % Subject data, checking whetehr the data recording is reliable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spring force pulse = 200ms, 12N

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3913_3921.mat', 'data');
data = reshape(data(1,1,:,:,:,:), size(data, 1+[2 3 4 5]));pertT_num = size(data,4);

color_arr = colormap('lines');
close all;

% 1. Check the position measurement if repeatable in the release. Plot out
% the each raw data, the average, and the raw-avg
for dist_i = 1:size(data,2) % for each spring 
    for fce_i = 1:size(data,1)
%         fh(dist_i) = figure(); hold on;
        figure(); hold on;
        for pi = 8
            % 1. get the data tobe plotted
            
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            sec = [-1 2];
            freq = 500;
            pair_t = sec(1):1/freq:sec(2);
            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate 
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 1;
                if (ifplot)
                    clf; 
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
%                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
                axh(1) = subplot(2,1,1); hold on;
                plot(t_mean, dat_interest, 'Marker', '.', 'color', [0.8 0.8 0.8]);
                plot(t_mean, dat_mean, 'color', [0 0 0]);
                axh(2) = subplot(2,1,2); hold on;
                plot(t_mean, dat_interest - dat_mean,  'color', [0.5 0.5 0.5]);
            end
            linkaxes(axh, 'x');
            xlim(axh(1), [-0.2 0.6]);
        end
        title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        title(axh(2), 'difference with unperturbed average');
    end
    
end

%% % Subject data, checking the difference between perturb and un-perturbed ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spring force pulse = 200ms, 12N

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3913_3921.mat', 'data');
data = reshape(data(1,1,:,:,:,:), size(data, 1+[2 3 4 5]));pertT_num = size(data,4);
force_list = [15, 20, 25];
dist_list = [2.5, 5.0, 7.5];

color_arr = colormap('lines');
close all;

sec = [-0.5 1];
freq = 500;
pair_t = sec(1):1/freq:sec(2);
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3913_3921/vDiff', 'MPEG-4');
v.FrameRate = 1;
open(v);

% 1. Check the position measurement if repeatable in the release. Plot out
% the each raw data, the average, and the raw-avg
for pi = 1:8
    figure();
    for dist_i = 1:size(data,2) % for each spring
        axh(1,dist_i) = subplot(4,3,dist_i); % plot the perturb command
        celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
        idx_release = find(celltmp1{ti,pi}.ts == 5);
        t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
        idx_interest = find(t>=sec(1) & t<=sec(2));
        dat_org = celltmp1{ti,pi}.Fp(2,:);
        plot(t(idx_interest),dat_org(idx_interest));
        for fce_i = 1:size(data,1)
            axh(fce_i+1,dist_i) = subplot(4,3,(fce_i-1)*3 + dist_i + 3); hold on; grid on;
            title(['force ' num2str(force_list(frc_i)) ' distance ' num2str(dist_list(dist_i)) ]);

            % 1. get the data tobe plotted
            
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity

            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
                dat_org = celltmp1{ti,1}.v(2,idx_interest);
%                 dat_org = celltmp1{ti,1}.f(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.v(2,idx_interest);
%                 dat_org = celltmp1{ti,pi}.f(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
                %axh(1) = subplot(2,1,1); hold on;
                %plot(t_mean, dat_interest, 'Marker', '.', 'color', [0.8 0.8 0.8]);
                %plot(t_mean, dat_mean, 'color', [0 0 0]);
                %axh(2) = subplot(2,1,2); hold on;
                plot(t_mean, dat_interest - dat_mean, 'color', color_arr(4+dist_i,:));
            end
            %            linkaxes(axh, 'x');
            xlim([-0.2 0.6]);
%             ylim([-12 12]);
            ylim([-0.2 0.2]);
            %xlim(axh(1), [-0.2 0.6]);
        end
        %         title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        %         title(axh(2), 'difference with unperturbed average');
    end
    
    
    linkaxes(axh(:), 'x');
%     linkaxes(axh(:), 'xy');
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

%% % Spring data SUBTRACTION, checking the difference between perturb and un-perturbed ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spring force pulse = 200ms, 12N

% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3873_3884.mat', 'data');
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3925_3937.mat', 'data');
data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);
force_list = [15, 20, 25];
dist_list = [2.5, 5.0, 7.5];

color_arr = colormap('lines');
close all;
sec = [-1 2];
figure(); hold on;
freq = 500;
pair_t = sec(1):1/freq:sec(2);
% v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3873_3884/fDiff', 'MPEG-4');
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3925_3937/fDiff', 'MPEG-4');
v.FrameRate = 1;
open(v);

% 1. Check the position measurement if repeatable in the release. Plot out
% the each raw data, the average, and the raw-avg

for pi = 1:6
%     figure(); hold on;
    clf;
    for dist_i = 1:size(data,2) % for each spring
        axh(1,dist_i) = subplot(4,3,dist_i); % plot the perturb command 
        celltmp1 = reshape(data(1,dist_i,:,:),size(data,3),pertT_num);
        idx_release = find(celltmp1{1,1}.ts == 5);
        t = celltmp1{1,1}.t - celltmp1{1,1}.t(idx_release(1));
        idx_interest = find(t>=sec(1) & t<=sec(2));
        dat_org = celltmp1{1,pi}.Fp(2,:);
        plot(t(idx_interest),dat_org(idx_interest));
        ylabel('command force (N)');
        
        for fce_i = 1:size(data,1)
            % plot the command f 
            
            % plot the response f
            axh(fce_i+1,dist_i) = subplot(4,3,(fce_i-1)*3+dist_i + 3); 
            title(['force ' num2str(force_list(fce_i)) ' distance ' num2str(dist_list(dist_i)) ]);
            hold on;
            grid on;
            
            % 1. get the data tobe plotted
            
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
                dat_org = celltmp1{ti,1}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,1}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,pi}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
                %                 axh(1) = subplot(2,1,1); hold on;
                %                 plot(t_mean, dat_interest, 'Marker', '.', 'color', [0.8 0.8 0.8]);
                %                 plot(t_mean, dat_mean, 'color', [0 0 0]);
                %                 axh(2) = subplot(2,1,2); hold on;
                plot(t_mean, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
                
%                 ylim([-0.5 0.5]); ylabel('velocity (m/s)');% velocity
                ylim([-12 12]); ylabel('censored force (N)');% force
            end
        end
        %         linkaxes(axh, 'x');
        %         title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        %         title(axh(2), 'difference with unperturbed average');
    end
   linkaxes(axh(:), 'xy'); sgtitle('force subtraction'); % force
%    linkaxes(axh(:), 'x'); sgtitle('velocity subtraction');% velocity
   set(gcf, 'position', [0,0, 1080, 680]);
   fname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3925_3937/fDiffp' num2str(pi) '.png'];
   saveas(gcf, fname);
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);

%% % Spring data, Another way to plot the difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spring force pulse = 200ms, 12N

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3925_3937.mat', 'data');
data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);
fce_list = [15, 20, 25];
dist_list = [2.5, 5.0, 7.5];
pert_list = [200 400 600 800 1000];
color_arr = colormap('lines');
close all;
sec = [-1 2];
% sec = [-0.4 0.8];
figure(); hold on;
freq = 500;
pair_t = sec(1):1/freq:sec(2);

clf;
for dist_i = 1:size(data,2) % for each spring
    pi = 1;
%     axh(1,(fce_i-1)*3+dist_i) = subplot(8,9,(fce_i-1)*3+dist_i); % plot the perturb command
    %         celltmp1 = reshape(data(1,dist_i,:,:),size(data,3),pertT_num);
    celltmp1 = reshape(data(1,dist_i,1:5,:),5,pertT_num);   % only use first 5 trials
    idx_release = find(celltmp1{1,pi}.ts == 5);
    t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
    idx_interest = find(t>=sec(1) & t<=sec(2));
    dat_org = celltmp1{1,pi}.Fp(2,:);
%     plot(t(idx_interest),dat_org(idx_interest));
%     ylabel('command force (N)');
    pertIdx = find(abs(celltmp1{1,pi}.Fp(2,:)) > max(abs(celltmp1{1,pi}.Fp(2,:)) * 0.05));
    if (~isempty(pertIdx))
%         xline(t(pertIdx(1)));
    end
    
    for fce_i = 1:size(data,1)
        % plot the command f
        for pi = 2:6
            % plot the response f
            %             axh(fce_i+1,dist_i) = subplot(4,3,(fce_i-1)*3+dist_i + 3);
            axh(pi-1,(dist_i-1)*3+fce_i) = subplot(5,9,(pi-2)*9+(dist_i-1)*3+fce_i);
            %             title(['force ' num2str(force_list(fce_i)) ' distance ' num2str(dist_list(dist_i)) ]);
            hold on;
            grid on;
            
            % 1. get the data tobe plotted
            
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
                dat_org = celltmp1{ti,1}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,1}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                pertIdx = find(abs(celltmp1{ti,pi}.Fp(2,:)) > max(abs(celltmp1{ti,pi}.Fp(2,:)) * 0.05));
                t_pert = [];
                if (~isempty(pertIdx))
                    t_pert = t(pertIdx(1));
                else
                    disp(['F' num2str(fce_i) 'D' num2str(dist_i) 'p' num2str(pi)]);
                end
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,pi}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
%                 plot(t_mean, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
                plot(t_mean-t_pert, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
%                 ylim([-0.5 0.5]); ylabel('velocity (m/s)');% velocity
                %                 ylim([-12 12]);  ylabel('censored force (N)');% force
                if (~isempty(t_pert) && ti==1)
%                     xline(t_pert);
                    xline(0);
                end
                
            end
%             xlim([t_pert - 0.1, t_pert + 0.8]);
%             xlim([ - 0.1, + 0.6]);
            xlim([ - 0.1, +1.5]);
            % off the xlabelTick if not in the last row
            if(pi~=6) 
                set(gca,'xTickLabel', {''});
                if (pi == 2)
                    title(['force' num2str(fce_list(fce_i)) ' dist' num2str(dist_list(dist_i)) ] );
                end
            else 
                xlabel('t (s)');
            end
            % off the ylabelTick if not in the first column
            if (~(dist_i == 1 && fce_i == 1))
                set(gca,'yTickLabel', {''});
            else 
                ylabel(['pert' num2str(pert_list(pi-1)) 'ms']);
            end
        end
        %         linkaxes(axh, 'x');
        %         title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        %         title(axh(2), 'difference with unperturbed average');
    end

end
    set(gcf, 'position', [0,0, 1080, 680]);
    linkaxes(axh(:), 'y');
%     sgtitle('Velocity Subtraction over movement')% velocity
    sgtitle('Force Subtraction over movement')% velocity
    figname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3925_3937/fdiffpOverlap.png'];
    saveas(gcf, figname);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New data trying to avoid the messy before-perturb effect 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpDatarg2(3923);


sec = [-1 2];
figure(); hold on;
color_arr = colormap('lines');
freq = 500;
pair_t = sec(1):1/freq:sec(2);
pertT_num = size(data,5);
for dist_i = 1:3 % for each spring
    %         axh(1,dist_i) = subplot(4,3,dist_i); % plot the perturb command
    %         celltmp1 = reshape(data(1,dist_i,:,:),size(data,3),pertT_num);
    %         idx_release = find(celltmp1{1,1}.ts == 5);
    %         t = celltmp1{1,1}.t - celltmp1{1,1}.t(idx_release(1));
    %         idx_interest = find(t>=sec(1) & t<=sec(2));
    %         dat_org = celltmp1{1,pi}.Fp(2,:);
    %         plot(t(idx_interest),dat_org(idx_interest));
    
    for fce_i = 1:3
        celltmp1 = reshape(data(1,fce_i,dist_i,:,:),size(data,4),pertT_num);
        idx_release = find(celltmp1{1,1}.ts == 5);
        t = celltmp1{1,1}.t - celltmp1{1,1}.t(idx_release(1));
        idx_interest = find(t>=sec(1) & t<=sec(2));
        dat_org = celltmp1{1,pi}.Fp(2,:);
        % plot the command f
        for pi = 1:8
            figure(); hold on;
            clf;
            % plot the response f
            %             axh(fce_i+1,dist_i) = subplot(4,3,(fce_i-1)*3+dist_i + 3);
            title(['force ' num2str(force_list(fce_i)) ' distance ' num2str(dist_list(dist_i)) ]);
            hold on;
            grid on;
            
            % 1. get the data tobe plotted
            
            %             celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:(size(celltmp1,1)-2)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
                dat_org = celltmp1{ti,1}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,1}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:(size(celltmp1,1)-2)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,pi}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
                %                 axh(1) = subplot(2,1,1); hold on;
                %                 plot(t_mean, dat_interest, 'Marker', '.', 'color', [0.8 0.8 0.8]);
                %                 plot(t_mean, dat_mean, 'color', [0 0 0]);
                %                 axh(2) = subplot(2,1,2); hold on;
                plot(t_mean, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
%                 ylim([-0.5 0.5]); % velocity
                ylim([-12 12]); % force
            end
        end
        %         linkaxes(axh, 'x');
        %         title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        %         title(axh(2), 'difference with unperturbed average');
    end
    % %    linkaxes(axh(:), 'xy');  % force
    %    linkaxes(axh(:), 'x'); % velocity
    %    frame = getframe(gcf);
    %    writeVideo(v,frame);
end
% close(v);

%% New subject data trying to avoid the predictio neffect 
% plot the Fp, f, v in one plot to see the difference 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cpDatarg2(3923);
% load('data/processedData/ss3938_3939.mat', 'data'); % 12N perturbation, various time
load('data/processedData/ss3938_3949.mat', 'data'); 
sec = [-1 2];
figure(); hold on;
color_arr = colormap('lines');
close all;
freq = 500;
pair_t = sec(1):1/freq:sec(2);
pertT_num = size(data,5);
for dist_i = 1:3 % for each spring
    %         axh(1,dist_i) = subplot(4,3,dist_i); % plot the perturb command
    %         celltmp1 = reshape(data(1,dist_i,:,:),size(data,3),pertT_num);
    %         idx_release = find(celltmp1{1,1}.ts == 5);
    %         t = celltmp1{1,1}.t - celltmp1{1,1}.t(idx_release(1));
    %         idx_interest = find(t>=sec(1) & t<=sec(2));
    %         dat_org = celltmp1{1,pi}.Fp(2,:);
    %         plot(t(idx_interest),dat_org(idx_interest));
    
    for fce_i = 1:3
        celltmp1 = reshape(data(1,fce_i,dist_i,:,:),size(data,4),pertT_num);
        idx_release = find(celltmp1{1,1}.ts == 5);
        t = celltmp1{1,1}.t - celltmp1{1,1}.t(idx_release(1));
        idx_interest = find(t>=sec(1) & t<=sec(2));
        dat_org = celltmp1{1,pi}.Fp(2,:);
        % plot the command f
        for pi = 1:8
%             figure('Position',[0 0 100 400]); hold on;
            figure('Position',[0 0 600 600]); hold on;
            clf;
            % plot the response f
            %             axh(fce_i+1,dist_i) = subplot(4,3,(fce_i-1)*3+dist_i + 3);
            sgtitle(['force ' num2str(force_list(fce_i)) ' distance ' num2str(dist_list(dist_i)) ]);
            hold on;
            grid on;
            
            % 1. get the data tobe plotted
            
            %             celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            dat_meanf = zeros(1,length(pair_t));
            dat_mean_matf = [];
            dat_meanv = zeros(1,length(pair_t));
            dat_mean_matv = [];
            for ti = 1:(size(celltmp1,1)-2)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
                dat_org = celltmp1{ti,1}.f(2,idx_interest);
                dat_interestf = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                dat_org = celltmp1{ti,1}.v(2,idx_interest);
                dat_interestv = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_matf = [dat_mean_matf; dat_interestf];
                dat_mean_matv = [dat_mean_matv; dat_interestv];
            end
            dat_meanf = mean(dat_mean_matf);
            dat_meanv = mean(dat_mean_matv);
            
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:(size(celltmp1,1)-2)
                if isempty(celltmp1{ti,pi})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                
                axh(1) = subplot(3,1,1);  hold on;
                dat_org = celltmp1{ti,pi}.Fp(2,idx_interest);
                plot(t(idx_interest), dat_org);
                
                
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,pi}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
                %                 axh(1) = subplot(2,1,1); hold on;
                %                 plot(t_mean, dat_interest, 'Marker', '.', 'color', [0.8 0.8 0.8]);
                %                 plot(t_mean, dat_mean, 'color', [0 0 0]);
                %                 axh(2) = subplot(2,1,2); hold on;
                axh(2) = subplot(3,1,2); hold on;
                plot(t_mean, dat_interest - dat_meanf,  'color', color_arr(4+dist_i,:));
                
                axh(3) = subplot(3,1,3); hold on;
                dat_org = celltmp1{ti,pi}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                plot(t_mean, dat_interest - dat_meanv,  'color', color_arr(4+dist_i,:));
                
                linkaxes(axh, 'x');
                xlim([-0.1 0.6]);
                ylim([-0.5 0.5]); % velocity
                ylim(axh(2),[-12 12]); % force
            end
        end
        %         linkaxes(axh, 'x');
        %         title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        %         title(axh(2), 'difference with unperturbed average');
    end
    % %    linkaxes(axh(:), 'xy');  % force
    %    linkaxes(axh(:), 'x'); % velocity
    %    frame = getframe(gcf);
    %    writeVideo(v,frame);
end
% close(v);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The subject data, compare the perturbed one and unperturbed ones

fce_list = [15 20 25];
dist_list = [2.5 5.0 7.5];
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
color_arr = colormap('lines');
pertT_num = 1 + 7;     % 1 without pert, and 12 perturbation time, and 1 stoc pert
close all;
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3938_3949/v_compare.mp4', 'MPEG-4');
v.FrameRate = 1;
open(v);

fh(1) = figure();
hold on;
for pi = 1:(pertT_num)
    clf;
    for dist_i = 1:size(data,3) % for each spring
        % 1. plot the perturbed force in the first panel
        trial_num = size(data,4); % 4 for spring data, 5 for human data
        axh(1,dist_i) = subplot(4,3,dist_i); hold on;                     % plot PF
%         celltmp1 = reshape(data(1,1,dist_i,:,:),trial_num,pertT_num);
        celltmp1 = reshape(data(1,1,dist_i,1:5,:),5,pertT_num);
        idx_release = find(celltmp1{1,pi}.ts == 5);
        t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
        plot(t, celltmp1{1,pi}.Fp(2,:), 'color', color_arr(1,:));
        ylabel('command force (N)');
        pertIdx = find(abs(celltmp1{1,pi}.Fp(2,:)) > max(abs(celltmp1{1,pi}.Fp(2,:)) * 0.05));
        if (~isempty(pertIdx))
            xline(t(pertIdx(1)));
        end
        
        % 2. plot the perturbed velocity in the following panels
        for fce_i = 1:size(data,2)
            axh(dist_i+1,fce_i) = subplot(4,3,dist_i+(fce_i*3)); hold on;         % plot each response
            celltmp1 = reshape(data(1,fce_i,dist_i,:,:),trial_num,pertT_num);
            % 2.1 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                plot(t, celltmp1{ti,1}.v(2,:), 'color', [0.5 0.5 0.5]);
%                 plot(t, celltmp1{ti,1}.f(2,:), 'color', [0.5 0.5 0.5]);
%                 pertIdx = find(abs(celltmp1{ti,1}.Fp(2,:)) > max(abs(celltmp1{1,1}.Fp(2,:)) * 0.05));
%                 if (~isempty(pertIdx))
%                     xline(t(pertIdx(1)));
%                 end
            end
            % 2.2 plot the perturbed one, -perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi}) || pi == 1
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                plot(t, celltmp1{ti,pi}.v(2,:), 'color', color_arr(dist_i+4,:));
%                 plot(t, celltmp1{ti,pi}.f(2,:), 'color', color_arr(dist_i+4,:));
                pertIdx = find(abs(celltmp1{ti,pi}.Fp(2,:)) > max(abs(celltmp1{1,pi}.Fp(2,:)) * 0.05));
                if (~isempty(pertIdx))
                    xline(t(pertIdx(1)));
                end
            end
        end
    end
        % plot notes here: 
        linkaxes(axh(:), 'x'); xlim([-0.1 0.8]);
        
%         sgtitle_str = 'unperturbed and perturbed force'; label_str = 'cnesored force (N)';
%         linkaxes(axh(2:end,:), 'y'); ylim([-12 28]); %force 
        sgtitle_str = 'unperturbed and perturbed velocity'; label_str = 'velocity (m/s)';
        linkaxes(axh(2:end,:), 'y'); ylim([-0.3 0.6]); % velocity
        
        for fce_i = 1:3
            for dist_i = 1:3
            subplot(4,3,fce_i*3 + dist_i); 
            title_str = (['force' num2str(fce_list(fce_i)) 'N dist' num2str(dist_list(dist_i)) 'cm']); 
            ylabel('N');
            subplot(4,3,fce_i*3 + dist_i); title(title_str); ylabel(label_str);
            xlabel('time');
            end
        end
        set(gcf, 'position', [0,0, 1080, 680]);
        
        sgtitle(sgtitle_str);
        fname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3938_3949/v_comparep' num2str(pi) '.png'];
        saveas(gcf, fname);
        frame = getframe(gcf);
        sgtitle(sgtitle_str);
        writeVideo(v,frame);
end
close(v);  

%% % Subject data, checking the difference between perturb and un-perturbed ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spring force pulse = 200ms, 12N

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);
force_list = [15, 20, 25];
dist_list = [2.5, 5.0, 7.5];

color_arr = colormap('lines');
close all;
% sec = [-1 2];
sec = [-0.4 0.8];
figure(); hold on;
freq = 500;
pair_t = sec(1):1/freq:sec(2);
v = VideoWriter('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3938_3949/vDiff', 'MPEG-4');
v.FrameRate = 1;
open(v);

% 1. Check the position measurement if repeatable in the release. Plot out
% the each raw data, the average, and the raw-avg

% for pi = 1:6
for pi = 1:8
%     figure(); hold on;
    clf;
    for dist_i = 1:size(data,2) % for each spring
        axh(1,dist_i) = subplot(4,3,dist_i); % plot the perturb command 
%         celltmp1 = reshape(data(1,dist_i,:,:),size(data,3),pertT_num);
        celltmp1 = reshape(data(1,dist_i,1:5,:),5,pertT_num);   % only use first 5 trials
        idx_release = find(celltmp1{1,pi}.ts == 5);
        t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
        idx_interest = find(t>=sec(1) & t<=sec(2));
        dat_org = celltmp1{1,pi}.Fp(2,:);
        plot(t(idx_interest),dat_org(idx_interest));
        ylabel('command force (N)');
        pertIdx = find(abs(celltmp1{1,pi}.Fp(2,:)) > max(abs(celltmp1{1,pi}.Fp(2,:)) * 0.05));
        if (~isempty(pertIdx))
            xline(t(pertIdx(1)));
        end
        
        for fce_i = 1:size(data,1)
            % plot the command f 
            
            % plot the response f
            axh(fce_i+1,dist_i) = subplot(4,3,(fce_i-1)*3+dist_i + 3); 
            title(['force ' num2str(force_list(fce_i)) ' distance ' num2str(dist_list(dist_i)) ]);
            hold on;
            grid on;
            
            % 1. get the data tobe plotted
            
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
%                 dat_org = celltmp1{ti,1}.f(2,idx_interest);
                dat_org = celltmp1{ti,1}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                pertIdx = find(abs(celltmp1{ti,pi}.Fp(2,:)) > max(abs(celltmp1{ti,pi}.Fp(2,:)) * 0.05));
                t_pert = [];
                if (~isempty(pertIdx))
                    t_pert = t(pertIdx(1));
                end
        
                idx_interest = find(t>=sec(1) & t<=sec(2));
%                 dat_org = celltmp1{ti,pi}.f(2,idx_interest);
                dat_org = celltmp1{ti,pi}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
                %                 axh(1) = subplot(2,1,1); hold on;
                %                 plot(t_mean, dat_interest, 'Marker', '.', 'color', [0.8 0.8 0.8]);
                %                 plot(t_mean, dat_mean, 'color', [0 0 0]);
                %                 axh(2) = subplot(2,1,2); hold on;
                plot(t_mean, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
                ylim([-0.5 0.5]); ylabel('velocity (m/s)');% velocity 
%                 ylim([-12 12]);  ylabel('censored force (N)');% force
                if (~isempty(t_pert))
                    xline(t_pert);
                end
        
            end
        end
        %         linkaxes(axh, 'x');
        %         title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        %         title(axh(2), 'difference with unperturbed average');
    end
    set(gcf, 'position', [0,0, 1080, 680]);
%    linkaxes(axh(:), 'xy'); sgtitle('Force Subtraction') % force
   linkaxes(axh(:), 'x'); sgtitle('Velocity Subtraction')% velocity
    
    figname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3938_3949/vdiffp' num2str(pi) '.png'];
    saveas(gcf, figname);
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);

%% % Subject data, Another way to plot the difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spring force pulse = 200ms, 12N

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3938_3949.mat', 'data');
data = reshape(data(1,:,:,:,:,:), size(data, [2 3 4 5 6]));
data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);
force_list = [15, 20, 25];
dist_list = [2.5, 5.0, 7.5];
pert_list = [100 125 150 175 200 225 250];
color_arr = colormap('lines');
close all;
% sec = [-1 2];
sec = [-0.4 0.8];
figure(); hold on;
freq = 500;
pair_t = sec(1):1/freq:sec(2);

clf;
for dist_i = 1:size(data,2) % for each spring
    pi = 1; % ... could be a unproper name
%     axh(1,(fce_i-1)*3+dist_i) = subplot(8,9,(fce_i-1)*3+dist_i); % plot the perturb command
    %         celltmp1 = reshape(data(1,dist_i,:,:),size(data,3),pertT_num);
    celltmp1 = reshape(data(1,dist_i,1:5,:),5,pertT_num);   % only use first 5 trials
    idx_release = find(celltmp1{1,pi}.ts == 5);
    t = celltmp1{1,pi}.t - celltmp1{1,pi}.t(idx_release(1));
    idx_interest = find(t>=sec(1) & t<=sec(2));
    dat_org = celltmp1{1,pi}.Fp(2,:);
%     plot(t(idx_interest),dat_org(idx_interest));
%     ylabel('command force (N)');
    pertIdx = find(abs(celltmp1{1,pi}.Fp(2,:)) > max(abs(celltmp1{1,pi}.Fp(2,:)) * 0.05));
    if (~isempty(pertIdx))
%         xline(t(pertIdx(1)));
    end
    
    for fce_i = 1:size(data,1)
        % plot the command f
        for pi = 2:8
            % plot the response f
            %             axh(fce_i+1,dist_i) = subplot(4,3,(fce_i-1)*3+dist_i + 3);
            axh(pi-1,(dist_i-1)*3+fce_i) = subplot(7,9,(pi-2)*9+(dist_i-1)*3+fce_i);
            %             title(['force ' num2str(force_list(fce_i)) ' distance ' num2str(dist_list(dist_i)) ]);
            hold on;
            grid on;
            
            % 1. get the data tobe plotted
            
            celltmp1 = reshape(data(fce_i,dist_i,:,:),size(data,3),pertT_num);
            % Get the un-perturbed avg velocity
            dat_mean = zeros(1,length(pair_t));
            dat_mean_mat = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
                dat_org = celltmp1{ti,1}.f(2,idx_interest);
%                 dat_org = celltmp1{ti,1}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat = [dat_mean_mat; dat_interest];
            end
            dat_mean = mean(dat_mean_mat);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                pertIdx = find(abs(celltmp1{ti,pi}.Fp(2,:)) > max(abs(celltmp1{ti,pi}.Fp(2,:)) * 0.05));
                t_pert = [];
                if (~isempty(pertIdx))
                    t_pert = t(pertIdx(1));
                else
                    disp(['F' num2str(fce_i) 'D' num2str(dist_i) 'p' num2str(pi)]);
                end
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.f(2,idx_interest);
                %dat_org = celltmp1{ti,pi}.v(2,idx_interest);
                dat_interest = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                
%                 plot(t_mean, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
                plot(t_mean-t_pert, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
%                 ylim([-0.5 0.5]); ylabel('velocity (m/s)');% velocity
                %                 ylim([-12 12]);  ylabel('censored force (N)');% force
                if (~isempty(t_pert) && ti==1)
%                     xline(t_pert);
                    xline(0);
                end
                
            end
%             xlim([t_pert - 0.1, t_pert + 0.8]);
            xlim([ - 0.1, + 0.6]);
            % off the xlabelTick if not in the last row
            if(pi~=8) 
                set(gca,'xTickLabel', {''});
                if (pi == 2)
                    title(['force' num2str(force_list(fce_i)) ' dist' num2str(dist_list(dist_i)) ] );
                end
            else 
                xlabel('t (s)');
            end
            % off the ylabelTick if not in the first column
            if (~(dist_i == 1 && fce_i == 1))
                set(gca,'yTickLabel', {''});
            else 
                ylabel(['pert' num2str(pert_list(pi-1)) 'ms']);
            end
        end
        %         linkaxes(axh, 'x');
        %         title(axh(1), ['force' num2str(fce_i), 'dist' num2str(dist_i)]);
        %         title(axh(2), 'difference with unperturbed average');
    end

end
    set(gcf, 'position', [0,0, 1080, 680]);
    linkaxes(axh(:), 'y');
%     sgtitle('Velocity Subtraction over movement')% velocity
    sgtitle('Force Subtraction over movement')% velocity
    figname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/dataDescriptions/ss3938_3949/fdiffpOverlap.png'];
    saveas(gcf, figname);