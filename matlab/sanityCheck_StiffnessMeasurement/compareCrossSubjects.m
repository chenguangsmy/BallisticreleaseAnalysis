% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3486_3534.mat'); % Chenguang testing
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3486_3534.mat'); % James testing  
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3491_3503.mat'); % Himanshu testing
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3615_3618.mat'); % Micheal testing
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3620_3629.mat'); % Adam testing
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3641_3644.mat'); % Marco testing
% subj_names_all = {'Chenguang', 'James', 'Himanshu', 'Michael', 'Adam', 'Marco'};
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/prelimData_6subj_fine.mat'); 
subj_names_all = {'Springs', 'Chenguang', 'James', 'Himanshu', 'Michael', 'Adam', 'Marco'};
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/prelimData_S_6subj_fine.mat'); 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/prelimData_6subj_fine_fcealign.mat'); 
%% 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. plot properties during perturbation for each subject

Data = data(1:2,:,:,:,:,:);
sbj_idx = 1;%[1:6];%[1 2];%[1:7];  %%[1 4 5 6];
subj_names = subj_names_all(sbj_idx);
Data = data(sbj_idx,:,:,:,:,:);
clear axh
Freq = 500;
t_step = 1/500;
fh = figure(); 
colors = colormap('lines');
r = size(Data, 1); %r = 2;% subj
c = size(Data, 2); c = 1; % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;
epoc_type = 2;
plot_type = 2%4;  % 1 displacement
                % 2 force
                % 3 Fp
                % 4 vlocity
ifavg = false;                
%
if (~ifavg)
for si = 1:r % subj
    for ri = 1:c % direction
        for ci = 1:f%1:2 % target force
            %axh(ri, ci) = subplot(f,c,f*(ci-1) + ci);
%             axh(1, ci) = subplot(1,f,ci); grid on; hold on; % plot on columns
%            axh(ci, 1) = subplot(f,1,ci); grid on; hold on; % plot on rows
            axh(ci, si) = subplot(f,r,(r*(ci-1)+si)); grid on; hold on; % plot on rows
            if (ci==1)
                title(subj_names{si});
            end
            for di = 1:d % target distance
%                   axh(di) = subplot(d,1,di); grid on; hold on;
%                 axh(di) = subplot(1,1,1); grid on; hold on;
                for li = 2%1:p % perturbation
                    trial_num = length(Data(si,ri,ci,di,:,li));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{si,ri,ci,di,ti,li}))
                            continue;
                        end

                        switch plot_type
                            case 1
                                dat = Data{si,ri,ci,di,ti,li}.x(xyi,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{si,ri,ci,di,ti,li}.f(xyi,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{si,ri,ci,di,ti,li}.Fp(xyi,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{si,ri,ci,di,ti,li}.v(xyi,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{si,ri,ci,di,ti,li}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{si,ri,ci,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{si,ri,ci,di,ti,li}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm force';
                            case 8 % optotrak x
                                switch length(size(Data{si,ri,ci,di,ti,li}.ox))
                                    case 2
                                        dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx);
                                    case 3
                                        dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,1);
                                end
                                titlestr = 'opto-position';
                            case 9 % optotrak v
                                switch length(size(Data{si,ri,ci,di,ti,li}.ov))
                                    case 2
                                        dat = Data{si,ri,ci,di,ti,li}.ov(xyi,idx);
                                    case 3
                                        dat = Data{si,ri,ci,di,ti,li}.ov(xyi,idx,1);
                                end
                                titlestr = 'opto-velocity';
                            case 10 
                                dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,1);
                                titlestr = 'opto1-position';
                            case 11
                                dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,2);
                                titlestr = 'opto2-position';
                            case 12 
                                dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,3);
                                titlestr = 'opto3-position';
                        end


                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end

                        plot(time, dat, 'Color', colors(4*(li-1)+di, :));
%                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                    end
                end
            end
        end
    end
end

elseif (ifavg)
    
for si = 1:r % subj
    for ri = 1:c % direction
        for ci = 1:f%1:2 % target force
            %axh(ri, ci) = subplot(f,c,f*(ci-1) + ci);
%             axh(1, ci) = subplot(1,f,ci); grid on; hold on; % plot on columns
%            axh(ci, 1) = subplot(f,1,ci); grid on; hold on; % plot on rows
            axh(ci, si) = subplot(f,r,(r*(ci-1)+si)); grid on; hold on; % plot on rows
            if (ci==1)
                title(subj_names{si});
            end
            for di = 1:d % target distance
%                   axh(di) = subplot(d,1,di); grid on; hold on;
%                 axh(di) = subplot(1,1,1); grid on; hold on;
                for li = 2%1:p % perturbation
                    trial_num = length(Data(si,ri,ci,di,:,li));
                    dat_mat = [];
                    dat_avg = [];   
                    dat_std = [];
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{si,ri,ci,di,ti,li}))
                            continue;
                        end
                        switch epoc_type
                            case 1
                                idx = find(Data{si,ri,ci,di,ti,li}.Fp(2,:)~=0);
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                            case 2
                                idx = find(Data{si,ri,ci,di,ti,li}.ts==5 | Data{si,ri,ci,di,ti,li}.ts==6);
                                if (length(idx)<300)
                                    continue
                                else
                                    idx = (idx-50):idx(300);
                                end
                                % limit idx length
                        end
                        time = t_step*(idx-idx(1));
                        switch plot_type
                            case 2
                                dat = Data{si,ri,ci,di,ti,li}.f(2,idx);
                                titlestr = 'force';
                            case 1
                                dat = Data{si,ri,ci,di,ti,li}.x(2,idx);
                                dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 3
                                dat = Data{si,ri,ci,di,ti,li}.Fp(2,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{si,ri,ci,di,ti,li}.v(2,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{si,ri,ci,di,ti,li}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{si,ri,ci,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{si,ri,ci,di,ti,li}.f(:,idx);
                                dat_submean = dat;% - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm force';
                            case 8 % optotrak x
                                switch length(size(Data{si,ri,ci,di,ti,li}.ox))
                                    case 2
                                        dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx);
                                    case 3
                                        dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,1);
                                end
                                titlestr = 'opto-position';
                            case 9 % optotrak v
                                switch length(size(Data{si,ri,ci,di,ti,li}.ov))
                                    case 2
                                        dat = Data{si,ri,ci,di,ti,li}.ov(xyi,idx);
                                    case 3
                                        dat = Data{si,ri,ci,di,ti,li}.ov(xyi,idx,1);
                                end
                                titlestr = 'opto-velocity';
                            case 10 
                                dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,1);
                                titlestr = 'opto1-position';
                            case 11
                                dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,2);
                                titlestr = 'opto2-position';
                            case 12 
                                dat = Data{si,ri,ci,di,ti,li}.ox(xyi,idx,3);
                                titlestr = 'opto3-position';
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        
                          dat_mat = [dat_mat; dat];
                          %plot(time, dat, 'Color', colors(4*(li-1)+di, :));
%                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                        %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
                    end
                    dat_avg = nanmean(dat_mat);
                    dat_std = nanstd(dat_mat);
                    % plot the center
                    plot(time, dat_avg, 'Color', colors(4*(li-1)+di, :), 'LineWidth', 3);
                    patch_x = [time time(end:-1:1)];
                    patch_y = [dat_avg+dat_std dat_avg(end:-1:1)-dat_std(end:-1:1)];
                    patch('XData', patch_x, 'YData', patch_y, ...
                        'FaceColor', colors(4*(li-1)+di, :), ...
                        'FaceAlpha', 0.3);
                    
                    %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
%                    %plot(figure(11), time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :), '--');
                end
            end
        end
    end
end

end


% 1.2, the procedure of adding axes
linkaxes(axh, 'x');
xlim([0 0.7])
if (epoc_type==1 && plot_type ==1)
    linkaxes(axh(:), 'y');
    ylim([-0.008 0.005]);
    sgtitle('Position');
elseif (epoc_type==1 && plot_type ==2)
    %ylim([-30, 0])
    linkaxes(axh(1,:),'y'); ylim(axh(1,1), [12 22]);
    linkaxes(axh(2,:),'y'); ylim(axh(2,1), [17 27]);
    linkaxes(axh(3,:),'y'); ylim(axh(3,1), [22 32]);
    %ylim([12 35])
    sgtitle('Sensored Force');
elseif (epoc_type==1 && plot_type ==4)
    %ylim([-0.2 0.2]);
    ylim([-0.1 0.1]);
    linkaxes(axh, 'y');
    sgtitle('Velocity');
elseif (epoc_type==2 && plot_type ==1)
    linkaxes(axh(1,:),'y'); 
    linkaxes(axh(2,:),'y');
    linkaxes(axh(3,:),'y');
    ylim([-0.01 0.18]);
    sgtitle('Position');
elseif    (epoc_type == 2 && plot_type == 2)
    %     ylim([-30 5]); % for force pulse
    linkaxes(axh(1,:),'y'); ylim(axh(1,1), [-5 20]);
    linkaxes(axh(2,:),'y'); ylim(axh(2,1), [-5 25]);
    linkaxes(axh(3,:),'y'); ylim(axh(3,1), [-5 30]);
    sgtitle('Sensored Force');
elseif (epoc_type==2 && plot_type ==4)
    %ylim([-0.9 0.9]);
    ylim([-0.1 0.5]);
    linkaxes(axh, 'y');
    sgtitle('Velocity');
end
%set(fh, 'Position', [183   277   180   420]);
%sgtitle(titlestr);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. plot overlap values in panels 

Data = data(1:2,:,:,:,:,:);
sbj_idx = [1 2];%[1:7]; %[1:6]; %%[1 4 5 6];
subj_names = subj_names_all(sbj_idx);
Data = data(sbj_idx,:,:,:,:,:);
clear axh
Freq = 500;
t_step = 1/500;
fh = figure(); 
colors = colormap('lines');
r = size(Data, 1); %r = 2;% subj
c = size(Data, 2); c = 1; % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;
epoc_type = 1;
plot_type_list = [1, 4, 2];
force_list = [15 20 25];
%plot_type = 2%4;  % 1 displacement
                % 2 force
                % 3 Fp
                % 4 vlocity
p = length(plot_type_list);
ifavg = false;             
ci = 1;
%
if (~ifavg)
%for ci = 1:f%1:2 % target force
for si = 1:r % subj
    for ri = 1:c % direction
        for pi = 1:p
            plot_type = plot_type_list(pi);
            axh(pi, si) = subplot(p,r,(r*(pi-1)+si)); grid on; hold on; % plot on rows
            if (pi==1)
                title(subj_names{si});
            end
            for di = 1:d % target distance
%                   axh(di) = subplot(d,1,di); grid on; hold on;
%                 axh(di) = subplot(1,1,1); grid on; hold on;
                for li = 2%1:p % perturbation
                    trial_num = length(Data(si,ri,ci,di,:,li));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{si,ri,ci,di,ti,li}))
                            continue;
                        end
                        switch epoc_type
                            case 1
                                idx = find(Data{si,ri,ci,di,ti,li}.Fp(2,:)~=0);
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                            case 2
                                idx = find(Data{si,ri,ci,di,ti,li}.ts==5 | Data{si,ri,ci,di,ti,li}.ts==6);
                                idx = (idx-100):idx(end);
%                                 idx = (idx(1)):idx(end);
                        end
                        time = t_step*(idx-idx(1));
                        switch plot_type
                            case 2
                                dat = Data{si,ri,ci,di,ti,li}.f(2,idx);
                                titlestr = 'force';
                            case 1
                                dat = Data{si,ri,ci,di,ti,li}.x(2,idx);
                                dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 3
                                dat = Data{si,ri,ci,di,ti,li}.Fp(2,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{si,ri,ci,di,ti,li}.v(2,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{si,ri,ci,di,ti,li}.tq(4,idx);
                                titlestr = 'torque4';
                            case 6
                                dat = Data{si,ri,ci,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{si,ri,ci,di,ti,li}.f(:,idx);
                                dat_submean = dat;% - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm force';
                                
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end

                        plot(time, dat, 'Color', colors(4*(li-1)+di, :));
%                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                    end
                end
            end
        end
    end
end

elseif (ifavg)
    % not changed this part to the time aligned yet...
for si = 1:r % subj
    for ri = 1:c % direction
        %for ci = 1:f%1:2 % target force
        for pi = 1:p%1:2 % target force
            plot_type = plot_type_list(pi);
            axh(pi, si) = subplot(p,r,(r*(pi-1)+si)); grid on; hold on; % plot on rows
            if (pi==1)
                title(subj_names{si});
            end
            for di = 1:d % target distance
%                   axh(di) = subplot(d,1,di); grid on; hold on;
%                 axh(di) = subplot(1,1,1); grid on; hold on;
                for li = 2%1:p % perturbation
                    trial_num = length(Data(si,ri,ci,di,:,li));
                    dat_mat = [];
                    dat_avg = [];   
                    dat_std = [];
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{si,ri,ci,di,ti,li}))
                            continue;
                        end
                        switch epoc_type
                            case 1
                                idx = find(Data{si,ri,ci,di,ti,li}.Fp(2,:)~=0);
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                            case 2
                                idx = find(Data{si,ri,ci,di,ti,li}.ts==5 | Data{si,ri,ci,di,ti,li}.ts==6);
                                if (length(idx)<300)
                                    continue
                                else
                                    idx = (idx-50):idx(300);
                                end
                                % limit idx length
                        end
                        time = t_step*(idx-idx(1));
                        switch plot_type
                            case 2
                                dat = Data{si,ri,ci,di,ti,li}.f(2,idx);
                                titlestr = 'force';
                            case 1
                                dat = Data{si,ri,ci,di,ti,li}.x(2,idx);
                                dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 3
                                dat = Data{si,ri,ci,di,ti,li}.Fp(2,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{si,ri,ci,di,ti,li}.v(2,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{si,ri,ci,di,ti,li}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{si,ri,ci,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{si,ri,ci,di,ti,li}.f(:,idx);
                                dat_submean = dat;% - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm force';
                                
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        
                          dat_mat = [dat_mat; dat];
                          %plot(time, dat, 'Color', colors(4*(li-1)+di, :));
%                         plot(time, smooth(dat, 20), 'Color', colors(4*(li-1)+di, :));
                        %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
                    end
                    dat_avg = nanmean(dat_mat);
                    dat_std = nanstd(dat_mat);
                    % plot the center
                    plot(time, dat_avg, 'Color', colors(4*(li-1)+di, :), 'LineWidth', 3);
                    patch_x = [time time(end:-1:1)];
                    patch_y = [dat_avg+dat_std dat_avg(end:-1:1)-dat_std(end:-1:1)];
                    patch('XData', patch_x, 'YData', patch_y, ...
                        'FaceColor', colors(4*(li-1)+di, :), ...
                        'FaceAlpha', 0.3);
                    title(titlestr);
                    %plot(time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :));
%                    %plot(figure(11), time, mean(dat_avg), 'Color', colors(4*(li-1)+di, :), '--');
                end
            end
        end
    end
end

end


% 1.2, the procedure of adding axes

if (epoc_type==1 && plot_type ==1)
    linkaxes(axh(:), 'y');
elseif (epoc_type==1 && sum(plot_type_list ==2))
    %ylim([-30, 0])
    linkaxes(axh(1,:),'y'); %ylim(axh(1,1), [12 22]);
    linkaxes(axh(2,:),'y'); %ylim(axh(2,1), [17 27]);
    linkaxes(axh(3,:),'y'); ylim(axh(3,1), [12 22]+(ci-1)*5);
elseif (epoc_type==1 && plot_type ==4)
    %ylim([-0.2 0.2]);
    %ylim([-0.1 0.1]);
    linkaxes(axh, 'y');
elseif (epoc_type==2 && plot_type ==1)
    linkaxes(axh(1,:),'y'); 
    linkaxes(axh(2,:),'y');
    linkaxes(axh(3,:),'y');
    %ylim([-0.01 0.18]);
elseif    (epoc_type == 2 && plot_type == 2)
    %     ylim([-30 5]); % for force pulse
    linkaxes(axh(1,:),'y'); %ylim(axh(1,1), [-5 20]);
    linkaxes(axh(2,:),'y'); %ylim(axh(2,1), [-5 25]);
    linkaxes(axh(3,:),'y'); %ylim(axh(3,1), [-5 30]);
elseif (epoc_type==2 && plot_type ==4)
    %ylim([-0.9 0.9]);
    %ylim([-0.1 0.5]);
    linkaxes(axh, 'y');
end

sgtitle(['force' num2str(force_list(ci)) 'N']);
linkaxes(axh, 'x');
xlim([0 0.7])


%% Compare the spring wih subjects...

close all;

load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/ProcessedData/ss3334_3344_alter.mat');
data_spring = data;
clear data
dexSpring = 1; % [336 N/m]
dexForce_spring = 1:3; % [15N, 20N, 25N]
dexDistance_spring = 1:3; % [will vary based on spring constant]
%             Spring [4.7, 6.25, 7.81]
depMeasures_spring = crossConditionAnalysis(data_spring,dexSpring,dexForce_spring,dexDistance_spring,'spring');

fh = figure(1);
for dex_subplot = 1:3
axh(dex_subplot) = subplot(1,3,dex_subplot);
xticklabels({'160', '320', '640'});
xlabel('stiffness (N/m)');
% plot([2 8], [160 160], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
% plot([2 8], [320 320], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
% plot([2 8], [640 640], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
end

%       160N/m      320N/m      640N/m
x0 = [  0.0938      0.0469      0.0234          % 15N
        0.1250      0.0625      0.0312          % 20N
        0.1562      0.0781      0.0391];        % 25N
x0 = x0*100;
x_axis_0 = [2.5     5       7.5];
for si = 1:3
    for fi = 1:3
        plot(x_axis_0(si)+[-1 1], [x0(fi,si) x0(fi,si)], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
    end
end
ylim([0 18]);

%% Compare across the subjects
subj_names_all = {'Chenguang', 'James', 'Himanshu', 'Michael', 'Adam', 'Marco'};
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/prelimData_6subj_fine.mat'); 
% save in the measurements 
subj_num = size(data,1);
for subj_i = 1:size(data,1)
        data_human = reshape(data(subj_i,1,:,:,:,:),1,3,3,15,3);
        dexSubject = 1; % [Chenguang]
        dexForce = 1:3; % [15N, 20N, 25N]
        dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
        depMeasures_human(subj_i) = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
end
%%%%%% plot them out 
%% 1. release -measurement
fceVal = {'15N', '20N', '25N'};
figure('Position',[300 314 900 320]);
for subj_i = 1:size(data,1)
    
    this = depMeasures_human(subj_i);
    if (subj_i == 2)
        this.k_hat_release(1,1,2,1) = nan;
        this.k_hat_release_VAF(1,1,2,1) = nan;
    end
%     k_pulse_mean = nanmean(this.k_hat_pulse(:,:,:,1:end),4);
%     k_pulse_std = nanstd(this.k_hat_pulse(:,:,:,1:end),0,4);
%     x0_pulse_mean = nanmean(this.x0_hat_pulse(:,:,:,1:end),4);
%     x0_pulse_std = nanstd(this.x0_hat_pulse(:,:,:,1:end),0,4);
%     
%     k_pulse_VAF_mean = nanmean(this.k_hat_pulse_VAF(:,:,:,1:end),4);
%     k_pulse_VAF_std = nanstd(this.k_hat_pulse_VAF(:,:,:,1:end),0,4);
%     
%     this.k_nonNan = squeeze(sum(~isnan(this.k_hat_pulse),4));
    
    k_release_mean = nanmean(this.k_hat_release(:,:,:,:),4);
    k_release_std = nanstd(this.k_hat_release(:,:,:,:),0,4);
    
    k_release_VAF_mean = nanmean(this.k_hat_release_VAF(:,:,:,1:end),4);
    k_release_VAF_std = nanstd(this.k_hat_release_VAF(:,:,:,1:end),0,4);
    
%     k_stocastic_mean = nanmean(this.k_hat_stocastic(:,:,:,:),4);
%     k_stocastic_std = nanstd(this.k_hat_stocastic(:,:,:,:),0,4);
%     
%     OC_hat_stocastic_mean = nanmean(this.OC_hat_stocastic(:,:,:,:),4);
%     OC_hat_stocastic_std = nanstd(this.OC_hat_stocastic(:,:,:,:),0,4);
    
    xRange = [1.5 8.5];
    distVal = [2.5 5.0 7.5];
    % 2. pulse -measurement
    % Stiffness Plot
    axh(subj_i) = subplot(1,subj_num,subj_i); hold on;
    
    title(subj_names_all{subj_i});
    if (subj_i == 1)
        ylabel('Stiffness (N/m)');
    else 
        yticks([500 1000]);
        yticklabels('');
    end
    ylim([0 1500]);
    xlabel('Distance');
    xticks(distVal); xlim(xRange);
    
    set(gca,'fontsize',16);grid on;
    
    for i = this.dexDirection
        errorbar(distVal(this.dexDistance),squeeze(k_release_mean(1,i,:)),squeeze(k_release_std(1,i,:)),'linewidth',2.5); hold on;
    end
    legend(fceVal);
end
sgtitle('release');

%% pulse-measurement 
fceVal = {'15N', '20N', '25N'};
figure('Position',[300 314 900 320]);
for subj_i = 1:size(data,1)
    
    this = depMeasures_human(subj_i);
    if (subj_i == 5)
        this.k_hat_pulse(1,3,1,11) = nan;
        this.k_hat_pulse_VAF(1,3,1,11) = nan;
        
        this.k_hat_pulse(1,3,3,4) = nan;
        this.k_hat_pulse_VAF(1,3,3,4) = nan;
    end
    k_pulse_mean = nanmean(this.k_hat_pulse(:,:,:,1:end),4);
    k_pulse_std = nanstd(this.k_hat_pulse(:,:,:,1:end),0,4);
    x0_pulse_mean = nanmean(this.x0_hat_pulse(:,:,:,1:end),4);
    x0_pulse_std = nanstd(this.x0_hat_pulse(:,:,:,1:end),0,4);
    
    k_pulse_VAF_mean = nanmean(this.k_hat_pulse_VAF(:,:,:,1:end),4);
    k_pulse_VAF_std = nanstd(this.k_hat_pulse_VAF(:,:,:,1:end),0,4);
    
    this.k_nonNan = squeeze(sum(~isnan(this.k_hat_pulse),4));
%     
%     k_release_mean = nanmean(this.k_hat_release(:,:,:,:),4);
%     k_release_std = nanstd(this.k_hat_release(:,:,:,:),0,4);
%     
%     k_release_VAF_mean = nanmean(this.k_hat_release_VAF(:,:,:,1:end),4);
%     k_release_VAF_std = nanstd(this.k_hat_release_VAF(:,:,:,1:end),0,4);
    
%     k_stocastic_mean = nanmean(this.k_hat_stocastic(:,:,:,:),4);
%     k_stocastic_std = nanstd(this.k_hat_stocastic(:,:,:,:),0,4);
%     
%     OC_hat_stocastic_mean = nanmean(this.OC_hat_stocastic(:,:,:,:),4);
%     OC_hat_stocastic_std = nanstd(this.OC_hat_stocastic(:,:,:,:),0,4);
    
    xRange = [1.5 8.5];
    distVal = [2.5 5.0 7.5];
    % 2. pulse -measurement
    % Stiffness Plot
    axh(subj_i) = subplot(1,subj_num,subj_i); hold on;
    
    title(subj_names_all{subj_i});
    if (subj_i == 1)
        ylabel('Stiffness (N/m)');
    else 
        yticks([500 1000]);
        yticklabels('');
    end
    ylim([0 1500]);
    xlabel('Distance');
    xticks(distVal); xlim(xRange);
    
    set(gca,'fontsize',16);grid on;
    
    for i = this.dexDirection
        errorbar(distVal(this.dexDistance),squeeze(k_pulse_mean(1,i,:)),squeeze(k_pulse_std(1,i,:)),'linewidth',2.5); hold on;
    end
    legend(fceVal);
end
sgtitle('pulse');

%% Question1: is the realease measurement measure correctly?
%   ALTERNATIVE1: USE THE 'V_MAX' ESTIMATE X0
