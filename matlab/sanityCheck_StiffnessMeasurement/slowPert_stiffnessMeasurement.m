% use the slow perturbation as the stiffness measurement. 
% the step last 1s, and the raising phase last for 0.5s. 
% use the \delta F/\delta x as the stiffness measurement. 
% as this measureemnt independent from the mass. This could be a easy
% measurement and easy to understand way.  

% ss3546 = SessionScan(3546);
% celltmp = ss3546.export_as_formatted_4(1);
% celltmp1 = celltmp(:,[1 4 3]); % 4 is the slow perturbation
% data(1,1,1,1,:,:) = celltmp1; 
load('data/processedData/ss3564_3569.mat', 'data');
Data = data;

tar_ALL = [2.5 5.0 7.5]
clear axhd
fh = figure(); 
colors = colormap('lines');
r = size(Data, 1); % subj
c = size(Data, 2); % direction
f = size(Data, 3); % force 
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_plot_last = 500;
if_subtract = 0;
epoc_type = 1;
plot_type = 2; % 1displacement 
stiffness = zeros(size(data));
t_step = 1/500;
fh(1) = figure(); % plot all 
fh(2) = figure(); % plot only perturb
axh = subplot(1,1,1); 
zone_1 = [-200:0];   % perturbation initial
zone_2 = [-200:0];  % perturbation ending
for ri = 1:r % subj
    for ci = 1:c % direction
        if (ci == 2 || ci == 4)
            tar_all = -tar_ALL;
        else 
            tar_all = tar_ALL;
        end
        for fi = 1:f
        for di = 1:d % target distance
            for li = 2%1:2%1:p % perturbation
                trial_num = length(Data(ri,ci,fi,di,:,li));
                for ti = 1:trial_num % each trial
                    if (isempty(Data{ri,ci,fi,di,ti,li}))
                        continue;
                    end
                    
                    switch epoc_type
                        case 1
                            idx = find(Data{ri,ci,fi,di,ti,li}.Fp(2,:)~=0);
                            %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                            idx_plot = [idx(1)+ [-200:-1] idx idx(end)+(1:idx_plot_last)]; % may error as the pert not long enough
                            if li == 1
                                display('ERROR: should use li == 2!!!');
                            end
                        case 2
                            idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
                            idx = (idx(1)):idx(end);
                    end
                    time = t_step*(idx-idx(1));
                    timep1 = t_step*(idx_plot - idx_plot(1));
                    dat_pert = Data{ri,ci,fi,di,ti,li}.Fp(2,idx);
                    dat_x = Data{ri,ci,fi,di,ti,li}.x(2,idx);
                    dat_v = Data{ri,ci,fi,di,ti,li}.v(2,idx);
                    if (ci == 2 || ci == 4)
                        dat_x = -dat_x;
                        dat_v = -dat_v;
                    end
                    delta_x = mean(Data{ri,ci,fi,di,ti,li}.x(2,idx(end)+zone_2)) - ...
                        mean(Data{ri,ci,fi,di,ti,li}.x(2,idx(1)+zone_1));
                    delta_f = mean(Data{ri,ci,fi,di,ti,li}.f(2,idx(end)+zone_2)) - ...
                        mean(Data{ri,ci,fi,di,ti,li}.f(2,idx(1)+zone_1));
                    stiffness(ri,ci,fi,di,ti,li) = delta_f/delta_x;
                    delta_x
                    delta_f
                    delta_f/delta_x
                    figure(fh(1)); % plot all data here
                    p = Data{ri,ci,fi,di,ti,li}.Fp(2,idx_plot);
                    x = Data{ri,ci,fi,di,ti,li}.x(2,idx_plot);
                    f = Data{ri,ci,fi,di,ti,li}.f(2,idx_plot);
                    subplot(3,1,1); hold on;
                    plot(timep1, p);
                    subplot(3,1,2); hold on;
                    plot(timep1, x);
                    subplot(3,1,3); hold on;
                    plot(timep1, f);
                    figure(fh(2)); % plot the selected data here
                    [~, ia1, ib1] = intersect((idx(1)+zone_1), idx_plot); 
                    [~, ia1, ib2] = intersect((idx(end)+zone_2), idx_plot); 
                    t1 = timep1(ib1);
                    p1 = Data{ri,ci,fi,di,ti,li}.Fp(2,idx(1)+zone_1);
                    x1 = Data{ri,ci,fi,di,ti,li}.x(2,idx(1)+zone_1);
                    f1 = Data{ri,ci,fi,di,ti,li}.f(2,idx(1)+zone_1);
                    t2 = timep1(ib2);
                    p2 = Data{ri,ci,fi,di,ti,li}.Fp(2,idx(end)+zone_2);
                    x2 = Data{ri,ci,fi,di,ti,li}.x(2,idx(end)+zone_2);
                    f2 = Data{ri,ci,fi,di,ti,li}.f(2,idx(end)+zone_2);
                    axh2(1,1) = subplot(3,2,1); hold on;
                    plot(t1, p1);
                    axh2(1,2) = subplot(3,2,2); hold on;
                    plot(t2, p2);
                    axh2(2,1) = subplot(3,2,3); hold on;
                    plot(t1, x1);
                    axh2(2,2) = subplot(3,2,4); hold on;
                    plot(t2, x2);
                    axh2(3,1) = subplot(3,2,5); hold on;
                    plot(t1, f1);
                    axh2(3,2) = subplot(3,2,6); hold on;
                    plot(t2, f2);
                    %linkaxes(axh2(1,:), 'y');
                    %linkaxes(axh2(2,:), 'y');
                    %linkaxes(axh2(3,:), 'y');
                    switch plot_type
                        case 2
                            dat = Data{ri,ci,fi,di,ti,li}.f(2,idx);
                            titlestr = 'force';
                        case 1
                            dat = Data{ri,ci,fi,di,ti,li}.x(2,idx);
                            dat = dat - dat(1);
                            titlestr = 'displacement';
                        case 3
                            dat = Data{ri,ci,fi,di,ti,li}.Fp(2,idx);
                            titlestr = 'Fp';
                        case 4
                            dat = Data{ri,ci,fi,di,ti,li}.v(2,idx);
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
%                     axh(1) = subplot(3,1,1);
%                     plot(time, dat_pert);
%                     title('pert_F');
%                     axh(2) = subplot(3,1,2); 
%                     plot(time, dat);
%                     title('F')
%                     axh(3) = subplot(3,1,3);
%                     title('x');
%                     plot(time, dat_x);
%                     linkaxes(axh, 'x');
%                     plot(dat_v, dat_x, 'Color', colors(4*(li-1)+di, :));
%                     plot(dat_x, dat_v, 'Color', colors(4*(li-1)+di, :));
                end
            end
        end
        end
    end
end

% xlim([-0.1 0.1])
linkaxes(axh, 'xy');
%sgtitle(titlestr);

%% show as a line
x = 1:length(stiffness(1,1,1,1,:,2));
y = reshape(stiffness(1,1,1,1,:,2), 1, length(x));
figure();
plot(x, -y, '-o', 'MarkerSize', 5);
