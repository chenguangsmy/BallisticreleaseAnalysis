%% Check the data of the relationship of XXX and velocity
% ########################### 5D data here ###########################
%load('data/processedData/ss3257_3261.mat', 'data')
% load('data/processedData/ss3307_3314.mat', 'data') % james testing at 11-03-2021
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3477_3479.mat', 'data')  % Hongwei test it with "blind" feedback
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3471_3474.mat', 'data')  % Chenguang test it with "blind" feedback 
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3480_3482.mat', 'data')
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3334_3336.mat', 'data')  % springs
%load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3334_3344.mat', 'data')
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3431_3440.mat','data') %himanshu with 6D data 
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3353_3417.mat','data') % chenguang did with 6D data
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3353_3369', 'data') % chenguang did 4 directions data, less sessions
%load('data/processedData/ss3345.mat', 'data')
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/prelimData_4subj_fine.mat', 'data') % all subjects
Data = data;
Data = data(1,4,:,:,:,:)
%Data = reshape(data(1,3,2,:,:,:), 1, 1, 3, 15, 3);
tar_ALL = [0.025 0.05, 0.075];
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(); 
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
plot_type = 1; % 1displacement 
axh(ri, ci) = subplot(1,1,1);grid on;hold on;
for ri = 1:r % subj
    for ci = 1:c % direction
        if (ci == 2 || ci == 4)
            tar_all = -tar_ALL;
        else 
            tar_all = tar_ALL;
        end
        
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        %axh(ri, ci) = subplot(c,r,r*(ci-1) + ri);grid on;hold on;
        for fi = 1:f
        for di = 1:d % target distance
        %for di = 3 % target distance
            for li = 1:2%1:p % perturbation
                trial_num = length(Data(ri,ci,fi,di,:,li));
                for ti = 1:trial_num % each trial
                    if (isempty(Data{ri,ci,fi,di,ti,li}))
                        continue;
                    end
                    
                    switch epoc_type
                        case 1
                            idx = find(Data{ri,ci,fi,di,ti,li}.Fp(2,:)~=0);
                            idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                            if li == 1
                                display('ERROR: should use li == 2!!!');
                            end
                        case 2
                            idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
                            %idx = (idx-9):idx(end);
                            idx = (idx(1)):idx(end);
                    end
                    %plot(Data{ri,ci,di,ti,li}.Fp(2,:));
                    %idx = find(Data{ri,ci,di,ti,li}.Fp(2,:)~=0);
                    %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                    %idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                    %idx = (idx-9):idx(end);
                    %idx = (idx(1)):(idx(end)+100);
                    time = t_step*(idx-idx(1));
                    dat_x = Data{ri,ci,fi,di,ti,li}.x(2,idx);
                    dat_x = dat_x - dat_x(1);
                    dat_x_sub = dat_x - tar_all(di);
                    %dat_x_sub = dat_x - dat_x(end);
                    dat_x = dat_x_sub;
                    
                    dat_v = Data{ri,ci,fi,di,ti,li}.v(2,idx);
                    dat = dat - dat(1);
                    %time = idx-idx(1);
                    if (ci == 2 || ci == 4)
                        dat_x = -dat_x;
                        dat_v = -dat_v;
                    end
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
                    plot(time, dat, 'Color', colors(4*(li-1)+di, :));
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