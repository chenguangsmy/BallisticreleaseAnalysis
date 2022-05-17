% Ballistic Release right

% export data and see how it look like 
%%% Chenguang tested with right directional ballistic release 

% ss_num =...
%         [   4130    4129    4131
%             4133    4134    4132
%             4137    4135    4136];

% ss_num =... %update with the new task design and 3 markers % wrongData!!!
%          [   4155    4156    4147
%              4148    4149    4150
%              4151    4152    4153];

ss_num = [  4169  4170  4171
            4175  4176  4174
            4179  4177  4178];
% Run the sessionScanOPT to make sure OPT was converted 
for i = 1:length(ss_num(:))
    ss_list = ss_num(:);
    SessionScanOPT(ss_list(i));
%     SessionScan(ss_list(i));
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% export the data into specified format 
% subj - direction - force - displacement - trial - pert 
% pert: 1(0). no perturbation; 
%       2(1). perturb before release;
%       3(2). perturb after release (during movement)
%       4(3). perturb at hold period 

%%% This part of data, when the ss_num is a matrix
data = cell(1, 1, 3, 3, 15, 4);

for fce_i = 1:size(ss_num,1)
    for tar_i = 1:size(ss_num,2) % step perts
        ss_tmp = SessionScan(ss_num(fce_i, tar_i));
        celltmp1 = ss_tmp.export_as_formatted(0);
        celltmp = reshape(celltmp1(1,1,:,:), size(celltmp1, [3 4]));
        trials_num = size(celltmp,1);
        % seperate perturbations 
        trial_pert_idx = zeros(size(celltmp,1),1); 
        for trial_i = 1:size(celltmp,1)
            if isempty(celltmp{trial_i,2})
                continue;
            end
            trialtmp = celltmp{trial_i,2};
            [~,pert_max_idx] = max(abs(trialtmp.Fp(1,:)));
            ifplot = 1;
            if (ifplot)
                clf;
                axh(1) = subplot(2,1,1);
                plot(trialtmp.t, trialtmp.ts);
                axh(2) = subplot(2,1,2);
                plot(trialtmp.t, trialtmp.Fp(1,:));
                linkaxes(axh, 'x');
            end
            switch (trialtmp.ts(pert_max_idx))
                case 4  % pert at force exertion
                    trial_pert_idx(trial_i) = 1;
                case 5  % pert at release
                    trial_pert_idx(trial_i) = 2;
                case 6  % pert at hold
                    trial_pert_idx(trial_i) = 3;
            end
        end
        

    % tidy them according to the perturb conditions 
         % unperturbed trials 
        if trials_num>15
            data(1,1, fce_i, tar_i,:,1) = celltmp(1:15,1);
            
        else
            data(1,1, fce_i, tar_i,1:trials_num,1) = celltmp(:,1);
        end
        % perturbed trial
        for pert_type = 1:3 
            pert_trial_idx = trial_pert_idx == pert_type;
            pert_trials_num = sum(pert_trial_idx);
            if pert_trials_num > 5
                pert_trials_num = 5;
                pert_trial_idx = find(pert_trial_idx);
                pert_trial_idx = pert_trial_idx(1:5);
            end
            data(1,1, fce_i, tar_i,1:pert_trials_num, pert_type+1) = celltmp(pert_trial_idx,2);
        end
    end
end

% save('data/processedData/ss4129_4137.mat', 'data')
% save('data/processedData/ss4147_4156.mat', 'data')
save('data/processedData/ss4169_4179.mat', 'data')
%% 
% concatinate two subjects 
clear; 
load('data/processedData/ss4129_4137.mat', 'data');
data1 = data;
load('data/processedData/ss4169_4179.mat', 'data');
data2 = data; 
data(1,:,:,:,:,:) = data1;
data(2,:,:,:,:,:) = data2;
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4129_4179.mat', 'data')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% export the data into specified format 
% subj - direction - force - displacement - trial - pert 
% pert: 1(0). no perturbation; 
%       2(1). perturb before release;
%       3(2). perturb after release (during movement)
%       4(3). perturb at hold period 
%%% This part of data, when the ss_num is a matrix cell

data = cell(1, 1, 3, 3, 15, 4);

ss_num =... %update with the new task design and 3 markers 
         {   4159    4160    4161
             4163    4164    4162
             4168    [4165 4166]    4167};

for fce_i = 1:size(ss_num,1)
    for tar_i = 1:size(ss_num,2) % step perts
        celltmp_cat = {};
        for session_i = 1:length(ss_num{fce_i, tar_i})
            ss_tmp = SessionScan(ss_num{fce_i, tar_i}(session_i));
            celltmp1 = ss_tmp.export_as_formatted(1);
            celltmp = reshape(celltmp1(1,1,:,:), size(celltmp1, [3 4]));
            celltmp_cat = [celltmp_cat; celltmp];
        end
        celltmp = celltmp_cat;
        trials_num = size(celltmp,1);
        % seperate perturbations 
        trial_pert_idx = zeros(size(celltmp,1),1); 
        for trial_i = 1:size(celltmp,1)
            if isempty(celltmp{trial_i,2})
                continue;
            end
            trialtmp = celltmp{trial_i,2};
            [~,pert_max_idx] = max(abs(trialtmp.Fp(1,:)));
            ifplot = 1;
            if (ifplot)
                clf;
                axh(1) = subplot(2,1,1)
                plot(trialtmp.t, trialtmp.ts);
                axh(2) = subplot(2,1,2)
                plot(trialtmp.t, trialtmp.Fp(1,:));
                linkaxes(axh, 'x');
            end
            switch (trialtmp.ts(pert_max_idx))
                case 4  % pert at force exertion
                    trial_pert_idx(trial_i) = 1;
                case 5  % pert at release
                    trial_pert_idx(trial_i) = 2;
                case 6  % pert at hold
                    trial_pert_idx(trial_i) = 3;
            end
        end
        

    % tidy them according to the perturb conditions 
         % unperturbed trials 
        if trials_num>15
            data(1,1, fce_i, tar_i,:,1) = celltmp(1:15,1);
            
        else
            data(1,1, fce_i, tar_i,1:trials_num,1) = celltmp(:,1);
        end
        % perturbed trial
        for pert_type = 1:3 
            pert_trial_idx = trial_pert_idx == pert_type;
            pert_trials_num = sum(pert_trial_idx);
            if pert_trials_num > 5
                pert_trials_num = 5;
                pert_trial_idx = find(pert_trial_idx);
                pert_trial_idx = pert_trial_idx(1:5);
            end
            data(1,1, fce_i, tar_i,1:pert_trials_num, pert_type+1) = celltmp(pert_trial_idx,2);
        end
    end
end

save('data/processedData/ss4159_4168.mat', 'data')

%% 
% Check the data, plot out the figure

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display here 
dist_list = [2.5 5.0 7.5];
fce_list = [15 20 25];
% load('data/processedData/ss4129_4137.mat', 'data');  % 6N perturbation on
% x, require 800ms
% load('data/processedData/ss4147_4156.mat', 'data');   % adapted task & 3 cursors
% load('data/processedData/ss4159_4168.mat', 'data');     % require 500ms
% load('data/processedData/ss4169_4179.mat', 'data');     % Himanshu
load('data/processedData/ss4129_4179.mat', 'data');     % Himanshu
Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 
colors = colormap('lines');
close(fh);
r = size(Data, 1); % subj
d = size(Data, 2); % direction
f = size(Data, 3); % force
l = size(Data, 4); % distance
t = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;

epoc_type = 2;  % 1 perturb
                % 2 release
plot_type = 2;  % 1 displacement
                % 2 force 
                % 3 feedforward force
                % 4 velocity
                % 5 torque
                % 6 norm displacement
                % 7 norm force
                % 8 opt displacement
                % 9 opt velocity
                % 10 opt 1
                % 11 opt 2
                % 12 opt 3
pert_type = 1; % choose option [2 3 4]
axh = zeros(d,r);
xyi = 1;        % x, y

fh = figure();
% pert_type = 2; colors = colors(4:end,:)
for ri = 1%1:r % subj
    for di = 1%:d % direction
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        %         axh(ri, di) = subplot(d,r,r*(di-1) + ri);grid on;hold on;
        %         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
        for fi = 1:f % target force
            for li = 1:3 % target distance
                axh(fi, li) = subplot(f,l,f*(li-1) + fi);grid on;hold on; %?
                for pi = pert_type%1:p % perturbation
                    trial_num = length(Data(ri,di,fi,li,:,pi));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,di,fi,li,ti,pi}))
                            continue;
                        end
                        
                        switch epoc_type
                            case 1
%                                 idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0 & Data{ri,di,fi,li,ti,pi}.ts==4);  % pert at y
                                idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0);
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if pi == 1
                                    disp('ERROR: should use pi == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,di,fi,li,ti,pi}.ts==5 | Data{ri,di,fi,li,ti,pi}.ts==6);
                                idx = (idx(1)-100):idx(end);
                                %idx = (idx(1)):idx(end);
                        end
                        %plot(Data{ri,ci,di,ti,li}.Fp(xyi,:));
                        %idx = find(Data{ri,ci,di,ti,li}.Fp(xyi,:)~=0);
                        %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                        %idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                        %idx = (idx-9):idx(end);
                        %idx = (idx(1)):(idx(end)+100);
                        time = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,di,fi,li,ti,pi}.x(xyi,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,di,fi,li,ti,pi}.f(xyi,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,di,fi,li,ti,pi}.Fp(xyi,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,di,fi,li,ti,pi}.v(xyi,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,di,fi,li,ti,pi}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,di,fi,li,ti,pi}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,di,fi,li,ti,pi}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm force';
                            case 8 % optotrak x
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ox))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                end
                                titlestr = 'opto-position';
                            case 9 % optotrak v
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ov))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx,1);
                                end
                                titlestr = 'opto-velocity';
                            case 10 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                titlestr = 'opto1-position';
                            case 11
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,2);
                                titlestr = 'opto2-position';
                            case 12 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,3);
                                titlestr = 'opto3-position';
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4+li, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
                title(['fce' num2str(fce_list(fi)) 'dist' num2str(dist_list(li))]);
            end
        end
    end
end
linkaxes(axh, 'xy');
xlim([0 1.0])
% xlim([0 0.5]);
% xlim([0 2])
sgtitle(titlestr);

%% do the dF/dx during pulse 
load('data/processedData/ss4129_4137.mat', 'data');  % 6N perturbation on x, 
% load('data/processedData/ss4159_4168_d.mat', 'data'); % ones contain df and dx
% load('data/processedData/ss4169_4179_d.mat', 'data');     % Himanshu
Data = data;
% for each perturbed trial 
sudoK = nan(3,3,5);
pert_i = 2; % 2, 3, or 4
for fce_i = 1:3
    for dist_i = 1:3
        for trial_i = 1:5
            if isempty(Data{1,1,fce_i,dist_i,trial_i,pert_i})
                continue
            end
            % write a function on the df/dx
            trialtmp = Data{1,1,fce_i,dist_i,trial_i,pert_i}; 
            if (pert_i == 3)
                psudoStiffness = getDf_dev_dx(trialtmp,1,0,1);
            else
                psudoStiffness = getDf_dev_dx(trialtmp,1,1,0);
            end
            if psudoStiffness<0 
                psudoStiffness = nan;
            end
            sudoK(fce_i,dist_i,trial_i) = psudoStiffness;
                % sudoK(3,3,4) = sudoK(3,3,5); % temperary code for run
        end
        % find the pert time 
        % pert_t = 1:
        % corresponding to the force peak and displacement peak 
    end
end

% see if psudoK change with conditions 
sudoK_mean = mean(sudoK,3, 'omitnan');
sudoK_std = std(sudoK, 0, 3, 'omitnan');
fh_sudoK = figure(); 
hold on;

for fce_i = 1:3
%     plot(dist_list, sudoK_mean(fce_i,:), ...
%         'LineWidth', 3);
    errorbar(dist_list, sudoK_mean(fce_i,:), sudoK_std(fce_i,:),...
        'LineWidth', 3);
end
xlim([1 9]);
xlabel('positio (cm)'); 
ylabel('psudo K (N/m)');
legend('15N', '20N', '25N');

% Try the 2-way anova
sudoK_prac = sudoK(:,:,1:4);
sudoK_pracMark = zeros(3,3,4);
for fcei = 1:3
    for disti = 1:3
        sudoK_pracMark(fcei,disti,:) = (fcei-1)*3+disti; %integer from 1 to 9
    end
end

% try use a something represent the real data (to make sure in right anova2 format) 
sudoK_2d = [];
for fce_i = 1:3
%     sudoK_2d = [sudoK_2d; reshape(sudoK_pracMark(fce_i,:,:), 3, 4)'];
    sudoK_2d = [sudoK_2d; reshape(sudoK_prac(fce_i,:,:), 3, 4)'];
end
nmbtrials = 4;
[~,~,stats] = anova2(sudoK_2d,nmbtrials);


%% do the dF/dx during the pulse, pulse during release, save it in the data 

% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4129_4137.mat', 'data');
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4159_4168.mat', 'data');
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4169_4179.mat', 'data');
data = reshape(data(1,:,:,:,:,:), size(data, [2 3 4 5 6]));
data = reshape(data(1,:,:,:,:), size(data, [2 3 4 5]));pertT_num = size(data,4);
force_list = [15, 20, 25];
dist_list = [2.5, 5.0, 7.5];
pert_list = [0 -1 1 2]; % 0, no eprt; -1, pert before release; 1, pert after release; 2, hold 
color_arr = colormap('lines');
close all;
% sec = [-1 2];
sec = [-10 10];
% sec = [-0.4 0.8];
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
    dat_org = celltmp1{1,pi}.Fp(1,:);
%     plot(t(idx_interest),dat_org(idx_interest));
%     ylabel('command force (N)');
    pertIdx = find(abs(celltmp1{1,pi}.Fp(1,:)) > max(abs(celltmp1{1,pi}.Fp(1,:)) * 0.05));
    if (~isempty(pertIdx))
%         xline(t(pertIdx(1)));
    end
    
    for fce_i = 1:size(data,1)
        % plot the command f
        for pi = 3
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
            dat_mean_mat1 = [];
            dat_mean_mat2 = [];
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,1})
                    continue;
                end
                idx_release = find(celltmp1{ti,1}.ts == 5);
                t = celltmp1{ti,1}.t - celltmp1{ti,1}.t(idx_release(1));
                
                
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                % interpolate
                dat_org = celltmp1{ti,1}.f(1,idx_interest);
                dat_interest1 = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');

                dat_org = celltmp1{ti,1}.ox(1,idx_interest);
                dat_interest2 = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');


                ifplot = 0;
                if (ifplot)
                    clf;
                    hold on;
                    plot(t(idx_interest), dat_org, 'Marker', '.');
                    plot(pair_t, dat_interest1, 'Marker','.');
                    legend('org', 'itp');
                end
                %                 v_mean_mat = [v_mean_mat; celltmp1{ti,1}.v(2, idx_aftrelease)];
                dat_mean_mat1 = [dat_mean_mat1; dat_interest1];
                dat_mean_mat2 = [dat_mean_mat2; dat_interest2];
            end
            dat_mean1 = mean(dat_mean_mat1);
            dat_mean2 = mean(dat_mean_mat2);
            t_mean =  pair_t; % 500 data points
            
            % 2 plot the original one, non-perturbed
            for ti = 1:size(celltmp1,1)
                if isempty(celltmp1{ti,pi})
                    continue;
                end
                idx_release = find(celltmp1{ti,pi}.ts == 5);
                t = celltmp1{ti,pi}.t - celltmp1{ti,pi}.t(idx_release(1));
                pertIdx = find(abs(celltmp1{ti,pi}.Fp(1,:)) > max(abs(celltmp1{ti,pi}.Fp(1,:)) * 0.05));
                t_pert = [];
                if (~isempty(pertIdx))
                    t_pert = t(pertIdx(1));
                else
                    disp(['F' num2str(fce_i) 'D' num2str(dist_i) 'p' num2str(pi)]);
                end
                
                idx_interest = find(t>=sec(1) & t<=sec(2));
                dat_org = celltmp1{ti,pi}.f(1,idx_interest);
                dat_interest1 = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');
                dat_org = celltmp1{ti,pi}.ox(1,idx_interest);
                dat_interest2 = interp1(t(idx_interest), dat_org, pair_t, 'linear', 'extrap');

%                 plot(t_mean, dat_interest - dat_mean,  'color', color_arr(4+dist_i,:));
                clf;
                axh(1) = subplot(2,1,1); 
                plot(t_mean-t_pert, dat_interest1 - dat_mean1,  'color', color_arr(4+dist_i,:));
                title('df');
                axh(2) = subplot(2,1,2);
                plot(t_mean-t_pert, dat_interest2 - dat_mean2,  'color', color_arr(4+dist_i,:));
                linkaxes(axh, 'x');
                title('dx');

                dat_save.df = interp1(pair_t, dat_interest1 - dat_mean1, t, 'linear', 'extrap');
                dat_save.dx = interp1(pair_t, dat_interest2 - dat_mean2, t, 'linear', 'extrap');

                data{fce_i, dist_i,ti,pi}.df = dat_save.df;
                data{fce_i, dist_i,ti,pi}.dx = dat_save.dx;

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
Data(1,1,:,:,:,:) = data; 
data = Data; 
% save it back 
% save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4129_4137.mat', 'data')
% save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4159_4168_d.mat', 'data');
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4169_4179_d.mat', 'data');

%% % 
% plot arm configuration in one single trial 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4159_4168.mat', 'data');
fig_name = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4159_4168/'...
    'arm3d_demo_pert1.gif'];
fh = figure(); 
for fce_i = 3
    for dist_i = 3
        for pert_i = 2
            for trial_i = 1
                switch pert_i
                    case 1
                        timezone = [-0.5 1.5];
                    case 2
                        timezone = [-1.0 1.0];
                    case 3
                        timezone = [-0.5 1.5];
                    case 4
                        timezone = [0.0 2.0];
                end
                for epoci = 1:100%??? how to pair this with time?
                    t_near = timezone(1) + epoci*range(timezone)/100;
                    trialtmp = data{1,1,fce_i,dist_i,trial_i,pert_i};
                    % on plot the 3 markers position with 2 lines
                    [~,t_idx] = min(abs(trialtmp.t-t_near));
                    pt1 = reshape(trialtmp.ox(1:3,t_idx,1), 1, 3);
                    pt2 = reshape(trialtmp.ox(1:3,t_idx,2), 1, 3);
                    pt3 = reshape(trialtmp.ox(1:3,t_idx,3), 1, 3);
                    
                    % plot
                    clf; hold on;
                    plot3([pt1(1) pt2(1)], [pt1(2) pt2(2)], [pt1(3) pt2(3)], 'LineWidth', 3);
                    plot3([pt2(1) pt3(1)], [pt2(2) pt3(2)], [pt2(3) pt3(3)], 'LineWidth', 3);
                    xlim([-0.5 0]);
                    ylim([-0.8 -0.5]);

                    exportgraphics(gcf,fig_name,'Append',true);
                    
                end
            end
        end
    end
end