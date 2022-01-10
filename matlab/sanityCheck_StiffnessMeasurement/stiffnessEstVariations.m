% test for the Stiffness measurement variations. 
% Look through the data of pulling back, the stiffness variates a lot (the
% pulse method have big variance). 

% data:  'ss3585_3596.mat' : ----- Chenguang pull back with hand fully extended
%

% Address following questions: 
% 1. The difference between perturb stiffness and release stiffness; 
% 2. The pulse method (perturb stiffness) is much higher than the release
% ones. 
% 3. The pulse method has very high variance 

% The following are explainations (using plot) to the previous questions
% Q1. The difference between perturb stiffness and release stiffness:
%   S1. See if there is difference in the release stiffness in release case
%   and have perturb case.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. 
% Data: 'ss3486_3534.mat': dir1234: fblr
% data I have: 
%               Chenguang-front     ?   'ss3353_3417'   ? 1
%               James-front         ?   'ss3307_3314_6D'? 2
%               Hongwei-front       ?   'ss3476_3479'   ? 3
%               Himanshu-front      ?   'ss3491_3493'   ? 4 
%               Chenguang-left      ?   'ss3353_3417'   ? 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% general variables
file_names = {'ss3353_3417', 'ss3307_3314_6D', 'ss3476_3479', 'ss3491_3503', 'ss3353_3417' };
dir_idx    = [1 1 1 1 3];
force_max = [3 3 1 3 3];
force = [15 20 25]; %N
dist =  [2.5 5 7.5];
names = {'Chenguang-F', 'James-F', 'Hongwei-F', 'Himanshu-F', 'Chenguang-L'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1.1 measurement  of release-case stiffness 
clear all; close all; clc;
% should use this one for cg: ss3564_3569
file_names = {'ss3353_3417', 'ss3307_3314_6D', 'ss3476_3479', 'ss3491_3503', 'ss3353_3417' };
dir_idx    = [1 1 1 1 3];
force_max = [3 3 1 3 3];
% save the data in a .mat 
rel_stiffness_mat = cell(5,2); 
pert_stiffness_mat  = cell(5,2); 
for dexSubj = 1:5
    load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/' file_names{dexSubj} '.mat']);
    data_human = reshape(data(1,dir_idx(dexSubj),:,:,:,:),1,3,3,15,3); % front here 
    clear data
    dexSubject = 1; % [Chenguang]
    dexForce = 1:force_max(dexSubj); % [15N, 20N, 25N]
    dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
    try 
        depMeasures_human1 = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
    catch 
        continue % go to the next subject 
    end
    rel_stiffness_mat{dexSubj,1} = depMeasures_human1.k_hat_release;
    pert_stiffness_mat{dexSubj,1} = depMeasures_human1.k_hat_pulse;
    % 1.2 measurement of the perturb-case stiffness 
    data_human2 = data_human; 
    data_human2(1,:,:,:,1) = data_human(1,:,:,:,2); 
    dexSubject = 1; % [Chenguang]
    dexForce = 1:force_max(dexSubj); % [15N, 20N, 25N]
    dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
    depMeasures_human2 = crossConditionAnalysis(data_human2, dexSubject, dexForce, dexDistance,'human');
    rel_stiffness_mat{dexSubj,2} = depMeasures_human2.k_hat_release;
    pert_stiffness_mat{dexSubj,1} = depMeasures_human1.k_hat_pulse;
    close all
end
%% 1.2 epxort data into datatmp/ReleaseStiffnesssCompare 
save('datatmp/ReleaseStiffness.mat', 'pert_stiffness_mat', '-append');
%save('datatmp/ReleaseStiffness.mat', 'rel_stiffness_mat', '-append');
%% 1.3 plot data: 
load('datatmp/ReleaseStiffness.mat', 'rel_stiffness_mat');

% each data points 
force = [15, 20, 25];
dist  = [2.5, 5.0, 7.5]; 
stiffness_val_theoratical = repmat(force', 1, 3) ./repmat(dist, 3, 1) * 100
for dexSubj = 1%1:5
    figure();
    %subplot(1,5,dex_subj); % subject
    depMeasures_human1.k_hat_release = rel_stiffness_mat{dexSubj,1}; 
    depMeasures_human2.k_hat_release = rel_stiffness_mat{dexSubj,2};
    if isempty(depMeasures_human1.k_hat_release)
        continue;
    end
    for dex_f = 1:force_max(dexSubj)
        for dex_d = 1:3
            axh(dex_f, dex_d) = subplot(3,3,(dex_f-1)*3+dex_d); hold on; grid on;
            trials_n1 = size(depMeasures_human1.k_hat_release,4);
            trials_n2 = size(depMeasures_human2.k_hat_release,4);
            plot(1:trials_n1, reshape(depMeasures_human1.k_hat_release(1,dex_f,dex_d,:), 1, trials_n1), 'x'); % release stiffness in no pert trials
            plot(1:trials_n2, reshape(depMeasures_human2.k_hat_release(1,dex_f,dex_d,:), 1, trials_n2), 'o'); % release stiffness in pert trials
             plot(1:15, ones(1,15)*stiffness_val_theoratical(dex_f, dex_d));
            if (dex_f == 3 && dex_d == 3)
                legend('no pert trials', 'perturbed trials', 'theoratical value');
            end
            title(['F=' num2str(force(dex_f)) 'N, dist=' num2str(dist(dex_d)) 'cm']);
            xlabel('trial'); ylabel('stiffness (N/m)');
            
        end
    end
    linkaxes(axh(:), 'y')
    %sgtitle(names{dexSubj});
    sgtitle('release measurement on have pert/no pert conditions');
end

% summarize figure
% [15, 20, 25]N
% [2.5, 5. 7.5]cm
%% 
names = {'Chenguang-F', 'James-F', 'Himanshu-F', 'Hongwei-F', 'Chenguang-L'};

col_i = colormap('lines');
dex_subj_avail = [1 2 4 5];
for dexSubj = 1:4
    %axh1(dex_subj) = subplot(1,5,dex_subj); % subject 
    axh1(dexSubj) = subplot(1,3,dexSubj); % subject 
    title(names{dex_subj_avail(dexSubj)});
    hold on;
    depMeasures_human1.k_hat_release = rel_stiffness_mat{dex_subj_avail(dexSubj),1};
    depMeasures_human2.k_hat_release = rel_stiffness_mat{dex_subj_avail(dexSubj),2};
    if isempty(depMeasures_human1.k_hat_release)
        continue;
    end
    k_hat_release1_mean = zeros(3,3);
    k_hat_release1_std  = zeros(3,3);
    k_hat_release2_mean = zeros(3,3);
    k_hat_release2_std  = zeros(3,3);
    for dex_f = 1:3
        for dex_d = 1:3
            k_hat_release1_mean(dex_f, dex_d) = nanmean(reshape(depMeasures_human1.k_hat_release(1,dex_f,dex_d,:), 1, 15));
            k_hat_release1_std(dex_f, dex_d)  = nanstd(reshape(depMeasures_human1.k_hat_release(1,dex_f,dex_d,:), 1, 15));
            k_hat_release2_mean(dex_f, dex_d) = nanmean(reshape(depMeasures_human2.k_hat_release(1,dex_f,dex_d,:), 1, 15));
            k_hat_release2_std(dex_f, dex_d)  = nanstd(reshape(depMeasures_human2.k_hat_release(1,dex_f,dex_d,:), 1, 15));
        end
        %plot(k_hat_release1_mean(dex_f,:),'color', col_i(dex_f, :));
        %plot(k_hat_release2_mean(dex_f,:),'--','color', col_i(dex_f, :));
        errorbar(k_hat_release1_mean(dex_f,:), 0.001*k_hat_release1_std(dex_f, :),...
            'Marker', 'x', 'MarkerSize', 5, ...
            'LineWidth', 2, ...
            'color', col_i(dex_f, :));
        errorbar(k_hat_release2_mean(dex_f,:), 0.001*k_hat_release2_std(dex_f, :),...
            'Marker', 'o', 'MarkerSize', 5, ...
            'LineWidth', 2, ...
            'color', 1-(1-col_i(dex_f+7, :))/2);
    end
    xlim([0 4]);
end
linkaxes(axh1(:), 'y')
ylabel('release stiffness N/m');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% 2. Stiffness measurement across time. 
% 2.1 Get Delta_time 
clear all; close all; clc;
file_names = {'ss3353_3417', 'ss3307_3314_6D', 'ss3476_3479', 'ss3491_3503', 'ss3353_3417' };
dir_idx    = [1 1 1 1 3];
% save the data in a .mat 
Delta_t_mat = cell(5,1); 
for dexSubj = 1:5
    load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/' file_names{dexSubj} '.mat']);
    data_human = reshape(data(1,dir_idx(dexSubj),:,:,:,:),1,3,3,15,3); % front here 
    Delta_t = zeros(3,3,15); % force, distance, trials
    clear data
    dexSubject = 1; % [Chenguang]
    for dexForce = 1:3 % [15N, 20N, 25N]
        for dexDistance = 1:3 % [2.5cm, 5cm, 7.5cm]
            for dexTrial = 1:15
                trial_tmp = data_human{1,dexForce,dexDistance,dexTrial,2};
                if isempty(trial_tmp)
                    continue;
                end
                
                pert_idx = find(trial_tmp.Fp(2,:) ~=0);
                pert_time = trial_tmp.t(pert_idx(1)); 
                rel_idx = find(trial_tmp.ts == 5); 
                rel_time = trial_tmp.t(rel_idx(1)); 
                Delta_t(dexForce, dexDistance, dexTrial) = rel_time - pert_time;
                
                clf; 
                axh1(1) = subplot(3,1,1); hold on;
                plot(trial_tmp.t, trial_tmp.Fp);
                axh1(2) = subplot(3,1,2); hold on;
                plot(trial_tmp.t, trial_tmp.ts);
                axh1(3) = subplot(3,1,3); hold on;
                plot(trial_tmp.t, trial_tmp.x(2,:));
                plot(pert_time, trial_tmp.x(2,pert_idx(1)),'o', 'MarkerSize', 5);
                plot(rel_time, trial_tmp.x(2,rel_idx(1)),'x', 'MarkerSize', 5);
                linkaxes(axh1, 'x');
            end
        end
    end
    Delta_t_mat{dexSubj} = Delta_t;
end
save('datatmp/ReleaseStiffness.mat', 'Delta_t_mat', '-append');

%% 2.2 Pair Delta_time with the actual stiffness 
% 2.2.1. plot out one condition of one subject 
clear; 
load('datatmp/ReleaseStiffness.mat', 'Delta_t_mat', 'rel_stiffness_mat');
names = {'Chenguang-F', 'James-F', 'Hongwei-F', 'Himanshu-F', 'Chenguang-L'};
force = [15, 20, 25];
dist  = [2.5, 5.0, 7.5]; 
force_max = [3 3 1 3 3];
% rel_stiffness_mat(1,:) = {un_pert_releaseStf, pert_releaseStf}
% Delta_t_mat(1,:) =  t_release - t_pert
for dexSubj = [1 2 3 4 5]
    figure();
    
    for dex_f = 1:force_max(dexSubj)
        for dex_d = 1:3
            axh2 = subplot(3,3,(dex_f-1)*3+dex_d);
            
            t = reshape(Delta_t_mat{dexSubj}(dex_f,dex_d,:), 1, 15);
            k_hat_release = reshape(rel_stiffness_mat{dexSubj,2}(1,dex_f,dex_d,:), 1, 15);
            plot(t, k_hat_release, 'o', 'MarkerSize', 5);
            refline
            title(['F=' num2str(force(dex_f)) 'N, dist=' num2str(dist(dex_d)) 'cm']);
            xlabel('\Delta t'); ylabel('release stiffness');
        end
    end
    suptitle(names{dexSubj});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% 3. TRY ANDY'S IDEA THAT GIVE TRIALS HAVE ONLY PERTURBATION BUT NOT
% RELEASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%file_names = {'ss3602_3611' };
file_names = {'ss3602_3611_adp' };
% save the data in a .mat 
pert_stiffness_mat  = cell(1,2); 

    load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/' file_names{1} '.mat']);
    % the trials that have both perturbation and release 
    data_pr = data(:,:,:,:,:,[1:3]);
    % the trials have only perturbation and no release 
    data_pn = data(:,:,:,:,:,[1,4,3]);
    
    data_human_pr = reshape(data_pr(1,1,:,:,:,:),1,3,3,15,3); % front here 
    data_human_pn = reshape(data_pn(1,1,:,:,:,:),1,3,3,15,3); % front here 
    clear data_pr data_pn data
    dexSubject = 1; % [Chenguang]
    dexForce = 1:3; % [15N, 20N, 25N]
    dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
    try 
        depMeasures_human_pr = crossConditionAnalysis(data_human_pr, dexSubject, dexForce, dexDistance,'human');
        depMeasures_human_pn = crossConditionAnalysis(data_human_pn, dexSubject, dexForce, dexDistance,'human');
        
    catch 

    end
    
    % save the data 
    %save('datatmp/PerturbStiffness.mat', 'depMeasures_human_pr', 'depMeasures_human_pn', '-append');

%% 3.2 plot out the measured stiffness difference
    % load the data 
    %load('datatmp/PerturbStiffness.mat', 'depMeasures_human_pr', 'depMeasures_human_pn');
    
    
    dexForce = 1:3; % [15N, 20N, 25N]
    dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
    pert_stiffness_mat{1,1} = depMeasures_human_pr.k_hat_pulse; % have release
    pert_stiffness_mat{1,2} = depMeasures_human_pn.k_hat_pulse; % no release

    stiffness_val_theoratical = repmat(force', 1, 3) ./repmat(dist, 3, 1) * 100
    figure(); 
    for dex_f = 1:3
        for dex_d = 1:3
            axh2(dex_f,dex_d) = subplot(3,3,(dex_f-1)*3+dex_d); hold on; grid on;
            
            k_hat_pert1 = reshape(pert_stiffness_mat{1}(1,dex_f,dex_d,:), 1, 15);
            if (dex_f == 3 && dex_d == 1) % for adp, failure trials
                k_hat_pert1 = reshape(pert_stiffness_mat{1}(1,dex_f,dex_d,1:10), 1, 10);
            end
            k_hat_pert2 = reshape(pert_stiffness_mat{2}(1,dex_f,dex_d,:), 1, 15);
            plot(k_hat_pert1, 'o', 'MarkerSize', 5); % have release
            plot(k_hat_pert2, 'x', 'MarkerSize', 5); % no release
            plot(1:15, ones(1,15)*stiffness_val_theoratical(dex_f, dex_d));
            title(['F=' num2str(force(dex_f)) 'N, dist=' num2str(dist(dex_d)) 'cm']);
            xlabel('trial seq'); ylabel('pert stiffness'); 
            if (dex_f == 3 && dex_d == 3)
                legend('have release', 'no release', 'theoratical value');
            end
        end
    end
    linkaxes(axh2(:), 'y');
    suptitle('perturb stiffness in release/no release trials');
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try release measurement with different data length
%% 1.1 measurement  of release-case stiffness 
clear all; close all; clc;
% should use this one for cg: ss3564_3569
file_names = {'ss3353_3417'};
dir_idx    = [1];
force_max = [3];
% save the data in a .mat 
rel_stiffness_mat = cell(1,3);      % 3 durations  
pert_stiffness_mat  = cell(5,2); 
dexSubj = 1;
for release_timedur = 1:4 % 300, 200, 100, 50
    load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/' file_names{dexSubj} '.mat']);
    data_human = reshape(data(1,dir_idx(dexSubj),:,:,:,:),1,3,3,15,3); % front here 
    clear data
    dexSubject = 1; % [Chenguang]
    dexForce = 1:force_max(dexSubj); % [15N, 20N, 25N]
    dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
    try 
        depMeasures_human1 = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
    catch 
        continue % go to the next subject 
    end
    rel_stiffness_mat{release_timedur} = depMeasures_human1.k_hat_release;
    %pert_stiffness_mat{dexSubj,1} = depMeasures_human1.k_hat_pulse;
%    close all
end

