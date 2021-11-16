% This scripts make packages of data out from the raw data 

% data format contains: 
%   1. Subject with one force level (not chromnically);
%       data(subject, direction, distance, trial, perturbation)
%   2. Subject with multiple force level: 
%       data(subject, direction, force, distance, trial, perturbation)
%   3. Parameter selection data, subject/springs with different pert magnitude.
%       
%
%   


%% %%%%%%%%%%%%% DATA FILES ARE UNDER FOLLOWING: %%%%%%%%%%
%% %%%% 1. Subject with one force level %%%%%%%%%%
%% experiment first test with different force levels and displacements. 
% This data comes from assorted trials and errors. (range from 24xx to 26xx)
% |sessions: | 9N*7.5cm | 9N*5cm   | 21N*7.5cm| 21N*5cm  |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2462     | 2460     | 2467     | 2464     | 
% |back      | 2487     | 2484     | 2492     | 2489     |
% |sessions: | 9N*10cm  | 9N*5cm   | 21N*10cm | 21N*5cm  |
% |left      | 2657     | 2655     | 2658     | 2661     |
% |right     | 2674     | 2672     | 2675     | 2676     |
data = cell(1,4,2, 15, 3);


%% experiment test with fixed force (15N) but different targets
%  data contains both pulse and stochastic perturbations
%  These datas are not aligned with the blackrock, thus less use-able.

% |sessions: | 15N*2.5cm| 15N*5cm  | 15N*10cm | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2858     | 2819     | 2820     | 2841     |
% |back      | 2865     | 2826     | 2827     | 2842     |
% |left      | 2829     | 2830     | 2831     | 2839     |
% |right     | 2832     | 2833     | 2834     | 2840     |


% ss_num = 2949;

% perturbation type: 0-no pert, 1-pulse pert, 2-stoc pert
% % ss_num = [  2858    2819    2820    2841
% %             2868    2826    2827    2842
% %             2829    2830    2831    2839
% %             2832    2833    2834    2840];
ss_num = [2985 2985 2987];%2962;
for dir_i = 1%:size(ss_num,1)
   for tar_i = 1:3 % step perts
       
        %dir_i = 4
        %tar_i = 2
        
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted(1);
       
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,dir_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
   end

        
    %ss_tmp = SessionScan(ss_num(dir_i, 4)); % stoc pert
    ss_tmp = SessionScan(2950); % stoc pert, 5cm 12N
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_ii = 1%:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        if trials_num>15
            data(1,dir_i,tar_ii,1:15,3) = celltmp(tar_ii,1:15,3);
        else
            data(1,dir_i,tar_ii,1:trials_num,3) = celltmp(tar_ii,:,3);
        end
    end
end
save('data/processedData/ss2985_2987.mat', 'data')

%% export a data from formatted requirements 
% James test during his first visit, having data with WAM part, squre wave 
% 15N. 
% This data did not align with the blackrock. Not good use.
% format like:
% |sessions: | 15N*2.5cm| 15N*5cm  | 15N*10cm | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2872     | 2873     | 2874     | 2876     |

ss_num = [  2872 2873 2874 2876 ];

for dir_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,dir_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(dir_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_ii = 1:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        if trials_num>15
            data(1,dir_i,tar_ii,1:15,3) = celltmp(tar_ii,1:15,3);
        else
            data(1,dir_i,tar_ii,1:trials_num,3) = celltmp(tar_ii,:,3);
        end
    end
end
save('data/processedData/ss2872_2876.mat', 'data')



%% export a data from formatted requirements:
% Date: 2021-11-01 Just before experiment tests, use one force level and
% different targets displacements.
% format like:
% |sessions: | 15N*2.5cm| 15N*5cm  | 15N*10cm | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 3259     | 3258     | 3257     | 3261     |

ss_num = [  3259 3258 3257 3261 ];

for dir_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,dir_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(dir_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_ii = 1:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        if trials_num>15
            data(1,dir_i,tar_ii,1:15,3) = celltmp(tar_ii,1:15,3);
        else
            data(1,dir_i,tar_ii,1:trials_num,3) = celltmp(tar_ii,:,3);
        end
    end
end
save('data/processedData/ss3257_3261.mat', 'data')


%% %%%% 2. Subject with multiple force levels %%%%%%%%%%
% export a data from formatted requirements 
% Date: 2021-11-03 SUBJECT JAMES TOOKS THESE EXPERIMENTS, trying different
% force level with 1 direction (the total last about 1h20min).
% format like:
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |15N       | 3306     | 3307     | 3308     | 3317     |
% |20N       | 3309     | 3310     | 3311     | 3318     |
% |25N       | 3312     | 3313     | 3314     | 3319     |
% export data defines as: 
% data = cell(n_subject, n_dir, n_fce, n_trials, 3); 
% whereas 3 means: 1) release; 2) step pert; 3) stoc pert;


ss_num = [  3306 3307 3308 3317;
            3309 3310 3311 3318;
            3312 3313 3314 3319];

for dir_i = 1:size(ss_num,1)
    for fce_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(dir_i, fce_i));
        celltmp = ss_tmp.export_as_formatted(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,dir_i,fce_i,:,:) = celltmp(1:15,:);
        else
            data(1,dir_i,fce_i,1:trials_num,:) = celltmp(:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(dir_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for fce_ii = 1:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        if trials_num>15
            data(1,dir_i,fce_ii,1:15,3) = celltmp(fce_ii,1:15,3);
        else
            data(1,dir_i,fce_ii,1:trials_num,3) = celltmp(fce_ii,:,3);
        end
    end
end
save('data/processedData/ss3307_3314.mat', 'data')

%% export a data from formatted requirements 
% Date: 2021-11-04 USE SPRINGS TO DO THE SAME THING with ss3307_ss3314
% make sure that:   1. The force exerted on robot is the same with subject's level; 
%                   2. The robot before-release position is on the same level with subject;
% Thus, I fixed the xr0 (robot nominal position) depend on the force, and
% varies the xs0 (spring nominal position) and Ks (springs) 

% format like:
% |sessions:| K160     | K320      | K640      | K0     |
% | ------- | -------- | --------- | --------- | ------ |
% | 15N     | 3341     | 3335      | 3340      |    -   |
% | 20N     | 3343     | 3334      | 3338      |    -   |
% | 25N     | 3344     | 3336      | 3337      |    -   |
% | 0N      |     -    |    -      |    -      | 3345   |

% To make sure I only need to adjust the spring combinations equilibrium 
% position once, I combined the pulse and stoc perturbation with the same 
% springs and equilibirum positions in one session. 

ss_num = [  3341        3335        3340
            3343        3334        3338
            3344        3336        3337];
for dir_i = 1:size(ss_num,1)
    for tar_i = 1:size(ss_num,2) % step perts
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,dir_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

end
save('data/processedData/ss3334_3344.mat', 'data')


ss_num = [3345];
for dir_i = 1:size(ss_num,1)
    for tar_i = 1:size(ss_num,2) % step perts
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,dir_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

end
save('data/processedData/ss3345.mat', 'data')

%% format like:
% |sessions:| K160     | K320      | 
% | ------- | -------- | --------- | 
% | 15N     | 3419     | 3418      | 
ss_num = [3418, 3419];
for dir_i = 1:size(ss_num,1)
    for tar_i = 1:size(ss_num,2) % different springs
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,dir_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

end
save('data/processedData/ss3318_3319.mat', 'data')

%% export a data from formatted requirements 
% Date: 2021-11-06 CHENGUNG DO MULTIPLE DIRECTIONS
% Have some data for better analysis. The perturbation is -6N pulse.  
% The data contains direction, force, target. 
% TO MAKE SURE ITS ABLE TO HAVE SUCH DATA, CHANGE THE FORMAT INTO: 
% data = cell(n_subj, n_dir, n_frc, n_dist, n_trials, n_pert);
% format like:
% FRONT
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | --------| --------  | -------- | -------- | -------- |
% |15N      | 3354      | 3355     | 3353     | 3415     |
% |20N      | 3359      | 3356     | 3358     | 3416     |
% |25N      | 3361      | 3363     | 3364     | 3417     |
% BACK
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | --------| --------  | -------- | -------- | -------- |
% |15N      | 3372      | 3373     | 3371     | 3382     |
% |20N      | 3376      | 3374     | 3375     | 3413     |
% |25N      | 3377      | 3378     | 3379     | 3414     |
% LEFT
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | --------| --------  | -------- | -------- | -------- |
% |15N      | 3400      | 3401     | 3399     | 3408     |
% |20N      | 3404      | 3402     | 3403     | 3409     |
% |25N      | 3405      | 3406     | 3407     | 3410     |
% RIGHT
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | --------| --------  | -------- | -------- | -------- |
% |15N      | 3385      | 3386     | 3384     | 3395     |
% |20N      | 3390      | 3387     | 3388     | 3411     |
% |25N      | 3391      | 3393     | 3394     | 3412     |
clear; close all; clc; 
ss_num = zeros(3, 3, 4);
data = cell(1, 4, 3, 3, 15, 3);
%front
ss_num(1,:,:) =...
        [   3354        3355        3353        3415
            3359        3356        3358        3416
            3361        3363        3364        3417 ];
%back
ss_num(2,:,:) =...
         [  3372        3373        3371        3382     
            3376        3374        3375        3413     
            3377        3378        3379        3414 ];

% left
ss_num(3,:,:) =...
         [  3400        3401        3399        3408      
            3404        3402        3403        3409      
            3405        3406        3407        3410 ];
%right
ss_num(4,:,:) =...
           [3385        3386        3384        3395     
            3390        3387        3388        3411     
            3391        3393        3394        3412 ];

for dir_i = 1:size(ss_num,1)
    for fce_i = 1:size(ss_num,2) 
        for tar_i = 1:3 % step perts
            ss_tmp = SessionScan(ss_num(dir_i, fce_i, tar_i));
            celltmp = ss_tmp.export_as_formatted(1);

            trials_num = size(celltmp,1);
            if trials_num>15
                data(1,dir_i, fce_i, tar_i,:,:) = celltmp(1:15,:);
            else
                data(1,dir_i, fce_i, tar_i,1:trials_num,:) = celltmp(:,:);
            end
        end
    end
        
    for fce_i = 1:size(ss_num,2) 
        ss_tmp = SessionScan(ss_num(dir_i, fce_i, 4));
        celltmp = ss_tmp.export_as_formatted(1);
        for tar_ii = 1:3 % as a session has 3 length
            trials_num = size(celltmp,2);
            if trials_num>15
                data(1,dir_i, fce_i, tar_ii,1:15,3) = celltmp(tar_ii,1:15,3);
            else
                data(1,dir_i, fce_i, tar_ii,1:trials_num,3) = celltmp(tar_ii,:,3);
            end
        end
    end
end
%save('data/processedData/ss3353_3417.mat', 'data')
%save('data/processedData/ss3345.mat', 'data')

% check the data
trial_i = 1
axh(1) = subplot(2,1,1);  plot(data{1,1,1,trial_i,1}.t, data{1,1,1,trial_i,1}.f);
axh(2) = subplot(2,1,2);  plot(data{1,1,1,trial_i,1}.t, data{1,1,1,trial_i,1}.ts);
linkaxes(axh, 'x');

%% %%%% 3. Subject/Springs with multiple perturb parameters %%%%%%%%%%

% ... This part 
% trying to find the optimal parameters for the perturbation. 
% | sessions | 2999                     | 3000                      | 3001                  | 
% | -------- | ------------------------ | ------------------------- | --------------------- |
% | tsigma   | 0.04                     | 0.02                      | 0.01                  |
% |magnitude | 0.1 0.2 0.5 1 2 4 6 8    | 0.1 0.2 0.5 1 2 4 6 8     | 0.1 0.2 0.5 1 2 4 6 8 |

% WANTED FORMAT:
% data = cell(3, 8);                    % size(tsigma) - by - size(magnitude)
%ss_num = [  2999 3000 3001 ];          % mag_list = [0.1 0.2 0.5 1 2 4 6 8];
%ss_num = [3003 3006 3005];             % mag_list = [1 2 4 6 8]
%ss_num = [3062 3063 3064];             % mag_list = -[1 2]N, with 16N force offset
%ss_num = [3087]                        % ----------------test the strange force reason -------------
%ss_num = [3085]                        % ----------------test the strange force reason -------------
%ss_num = 3101;                         % mag = 8N, perturbation with Kx = 700N/m
ss_num = 3176;                          %
%ss_num = 3103;
%ss_num = 3104;
%ss_num = [3047];                       % mag = [+2N], sigma = 0.04s 10cm 15N ballistic release
% mag_list = [ 1 2 4 6 8];
mag_list = [-4 4];
%mag_list = [20];
%mag_list = [-2 2];
data = cell(length(ss_num), length(mag_list)); 
for ss_i = 1:length(ss_num)
        ss_tmp = SessionScan(ss_num(ss_i));
        celltmp = ss_tmp.export_as_formatted(1);
        % find idx in this cell tmp
        pert_mag_C = celltmp(:,2); 
        pert_mags = [];
        for ci = 1:length(pert_mag_C)
            if (isempty(pert_mag_C{ci}))
                pert_mags = [pert_mags -1];
                continue;
            end
            %pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))];
            pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))*setdiff(unique(sign(pert_mag_C{ci}.Fp(2,:))), 0)];
        end
        % make sure the pert_mags are categorical:
        pert_mags(pert_mags<1) = round(pert_mags(pert_mags<1)*10)/10;
        pert_mags(pert_mags>=1) = round(pert_mags(pert_mags>=1));
        % get the index of each magnitude
        tridx = cell(1, length(mag_list)); % trials in each session;
        for mag_idx = 1:length(mag_list)
            magtmp = mag_list(mag_idx);
            tridx{mag_idx} = find(pert_mags == magtmp);
        end
        % package
        for mag_idx = 1:length(mag_list)
            data{ss_i,mag_idx} = pert_mag_C(tridx{mag_idx});
        end
        
end
%save('data/processedData/ss2999_3001.mat', 'data')
%save('data/processedData/ss3003_3005.mat', 'data')
%save('data/processedData/ss3062_3064.mat', 'data')
%save('data/processedData/ss3085.mat', 'data')
%save('data/processedData/ss3103.mat', 'data')
%save('data/processedData/ss3104.mat', 'data')
save('data/processedData/ss3176.mat', 'data')

%% export according to the magnitude of perturbation 
% | target      | 10cm          | 5cm             |  2.5cm                |
% | ----------- | ------------- | --------------- | --------------------- |
% | +8N pulse   |  3090         | 3092            | 3095                  | 
% | -8N pulse   |  3089         | 3093            | 3094                  |
% sigma is always 0.04s

% wanted format:
% data = cell(2, 3); % size(tsigma) - by - size(magnitude)
ss_num = [3090 3092 3095; 
          3089 3093 3094]

mag_list = [-8 8];
data = cell(size(ss_num,2), length(mag_list)); 
for ss_i = 1:length(ss_num(:))
    close all;
        ss_tmp = SessionScan(ss_num(ss_i));
        celltmp = ss_tmp.export_as_formatted(1);
        % find idx in this cell tmp
        pert_mag_C = celltmp(:,2); 
        pert_mags = [];
        for ci = 1:length(pert_mag_C)
            if (isempty(pert_mag_C{ci}))
                pert_mags = [pert_mags -1];
                continue;
            end
            %pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))];
            pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))*setdiff(unique(sign(pert_mag_C{ci}.Fp(2,:))), 0)];
        end
        % make sure the pert_mags are categorical:
        pert_mags(pert_mags<1) = round(pert_mags(pert_mags<1)*10)/10;
        pert_mags(pert_mags>=1) = round(pert_mags(pert_mags>=1));
        % package            
        data{floor((ss_i+1)/2),find(pert_mags(1) == mag_list)} = pert_mag_C;

end
save('data/processedData/ss3089_3094.mat', 'data')

%% 
%% export according to the magnitude of perturbation, try only +-2N
% | target      | 10cm          | 5cm             |  2.5cm                |
% | ----------- | ------------- |  -------------- |  -------------------- |
% | +2N pulse   |  3110         | 3111            | 3115                  | 
% | -2N pulse   |  3108         | 3113            | 3114                  |
% sigma is always 0.16s (which means the whole pulse last about 1s

% | target      | 10cm          | 5cm             |  2.5cm                |
% | ----------- | ------------- |  -------------- |  -------------------- |
% | +2N pulse   |  3116         | 3119            | 3120                  | 
% | -2N pulse   |  3117         | 3118            | 3122                  |
% sigma is always 0.04s (which means the whole pulse last about 1s


% wanted format:
% data = cell(2, 3); % size(tsigma) - by - size(magnitude)
%%%%%%%%%%%%%%%% Where Kr=100N/m %%%%%%%%%%%%%%%%%
ss_num = [3110 3111 3115; 3108 3113 3114]; % Kr=100; the combination of �2N sigma = 0.16s 
% ss_num = [3116 3119 3120; 3117 3118 3122]; % Kr=100; the combination of �2N sigma = 0.04s 
% ss_num = [3124 3125 3128; 3123 3126 3127]; % Kr=100; the combination of �4N sigma = 0.04s

mag_list = [-2 2];
%mag_list = [-4 4];
%mag_list = [-6 6];
data = cell(size(ss_num,2), length(mag_list)); 
for ss_i = 1:length(ss_num(:))
    close all;
        ss_tmp = SessionScan(ss_num(ss_i));
        celltmp = ss_tmp.export_as_formatted(1);
        % find idx in this cell tmp
        pert_mag_C = celltmp(:,2); 
        pert_mags = [];
        for ci = 1:length(pert_mag_C)
            if (isempty(pert_mag_C{ci}))
                pert_mags = [pert_mags -1];
                continue;
            end
            %pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))];
            pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))*setdiff(unique(sign(pert_mag_C{ci}.Fp(2,:))), 0)];
        end
        % make sure the pert_mags are categorical:
        pert_mags(pert_mags<1) = round(pert_mags(pert_mags<1)*10)/10;
        pert_mags(pert_mags>=1) = round(pert_mags(pert_mags>=1));
        % package            
        data{floor((ss_i+1)/2),find(pert_mags(1) == mag_list)} = pert_mag_C;

end

save('data/processedData/ss3108_3114.mat', 'data')
%save('data/processedData/ss3116_3122.mat', 'data')
%save('data/processedData/ss3123_3128.mat', 'data')
%save('data/processedData/ss3129_3134.mat', 'data')
%save('data/processedData/ss3147_3154.mat', 'data')

%% export the parameter scanned data. 
% The subject cg was performing all these data on Oct 17 2021. 
ss_num1 = [3147 3150 3151; 3148 3149 3154]; % Kr=300; the combination of �2N sigma = 0.04s
ss_num2 = [3142 3143 3146; 3141 3144 3145]; % Kr=300; the combination of �2N sigma = 0.16s
ss_num3 = [3129 3132 3133; 3130 3131 3134]; % Kr=300; the combination of �4N sigma = 0.04s
ss_num4 = [3136 3137 3140; 3135 3138 3139]; % Kr=300; the combination of �4N sigma = 0.16s
ss_num5 = [3156 3157 3160; 3155 3158 3159]; % Kr=300; the combination of �6N sigma = 0.04s
ss_num6 = [3163 3166 3167; 3164 3165 3168]; % Kr=300; the combination of �6N sigma = 0.16s

ss_num_mat = {ss_num1  ss_num3  ss_num5; ss_num2 ss_num4 ss_num6}'
mags_idx = [2 4 6];
% Data-type: pert_mag * pert_duration * target_dist * pert_pos/neg
Data = cell(3, 2, 3, 2);
for matri = 1:3
    for matci = 1:2
        ss_num = ss_num_mat{matri, matci};
        fname = sprintf('ss%d_%d.mat', min(ss_num(:)), max(ss_num(:)));
        mag_list = [-mags_idx(matri), mags_idx(matri)];
        data = cell(size(ss_num,2), length(mag_list));
        for ss_i = 1:length(ss_num(:))
            close all;
            ss_tmp = SessionScan(ss_num(ss_i));
            celltmp = ss_tmp.export_as_formatted(1);
            % find idx in this cell tmp
            pert_mag_C = celltmp(:,2);
            pert_mags = [];
            for ci = 1:length(pert_mag_C)
                if (isempty(pert_mag_C{ci}))
                    pert_mags = [pert_mags -1];
                    continue;
                end
                %pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))];
                pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))*setdiff(unique(sign(pert_mag_C{ci}.Fp(2,:))), 0)];
            end
            % make sure the pert_mags are categorical:
            pert_mags(pert_mags<1) = round(pert_mags(pert_mags<1)*10)/10;
            pert_mags(pert_mags>=1) = round(pert_mags(pert_mags>=1));
            % package
            data{floor((ss_i+1)/2),find(pert_mags(1) == mag_list)} = pert_mag_C%
        end
        Data(matri, matci, :, :) = data;
    end
end
save('data/processedData/PertParamSelection_gaussian_cg.mat', 'Data')

%% export the method validation data with Spring tests. 
% After the experiment, James want to use spring data to test the new model
% with deconvolve. 

% with all of them have the [-4 +4 -6 +6] gaussian pert, the sigma of
% gaussian is 0.04s
%ss_num = [3180 3182 3184]; 
ss_num = [3192 3189 3191 3201]; 

mags_idx = [4 6];
% Data-type: pert_mag * pert_duration * spring_stiffness * pert_pos/neg
Data = cell(2, 1, 3, 2);
for matri = 1:4
    ss_tmp = SessionScan(ss_num(matri));
    celltmp = ss_tmp.export_as_formatted(1);
    for mag_i = 1:2
        fname = sprintf('ss%d_%d.mat', min(ss_num(:)), max(ss_num(:)));
        mag_list = [-mags_idx(mag_i), mags_idx(mag_i)];
        data = cell(1, length(mag_list));
            % find idx in this cell tmp
        pert_mag_C = celltmp(:,2);
        pert_mags = [];
        for ci = 1:length(pert_mag_C)
            if (isempty(pert_mag_C{ci}))
                pert_mags = [pert_mags -1];
                continue;
            end
            %pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))];
            pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))*setdiff(unique(sign(pert_mag_C{ci}.Fp(2,:))), 0)];
        end
            % make sure the pert_mags are categorical:
        pert_mags(pert_mags<1) = round(pert_mags(pert_mags<1)*10)/10;
        pert_mags(pert_mags>=1) = round(pert_mags(pert_mags>=1));
        % package
        data{1,1} = pert_mag_C(pert_mags == mag_list(1))% pos and neg
        data{1,2} = pert_mag_C(pert_mags == mag_list(2))
        Data(mag_i, 1, matri, :) = data;
    end
end
save('data/processedData/SpringGaussian3.mat', 'Data')
%% 
% plot to check release 
dir_i = 1;
tar_i = 1;
sub_i = 1;
trials_num = 14;
dattmp = reshape(data(sub_i,dir_i,tar_i,:,1), trials_num, 1);
figure(); 
axh(1) = subplot(2,1,1); hold on;
axh(2) = subplot(2,1,2); hold on;
linkaxes(axh, 'x');
for i = 1:trials_num
    if isempty(dattmp{i})
        continue;
    end
    idx1 = find(dattmp{i}.mvst ~=0);
    time = dattmp{i}.t - dattmp{i}.t(idx1(1));
    
    subplot(axh(1));
    plot(data{1,1,1,i,1}.x(2,:));
%    plot(dattmp{i}.t, dattmp{i}.x(2,:));
    %plot(time, dattmp{i}.x(2,:));
    subplot(axh(2));
    plot(data{1,1,1,i,1}.f(2,:));
%    plot(dattmp{i}.t, dattmp{i}.f(2,:));
    %plot(time, dattmp{i}.f(2,:));
end

% plot to check pert
for i = 1:15
    if isempty(dattmp{i})
        continue;
    end
    idx1 = find(dattmp{i}.mvst ~=0);
    time = dattmp{i}.t - dattmp{i}.t(idx1(1));
    
    subplot(axh(1));
%    plot(dattmp{i}.t, dattmp{i}.x(2,:));
    plot(time, dattmp{i}.x(2,:));
    subplot(axh(2));
%    plot(dattmp{i}.t, dattmp{i}.f(2,:));
    plot(time, dattmp{i}.f(2,:));
end

% plot to check step-pert
dattmp = reshape(data(sub_i,2,1,:,2), 15, 1);
figure(); hold on;
for i = 1:15
    if isempty(dattmp{i})
        continue;
    end
    % pert-time
    idx1 = find(dattmp{i}.Fp(2,:) ~=0);
    time = dattmp{i}.t - dattmp{i}.t(idx1(1));
    % shift time
    
    % plot
    plot(time, dattmp{i}.x(2,:));
    
end  


%% export the optotrak data with Spring tests. 
% Andy wants to see if the 'hystesis effect' still exists in the new data. 

% with all of them have the [-4 +4 -6 +6] gaussian pert, the sigma of
% gaussian is 0.04s
% The perturbations are intermediate, rather than Jammed together. 
%ss_num = [3180 3182 3184]; 
%ss_num = [3216 3213 3218]; % 3213 has the problem of only one type of pert (-6N)
ss_num = [3216 3218];

mags_idx = [4 6];
% Data-type: pert_mag * pert_duration * spring_stiffness * pert_pos/neg
%Data = cell(2, 1, 3, 2);
Data = cell(2, 1, length(ss_num), 2);
for matri = 1:length(ss_num)
    ss_tmp = SessionScan(ss_num(matri));
    celltmp = ss_tmp.export_as_formatted(1);
    for mag_i = 1:2
        fname = sprintf('ss%d_%d.mat', min(ss_num(:)), max(ss_num(:)));
        mag_list = [-mags_idx(mag_i), mags_idx(mag_i)];
        data = cell(1, length(mag_list));
            % find idx in this cell tmp
        pert_mag_C = celltmp(:,2);
        pert_mags = [];
        for ci = 1:length(pert_mag_C)
            if (isempty(pert_mag_C{ci}))
                pert_mags = [pert_mags -1];
                continue;
            end
            %pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))];
            pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))*setdiff(unique(sign(pert_mag_C{ci}.Fp(2,:))), 0)]
        end
            % make sure the pert_mags are categorical:
        pert_mags(pert_mags<1) = round(pert_mags(pert_mags<1)*10)/10;
        pert_mags(pert_mags>=1) = round(pert_mags(pert_mags>=1));
        % package
        data{1,1} = pert_mag_C(pert_mags == mag_list(1))% pos and neg
        data{1,2} = pert_mag_C(pert_mags == mag_list(2))
        Data(mag_i, 1, matri, :) = data;
    end
end
save('data/processedData/SpringGaussianFixXs0_1.mat', 'Data')

%% %%%%%%%%%%%%%% THIS SECTION CONTAINS THE HANDY-DANDY PLOT CODES %%%%%%%%%%%%%%%%%%%

%% brief try to debug the export
celltmp = ss3185.export_as_formatted(1)
%%
%  quick look
dattmp = ss2833;
             figure;
             ax1 = subplot(3,1,1); plot(ss2833.wam_t, ss2833.wam.jp(:,2));
             ax2 = subplot(3,1,2); plot(ss2833.force_t, ss2833.force_h(2,:));
             ax3 = subplot(3,1,3); plot(ss2833.wam_t, ss2833.wam.cf(:,2));
             linkaxes([ax1,ax2, ax3],'x');
             
%% 

dataTmp = celltmp{6,2};
            figure;
            ax1 = subplot(3,1,1); plot(dataTmp.t, dataTmp.x(2,:)); 
            ylabel('x');
            ax2 = subplot(3,1,2); plot(dataTmp.t, dataTmp.f(2,:));
            ylabel('f');
            ax3 = subplot(3,1,3); plot(dataTmp.t, dataTmp.Fp(2,:));
            ylabel('Fp');
            linkaxes([ax1,ax2,ax3],'x');
            
%%
trial = 11;
f = data{1,1,1,trial,2}.f(2,:);
x = data{1,1,1,trial,2}.x(2,:);
Fp = data{1,1,1,trial,2}.Fp(2,:);
v = data{1,1,1,trial,2}.v(2,:);

%just plot
figure();
axh(1) = subplot(3,1,1); plot(f);
axh(2) = subplot(3,1,2); plot(x); 
axh(3) = subplot(3,1,3); plot(Fp);
linkaxes(axh, 'x');

%% 
figure()
axh(1) = subplot(2,1,1); plot(x);
axh(2) = subplot(2,1,2); plot(smooth(diff(diff(smooth(x)))));
linkaxes(axh, 'x');

%% 
obj = ss2975;
figure();
axh(1) = subplot(3,1,1); plot(obj.wam_t,obj.wam.cf);
axh(2) = subplot(3,1,2); plot(obj.wam_t,obj.wam.tp(:,2)); 
axh(3) = subplot(3,1,3); plot(obj.wam_t,obj.wam.state); 

linkaxes(axh, 'x');