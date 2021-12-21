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

%ss_num = [  3259 3258 3257 3261 ];
%ss_num = [3473 3471 3472 3474];
%ss_num = [3478 3476 3477 3479];    % hongwei with block visual feedback
%ss_num = [3480 3481 3482 3479];    % hongwei with visual feedback 
ss_num = [3473 3471 3472 3474];     % Chenguang with block visual feedback

%cpDatarg2(setdiff(ss_num, [3481 3479]));
for dir_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted(1);
        
        trials_num = size(celltmp,1);
    %    if trials_num>15
    %        data(1,dir_i,tar_i,:,:) = celltmp(1:15,:);
    %    else
            data(1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
     %   end
    end

    ss_tmp = SessionScan(ss_num(dir_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_ii = 1:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        %if trials_num>15
        %    data(1,dir_i,tar_ii,1:15,3) = celltmp(tar_ii,1:15,3);
        %else
            data(1,dir_i,tar_ii,1:trials_num,3) = celltmp(tar_ii,:,3);
        %end
    end
end
%save('data/processedData/ss3257_3261.mat', 'data')
%save('data/processedData/ss3471_3474.mat', 'data')
%save('data/processedData/ss3477_3479.mat', 'data')
%save('data/processedData/ss3480_3482.mat', 'data')
save('data/processedData/ss3471_3474.mat', 'data')


%% %%%% 2. Subject with multiple force levels %%%%%%%%%%
%% export a data from formatted requirements 
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

clear;
ss_num = [      3306    3307    3308    3317;
                3309    3310    3311    3318;
                3312    3313    3314    3319];
data = cell(1,4,3,3,15,3);
for fce_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(fce_i, tar_i));
        celltmp = ss_tmp.export_as_formatted(1);
        
        trials_num = size(celltmp,3);
        if trials_num>15
            data(1,1,fce_i,tar_i,:,:) = celltmp(1,1,1:15,:);
        else
            data(1,1,fce_i,tar_i,1:trials_num,:) = celltmp(1,1,:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(fce_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_i = 1:3 % as a session has 3 length
        trials_num = size(celltmp,3);
        if trials_num>15
            data(1,1,fce_i,tar_i,1:15,3) = celltmp(1,tar_i,1:15,3);
        else
            data(1,1,fce_i,tar_i,1:trials_num,3) = celltmp(1,tar_i,:,3);
        end
    end
end
save('data/processedData/ss3307_3314.mat', 'data')

%% export a data from formatted requirements 
% Date: 2021-11-04 USE SPRINGS TO DO THE SAME THING with ss3307_ss3314
% SPRING TESTING
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

% backward movements: (make James happier, do not change to order of
% subject)
% format like:
% |sessions:| K160     | K320      | K640      | K0     |
% | ------- | -------- | --------- | --------- | ------ |
% | 15N     | 3683     | 3688      | 3689      |    -   |
% | 20N     | 3684     | 3687      | 3690      |    -   |
% | 25N     | 3685     | 3686      | 3691      |    -   |
% | 0N      |     -    |    -      |    -      | 3682   |

% To make sure I only need to adjust the spring combinations equilibrium 
% position once, I combined the pulse and stoc perturbation with the same 
% springs and equilibirum positions in one session. 

% ss_num = [  3341        3335        3340
%             3343        3334        3338
%             3344        3336        3337];

% ss_num = [  3683        3688        3689
%             3684        3687        3690
%             3685        3686        3691];

ss_num = [  3725        3722        3689
            3684        3723        3690
            3685        3724        3691];

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
% save('data/processedData/ss3334_3344.mat', 'data')
save('data/processedData/ss3683_3691.mat', 'data')

%% Do the spring test for mutiple directions... (for x & y direction movement)
% use springs 320N/m, F0 = 25N (maximize the spring exerting force to
% enable under-damped movement after release.)  

% format like:
% |sessions:| K320,25N  |
% | ------- | --------- |
% | front   | 3763      |
% | back    | 3764      |
% | left    | 3761      |
% | right   | 3759      |  

clear;
ss_num = [      3763    3764    3761    3759];
data = cell(4,1,1,15,3); % direction-force-stiffness-trial-perturbation
for dir_i = 1:4
        ss_tmp = SessionScan(ss_num(dir_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(dir_i,1,1,:,:) = celltmp(1:15,:);
        else
            data(dir_i,1,1,1:trials_num,:) = celltmp(:,:);
        end


end
save('data/processedData/ss3759_3764.mat', 'data')



%% ss_num = [3345];
ss_num = [3682];
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
% save('data/processedData/ss3345.mat', 'data')
save('data/processedData/ss3682.mat', 'data')


%% Convert spring using the same order with the subjects

ss_num = [  3340        3335        3341
            3338        3334        3343
            3337        3336        3344];
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
save('data/processedData/ss3334_3344_alter.mat', 'data')

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

%% format like:
% still test with the 8N pulse perturbation and 2N stochastic perturbation.
% 8N pulse: 1) sigma = 0.04s, 2) sigma = 0.02s
% |sessions:        | K160     | K320      | 
% | -------         | -------- | --------- | 
% | sigma0.04       | 3466     | 3467      | 
% | sigma0.02       | 3465     | 3468      | ??? 3464? 3465???
ss_num = [3466, 3467; 3465, 3468];
cpDatarg2(ss_num);
for dir_i = 1:size(ss_num,1)
    for tar_i = 1:size(ss_num,2) % different springs
        ss_tmp = SessionScan(ss_num(dir_i, tar_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,1,dir_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,1,dir_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

end
save('data/processedData/ss3465_3468.mat', 'data')

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
% % % %%% Chenguang tested with the arm sling off...
% % % %front
% % % ss_num(1,:,:) =...
% % %         [   3354        3355        3353        3415
% % %             3359        3356        3358        3416
% % %             3361        3363        3364        3417 ];
% % % %back
% % % ss_num(2,:,:) =...
% % %          [  3372        3373        3371        3382     
% % %             3376        3374        3375        3413     
% % %             3377        3378        3379        3414 ];
% % % 
% % % % left
% % % ss_num(3,:,:) =...
% % %          [  3400        3401        3399        3408      
% % %             3404        3402        3403        3409      
% % %             3405        3406        3407        3410 ];
% % % %right
% % % ss_num(4,:,:) =...
% % %            [3385        3386        3384        3395     
% % %             3390        3387        3388        3411     
% % %             3391        3393        3394        3412 ];
%%% Chenguang tested with the arm sling on... 
ss_num(1,:,:) =...
        [   3487   	3486  	3488    3552
            3495   	3494  	3489    3550
            3507   	3506  	3505    3551];
%back
ss_num(2,:,:) =...
         [ 	3510     3509    3508   3553       
            3513     3512    3511   3554       
            3514     3516    3515   3555];

% left
ss_num(3,:,:) =...
         [	3518     3517    3519   3556       
            3522     3521    3520   3557       
            3523     3525    3535   3558];
%right
ss_num(4,:,:) =...
         [  3527     3526    3528   3559       
            3531     3530    3529   3562       
            3532     3534    3533   3563];

for dir_i = 1:size(ss_num,1)%1:size(ss_num,1)
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
% save('data/processedData/ss3486_3534.mat', 'data')
%save('data/processedData/ss3345.mat', 'data')

% % check the data
% trial_i = 1
% axh(1) = subplot(2,1,1);  plot(data{1,1,1,trial_i,1}.t, data{1,1,1,trial_i,1}.f);
% axh(2) = subplot(2,1,2);  plot(data{1,1,1,trial_i,1}.t, data{1,1,1,trial_i,1}.ts);
% linkaxes(axh, 'x');

%% export a data from formatted requirements 
%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PART OF THE CODE HERE! %%%%%%%%%%%%%%%%%%%%%
% Date: 2021-11-10 SUBJECT HIMANSHU TRYS THESE EXPERIMENTS, 
% With the same force level, 3 distance, and with different task settings.
% TASK SETTING  1. Allow overshoot and stop earlier; 
%               2. Do not allow the trials have >1cm overshoot, and trials
%                  of movement too long will fail (0.4s). 
% EXPORT FORMAT: will fill these two settings in different force levels. 
% format like:
% 1. ALLOW OVERSHOOT
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |15N       | 3433     | 3432     | 3431     | 3434     |
% 2. FORBID OVERSHOOT
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |15N       | 3438     | 3437     | 3436     | 3440     |

% export data defines as: 
% data = cell(n_subject, n_dir, n_fce, n_tar, n_trials, 3); 
% whereas 3 means: 1) release; 2) step pert; 3) stoc pert;

% Date: 2021-11-10 SUBJECT HONGWEI TRYS THESE EXPERIMENT
% With the same force level, 3 distance, and with different task settings.
% TASK SETTING  1. Subject cannot see the cursor during the movement; 
%               2. Subject is able to visualize the dursor during movement 

% 1. BLOCK visual feedback
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |15N       | 3478     | 3476     | 3477     | 3479     |
% 2. allow visual feedback
% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |15N       | 3480     | 3481     | 3482     | ----     |
% note: ss3480 hase just not enough trials

clear; clc; close all;
% % HIMANSHU
% ss_num = [  3433  3432	3431	3434;
%             3438	3437 	3436	3440];
% HONGWEI
% ss_num = [3478 3476 3477 3479];
% Chenguang with harder condition 
%ss_num = [3487 3486 3488 3415
%          3495 3494 3489 3416];

% % HIMANSHU AT HARDER CONDITION, REALLY GOOD LOOKING DATA 
% ss_num = [3492 3491 3493 3434];
% % HIMANSHU AT HARDER CONDITION, 3 FORCE LEVELS!!!
% ss_num = [  3499   3498  3500  3597
%             3492   3491  3497  3599
%             3502   3501  3503  3600];

% DELIN
%ss_num = [3581 3580 3582 3434];
%ss_num = [3581 3580 3584 3434];

% % MICHAEL WITH 20N force level. 
% 
% ss_num = [  3629    3630    3661    3632
%             3650    3648    3665    3651
%             3627    3664    3662    3628];

% %ADAM WITH 3 FORCE LEVELS;
% ss_num = [  3621    3620    3622    3667
%             3637    3635    3636    3638
%             3668    3670    3669    3672]; 
        
% MARCO WITH 3 FORCE LEVELS: 
ss_num = [  3643    3641    3642    3644
            3647    3645    3646    3654
            3656    3655    3657    3658];

% CHENGUANG TRY BACKWARD FULLY EXTENDED
% ss_num = [  3587    3586    3585    3596
%             3588    3590    3589    3595
%             3592    3591    3593    3594];
data = cell(1, 4, 3, 3, 15, 3);

idx_sf = 1; % only choose suceed trials 
for fce_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        fprintf('ss: %d', ss_num(fce_i, tar_i));
        ss_tmp = SessionScan(ss_num(fce_i, tar_i));
        ss_tmp.displayBlockCondition();
        celltmp = ss_tmp.export_as_formatted(1);
        
        trials_num = size(celltmp,3);

        if trials_num>15
            data(1,1,fce_i,tar_i,:,:) = celltmp(idx_sf,1,1:15,:);
        else
            data(1,1,fce_i,tar_i,1:trials_num,:) = celltmp(idx_sf,1,:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(fce_i, 4)); % stoc pert
    fprintf('ss: %d', ss_num(fce_i, 4));
    ss_tmp.displayBlockCondition();
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_i = 1:3 % as a session has 3 length
        trials_num = size(celltmp,3);
        data(1,1,fce_i,tar_i,1:5,3) = celltmp(idx_sf,tar_i,1:5,3);
%         if trials_num>15
%             data(1,1,fce_i,tar_i,1:15,3) = celltmp(idx_sf,tar_i,1:15,3);
%         else
%             data(1,1,fce_i,tar_i,1:trials_num,3) = celltmp(idx_sf,tar_i,:,3);
%         end
    end
end
%save('data/processedData/ss3431_3440.mat', 'data')
%save('data/processedData/ss3476_3479.mat', 'data')
%save('data/processedData/ss3486_3495.mat', 'data')
%save('data/processedData/ss3491_3493.mat', 'data')
%save('data/processedData/ss3580_3582.mat', 'data')
%save('data/processedData/ss3585_3596.mat', 'data')
% save('data/processedData/ss3491_3503.mat', 'data')    % HIMANSHU'S NEW TASK CONDITION
% save('data/processedData/ss3615_3618.mat', 'data')      % MICHAEL TRY "TOLERANT" TASK CONDITION
% save('data/processedData/ss3620_3629.mat', 'data')       % ADAM TRY "TOLERANT" TASK CONDITION
save('data/processedData/ss3641_3644.mat', 'data')       % MARCO TRY "TOLERANT" TASK CONDITION


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try Andy's idea on no release only pulse 
clear; clc; close all;

% % CHENGUANG TRIED AT NO RELEASE CONDITIONS. EACH TRIAL HAS 15 TRIALS EACH
% THAT: 1. NO PERTURB; 2. PERT+RELEASE; 3. NO RELEASE

% |sessions: | 2.5cm    | 5cm      | 7.5cm    | StocPert |
% | --------| --------  | -------- | -------- | -------- |
% |15N      | 3604      | 3612     | 3613     |   -      |
% |20N      | 3607      | 3606     | 3608     |   -      |
% |25N      | 3611      | 3610     | 3609     |   -      |

ss_num = [  3604        3612        3613          3600
            3607        3606        3608          3600
            3611        3610        3609          3600];

data = cell(1, 4, 3, 3, 15, 4);     % 1. subj
                                    % 2. directions
                                    % 3. Force levels
                                    % 4. Distance levels
                                    % 5. Trial count
                                    % 6. Task conditions:   1. no pulse
                                    %                       2. have pulse
                                    %                       3. stoc
                                    %                       4. release

for fce_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(fce_i, tar_i));
        celltmp = ss_tmp.export_as_formatted_5(1);
            if (fce_i == 3 && tar_i == 1) 
                celltmp = ss_tmp.export_as_formatted_5_failedTrials(1);
            end
        % the release are not so important, and these two sessions have
        % XM.module issue (that do not try next trial).
        if sum(ss_tmp.ssnum == [3603 3602 3612 3613])
            celltmp(15,1) = celltmp(1,1);
            if sum(ss_tmp.ssnum == [3602 3613])
                celltmp(14,1) = celltmp(1,1);
            end
        end
        
        trials_num = size(celltmp,1);
        if trials_num>15
            data(1,1,fce_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,1,fce_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(fce_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_i = 1:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        if trials_num>15
            data(1,1,fce_i,tar_i,1:15,3) = celltmp(tar_i,1:15,3);
        else
            data(1,1,fce_i,tar_i,1:trials_num,3) = celltmp(tar_i,:,3);
        end
    end
end
save('data/processedData/ss3602_3611_adp.mat', 'data');

%% also tidy up himanshu's data in our data 
% subj: Chenguang, James, Himanshu
clear; 
data1 = cell(4, 4, 3, 3, 15, 3);
load('data/processedData/ss3353_3417.mat', 'data');     %chenguang
data1(1,:,:,:,:,:) = data(1,:,:,:,:,:);
clear data;
load('data/processedData/ss3307_3314_6D.mat', 'data');  % James
data1(2,:,:,:,:,:) = data(1,:,:,:,:,:);
clear data; 
load('data/processedData/ss3431_3440.mat', 'data');     % Himanshu
%data1(3,:,1,:,:,:) = data(1,:,1,:,:,:); % the data have too much submovements
%data1(4,:,1,:,:,:) = data(1,:,2,:,:,:);
data1(3,:,1,:,:,:) = data(1,:,2,:,:,:);
load('data/processedData/ss3476_3479.mat', 'data');     %HONGWEI
data1(4,:,1,:,:,:) = data(1,:,1,:,:,:); 
clear data; 

data = data1;
%save('data/processedData/prelimData_3subj.mat', 'data')
%save('data/processedData/prelimData_3subj_restict.mat', 'data')
save('data/processedData/prelimData_4subj_fine.mat', 'data')

%% tidy up 6 subjects' data in our data
clear;
fnames = {
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3486_3534.mat'; % Chenguang testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3307_3314.mat'; % James testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3491_3503.mat'; % Himanshu testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3615_3618.mat'; % Micheal testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3620_3629.mat'; % Adam testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3641_3644.mat'; % Marco testing
    };
data1 = cell(6, 4, 3, 3, 15, 3);
for idx_subj = 1:6
    clear data;
    load(fnames{idx_subj});
    data1(idx_subj,:,:,:,:,:) = data(1,:,:,:,:,:);
end
data = data1;
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/prelimData_6subj_fine.mat', 'data', '-v7.3')

%% spring + 6 subjects
clear;
fnames = {
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/ProcessedData/ss3334_3344_alter.mat'; % Springs
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3486_3534.mat'; % Chenguang testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3307_3314.mat'; % James testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3491_3503.mat'; % Himanshu testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3615_3618.mat'; % Micheal testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3620_3629.mat'; % Adam testing
    '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3641_3644.mat'; % Marco testing
    };
data1 = cell(7, 4, 3, 3, 15, 3);
for idx_subj = 1:7
    clear data;
    load(fnames{idx_subj});
    if idx_subj==1
        data1(idx_subj,1,:,:,:,:) = data(1,:,:,:,:,:);
    else
        data1(idx_subj,:,:,:,:,:) = data(1,:,:,:,:,:);
    end
end
data = data1;
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/prelimData_S_6subj_fine.mat', 'data', '-v7.3')

%% tidy up Himanshu's data as a comparation 
clear; 
data1 = cell(3, 4, 3, 3, 15, 3);
load('data/processedData/ss3431_3440.mat', 'data');     % Himanshu
data1(1,:,1,:,:,:) = data(1,:,1,:,:,:);     % the data have too much submovements
data1(2,:,1,:,:,:) = data(1,:,2,:,:,:);     % limit submovemnets by displacement
load('data/processedData/ss3491_3493.mat', 'data');     % Himanshu with change condition
data1(3,:,1,:,:,:) = data(1,:,1,:,:,:);     % limit submovements by velocity
clear data; 

data = data1;
%save('data/processedData/prelimData_3subj.mat', 'data')
%save('data/processedData/prelimData_3subj_restict.mat', 'data')
save('data/processedData/prelimData_subjHA_compare.mat', 'data')

%% Chenguang test with the step perturbation (long-step )

ss_num = [  3566 3565 3564  3434
            3567 3569 3568  3434
            3573 3572 3574  3434];
data = cell(1, 4, 3, 3, 15, 3);

for fce_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(fce_i, tar_i));
        celltmp = ss_tmp.export_as_formatted_4(1);
        
        trials_num = size(celltmp,1);
        celltmp1 = celltmp(:,[1,4,3]);
        celltmp = celltmp1;
        if trials_num>15
            data(1,1,fce_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,1,fce_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(fce_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_i = 1:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        if trials_num>15
            data(1,1,fce_i,tar_i,1:15,3) = celltmp(tar_i,1:15,3);
        else
            data(1,1,fce_i,tar_i,1:trials_num,3) = celltmp(tar_i,:,3);
        end
    end
end
save('data/processedData/ss3564_3569.mat', 'data')

%% Chenguang test backward mvoement with big perturbation
ss_num = [  3577 3579 3578 3434];
data = cell(1, 4, 3, 3, 15, 3);

for fce_i = 1:size(ss_num,1)
    for tar_i = 1:3 % step perts
        ss_tmp = SessionScan(ss_num(fce_i, tar_i));
        celltmp = ss_tmp.export_as_formatted_4(1);
        
        trials_num = size(celltmp,1);
        celltmp1 = celltmp(:,[1,4,3]);
        celltmp = celltmp1;
        if trials_num>15
            data(1,1,fce_i,tar_i,:,:) = celltmp(1:15,:);
        else
            data(1,1,fce_i,tar_i,1:trials_num,:) = celltmp(:,:);
        end
    end

    ss_tmp = SessionScan(ss_num(fce_i, 4)); % stoc pert
    celltmp = ss_tmp.export_as_formatted(1);
    for tar_i = 1:3 % as a session has 3 length
        trials_num = size(celltmp,2);
        if trials_num>15
            data(1,1,fce_i,tar_i,1:15,3) = celltmp(tar_i,1:15,3);
        else
            data(1,1,fce_i,tar_i,1:trials_num,3) = celltmp(tar_i,:,3);
        end
    end
end
save('data/processedData/ss3575_3578.mat', 'data')


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
ss_num = [3110 3111 3115; 3108 3113 3114]; % Kr=100; the combination of ±2N sigma = 0.16s 
% ss_num = [3116 3119 3120; 3117 3118 3122]; % Kr=100; the combination of ±2N sigma = 0.04s 
% ss_num = [3124 3125 3128; 3123 3126 3127]; % Kr=100; the combination of ±4N sigma = 0.04s

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
ss_num1 = [3147 3150 3151; 3148 3149 3154]; % Kr=300; the combination of ±2N sigma = 0.04s
ss_num2 = [3142 3143 3146; 3141 3144 3145]; % Kr=300; the combination of ±2N sigma = 0.16s
ss_num3 = [3129 3132 3133; 3130 3131 3134]; % Kr=300; the combination of ±4N sigma = 0.04s
ss_num4 = [3136 3137 3140; 3135 3138 3139]; % Kr=300; the combination of ±4N sigma = 0.16s
ss_num5 = [3156 3157 3160; 3155 3158 3159]; % Kr=300; the combination of ±6N sigma = 0.04s
ss_num6 = [3163 3166 3167; 3164 3165 3168]; % Kr=300; the combination of ±6N sigma = 0.16s
% Kr=300; the combination of ±8N sigma = 0.02s

ss_num_mat = {  ss_num1     ss_num3     ss_num5;
                ss_num2     ss_num4     ss_num6}';
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
                plot(pert_mag_C{ci}.Fp(2,:));
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

%% export the parameter scanned data, with more parameter selection
% The subject cg was performing all these data on Oct 17 2021 + Nov 10 2021 
ss_num1 = [3147 3150 3151; 3148 3149 3154]; % Kr=300; the combination of ±2N sigma = 0.04s
ss_num2 = [3142 3143 3146; 3141 3144 3145]; % Kr=300; the combination of ±2N sigma = 0.16s
ss_num3 = [3129 3132 3133; 3130 3131 3134]; % Kr=300; the combination of ±4N sigma = 0.04s
ss_num4 = [3136 3137 3140; 3135 3138 3139]; % Kr=300; the combination of ±4N sigma = 0.16s
ss_num5 = [3156 3157 3160; 3155 3158 3159]; % Kr=300; the combination of ±6N sigma = 0.04s
ss_num6 = [3163 3166 3167; 3164 3165 3168]; % Kr=300; the combination of ±6N sigma = 0.16s
ss_num7 = [3445 3447 3448; 3451 3450 3449]; % Kr=300; the combination of ±8N sigma = 0.04s
ss_num8 = [3458 3459 3460; 3463 3462 3461]; % Kr=300; the combination of ±8N sigma = 0.02s
ss_num_blk = [];
% Kr=300; the combination of ±8N sigma = 0.02s

ss_num_mat = {  ss_num_blk  ss_num_blk  ss_num_blk  ss_num8;        % ... sigma = 0.02s
                ss_num1     ss_num3     ss_num5     ss_num7;        % ... sigma = 0.04s 
                ss_num2     ss_num4     ss_num6     ss_num_blk}';   % ... sigma = 0.16s
mags_idx = [2 4 6 8];
% Data-type: pert_mag * pert_duration * target_dist * pert_pos/neg
r = size(ss_num_mat, 1);
c = size(ss_num_mat, 1);
Data = cell(r, c, 3, 2);
for matri = 1:r
    for matci = 1:c
        ss_num = ss_num_mat{matri, matci};
        if isempty(ss_num) % in case of no such data
            continue;
        end
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
                plot(pert_mag_C{ci}.Fp(2,:));
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
save('data/processedData/PertParamSelection_gaussian_cg_4by3.mat', 'Data')

%% export the method validation data with Spring tests. 
% After the experiment, James want to use spring data to test the new model
% with deconvolve. 

% with all of them have the [-4 +4 -6 +6] gaussian pert, the sigma of
% gaussian is 0.04s
%ss_num = [3180 3182 3184]; 
% | ss_num  | descriptions                              | 
% | ------- | ----------------------------------------- |
% | 3192    | peak±4N,±6N, sigma = 0.04s Ks = 160N/m    |
% | 3191    | peak±4N,±6N, sigma = 0.04s Ks = 320N/m    |
% | 3189    | peak±4N,±6N, sigma = 0.04s Ks = 640N/m    |
% | 3201    | no spring hooked                          |

% Data-type: pert_mag * pert_duration * spring_stiffness * pert_pos/neg
Data = cell(2, 1, 3, 2);
ss_num = [3192 3189 3191 3201]; 
mags_idx = [4 6];

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
%save('data/processedData/SpringGaussian3.mat', 'Data')

%% export the method validation data with Spring tests. 
% After the experiment, James want to use spring data to test the new model
% with deconvolve, try to use the 8N to test 

% | ss_num  | descriptions                              | 
% | ------- | ----------------------------------------- |
% | 3466    | peak -8N, sigma = 0.04s Ks = 160N/m       |
% | 3467    | peak -8N, sigma = 0.04s Ks = 320N/m       |
% | 3464    | peak -8N, sigma = 0.02s Ks = 160N/m       |3465?
% | 3468    | peak -8N, sigma = 0.02s Ks = 320N/m       |

% Data-type: pert_mag * pert_duration * spring_stiffness * pert_pos/neg
Data = cell(1, 2, 2, 1);
ss_num = [  3466    3467;
            3465    3468]; 
mags_idx = [ 8 ];

for dur_i = 1:2
    for stf_i = 1:2
        ss_tmp = SessionScan(ss_num(dur_i, stf_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
 
        mag_list = [-mags_idx(1)];
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
        Data(1, dur_i, stf_i, :) = data;
    end
end
%save('data/processedData/SpringGaussian_8N.mat', 'Data')
%save('data/processedData/SpringGaussian3.mat', 'Data')
%% 

%% export the method validation data with Spring tests. 
% JAMES said he want the 'standard format'. I'm guessing it is the subject
% format with 6D

% | ss_num  | descriptions                              | 
% | ------- | ----------------------------------------- |
% | 3466    | peak -8N, sigma = 0.04s Ks = 160N/m       |
% | 3467    | peak -8N, sigma = 0.04s Ks = 320N/m       |
% | 3464    | peak -8N, sigma = 0.02s Ks = 160N/m       |3465?
% | 3468    | peak -8N, sigma = 0.02s Ks = 320N/m       |

% Data-type: pert_mag * pert_duration * spring_stiffness * pert_pos/neg
Data = cell(2, 1, 2, 2, 15, 3);
ss_num = [  3466    3467;
            3465    3468]; 
mags_idx = [ 8 ];

for dur_i = 1:2
    for stf_i = 1:2
        ss_tmp = SessionScan(ss_num(dur_i, stf_i));
        celltmp = ss_tmp.export_as_formatted_hybridss(1);
 
        mag_list = [-mags_idx(1)];
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
        if length(pert_mag_C) > 15 % could be more than 15 trials 
            Data(1, 1, dur_i, stf_i, 1:15,1) = celltmp(1:15, 1);
            Data(1, 1, dur_i, stf_i, 1:15,2) = celltmp(1:15, 2);
            Data(1, 1, dur_i, stf_i, 1:15,3) = celltmp(1:15, 3);
        end
        
            
        %data{1,1} = pert_mag_C(pert_mags == mag_list(1))% pos and neg
        %Data(1, dur_i, stf_i, :) = data;
    end
end
%Data(2,1,1,:,:,3) = Data(1,1,1,:,:,3);
save('data/processedData/SpringGaussian_8N_format.mat', 'Data')
%save('data/processedData/SpringGaussian3.mat', 'Data')

%% export data with only 3N stoc pert
Data = cell(1, 1, 1, 1, 15, 3);
ss_num = [  3504]; 
mags_idx = [ 0 ];

ss_tmp = SessionScan(ss_num(1));
celltmp = ss_tmp.export_as_formatted_hybridss(1);
Data(1, 1, 1, 1, 1:9,3) = celltmp(1:9, 3);
save('data/processedData/SpringStoc_3N.mat', 'Data')
%save('data/processedData/SpringGaussian3.mat', 'Data')



%% plot to check release 
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