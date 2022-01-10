% trialCountInTrainingTesting
% count how many trials are needed in the training and testing. 
% So that I can have an idea on how much trial is needed both in testing
% and training. 
% Compare among the subjects did in the first direction testing 

% subject list: Chenguang, James, Himanshu, Michael, Marco, Adam

% Will compare among sessions did for each subject! 
clear all; 

ss_num(1,:,:) =...                          % Chenguang
        [   3487   	3486  	3488    3552
            3495   	3494  	3489    3550
            3507   	3506  	3505    3551];
        
ss_num(2,:,:) =...                          % James
        [   3306    3307    3308    3317;
            3309    3310    3311    3318;
            3312    3313    3314    3319];
    
ss_num(3,:,:) =...                          % Himanshu
        [   3499   3498  3500  3597
            3492   3491  3497  3599
            3502   3501  3503  3600];
% ss_num = [  3433  3432	3431        3434;  % himanshu in training
%             3438	3437 	3436,3493	3440];

ss_num(4,:,:) =...
        [   3629    3630    3661    3632    % Michael
            3650    3648    3665    3651
            3627    3664    3662    3628]; 
% ss_num = [ -      -       3631         -;        % Michael in training 
%           3617    3615    3616,3649    3618
%           -       3624    3626    -   ];
        
ss_num(5,:,:) =...
        [   3643    3641    3642    3644    % Marco
            3647    3645    3646    3654
            3656    3655    3657    3658];
        % Marco in training
        
ss_num(6,:,:) =...
        [   3621    3620    3622    3667    % Adam
            3637    3635    3636    3638
            3668    3670    3669    3672]; 
        % Adam in training

        
%% 
% 1. Sum up all the sessions and find the:  1. sucess trials, 
%                                           2. total trials, 
%                                           3. sucessful rate 
% for each subject 
sucess_rate = zeros(3,6,3,4);

for subj_i = 1:6
    for ri = 1:3
        for ci = 1:4
            ss_num(subj_i,ri,ci)
            sstmp = SessionScan(ss_num(subj_i,ri,ci));
            sucess_rate(1,subj_i,ri,ci) = sum([sstmp.trials.outcome]);
            sucess_rate(2,subj_i,ri,ci) = sstmp.trials_num;
            sucess_rate(3,subj_i,ri,ci) = sucess_rate(1,subj_i,ri,ci)/sucess_rate(2,subj_i,ri,ci);
        end
    end
end

save('trialCountInTrainingTesting_sucessRate.mat', 'sucess_rate');

%% 
sucess_Rate = zeros(3,6);
for subj_i = 1:6
    for ri = 1:3
        for ci = 1:3
            sucess_Rate(1,subj_i) = sucess_Rate(1,subj_i) + sucess_rate(1,subj_i,ri,ci);
            sucess_Rate(2,subj_i) = sucess_Rate(2,subj_i) + sucess_rate(2,subj_i,ri,ci);
        end
    end
    sucess_Rate(3,subj_i) = sucess_Rate(1,subj_i) / sucess_Rate(2,subj_i);
end


% 2. Find how many trials are taken in the subjects previously "trained" 
