%% par1. export data into formatted 

clear; clc; close all;

% % HIMANSHU 
% ss_num = [  3499   3498  3500  3597
%             3492   3491  3497  3599
%             3502   3501  3503  3600];

% % MICHAEL 
% 
% ss_num = [  3629    3630    3661    3632
%             3650    3648    3665    3651
%             3627    3664    3662    3628];

% %ADAM 
% ss_num = [  3621    3620    3622    3667
%             3637    3635    3636    3638
%             3668    3670    3669    3672]; 
        
% MARCO WITH 3 FORCE LEVELS: 
ss_num = [  3643    3641    3642    3644
            3647    3645    3646    3654
            3656    3655    3657    3658];


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

% save('data/processedData/ss3491_3503.mat', 'data')    % HIMANSHU'S NEW TASK CONDITION
% save('data/processedData/ss3615_3618.mat', 'data')      % MICHAEL TRY "TOLERANT" TASK CONDITION
% save('data/processedData/ss3620_3629.mat', 'data')       % ADAM TRY "TOLERANT" TASK CONDITION
save('data/processedData/ss3641_3644.mat', 'data')       % MARCO TRY "TOLERANT" TASK CONDITION

%% part2. check the data, if Force data is continuous

% seems bad trials here: 
% [s3648, t8], [s3641, t54], [s3647, t15], [s3645, t26]

ss_bad_list = [3650, 3641, 3647, 3645];
trial_bad_list = {[8], [54], [15], [26]};

for ss_i = 1:length(ss_bad_list)
    sstmp = SessionScan(ss_bad_list(ss_i));
    ss_bad_list(ss_i)
    celltmp = sstmp.export_as_formatted(1);
    
    fh = sstmp.plotTaskEndpointForce;
    axh(1) = sstmp.plotAddTrialMark(fh.Children(1));
    axh(2) = sstmp.plotAddTrialMark(fh.Children(3));
    linkaxes(axh, 'x')
end
