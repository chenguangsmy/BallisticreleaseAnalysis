% make packages of data

% data(subject, direction, distance, trial, perturbation)
% find max trials... 
% HOW TO MAKE SURE EVERY PARTS HAS 15 TRIALS?
data = cell(1,4,2, 15, 3);

% direction: front | back | left | right | 
% distance:  5     | 10   |
% only use 21 N 


% |sessions: | 9N*10cm  | 9N*5cm   | 21N*10cm | 21N*5cm  |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2462     | 2460     | 2467     | 2464     | 
% |back      | 2487     | 2484     | 2492     | 2489     |
% |left      | 2657     | 2655     | 2658     | 2661     |
% |right     | 2674     | 2672     | 2675     | 2676     |

% subject: chenguang
ss_num =    [2464, 2489, 2661, 2676;
            2467 2492 2568 2675];

% |sessions: | 15N*2.5cm| 15N*5cm  | 15N*10cm | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2818     | 2819     | 2820     | 2841     |
% |back      | 2825     | 2826     | 2827     | 2842     |
% |left      | 2829     | 2830     | 2831     | 2839     |
% |right     | 2832     | 2833     | 2834     | 2840     |
% perturbation type: 0-no pert, 1-pulse pert, 2-stoc pert
ss_num = [  2818    2819    2820    2841
            2825    2826    2827    2842
            2829    2830    2831    2839
            2832    2833    2834    2840];

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
save('data/processedData/ss2818_2842.mat', 'data')

%% 
% plot to check release 
dattmp = reshape(data(1,1,1,:,1), 15, 1);
figure(); 
axh(1) = subplot(2,1,1); hold on;
axh(2) = subplot(2,1,2); hold on;
linkaxes(axh, 'x');
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

% plot to check pert
dattmp = reshape(data(1,1,1,:,2), 15, 1);
figure(); hold on;
for i = 1:15
    if isempty(dattmp{i})
        continue;
    end
    % pert-time
    idx1 = find(dattmp{i}.Fp ~=0);
    time = dattmp{i}.t - dattmp{i}.t(idx1(1));
    % shift time
    
    % plot
    plot(time, dattmp{i}.x(2,:));
end  

%% move file 
srcdir1 = '/Volumes/rg2/data/KingKong/Raw/';
srcdir2 = '/Volumes/rg2/data/KingKong/Formatted/';
dstdir  = 'data/';
ss_num = [2839 2840 2841 2842];
for ss_i = 1:length(ss_num)
    fname1 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongFT%05d.csv'], ss_num(ss_i), ss_num(ss_i))
    fname2 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAM%05d.csv'], ss_num(ss_i), ss_num(ss_i))
    fname3 = sprintf([srcdir2 'KingKong.%05d.mat'], ss_num(ss_i))
    dstn1  = sprintf([dstdir 'KingKongFT%05d.csv'], ss_num(ss_i))
    dstn2  = sprintf([dstdir 'KingKongWAM%05d.csv'], ss_num(ss_i))
    dstn3  = sprintf([dstdir 'KingKong.%05d.mat'], ss_num(ss_i))
    copyfile(fname1, dstn1);
    copyfile(fname2, dstn2);
    copyfile(fname3, dstn3);
    %copyfile(dstn1, fname1);
    %copyfile(dstn2, fname2);
    %copyfile(dstn3, fname3);
end