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

%% move file 
srcdir1 = '/Volumes/rg2/data/KingKong/Raw/';
srcdir2 = '/Volumes/rg2/data/KingKong/Formatted/';
srcdir3 = '/Volumes/rg2/data/KingKong/Intermediate/';

dstdir  = 'data/';
dstdir2 = 'data/Intermediate/'
% ss_num = [ 2858    2819    2820    2841
%             2868    2826    2827    2842
%             2829    2830    2831    2839
%             2832    2833    2834    2840 ];
ss_num = [3006]
ss_num = ss_num(:);
for ss_i = 1:length(ss_num)
    fname1 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongFT%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname2 = sprintf([srcdir1 'KingKong.DK.%05d/KingKongWAM%05d.csv'], ss_num(ss_i), ss_num(ss_i));
    fname3 = sprintf([srcdir2 'KingKong.%05d.mat'], ss_num(ss_i));
    fname4 = sprintf([srcdir3 'KingKong.%05d.mat'], ss_num(ss_i));
    dstn1  = sprintf([dstdir 'KingKongFT%05d.csv'], ss_num(ss_i));
    dstn2  = sprintf([dstdir 'KingKongWAM%05d.csv'], ss_num(ss_i));
    dstn3  = sprintf([dstdir 'KingKong.%05d.mat'], ss_num(ss_i));
    dstn4  = sprintf([dstdir2 'KingKong.%05d.mat'], ss_num(ss_i));
    copyfile(fname1, dstn1);
    copyfile(fname2, dstn2);
    copyfile(fname3, dstn3);
    copyfile(fname4, dstn4);
    %copyfile(dstn1, fname1);
    %copyfile(dstn2, fname2);
    %copyfile(dstn3, fname3);
end

%%
% |sessions: | 15N*2.5cm| 15N*5cm  | 15N*10cm | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2818     | 2819     | 2820     | 2841     |
% |back      | 2825     | 2826     | 2827     | 2842     |
% |left      | 2829     | 2830     | 2831     | 2839     |
% |right     | 2832     | 2833     | 2834     | 2840     |
% perturbation type: 0-no pert, 1-pulse pert, 2-stoc pert
% % ss_num = [  2858    2819    2820    2841
% %             2868    2826    2827    2842
% %             2829    2830    2831    2839
% %             2832    2833    2834    2840];
% ss_num = 2949;
ss_num = [2985 2985 2987]%2962;
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
    ss_tmp = SessionScan(2950); % stoc pert
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

%%
% format like:
% |sessions: | 15N*2.5cm| 15N*5cm  | 15N*10cm | StocPert |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2818     | 2819     | 2820     | 2841     |
% |back      | 2825     | 2826     | 2827     | 2842     |
% |left      | 2829     | 2830     | 2831     | 2839     |
% |right     | 2832     | 2833     | 2834     | 2840     |

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

%% export according to the magnitude of perturbation 
% | sessions | 2999                     | 3000                      | 3001                  | 
% | tsigma   | 0.04                     | 0.02                      | 0.01                  |
% |magnitude | 0.1 0.2 0.5 1 2 4 6 8    | 0.1 0.2 0.5 1 2 4 6 8     | 0.1 0.2 0.5 1 2 4 6 8 |

% wanted format:
% data = cell(3, 8); % size(tsigma) - by - size(magnitude)
%ss_num = [  2999 3000 3001 ];
%mag_list = [0.1 0.2 0.5 1 2 4 6 8];
ss_num = [3003 3006 3005];
mag_list = [ 1 2 4 6 8];
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
            pert_mags = [pert_mags max(abs(pert_mag_C{ci}.Fp(2,:)))];
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

%%
%  quick look
dattmp = ss2833;
             figure;
             ax1 = subplot(3,1,1); plot(ss2833.wam_t, ss2833.wam.jp(:,2));
             ax2 = subplot(3,1,2); plot(ss2833.force_t, ss2833.force_h(2,:));
             ax3 = subplot(3,1,3); plot(ss2833.wam_t, ss2833.wam.cf(:,2));
             linkaxes([ax1,ax2, ax3],'x');
             
%% 

dataTmp = data{1,4,2,8,2};
            figure;
            ax1 = subplot(3,1,1); plot(dataTmp.t, dataTmp.x(2,:));
            ax2 = subplot(3,1,2); plot(dataTmp.t, dataTmp.f(2,:));
            ax3 = subplot(3,1,3); plot(dataTmp.t, dataTmp.Fp(2,:));
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