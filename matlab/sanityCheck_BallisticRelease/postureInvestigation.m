% Finding postures that can do 2D movements. 

% |      stoc force pert             |   3N  |   5N      |
% |----------------------------------|-------|-----------|
% |posture 1. jp = [0 -1.2 0 2.249]; | 4010  |  4011     |
% |posture 2. jp = [0 0 0 1.5708];   | 4013  |  4012     |

ss_mat = [4010, 4011; 
          4013, 4012];

data = cell(1,1,2,2,15,3); % 
for pos_i = 1:2
    for fce_i = 1:2
        ss_tmp = SessionScan(ss_mat(pos_i,fce_i)); 
        celltmp = ss_tmp.export_as_formatted_hybridss();
        
        try
            data(1,1,fce_i, pos_i, 1:15, 3) = celltmp(1:15,3);
        catch 
            trial_num = size(celltmp,1);
            data(1,1,fce_i, pos_i, 1:trial_num, 3) = celltmp(1:trial_num,3);
        end
    end
end
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4010_4013.mat', 'data'); 

%% load data and plot it out 
clear; close all;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4010_4013.mat', 'data'); 
fce_list = [3, 5];
pos_list = [1, 2];
for pos_i = 1:2
    for fce_i = 1:2
        fh(pos_i, fce_i) = figure('unit', 'inch', 'position', [0 0 10 6]); 
        trials_num = 5;
        for trial_i = 1:trials_num
            trialtmp = data{1,1,fce_i, pos_i, trial_i, 3};
            subplot(trials_num,2,(trial_i-1)*2 + 1); 
            plot(trialtmp.t, trialtmp.Fp); 
            subplot(trials_num,2,(trial_i-1)*2 + 2);
            plot(trialtmp.t, trialtmp.v); 
            if (trial_i == 1)
                legend('x', 'y', 'z');
            end
            ylim([-0.4 0.4])
%             xticklabels({});
        end
        sgtitle(['posture' num2str(pos_list(pos_i)), 'pert fce' num2str(fce_list(fce_i)) 'N']);
    end
end

%% 3-min data
clear
ss_mat = [4017, 4022, 4019];

data = cell(1,1,1,3,1,3); % 
t_range = [2974.1704, 3153.168
    2140.484, 2299.484
    3999.62, 4178.62];
            
for pos_i = 1:3

% for trial 2
% some procedure to deal with session 4022
sstmp = SessionScan(ss_mat(pos_i));
clf; plot(sstmp.data.t, sstmp.data.Fp)
t_idx       = sstmp.data.t>t_range(pos_i,1) & sstmp.data.t<t_range(pos_i,2);
dattmp.t    = sstmp.data.t(t_idx);
dattmp.x    = sstmp.data.x(:,t_idx);
dattmp.v    = sstmp.data.v(:,t_idx);
dattmp.f    = sstmp.data.f(:,t_idx);
dattmp.ftq  = sstmp.data.ftq(:,t_idx);
dattmp.Fp   = sstmp.data.Fp(:,t_idx);
dattmp.ts   = sstmp.data.ts(t_idx);
dattmp.tq   = sstmp.data.tq(:,t_idx);
dattmp.mvst = zeros(size(dattmp.t));
data{1,1,1,pos_i,1,3} = dattmp;

end
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4017_4019.mat', 'data'); 

% load data and plot it out 
clear; close all;
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4017_4019.mat', 'data'); 
fce_list = [5];
pos_list = [1, 2, 3];

fh = figure('unit', 'inch', 'position', [0 0 10 6]); 
for pos_i = 1:3
    subplot(3,1,pos_i);
    fce_i = 1;
        trials_num = 1;
        for trial_i = 1:trials_num
            trialtmp = data{1,1,fce_i, pos_i, trial_i, 3};
%             subplot(trials_num,2,(trial_i-1)*2 + 1); 
%             plot(trialtmp.t, trialtmp.Fp); 
%             subplot(trials_num,2,(trial_i-1)*2 + 2);
            plot(trialtmp.t, trialtmp.x - trialtmp.x(:,1)); 
            if (pos_i == 1)
                legend('x', 'y', 'z');
            end
%             ylim([-0.4 0.4])
%             xticklabels({});
        end
        sgtitle(['posture' num2str(pos_list(pos_i)), 'pert fce' num2str(fce_list(fce_i)) 'N']);

end

% cg's note: I should mannually process the data in ss4021/4022 to get the
% proper data out.


