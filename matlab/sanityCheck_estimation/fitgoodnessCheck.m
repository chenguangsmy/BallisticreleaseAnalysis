% Check the fit goodness (r2) by conditions. 
% Using Federico's fitting, and check the reason why some conditions have
% less fitting goodness than others 

% plot the fitting on different conditions
clear; clc; close all;
% load('../data/processedData/ss4379_4438_OPTupdate_results.mat');  subj = [1:6];
% load('../data/processedData/ss4446_4467_results.mat');            subj = [7:9];
% load('../data/processedData/ss4472_4495_results.mat');            subj = [10:12];
% load('../data/processedData/ss4500_4524_results.mat');            subj = [13:15];
% load('../data/processedData/ss4530_4563_results.mat');            subj = [16:18];
% load('../data/processedData/ss4573_4587_rsults.mat');             subj = [19 20];

load('../data/processedData/all_subj_results.mat', 'results');      subj = 1:20;

for i = 1:20
    subj_i = subj(i);
    for dir_i = 1:4
        results{subj_i, dir_i};
        fh = figure(); 
        fce     = [15 20 25];
        disp    = [0.025 0.05 0.075];

        for fce_i = 1:3
            for disp_i = 1:3
                subplot(3,3, disp_i+(fce_i-1)*3);
                hold on;
                t       = results{subj_i, 1}.FD_UP{fce_i,disp_i}{1};
                f       = results{subj_i, 1}.FD_UP{fce_i,disp_i}{2};
                x       = results{subj_i, 1}.FD_UP{fce_i,disp_i}{3};
                x_pred  = results{subj_i, 1}.FD_UP{fce_i,disp_i}{4};
                goodfit = results{subj_i, 1}.FIT_up(fce_i, disp_i);

                plot(t, x, 'k:');
                plot(t, x_pred, '-', 'color', [0.8 0.8 0.8]);

                %         title([num2str(fce(fce_i)), 'N ', num2str(disp(disp_i)), 'cm score:', num2str(goodfit)]);
                title(['score:', num2str(goodfit)]);
            end

        end
        sgtitle(['subj' num2str(subj_i), 'dir' num2str(dir_i)]);

        if min(results{subj_i, 1}.FIT_up(:)) > 85
            close(fh)
        end
    end
end