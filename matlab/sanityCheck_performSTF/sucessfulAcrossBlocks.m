%% check blocks with accum sucessful rate 
%% BLOCK1. plot examplary subject 

load('test_data/exampleData_subj8dir4.mat')

successful = [headers.trialHeader.outcome];
block_num  = [headers.trialHeader.bNo];
block_chg = [0 diff(block_num)];
block_chgdex=find(block_chg ~=0);

sucessful_movAvg = smooth(successful, 5); 


%
% 1. a plot of successful rate 
fh = figure(); hold on;
set(gca, 'fontsize', 20);
set(gca, 'linewidth', 2);
plot(sucessful_movAvg, 'linewidth', 1); %hold on;
% plot(successful, '*'); hold on;
xline(block_chgdex, 'linewidth', 0.5)
xlabel('trials'); 
ylabel('successful rate');
title('moving avg successfulness');
legend('successful rate', 'change task condition');

saveas(fh, 'beh_figures/sucessful.movingAvgSuccR.subj8.dir1.allcond.png')

%% 
% check the trial count before succeed 
clear; close all; clc; 
subj_list=[2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 20 21]; 
fail_trials_bef = cell(18,4,10);
dir_label = {'+X', '-X', '+Y', '-Y'};
for subj_i = 1:18
    for dir_i  = 1:4
        fname = sprintf('exampleData_subj%ddir%d.mat', subj_list(subj_i), dir_i);
        fdir  = 'test_data';
        load([fdir, '/', fname], 'headers');

        trials  = 1:length(headers.trialHeader);
        outcome = [headers.trialHeader.outcome];
        bNo     = [headers.trialHeader.bNo];

        for b_i = 1:10
            outcome_tmp = logical(outcome(bNo == b_i));
            outcomef_cusum = cumsum(double(outcome_tmp == 0));

            trials_fail_cum = outcomef_cusum(outcome_tmp);
            trials_fail_bef = diff([0 trials_fail_cum]);

            if (length(trials_fail_bef) >= 9)
                fail_trials_bef{subj_i,dir_i,b_i} = trials_fail_bef(1:9);
            else % when there is less than 9 trials 
                trials_fail_bef_wn = nan(1,9);
                trials_fail_bef_wn(1:length(trials_fail_bef)) = trials_fail_bef;
                fail_trials_bef{subj_i,dir_i,b_i} = trials_fail_bef_wn;
            end

        end
    end
end

% plot out 
fh = figure(); 
% 1. box plot
for dir_i = 1:4
    axh(dir_i) =  subplot(4,1,dir_i); hold on;
    set(axh(dir_i), 'linewidth', 1);
    set(axh(dir_i), 'fontsize', 15);

    if (dir_i == 4) 
        xlabel('successful trial index')
    end 
    ylabel('count')
    title(dir_label(dir_i));

end

for dir_i = 1:4
    all_failedtrials_bef = []; 
for b_i = 1:10 
    block_failedtrials_bef = [];
    for subj_i = 1:18
        block_failedtrials_bef = [block_failedtrials_bef; fail_trials_bef{subj_i,dir_i,b_i}];
    end
    all_failedtrials_bef = [all_failedtrials_bef, block_failedtrials_bef];
end
    subplot(axh(dir_i));
    boxplot(all_failedtrials_bef)

    xticks([10:10:90]);
    xticklabels({'10', '20', '30', '40', '50', '60', '70', '80', '90'});
end
sgtitle('trial count prior successful')

saveas(fh, 'beh_figures/sucessful.trialscountpriorsucess.4dirs.allcond.png');

%% overlap every block (pool subjects, directions, and conditions) 
clear fh axh lnh

all_failedtrials_bef = [];
all_failedtrials_bef_x=[];
for dir_i = 1:4
    for b_i = 1:10
        block_failedtrials_bef = [];
        for subj_i = 1:18
            block_failedtrials_bef = [block_failedtrials_bef; fail_trials_bef{subj_i,dir_i,b_i}];
        end

        x_ = repmat(1:9, 18, 1)
        all_failedtrials_bef = [all_failedtrials_bef; block_failedtrials_bef ];
        all_failedtrials_bef_x = [all_failedtrials_bef_x; x_];

    end
end

all_failedtrials_bef   = all_failedtrials_bef(:);
all_failedtrials_bef_x = all_failedtrials_bef_x(:);

all_failedtrials_bef_dex = ~isnan(all_failedtrials_bef);
all_failedtrials_bef   = all_failedtrials_bef(all_failedtrials_bef_dex);
all_failedtrials_bef_x = all_failedtrials_bef_x(all_failedtrials_bef_dex);

fh = figure('unit', 'inch', 'position', [0 0 8 3]);
axh(1) = subplot(1,2,1); hold on; 
set(axh(1), 'fontsize', 15);
set(axh(1), 'linewidth', 1);
boxplot(all_failedtrials_bef, all_failedtrials_bef_x, ...
    'OutlierSize', 2, 'Jitter', 0.3);
xlabel('trials'); ylabel('failed trial count');
title('bar plot');

axh(2) = subplot(1,2,2); hold on;
set(axh(2), 'fontsize', 15);
set(axh(2), 'linewidth', 1);
scatter(all_failedtrials_bef_x(:), all_failedtrials_bef(:), 20);
xlabel('trials'); ylabel('failed trial count');
lnh = refline();
set(lnh, 'linewidth', 2);
title('scatters');

linkaxes(axh, 'y')

% try bar plot 
saveas(fh, 'beh_figures/sucessful.trialscountpriorsucess.pooled_subj_dir_cond.png');
ylim(axh(1), [0 5]);
saveas(fh, 'beh_figures/sucessful.trialscountpriorsucess.pooled_subj_dir_cond.zoomin.png');

% %%  2. many scatter plot
close all; 
hline_k = zeros(4,10);
hline_x = zeros(4,10,2); 
hline_y = zeros(4,10,2); 
for dir_i = 1:4
figure('unit', 'inch', 'position', [0 0 3 14], 'visible', 'off');
for b_i = 1:10 
    block_failedtrials_bef = [];
    for subj_i = 1:18
        block_failedtrials_bef = [block_failedtrials_bef; fail_trials_bef{subj_i,dir_i,b_i}];
    end
    subplot(10,1,b_i);
    x_ = repmat(1:9, 18, 1);
    plot(x_(:), block_failedtrials_bef(:), 'o')
    hline = refline;
    hline_k(dir_i, b_i) = (hline.YData(2) - hline.YData(1)) / 10; 
%     if hline.YData(2) - hline.YData(1) > 0
%         set(hline,'Color', 'r')
%     end
    hline_x(dir_i,b_i,:) = hline.XData;
    hline_y(dir_i,b_i,:) = hline.YData;

    set(hline, 'linewidth', 3)
end

    sgtitle(['dir' num2str(dir_i)])
end


%% many box plot 
for dir_i = 1:4
    fh = figure('unit', 'inch', 'position', [0 0 3 14]);
    all_failedtrials_bef = [];
    for b_i = 1:10
        block_failedtrials_bef = [];
        for subj_i = 1:18
            block_failedtrials_bef = [block_failedtrials_bef; fail_trials_bef{subj_i,dir_i,b_i}];
        end
        subplot(10,1,b_i); hold on;
        x_ = repmat(1:9, 18, 1);
        boxplot(block_failedtrials_bef(:), x_(:), 'PlotStyle', 'compact')

        line_x = reshape(hline_x(dir_i,b_i,:),1,2);
        line_y = reshape(hline_y(dir_i,b_i,:),1,2);

        if (diff(line_y) < 0)
            line(line_x, line_y, 'linewidth', 3, 'color', 'b');
        else
            line(line_x, line_y, 'linewidth', 3, 'color', 'r');
        end

        if dir_i == 1
            ylabel(['block' num2str(b_i)]);
        else
            ylabel('');
        end

    end
    if b_i == 10
        xlabel('trial index')
    end
    sgtitle(['dir' num2str(dir_i)])

    saveas(fh, ['beh_figures/sucessful.trialscountpriorsucess.pooled_subj_dir' num2str(dir_i) '.png']);
end




%%  
% check sucessful rate at release time 
clear; close all; clc;
load('test_data/exampleData_subj8dir1.mat', 'headers')
headers_all{1} = headers;
load('test_data/exampleData_subj9dir1.mat', 'headers')
headers_all{2} = headers;
load('test_data/exampleData_subj10dir1.mat', 'headers')
headers_all{3} = headers;
load('test_data/exampleData_subj12dir1.mat', 'headers')
headers_all{4} = headers;

figure(); hold on;
for subj_dex = 1:4
    headers = headers_all{subj_dex};

    successful = [headers.trialHeader.outcome];
    block_num  = [headers.trialHeader.bNo];
    block_chg = [0 diff(block_num)];
    block_chgdex=find(block_chg ~=0);

    sucessful_movAvg = smooth(successful, 5);


    %
    % 1. a plot of successful rate
    subplot(4,1,subj_dex)
    set(gca, 'fontsize', 20);
    set(gca, 'linewidth', 2);
    plot(sucessful_movAvg, 'linewidth', 1); %hold on;
    % plot(successful, '*'); hold on;
    xline(block_chgdex, 'linewidth', 0.5)
    xlabel('trials');
    ylabel('successful rate');
    legend('successful rate', 'change task condition');
end
sgtitle('moving avg successfulness');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cg: move this to other files...
% 2. check a data as release 
load('exampleData_subj8dir1.mat')

trials_num = length(trials); 
t_range = [-0.1 1.5];
fh = figure();
for trial_dex = 1:length(trials)

    t_orig = trials(trial_dex).data.t_shift;
    t_dex  = t_orig > t_range(1) & t_orig < t_range(2); 
    t = t_orig(t_dex);
    x = trials(trial_dex).data.x(1,t_dex);
%     x = trials(trial_dex).data.ox(1,t_dex,1); % caucious using this as
%     not all was recorded
    x_hold = mean(x(t>-0.1 & t<0));
    x = x - x_hold;
    F = trials(trial_dex).data.f(1,t_dex);
    
    
    clf;
    axh(1) = subplot(2,1,1);
    axh(2) = subplot(2,1,2); 
    if (trials(trial_dex).outcome == 1)
        plot(axh(1), t, x, 'b');
        plot(axh(2), t, F, 'b');
        sgtitle(num2str(trials(trial_dex).tNo)); 
    else % fail
        plot(axh(1), t, x, 'r');
        plot(axh(2), t, F, 'r');
        sgtitle({num2str(trials(trial_dex).tNo); trials(trial_dex).failReason}); 
    end

    saveas(fh, ['trials_fig_long/trial' num2str(trial_dex) '.png'])
    
end
% save('exampleData_subj8dir1.mat', 'trials')

%

%% plot estimations against trials 
x_tars = trials(:).tarL;
F_tars = trials(:).tarF;
k_tars = [trials(:).tarF]./[trials(:).tarL]; 
% 
k_ests = [trials(:).estK];
b_ests = [trials(:).estB];
m_ests = [trials(:).estM]; 
fit_ests = [trials(:).estFIT]; 

outcome= [trials(:).outcome];

% figure 1 plot k
figure(1); 
trials_all = 1:length(trials); hold on; 
plot(trials_all, k_tars); 
plot(trials_all, k_ests); 
plot(trials_all(outcome), k_ests(outcome), 'g*');
plot(trials_all(~outcome), k_ests(~outcome), 'r*');

%% figure 2 plot FIT
figure(2); hold on; 
trials_all = 1:length(trials); 
% plot(trials_all, k_tars); 
plot(trials_all(outcome), fit_ests(outcome), 'g.');
plot(trials_all(~outcome), fit_ests(~outcome), 'r.');
%% figure 3 plot m
figure(3); 
trials_all = 1:length(trials); hold on; 
% plot(trials_all, k_tars); 
plot(trials_all(outcome), m_ests(outcome), 'g.');
plot(trials_all(~outcome), m_ests(~outcome), 'r.');

%%%
% Run a few other subjects and directions. Discuss the stiffness meeting. 