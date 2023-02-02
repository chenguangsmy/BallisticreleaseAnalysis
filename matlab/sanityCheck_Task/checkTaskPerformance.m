% Check task performance around different subjects... 
%% save all of them 
ss_num = {  ...
                [4401 4402 4404 4405 4426 4427] ...     % BH conducting
                [4408 4409 4412 4414] ...               % MR conducting
                [4418 4419 4422 4421] ...               % NN conducting
                [4432 4434 4437 4438] ...               % HM conducting
                [4379 4380 4382 4383] ...               % CZ testing
                [4385 4386 4387 4388] ...               % HA testing
                [4446 4448 4450 4451]  ...              % FM conducting
                [4455 4456 4458 4459]  ...              % FT conduting
                [4463 4464 4466 4467]  ...              % QX conducting
                [4472 4473 4476 4477] ...               % VC
                [4481 4482 4483 4485 4486] ...          % DS
                [4491 4492 4494 4495] ...               % BW
                [4500 4501 4503 4504] ...               % AS
                [4512 4513 4515 4516] ...               % XZ
                [4520 4521 4523 4524] ...               % ZC
                [4530 4531 4533 4534] ...               % KO
                [4542 4543 4545 4546] ...               % SL
                [4558 4560 4562 4563] ...               % AK %                 [4568 4569] ...                         % AR % unfinished
                [4573 4574 4576 4577] ...               % RL
};  


for subj_i = 1:length(ss_num) % 15
    ssstmp = SessionsScan(ss_num(subj_i))
    perf(subj_i) = ssstmp.getTrialPerformance()
end

save('subj_perf.mat', 'perf');

%% load all of them 
clear; clc; close all; 

load('subj_perf.mat', 'perf');

% 1. The sucessful rate by directions...
% plot the sucessful rate finished by time 
subj_exclude = 1;
subj_idx = setdiff(1:length(perf), subj_exclude);
subj_num = length(perf); 
time_idx = [1 3 2 4]; 
rate_all_dir = zeros(subj_num, 4)
for subj_i = 1:subj_num 
    rate_all_dir(subj_i,:) = perf{subj_i}.rate_eachd(time_idx);
end
time_seq = repmat([1 2 3 4], subj_num , 1);
scatter(time_seq(:), rate_all_dir(:), 'b'); 

% do the bar plot with scattered on... 
figure(); hold on;
bar(1:4, mean(rate_all_dir(:,:)));
errorbar(mean(rate_all_dir(:,:)), std(rate_all_dir(:,:)), 'LineWidth',3, 'Color',[0 0 0]); 
time_randomnized_mat = time_seq(:,:)+(rand(subj_num,4)-0.5)/2;
scatter(time_randomnized_mat(:), rate_all_dir(:), 'b'); 
xticks([1:4]);
xticklabels({'1 (right)', '2 (left)', '3 (front)', '4 (back)'});
xlabel('The sequence of blocks'); 
ylabel('The sucessful rate');
legend('sucessful rate avg','sucessful rate std','single subject');
title('performance rate across directions'); 

% 2. The sucessful rate by each conditions 
figure();
for dir = 1:4
    subplot(1,4,dir);  hold on;
rate_all_cond = zeros(subj_num, 9);
for subj_i = setdiff(1:subj_num, subj_exclude)
    rate_all_cond(subj_i,:) = [reshape(perf{subj_i}.rate_eachc(1,dir,1,:),1,3) ...
                              reshape(perf{subj_i}.rate_eachc(1,dir,2,:),1,3) ...
                              reshape(perf{subj_i}.rate_eachc(1,dir,3,:),1,3) ];
end
% cond_seq = repmat([1:9], subj_num , 1);
cond_seq = repmat([7 3 1 8 5 2 9 6 4], subj_num , 1); % according to stiffness levels
bar(1:9, mean(rate_all_cond(subj_idx,:)));
errorbar(mean(rate_all_cond(subj_idx,:)), std(rate_all_cond(subj_idx,:)), 'LineWidth',3, 'Color',[0 0 0]); 
cond_randomnized_mat = cond_seq(:,:)+(rand(subj_num,9)-0.5)/2;
rate_all_cond_idx = rate_all_cond(subj_idx,:);
cond_randomnized_mat_idx = cond_randomnized_mat(subj_idx,:);
scatter(cond_randomnized_mat_idx(:), rate_all_cond_idx(:), 'b'); 
xticks([1:9]);
ylim([0 1]);
% xticklabels({'15N 2.5cm', '15N 5.0cm', '15N 7.5cm', '20N 2.5cm', '20N 5.0cm', '20N 7.5cm', '25N 2.5cm', '25N 5.0cm', '25N 7.5cm'});
xticklabels({'15N 7.5cm', '20N 7.5cm', '15N 5.0cm', '25N 7.5cm', '20N 5.0cm', '25N 5.0cm', '15N 2.5cm', '20N 2.5cm', '25N 2.5cm'});
xlabel('The sequence of blocks'); 
ylabel('The sucessful rate');
legend('sucessful rate avg','sucessful rate std','single subject');
title({'performance rate across condition', ...
    ['direction' num2str(dir)]}); 
end
sgtitle('finish rate across conditions');

%% get the sucessful rate 
perf_rate = zeros(1, size(perf,2));
for subj_i = 1:length(perf_rate)
    perf_rate(subj_i) = perf{subj_i}.rate_all
end

perf_idx = setdiff(1:length(perf_rate), [1]); % only removal subject 1
fprintf('successful rate: mean:')
mean(perf_rate(perf_idx)) 
fprintf('std:')
std(perf_rate(perf_idx))
fprintf('min:')
min(perf_rate(perf_idx))
fprintf('max:')
max(perf_rate(perf_idx))

%% get the duration of experiment
perf_dur =  zeros(1, size(perf,2));
for subj_i = 1:length(perf_dur)
    perf_dur(subj_i) = perf{subj_i}.time_all
end

perf_idx = setdiff(1:length(perf_dur), [1]); % only removal subject 1
fprintf('Duration, mean:')
mean(perf_dur(perf_idx)) 
fprintf('std:')
std(perf_dur(perf_idx))
fprintf('min:')
min(perf_dur(perf_idx))
fprintf('max:')
max(perf_dur(perf_idx))

