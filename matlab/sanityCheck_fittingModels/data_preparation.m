clear; close all; clc;
% load('exampleData_subj8dir1.mat')
load('../sanityCheck_performSTF/test_data/exampleData_subj3dir1.mat')

%%
dex_trial_sel = 9; 


% data_f = trials(dex_trial_sel).data.f;
% data_ox = reshape(trials(dex_trial_sel).data.ox(:,:,1), ...
%     size(trials(dex_trial_sel).data.ox, [1 2]));
% data_t = trials(dex_trial_sel).data.t_shift;

data_f = recordedData(dex_trial_sel).data.f;
data_ox = reshape(recordedData(dex_trial_sel).data.ox(:,:,1), ...
    size(recordedData(dex_trial_sel).data.ox, [1 2]));
data_t = recordedData(dex_trial_sel).data.t_shift;


%
%figure(); 
clf; 
axh(1) = subplot(2,1,1); 
set(axh(1), 'linewidth', 1, 'fontsize', 12); hold on;
plot(data_t, data_ox(1,:)); 
axh(2) = subplot(2,1,2); 
set(axh(2), 'linewidth', 1, 'fontsize', 12); hold on;
plot(data_t, data_f(1,:));

linkaxes(axh, 'x');
xlim([-0.1 0.7])


%% 
t_zone = [-0.1 0.7];
dex_time_sel = data_t > t_zone(1) & data_t <= t_zone(2);
data_eg.t = data_t(dex_time_sel);
data_eg.F = data_f(1,dex_time_sel);
data_eg.x = data_ox(1,dex_time_sel);

plot(data_eg.t, data_eg.x)
plot(data_eg.t, data_eg.F)

save('data_tmp.mat', 'data_eg')

%% get data_tmp_batch.mat
% sel_idx  = [trials.tarF] == 20 & [trials.tarL] == 0.05 & [trials.outcome] == 1;
sel_idx  = [headers.trialHeader.tarF] == 20 & [headers.trialHeader.tarL] == 0.05 & [headers.trialHeader.outcome] == 1;
i = 0;
figure(); clf;
axh(1) = subplot(2,1,1); 
set(axh(1), 'linewidth', 1, 'fontsize', 12); hold on;
axh(2) = subplot(2,1,2); 
set(axh(2), 'linewidth', 1, 'fontsize', 12); hold on;

for sel_i = find(sel_idx)
    i = i+1;


    dex_trial_sel = sel_i; 


    bdata_f{i} = recordedData(dex_trial_sel).data.f;
    bdata_ox{i} = reshape(recordedData(dex_trial_sel).data.ox(:,:,1), ...
        size(recordedData(dex_trial_sel).data.ox, [1 2]));
    bdata_t{i} = recordedData(dex_trial_sel).data.t_shift;


    t_zone = [-0.1 0.7];
    dex_time_sel = bdata_t{i} > t_zone(1) & bdata_t{i} <= t_zone(2);
    bdata_eg{i}.t = bdata_t{i}(dex_time_sel);
    bdata_eg{i}.F = bdata_f{i}(1,dex_time_sel);
    bdata_eg{i}.x = bdata_ox{i}(1,dex_time_sel);


    plot(axh(1), bdata_t{i}, bdata_ox{i}(1,:)); 
    plot(axh(2), bdata_t{i}, bdata_f{i}(1,:));

    
end

save('bdata_tmp.mat', 'bdata_eg')

linkaxes(axh, 'x');
xlim([-0.1 0.7])