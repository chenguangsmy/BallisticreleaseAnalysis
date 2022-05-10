% This file check the ballistic-release on 4 directional movements. 
% It takes the both OPT and WAM positions 

% sanity check contents: 
% 1. Compare the y directional OPTOTRAK position and the WAM position 
% 2. Compare the y directional force vs x directional force 







% |sessions:| 2.5cm    | 5cm        | 7.5cm     | 
% | ------- | -------- | ---------  | --------- | 
% | 15N     | 4062     | 4061       | 4060      | 
% | 20N     | 4063     | 4065       | 4064      | 
% | 25N     | 4067     | 4066       | 4068      | 

% four example sesions in the ballistic-release testing. 

%%
% 1. Compare the y directional OPTOTRAK position and the WAM position 
% (ss4060)
sstmp = SessionScan(4060);
celltmp = sstmp.export_as_formatted(1);
celltmp1 = reshape(celltmp(1,1,:,2), size(celltmp,3),1);

%% see the position at the time of release 
figure('unit', 'inch', 'position', [0 0 3 8], ...
        'name', 'release_displacement'); 
hold on;
for trial_i = 1:13 %-valid trials
    axh(1) = subplot(3,1,1);
    hold on;
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.x(1,:), 'b');
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.ox(1,:), 'r');
    ylabel('x (m)');
    
    axh(2) = subplot(3,1,2);
    hold on;
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.x(2,:), 'b');
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.ox(2,:), 'r');
    ylabel('y (m)');
    
    axh(3) = subplot(3,1,3);
    hold on;
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.x(3,:), 'b');
    plot(celltmp1{trial_i}.t, celltmp1{trial_i}.ox(3,:), 'r');
    ylabel('z (m)');
end
linkaxes(axh(1:3), 'x');
xlabel('t (s)');
xlim([-0.1 0.5])
sgtitle('WAM(x) vs OPT (x)')

% see the position at the time of perturb
figure('unit', 'inch', 'position', [0 0 3 8], ...
        'name', 'pert_displacement'); 
hold on;
for trial_i = 1:13 %-valid trials
    pert_idx = find(celltmp1{trial_i}.Fp(2,:) ~= 0);
    t_shift = celltmp1{trial_i}.t - celltmp1{trial_i}.t(pert_idx(1));
    axh(1) = subplot(3,1,1);
    hold on;
    plot(t_shift, celltmp1{trial_i}.x(1,:), 'b');
    plot(t_shift, celltmp1{trial_i}.ox(1,:), 'r');
    ylabel('x (m)');
    
    axh(2) = subplot(3,1,2);
    hold on;
    plot(t_shift, celltmp1{trial_i}.x(2,:), 'b');
    plot(t_shift, celltmp1{trial_i}.ox(2,:), 'r');
    ylabel('y (m)');
    
    axh(3) = subplot(3,1,3);
    hold on;
    plot(t_shift, celltmp1{trial_i}.x(3,:), 'b');
    plot(t_shift, celltmp1{trial_i}.ox(3,:), 'r');
    ylabel('z (m)');
end
linkaxes(axh(1:3), 'x');
xlabel('t (s)');
xlim([-0.1 0.5])
sgtitle('WAM(x) vs OPT (x)')

% see the position overlap at the time of perturb
figure('unit', 'inch', 'position', [0 0 3 8], ...
        'name', 'pert_displacement'); 
hold on;
for trial_i = 1:13 %-valid trials
    pert_idx = find(celltmp1{trial_i}.Fp(2,:) ~= 0);
    t_shift = celltmp1{trial_i}.t - celltmp1{trial_i}.t(pert_idx(1));
    x_shift = celltmp1{trial_i}.x(:,:) - celltmp1{trial_i}.x(:,pert_idx(1));
    ox_shift = celltmp1{trial_i}.ox(:,:) - celltmp1{trial_i}.ox(:,pert_idx(1));
    axh(1) = subplot(3,1,1);
    hold on;
    plot(t_shift, x_shift(1,:), 'b');
    plot(t_shift, ox_shift(1,:), 'r');
    ylabel('x (m)');
    
    axh(2) = subplot(3,1,2);
    hold on;
    plot(t_shift, x_shift(2,:), 'b');
    plot(t_shift, ox_shift(2,:), 'r');
    ylabel('y (m)');
    
    axh(3) = subplot(3,1,3);
    hold on;
    plot(t_shift, x_shift(3,:), 'b');
    plot(t_shift, ox_shift(3,:), 'r');
    ylabel('z (m)');
end
linkaxes(axh(1:3), 'x');
xlabel('t (s)');
xlim([-0.1 0.5])
sgtitle('WAM(x) vs OPT (x), subtracted 0')

%%
% 2. Compare the y directional force vs x directional force 
% (ss4065 v.s. ss4069)

