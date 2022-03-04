% See the release parameters without being perturb. 
%% Export data 
clear; close all; clc;
ss_num = {  [3938,3939],    [3944,3945],        3949 
            3940,           3943,        [3947,3948]
            3941,           3942,       3946} 
pert_f = 12; % only use 12N pert
pertT_num = 1 + 7;     % 1 without pert, and 7 perturbation time
freq = 500;
t_lim = [-0.1 0.9];
t = [t_lim(1) : 1/freq : t_lim(2)];
data = cell(1, size(ss_num,1), size(ss_num,2), 5, pertT_num); % The stochastic ones are not attached at the end 
dataMat = cell(3,3);
for frc_i = 1:size(ss_num,1) % actually force 
    for dist_i = 1:size(ss_num,2) % step perts
        % if multiple sessions in it 
%         ss_tmp = SessionScan();
        celltmp = cell(200,3);
        cell_idx_to = [1 1 1];
        for si = 1:length(ss_num{frc_i, dist_i})
            ss_tmp = SessionScan(ss_num{frc_i, dist_i}(si));
            cell_idx_from = [0 0 0];
%             celltmptmp = ss_tmp.export_as_formatted_hybridss(1);
            celltmptmp = ss_tmp.export_as_formatted_hybridss();
            for trial_i = 1:size(celltmptmp,1)
                % get the index from the 'from' cells
                if ~isempty(celltmptmp{trial_i,1})
                    cell_idx_from(1) = cell_idx_from(1) + 1;
                else
                    % do nothing 
                end
            end
                
            % check how many cells have data inside.
             celltmp(cell_idx_to(1)+[0:(cell_idx_from(1)-1)],1) = ...
                 celltmptmp(1:cell_idx_from(1),1);
            % add the index so that next session copy from here  
            cell_idx_to(1) = cell_idx_to(1) + cell_idx_from(1);
        end
        % interp each tial data into 1 1-by-500 array 
             % oversteak the same condition 
            % see how many trials leftover
            trial_num = 0;
            for trial_i = 1:size(celltmp,1)
                if ~isempty(celltmp{trial_i,1})
                    trial_num = trial_num + 1;
                end
            end

       dat_x = zeros(trial_num,length(t));
       dat_v = zeros(trial_num, length(t));
       dat_f = zeros(trial_num, length(t));
       for trial_i = 1:trial_num
           idx_release = find(celltmp{trial_i,1}.ts == 5);
           t_shift = celltmp{trial_i,1}.t - celltmp{trial_i,1}.t(idx_release(1));
           % interp data here
           dat_x(trial_i,:) = interp1(t_shift, celltmp{trial_i,1}.x(2,:), t, 'Linear', 'extrap');
           dat_v(trial_i,:) = interp1(t_shift, celltmp{trial_i,1}.v(2,:), t, 'Linear', 'extrap');
           dat_f(trial_i,:) = interp1(t_shift, celltmp{trial_i,1}.f(2,:), t, 'Linear', 'extrap');
           ifplot = 1;
           if(ifplot)
               clf; 
               axh(1) = subplot(3,1,1); hold on;
               plot(t_shift, celltmp{trial_i,1}.x(2,:), 'b.');
               plot(t, dat_x(trial_i,:), 'r*');
               title(['trial' num2str(trial_i)]);
               axh(2) = subplot(3,1,2); hold on;
               plot(t_shift, celltmp{trial_i,1}.v(2,:), 'b.');
               plot(t, dat_v(trial_i,:), 'r*');
               title(['trial' num2str(trial_i)]);
               axh(3) = subplot(3,1,3); hold on;
               plot(t_shift, celltmp{trial_i,1}.f(2,:), 'b.');
               plot(t, dat_f(trial_i,:), 'r*');
               title(['trial' num2str(trial_i)]);
               linkaxes(axh, 'x');
           end
       end
       dataMat{frc_i,dist_i}.x = dat_x;
       dataMat{frc_i,dist_i}.v = dat_v;
       dataMat{frc_i,dist_i}.f = dat_f;
    end
end
save('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/datatmp/Unperturbeddat.mat', 'dataMat');

%% Load the data 

% Pile the data up into a matrix 
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/datatmp/Unperturbeddat.mat', 'dataMat');
Dat_x = [];
for fce_i = 1:3
    for dist_i = 1:3
        dat_x = dataMat{fce_i,dist_i}.v;
        Dat_x = [Dat_x; dat_x];
    end
end


dat_x = Dat_x(:,50:10:400);
dat_x_center = mean(dat_x);
plot((dat_x - dat_x_center)');

dat_x_centered = dat_x - dat_x_center;

[coeff,score,latent] = pca(dat_x_centered);

plot(score(:,1), score(:,2), '.');

% Find the nearest neighbor using the data 
Dat_x = [];
Dat_v = [];
Dat_f = [];
for fce_i = 1:3
    for dist_i = 1:3
        dat_x = dataMat{fce_i,dist_i}.x;
        dat_v = dataMat{fce_i,dist_i}.v;
        dat_f = dataMat{fce_i,dist_i}.f;
        Dat_x = [Dat_x; dat_x];
        Dat_v = [Dat_v; dat_v];
        Dat_f = [Dat_f; dat_f];
    end
end
fh = figure(1); % see the raw data 
subplot(1,3,1); plot(Dat_x'); 
subplot(1,3,2); plot(Dat_v'); 
subplot(1,3,3); plot(Dat_f');

figure(2); % do the z-center 
Dat_x0 = Dat_x - mean(Dat_x);
Dat_x1 = Dat_x0./(ones(size(Dat_x0,1), 1) * std(Dat_x0));
% Perform the data into this matrix 

% visualize its distribution use dimention reduction 


% given a part of one trial, predict the other trials using existing data