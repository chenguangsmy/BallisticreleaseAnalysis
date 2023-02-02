% pulseParameterSelection 

% The pulse parameter was selected carelessly prior to 2022-05-24. We
% choose 300ms wide 6N pulse is because it is small enough (peak force less
% than 30N), and it do not change the displacement measurement too much.  
% However, as we have the optotrak running, we could have totally
% independent displacement and force measureremnt, which makes the thinner
% pulse possible.  

% 1. Export the data 
ss_num = [4229 4230 4231];
data = cell(3, ...                                       % pulse format
            4, ...                                       % directions
            3, ...                                       % fce
            3, ...                                       % dist
            7, ...                                       % trials
            2);                                          % 1. release, 2. pulse

trialout_num = 7;

obj.export_cond.direction = 1; 
obj.export_cond.fce = [15 20 25];
obj.export_cond.disp = [0.025 0.05 0.075];

for subj_i = 1:3
    sstmp = SessionScan(ss_num(subj_i));
    cond.fce = [sstmp.trials.tarF];
    cond.disp = [sstmp.trials.tarL];
    cond.sf = [sstmp.trials.outcome];
    cond.pert = [sstmp.trials.ifpert] + 1;
    for direction_i = 1
        for fce_i = 1:3
            for disp_i = 1:3
                for pert_i = 1:2
                    trialMask = [cond.fce == obj.export_cond.fce(fce_i) & ...
                        cond.disp == obj.export_cond.disp(disp_i) & ...
                        cond.pert == pert_i & ...
                        cond.sf == 1
                        ];


                    if (sum(trialMask)) == 0 % no trials
                        continue;
                    end

                    trial_idx = find(trialMask);
                    % if trial is enough, get trials,
                    if (sum(trialMask))>=trialout_num
                        trial_idx = find(trialMask);
                        trial_idx = trial_idx(1:trialout_num);
                        % if trial is not enough, get more trials
                        % from repeating
                    else
                        trials_qualify_num = sum(trialMask);
                        trials_lack = trialout_num - trials_qualify_num;

                        for trials_lacki = 1:trials_lack
                            trial_idx(trials_qualify_num+trials_lacki) = ...
                                trial_idx(trials_lacki);
                            disp(['put trial' num2str(trials_lacki) 'in slot' num2str(trials_qualify_num+trials_lacki)]);
                        end
                    end

                    for trial_idx_dest = 1:length(trial_idx)
                        trial_idx_from = trial_idx(trial_idx_dest);
                        data{subj_i,direction_i,fce_i,disp_i,trial_idx_dest,pert_i} = ...
                            sstmp.trials(trial_idx_from).export_as_formatted(); % need edition.
                        %                                         obj.trials(trial_idx_from).export_as_formatted(1); % need edition.
                    end
                end
            end
        end
    end
end

% save 
params.height = [6, 12, 24];
params.width = [300 150 75]; % ms
save('ss4229_4231_pulseParameterCompare', 'data', 'params', '-append');

% shift to the different pulse format 
data1 = cell(1, ...                                       % pulse format
            4, ...                                       % directions
            3, ...                                       % fce
            3, ...                                       % dist
            7, ...                                       % trials
            4);                                          % 1. release, 2. pulse1 3. pulse 2 4. pulse 3

data1(1,:,:,:,:,1) = data(1,:,:,:,:,1);
data1(1,:,:,:,:,2) = data(1,:,:,:,:,2);
data1(1,:,:,:,:,3) = data(2,:,:,:,:,2);
data1(1,:,:,:,:,4) = data(3,:,:,:,:,2); 
data = data1;
save('ss4229_4231_pulseParameterCompare_', 'data');


%% 2. show the data

dist_list = [2.5 5.0 7.5];
fce_list = [15 20 25];

% load('pulseParameterCompare.mat');       % parameter selection
Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 
colors = colormap('lines');
close(fh);
r = size(Data, 1); % subj
d = size(Data, 2); % direction
f = size(Data, 3); % force
l = size(Data, 4); % distance
t = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 1;

%%%% plot the force
epoc_type = 1;  % 1 perturb
                % 2 release
plot_type = 8;  % 1 displacement
                % 2 force 
                % 3 feedforward force
                % 4 velocity
                % 5 torque
                % 6 norm displacement
                % 7 norm force
                % 8 opt displacement
                % 9 opt velocity
                % 10 opt 1
                % 11 opt 2
                % 12 opt 3

pert_type = 2; % choose option [2 3 4]
axh = zeros(d,r);
xyi = 1;        % x, y

fh = figure(); hold on;
% pert_type = 2; colors = colors(4:end,:)
for ri = 1:r % perturb type
    for di = 1%:d % direction
%         fh = figure('Name', ['direction' num2str(di)]);
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        %         axh(ri, di) = subplot(d,r,r*(di-1) + ri);grid on;hold on;
        %         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
        for fi = 1:f % target force
            for li = 1:3 % target distance
                axh(fi, li) = subplot(f,l,f*(li-1) + fi);grid on;hold on; %?
                for pi = 2:p % perturbation
                    trial_num = length(Data(ri,di,fi,li,:,pi));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,di,fi,li,ti,pi}))
                            continue;
                        end
                        
                        switch epoc_type
                            case 1
%                                 idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0 & Data{ri,di,fi,li,ti,pi}.ts==4);  % pert at y
                                idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0);
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if pi == 1
                                    disp('ERROR: should use pi == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,di,fi,li,ti,pi}.ts==5 | Data{ri,di,fi,li,ti,pi}.ts==6);
                                idx = (idx(1)-100):idx(end);
                                %idx = (idx(1)):idx(end);
                        end
                        %plot(Data{ri,ci,di,ti,li}.Fp(xyi,:));
                        %idx = find(Data{ri,ci,di,ti,li}.Fp(xyi,:)~=0);
                        %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                        %idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                        %idx = (idx-9):idx(end);
                        %idx = (idx(1)):(idx(end)+100);
                        time = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,di,fi,li,ti,pi}.x(xyi,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,di,fi,li,ti,pi}.f(xyi,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,di,fi,li,ti,pi}.Fp(xyi,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,di,fi,li,ti,pi}.v(xyi,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,di,fi,li,ti,pi}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,di,fi,li,ti,pi}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,di,fi,li,ti,pi}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm force';
                            case 8 % optotrak x
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ox))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                end
                                titlestr = 'opto-position';
                            case 9 % optotrak v
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ov))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx,1);
                                end
                                titlestr = 'opto-velocity';
                            case 10 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                titlestr = 'opto1-position';
                            case 11
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,2);
                                titlestr = 'opto2-position';
                            case 12 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,3);
                                titlestr = 'opto3-position';
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:5));
                        end
%                         plot(time, dat, 'Color', colors(4+li, :));
%                         plot(time, dat, '.', 'Color', colors(4+li + 4*(ri-1), :));
                        plot(time, dat, '.', 'Color', colors(ri, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
                title(['fce' num2str(fce_list(fi)) 'dist' num2str(dist_list(li))]);

            end
        end
        linkaxes(axh, 'xy');
        xlim([0 1.0])
    end
end
linkaxes(axh, 'xy');
xlim([0 1.0])
if (if_subtract)
    ylim([-8 25]);
else
    ylim([6 50]);
end
% xlim([0 0.5]);
% xlim([0 2])
sgtitle(titlestr);



%% %% plot the displacement
epoc_type = 1;  % 1 perturb
                % 2 release
plot_type = 8;  % 1 displacement
                % 2 force 
                % 3 feedforward force
                % 4 velocity
                % 5 torque
                % 6 norm displacement
                % 7 norm force
                % 8 opt displacement
                % 9 opt velocity
                % 10 opt 1
                % 11 opt 2
                % 12 opt 3

pert_type = 2; % choose option [2 3 4]
axh = zeros(d,r);
xyi = 1;        % x, y

fh = figure(); hold on;
% pert_type = 2; colors = colors(4:end,:)
for ri = 1:r % perturb type
    for di = 1%:d % direction
%         fh = figure('Name', ['direction' num2str(di)]);
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        %         axh(ri, di) = subplot(d,r,r*(di-1) + ri);grid on;hold on;
        %         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
        for fi = 1:f % target force
            for li = 1:3 % target distance
                axh(fi, li) = subplot(f,l,f*(li-1) + fi);grid on;hold on; %?
                for pi = pert_type%1:p % perturbation
                    trial_num = length(Data(ri,di,fi,li,:,pi));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,di,fi,li,ti,pi}))
                            continue;
                        end
                        
                        switch epoc_type
                            case 1
%                                 idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0 & Data{ri,di,fi,li,ti,pi}.ts==4);  % pert at y
                                idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0);
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if pi == 1
                                    disp('ERROR: should use pi == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,di,fi,li,ti,pi}.ts==5 | Data{ri,di,fi,li,ti,pi}.ts==6);
                                idx = (idx(1)-100):idx(end);
                                %idx = (idx(1)):idx(end);
                        end
                        %plot(Data{ri,ci,di,ti,li}.Fp(xyi,:));
                        %idx = find(Data{ri,ci,di,ti,li}.Fp(xyi,:)~=0);
                        %idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                        %idx = find(Data{ri,ci,di,ti,li}.ts==5 | Data{ri,ci,di,ti,li}.ts==6);
                        %idx = (idx-9):idx(end);
                        %idx = (idx(1)):(idx(end)+100);
                        time = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,di,fi,li,ti,pi}.x(xyi,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,di,fi,li,ti,pi}.f(xyi,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,di,fi,li,ti,pi}.Fp(xyi,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,di,fi,li,ti,pi}.v(xyi,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,di,fi,li,ti,pi}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,di,fi,li,ti,pi}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,di,fi,li,ti,pi}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(xyi,:));
                                titlestr = 'norm force';
                            case 8 % optotrak x
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ox))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                end
                                titlestr = 'opto-position';
                            case 9 % optotrak v
                                switch length(size(Data{ri,di,fi,li,ti,pi}.ov))
                                    case 2
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx);
                                    case 3
                                        dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx,1);
                                end
                                titlestr = 'opto-velocity';
                            case 10 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                                titlestr = 'opto1-position';
                            case 11
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,2);
                                titlestr = 'opto2-position';
                            case 12 
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,3);
                                titlestr = 'opto3-position';
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:5));
                        end
%                         plot(time, dat, 'Color', colors(4+li, :));
%                         plot(time, dat, '.', 'Color', colors(4+li + 4*(ri-1), :));
                        plot(time, dat, '.', 'Color', colors(ri, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
                title(['fce' num2str(fce_list(fi)) 'dist' num2str(dist_list(li))]);

            end
        end
        linkaxes(axh, 'xy');
        xlim([0 1.0])
    end
end
linkaxes(axh, 'xy');
xlim([0 1.0])
if (if_subtract)
    ylim([-7 4]*1e-3);
else
    ylim([-0.492 -0.473]);
end
% xlim([0 0.5]);
% xlim([0 2])
sgtitle(titlestr);


%%

% code preparing to run Federico's system identification 
load('ss4229_4231_pulseParameterCompare', 'data');  % just release, Chenguang

% turn data into 1 subject with 3 pulses 
% data1 = cell([1,size(data,2:5),4]);
% data1(1,:,:,:,:,1:2) = data(1,:,:,:,:,:);
% data1(1,:,:,:,:,3) = data(2,:,:,:,:,2);
% data1(1,:,:,:,:,4) = data(3,:,:,:,:,2);
% data = data1;
% size(data)


Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 

colors = colormap('lines');
r = size(Data, 1); % subj
c = size(Data, 2); % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;
epoc_type = 2;
plot_type = 1;          % 1 displacement
                        % 2 force
                        % 3 force command
                        % 4 velocity
                        % 5 torque of 3rd joint
                        % 6 displacement (vector length in cartesian space)
                        % 7 force (vector length in cartesian space)
axh = zeros(f,r);



for ri = 1:r % subj
    for ci = 1:c
        for fi = 1:f % fce
            for di = 1:d % target distance
                axh(ri, fi) = subplot(d,f,d*(fi-1) + di);grid on;hold on;
                %for di = 3 % target distance
                for li = 1:p % perturbation
                    trial_num = length(Data(ri,ci,fi,di,:,li));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,ci,fi,di,ti,li}))
                            continue;
                        end
                        
                        % chenguang: cannot run in current setting, do the
                        % folloing alteration: 2022-05-25
                        if (li == 1)
                            epoc_type = 2;
                        else
                            epoc_type = 1;
                        end
                        %
                        switch epoc_type
                            case 1
                                idx = find(Data{ri,ci,fi,di,ti,li}.Fp(1,:)~=0 & ...
                                    Data{ri,ci,fi,di,ti,li}.ts==4);  % pert at y
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if li == 1
                                    disp('ERROR: should use li == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
                                idx = (idx(1)-250):idx(end);
                                %idx = (idx(1)):idx(end);
                        end

                        time = t_step*(idx-idx(1));
                        idx_t{ri,ci,fi,di,ti,li} = idx;
                        time_t{ri,ci,fi,di,ti,li} = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,ci,fi,di,ti,li}.ox(1,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,ci,fi,di,ti,li}.f(1,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,ci,fi,di,ti,li}.Fp(1,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,ci,fi,di,ti,li}.v(1,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,ci,fi,di,ti,li}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,ci,fi,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,ci,fi,di,ti,li}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm force';

                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
            end
        end
    end
end
%xlim([0 0.7])
linkaxes(axh, 'xy');
%xlim([0 0.5]);
xlim([0 2])
sgtitle(titlestr);


% 3. Run the code to get the estimation matric
results1 = sys_id_pitt(Data,idx_t,time_t,1);
results2 = sys_id_pitt(Data,idx_t,time_t,2);
results3 = sys_id_pitt(Data,idx_t,time_t,3);

save('ss4229_4231_pulseParameterCompare.mat', 'results1', 'results2', 'results3', '-append');

% 4. Compare between matrics
