% After setting everything up in the experiment, (2022-04-20), 
% Chenguang is able to export the data to tidied up format. 

% Data is 6-D form

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For short-time processing, process each subject & each direction 
% individually
ss_list = [ 4076    4075    4074
            4077    4079    4078
            4081    4080    4082];
        
dir_i = 1;
subj_i = 1;
    
data = cell(1, 4, 3, 3, 15, 3);
fname = ['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss' ...
        num2str(min(ss_list(:))) '_' num2str(max(ss_list(:))) '.mat']

idx_sf = 1; % only choose suceed trials 
for fce_i = 1:size(ss_list,1)
    for tar_i = 1:3 % step perts
        if (isnan(ss_list(fce_i,tar_i)))
            continue;
        end
        fprintf('ss: %d', ss_list(fce_i, tar_i));
        ss_tmp = SessionScan(ss_list(fce_i, tar_i));
        ss_tmp.displayBlockCondition();
        celltmp = ss_tmp.export_as_formatted(0);    % 1 for show trials
        
        trials_num = size(celltmp,3);

        if trials_num>15
            data(subj_i,dir_i,fce_i,tar_i,:,:) = celltmp(idx_sf,1,1:15,:);
        else
            data(subj_i,dir_i,fce_i,tar_i,1:trials_num,:) = celltmp(idx_sf,1,:,:);
        end
    end

%     ss_tmp = SessionScan(ss_list(fce_i, 4)); % stoc pert
%     fprintf('ss: %d', ss_list(fce_i, 4));
%     ss_tmp.displayBlockCondition();
%     celltmp = ss_tmp.export_as_formatted(1);
%     for tar_i = 1:3 % as a session has 3 length
%         trials_num = size(celltmp,3);
%         data(1,1,fce_i,tar_i,1:5,3) = celltmp(idx_sf,tar_i,1:5,3);
% %         if trials_num>15
% %             data(1,1,fce_i,tar_i,1:15,3) = celltmp(idx_sf,tar_i,1:15,3);
% %         else
% %             data(1,1,fce_i,tar_i,1:trials_num,3) = celltmp(idx_sf,tar_i,:,3);
% %         end
%     end
end

save(fname, 'data')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display here 
load('data/processedData/ss4074_4082.mat', 'data');  % 6N perturbation on x, 
Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 
colors = colormap('lines');
r = size(Data, 1); % subj
d = size(Data, 2); % direction
f = size(Data, 3); % force
l = size(Data, 4); % distance
t = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 1;

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
axh = zeros(d,r);
xyi = 2;        % x, y

figure();
for ri = 1:r % subj
    for di = 1%:d % direction
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        %         axh(ri, di) = subplot(d,r,r*(di-1) + ri);grid on;hold on;
        %         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
        for fi = 1:f % target force
            for li = 1:3 % target distance
                axh(fi, li) = subplot(f,l,f*(li-1) + fi);grid on;hold on; %?
                for pi = 2%1:p % perturbation
                    trial_num = length(Data(ri,di,fi,li,:,pi));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,di,fi,li,ti,pi}))
                            continue;
                        end
                        
                        switch epoc_type
                            case 1
                                idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0 & ...
                                    Data{ri,di,fi,li,ti,pi}.ts==4);  % pert at y
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if pi == 1
                                    disp('ERROR: should use li == 2!!!');
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
                            case 8 % optotrak y
                                dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx);
                                titlestr = 'opto-position';
                                
                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4*(pi-1)+fi, :));
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
