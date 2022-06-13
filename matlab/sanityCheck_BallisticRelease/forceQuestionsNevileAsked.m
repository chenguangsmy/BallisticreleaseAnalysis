% Nevile was asking why the force after released at 200~400ms. There is a
% small blip at that time and I suspect that is just due to the static
% friction (I'll show it with align them with the velocity as well as the
% position). 


% Figure 1, an illustration of the 'force change' consistent with the zero
% velocity

dist_list = [2.5 5.0 7.5];
fce_list = [15 20 25];

% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4216_4226.mat', 'data');       % chenguang & Himanshu
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
if_subtract = 0;

epoc_type = 2;  % 1 perturb
                % 2 release
plot_type = 9;  % 1 displacement
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
pert_type = 1; % choose option [2 3 4]
axh = zeros(d,r);
xyi = 1;        % x, y

% fh = figure();
% pert_type = 2; colors = colors(4:end,:)
for ri = 1:2 % subj
    for di = 1:1 % direction
        fh = figure('Name', ['direction' num2str(di)]);
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci);grid on;hold on;
        %         axh(ri, di) = subplot(d,r,r*(di-1) + ri);grid on;hold on;
        %         axh(ri, ci) = subplot(1,1,1);grid on;hold on;
        for fi = 3%1:1 % target force
            for li = 2%1:1 % target distance
%                 axh(fi, li) = subplot(f,l,f*(li-1) + fi);grid on;hold on; %?
                for pi = 1%1:p % perturbation
                    trial_num = length(Data(ri,di,fi,li,:,pi));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,di,fi,li,ti,pi}))
                            continue;
                        end
                        
                        switch epoc_type
                            case 1
%                               idx = find(Data{ri,di,fi,li,ti,pi}.Fp(xyi,:)~=0 & Data{ri,di,fi,li,ti,pi}.ts==4);  % pert at y
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

                        time = t_step*(idx-idx(1));
                        axh(1) = subplot(3,1,1); hold on; grid on;
                            dat = Data{ri,di,fi,li,ti,pi}.f(xyi,idx);
                            plot(time, dat, 'Color', colors(4+li, :));
                            title('force');
                        axh(2) = subplot(3,1,2); hold on; grid on;
                            dat = Data{ri,di,fi,li,ti,pi}.ov(xyi,idx,1);
                            plot(time, dat, 'Color', colors(4+li, :));
                            title('velocity');
                        axh(3) = subplot(3,1,3); hold on; grid on;
                            dat = Data{ri,di,fi,li,ti,pi}.ox(xyi,idx,1);
                            plot(time, dat, 'Color', colors(4+li, :));
                            title('displacement');
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4+li, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
                title(['fce' num2str(fce_list(fi)) 'dist' num2str(dist_list(li))]);

            end
        end
        linkaxes(axh, 'x');
        xlim([0 1.0])
    end
end
linkaxes(axh, 'x');
xlim([0 1.0])
% xlim([0 0.5]);
% xlim([0 2])
% sgtitle(titlestr);




% Figure 2, an illustration of the time consistent. 