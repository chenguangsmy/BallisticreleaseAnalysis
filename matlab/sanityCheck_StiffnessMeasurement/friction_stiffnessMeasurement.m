% use smallest and sloest movement as possible;
% use the before-perturb and after-perturb force and displacement as
% measurement

% 160N/m
%load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3104.mat','data') 
% 320N/m
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss3103.mat', 'data')

Freq = 500;
t_step = 1/500;
fh = figure(); 
colors = colormap('lines');
r = size(data, 1);
c = 3; % F_cmd, F_csd vs x
l = size(data, 2);
idx_last = 200;

k_est = [];
for ri = 1:r % sigmas
    for ci = 1:c % cmd vs csd
        %axh(ri, ci) = subplot(r,c,c*(ri-1) + ci); 
        axh(ri, ci) = subplot(c,r,r*(ci-1) + ri); 
        hold on;
        for mi = 1:l % magnitude of the lines 
            trial_num = length(data{ri,mi});
            for ti = 1:trial_num % each trial
                if (isempty(data{ri,mi}{ti}))
                    continue;
                end
                idx = find(data{ri,mi}{ti}.Fp(2,:)~=0);
                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                time = t_step*(idx-idx(1));
                switch ci
                    case 1
                        dat = -data{ri,mi}{ti}.Fp(2,idx);
                    case 2
                        dat = data{ri,mi}{ti}.f(2,idx);
                    case 3
                        dat = data{ri,mi}{ti}.x(2,idx);
                end
                plot(time, dat, 'Color', colors(mi, :));
                
                if ci == 2 % measuring stiffness
                    idx_bef = idx(1:100);        
                    idx_aft = idx(end-199:end); 
                    plot(time([1:100]), dat([1:100]), 'r.');
                    plot(time(end-199:end), dat(end-199:end), 'r.'); 
                    dF = mean(data{ri,mi}{ti}.f(2,idx_aft)) - mean(data{ri,mi}{ti}.f(2,idx_bef));
                    dx = mean(data{ri,mi}{ti}.x(2,idx_aft)) - mean(data{ri,mi}{ti}.x(2,idx_bef));
                    k  = dF/dx
                    k_est = [k_est k];
                end
            end
        end
    end
end
linkaxes(axh(:), 'x');
linkaxes(axh(:, 1), 'y');
ylim(axh(1,2), 4+[-6 6]);
linkaxes(axh(:, 2), 'y');
%ylim(axh(1,2), 15+[-1 8]);
for ri = 1:r
ylabel(axh(ri,1), 'F_{cmd} (N)');
ylabel(axh(ri,2), 'F_{csd} (N)');
ylabel(axh(ri,3), 'x_{csd} (m)');
end 
for ci = 1:c
xlabel(axh(end,ci), 'time (s)');
end
for ri = 1:r
    for ci = 1:c
        grid(axh(ri, ci));
    end
end
sgtitle('multiple bell-shaped perturbations'); 


