% tidy up this model 

close all; clc; clear;
load('bdata_tmp.mat')

fit_all = zeros(length(bdata_eg), 7);

IFPLOT = 0; 
for trial_i = 1:length(bdata_eg)
clearvars x_pred_xc

Ts = mean(diff(bdata_eg{trial_i}.t));
bdata_eg{trial_i}.x = bdata_eg{trial_i}.x - bdata_eg{trial_i}.x(1);
bdata_eg{trial_i}.F_ivt = -(bdata_eg{trial_i}.F - mean(bdata_eg{trial_i}.F(1:200)));

if (IFPLOT)
% plot a force and displacement figure 
figure('unit', 'inch', 'position', [0 0 4 2]); 
axh = gca;
set(axh(1), 'fontsize', 14, 'linewidth', 1); hold on; 
yyaxis(axh, 'left')
% plot(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F_ivt, 'linewidth', 2);
plot(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, 'linewidth', 2);
ylabel('force (N)');
yyaxis(axh, 'right')
plot(bdata_eg{trial_i}.t, bdata_eg{trial_i}.x , 'linewidth', 2);
ylabel('displacement (m)');
xlabel('time (s)');
title('Raw Data');
% write a model and test it
end

%%  model the data to have K, B, and M
data_est_UP = iddata(bdata_eg{trial_i}.x',bdata_eg{trial_i}.F_ivt',Ts);
sysUP = tfest(data_est_UP,2,0);
opt = predictOptions('InitialCondition','z');
[yp,~,~] = predict(sysUP,data_est_UP,0,opt);
bdata_eg{trial_i}.x_pred = yp.OutputData;
[NUM_UP,DEN_UP] = tfdata(sysUP);
K_est   = DEN_UP{1}(3)/NUM_UP{1}(3);
B_est   = DEN_UP{1}(2)/NUM_UP{1}(3);
M_est   = DEN_UP{1}(1)/NUM_UP{1}(3);
% FIT_up_s = sysUP.Report.Fit.FitPercent;

fit_tmp = sqrt(sum((yp.OutputData - bdata_eg{trial_i}.x').^2)/length(yp.OutputData)) / ...
             sqrt(sum((bdata_eg{trial_i}.x' - mean(bdata_eg{trial_i}.x)).^2)/length(yp.OutputData));
FIT_up_ = 100 - fit_tmp*100;

up_omegan_s = sqrt(K_est/M_est);
up_dampr_s = B_est/(2*sqrt(K_est*M_est));
fitpercent(1) = FIT_up_;

% sgtitle('Compare Dynamical Models');


% %% Use the model and 2nd order ode to predict displacement 
%% model_1: ddx = 1/m*(F-k(x-x0) - bx)

md1.x0 = bdata_eg{trial_i}.x(end); 
md1.k  = K_est;
md1.b  = B_est;
md1.m  = M_est;

freq = mean(1./diff(bdata_eg{trial_i}.t));


f = @(t,x) [x(2); -(md1.b/md1.m)*x(2) + (md1.k/md1.m)*(md1.x0 - x(1)) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m];

[tt,xx] = ode45(f,bdata_eg{trial_i}.t(1):1/freq:bdata_eg{trial_i}.t(end),[0.0,0]);

fit_tmp = sqrt(sum((xx(:,1) - bdata_eg{trial_i}.x').^2)/length(xx(:,1))) / ...
             sqrt(sum((bdata_eg{trial_i}.x' - mean(bdata_eg{trial_i}.x)).^2)/length(xx(:,1)));
FIT_up_ = 100 - fit_tmp*100;

fit_all(trial_i, 1) = FIT_up_;
% plot 
if (IFPLOT)
figure('unit', 'inch', 'position', [0 0 4 3]); 
axh = gca;
set(axh, 'fontsize', 14, 'linewidth', 1); hold on;
plot(bdata_eg{trial_i}.t, bdata_eg{trial_i}.x, 'b', 'linewidth', 2)
plot(bdata_eg{trial_i}.t, bdata_eg{trial_i}.x_pred, 'r', 'linewidth', 2);
plot(tt, xx(:,1), 'r', 'linewidth', 2, 'LineStyle', '--'); 

ylabel('displacement (m)');
xlabel('time (s)');
legend('data', '${x(s)}/{F(s)} = 1/(ms^2 + bs + k)$', '$F = m\ddot{x} + b\dot{x} + k(x-x_0)$', 'Interpreter', 'latex'); 
% title(['fit ' num2str(fitpercent)]);
title(['prediction of frequency and time domains']);
end

%% model_2: ddx = 1/m*(F-k(x-x0) - bx) % ... need change ... assuming x0 changes after release;
% k1 = 300; % x0_1 = F/k1, x0_2 = x(end)
% ! use t_shift == 0 here is meaning less, but that means subject changes x0 instantaneously. 
t_shifts = [0.05 0.15 0.30]

x_pred = zeros(length(t_shifts), length(bdata_eg{trial_i}.x)); 
fit_scores_xc = zeros(length(t_shifts), 1); 

for t_i = 1:3
x0_ref  = zeros(size(bdata_eg{trial_i}.t));
k_ref   = zeros(size(bdata_eg{trial_i}.t));
% t_shift = 0.05; % assume a transcortical, visual feedback 
t_shift = t_shifts(t_i)
md1.k_dumb = 200; % a low value, not necessary to be the target stiffness 

x0_ref(bdata_eg{trial_i}.t<t_shift) = bdata_eg{trial_i}.F(1) / md1.k_dumb;
x0_ref(bdata_eg{trial_i}.t>=t_shift) = md1.x0;

k_ref(bdata_eg{trial_i}.t<0) = md1.k_dumb; 
k_ref(bdata_eg{trial_i}.t>=0) = md1.k_dumb; 
% f = @(t,x) [x(2); -(md1.b/md1.m)*x(2) + (md1.k/md1.m)*(md1.x0 - x(1)) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m];
% f1 = @(t,x) [x(2); -(md1.b/md1.m)*x(2) - (md1.k/md1.m)*x(1) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m + (md1.k/md1.m)* interp1(bdata_eg{trial_i}.t, x0_ref, t)];
f1 = @(t,x) [x(2); -(md1.b/md1.m)*x(2) - (md1.k_dumb/md1.m)*x(1) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m + (md1.k_dumb/md1.m)* interp1(bdata_eg{trial_i}.t, x0_ref, t)];

[tt,xx] = ode45(f1,bdata_eg{trial_i}.t(1):1/freq:bdata_eg{trial_i}.t(end),[0.0,0]);
x_pred_xc(t_i,:) = xx(:,1)';

fit_tmp = sqrt(sum((xx(:,1) - bdata_eg{trial_i}.x').^2)/length(xx(:,1))) / ...
             sqrt(sum((bdata_eg{trial_i}.x' - mean(bdata_eg{trial_i}.x)).^2)/length(xx(:,1)));
FIT_score = 100 - fit_tmp*100;
fit_scores_xc(t_i) = FIT_score;
end
xticks_xc(1) = x0_ref(1);
xticks_xc(2) = x0_ref(end);
% % find a better color-code
% figure('unit', 'inch', 'position', [0 0 5 4]); 
% axh = gca; 
% set(axh, 'fontsize', 14, 'linewidth', 1); hold on;
% plot(tt, bdata_eg{trial_i}.x, 'linewidth', 2,   'color',   [0.5 0.5 0.5], 'LineStyle','--'); 
% plot(tt, x_pred_xc(1,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*2/5); 
% plot(tt, x_pred_xc(2,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*4/5); 
% plot(tt, x_pred_xc(3,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*6/5); 
% % plot(tt, x_pred_xc(4,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*8/5); 
% yline([x0_ref(1), x0_ref(end)]);
% % legend('origin', '$\tau = 0$', '$\tau = 0.05$', '$\tau = 0.15$', '$\tau = 0.30$', 'Interpreter', 'latex'); 
% 
% legend('origin', '$\tau = 0.05$', '$\tau = 0.15$', '$\tau = 0.30$', 'Interpreter', 'latex'); 
% title('Subjects Changes x_0 after Release');
% 
% % yticks(sort([0 0.05 0.10 x0_ref(1) x0_ref(end)]));
% % yticklabels({'0', '0.05', 'x_0_Phd', 'x_0_Fhd', '0.10'})
% yticks(sort([0 x0_ref(1) x0_ref(end)]));
% yticklabels({'0', 'x_0_{Phd}', 'x_0_{Fhd}'})
% % legend('data', '${x(s)}/{F(s)} = 1/(ms^2 + bs + k)$', '$F = m\ddot{x} + b\dot{x} + k(x-x_0)$', 'Interpreter', 'latex'); 


%%% model_3: ddx = 1/m*(F-k(x-x0) - bx) % ... need change ... assuming k changes after release; 
% k1 = 300; % k(1) = 
t_shifts = [0.05 0.15 0.30];

x_pred_kc = zeros(length(t_shifts), length(bdata_eg{trial_i}.x)); 
fit_scores_kc = zeros(length(t_shifts), 1); 


for t_i = 1:3
x0_ref  = zeros(size(bdata_eg{trial_i}.t)); 
k_ref   = zeros(size(bdata_eg{trial_i}.t)); 
% t_shift = 0.150; % assume a transcortical, visual feedback 
t_shift = t_shifts(t_i)
md1.k_dumb = 200; % a low value, not necessary to be the target stiffness 

x0_ref(bdata_eg{trial_i}.t<t_shift) = md1.x0;
x0_ref(bdata_eg{trial_i}.t>=t_shift) = md1.x0;

k_ref(bdata_eg{trial_i}.t<t_shift) = md1.k; 
k_ref(bdata_eg{trial_i}.t>=t_shift) = md1.k / 2; 
% f = @(t,x) [x(2); -(md1.b/md1.m)*x(2) + (md1.k/md1.m)*(md1.x0 - x(1)) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m];
% f1 = @(t,x) [x(2); -(md1.b/md1.m)*x(2) - (md1.k/md1.m)*x(1) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m + (md1.k/md1.m)* interp1(bdata_eg{trial_i}.t, x0_ref, t)];
% f1 = @(t,x) [x(2); -(md1.b/md1.m)*x(2) - (md1.k_dumb/md1.m)*x(1) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m + (md1.k_dumb/md1.m)* interp1(bdata_eg{trial_i}.t, x0_ref, t)];
f1 = @(t,x) [x(2); -(md1.b/md1.m)*x(2) - (interp1(bdata_eg{trial_i}.t, k_ref, t)/md1.m)*x(1) - interp1(bdata_eg{trial_i}.t, bdata_eg{trial_i}.F, t)/md1.m + (interp1(bdata_eg{trial_i}.t, k_ref, t)/md1.m)* md1.x0];

[tt,xx] = ode45(f1,bdata_eg{trial_i}.t(1):1/freq:bdata_eg{trial_i}.t(end),[0.0,0]);

x_pred_kc(t_i,:) = xx(:,1)';

fit_tmp = sqrt(sum((xx(:,1) - bdata_eg{trial_i}.x').^2)/length(xx(:,1))) / ...
             sqrt(sum((bdata_eg{trial_i}.x' - mean(bdata_eg{trial_i}.x)).^2)/length(xx(:,1)));
FIT_score = 100 - fit_tmp*100;
fit_scores_kc(t_i) = FIT_score;
end 

if(IFPLOT)
figure('unit', 'inch', 'position', [0 0 8 4]); 
axh(1) = subplot(1,2,1); 
set(axh(1), 'fontsize', 14, 'linewidth', 1); hold on;
plot(tt, bdata_eg{trial_i}.x, 'linewidth', 2,   'color',   [0.5 0.5 0.5], 'LineStyle','--'); 
plot(tt, x_pred_kc(1,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*2/5); 
plot(tt, x_pred_kc(2,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*4/5); 
plot(tt, x_pred_kc(3,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*6/5); 

legend('origin', '$\tau = 0.05$', '$\tau = 0.15$', '$\tau = 0.30$', 'Interpreter', 'latex'); 
title('Endpoint Stiffness Change');
% yticks(sort([0 x0_ref(1) x0_ref(end)]));
% yticklabels({'0', 'x_0_{Phd}', 'x_0_{Fhd}'})
fit_scores_kc
% find a better color-code

axh(2) = subplot(1,2,2); 
set(axh(2), 'fontsize', 14, 'linewidth', 1); hold on;
plot(tt, bdata_eg{trial_i}.x, 'linewidth', 2,   'color',   [0.5 0.5 0.5], 'LineStyle','--'); 
plot(tt, x_pred_xc(1,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*2/5); 
plot(tt, x_pred_xc(2,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*4/5); 
plot(tt, x_pred_xc(3,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*6/5); 
% plot(tt, x_pred_xc(4,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*8/5); 
yline([xticks_xc(1), xticks_xc(2)]);
% legend('origin', '$\tau = 0$', '$\tau = 0.05$', '$\tau = 0.15$', '$\tau = 0.30$', 'Interpreter', 'latex'); 

legend('origin', '$\tau = 0.05$', '$\tau = 0.15$', '$\tau = 0.30$', 'Interpreter', 'latex'); 
title('Equilibrium Position Change');

linkaxes(axh, 'xy')
% yticks(sort([0 0.05 0.10 x0_ref(1) x0_ref(end)]));
% yticklabels({'0', '0.05', 'x_0_Phd', 'x_0_Fhd', '0.10'})
yticks(sort([0 xticks_xc(1) xticks_xc(2)]));
yticklabels({'0', 'x_0_{Phd}', 'x_0_{Fhd}'})
end
% legend('data', '${x(s)}/{F(s)} = 1/(ms^2 + bs + k)$', '$F = m\ddot{x} + b\dot{x} + k(x-x_0)$', 'Interpreter', 'latex'); 
fit_scores_xc

fit_all(trial_i, 2:7) = [fit_scores_kc' fit_scores_xc'];
end

%% box plot 
color_map = colormap('lines');

figure('unit', 'inch', 'position', [0 0 8 4]); 
axh(1) = subplot(1,2,1); hold on; 
axh(2) = subplot(1,2,2); hold on; 
boxplot(axh(1),fit_all(:,[1 2 3 4]), 'Colors', [color_map(1,:); color_map(2,:); color_map(3,:); color_map(4,:); ], 'PlotStyle','compact')
boxplot(axh(2),fit_all(:,[1 5 6 7]), 'Colors', [color_map(1,:); color_map(2,:); color_map(3,:); color_map(4,:); ], 'PlotStyle','compact')

set(axh(1), 'fontsize', 14, 'linewidth', 1); hold on;
set(axh(2), 'fontsize', 14, 'linewidth', 1); hold on;


xticks(axh(1), [1 2 3 4])
xticklabels(axh(1), {'TI', 'Spinal', 'Cortical Sensory' 'Cortical Visual'})
ylabel(axh(1), 'Fitting (%)')
title(axh(1), 'Endpoint Stiffness Change')

xticks(axh(2), [1 2 3 4])
xticklabels(axh(2), {'TI', 'Spinal', 'Cortical Sensory' 'Cortical Visual'})
ylabel(axh(2), 'Fitting (%)')
title(axh(2), 'Equilibrium Position Change')

pvals = zeros(1,6);
for test_i = 1:6
    [H, P] = ttest2(fit_all(:,1), fit_all(:,test_i+1));
    pvals(test_i) = P;
end

linkaxes(axh, 'y');

pvals

% add sig ** 
ylim([-20 150])
% left 
subplot(axh(1))
for i = 1:3
    if pvals(i) < 0.05
        pmark = '*'
    else
        continue
    end
    line([1 i+1], [100 100]+(i-1)*15, 'linewidth', 1, 'color', 'k')
    text(mean([1 i+1]), 100 + (i-1)* 15 + 5, '*', 'color', 'k', 'FontSize', 15)
end

% right 
subplot(axh(2))
for i = 1:3
    if pvals(3+i) < 0.05
        pmark = '*'
    else
        continue
    end
    line([1 i+1], [100 100]+(i-1)*15, 'linewidth', 1, 'color', 'k')
    text(mean([1 i+1]), 100 + (i-1)* 15 + 5, '*', 'color', 'k', 'FontSize', 15)
end

%% 
% across trials 
