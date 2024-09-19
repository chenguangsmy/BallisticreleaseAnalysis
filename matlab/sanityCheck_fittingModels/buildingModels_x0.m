% tidy up this model 

close all; clc; clear;
color_map = colormap('lines'); close all; 
load('data_tmp.mat')

Ts = mean(diff(data_eg.t));
data_eg.x = data_eg.x - data_eg.x(1);
data_eg.F_ivt = -(data_eg.F - mean(data_eg.F(1:200)));

% plot a force and displacement figure 
figure('unit', 'inch', 'position', [0 0 4 2]); 
axh = gca;
set(axh(1), 'fontsize', 14, 'linewidth', 1); hold on; 
yyaxis(axh, 'left')
plot(data_eg.t, data_eg.F_ivt, 'linewidth', 2);
%plot(data_eg.t, data_eg.F, 'linewidth', 2);
ylabel('Force (N)');
yyaxis(axh, 'right')
plot(data_eg.t, data_eg.x , 'linewidth', 2);
ylabel('Displacement (m)');
xlabel('Time (s)');
title('Raw Data');
% write a model and test it

%%  model the data to have K, B, and M
data_est_UP = iddata(data_eg.x',data_eg.F_ivt',Ts);
sysUP = tfest(data_est_UP,2,0);
opt = predictOptions('InitialCondition','z');
[yp,~,~] = predict(sysUP,data_est_UP,0,opt);
data_eg.x_pred = yp.OutputData;
[NUM_UP,DEN_UP] = tfdata(sysUP);
K_est   = DEN_UP{1}(3)/NUM_UP{1}(3);
B_est   = DEN_UP{1}(2)/NUM_UP{1}(3);
M_est   = DEN_UP{1}(1)/NUM_UP{1}(3);
% FIT_up_s = sysUP.Report.Fit.FitPercent;

fit_tmp = sqrt(sum((yp.OutputData - data_eg.x').^2)/length(yp.OutputData)) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(yp.OutputData));
FIT_up_ = 100 - fit_tmp*100;

up_omegan_s = sqrt(K_est/M_est);
up_dampr_s = B_est/(2*sqrt(K_est*M_est));
fitpercent(1) = FIT_up_;


% %% Use the model and 2nd order ode to predict displacement 
%% model_1: ddx = 1/m*(F-k(x-x0) - bx)

md1.x0 = data_eg.x(end); 
md1.k  = K_est;
md1.b  = B_est;
md1.m  = M_est;

freq = mean(1./diff(data_eg.t));


f = @(t,x) [x(2); -(md1.b/md1.m)*x(2) + (md1.k/md1.m)*(md1.x0 - x(1)) - interp1(data_eg.t, data_eg.F, t)/md1.m];

[tt,xx] = ode45(f,data_eg.t(1):1/freq:data_eg.t(end),[0.0,0]);

fit_tmp = sqrt(sum((xx(:,1) - data_eg.x').^2)/length(xx(:,1))) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(xx(:,1)));
FIT_up_ = 100 - fit_tmp*100;


% plot 
figure('unit', 'inch', 'position', [0 0 4 3]); 
axh = gca;
set(axh, 'fontsize', 14, 'linewidth', 1); hold on;
plot(data_eg.t, data_eg.x, 'b', 'linewidth', 2, 'Color', color_map(1,:), 'LineStyle', '--')
plot(data_eg.t, data_eg.x_pred, 'r', 'linewidth', 2, 'Color', color_map(6,:))
plot(tt, xx(:,1), 'r', 'linewidth', 2,  'Color', color_map(1,:))

ylabel('displacement (m)');
xlabel('time (s)');
legend('data', '${x(s)}/{F(s)} = 1/(ms^2 + bs + k)$', '$F = m\ddot{x} + b\dot{x} + k(x-x_0)$', 'Interpreter', 'latex'); 
% title(['fit ' num2str(fitpercent)]);
title(['prediction of frequency and time domains']);


%% model_2: ddx = 1/m*(F-k(x-x0) - bx) % ... need change ... assuming x0 changes after release;
% k1 = 300; % x0_1 = F/k1, x0_2 = x(end)
% ! use t_shift == 0 here is meaning less, but that means subject changes x0 instantaneously. 
t_shifts = [0.05 0.15 0.30]

x_pred = zeros(length(t_shifts), length(data_eg.x)); 
fit_scores_xc = zeros(length(t_shifts), 1); 

for t_i = 1:3
x0_ref  = zeros(size(data_eg.t));
k_ref   = zeros(size(data_eg.t));
% t_shift = 0.05; % assume a transcortical, visual feedback 
t_shift = t_shifts(t_i)
md1.k_dumb = 200; % a low value, not necessary to be the target stiffness 

x0_ref(data_eg.t<t_shift) = data_eg.F(1) / md1.k_dumb;
x0_ref(data_eg.t>=t_shift) = md1.x0;

k_ref(data_eg.t<0) = md1.k_dumb; 
k_ref(data_eg.t>=0) = md1.k_dumb; 
% 
f1 = @(t,x) [x(2); -(md1.b/md1.m)*x(2) - (md1.k_dumb/md1.m)*x(1) - interp1(data_eg.t, data_eg.F, t)/md1.m + (md1.k_dumb/md1.m)* interp1(data_eg.t, x0_ref, t)];

[tt,xx] = ode45(f1,data_eg.t(1):1/freq:data_eg.t(end),[0.0,0]);
x_pred_xc(t_i,:) = xx(:,1)';

fit_tmp = sqrt(sum((xx(:,1) - data_eg.x').^2)/length(xx(:,1))) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(xx(:,1)));
FIT_score = 100 - fit_tmp*100;
fit_scores_xc(t_i) = FIT_score;
end
xticks_xc(1) = x0_ref(1);
xticks_xc(2) = x0_ref(end);
% % find a better color-code


%%% model_3: ddx = 1/m*(F-k(x-x0) - bx) % ... need change ... assuming k changes after release; 
% k1 = 300; % k(1) = 
t_shifts = [0.05 0.15 0.30];

x_pred_kc = zeros(length(t_shifts), length(data_eg.x)); 
fit_scores_kc = zeros(length(t_shifts), 1); 


for t_i = 1:3
x0_ref  = zeros(size(data_eg.t)); 
k_ref   = zeros(size(data_eg.t)); 
% t_shift = 0.150; % assume a transcortical, visual feedback 
t_shift = t_shifts(t_i)
md1.k_dumb = 200; % a low value, not necessary to be the target stiffness 

x0_ref(data_eg.t<t_shift) = md1.x0;
x0_ref(data_eg.t>=t_shift) = md1.x0;

k_ref(data_eg.t<t_shift) = md1.k; 
k_ref(data_eg.t>=t_shift) = md1.k / 2; 
% 
f1 = @(t,x) [x(2); -(md1.b/md1.m)*x(2) - (interp1(data_eg.t, k_ref, t)/md1.m)*x(1) - interp1(data_eg.t, data_eg.F, t)/md1.m + (interp1(data_eg.t, k_ref, t)/md1.m)* md1.x0];

[tt,xx] = ode45(f1,data_eg.t(1):1/freq:data_eg.t(end),[0.0,0]);

x_pred_kc(t_i,:) = xx(:,1)';

fit_tmp = sqrt(sum((xx(:,1) - data_eg.x').^2)/length(xx(:,1))) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(xx(:,1)));
FIT_score = 100 - fit_tmp*100;
fit_scores_kc(t_i) = FIT_score;
end 

figure('unit', 'inch', 'position', [0 0 8 4]); 
axh(1) = subplot(1,2,1); 
set(axh(1), 'fontsize', 14, 'linewidth', 1); hold on;
plot(tt, data_eg.x, 'linewidth', 2,   'color', color_map(1,:), 'LineStyle','--'); 
plot(tt, x_pred_kc(1,:), 'linewidth', 2, 'color', color_map(2,:)); 
plot(tt, x_pred_kc(2,:), 'linewidth', 2, 'color', color_map(3,:)); 
plot(tt, x_pred_kc(3,:), 'linewidth', 2, 'color', color_map(4,:)); 

% legend('data', '$\tau = 0.05$', '$\tau = 0.15$', '$\tau = 0.30$', 'Interpreter', 'latex'); 
legend('data', 'spinal', 'cortical sensory', 'cortical visual'); 
title('Endpoint Stiffness Change');
ylabel('Displacement (m)');
xlabel('Time (s)');
% yticks(sort([0 x0_ref(1) x0_ref(end)]));
% yticklabels({'0', 'x_0_{Phd}', 'x_0_{Fhd}'})
fit_scores_kc
% find a better color-code

axh(2) = subplot(1,2,2); 
set(axh(2), 'fontsize', 14, 'linewidth', 1); hold on;
plot(tt, data_eg.x, 'linewidth', 2,   'color',   color_map(1,:), 'LineStyle','--'); 
plot(tt, x_pred_xc(1,:), 'linewidth', 2.5, 'color', color_map(2,:), 'LineStyle','-.'); 
plot(tt, x_pred_xc(2,:), 'linewidth', 2.5, 'color', color_map(3,:), 'LineStyle','-.'); 
plot(tt, x_pred_xc(3,:), 'linewidth', 2.5, 'color', color_map(4,:), 'LineStyle','-.'); 
% plot(tt, x_pred_xc(4,:), 'linewidth', 2, 'color', [0.5 0.5 0.5]*8/5); 
yline([xticks_xc(1), xticks_xc(2)]);
% legend('origin', '$\tau = 0$', '$\tau = 0.05$', '$\tau = 0.15$', '$\tau = 0.30$', 'Interpreter', 'latex'); 
legend('data', 'spinal', 'cortical sensory', 'cortical visual');
title('Equilibrium Position Change');
xlabel('Time (s)');

linkaxes(axh, 'xy')
% yticks(sort([0 0.05 0.10 x0_ref(1) x0_ref(end)]));
% yticklabels({'0', '0.05', 'x_0_Phd', 'x_0_Fhd', '0.10'})
yticks(sort([0 xticks_xc(1) xticks_xc(2)]));
yticklabels({'0', 'x_0_{Phd}', 'x_0_{Fhd}'})
% legend('data', '${x(s)}/{F(s)} = 1/(ms^2 + bs + k)$', '$F = m\ddot{x} + b\dot{x} + k(x-x_0)$', 'Interpreter', 'latex'); 
fit_scores_xc


%% 
% across trials 
