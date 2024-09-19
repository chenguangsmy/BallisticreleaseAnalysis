% tidy up this file 

close all; clc; clear;
load('data_tmp.mat')

Ts = mean(diff(data_eg.t));
disp_interp_t = data_eg.x - data_eg.x(1);
data_eg.x = data_eg.x - data_eg.x(1);
force_interp_t = -(data_eg.F - mean(data_eg.F(1:200)));

color_map = colormap('lines');

figure('unit', 'inch', 'position', [0 0 4 2]); 
axh = gca;
set(axh(1), 'fontsize', 14, 'linewidth', 1); hold on; 
yyaxis(axh, 'left')
set(gca, 'YColor', color_map(2,:))
plot(data_eg.t, force_interp_t, 'linewidth', 2, 'Color',color_map(2,:));
ylabel('force (N)');
yyaxis(axh, 'right')
set(gca, 'YColor', color_map(1,:))
plot(data_eg.t, data_eg.x , 'linewidth', 2, 'Color',color_map(1,:));
ylabel('displacement (m)');
xlabel('time (s)');
title('Raw Data');
% write a model and test it
i = 1

%% model 1
data_est_UP = iddata(disp_interp_t(i,:)',force_interp_t(i,:)',Ts);
sysUP = tfest(data_est_UP,0,0);
[NUM_UP,DEN_UP] = tfdata(sysUP);
K_est_up(i) = DEN_UP{1}(1);
opt = predictOptions('InitialCondition','z');
[yp,~,~] = predict(sysUP,data_est_UP,0,opt);
preddisp_interp_t(1,:) = yp.OutputData;
% fitpercent(1) = sysUP.Report.Fit.FitPercent;
fit_tmp = sqrt(sum((yp.OutputData - data_eg.x').^2)/length(yp.OutputData)) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(yp.OutputData));
FIT_up_(i) = 100 - fit_tmp*100;
fitpercent(1) = FIT_up_(i);

%% model 2
data_est_UP = iddata(disp_interp_t(i,:)',force_interp_t(i,:)',Ts);
sysUP = tfest(data_est_UP,1,0);
[NUM_UP,DEN_UP] = tfdata(sysUP);
K_est_up(i) = DEN_UP{1}(1);
B_est_up(i) = DEN_UP{1}(2);
opt = predictOptions('InitialCondition','z');
[yp,~,~] = predict(sysUP,data_est_UP,0,opt);
preddisp_interp_t(2,:) = yp.OutputData;

sysUP = tfest(data_est_UP,1,0);
% fitpercent(2) = sysUP.Report.Fit.FitPercent;
fit_tmp = sqrt(sum((yp.OutputData - data_eg.x').^2)/length(yp.OutputData)) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(yp.OutputData));
FIT_up_(i) = 100 - fit_tmp*100;
fitpercent(2) = FIT_up_(i);

%% model 3
data_est_UP = iddata(disp_interp_t(i,:)',force_interp_t(i,:)',Ts);
sysUP = tfest(data_est_UP,2,0);
opt = predictOptions('InitialCondition','z');
[yp,~,~] = predict(sysUP,data_est_UP,0,opt);
preddisp_interp_t(3,:) = yp.OutputData;
[NUM_UP,DEN_UP] = tfdata(sysUP);
K_est_up(i) = DEN_UP{1}(3)/NUM_UP{1}(3);
B_est_up_s(i) = DEN_UP{1}(2)/NUM_UP{1}(3);
M_est_up_s(i) = DEN_UP{1}(1)/NUM_UP{1}(3);
% FIT_up_s(i) = sysUP.Report.Fit.FitPercent;

fit_tmp = sqrt(sum((yp.OutputData - data_eg.x').^2)/length(yp.OutputData)) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(yp.OutputData));
FIT_up_(i) = 100 - fit_tmp*100;

up_omegan_s(i) = sqrt(K_est_up(i)/M_est_up_s(i));
up_dampr_s(i) = B_est_up_s(i)/(2*sqrt(K_est_up(i)*M_est_up_s(i)));
% fitpercent(3) = sysUP.Report.Fit.FitPercent;
fitpercent(3) = FIT_up_(i);


%%
figure('unit', 'inch', 'position', [0 0 4 6]); 
axh(1) = subplot(3,1,1); 
set(axh(1), 'fontsize', 14, 'linewidth', 1); hold on;
plot(data_eg.t, disp_interp_t(i,:), 'b', 'linewidth', 2, 'Color', color_map(1,:), 'LineStyle','--')
plot(data_eg.t, preddisp_interp_t(1,:), 'r', 'linewidth', 2, 'Color', color_map(1,:), 'LineStyle','-');
legend('data', '$F=kx$',  'Interpreter', 'latex'); 
% ylabel('displacement (m)');
% xlabel('time (s)');
% title({'Model Prediction'; ['F=kx, fit' num2str(sysUP_s.Report.Fit.FitPercent)]});
title(['fit ' num2str(fitpercent(1))]);

axh(2) = subplot(3,1,2); 
set(axh(2), 'fontsize', 14, 'linewidth', 1); hold on;
plot(data_eg.t, disp_interp_t(i,:), 'b', 'linewidth', 2, 'Color', color_map(1,:), 'LineStyle','--')
plot(data_eg.t, preddisp_interp_t(2,:), 'r', 'linewidth', 2, 'Color', color_map(1,:), 'LineStyle','-');
% legend('data', 'prediction'); 
legend('data', '$F = b\dot{x} + kx$',  'Interpreter', 'latex'); 
ylabel('displacement (m)');
% xlabel('time (s)');
% title({'Model Prediction'; ['$F = b\dot{x} + kx$  fit' num2str(sysUP_s.Report.Fit.FitPercent)]}, 'Interpreter', 'latex')
% title('$F = b\dot{x} + kx$', 'Interpreter', 'latex')
title(['fit ' num2str(fitpercent(2))]);

axh(3) = subplot(3,1,3); 
set(axh(3), 'fontsize', 14, 'linewidth', 1); hold on;
plot(data_eg.t, disp_interp_t(i,:), 'b', 'linewidth', 2, 'Color', color_map(1,:), 'LineStyle','--')
plot(data_eg.t, preddisp_interp_t(3,:), 'r', 'linewidth', 2, 'Color', color_map(1,:), 'LineStyle','-');
% legend('data', 'prediction'); 
legend('data', '$F = m\ddot{x} + b\dot{x} + kx$',  'Interpreter', 'latex'); 
% ylabel('displacement (m)');
xlabel('time (s)');
% title({'Model Prediction'; ['$F = m\ddot{x} + b\dot{x} + kx$  fit' num2str(sysUP_s.Report.Fit.FitPercent)]}, 'Interpreter', 'latex')
% title('$F = m\ddot{x} + b\dot{x} + kx$', 'Interpreter', 'latex')
title(['fit ' num2str(fitpercent(3))]);

sgtitle('Compare Dynamic Models');


% % % %% model 4
% % % % clc; close all;
% % % % a delayed model with state-space
% % % % Define the parameters
% % % m = M_est_up_s;  % Mass (kg)
% % % b = B_est_up_s;  % Damping coefficient (Ns/m)
% % % k = K_est_up;  % Spring constant (N/m)
% % % % tau = -0.05; % Time delay (s)
% % % tau = -0.05;
% % % 
% % % params.m = m;
% % % params.b = b;
% % % params.k = k;
% % % params.tau = tau;
% % % save('data_tmp.mat', 'params', '-append');
% % % 
% % % % Define the state-space matrices without delay
% % % A = [0 1; -k/m -b/m];
% % % B = [0; 1/m];
% % % C = [1 0];
% % % D = 0;
% % % 
% % % % Define the delayed state matrix A_d
% % % A_d = [0 0; -k/m/4 0];
% % % % A_d = [0 0; 0 0];
% % % B_d = [0; 0];
% % % C_d = [0 0];
% % % D_d = 0; 
% % % 
% % % % Create the state-space model (without delay)
% % % sys = ss(A, B, C, D);
% % % 
% % % % Introduce the time delay in the system using `pade` approximation
% % % sys_delay = ss(A_d, B_d, C, D);
% % % 
% % % DelayT = struct('delay', -tau, 'a', A_d, 'b', B_d, 'c', C_d, 'd', D_d);
% % % 
% % % sys_total = delayss(A, B, C, D, DelayT)
% % % % Combine both systems (using a linear time-delay system model)
% % % 
% % % % Display the state-space model
% % % % disp('State-Space Model with Time Delay:')
% % % % sys_total
% % % 
% % % % step(sys_total)
% % % 
% % % %% Simulate the response to an input force
% % % figure;
% % % % step(sys_total);
% % % % [yp,~,~] = predict(sys_total,data_est_UPs,0,opt);
% % % 
% % % output = lsim(sys_total, force_interp_t, data_eg.t);
% % % % output = lsim(sys_delay, force_interp_t, data_eg.t);
% % % plot(output)
% % % title('Step Response of the System with Time Delay');
% % % grid on;
% % % 
% % % fold = 2;
% % % force_interp_t_long = [force_interp_t repmat(force_interp_t(end),1, (fold-1)*length(force_interp_t))];
% % % data_eg.t_long = [data_eg.t data_eg.t - data_eg.t(1) + data_eg.t(end) +  data_eg.t(2) - data_eg.t(1)]; % think how to make fold-able
% % % output = lsim(sys_total, force_interp_t_long, data_eg.t_long);
% % % % output = lsim(sys_delay, force_interp_t, data_eg.t);
% % % plot(data_eg.t_long, output)
% % % title('Step Response of the System with Time Delay');
% % % grid on;
% % % 
% % % %%
