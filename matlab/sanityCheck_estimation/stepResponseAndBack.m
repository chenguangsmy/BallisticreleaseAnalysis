% A time-variant impedance 

% build a model with the variable impedance 
clear, clc, close all

F = 15;

t = 0:0.001:1;

m_n = 2; %[kg]
b_n = 20; %[Ns/m]
k_n = 500; %[N/m]

omega_n = sqrt(k_n/m_n); %[rad/s]
zita = b_n/(2*sqrt(k_n*m_n));
omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

x_n = (F/k_n)*(1-exp(-zita*omega_n*t)/(sqrt(1-zita^2)).*sin(omega_d*t+acos(zita)));
f_n = (x_n(end) - x_n)*k_n - [diff(x_n)/(1e-3), 0]*b_n + [diff(x_n,2)/1e-6, 0, 0]*m_n;
% f_n = ones(size(x_n))*F;
% f_n(1) = 0;

% something wrong with the analytical displacement? 
% solve the displacement, given the step function and the force input. 
f_input = ones(size(x_n))*F;
f_input(1) = 0; 
s = tf('s');
G = 1/(m_n*s^2+b_n*s+k_n);
x_pred = lsim(G,f_input,t);

% check the force, 
f_pred = (x_pred(end) - x_pred')*k_n + [diff(x_pred)'/(1e-3), 0]*b_n + [diff(x_pred,2)'/1e-6, 0, 0]*m_n;
plot(f_pred);

figure(); 
subplot(2,1,1); 
plot(t, f_n);
title('force');
subplot(2,1,2); 
plot(t, x_n); 
title('displacement');
sgtitle(' TI x and f analytical');

% use the tf to estimate the impedance 
f_input = f_n;
Ts = 1000;
data_est_UPs = iddata(x_n',f_input',1/Ts);
sysUP_s = tfest(data_est_UPs,2,0);
opt = predictOptions('InitialCondition','z');
[xp,~,~] = predict(sysUP_s,data_est_UPs,0,opt);
preddisp_interp_t = xp.OutputData;
[NUM_UPs,DEN_UPs] = tfdata(sysUP_s);
K_est_up_s = DEN_UPs{1}(3)/NUM_UPs{1}(3);
B_est_up_s = DEN_UPs{1}(2)/NUM_UPs{1}(3);
M_est_up_s = DEN_UPs{1}(1)/NUM_UPs{1}(3);
FIT_up_s = sysUP_s.Report.Fit.FitPercent;
