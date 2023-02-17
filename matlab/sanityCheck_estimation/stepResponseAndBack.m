% Step response and calculate its force back
% Question: Why cannot I get the original force?

%% build a model with the actual step response
clear, clc,
% close all

F = 15;

freq = 10000;
t = 0:1/freq:1;

m_n = 2; %[kg]
b_n = 20; %[Ns/m]
k_n = 500; %[N/m]

omega_n = sqrt(k_n/m_n); %[rad/s]
zita = b_n/(2*sqrt(k_n*m_n));
omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

% solve the displacement, given the step function and the force input. 
f_input = ones(size(t))*F;
f_input(1) = 0; 
s = tf('s');
G = 1/(m_n*s^2+b_n*s+k_n);
x_pred = lsim(G,f_input,t);

% check the force, 
f_pred = (x_pred'-x_pred(1))*k_n + [diff(x_pred)'/(1/freq), 0]*b_n + [diff(x_pred,2)'/(1/freq^2), 0, 0]*m_n;
figure();
subplot(2,1,1);
hold on;
plot(t, f_input);
plot(t, f_pred);
legend('input force', 'reconstructed force');
subplot(2,1,2);
plot(t, x_pred);
sgtitle(['freq' num2str(freq)]);

% conclusion: the step-response is able to predict step-force... 

%% build a model with the theoratical impulse reseponse 
x_pred_ir = impulse(G,t);
figure(); hold on;
plot(t, x_pred, 'k--');
plot(t, x_pred_ir, 'b.-');
legend('actual impulse', 'theoratical impulse')