clear, clc, close all
set(0, 'DefaultLineLineWidth', 2);
set(groot,'defaultAxesFontSize',16);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background
set(groot,'defaultAxesBox','on'); % box on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');


s = tf('s');
t = 0:0.001:0.5;
f_amp = 12; %[N]
t_off = 0.1; %[s]
f_std = 0.015; %[N]
F = f_amp*exp(-(t-t_off).^2/(2*f_std^2));


m_n = 2; %[kg]
b_n = 20; %[Ns/m]
k_n = 500; %[N/m]

omega_n = sqrt(k_n/m_n); %[rad/s]
zita = b_n/(2*sqrt(k_n*m_n));
omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

G_n = 1/(m_n*s^2+b_n*s+k_n);
x_n = lsim(G_n,F,t);

for ii = 1:2
    m = ((0.4*(ii-1))/1+0.8)*m_n;
    b = b_n;
    k = k_n;

    omega_n = sqrt(k/m); %[rad/s]
    zita = b/(2*sqrt(k*m));
    omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

    G = 1/(m*s^2+b*s+k);
    x_m(ii,:) = lsim(G,F,t);
end

for ii = 1:2
    m = m_n;
    b = ((0.4*(ii-1))/1+0.8)*b_n;
    k = k_n;

    omega_n = sqrt(k/m); %[rad/s]
    zita = b/(2*sqrt(k*m));
    omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

    G = 1/(m*s^2+b*s+k);
    x_b(ii,:) = lsim(G,F,t);
end

for ii = 1:2
    m = m_n;
    b = b_n;
    k = ((0.4*(ii-1))/1+0.8)*k_n;

    omega_n = sqrt(k/m); %[rad/s]
    zita = b/(2*sqrt(k*m));
    omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

    G = 1/(m*s^2+b*s+k);
    x_k(ii,:) = lsim(G,F,t);
end

figure(),
subplot(4,1,1)
plot(t,F,'k'), grid on
ylabel('Force [N]')
subplot(4,1,2),
fill([t,fliplr(t)],[x_m(1,:),fliplr(x_m(2,:))]*1000,'r','FaceAlpha',0.5), hold on
plot(t,x_n*1000,'k'), grid on
subplot(4,1,3),
fill([t,fliplr(t)],[x_b(1,:),fliplr(x_b(2,:))]*1000,'b','FaceAlpha',0.5), hold on
plot(t,x_n*1000,'k'), grid on
ylabel('Displacement [mm]')
xlabel('Time [s]')
subplot(4,1,4),
fill([t,fliplr(t)],[x_k(1,:),fliplr(x_k(2,:))]*1000,'g','FaceAlpha',0.5), hold on
plot(t,x_n*1000,'k'), grid on
xlabel('Time [s]')
