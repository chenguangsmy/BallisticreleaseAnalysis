clear, clc, close all
set(0, 'DefaultLineLineWidth', 2);
set(groot,'defaultAxesFontSize',16);
set(0,'defaultfigurecolor',[1 1 1]); % white figure background
set(groot,'defaultAxesBox','on'); % box on
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');


s = tf('s');

load('step_force_WAM.mat');

m_n = 2; %[kg]
b_n = 20; %[Ns/m]
k_n = 500; %[N/m]

omega_n = sqrt(k_n/m_n); %[rad/s]
zita = b_n/(2*sqrt(k_n*m_n));
omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

G_n = 1/(m_n*s^2+b_n*s+k_n);
x_n = lsim(G_n,-(F(:,2)-max(F(:,2))),F(:,1));

for ii = 1:2
    m = ((0.4*(ii-1))/1+0.8)*m_n;
    b = b_n;
    k = k_n;

    omega_n = sqrt(k/m); %[rad/s]
    zita = b/(2*sqrt(k*m));
    omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

    G = 1/(m*s^2+b*s+k);
    x_m(ii,:) = lsim(G,-(F(:,2)-max(F(:,2))),F(:,1));
end

for ii = 1:2
    m = m_n;
    b = ((0.4*(ii-1))/1+0.8)*b_n;
    k = k_n;

    omega_n = sqrt(k/m); %[rad/s]
    zita = b/(2*sqrt(k*m));
    omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

    G = 1/(m*s^2+b*s+k);
    x_b(ii,:) = lsim(G,-(F(:,2)-max(F(:,2))),F(:,1));
end

for ii = 1:2
    m = m_n;
    b = b_n;
    k = ((0.4*(ii-1))/1+0.8)*k_n;

    omega_n = sqrt(k/m); %[rad/s]
    zita = b/(2*sqrt(k*m));
    omega_d = omega_n*sqrt(1-zita^2); %[rad/s]

    G = 1/(m*s^2+b*s+k);
    x_k(ii,:) = lsim(G,-(F(:,2)-max(F(:,2))),F(:,1));
end

figure(),
subplot(4,1,1)
plot(F(:,1),-(F(:,2)-max(F(:,2))),'k'), grid on
ylabel('Force [N]')
subplot(4,1,2),
fill([F(:,1)',fliplr(F(:,1)')],[x_m(1,:),fliplr(x_m(2,:))]*1000,'r','FaceAlpha',0.5), hold on
plot(F(:,1),x_n*1000,'k'), grid on
% ylabel('Displacement [mm]')
subplot(4,1,3),
fill([F(:,1)',fliplr(F(:,1)')],[x_b(1,:),fliplr(x_b(2,:))]*1000,'b','FaceAlpha',0.5), hold on
plot(F(:,1),x_n*1000,'k'), grid on
ylabel('Displacement [mm]')
subplot(4,1,4),
fill([F(:,1)',fliplr(F(:,1)')],[x_k(1,:),fliplr(x_k(2,:))]*1000,'g','FaceAlpha',0.5), hold on
plot(F(:,1),x_n*1000,'k'), grid on
% ylabel('Displacement [mm]')
xlabel('Time [s]')

%% Analytical Step Response
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



%Variable Stiffness
a = 0.5;
toff = 0.2;
k = k_n*(1-a*exp(-(t-toff).^2/(2*0.05^2)));
figure(),
plot(t,k), grid on
xlabel('Time [s]')
ylabel('Stiffness [N/m]')

i = 1;
for toff = [0.03 0.06 0.1 0.2]
    %Differential Equation
    f = @(t,x) [x(2);-(b_n/m_n)*x(2)-((k_n*(1-a*exp(-(t-toff)^2/(2*0.03^2))))/m_n)*x(1)+F/m_n];

    [tt,xx] = ode45(f,0:0.001:1,[0,0]);

    x_kvar(:,i) = xx(:,1);
    i = i + 1;
end

figure(),
subplot(4,1,1)
plot(t,x_n*1000,'k'), hold on
plot(tt,x_kvar(:,1)*1000), hold on
plot([0.03 0.03],[0 50],'b'), hold on
ylabel('Disp. [mm]')
grid on
subplot(4,1,2)
plot(t,x_n*1000,'k'), hold on
plot(tt,x_kvar(:,2)*1000), hold on
plot([0.06 0.06],[0 50],'b'), hold on
grid on
ylabel('Disp. [mm]')
subplot(4,1,3)
plot(t,x_n*1000,'k'), hold on
plot(tt,x_kvar(:,3)*1000), hold on
plot([0.1 0.1],[0 50],'b'), hold on
grid on
ylabel('Disp. [mm]')
subplot(4,1,4)
plot(t,x_n*1000,'k'), hold on
plot(tt,x_kvar(:,4)*1000), hold on
plot([0.2 0.2],[0 50],'b'), hold on
grid on
xlabel('Time [s]')
ylabel('Disp. [mm]')