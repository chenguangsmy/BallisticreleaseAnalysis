%% A time-invariant impedance 
clear, clc, close all

F       = 15;
freq    = 1000;
m_n     = 2; %[kg]
b_n     = 20; %[Ns/m]
k_n     = 500; %[N/m]
omega_n = sqrt(k_n/m_n); %[rad/s]
zita    = b_n/(2*sqrt(k_n*m_n));
omega_d = omega_n*sqrt(1-zita^2); %[rad/s]


t   = 0:1/freq:2;
x_n = (F/k_n)*(1-exp(-zita*omega_n*t)/(sqrt(1-zita^2)).*sin(omega_d*t+acos(zita)));
f_n = (x_n-x_n(1))*k_n + [diff(x_n)/(1/freq), 0]*b_n + [diff(x_n,2)/(1/freq^2), 0, 0]*m_n;
f_n0= ones(size(x_n))*F; % force input
f_n0(1) = 0;

% solve the displacement, given the step function and the force input. 
f_input = ones(size(x_n))*F;
f_input(1) = 0; 
s = tf('s');
G = 1/(m_n*s^2+b_n*s+k_n);
x_pred = lsim(G,f_input,t);

% check the force, 
f_pred = (x_pred'-x_pred(1))*k_n + [diff(x_pred)'/(1/freq), 0]*b_n + [diff(x_pred,2)'/(1/freq^2), 0, 0]*m_n;

figure(); 
subplot(2,1,1); 
hold on;
plot(t, f_n0, 'k--', 'linewidth', 1);
plot(t, f_n, 'b', 'linewidth', 1); 
plot(t, f_pred, 'r', 'linewidth', 1);
ylim([0 20]);
legend('F_{step}', 'F_{analytical}', 'F_{emperical}')
title('force');
subplot(2,1,2); 
hold on;
plot(t, x_n, 'b', 'linewidth', 1); 
plot(t, x_pred, 'r', 'linewidth', 1);
legend('x_{analytical}', 'x_{emperical}')
title('displacement');
sgtitle('analytical f, x agree with emperical');

% use the tf to estimate the impedance 
f_input = f_n;
Ts = 1000;
data_est_UPs = iddata(x_n',f_input',1/Ts);
sysUP_s = tfest(data_est_UPs,2,0);
opt = predictOptions('InitialCondition','z');
[xp,~,~] = predict(sysUP_s,data_est_UPs,0,opt);
x_pred = xp.OutputData;
[NUM_UPs,DEN_UPs] = tfdata(sysUP_s);
K_est_up_s  = DEN_UPs{1}(3)/NUM_UPs{1}(3);
B_est_up_s  = DEN_UPs{1}(2)/NUM_UPs{1}(3);
M_est_up_s  = DEN_UPs{1}(1)/NUM_UPs{1}(3);
FIT_up_s    = sysUP_s.Report.Fit.FitPercent;
% estimation agrees with the setting 

%% A time-variant impedance 
clc; 

a = 0.5;
t = 0:1/freq:2;
toff = 0.20;
k = k_n*(1-a*exp(-(t-toff).^2/(2*0.05^2)));
figure(),
hold on;
plot(t,k), grid on
yline(mean(k));
ylim([0 max(k)]);
xlabel('Time [s]')
ylabel('Stiffness [N/m]')
title('examplary stiffness change'); 



i = 1;
for toff = [0.03 0.06 0.1 0.2]
    %Differential Equation
    f = @(t,x) [x(2);-(b_n/m_n)*x(2)-((k_n*(1-a*exp(-(t-toff)^2/(2*0.05^2))))/m_n)*x(1)+F/m_n];

    [tt,xx] = ode45(f,0:1/freq:2,[0,0]);

    x_kvar(:,i) = xx(:,1);
    k = k_n*(1-a*exp(-(t-toff).^2/(2*0.05^2)));
    fce(:,i) = (xx(:,1)' - xx(1,1)).*k + xx(:,2)'.*b_n + [diff(xx(:,2)')/(1e-3), 0]*m_n;

    i = i + 1;
end

x_pred = [];
% estimating each of the system using TI system 
for i = 1:4
%     f = fce(:,i);    
    x = x_kvar(:,i);
    f = ones(size(x))*F;    
    f(1) = 0; 
    f(end) = f(end-1);
    Ts = 1000;
    data_est_UPs = iddata(x,f,1/Ts);

    sysUP_s = tfest(data_est_UPs,2,0);
    opt = predictOptions('InitialCondition','z');
    [xp,~,~] = predict(sysUP_s,data_est_UPs,0,opt);
    x_pred(:,i) = xp.OutputData';
    [NUM_UPs,DEN_UPs] = tfdata(sysUP_s);
    K_est_up_s(i) = DEN_UPs{1}(3)/NUM_UPs{1}(3);
    B_est_up_s(i) = DEN_UPs{1}(2)/NUM_UPs{1}(3);
    M_est_up_s(i) = DEN_UPs{1}(1)/NUM_UPs{1}(3);
    FIT_up_s(i) = sysUP_s.Report.Fit.FitPercent;

%     figure(); 
%     subplot(2,1,1); 
%     hold on;
%     plot(t, f);
%     plot(t, fce(:,i));
%     subplot(2,1,2); 
%     hold on; 
%     plot(t, x);
%     plot(t, x_pred(:,i));
%     legend('actual pos', 'predicted pos');
end


figure(),
subplot(4,1,1);
hold on;
plot([0.03 0.03],[0 50],'--b');
plot(t,x_n*1000);
plot(tt,x_kvar(:,1)*1000);
plot(tt,x_pred(:,1)*1000, ':k');
legend('K peak time', 'TI model', 'TV model', 'TI fitting of TV');
ylabel('Disp. [mm]');
grid on
% title(['K' num2str(K_est_up_s(1)), ' fit' num2str(FIT_up_s(1))  '%']);
title('K peak at 0.03s');


subplot(4,1,2);
hold on;
plot([0.06 0.06],[0 50],'--b');
plot(t,x_n*1000);
plot(tt,x_kvar(:,2)*1000);
plot(tt,x_pred(:,2)*1000, ':k');
grid on
% legend('K peak time', 'TI model', 'TV model', 'TI fitting of TV');
ylabel('Disp. [mm]');
% title(['K' num2str(K_est_up_s(2)), ' fit' num2str(FIT_up_s(2))  '%']);
title('K peak at 0.06s')

subplot(4,1,3)
hold on;
plot([0.1 0.1],[0 50],'--b');
plot(t,x_n*1000);
plot(tt,x_kvar(:,3)*1000);
plot(tt,x_pred(:,3)*1000, ':k');
grid on
% legend('K peak time', 'TI model', 'TV model', 'TI fitting of TV');
ylabel('Disp. [mm]');
% title(['K' num2str(K_est_up_s(3)), ' fit' num2str(FIT_up_s(3))  '%']);
title('K peak at 0.10s')

subplot(4,1,4)
hold on;
plot([0.2 0.2],[0 50],'--b');
plot(t,x_n*1000);
plot(tt,x_kvar(:,4)*1000);
plot(tt,x_pred(:,4)*1000, ':k');
grid on
title(['K' num2str(K_est_up_s(4)), ' fit' num2str(FIT_up_s(4))  '%']);
xlabel('Time [s]');
ylabel('Disp. [mm]');
% legend('K peak time', 'TI model', 'TV model', 'TI fitting of TV');
title('K peak at 0.20s')
% sgtitle('displacement variates with reflex delay');

%% Time-variant impedance on multiple values

clc; 

K_pred_all = zeros(5, 4);
B_pred_all = zeros(5, 4);
M_pred_all = zeros(5, 4);
FIT_pred_all=zeros(5, 4);

k_list = [100 300 500 700 900];
toff_list =  [0.03 0.06 0.1 0.2];

for k_i = 1:5
    k_n = k_list(k_i);
    a = 0.5;
    t = 0:1/freq:2;
    toff = 0.2;
    k = k_n*(1-a*exp(-(t-toff).^2/(2*0.05^2)));


    i = 1;
    for toff = [0.03 0.06 0.1 0.2]
        %Differential Equation
        f = @(t,x) [x(2);-(b_n/m_n)*x(2)-((k_n*(1-a*exp(-(t-toff)^2/(2*0.05^2))))/m_n)*x(1)+F/m_n];

        [tt,xx] = ode45(f,0:1/freq:2,[0,0]);

        x_kvar(:,i) = xx(:,1);
        k = k_n*(1-a*exp(-(t-toff).^2/(2*0.05^2)));
        fce(:,i) = (xx(:,1)' - xx(1,1)).*k + xx(:,2)'.*b_n + [diff(xx(:,2)')/(1e-3), 0]*m_n;

        i = i + 1;
    end

    x_pred = [];
    % estimating each of the system using TI system
    for i = 1:4
        %     f = fce(:,i);
        x = x_kvar(:,i);
        f = ones(size(x))*F;
        f(1) = 0;
        f(end) = f(end-1);
        Ts = 1000;
        data_est_UPs = iddata(x,f,1/Ts);

        sysUP_s = tfest(data_est_UPs,2,0);
        opt = predictOptions('InitialCondition','z');
        [xp,~,~] = predict(sysUP_s,data_est_UPs,0,opt);
        x_pred(:,i) = xp.OutputData';
        [NUM_UPs,DEN_UPs] = tfdata(sysUP_s);
        K_est_up_s(i) = DEN_UPs{1}(3)/NUM_UPs{1}(3);
        B_est_up_s(i) = DEN_UPs{1}(2)/NUM_UPs{1}(3);
        M_est_up_s(i) = DEN_UPs{1}(1)/NUM_UPs{1}(3);
        FIT_up_s(i) = sysUP_s.Report.Fit.FitPercent;

        %     figure();
        %     subplot(2,1,1);
        %     hold on;
        %     plot(t, f);
        %     plot(t, fce(:,i));
        %     subplot(2,1,2);
        %     hold on;
        %     plot(t, x);
        %     plot(t, x_pred(:,i));
        %     legend('actual pos', 'predicted pos');
    end
    K_pred_all(k_i,:) = K_est_up_s;
    B_pred_all(k_i,:) = B_est_up_s;
    M_pred_all(k_i,:) = M_est_up_s;
    FIT_pred_all(k_i,:) = FIT_up_s;
end

K_pred_all_arr = K_pred_all(:);
FIT_pred_all_arr = FIT_pred_all(:);
% K_pred_all_arr(FIT_pred_all_arr<80) = nan;
K_pred_all_val = reshape(K_pred_all_arr,5,4);
plot(toff_list, K_pred_all_val, 'o-', 'linewidth', 2);
title('neglegable difference in model K and estimated K using TI model'); 
ylabel('K (N/m)');
xlabel('reflex peak time (s)');
legend('K = 100N/m', 'K = 300N/m', 'K = 500N/m', 'K = 700N/m', 'K = 900N/m');