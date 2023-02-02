% Assume a pulse was exerted on a second-order system 
% F(t) = Kx + Bx' + Mx''
% close all; 
% clear; 
clc

% build a transfer function
% wn = sqrt(k/m)
% z = B/2*sqrt(k*m)
K_parm = 250;
M_parm = 2; 
B_parm = 20;
wn = sqrt(K_parm/M_parm);
% z = B_parm/(2*sqrt(K_parm*M_parm));
% z = B_parm/(2*sqrt(K_parm*M_parm));
% [num,den] = ord2(wn,z)
% num = num; % unknown???
% den = den;
num = 1;
den = [M_parm B_parm K_parm];


G=tf(num,den)
% create a gaussian force pulse
Fs = 2000;
t = linspace(0, 1, Fs);
Ugus = gaussmf(t,[0.015 0.05])*15;
Ugus = yF; % get from another script, aimed to do a position control
y = lsim(G, Ugus, t);

% calculate the force from each component
yd = [diff(y)*Fs; 0];
ydd = [diff(yd)*Fs; 0];
F_mass = M_parm * ydd; 
F_damp =-B_parm * yd;
F_stff = K_parm * y;

plt_idx = 1:20:length(t);
figure
axh(1) = subplot(2,1,1);
hold on
grid 
plot(t(plt_idx), Ugus(plt_idx), 'k--')
plot(t(plt_idx), F_mass(plt_idx), 'k--o');
plot(t(plt_idx), F_damp(plt_idx), '.k');
plot(t(plt_idx), F_stff(plt_idx), 'k+', 'linewidth', 1);
% plot(t, +(F_mass'-F_damp'+F_stff'), '.b')
legend('F_{pulse}', 'F_{M}', 'F_{B}', 'F_{K}');
% legend('F_{pulse}', 'F_{M}', 'F_{B}', 'F_{K}', 'F_{all}');
title('Forces');
ylabel('F (N)');


axh(2) = subplot(2,1,2);
plot(t, y, '-k')
legend('x'); 
title('displacement');
grid
linkaxes(axh, 'x');
xlim([0 0.5]);
xlabel('time (s)');
ylabel('x (m)');
sgtitle('force and displacement');


%% two comparations of "Human high stiffness" and "monkey low stiffness"

K_parm_alt1 = 900;
M_parm_alt1 = 2; 
B_parm_alt1 = 20;

K_parm_alt2 = 300;
M_parm_alt2 = 0.2; 
B_parm_alt2 = 20;

num = 1;
den1 = [M_parm_alt1 B_parm_alt1 K_parm_alt1];
den2 = [M_parm_alt2 B_parm_alt2 K_parm_alt2];
G1 = tf(num,den1);
G2 = tf(num,den2);

% create a gaussian force pulse
Fs = 2000;
t = linspace(0, 1, Fs);
Ugus = gaussmf(t,[0.015 0.05])*15;
y1 = lsim(G1, Ugus, t);
y2 = lsim(G2, Ugus, t);

% calculate the force from each component
y1d = [diff(y1)*Fs; 0];
y1dd = [diff(y1d)*Fs; 0];
F_mass1 = M_parm_alt1 * y1dd; 
F_damp1 =-B_parm_alt1 * y1d;
F_stff1 = K_parm_alt1 * y1;

y2d = [diff(y2)*Fs; 0];
y2dd = [diff(y2d)*Fs; 0];
F_mass2 = M_parm_alt2 * y2dd; 
F_damp2 =-B_parm_alt2 * y2d;
F_stff2 = K_parm_alt2 * y2;

plt_idx = 1:30:length(t);
figure
axh(1) = subplot(3,1,1);
hold on
grid 
% plot(t(plt_idx), Ugus(plt_idx), 'k--')
% plot(t(plt_idx), F_mass(plt_idx), 'k--o');
% plot(t(plt_idx), F_mass1(plt_idx), 'g--o');
% plot(t(plt_idx), F_mass2(plt_idx), 'b--o');
plot(t(plt_idx), Ugus(plt_idx), 'k')
plot(t(plt_idx), F_mass(plt_idx), 'k-o');
plot(t(plt_idx), F_mass1(plt_idx), 'g-o');
plot(t(plt_idx), F_mass2(plt_idx), 'b-o');

legend('F_{pulse}', 'F_{M}', 'F_{M}alt1', 'F_{M}alt2');
title('Force_M');
ylabel('F (N)');

axh(2) = subplot(3,1,2);
hold on; 
grid
plot(t(plt_idx), F_stff(plt_idx), 'k--+', 'linewidth', 1);
plot(t(plt_idx), F_stff1(plt_idx), 'g--+', 'linewidth', 1);
plot(t(plt_idx), F_stff2(plt_idx), 'b--+', 'linewidth', 1);

% plot(t, +(F_mass'-F_damp'+F_stff'), '.b')
legend('F_{K}', 'F_{K}alt1', 'F_{K}alt2');
% legend('F_{pulse}', 'F_{M}', 'F_{B}', 'F_{K}', 'F_{all}');
title('Force_K');
ylabel('F (N)');


axh(3) = subplot(3,1,3);
hold on; 
plot(t, y, 'k-');
plot(t, y1, 'g-');
plot(t, y2, 'b-');
legend('x', 'x\_alt1', 'x\_alt2'); 
title('displacement');
grid
linkaxes(axh, 'x');
xlim([0 0.5]);
xlabel('time (s)');
ylabel('x (m)');
sgtitle('force and displacement');

