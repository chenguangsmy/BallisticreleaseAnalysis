% Create a math model that describe ballistic release: 
% System: 
% |                                        |
% |                                        |
% |                                        |
% |---Spring(r)------ Mass -----Spring(s)--|
% |                                        |
% |                                        |
% |                                        |

% position:     x = 0.1m, the position of the Mass
%               xs0 = 0.15m, the origin length of spring(s)
%               xr0 = 0.15m, the origin length of spring(r)
%               xw1 = 0.20m, the right wall position (left wall 0).
% Force:        F, the force exerted in the non-compressible connection. 
% acceleration: a = x``;
% stiffness:    K_r1 = 320N/m
%               K_r2 = 160N/m
%               K_s  = 320N/m
% Mass:         m = 2kg
% Perturbation: at time 0, the K_r1 change to K_r2

% Question:     Describe the force (F) change with time,
%                        the position (x) change with time,
%               how dF/dx correlate with spring stiffness 
% 
clf; clc; clear;
x0 = 0.1;
xs0 = 0.15;
xr0 = 0.15;
xw1 = 0.20;
Kr1 = 320;
Kr2 = 160;
Ks  = 320;
m = 0.1;

F1 = (Kr1*Ks)/(Kr1+Ks)*(xs0+xr0-xw1); % F(t=0-)
F2 = (Kr2*Ks)/(Kr2+Ks)*(xs0+xr0-xw1); % F(t=+inf)

dt = 1e-4;
t = 0;
x = x0;
v = 0;
X = [];
F = [];
Fr= [];
T = [];

i = 1;
t_stt = 0;
nsteps = 1/dt;
for i = 1:nsteps
%     time(i+1) = t+t_step;
%     F_r(i+1) = Kr2*(xr0-x(i)); 
%     F_s(i+1) = Ks*(xw1-xs0-(x(i))); %?
%     F(i+1)   = F_r(i)-F_s(i);
%     a(i+1)   = F(i)/m;
%     v(i+1)   = v(i)+a(i);
%     x(i+1)   = x(i)+v(i)*t_step;
%     i = i+1;
    % calc next position using current speed
    dx = v * dt;
    x_next = x + dx;
    x = min(max(x, 0), xw1);
    
    % calc current forces
    d1 = xr0 - x;            % dx of spring 1
    f1 = Kr2*d1;
    
    d2 = xw1 - xs0 - x;     % dx of spring 2
    f2 = Ks*d2;
    
    f = f1 + f2;
    acc = f / m;
    
    % calc next speed
    dv = acc * dt;
    v_next = v + dv;
    
    % for next time stamp
    x = x_next;
    v = v_next;
    
    X = [X x];
    F = [F f];
    Fr= [Fr f1];
    T = [T i*dt];
    
end
estK = diff(Fr) ./ diff(X);

figure;
subplot(3, 1, 1); plot(T, F); ylabel('F (N)'); hold on;
subplot(3, 1, 2); plot(T, X); ylabel('X (m)'); hold on;
subplot(3, 1, 3); plot(T(2:end), estK); ylabel('K (N/m)'); hold on;
