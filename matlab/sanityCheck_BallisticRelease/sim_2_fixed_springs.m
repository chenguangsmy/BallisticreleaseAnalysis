%% Two-spring and mass system simulation
%    analogy of arm-robot interaction
%
%   wall |-\/\/\/\/- O -/\/\/\/\| wall
%     x: 0 m                    1 m 

close all; clear; clc;

% spring 1
k1 = 10;    % N/m
l1 = 0.6;   % m

% spring 2
k2 = 20;    % N/m
l2 = 0.6;   % m

% initial position of the joint
x0 = 0.7;
xmax = 1;

% simulate dynamics
dt = 0.0001;
m  = 0.2;   % kg; mass of joint

t = 0;
x = x0;
v = 0;
X = [];
F = [];
T = [];
nsteps = 1/dt;
for i = 1:nsteps
    % calc next position using current speed
    dx = v * dt;
    x_next = x + dx;
    
    % calc current forces
    d1 = l1 - x;            % dx of spring 1
    f1 = k1*d1;
    
    d2 = xmax - x - l2;     % dx of spring 2
    f2 = k2*d2;
    
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
    T = [T i*dt];
end

estK = diff(F) ./ diff(X);

figure;
subplot(3, 1, 1); plot(T, F); ylabel('F (N)');
subplot(3, 1, 2); plot(T, X); ylabel('X (m)');
subplot(3, 1, 3); plot(T(2:end), estK); ylabel('K (N/m)');
