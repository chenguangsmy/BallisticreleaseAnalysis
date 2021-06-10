 clear all
 close all
 clc
 
 syms s
 syms m11 m12 m21 m22
 syms b11 b12 b21 b22
 syms k11 k12 k21 k22
 
 M = [m11 m12; m21 m22];
 B = [b11 b12; b21 b22];
 K = [k11 k12; k21 k22];
 
 Z = M*s^2+B*s+K;
 C = inv(Z);

[N,D] = numden(C);

coeffs(N(1,1),s)
coeffs(D(1,1),s)

%%

clear all
close all
clc

s = tf('s');
M = [1.7099, -0.2566; -0.2566, 2.1775];
B = [5.2510, -1.0215; -1.0215, 39.0782];
K = [105.0196, -20.4292; -20.4292, 781.5645];
    
sfrq = 100;
dt = 1/sfrq;
t = 0:dt:2;

Z = M*s^2+B*s+K;
% Z = M(2,2)*s^2+B(2,2)*s+K(2,2);
C_total = inv(Z);
C_c = C_total(1,1);
H_c = impulse(C_c,t);
C_d = c2d(C_c,dt,'matched');
H_d = impulse(C_d,t);
% H_d2 = impz(C_d.numerator{1},C_d.Denominator{1},length(t));
% figure;plot(t,H_d1,t,H_d2);

figure; plot(t, H_c,'-',t,H_d,'o');

[b,a] = prony(H_d,4,4);
C_hat_d = tf(b,sfrq*a,dt); % Why is this off by a factor of the sampling frequency??
C_hat_c = d2c(C_hat_d,'matched');
H_hat_d = impulse(C_hat_d,t); 
H_hat_c = impulse(C_hat_c,t);

figure; plot(t,H_c,'-',t,H_hat_c,'o','linewidth',2.5);
figure; plot(t,H_d,'-',t,H_hat_d,'o','linewidth',2.5);

syms s m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22
 
M = [m11 m12; m21 m22];
B = [b11 b12; b21 b22];
K = [k11 k12; k21 k22];

Z = M*s^2+B*s+K;
C = inv(Z);

[N,D] = numden(C);
n(:,1,1) = coeffs(N(1,1),s);
d(:,1,1) = coeffs(D(1,1),s);

n(:,1,2) = coeffs(N(1,2),s);
d(:,1,2) = coeffs(D(1,2),s);

n(:,2,1) = coeffs(N(2,1),s);
d(:,2,1) = coeffs(D(2,1),s);

n(:,2,2) = coeffs(N(2,2),s);
d(:,2,2) = coeffs(D(2,2),s);

equns = [n(:,1,1) == C_hat_c.numerator{1}(end-2:end)';...
         d(:,1,1) == C_hat_c.Denominator{1}(:);...
         ];

vars = [m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22];

[m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22] = solve(eqns,vars);


