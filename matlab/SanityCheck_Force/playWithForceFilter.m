clear all; close all; clc;
load('20220921_forcedata.mat', 'datfmt');

% 
dat = datfmt; 
dat.t = dat.t - dat.t(1);

Fs = unique(round(1./diff(dat.t)));

% show the raw 
figure('name', 'raw data');
plot(dat.t, dat.f);
xlabel('time (s)'); 
ylabel('Force (N)'); 

% see the power spectrum 
[P, F] = pspectrum(dat.f, dat.t);

figure('name', 'power spectrum'); hold on;
plot(F, P); 
xlabel('frequency (Hz)');
ylabel('power');

[P, F] = pspectrum(dat.f, Fs);

plot(F, P); 

xlim([1 50]);