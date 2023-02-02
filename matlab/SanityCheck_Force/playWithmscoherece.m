clear; 
Fs = 2000;
band_interest = [5:0.2:20];
t = 0:(1/Fs):1.2;

x_mat = [sin(2*pi*(11*t) - pi/4) ;
        sin(2*pi*(31*t));
        sin(2*pi*(79*t))];

x = sum(x_mat) + 0.5*randn(size(t));

y = sum(sin(2*pi*([11,101,113]'*t))) + 0.5*randn(size(t));

t = repmat(t, 1, 1);
x = repmat(x, 1, 1);
y = repmat(y, 1, 1);
% [Cxy, F] = mscohere(x,y,hamming(1000),500,band_interest,Fs)
[Cxy, F] = mscohere(x,y,hamming(1000),500,band_interest,Fs)

figure();
subplot(3,1,1); 
plot(t,x);

subplot(3,1,2);
plot(t,y);

subplot(3,1,3);
plot(F, Cxy); grid on;

figure();
Fs = 500;
[Pxy,F] = cpsd(x,y,hamming(100),90,band_interest,Fs);

Pxy(Cxy < 0.99) = 0;

plot(F,angle(Pxy)/pi)
title('Cross Spectrum Phase')
xlabel('Frequency (Hz)')
ylabel('Lag (\times\pi rad)')
grid

%% their example
clear; clc
rng default

Fs = 1000;
t = 0:1/Fs:1-1/Fs;

x = cos(2*pi*100*t) + sin(2*pi*200*t) + 0.5*randn(size(t));
y = 0.5*cos(2*pi*100*t - pi/4) + 0.35*sin(2*pi*200*t - pi/2) + 0.5*randn(size(t));

[Cxy,F] = mscohere(x,y,hamming(100),80,100,Fs);

plot(F,Cxy)
title('Magnitude-Squared Coherence')
xlabel('Frequency (Hz)')
grid

%%
[Pxy,F] = cpsd(x,y,hamming(100),80,100,Fs);

Pxy(Cxy < 0.2) = 0;

plot(F,angle(Pxy)/pi)
title('Cross Spectrum Phase')
xlabel('Frequency (Hz)')
ylabel('Lag (\times\pi rad)')
grid