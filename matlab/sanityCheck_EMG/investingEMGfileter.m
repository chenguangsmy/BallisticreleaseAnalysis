% investigating EMG filters

% let's say I'm collecting EMG already


% by reading it, execute: 
% sstmpemg = SessionScanEMG(4273) 
% and stop at the middle 


X = emg(4,:);


Fs = 2000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(X);             % Length of signal
t = (0:L-1)*T;        % Time vector

plot(1000*t(1:5000),X(1:5000))
title('Signal of EMG')
xlabel('t (milliseconds)')
ylabel('X(t)')


% do the FFT

% Y = fft(X);
Y = fft(X_filter);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% looks like the data have strong 

%%%%%%%%%%%%%%%%%%%%%%%

% build a notch filter
% wo = 60/(Fs/2);  
wo = 400/(Fs/2);  
bw = wo/35;
[num,dem] = iirnotch(wo,bw);

% fvtool(num,dem)

X_orig = X;
X_filter = filter(num, dem, X_orig); 
figure(); hold on;
plot(t, X_orig); 
plot(t, X_filter); 
legend('origin', 'filtered');

%%%%%%%%%%%%%%%%%%%%%%%%%
% notch filter with a series of band width 
wo_list = [60 200 400 500 600 800];
X_orig0 = X_orig;
for woi = 1:length(wo_list)
    wo = wo_list(woi)/(Fs/2);
    bw = wo/35;
    [num,dem] = iirnotch(wo,bw);
    X_orig1 = X_orig0;
    X_filter = filter(num, dem, X_orig1);
    X_orig0 = X_filter;
end
figure(); hold on;
plot(t, X_orig);
plot(t, X_filter);
legend('origin', 'filtered');
