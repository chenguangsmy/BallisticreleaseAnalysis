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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InvestigEMGfilter using the both MVF data or another EMG data with
% experiment 4310
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% emgtmp = SessionScanEMG(4315);
% sstmp = SessionScan(4315); 
emgtmp = SessionScanEMG(4310);
sstmp = SessionScan(4310); 
data.emg_t   = emgtmp.data.t;
data.emg_raw = emgtmp.data.emgraw;
data.emg_rec = emgtmp.data.emg;         % rectified
data.emg_evl = emgtmp.data.emgevl;
data.fce_t   = sstmp.data.t;
data.fce     = sstmp.data.f;  
% to see if the fce and emg are aligned 
ifplot = 1; 
if(ifplot)
    figure();

    axh(1) = subplot(5,1,1);
    plot(data.fce_t, data.fce); 

    axh(2) = subplot(5,1,2); hold on;
    plot(data.emg_t, data.emg_raw(1,:));
    plot(data.emg_t, data.emg_raw(2,:));

    axh(3) = subplot(5,1,3); hold on;
    plot(data.emg_t, data.emg_raw(3,:));
    plot(data.emg_t, data.emg_raw(4,:));

    axh(4) = subplot(5,1,4); hold on;
    plot(data.emg_t, data.emg_raw(5,:));
    plot(data.emg_t, data.emg_raw(6,:));

    axh(5) = subplot(5,1,5); hold on;
    plot(data.emg_t, data.emg_raw(7,:));
    plot(data.emg_t, data.emg_raw(8,:));

    linkaxes(axh, 'x'); 
end
save('data/emgFilter_mvf.mat', 'data');
% save('data/emgFilter_tsk.mat', 'data');

%% load the data and try different combinations of processing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('data/emgFilter_mvf.mat', 'data'); % mvf block
load('data/emgFilter_tsk.mat', 'data'); % task block
emg_t   = data.emg_t; 
emg_raw = data.emg_raw;
fce_t   = data.fce_t; 
fce     = data.fce; 

%% try the different of filter here 
EMGCutoff = 30; 
Fs = unique(round(1./diff(emg_t)));
clear emg_rct emg_filtered
emg_rct = abs(emg_raw - mean(emg_raw));
EMGSF = Fs; 
EMGOrder = 6;
[EMGb,EMGa] = butter(EMGOrder,EMGCutoff/(EMGSF/2));
emg_filtered = filtfilt(EMGb,EMGa,emg_rct')';

%% plot here 
ifplot = 1; 
clear axh 
if (ifplot)
    figure(); 
    axh(1) = subplot(5,1,1);
    plot(fce_t, fce(1,:)); 
    xlabel('t(s)'); ylabel('Force (N)');

    axh(2) = subplot(5,1,2); hold on;
    plot(emg_t, emg_rct(1,:)); 
    plot(emg_t, emg_filtered(1,:), 'linewidth', 2);
    plot(emg_t, -emg_rct(2,:));
    plot(emg_t, -emg_filtered(2,:), 'linewidth', 2);
    xlabel('t(s)'); ylabel('EMG (portion)');
    legend('ch1 rect', 'ch1 filtered', 'ch2 rect', 'ch2 filtered');

    axh(3) = subplot(5,1,3); hold on;
    plot(emg_t, emg_rct(3,:)); 
    plot(emg_t, emg_filtered(3,:), 'linewidth', 2);
    plot(emg_t, -emg_rct(4,:));
    plot(emg_t, -emg_filtered(4,:), 'linewidth', 2);
    xlabel('t(s)'); ylabel('EMG (portion)');
    legend('ch3 rect', 'ch3 filtered', 'ch4 rect', 'ch4 filtered');

    axh(4) = subplot(5,1,4); hold on;
    plot(emg_t, emg_rct(5,:)); 
    plot(emg_t, emg_filtered(5,:), 'linewidth', 2);
    plot(emg_t, -emg_rct(6,:));
    plot(emg_t, -emg_filtered(6,:), 'linewidth', 2);
    xlabel('t(s)'); ylabel('EMG (portion)');
    legend('ch5 rect', 'ch5 filtered', 'ch6 rect', 'ch6 filtered');

    axh(5) = subplot(5,1,5); hold on;
    plot(emg_t, emg_rct(7,:)); 
    plot(emg_t, emg_filtered(7,:), 'linewidth', 2);
    plot(emg_t, -emg_rct(8,:));
    plot(emg_t, -emg_filtered(8,:), 'linewidth', 2);
    xlabel('t(s)'); ylabel('EMG (portion)');
    legend('ch7 rect', 'ch7 filtered', 'ch8 rect', 'ch8 filtered');

    linkaxes(axh, 'x');
end
sgtitle(['filter with freq ' num2str(EMGCutoff) 'Hz']); 