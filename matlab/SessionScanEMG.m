classdef SessionScanEMG
    %SESSIONSCANEMG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ss_num
        anin_idx = 1:8
        anin_idx_offset = 0
        dat             % intermediate
        data            % full
        unit
        time
        brtime          % blackrock time
        freq
        fname_raw
        
        CHANNELS_NUM = 128;
        AINPUTS_NUM  = 16;
        DATA_PER_SHORT = 10;
        
        emg_means = [...
                    51.0777
                    23.6027
                    41.4865
                    80.8630
                    16.5060
                    46.2156
                    36.5228
                    23.6265];
        emg_stds = 1.0e+03 * [ ...
                    1.7848
                    3.1263
                    0.3451
                    1.2234
                    4.0816
                    1.3261
                    0.8344
                    1.2729];


    end
    
    methods
        function obj = SessionScanEMG(ss_num)
            %SESSIONSCANEMG Construct an instance of this class
            %   Detailed explanation goes here
            %SESSIONSCANEMG scan from the intermediate mat file (which 
            %   contains the EMG continuous sitnal), or the raw file and 
            %   get a variable.
            %   Current setting is from the `intermediate` file
            %   header read from %DT sequence% described above.

            % 1. read data 
            %    Read either from `IntermediateFile` 
            %                  or `ASeperateFile`
            obj.ss_num = ss_num;
            fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate';
            fname = sprintf('KingKong.%05d.mat', ss_num);
            filename_intermediate = [fdir '/' fname];
            % ... read the intermediate file
            load(filename_intermediate, 'Data'); 
            DataInt = Data;
            Dat = DataInt.QL.Data.RAW_CTSDATA.data;
            N = size(Dat,2);
            dat_all = reshape(Dat, obj.CHANNELS_NUM+obj.AINPUTS_NUM, ...
                N*obj.DATA_PER_SHORT);
            obj.dat = dat_all(obj.anin_idx+obj.anin_idx_offset,:);
            time_msg = DataInt.QL.Data.RAW_CTSDATA.source_timestamp; 
            % get brtime
            obj.brtime = obj.getBRtimefromIntermediate(time_msg);
            % get RTMA time (task time)
            obj.time   = obj.getRTMAtimeAligned(); % to write
            fnameraw = sprintf('KingKongEMG.%05d.csv', ss_num);
            fdiremg = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
            obj.fname_raw = [fdiremg '/' fnameraw];
            try
                obj = obj.readRawData();
                rawdata_flag = 1;
            catch 
                disp('SessionScanEMG: no raw EMG readed');
                obj.data.emg = zeros(8,0);
                rawdata_flag = 0;
            end

            % if no raw data, fill the raw data with the intermediate data
            if (rawdata_flag == 0)
                obj.data.emg = double(obj.dat);
                obj.data.t = obj.brtime';
                obj.freq = round(mean(1./diff(obj.data.t)));
            end

            % magnify the channels according to the magnification. 
            amps = readConfigChannelAmp(obj.ss_num);
            obj = obj.convertTomV(amps);
            % preprocess data (filter, normalize and take the envolope)-
            obj = obj.preprocessRawData(1);

            % rotate the channels for the config file 
            chmap = readConfigChannelMap(obj.ss_num);
            obj.data.emg = obj.data.emg(chmap, :);
        end
        
        function time_highsampled = getBRtimefromIntermediate(obj, time_msg)
            % time_highsampled = getBRtimefromIntermediate(time_msg)
            % Out put the high_sampled time from the data read from
            % intermediate file. 
            % as the intermediate 
            ifdisp = 0;
            if (ifdisp)
                disp([num2str(sum(isnan(time_msg))) 'timepoints were nan']);
            end
            
            timeidx_msg = (1:length(time_msg)) * obj.DATA_PER_SHORT;
            timeidx_raw = 1:timeidx_msg(end);
            time_highsampled = interp1(timeidx_msg, time_msg, timeidx_raw, 'linear', 'extrap');
        end

        function time = getRTMAtimeAligned(obj)
            % align the task time with the current data
            % so that I can directly plot the current data time with data
            disp('NOT FINISHED GETRTMATIMEALIGNED YET! ERROR :)))');
            time = [];
        end
        
        function obj = convertTomV(obj, amp)
            % convert the EMG data to mV measured, un-gain the OCTOPUS and
            % the BLACKROCK
            op_amp = amp(:,1);
            bk_amp = amp(:,2);
            bk_range = 5000 ./ bk_amp;
            sampled_range = 32767;
            % 1. calculate the voltage to BLACKROCK input 
            emgBkmV = obj.data.emg/sampled_range.*bk_range;
            % 2. calculate the voltage to OCTOPUS input 
            emgOcmV = emgBkmV./op_amp;

            obj.data.emg = emgOcmV;
        end

        function obj = readRawData(obj)
            %data = readRawData(obj)
            %   Read raw data from the *.csv file
            data_tmp = readtable(obj.fname_raw);
            t_tmp = [data_tmp{:,1}]';
            idx_tmp = 1:length(t_tmp); 
            [a, ia, ic] = unique(t_tmp);
            time_first_idx = find([0; diff(ic) == 1]);
            time_first = t_tmp(time_first_idx); 
            ifplot = 0;
            if (ifplot)
                clf; hold on;
                plot(idx_tmp, t_tmp, 'b.');
                plot(time_first_idx, time_first, 'ro');
            end
            t = interp1(time_first_idx, time_first, idx_tmp, 'linear', 'extrap')'; 
            if (ifplot)
                clf; hold on;
                plot(idx_tmp, t_tmp, 'bo');
                plot(idx_tmp, t, 'r.');
                legend('origin', 'introp');
            end
            obj.data.t = t;
            obj.freq = [unique(data_tmp{:,2})];
            obj.data.emg = [data_tmp{:,3:10}'];
        end
        
        function obj = preprocessRawData(obj, ifplot)

            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            emg = obj.data.emg;
%             emg_processed = zeros(size(emg));
%             emg_processed0= zeros(size(emg));
            t = obj.data.t;

            % remove line noise from the raw data 
            [t_filter, emg_filter] = removeLineNoise(t, emg);
            
            [t_filter, emg_filter] = removeMotionNoise(t_filter, emg_filter);
            %
            t = t_filter; 
            emg = emg_filter;
            emg_processed = zeros(size(emg));
            emg_processed0= zeros(size(emg));
                
            for chi = 1:8 % iterate through channels
                % 1. do the 100Hz low pass filter
                fs = obj.freq; % data frequency
                lpf1 = 100; % were60 
%                 emg1 = bandstop(emg(ch_i,:),[118.5 121.5],fs);
                emg_stp1 = highpass(emg(chi,:), lpf1, fs);
%                 emg_stp1 = emg(chi,:);
                % 2. mean-centered and scaled by standard deviation...
%                 emg_stp2 = (emg_stp1 - obj.emg_means(chi)) / obj.emg_stds(chi);
                emg_stp2 = detrend(emg_stp1);
                % 3. squared, and do low-pass filter of 30Hz
                lpf2 = 30; %30
%                 emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
                emg_stp3 = emg_stp2.^2;
                % 3. squre root transform and *2
                emg_stp4 = sqrt(emg_stp3)*2;                                                   % bad name, need chagne
%                 emg_stp4 =  abs(emg(chi,:)/32767*5000/2000); % convert to mV % already mV since 2022-06-09
                [emg_stp5, ~] = envelope(real(emg_stp4), 50, 'peak'); % check what it used to be?
%                 [emg_stp5, ~] = envelope(real(emg_stp4), 150, 'peak'); % check what it used to be?
%                 ifplot = 1;
                if (ifplot)
                    clf; hold on;
                    plot(t, emg(chi,:), 'Color', [0.7 0.7 0.7], 'Marker', '.')
                    plot(t, emg_stp1, 'r', 'Marker', '.');
                    plot(t, emg_stp2, 'b', 'Marker', '.');
                    plot(t, emg_stp3, 'g', 'Marker', '.');
                    plot(t, emg_stp4, 'm', 'Marker', '.');
                    plot(t, emg_stp5, 'c', 'Marker', '.');
                    ylim([-0, 0.2]);
                    title(['chennel ' num2str(chi)]);
                    xlabel('time (s)' );
                    ylabel('magnitude (mV)');
                    legend('raw', 'highpass', 'z-score', 'lowpass', 'sqrt transform', 'envolope');
                end
                emg_processed0(chi,:)=emg_stp4;
                emg_processed(chi,:) = emg_stp5;%emg_stp4;
                
            end
            obj.data.emg = emg_processed;
            obj.data.t = t;
%             ifplot = 1;
            if (ifplot)
                clf;
                % plot

                axh(2) = subplot(5,1,2); hold on;
%                 plot(t,emg_processed0(1:2,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(1:2,:)); %title('EMG12'); %ylabel('N');
                legend('1', '2');
                
                axh(3) = subplot(5,1,3); hold on;
%                 plot(t,emg_processed0(3:4,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(3:4,:)); %title('EMG34'); %ylabel('N');
                legend('3', '4');
                
                % hold on; plot(t,emg_processed(3,:));
                % plot(t,emg_processed0(3,:)); %title('EMG34'); %ylabel('N');
                
                axh(4) = subplot(5,1,4); hold on;
%                 plot(t,emg_processed0(5:6,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(5:6,:)); %title('EMG56'); %ylabel('N');
                legend('5', '6');
                
                axh(5) = subplot(5,1,5); hold on;
%                 plot(t,emg_processed0(7:8,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(7:8,:)); %title('EMG78'); %ylabel('N');
                legend('7', '8');
                
                linkaxes(axh, 'x');
                linkaxes(axh(2:end), 'y');
%                 ylim(axh(5), [0, 10]);
                ylim(axh(5), [-0.2 0.2]);
                
                linkaxes(axh, 'x');
            end
        end
        
        function [t_list] = scanThresholdCrossing(obj, threshold)
            % [t_list, ch_list] = scanThresholdCrossing(threshold)
            % return of a list of time and channel number as threshold 
            % time is a channel-num array 

            % threshold is measured in mV
            t_list = cell(8,1);
            for ch_i = 1:8
                t_list_idx = abs(obj.data.emg(ch_i,:)) > threshold;
                t_list{ch_i} = obj.data.t(t_list_idx);
            end
        end


        function axh = plotRawEMGData(obj)
            % axh = plotRawEMGData(obj)
            % plot the raw data of EMG
            % with BlackRock time, and the raw data
            axh = figure(); 
            axesh = zeros(1,8);
            for axi = 1:8
                axesh(axi) = subplot(8,1,axi);
                plot(obj.brtime, obj.dat(axi,:));
                ylabel(['ch' num2str(axi)] );
                switch axi
                    case 1
                        title('raw EMG data');
                    case 8
                        xlabel('time (BlackRock Sys)');
                end
            end
            linkaxes(axesh, 'x');
        end

        function fh = plotRawData(obj)
            % fh = plotRawData(obj)
            % plot the raw data in each subplot
            fh = figure(); 
            for i = 1:8
                axh(i) = subplot(8,1,i); grid on
                plot(obj.data.t, obj.data.emg(i,:), 'Marker' , '.');
            end
            linkaxes(axh, 'xy')
            sgtitle('EMG data');
        end
    end
end

function chmap = readConfigChannelMap(ssnum)
% read the channel-map, channel map to the muscles
filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetEMGChannelMap.conf';
fid = fopen(filename);
C = textscan(fid, '%s\n','CommentStyle','#');
fclose(fid);
for li = 1:size(C{1}, 1)
    str = C{1}{li};
    freadtmp = textscan(str,'%d,');
    ss_num = freadtmp{1}(1);
    if ss_num ~= ssnum
        continue;
    end
    chmap_str = freadtmp{1}(2:end);
end
if exist('chmap_str', 'var')    % defined by the config file
    chmap = double(chmap_str)';
    clear chmap_str;
else % default value
    chmap = 1:8;    % The right order
end
end

function amp = readConfigChannelAmp(ssnum)
% find the amplifier gain on each channel
% return value:
% amp = [O_gain for each channel (1~8), B_gain for each channel (1~8)];
filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetEMGAmpData.conf';
fid = fopen(filename);
C = textscan(fid, '%s\n','CommentStyle','#');
fclose(fid);
for li = 1:size(C{1}, 1)
    str = C{1}{li};
    freadtmp = textscan(str,'%d,');
    ss_num = freadtmp{1}(1);
    if ss_num ~= ssnum
        continue;
    end
    amp_str = freadtmp{1}(2:end);
end
if exist('amp_str', 'var')
    amp = double(amp_str)';
    amp = reshape(amp, 8, 2);
    clear amp_str;
else
    amp = [500*ones(8,1), 10*ones(8,1)];    % The right order
end
end

function [t_filter, emg_filter] = removeLineNoise(t, emg, ifplot)
% might need de-trend (cg)
if (~exist('ifplot', 'var'))
    ifplot = 0;
else

end

% remove line noise from the filter
Fs_raw = mean(1./diff(t));

if (Fs_raw == 2000)
    L = (max(t) - min(t))/(1./round(Fs_raw));
else
    L = (max(t) - min(t))/(1./2000);
end
t_resample = linspace(min(t), max(t), L);
Fs = mean(1./diff(t_resample));
t_filter = t_resample;
emg_filter = zeros(8, length(t_filter));

dat = emg(:,:);
% resample with the ideal frequency
emg_resample = interp1(t, dat', t_resample', 'spline');
emg_resample = emg_resample';

%     ifplot = 1;
if (ifplot)
    figure();
end

for ch_i = 1:size(emg,1) % each channel
%     figure('name', ['channel' num2str(ch_i)])
    % remove the line noise
    wo = 60/(Fs/2);
    bw = wo/35;
    [num,dem] = iirnotch(wo,bw);

    X_orig = emg_resample(ch_i,:);
    X_filter = filter(num, dem, X_orig);
    emg_filter(ch_i,:) = X_filter;

    % use fft to see whether there is high line noise
    % original data
    Y1 = fft(X_orig);
    P12 = abs(Y1/L);                  % two-side spectrum
    P11 = P12(1:L/2+1);               % single-side spectrum
    P11(2:end-1) = 2*P11(2:end-1);
    f1 = Fs*(0:(L/2))/L;
    % filtered data
    Y2 = fft(X_filter);
    P22 = abs(Y2/L);                  % two-side spectrum
    P21 = P22(1:L/2+1);               % single-side spectrum
    P21(2:end-1) = 2*P21(2:end-1);
    f2 = Fs*(0:(L/2))/L;

    % plot them out
    if (ifplot)
        clf;
        subplot(2,1,1); hold on;
        plot(t_resample, X_orig);
        plot(t_resample, X_filter);
        legend('origin', 'filtered', 'lowpassed');
        title('raw and filtered signal of X(t)')
        xlabel('t');
        ylabel('signal value (mV)'); % assume its mV

        subplot(2,1,2); hold on;
        % spectrum
        plot(f1,P11);
        plot(f2,P21);
        legend('origin', 'filtered');
        title('Single-Sided Amplitude Spectrum of X(t)');
        xlabel('f (Hz)');
        ylabel('|P1(f)|');
        % return the data
    end
end

end

function [t_filter, emg_filter] = removeMotionNoise(t_filter, emg_filter, ifplot)
% remove motion noise with crazy values
if (~exist('ifplot', 'var'))
    ifplot = 0;
else

end

for ch_i = 1:size(emg_filter)
    emg_filter1(ch_i,:) = filloutliers(emg_filter(ch_i,:), 'clip', 'movmedian', [200 200]);
    t_outlair_idx = abs(emg_filter(ch_i,:) - emg_filter1(ch_i,:)) > 3*std(emg_filter1(ch_i,:));
    emg_filter2(ch_i,:) = emg_filter(ch_i,:);
    emg_filter2(ch_i,t_outlair_idx) = emg_filter1(ch_i,t_outlair_idx);


    if (ifplot)
        clf; hold on;
        plot(t_filter, emg_filter(ch_i,:));
        plot(t_filter, emg_filter2(ch_i,:));
        legend('before filter', 'after filter');
        title('filter high magnitude outlairs');
    end
end
emg_filter = emg_filter2;
end


