classdef SessionScanEMG
    %SESSIONSCANEMG Summary of this class goes here
    %   Detailed explanation goes here
    
    %%%%%%%%%%%%%%%%
    % todo: have a bandstop notch filter
    %%%%%%%%%%%%%%%
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
        fname_rawm      % for .mat file
        fname_fmt
        CHANNELS_NUM = 128;
        AINPUTS_NUM  = 16;
        DATA_PER_SHORT = 10;
        mvf_shrink = true;
        amp_norm = zeros(1,8);
        if_calibSS
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
        function obj = SessionScanEMG(ss_num, ifReadAgain)
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

            if (~exist('ifReadAgain', 'var'))
                ifReadAgain = 1;
            end
            obj.ss_num = ss_num;
            fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate';
            fdiremg = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
            fnameraw = sprintf('KingKongEMG.%05d.csv', ss_num);
            fnamerawm= sprintf('KingKongEMG.%05d.mat', ss_num);
            fnamefmt = sprintf('KingKongEMGFMT.%05d.mat', ss_num); % formatted
            obj.fname_raw = [fdiremg '/' fnameraw];
            obj.fname_rawm = [fdiremg '/' fnamerawm];
            obj.fname_fmt = [fdiremg '/' fnamefmt];
            fname = sprintf('KingKong.%05d.mat', ss_num);
            filename_intermediate = [fdir '/' fname];
            % ... read the intermediate file
            
            load(filename_intermediate, 'Data'); 
            DataInt = Data;
            Dat = DataInt.QL.Data.RAW_CTSDATA.data;
            N = size(Dat,2);
            dat_all = reshape(Dat, obj.CHANNELS_NUM+obj.AINPUTS_NUM, ...
                N*obj.DATA_PER_SHORT);
            % message runs at 100Hz, and get 10 data points per msg
            obj.dat = dat_all(obj.anin_idx+obj.anin_idx_offset,:);
            time_msg = DataInt.QL.Data.RAW_CTSDATA.source_timestamp; 
            % get brtime
            obj.brtime = obj.getBRtimefromIntermediate(time_msg);
            % get RTMA time (task time)
            obj.time   = obj.getRTMAtimeAligned(); % to write
            try
                obj = obj.readRawData();
                rawdata_matflag = 1;
            catch 
                obj = obj.readRawDataM(); 
%                 disp('SessionScanEMG: no raw EMG readed');
                obj.data.emg = zeros(8,0);
                rawdata_matflag = 0;
            end

            % if no raw data, fill the raw data with the intermediate data
            if (rawdata_matflag == 0)
                obj.data.emg = double(obj.dat);
                obj.data.t = obj.brtime';
                obj.freq = round(mean(1./diff(obj.data.t)));
            end

            obj.if_calibSS = findCalibSS(obj.ss_num);

            if (ifReadAgain || ~exist(obj.fname_fmt, 'file') || obj.if_calibSS) % If want to re-perform pre-process
                
                % convert the channel to the right order 
                chmap = readConfigChannelMap(obj.ss_num);
                obj.data.emg = obj.data.emg(chmap, :);
                % magnify the channels according to the magnification.
                amps = readConfigChannelAmp(obj.ss_num);
                obj = obj.convertTomV(amps);
                
                % do the frequency-based processing 
                % preprocess data (filter, and take the envolope)-
                obj = obj.preprocessRawData(1);
                
                % convert chanel to the MVF
                [obj.amp_norm, exist_flag] = readConfigChannelMVF(obj.ss_num);
                obj.mvf_shrink = exist_flag;
                obj = obj.normalizebyMVF(obj.amp_norm); % already

                data = obj.data;
                save(obj.fname_fmt, 'data');
            else    % directly read from existing data
                data = load(obj.fname_fmt);
                obj.data = data.data;
            end
            % save data 

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
            bk_amp = ones(size(amp(:,2))); % Chenguang's test 
%             anl_range = 5000./bk_amp;
            anl_range = 5000 * ones(size(bk_amp));
            digital_range = 32767;
            % 1. calculate the voltage to BLACKROCK input 
            emgBkmV = obj.data.emg/digital_range.*anl_range;
            % 2. calculate the voltage to OCTOPUS input 
            emgOcmV = emgBkmV./op_amp;

            obj.data.emg = emgOcmV;
        end

        function obj = normalizebyMVF(obj, amp)
            % convert the EMG data to % of MVF measured
            % devided by the MVF
            % 1. calculate the voltage to BLACKROCK input
            amp_mat = repmat((1./amp)',1,size(obj.data.emg,2));
            emg_pct_mvf = obj.data.emgevl.*amp_mat;

            obj.data.emgevl = emg_pct_mvf;
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
        
        function obj = readRawDataM(obj)
            %data = readRawData(obj)
            %   Read raw data from the *.mat file
            data_tmp = load(obj.fname_rawm);
            % concatinate each file
            varnames = fieldnames(data_tmp);
            num_trials = length(varnames)/2;
            t_tmp = [];
            d_tmp = [];
            for trial_i = 1:(num_trials)
                t_tmp = [t_tmp; data_tmp.(['time' num2str(trial_i-1)])];
                d_tmp = [d_tmp; data_tmp.(['dat' num2str(trial_i-1)])];
            end

            %%% A sanity check whether the data has some lost 
            sample_freq = 2000;
            t_tmpidx = t_tmp == max(t_tmp);
            t_tmp_range = max(t_tmp) - min(t_tmp) + (sum(t_tmpidx)-1)*1/sample_freq; 
            t_sample_should = t_tmp_range*sample_freq; 
            t_sample_prac = length(t_tmp);
            fprintf('>> EMG sanity check: Sample Theory: %d, Pract: %d, miss: %d points in %.2f(s)', ...
                round(t_sample_should), t_sample_prac, round(t_sample_should-t_sample_prac), t_tmp_range);
%             t_tmp = [data_tmp{:,1}]';
            ifplot = 1;
            if (ifplot)
                clf; 
                hold on;
                plot(t_tmp(1)+((1:length(t_tmp))/2000), t_tmp(1)+((1:length(t_tmp))/2000), 'b.'); 
                plot(t_tmp(1)+((1:length(t_tmp))/2000), t_tmp, 'r.'); 
                xlabel('time start from 1st data (s)');
                ylabel('time when data recorded (s)');
            end
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
%             t = interp1(time_first_idx, time_first, idx_tmp, 'linear', 'extrap')'; 
%             t = (idx_tmp-1)*(1/sample_freq) + t_tmp(1); 
              t_first = min(t_tmp); 
              t_last = max(t_tmp) + (sum(t_tmpidx)-1)*1/sample_freq; 
              t = linspace(t_first,t_last,length(t_tmp));

            if (ifplot)
                clf; hold on;
                plot(idx_tmp, t_tmp, 'bo');
                plot(idx_tmp, t, 'r.');
                legend('origin', 'introp');
            end
            obj.data.t = t;
            obj.freq = mean(1./diff(t));
%             obj.freq = [unique(data_tmp{:,2})];
            obj.data.emg = d_tmp;
        end

        function obj = preprocessRawData(obj, ifplot)

            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            emg_raw = obj.data.emg;
%             emg_processed = zeros(size(emg));
%             emg_processed0= zeros(size(emg));
            t_raw = obj.data.t;
            % do something deal with the time
            % see the power spectrum of the current signal 
            figure(); 
            axh = subplot(1,1,1);
%             axh = obj.powerspectrum(emg_raw, t, axh); % data, time, and axis

            % try to remove the flat shifting data
            emg_raw1 = zeros(size(emg_raw));
            for ch_i = 1:8 
                emg_raw1(ch_i,:) = detrend(emg_raw(ch_i,:));
%                 emg_raw1(ch_i,:) = emg_raw1(ch_i,:) - mean(emg_raw1(ch_i,:));
            end
            % remove line noise from the raw data 

            % look the power spectrum, and spectrugram for each of the
            % channel (see if there is noise throughout the session) 
            for chi = 1:8 % what is the dimension of p_tmp and f_tmp???
                [p_tmp(chi,:), f_tmp(chi,:)] = pspectrum(emg_raw1(chi,:), 2000);
            end
            if (ifplot) % the raw data 
                figure('name', ['raw spectrum of Session' num2str(obj.ss_num)])
                hold on; 
                for chi = 1:8
                    plot(f_tmp(chi,:), log(p_tmp(chi,:)));
                end
                legend('ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8');
            end

            [t_filter1, emg_filter1] = obj.removeLineNoise(t_raw, emg_raw1);

            if (ifplot) % the raw data 
            for chi = 1:8 % what is the dimension of p_tmp and f_tmp???
                [p_tmp(chi,:), f_tmp(chi,:)] = pspectrum(emg_filter1(chi,:), 2000);
            end

            figure('name', ['line noise removal of spectrum of Session' num2str(obj.ss_num)])
            hold on; 
            for chi = 1:8
                plot(f_tmp(chi,:), log(p_tmp(chi,:)));
            end
            legend('ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8');
            
            end


            [t_filter2,emg_filter2] = removeOptNoise(t_filter1, emg_filter1, p_tmp, f_tmp, ifplot);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             for ch_i = 1:8
% %             % See the power spectrum of the raw data
%              figure(ch_i)
%              x = emg_raw(ch_i,:);
% %             x = emg_filter(7,:);
%              Fs = 1000;
%              N = length(x(1,:));
% % %             xdft = fft(emgtmp);
%              xdft = fft(x);
%              xdft = xdft(:,1:N/2+1);
%              psdx = (1/(Fs*N)) * abs(xdft).^2;
%              psdx(2:end-1) = 2*psdx(2:end-1);
%              freq = 0:Fs/length(x):Fs/2;
% % 
%              plot(freq,10*log10(psdx));
%              grid on
%              title('Periodogram Using FFT')
%              xlabel('Frequency (Hz)')
%              ylabel('Power/Frequency (dB/Hz)')
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%             % do a high-pass filter
%             x = hiaghpass(emg(1,:)', 200, Fs);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [t_filter, emg_filter] = removeMotionNoise(t_filter2, emg_filter2);
            %
            t = t_filter; 
            emg_int = emg_filter;
            emg_processed = zeros(size(emg_int));
            emg_processed0= zeros(size(emg_int));
                
            for chi = 1:8 % iterate through channels
                % 1. do the 100Hz low pass filter
                fs = obj.freq; % data frequency
                lpf1 = 20; % were60 
%                 emg1 = bandstop(emg(ch_i,:),[118.5 121.5],fs);
                emg_stp1 = highpass(emg_int(chi,:), lpf1, fs);
%                 emg_stp1 = emg(chi,:);
                % 2. mean-centered and scaled by standard deviation...
%                 emg_stp2 = (emg_stp1 - obj.emg_means(chi)) / obj.emg_stds(chi);
%                 emg_stp2 = detrend(emg_stp1);
                emg_stp2 = emg_stp1; % already done previously...
                % 3. squared, and do low-pass filter of 30Hz
%                 lpf2 = 30; %30
%                 emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
%                 emg_stp3 = emg_stp2.^2;
                % 3. squre root transform and *2
%                 emg_stp4 = real(sqrt(emg_stp2.^2));                                                   % bad name, need chagne
                emg_stp4 = abs(emg_stp2);
%                 emg_stp4 =  abs(emg(chi,:)/32767*5000/2000); % convert to mV % already mV since 2022-06-09
%                 [emg_evl, ~] = envelope(real(emg_stp4), 50, 'peak'); % check what it used to be?
%                 [emg_evl, ~] = envelope(real(emg_stp4), 10, 'peak'); % check what it used to be?
%                 [emg_evl, ~] = envelope(emg_stp4, 50, 'rms'); 
%                 [emg_evl, ~] = envelope(emg_stp4, 300, 'rms'); 
                [emg_evl, ~] = envelope(emg_stp2, 50, 'rms'); 
%                 if (~obj.if_calibSS || 1 ) % 1 when checking calibrationresult
%                     emg_evl = emg_evl/obj.amp_norm(chi)*100; % norm by maximum
%                     emg_evl = emg_evl; % norm by maximum
%                 end
%                 emg_evl = emg_evl/obj.amp_norm(chi); % norm by maximum
%                 [emg_evl, ~] = envelope(emg_stp4, 500, 'analytic'); 
%                 [emg_evlu, emg_evlp] = envelope(emg_stp4, 300, 'rms'); 
                emg_stp5 = emg_stp4; % do not do anything (rectify and envolope)
%                 [emg_stp5, ~] = envelope(real(emg_stp4), 150, 'peak'); % check what it used to be?
%                 ifplot = 1;
%                 if (ifplot)
%                     clf; hold on;
%                     plot(t, emg_int(chi,:), 'Color', [0.7 0.7 0.7], 'Marker', '.')
%                     plot(t, emg_stp1, 'r', 'Marker', '.');
%                     plot(t, emg_stp2, 'b', 'Marker', '.');
%                     plot(t, emg_stp3, 'g', 'Marker', '.');
%                     plot(t, emg_stp4, 'm', 'Marker', '.');
%                     plot(t, emg_evl, 'c', 'Marker', '.');
%                     ylim([-0, 0.2]);
%                     title(['chennel ' num2str(chi)]);
%                     xlabel('time (s)' );
%                     ylabel('magnitude (mV)');
%                     legend('raw', 'highpass', 'z-score', 'lowpass', 'sqrt transform', 'envolope');
%                 end
                ifplot = 1;
                if (ifplot)

                    clear axh
                    figure(chi);
                    axh(1) = subplot(2,1,1); hold on;
                    % plot the time-domain signal 
                    plot(t_raw, emg_raw(chi,:), ':');   % raw
                    plot(t, emg_int(chi,:), '--');       % remove line noise
                    plot(t, emg_stp1, '-');              % highpass
                    plot(t, emg_stp4, 'Color', [0.5 0.5 0.5]);              % rectify
                    plot(t, emg_evl, 'LineWidth',1);               % envelope
%                     plot(t, emg_evlu, 'k.');
%                     plot(t, emg_evlp, 'k.');
                    legend('raw', 'line noise removed', 'highpass', 'rectified', 'envelop');
                    xlabel('time (s)');
                    if (obj.mvf_shrink)
                        ylabel('intensity (%)');
                        ylim([-10 200]);
                    else
                        ylabel('value (mV)');
                    end
                    axh(2) = subplot(2,1,2); hold on;
                    % plot the frequency-domain power spectrum
                    
%                     obj.powerspectrum(emg_raw(chi,:),t,axh(2)); 
%                     obj.powerspectrum(emg_int(chi,:),t,axh(2)); 
%                     obj.powerspectrum(emg_stp1,t,axh(2)); 
%                     obj.powerspectrum(emg_stp4,t,axh(2)); 
%                     obj.powerspectrum(emg_evl,t,axh(2)); 

                    [p1, f1] = pspectrum(emg_raw(chi,:), 2000);
                    [p2, f2] = pspectrum(emg_int(chi,:), 2000);
                    [p3, f3] = pspectrum(emg_stp1, 2000);
                    [p4, f4] = pspectrum(emg_stp4, 2000);
                    [p5, f5] = pspectrum(emg_evl, 2000);
                    plot(f1, log(p1), 'r');
                    plot(f2, log(p2), 'g');
                    plot(f3, log(p3), 'b');
                    plot(f4, log(p4), 'c');
                    plot(f5, log(p5), 'm');

%                     obj.powerspectrum(emg_evlu-emg_evlp,t,axh(2)); 
                    xlim([0 1000]); % Fs==2000;
                    legend('raw', 'line noise removed', 'highpass', 'rectified', 'envelop');
                    sgtitle(['emg data for ch' num2str(chi)]);
                end
                
                emg_processed0(chi,:)=emg_stp4;
                emg_processed(chi,:) = emg_stp5;%emg_stp4;
                emg_processedenv(chi,:) = emg_evl;
                
            end
            close all;
            obj.data.emg = emg_processed;
            obj.data.emgevl=emg_processedenv;
            obj.data.t = t;
            ifplot = 0;
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

        function [t_filter, emg_filter] = removeLineNoise(obj, t, emg, ifplot)
            % might need de-trend (cg)
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            else

            end

            % The octives of the line nosie, below 1000Hz
            freq_interested = 60:60:1000;
            Fs = 2000;
            t_filter = t;
            emg_filter = zeros(8, size(t_filter,1));
            [tmp, ~] = pspectrum(emg(1,:), 2000); % only care about the size
            p_tmp_raw = zeros(8,size(tmp,1));
            f_tmp_raw = zeros(8,size(tmp,1));

            for ch_i = 1:size(emg,1)
                [p_tmp, f_tmp_raw(ch_i,:)] = pspectrum(emg(ch_i,:), 2000);
                p_tmp_raw(ch_i,:) = pow2db(p_tmp);
            end

            % 1. Show a spectrum at the begining of the peak finding.
            ...ifplot = 0;
            if(ifplot)
                figure('name', 'Sanity Check inSessionScanEMG::removeLineNoise');
                hold on;
                for ch_i = 1:size(emg,1)
                    plot(f_tmp_raw(ch_i,:), p_tmp_raw(ch_i,:));
                end
                xlabel('freuency (Hz)');
                ylabel('intensity log(P) (dB)'); % ... check! not sure about the unit!!!
                legend('ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8');
            end

            % 2. find peak and see if it is intersected with the line frequency
            locs_cell = cell(8,1);
            for ch_i = 1:size(emg,1) % each channel
                % find peaks in the spectrum
                pwdb_tmp = p_tmp_raw(ch_i,:);
                p_val_subtracted = pwdb_tmp - smooth(pwdb_tmp)';
                [~, locs] = findpeaks(p_val_subtracted, 'MinPeakHeight', 0.1);
                f_peak_loc = round(f_tmp_raw(ch_i,locs));

                if (~isempty(intersect(f_peak_loc,freq_interested)))
                    locs_cell{ch_i} = intersect(f_peak_loc,freq_interested);
                    % report and remove
                    disp('there is A line noise!');
                    freq_list = intersect(f_peak_loc,freq_interested);
                    % build up a filter based on the frequencies
                    filter_arr = [];
                    for f0_i = 1:length(freq_list)
                        f0 = freq_list(f0_i);
                        f = fdesign.notch('N,F0,BW',2,f0,3,Fs);
                        h(f0_i) = design(f);

                        if (f0_i) == 1
                            filter_arr = [filter_arr 'h(' num2str(f0_i) ')'];
                        else
                            filter_arr = [filter_arr ',h(' num2str(f0_i) ')'];
                        end
                    end
                    % concatinate the filter 
                    eval(['hd = dfilt.cascade(' filter_arr ');']);
                    emg_tmp = emg(ch_i,:);
                    emg_out = filter(hd, emg_tmp);
                    emg_filter(ch_i,:) = emg_out;
                else
                    % no line noise detected
                    emg_filter(ch_i,:) = emg(ch_i,:);
                end

                ifplot = 0;
                if (ifplot)
                    figure(); hold on;
                    plot(f_tmp_raw(ch_i,:),p_tmp_raw(ch_i,:));
                    plot(f_tmp_raw(ch_i,locs), pwdb_tmp(locs), 'rs');
                    xlabel('frequency (Hz)');
                    ylabel('relative noise intensity (dB)' );
                    title(['channel ' num2str(ch_i)]);
                    legend('raw data', 'captured peak');
                end

            end
            % calculate the filtered power spectrum
            p_tmp_filtered = zeros(size(p_tmp_raw));
            f_tmp_filtered = zeros(size(f_tmp_raw));
            for ch_i = 1:size(emg,1)
                [p_tmp_ftd, f_tmp_filtered(ch_i,:)] = pspectrum(emg_filter(ch_i,:),Fs);
                p_tmp_filtered(ch_i,:) = pow2db(p_tmp_ftd);
            end

            % after noise removal 
            ...ifplot = 0;
            if (ifplot)
                figure('name', 'after lineNoise Removal');
                hold on;
                for ch_i = 1:size(emg,1)
                    [p_tmp, f_tmp] = pspectrum(emg_filter(ch_i,:),Fs);
                    plot(f_tmp, pow2db(p_tmp));
                    xlabel('frequency (Hz)');
                    ylabel('?intensity (dB)');
                end
                xline(freq_interested);
                legend('ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8');
            end

            ifplot = 1; % the summary plot of the processing
            if (ifplot)
                figure('position', [0 0 552 1000]);
                % a 8-by-1 figure shows the power intensity of all the band (to avoid the
                % overleat)
                clear lnh;
                for ch_i = 1:8 % plot the origin power spectrum and filtered power spectrum
                    subplot(8,1,ch_i); hold on;
                    lnh(1) = plot(f_tmp_raw(ch_i,:), p_tmp_raw(ch_i,:), 'b-');
                    lnh(2) = plot(f_tmp_filtered(ch_i,:), p_tmp_filtered(ch_i,:), 'r-');
                    xlabel('frequency (Hz)');
                    ylabel('dB');
                    title(['channel ' num2str(ch_i)]);

                    if (~isempty(locs_cell{ch_i}))
                        [~, loc] = min(abs(f_tmp_raw(ch_i,:) - locs_cell{ch_i}'),[],2);
                        lnh(3) = plot(f_tmp_raw(ch_i,loc), p_tmp_raw(ch_i,loc), 'rs');
                    end
                end
                switch length(lnh)
                    case 2 % not detected peak
                        legend(lnh, {'before', 'after'});
                    case 3 % more than 1 peak(s)
                        legend(lnh, {'before', 'after', 'detected peak'});
                end
                sgtitle(['Session' num2str(obj.ss_num) ' Power spectrum of Line Noise Removal']);
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

        function axh = powerspectrum(obj, data, t, axh)
            % obj.powerspectrum(data, t, axh) 
            % plot the specific data power spectrum in the axis axh 
            % data format: m-by-n, m channels, and n datapoint in each channel
            % For multiple channels, plot each of them individually. 
            % t specifies time of data point collection (give an idea of frequency)  
            FS_raw = 1./diff(t);
            Fs = unique(round(FS_raw)); 
            ifinterp = 0;
            if (ifinterp)
                t_fmt = t(1):1/Fs:t(end); 
                dat_fmt = interp1(t, data', t_fmt'); 
                dat_fmt = dat_fmt';
            else
                dat_fmt = data;
                t_fmt = t; 
            end
            chs = size(data, 1);

            %%%%%%%%%%%%% intermediate plot
            ifplot = 0; % intermediate step getting power spectrum figure
            if (ifplot)
                
                for ch_i = 1:chs
                    axhtmp(ch_i) = subplot(chs,1,ch_i); 
                    hold on;
                    plot(t, data(ch_i,:), 'b.'); 
                    plot(t_fmt, dat_fmt(ch_i,:), 'r');
                    if (ch_i==1)
                        legend('raw', 'interp');
                    end
                    ylabel('EMG magnitude');
                    xlabel('t (s)'); 
                end
                linkaxes(axhtmp, 'x');
                sgtitle('EMG power-spectrum: re-sample sanity check')
            end
            %%%%%%%%%%%%%%
    
            % for each channel in the data
%                 axh = subplot(1,1,1); % consider?
            if (~exist('axh', 'var') || isempty(axh)); hold on;
                axh = subplot(1,1,1); % consider?
            else
                axh = subplot(axh); hold on;
            end

            % plot for each of the channels
            for ch_i = 1:chs
                x = dat_fmt(ch_i,:);
                N = length(x(1,:));

                xdft = fft(x);
                xdft = xdft(:,1:N/2+1);
                psdx = (1/(Fs*N)) * abs(xdft).^2;
                psdx(2:end-1) = 2*psdx(2:end-1);
                freq = 0:Fs/length(x):Fs/2;
                %
                plot(freq,10*log10(psdx));

            end
            legend
            grid on
            
            title('Periodogram Using FFT');
            ylabel('Power/Frequency (dB/Hz)');
            xlabel('Frequency (Hz)');


            
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
        
        function axh = plotEVLEMGData(obj)
            % axh = plotRawEMGData(obj)
            % plot the raw data of EMG
            % with BlackRock time, and the raw data
            axh = figure(); 
            axesh = zeros(1,8);
            for axi = 1:8
                axesh(axi) = subplot(8,1,axi);
                plot(obj.data.t, obj.data.emgevl(axi,:));
                ylabel(['ch' num2str(axi)] );
                switch axi
                    case 1
                        title('envloped EMG data');
                    case 8
                        xlabel('time (BlackRock Sys)');
                end
            end
            linkaxes(axesh, 'x');
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

function [amp, exist_flag] = readConfigChannelMVF(ssnum)
% find the MVF gain on each channel
% return value:
% amp = [O_gain for each channel (1~8), B_gain for each channel (1~8)];
filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetEMGMVFData.conf';
fid = fopen(filename);
C = textscan(fid, '%s\n','CommentStyle','#');
fclose(fid);
for li = 1:size(C{1}, 1)
    str = C{1}{li};
    freadtmp = textscan(str,'%f,');
    ss_num = freadtmp{1}(1);
    if ss_num ~= ssnum
        continue;
    end
    amp_str = freadtmp{1}(2:end);
end
if exist('amp_str', 'var')
    amp = double(amp_str)';
    clear amp_str;
    exist_flag = true;
else
    amp = 1*ones(1,8);    % The right order
    exist_flag = false; 
end
end

function [t_filter, emg_filter] = removeOptNoise(t, emg, p_tmp, f_tmp, ifplot)
% p_tmp is the matrix of power spectrum ;
% f_tmp is the matrix of the frequency of power spectrum; 
if (~exist('ifplot', 'var'))
    ifplot = 0;
else
    ifplot = 1; 
end
if (ifplot)
    figure(); hold on;
for chi = 1:8 
    pspectrum(emg(chi,:),2000); 
end
title('spectrum of the detrend data');
legend('ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8');
end

emg_filter = emg;
t_filter = t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A array with session number at the begining and the channels later 
chs = [6 8]; % should be session dependent 

freq_rmv = [140 160 240 260 340 360 440 460 540 560 640 660 740 760 840 860 940 960];
Fs = 2000;
ifplot = 1;
% if(ifplot)
%     figure();
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ch_i = chs 
    emg_tmp = emg(ch_i,:);
    for f0_i = 1:length(freq_rmv)
        f0 = freq_rmv(f0_i);
        f = fdesign.notch('N,F0,BW',2,f0,3,Fs);
        h(f0_i) = design(f); 
%         if f0_i == 1
%             hd = h(1);
%         else
%             hd = addstage(hd, h(f0_i));
%         end
    end
    % change the filter name, 
    % also change the filter writing way.

    hd = dfilt.cascade(h(1), h(2), h(3), h(4), h(5), h(6), h(7), ...
        h(8), h(9), h(10), h(11), h(12), h(13), h(14), h(15), h(16), h(17), h(18));

%     hfvt = fvtool(hd, 'Color', 'white'); % plot the figure of filter
    % cascade the filter 
    emg_out = filter(hd, emg_tmp); 
    
    if(ifplot)
%         clf; 
        hold on; 
        [p1, f1] = pspectrum(emg_tmp, 2000);
        [p2, f2] = pspectrum(emg_out, 2000);
        plot(f1, log(p1), 'b');
        plot(f2, log(p2), 'g-'); 
        legend('raw', 'filtered');
        title(['channel' num2str(ch_i)]);
    end
    
    emg_filter(ch_i,:) = emg_out;
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

function if_calibSS = findCalibSS(ss_num)
    % remove this into some .conf file in the future
    ss_list = [4315, 4312 ...
        4324, 4327 ...
        4335, 4338 ...
        4349, 4354 ...
        ];
    if sum(ss_num==ss_list)
        if_calibSS = true;
    else
        if_calibSS = false;
    end
end
