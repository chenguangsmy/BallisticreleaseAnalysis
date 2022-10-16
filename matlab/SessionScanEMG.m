classdef SessionScanEMG
    %SESSIONSCANEMG Summary of this class goes here
    %   Detailed explanation goes here
    %SESSIONSCANEMG scan from the intermediate mat file (which
    %   contains the EMG continuous sitnal), or the raw file and
    %   get a variable.
    %   SessionScanEMG(ss_num, ifReadAgain)
    %   Current setting is from the `intermediate` file
    %   header read from %DT sequence% described above.


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
                ifReadAgain = 0;
            end

            obj.ss_num = ss_num;
            obj.if_calibSS = findCalibSS(obj.ss_num);
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

            if (ifReadAgain || ~exist(obj.fname_fmt, 'file') || obj.if_calibSS) % If want to re-perform pre-process

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

                % convert the channel to the right order
                chmap = readConfigChannelMap(obj.ss_num);
                obj.data.emg = obj.data.emg(chmap, :);
                % magnify the channels according to the magnification.
                amps = readConfigChannelAmp(obj.ss_num);
                obj = obj.convertTomV(amps);

                [obj.amp_norm, exist_flag] = readConfigChannelMVF(obj.ss_num);

                % visualize spectrumgram to see time-related frequency
                % signal
                ifplot = 0;
                if(ifplot)
                    for ch_i = 1:8
                        figure(ch_i);
                        % plot the spectrum regard to time
                        spectrogram(obj.data.emg(ch_i,:), 5*2e3, 4*2e3, 2^13);
                        title(['channel ' num2str(ch_i)]);
                    end
                end

                % do the frequency-based processing
                % preprocess data (filter, and take the envolope)-
% %                 obj = obj.preprocessRawData(1); % do not pre-process as Federico asked so 

                % convert chanel to the MVF
% %                 obj.mvf_shrink = exist_flag;
% %                 obj = obj.normalizebyMVF(obj.amp_norm); % already

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
            emgevl_pct_mvf = obj.data.emgevl.*amp_mat;
            obj.data.emgevl = emgevl_pct_mvf;

            amp_mat = repmat((1./amp)',1,size(obj.data.emg,2));
            emg_pct_mvf = obj.data.emg.*amp_mat;
            obj.data.emg = emg_pct_mvf;

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
            %
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
            ifplot = 0;
            if (ifplot)
                clf; 
                hold on;
%                 plot(t_tmp(1)+((1:length(t_tmp))/2000), t_tmp(1)+((1:length(t_tmp))/2000), 'b.'); 
                plot(t_tmp(end)+((-length(t_tmp):-1)/2000), t_tmp(end)+((-length(t_tmp):-1)/2000), 'b.'); 
                plot(t_tmp(1)+((1:length(t_tmp))/2000), t_tmp, 'r.'); 
                xlabel('time start from last data (s)');
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
%               t_first = min(t_tmp); % should be some value, 
              t_last = max(t_tmp) + (sum(t_tmpidx)-1)*1/sample_freq; 
              t_first = t_last - (length(t_tmp)-1) * 1/sample_freq;
              t = linspace(t_first,t_last,length(t_tmp));
            ifplot = 0; 
            if (ifplot)
                clf; hold on;
                plot(idx_tmp, t_tmp, 'bo');
                plot(idx_tmp, t, 'r.');
                legend('origin', 'introp');
                xlabel('index of time');
                ylabel('actual time(s)');
                title('time acquire from *.mat and reconstruct');
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
                emg_raw1(ch_i,:) = detrend(emg_raw(ch_i,:)); % NOT nesessary for the constant direction...
%                 emg_raw1(ch_i,:) = emg_raw1(ch_i,:) - mean(emg_raw1(ch_i,:));
            end
            % remove line noise from the raw data 

            % look the power spectrum, and spectrugram for each of the
            % channel (see if there is noise throughout the session) 
            ifcheckpointOn = 1;
            if (ifcheckpointOn)
                for chi = 1:8
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
            end

            [t_filter1, emg_filter1] = obj.removeLineNoise(t_raw, emg_raw1);

            ifcheckpointOn = 0;
            if (ifcheckpointOn)
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


            [t_filter2,emg_filter2] = obj.removeOPTONoise(t_filter1, emg_filter1, ifplot);

            ifMotionNoiseRemoval = 1;
            if (ifMotionNoiseRemoval)
                [t_filter3, emg_filter3] = obj.removeMotionNoise(t_filter2, emg_filter2);
            else
                t_filter3 = t_filter2;
                emg_filter3 = emg_filter2;
            end

            ifCheckpointOn = 0;
            if (ifCheckpointOn) % the raw data
                for chi = 1:8 % what is the dimension of p_tmp and f_tmp???
                    [p_tmp(chi,:), f_tmp(chi,:)] = pspectrum(emg_filter2(chi,:), 2000);
                end

                figure('name', ['motion noise removal of spectrum of Session' num2str(obj.ss_num)])
                hold on;
                for chi = 1:8
                    plot(f_tmp(chi,:), log(p_tmp(chi,:)));
                end
                legend('ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8');
            end



            t = t_filter3; 
            emg_int = emg_filter3;
            emg_processed = zeros(size(emg_int));
            emg_processed0= zeros(size(emg_int));
                
            tic 
% %             for chi = 1:8 % iterate through channels
                % 1. do the 100Hz low pass filter
                fs = obj.freq; % data frequency
                lpf1 = 10; % were60 
%                 lpf1 = 100; % were60
                % prior to 2022-09-20 set the freq to be 100
%                 emg1 = bandstop(emg(ch_i,:),[118.5 121.5],fs);
%                 emg_stp1 = highpass(emg_int(chi,:), lpf1, fs); % highpass treated each column independently...
                emg_stp1 = highpass(emg_int', lpf1, fs, 'Impulseresponse', 'iir', 'Steepness',0.5); % matrix operation
%                 emg_stp1 = emg(chi,:);
                % 2. mean-centered and scaled by standard deviation...
%                 emg_stp2 = (emg_stp1 - obj.emg_means(chi)) / obj.emg_stds(chi);
%                 emg_stp2 = detrend(emg_stp1);
%                 emg_stp2 = emg_stp1; % already done previously...
                % 3. squared, and do low-pass filter of 30Hz
%                 lpf2 = 30; %30
%                 emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
%                 emg_stp3 = emg_stp2.^2;
                % 3. squre root transform and *2
%                 emg_stp4 = real(sqrt(emg_stp2.^2));                                                   % bad name, need chagne
%                 emg_stp4 = abs(emg_stp2);
                % prior to 2022-09-20, have abs
%                 emg_stp4 = emg_stp2;
                emg_stp4 = abs(emg_stp1);
%                 emg_stp4 =  abs(emg(chi,:)/32767*5000/2000); % convert to mV % already mV since 2022-06-09
%                 [emg_evl, ~] = envelope(real(emg_stp4), 50, 'peak'); % check what it used to be?
%                 [emg_evl, ~] = envelope(real(emg_stp4), 10, 'peak'); % check what it used to be?
%                 [emg_evl, ~] = envelope(emg_stp4, 50, 'rms'); 
%                 [emg_evl, ~] = envelope(emg_stp4, 300, 'rms'); 
%                 [emg_evl, ~] = envelope(emg_stp2, 50, 'rms'); 
%                 [emg_evl, ~] = envelope(emg_stp2, 200, 'rms');
                [emg_evl1, ~] = envelope(emg_stp1, 50, 'rms'); % envelope treated each column independently...
%                 emg_evl1 = emg_evl1';
%               previous RMS filter with RMS method 
%                 [emg_evl1, ~] = envelope(emg_stp2, 50, 'analytic');
%                 if (~obj.if_calibSS || 1 ) % 1 when checking calibrationresult
%                     emg_evl = emg_evl/obj.amp_norm(chi)*100; % norm by maximum
%                     emg_evl = emg_evl; % norm by maximum
%                 end

%%%%%%%%%%%%%%%%%%%%%%%%%
                % implement a buterworth filter
                Freq = fs;
                NyqFreq = Freq/2; 

%                 LP = 0.5;
%                 Wp = LP/NyqFreq; 
%                 [F,E] = butter(4, Wp, 'low'); % butter worthfilter with 0.5Hz
%                 emg_evl_bt05 = filtfilt(F,E,abs(emg_stp2)); 

                LP = 2;
                Wp = LP/NyqFreq; 
                [F,E] = butter(4, Wp, 'low'); % butter worthfilter with 2Hz
%                 emg_evl_bt2 = filtfilt(F,E,abs(emg_stp2)); 
                emg_rms_bt2 = filtfilt(F,E,emg_evl1); 

                LP = 20;
                Wp = LP/NyqFreq; 
                [F,E] = butter(2, Wp, 'low'); % butter worthfilter with 10Hz
                emg_rms_bt10 = filtfilt(F,E,emg_evl1); 

% 
%                 LP = 5;
%                 Wp = LP/NyqFreq; 
%                 [F,E] = butter(4, Wp, 'low'); % butter worthfilter with 10Hz

%                 F = [
% 
%                 1.0000    2.0000    1.0000    1.0000   -1.9949    0.9959
%                 1.0000    2.0000    1.0000    1.0000   -1.9868    0.9879
%                 1.0000    2.0000    1.0000    1.0000   -1.9790    0.9801
%                 1.0000    2.0000    1.0000    1.0000   -1.9716    0.9726
%                 1.0000    2.0000    1.0000    1.0000   -1.9646    0.9657
%                 1.0000    2.0000    1.0000    1.0000   -1.9582    0.9593
%                 1.0000    2.0000    1.0000    1.0000   -1.9525    0.9536
%                 1.0000    2.0000    1.0000    1.0000   -1.9476    0.9486
%                 1.0000    2.0000    1.0000    1.0000   -1.9434    0.9445
%                 1.0000    2.0000    1.0000    1.0000   -1.9401    0.9412
%                 1.0000    2.0000    1.0000    1.0000   -1.9378    0.9388
%                 1.0000    2.0000    1.0000    1.0000   -1.9363    0.9373
%                 1.0000    1.0000         0    1.0000   -0.9679         0];
% 
%                 E = [
% 
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0003
%                 0.0160
%                 1.0000];

            SOS =[
                1.000000000000000   2.000000000000000   1.000000000000000   1.000000000000000  -1.993884619309987   0.994421634925319
                1.000000000000000   2.000000000000000   1.000000000000000   1.000000000000000  -1.983143383413954   0.983677506077834
                1.000000000000000   2.000000000000000   1.000000000000000   1.000000000000000  -1.973448250544566   0.973979762005324
                1.000000000000000   2.000000000000000   1.000000000000000   1.000000000000000  -1.965326206532546   0.965855530472311
                1.000000000000000   2.000000000000000   1.000000000000000   1.000000000000000  -1.959206654373552   0.959734330126138
                1.000000000000000   2.000000000000000   1.000000000000000   1.000000000000000  -1.955405908813222   0.955932560905914
                1.000000000000000   1.000000000000000                   0   1.000000000000000  -0.977058660983481                   0];
            G = [0.000134253903833
                0.000133530665970
                0.000132877865190
                0.000132330984941
                0.000131918938146
                0.000131663023173
                0.011470669508260
                1.000000000000000];

%                 emg_evl_bt10 = filtfilt(SOS,G,abs(emg_stp4)); 

%                 emg_evl_bt10 = filtfilt(SOS,G,abs(emg_evl1)); 

%%%%%%%%%%%%%%%%%%%%%%%%%
%                 emg_evl = emg_evl/obj.amp_norm(chi); % norm by maximum
%                 emg_evl = emg_evl_bt05;
%                     emg_evl = emg_evl_bt10;
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
%                 ifplot = 1;
                if (ifplot)

                    clear axh
                    figure(chi);
                    axh(1) = subplot(2,1,1); hold on;
                    % plot the time-domain signal 
%                     plot(t_raw, emg_raw(chi,:), ':');   % raw
%                     plot(t, emg_int(chi,:), '--');       % remove line noise
%                     plot(t, emg_stp1, '-');              % highpass
                    plot(t, emg_stp4, 'Color', [0.5 0.5 0.5]);              % rectify
%                     plot(t, emg_evl, 'LineWidth',1);               % envelope
%                     plot(t, emg_evl_bt05, 'LineWidth',1);               % envelope
%                     plot(t, emg_evl1, 'LineWidth',2);               % envelope
%                     plot(t, emg_evl_bt10, 'LineWidth',1);               % envelope
%                     plot(t, emg_evlu, 'k.');
%                     plot(t, emg_evlp, 'k.');
                    plot(t, emg_evl1, 'LineWidth', 1);
%                     plot(t-t(1), emg_evl_bt2_rms, 'LineWidth', 1);
                    plot(t, emg_rms_bt2, 'LineWidth', 1);
%                     legend('raw', 'line noise removed', 'highpass', 'rectified', 'envelop');
                    legend('rectified', 'butterworth-filtered', 'rms-filtered', 'rms*btw-2Hz');
                    xlabel('time (s)');
                    if (obj.mvf_shrink) % this line should not be here!!!
%                         plot(t, emg_evl./obj.amp_norm(chi)*100, 'LineWidth', 1);
                        ylabel('intensity (portion)');
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
%                     [p5, f5] = pspectrum(emg_evl, 2000);
%                     [p6, f6] = pspectrum(emg_evl_bt05, 2000);
%                     [p6, f6] = pspectrum(emg_evl_bt10, 2000);
%                     plot(f1, log(p1), 'r');
%                     plot(f2, log(p2), 'g');
%                     plot(f3, log(p3), 'b');
                    plot(f4, log(p4), 'c');
%                     plot(f5, log(p5), 'm');
%                     plot(f6, log(p6), 'k');

%                     obj.powerspectrum(emg_evlu-emg_evlp,t,axh(2)); 
                    xlim([0 1000]); % Fs==2000;
%                     legend('raw', 'line noise removed', 'highpass', 'rectified', 'envelop', 'envelop-btw0.5Hz');
                    legend('raw', 'line noise removed', 'highpass', 'rectified', 'envelop', 'envelop-btw10Hz');
                    sgtitle(['emg data for ch' num2str(chi)]);
                end
                toc 
                
%                 ifplot = 0;
%                 if (ifplot)
%                     figure(); hold on;
%                     plot(t, emg_stp4, ':');   % raw
%                     plot(t, emg_evl_bt05, 'LineWidth',1);               % envelope
%                 end

%                 emg_processed0(chi,:)=emg_stp4;
%                 emg_processed(chi,:) = emg_stp5;%emg_stp4;
%                 emg_processedenv(chi,:) = emg_evl;
%                 emg_processedenv(chi,:) = emg_evl_bt05;
%                 emg_processedenv(chi,:) = emg_evl_bt10;
%                emg_processedenv(chi,:) = emg_rms_bt2;
                emg_processed = emg_stp4; % only rectified
                emg_processedenv = emg_rms_bt10;
% %             end
%             close all;
            obj.data.emg = emg_processed';
            obj.data.emgevl2=emg_rms_bt2';
            obj.data.emgevl10=emg_rms_bt10';
            obj.data.emgraw = emg_int;      % a data on the intermediate after the filter the line noise
            obj.data.emgevl = emg_rms_bt10';     % I prefer higher frequency.
            obj.data.t = t';
%             ifplot = 1;
            if (ifplot)
                clf;
                % plot

                axh(2) = subplot(5,1,2); hold on;
%                 plot(t,emg_processed0(1:2,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(:,1)); %title('EMG12'); %ylabel('N');
                plot(t,-emg_processed(:,2)); %title('EMG12'); %ylabel('N');
                plot(t,emg_processedenv(:,1), 'linewidth', 2);
                plot(t,-emg_processedenv(:,2), 'linewidth', 2);
                legend('1', '2', '1env', '2env');
                
                axh(3) = subplot(5,1,3); hold on;
%                 plot(t,emg_processed0(3:4,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(:,3)); %title('EMG12'); %ylabel('N');
                plot(t,-emg_processed(:,4)); %title('EMG12'); %ylabel('N');
                plot(t,emg_processedenv(:,3), 'linewidth', 2);
                plot(t,-emg_processedenv(:,4), 'linewidth', 2);
                legend('3', '4', '1env', '2env');
                
                % hold on; plot(t,emg_processed(3,:));
                % plot(t,emg_processed0(3,:)); %title('EMG34'); %ylabel('N');
                
                axh(4) = subplot(5,1,4); hold on;
%                 plot(t,emg_processed0(5:6,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(:,5)); %title('EMG12'); %ylabel('N');
                plot(t,-emg_processed(:,6)); %title('EMG12'); %ylabel('N');
                plot(t,emg_processedenv(:,5), 'linewidth', 2);
                plot(t,-emg_processedenv(:,6), 'linewidth', 2);
                legend('5', '6', '1env', '2env');
                
                axh(5) = subplot(5,1,5); hold on;
%                 plot(t,emg_processed0(7:8,:), '.'); %title('EMG12'); %ylabel('N');
                plot(t,emg_processed(:,7)); %title('EMG12'); %ylabel('N');
                plot(t,-emg_processed(:,8)); %title('EMG12'); %ylabel('N');
                plot(t,emg_processedenv(:,7), 'linewidth', 2);
                plot(t,-emg_processedenv(:,8), 'linewidth', 2);
                legend('7', '8', '1env', '2env');
                
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

%             ifplot = 1; % the summary plot of the processing
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

        function [t_filter, emg_filter] = removeOPTONoise(obj, t, emg, ifplot)
            % might need de-trend (cg)
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            else

            end

            % The octives of the line nosie, below 1000Hz
            freq_interested_center = 100:200:1000;
            freq_interested = freq_interested_center + [-20 20]';
            freq_interested = freq_interested(:);

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
                f_peak_loc = round(f_tmp_raw(ch_i,locs)/10)*10; % give a ±2 
%                 f_peak_loc = round(f_tmp_raw(ch_i,locs)/10)*10 % display so that I know...
                if (~isempty(intersect(f_peak_loc,freq_interested)))
                    locs_cell{ch_i} = intersect(f_peak_loc,freq_interested);
                    % report and remove
                    disp('there is A OPTO_MARKER noise!');
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

%             ifplot = 1; % the summary plot of the processing
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
                        loc_tmp = reshape(locs_cell{ch_i}, length(locs_cell{ch_i}), 1);
                        [~, loc] = min(abs(f_tmp_raw(ch_i,:) - loc_tmp),[],2);
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

        function [t_filter, emg_filter] = removeMotionNoise(obj, t_filter, emg_filter, ifplot)
            % remove motion noise with crazy values
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            else

            end

            chs = readConfigMotionNoiseChannel(obj.ss_num); 

            emg_filter1 = emg_filter; 
            emg_filter2 = emg_filter;
            for ch_i = 1:length(chs)
                ch = chs(ch_i);
                emg_filter1(ch,:) = filloutliers(emg_filter(ch,:), 'clip', 'movmedian', [200 200]);
                t_outlair_idx = abs(emg_filter(ch,:) - emg_filter1(ch,:)) > 3*std(emg_filter1(ch,:));
                emg_filter2(ch,:) = emg_filter(ch,:);
                emg_filter2(ch,t_outlair_idx) = emg_filter1(ch,t_outlair_idx);


                if (ifplot)
                    clf; hold on;
                    plot(t_filter, emg_filter(ch,:));
                    plot(t_filter, emg_filter2(ch,:));
                    legend('before filter', 'after filter');
                    title(['filter high magnitude outlairs for channel' num2str(ch)]);
                end
            end
            emg_filter = emg_filter2;
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

function chs = readConfigMotionNoiseChannel(ssnum)
% read the channel-map, channel map to the muscles
filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetEMGMotionNoise.conf';
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
    chs_str = freadtmp{1}(2:end);
end
if exist('chs_str', 'var')    % defined by the config file
    chs = double(chs_str)';
    clear chs_str;
else % default value
    chs = [];    % no channel contains noise
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
