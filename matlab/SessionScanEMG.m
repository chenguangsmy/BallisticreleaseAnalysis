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
            obj.ss_num = ss_num;
            fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate';
%             fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
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
            obj = obj.readRawData();
            
            % preprocess data
            obj = obj.preprocessRawData();
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
        
        function time = getRTMAtimeAligned(obj)
            % align the task time with the current data
            % so that I can directly plot the current data time with data
            disp('NOT FINISHED GETRTMATIMEALIGNED YET! ERROR :)))');
            time = [];
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
        
        function obj = preprocessRawData(obj)
            emg = obj.data.emg;
            emg_processed = zeros(size(emg));
            emg_processed0= zeros(size(emg));
            t = obj.data.t;
            for chi = 1:8 % iterate through channels
                % 1. do the 100Hz low pass filter
                fs = obj.freq; % data frequency
                lpf1 = 60;%100;
                % emg1 = bandstop(emg(ch_i,:),[118.5 121.5],fs);
%                 emg_stp1 = highpass(emg(chi,:), lpf1, fs);
                emg_stp1 = emg(chi,:);
                % 2. mean-centered and scaled by standard deviation...
                emg_stp2 = (emg_stp1 - obj.emg_means(chi)) / obj.emg_stds(chi);
                % 3. squared, and do low-pass filter of 30Hz
                lpf2 = 30; %30
%                 emg_stp3 = lowpass(emg_stp2.^2, lpf2, fs);
                emg_stp3 = emg_stp2.^2;
                % 3. squre root transform and *2
                emg_stp4 = sqrt(emg_stp3)*2;                                                   % bad name, need chagne
                emg_stp4 =  abs(emg(chi,:)/32767*5000/2000); % convert to mV
                [emg_stp5, ~] = envelope(real(emg_stp4), 50, 'peak'); % check what it used to be?
                ifplot = 1;
                if (ifplot)
                    clf; hold on;
                    plot(t, emg(chi,:), 'Color', [0.7 0.7 0.7], 'Marker', '.')
                    plot(t, emg_stp1, 'r', 'Marker', '.');
                    plot(t, emg_stp2, 'b', 'Marker', '.');
                    plot(t, emg_stp3, 'g', 'Marker', '.');
                    plot(t, emg_stp4, 'm', 'Marker', '.');
                    plot(t, emg_stp5, 'c', 'Marker', '.');
                    ylim([-5, 10]);
                    title(['chennel ' num2str(chi)]);
                    legend('raw', 'highpass', 'z-score', 'lowpass', 'sqrt transform', 'envolope');
                end
                emg_processed0(chi,:)=emg_stp4;
                emg_processed(chi,:) = emg_stp5;%emg_stp4;
                
            end
            obj.data.emg = emg_processed;
            ifplot = 1;
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

