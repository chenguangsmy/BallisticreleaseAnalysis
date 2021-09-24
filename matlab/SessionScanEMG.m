classdef SessionScanEMG
    %SESSIONSCANEMG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ss_num
        anin_idx = 1:8
        anin_idx_offset = 0
        dat
        unit
        time
        brtime          % blackrock time
        freq
        
        CHANNELS_NUM = 128;
        AINPUTS_NUM  = 16;
        DATA_PER_SHORT = 10;
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
            fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/intmediate';
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
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

