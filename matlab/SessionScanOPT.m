classdef SessionScanOPT
    %SESSIONSCANOPT read the data from the optotrack files. Either the
    %intermediate file or the raw data file 
    %   Read file from the optotrack message and the 
    
    properties
        ss_num
        marker_idx = 1:10
        data
        datah           % high resolution data
        localtime       % time from the recording start % useless if I had
        brtime          % blackrock time
        freq
%         %%%%%%%%%%%% for springs
%         X_opt0 =    [0.0460;-0.0235;-1.4282];
%         X_wam0 =    [-0.5254;0.4627;0.0012];
%         A =    [-0.0391    0.9870   -0.0894
%                 -0.9523   -0.0141    0.3063
%                  0.3232   -0.0270    0.9585];
%         %%%%%%%%%%%% for subject
%         X_opt0 =    [0.0076;-0.0297;-1.4656];
%         X_wam0 =    [-0.5306;0.4854;-0.0172];
%         A =    [-0.0175    0.9850   -0.0416
%                 -0.9492   -0.0218    0.3034
%                 0.3234   -0.0249    0.9513];
%         %%%%%%%%%%%% for horizontal mount WAM
%         X_opt0 =    [0.1739;0.0185;-1.3094];
%         X_wam0 =    [-0.4850;0.0020;0.5141];
%         A =    [    0.0253    0.8596    0.0297
%                    -0.8186    0.0128    0.2829
%                     0.3296   -0.0582    0.9329];
        %%%%%%%%%%%% for horizontal mount WAM, 'inmotion'
%         X_opt0 =    [-0.039926249838378;-0.120462853038858;-1.526446107588706];
%         X_wam0 =    [-0.484880470027251;-0.515017776586255;0.001370295584099];
%         A =    [  -0.037848635084157   1.006139173669417  -0.049302507925071
%                   -0.925253905845084   0.004383012560907   0.338479651807431
%                    0.382639418084085  -0.005081861415701   0.919003503546431]; 
        %%%%%%%%%%%%% after I touched the optotrack, after ss4296
        X_opt0 = [ -0.207025909363952;  -0.159552961960486;  -1.672379585727084];
        X_wam0 = [ -0.509554118480000;  -0.503134221202000;  -0.155901906928000];
        A = [   -0.046020953751129   1.000420282996839  -0.075813059226936
                -0.893388205043085   0.000259334777461   0.430759897297924
                 0.417832205169391  -0.015157199405532   0.850685498567284
        ];
                

        read_from_intermediate = 1;
    end
    
    methods
        function obj = SessionScanOPT(ss_num, fname_other)
            %SESSIONSCANOPT Construct an instance of this class by reading
            %all the data in.
            %   read the data of intermediate file and the *.csv file. The
            %   intermediate file do not have a hardware sync time, so it
            %   could be not accurate.
            %   Plan to do the hardware alignment latter (2021-10-22, cg)
            
            % take the ss3197 as an example
            %        ss_num = 3197
            obj.ss_num = ss_num;
            fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate';
            fname = sprintf('KingKong.%05d.mat', ss_num);
            filename_intermediate = [fdir '/' fname];
            
            if (~obj.read_from_intermediate)
                fdirh = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
                fnameh = sprintf('KingKongOPT.%05d.csv', ss_num);
                filename_h = [fdirh '/' fnameh];
            end
            % 1.... read the intermediate file
            load(filename_intermediate, 'Data');
            if_markerdata = isfield(Data.QL.Data, 'OPTOTRAK_MARKER_DATA');
            if (if_markerdata)
                recv_time = Data.QL.Headers.OPTOTRAK_MARKER_DATA.recv_time;
                sg_recv_time = Data.QL.Headers.SAMPLE_GENERATED.recv_time;
                bk_time = Data.QL.Data.SAMPLE_GENERATED.source_timestamp; % check...

                % assume the message arrives at equal time interval (which hard
                % be true)
                optmsg_time = interp1(sg_recv_time, bk_time, recv_time, 'linear', 'extrap');

                data.t = optmsg_time;
                data.rdt = double([Data.QL.Data.OPTOTRAK_MARKER_DATA.nMarkers]);
                data.x   = Data.QL.Data.OPTOTRAK_MARKER_DATA.x;
                data.y   = Data.QL.Data.OPTOTRAK_MARKER_DATA.y;
                data.z   = Data.QL.Data.OPTOTRAK_MARKER_DATA.z;
                data.inview = Data.QL.Data.OPTOTRAK_MARKER_DATA.inView;

                % as we now only use one marker
                data.x = data.x(obj.marker_idx, :);
                data.y = data.y(obj.marker_idx, :);
                data.z = data.z(obj.marker_idx, :);
                obj.data = data;
            end
            
            
            % 2.... read the raw data file
            %  .... or, read batch file from recorded
            
            if (~obj.read_from_intermediate)
                Datah = readmatrix(filename_h);
                if (exist('fname_other', 'var'))
                    Datah = readmatrix([fdirh '/' fname_other]);
                end
                % NOTE: the following code only support the 1 MARKER case !!!
                % MORE MARKER NEED TO ADAPT THIS CODE!
                idx_pos   = 3:5;
                idx_rdt  = 1;
                idx_nmarker_idx = 2;
                idx_time = 7;%
                
                %             datah.rdt = reshape(datah.rdt, 1, length(datah.rdt));
                
                %datah.t = interp1(data.rdt, data.t, datah.rdt, 'linear', 'extrap');
                datah.t = Datah(:,idx_time) + 0.0025; % Assume collection time.
                % use interp to make all the Marker on the same time
                markers_num = length(unique(Datah(:,idx_nmarker_idx)));
                datats_num = size(Datah,1);
                datah.t = Datah(Datah(:,idx_nmarker_idx) == 0, idx_time)'; % first column
                datah.rdt=Datah(Datah(:,idx_nmarker_idx) == 0, idx_rdt)';     % ... useless, need to adapt
                datah.bkt = interp1(data.rdt, data.t, datah.rdt, 'linear', 'extrap')';
                
                ifplot = 1;
                if (ifplot)
                    clf; hold on;
                    plot(data.rdt, data.t, '.'); % message time;
                    plot(datah.rdt, datah.bkt, '.'); % raw paired BK time
                    legend('message time', 'reconstructed bk time');
                end
                datah.x = nan(markers_num, datats_num/markers_num);
                datah.y = nan(markers_num, datats_num/markers_num);
                datah.z = nan(markers_num, datats_num/markers_num);
                
                % for the first marker
                %   HERE I ASSUME ALL MARKER WERE COLLECTED IN 1 TIME, WHICH
                %   COULD BE WRONG
                for marker_idx = 1:markers_num
                    datah.x(marker_idx,:) = Datah(Datah(:,idx_nmarker_idx) == (marker_idx-1), idx_pos(1))';
                    datah.y(marker_idx,:) = Datah(Datah(:,idx_nmarker_idx) == (marker_idx-1), idx_pos(2))';
                    datah.z(marker_idx,:) = Datah(Datah(:,idx_nmarker_idx) == (marker_idx-1), idx_pos(3))';
                end
                
                obj.datah = datah;
                obj.localtime = data.t;
                obj.brtime = obj.datah.t; % useless variable here
            else
                % if there is tidied up KingKong****OPT.mat, read it and
                % else, create it and save it as the KingKong*****OPT.mat
                fname_OPT = sprintf('KingKongOPT.%05d.mat', ss_num); %05d
                if exist(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' fname_OPT], 'file')
                    load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' fname_OPT]);
                    
                else
                    
                    % read from intermediate.mat with each trial
%                     num_trials = length(Data.QL.Data.OPTO_BUFFER_DATA); % could be inaccurate as some trials not record
                    
                    % also read TimeSync.mat so that I'm able to know time
                    fname = sprintf('KingKongTSync.%05d.mat', ss_num);
                    datadir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
                    dataT = load([datadir '/' fname]);
                    num_trials = max(dataT.data.eventsTrials); 
                    %time_start = min(dataT.data.eventsT(dataT.data.eventsL==5));
                    %time_end = max(dataT.data.eventsT(dataT.data.eventsL==5));
                    time_interval = diff(dataT.data.eventsT(dataT.data.eventsL==5));
                    def_interval = mode(time_interval);
                    % 1. Get the blackrock recorded time
                    datTrial_stt_end_t = nan(num_trials,2);
                    time_interTrial = cell(num_trials-1,1);
                    num_interTrial = zeros(num_trials-1,1); % no data, but still has time count
                    num_pulse = zeros(1,num_trials);
                    for trial_i = 1:num_trials
                        events_idx = (dataT.data.eventsTrials==trial_i & dataT.data.eventsL==5);
                        if (sum(dataT.data.eventsTrials==trial_i & dataT.data.eventsL==5)==0) % no data recorded
                            continue;
                        end
                        datTrial_stt_end_t(trial_i,1) = min(dataT.data.eventsT(events_idx));
                        datTrial_stt_end_t(trial_i,2) = max(dataT.data.eventsT(events_idx));
                        % stack all the data time
                        if (trial_i > 1)
                            time_interTrial{trial_i-1} = datTrial_stt_end_t(trial_i-1,2):def_interval:datTrial_stt_end_t(trial_i,1);
                            time_interTrial{trial_i-1} = time_interTrial{trial_i-1}(2:end-1); % remove the begining and the last
                            num_interTrial(trial_i-1) = length(time_interTrial{trial_i-1});
                        end
                       num_pulse(trial_i) = sum(events_idx);
                    end
                    
                    % 2. Get the data out
                    dataint_idx_all = Data.QL.Data.TRIAL_CONFIG.trial_no;
                    for trial_i = 1:length(dataint_idx_all)
                        trial_idx = dataint_idx_all(trial_i);
                        events_idx = (dataT.data.eventsTrials==trial_idx & ...
                            dataT.data.eventsL==5);
                        
                        if (sum(events_idx) == 0) | sum(trial_i==dataint_idx_all)==0
                            datah.t{trial_i} = [];
                            datah.x{trial_i} = [];
                            bufferdat{trial_idx} = []; % valid?
                            continue;
                        end
                        events_t = dataT.data.eventsT(events_idx);
                        [n_markers, buffertmp] = obj.getBufferDatafromRawReading(Data.QL.Data.OPTO_BUFFER_DATA{dataint_idx_all(trial_i)}); %
                        bufferdat{trial_idx} = buffertmp;
                        num_data(trial_idx) = length(bufferdat{trial_idx}{1});
                    end


                    % 3. align them
                    for trial_i = 1:length(dataint_idx_all)
                        trial_idx = dataint_idx_all(trial_i);
                        datah.x{trial_idx} = bufferdat{trial_idx};
                        if isempty(bufferdat{trial_idx})
                            continue;
                        end
%                       datah.t{trial_i} = 1:size(bufferdat{trial_i}{1},2); % don't mean anything, time at blackrock time
                        events_idx = (dataT.data.eventsTrials==trial_idx & ...
                            dataT.data.eventsL==5);
                        events_t = dataT.data.eventsT(events_idx);

                        datah.t{trial_idx} = events_t;
                        % need some way to get the datah.t here
                        
                        % if datah.x and datah.t are not same length
                        if (size(datah.x{trial_idx}{1},2) ~= length(events_t))
                            num_lost = length(events_t)-size(datah.x{trial_idx}{1},2);
                            if (num_lost) < 0 
                                disp('Data Error: pulse < data'); 
                                % after I check the code, I found the
                                % majority reason is the lost of pulse
                                % (might be electrical problem let it
                                % unrecognized), I need to change cable. 
                                
                                disp('doing correction');
                                dur_interpulse = mode(diff(events_t));% 
                                events_t_edt = events_t; 
                                pulse_insert_idx = find(diff(events_t_edt)>dur_interpulse*1.1); 
                                while (~isempty(pulse_insert_idx))
                                    % get how many pulses was not count 
                                    diff_time = diff(events_t_edt);
                                    dur_time = diff_time(pulse_insert_idx(1));
                                    insert_num = round(dur_time/dur_interpulse);
                                    % get the pulses needed to be inserted 
                                    insert_arr = events_t_edt(pulse_insert_idx(1)):dur_interpulse:events_t_edt(pulse_insert_idx(1)+1); 
                                    insert_arr = insert_arr(2:end);  
                                    if (insert_arr(end) == events_t_edt(pulse_insert_idx(1)+1))
                                        insert_arr = insert_arr(1:end-1);
                                    end
                                    % insert the pulses 
                                    events_t_edt = [events_t_edt(1:pulse_insert_idx(1)), ...
                                        insert_arr,...
                                        events_t_edt(pulse_insert_idx(1)+1:end)];
                                    ifplot = 1;
                                    if (ifplot)
                                        clf; 
                                        plot(diff(events_t_edt), '.');
                                    end
                                    pulse_insert_idx = find(diff(events_t_edt)>dur_interpulse*1.1);
                                end
                                events_t = events_t_edt;
                                num_lost = length(events_t)-size(datah.x{trial_idx}{1},2);
                            end
                            
                                disp(['num_pulse - num_data == ' num2str(num_lost)]);
                            
                            % deal with the most daunting part, assuming the data
                            % only lost at the start of the pulse (pulse start but
                            % recording did not start).
                            % discard/dispose the early pulses
%                             datah.t{trial_i} = events_t((num_lost+1):end);
                            datah.t{trial_idx} = events_t(1:end-num_lost); % discard/dispose the late pulses
                        end
                    end
                    
                    % 3. Concatinate each trial's data into a long datah
                    dataH.t = [];       % use H to represent the whole session
                    dataH.x = [];       % n_markers, initial as 10
%                     for marker_i = 1:10
%                         dataH.x{marker_i} = [];
%                     end
                    for trial_i = 1:length(dataint_idx_all)
                        trial_idx = dataint_idx_all(trial_i);
                        dataH.t = [dataH.t datah.t{trial_idx} ];
                        % 
                        ifplot = 1; 
                        if (ifplot)
                            plot(dataH.t, '.'); 
                        end
%                         try 
                            % convert markers array to a 3d array
%                             for marker_i = 1:10
%                                 dataH.x = [dataH.x datah.x{trial_i}{marker_i} ];
%                             end
                              % format a 3-d matrix for each trial
                              if isempty(datah.x{trial_idx})
                                  data_length = 0;
                              else
                                  data_length = size(datah.x{trial_idx}{1},2);
                              end
                              dataH_allm = nan(3,data_length,10);
                              if isempty(datah.x{trial_idx})
                                  % do nothing
                              else
                                  for marker_i = 1:3 % only 3 marker used
                                      dataH_allm(:,:,marker_i) = datah.x{trial_idx}{marker_i};
                                  end
                              end

                              dataH.x = cat(2,dataH.x, dataH_allm);
%                         catch % blank in datah.x{trial_i}
%                             dataH.x = [dataH.x datah.x{trial_i} ];
%                         end
                        
                        if (trial_idx < num_trials)
                            dataH.t = [dataH.t time_interTrial{trial_idx}];
%                             for marker_i = 1:10
%                                 dataH.x{marker_i} = [dataH.x{marker_i} nan(3, num_interTrial(trial_i))];
%                             end
                              dataH.x = cat(2,dataH.x, nan(3, num_interTrial(trial_idx),10));
                        end
                    end
                    
                    
                    ifplot = 1;
                    if (ifplot)
                        clf; hold on;
                        for marker_i = 1:3
                            subplot(3,1,marker_i);
                            plot(dataH.t, dataH.x(:,:,marker_i)', '.'); % message time;
                            legend('x', 'y', 'z');
                            title(['marker' num2str(marker_i)]);
                        end
                    end
                    
                    clear datah
                    for marker_i = 1:10
                        datah.x(marker_i,:) = dataH.x(1,:,marker_i);
                        datah.y(marker_i,:) = dataH.x(2,:,marker_i);
                        datah.z(marker_i,:) = dataH.x(3,:,marker_i);
                    end
                    datah.t = dataH.t;
                    % save to file after prrocessing
                    save(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' fname_OPT], ...
                        'datah');
                end
                
                
                obj.datah = datah;
                obj.localtime = datah.t;
                obj.brtime = datah.t; % useless variable here
                obj.datah.t = datah.t;
                
                
            end

            % diagnize the nan
            nan_idx_x = abs(obj.datah.x)>9;
            nan_idx_y = abs(obj.datah.y)>9;
            nan_idx_z = abs(obj.datah.z)>9;
            obj.datah.x(nan_idx_x) = nan;
            obj.datah.y(nan_idx_y) = nan;
            obj.datah.z(nan_idx_z) = nan;
            
            obj = obj.rotate();  % change into robot axis
            
        end

        
        function axh = plotRawEMGData(obj) % consider shift this function into the opt function
            % axh = plotRawEMGData(obj)
            % plot the raw data of EMG
            % with BlackRock time, and the raw data
            axh = figure(); 
            axesh = zeros(1,8);
            for axi = 1:8
                axesh(axi) = subplot(8,1,axi);
                plot(obj.brtime, (obj.dath.x));
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
        
        function fh = plotInView(obj) % plot the inview situation of the markers
            % axh = plotInView(obj) 
            % when inview, plot as 1
            
            fh = figure();
            axh = zeros(1,10); 
            for axi = 1:10
                axh(axi) = subplot(10,1,axi);
                plot(obj.data.t,obj.data.inview(axi,:), '.');
                ylim([-0.1 1.1]);
                ylabel('inview');
            end
            sgtitle('in view of all markers');
            xlabel('time'); 
        end
        
        function obj = rotate(obj,inputArg)
            %METHOD1 align the data into certain direction 
            for marker_i = 1:size(obj.datah.x,1)
                pos = [ obj.datah.x(marker_i,:);
                        obj.datah.y(marker_i,:);
                        obj.datah.z(marker_i,:);];
                pos1= obj.X_wam0 + obj.A*(pos-obj.X_opt0);
                obj.datah.x(marker_i,:) = pos1(1,:);
                obj.datah.y(marker_i,:) = pos1(2,:);
                obj.datah.z(marker_i,:) = pos1(3,:);
            end
            % do the shift and the linear transformation (only rotate). 
        end
        
        function [n_markers, bufferdat] = getBufferDatafromRawReading(obj, data_uint8) % Dealing with int data to the exact number
            % first 8 - number
            % The others, n_markers * [x, y, z] position
            
            dat_arr = zeros(1,length(data_uint8)/8);
            unit8_mat = repmat('0', [length(data_uint8)/8, 64]);
            for byte_i = 1:length(data_uint8)/8
%                 uint8_arr = [];
                offset = (byte_i-1)*8;
                for dat_pos_i = 1:8
                    %    uint8_arr = [uint8_arr dec2bin(typecast(int8(data_uint8(dat_pos_i)),'uint8'),8)]
%                     uint8_arr = [dec2bin(typecast(int8(data_uint8(offset+dat_pos_i)),'uint8'),8) uint8_arr];
                    unit8_mat(byte_i,((8-dat_pos_i)*8 + (1:8))) = dec2bin(typecast(int8(data_uint8(offset+dat_pos_i)),'uint8'),8);
                end
%                 q = quantizer('double');
% %                 B = bin2num(q, uint8_arr);
%                 B = bin2num(q, unit8_mat);
%                 dat_arr(byte_i) = B;
            end
                q = quantizer('double');
%                 B = bin2num(q, uint8_arr);
                B = bin2num(q, unit8_mat);
%                 dat_arr(byte_i) = B;
                dat_arr = B;


            n_markers = dat_arr(1); 
            pos = reshape(dat_arr(2:end), [3, length(dat_arr(2:end))/3]);
            plot(pos');     % position of all markers
            
            num_eachMarker = size(pos,2)/n_markers; 
%             if mod(num_eachMarker, n_markers) ~= 0        % what do you mean here? Useless
%                 disp('error on the marker_num or intermediate data');
%             end
            
            bufferdat = cell(n_markers,1);
            for marker_i = 1:n_markers
                idx_start = (marker_i-1)*num_eachMarker + 1;
                idx_end = marker_i*num_eachMarker;
                bufferdat{marker_i} = pos(:,idx_start:idx_end); % x, y, z
            end
            
        end
        
    end
end

