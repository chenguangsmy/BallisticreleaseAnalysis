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
        X_opt0 =    [0.0460;-0.0235;-1.4282];
        X_wam0 =    [-0.5254;0.4627;0.0012];
        A =    [-0.0391    0.9870   -0.0894
                -0.9523   -0.0141    0.3063
                 0.3232   -0.0270    0.9585];
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
            
            fdirh = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
            fnameh = sprintf('KingKongOPT.%05d.csv', ss_num);
            filename_h = [fdirh '/' fnameh];
            % 1.... read the intermediate file
            load(filename_intermediate, 'Data'); 
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
            
            
            % 2.... read the raw data file 
            % assume the raw data file has been put into right place
            
            Datah = readmatrix(filename_h);
            if (exist('fname_other', 'var'))
                Datah = readmatrix([fdirh '/' fname_other]);
            end
            % NOTE: the following code only support the 1 MARKER case !!!
            % MORE MARKER NEED TO ADAPT THIS CODE! 
            idx_pos   = 3:5;
            idx_rdt  = 1;
            idx_nmarker_idx = 2; 
            idx_time = 7;
            
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
    end
end

