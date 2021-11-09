classdef SessionScanOPT
    %SESSIONSCANOPT read the data from the optotrack files. Either the
    %intermediate file or the raw data file 
    %   Read file from the optotrack message and the 
    
    properties
        ss_num
        marker_idx = 1
        data
        datah           % high resolution data
        localtime       % time from the recording start % useless if I had
        brtime          % blackrock time
        freq
        
    end
    
    methods
        function obj = SessionScanOPT(ss_num)
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
            
            % as we now only use one marker
            data.x = data.x(obj.marker_idx, :);
            data.y = data.y(obj.marker_idx, :);
            data.z = data.z(obj.marker_idx, :);
            obj.data = data;
            
            % 2.... read the raw data file 
            % assume the raw data file has been put into right place
            
            Datah = readmatrix(filename_h);
            % NOTE: the following code only support the 1 MARKER case !!!
            % MORE MARKER NEED TO ADAPT THIS CODE! 
            idx_pos   = 3:5;
            idx_rdt  = 1;
            datah.rdt=Datah(:,idx_rdt);
            datah.rdt = reshape(datah.rdt, 1, length(datah.rdt));
            datah.t = interp1(data.rdt, data.t, datah.rdt, 'linear', 'extrap');
            datah.x = Datah(:,idx_pos(1))';
            datah.y = Datah(:,idx_pos(2))';
            datah.z = Datah(:,idx_pos(3))';
            
            obj.datah = datah;
            obj.localtime = data.t;
            obj.brtime = obj.datah.t; % useless variable here
            
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
        
        
        function obj = rotate(obj,inputArg)
            %METHOD1 align the data into certain direction 
            
            % do the shift and the linear transformation (only rotate). 
        end
    end
end

