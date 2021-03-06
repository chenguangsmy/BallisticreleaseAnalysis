classdef SessionScanFT
    %FTSEPERATEDAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        force_origin
        force_net
        force
        torque_origin
        FTrot_M = ...   % global: x-right, y-front, z-up, FT_base x-backup, y-frontup, z-left
            [0          0           cosd(180)
            -sind(45)   cosd(45)    0
            cosd(45)    sind(45)    0];
        RDT
        FT
        elapse
        force0
    end
    
    methods
        function obj = SessionScanFT(ss_num)
            %FTSEPERATEDAT Construct an instance of this class
            %   Detailed explanation goes here
%             fdir = '/Users/cleave/Documents/projPitt/Ballistic_release_data/FT.data';
            fdir = ['data/'];
            %fname = 'KingKongFT01865.csv';
            fname = sprintf('KingKongFT%05d.csv', ss_num);
            Data = readtable([fdir '/' fname]);
            obj.force_origin = [Data.Fx' + Data.Fx0'
                                Data.Fy' + Data.Fy0'
                                Data.Fz' + Data.Fz0'];
            obj.force_net    = [Data.Fx'
                                Data.Fy'
                                Data.Fz'];
            obj.torque_origin = [Data.Tx' + Data.Tx0'
                                Data.Ty'  + Data.Ty0'
                                Data.Tz'  + Data.Tz0'];
            obj.force0 = [  Data.Fx0';
                            Data.Fy0';
                            Data.Fz0'];
            obj.RDT = [Data.RDT];           % read-time sequence
            obj.FT = [Data.FT];             % FT sequence
            obj.elapse = [Data.elapse];     % read time elapse
            obj = forceFTconvert(obj);
            %plotForce(obj)

        end
        
        function plotForceOrigin(obj)
           figure();
           SAMPLE_R = 700; % plot from every 700 data points
           plot(obj.force_origin(:,1:SAMPLE_R:end)');
           legend('x','y','z');
           xlabel('read_timepints');
           ylabel('force / N'); 
           title('original force');
        end
        function plotForce(obj, sample_r)
           figure();
           if (nargin < 2)
            sample_r = 100; % plot from every 100 data points
           end
           plot(obj.force(:,1:sample_r:end)');
           legend('x','y','z');
           xlabel('read timepints');
           ylabel('force / N'); 
           title(['converted force every ' num2str(sample_r) ' pts']);
        end
        function plotForceOffset(obj, sample_r)
           figure();
           if (nargin < 2)
            sample_r = 700; % plot from every 700 data points (now using 700Hz)
           end
           plot(obj.force0(:,1:sample_r:end)');
           legend('x','y','z');
           xlabel('seconds s');
           ylabel('force / N'); 
           title(['force offset every ' num2str(sample_r) ' pts']);
        end
        function plotForceOffsetTrial(obj)
           figure();
           force_idx = find(diff(obj.force0(1,:),1,2)~=0);
           plot(obj.force0(force_idx,1:sample_r:end)');
           legend('x','y','z');
           xlabel('trial');
           ylabel('force / N'); 
           title(['force offset every trial']);
        end
        function obj = forceFTconvert(obj) % convert from select into world axis
            obj.force = obj.FTrot_M * obj.force_origin;
        end
        function plotElapse(obj)
           figure();
           plot(obj.elapse');
        end
        function plotDiffFt(obj)
            figure();
            plot(diff(obj.FT'));
        end
        % function align the data here?
        
        function plotForceTorque(obj, sample_r)
           figure();
           if (nargin < 2)
            sample_r = 5; % plot from every 100 data points
           end
           plot([obj.force(:,1:sample_r:end)', obj.torque_origin(:,1:sample_r:end)']); 
           legend('x','y','z', 'tx', 'ty', 'tz');
           xlabel('read timepints');
           ylabel('force / N'); 
           title('original force');
        end
    end
end

