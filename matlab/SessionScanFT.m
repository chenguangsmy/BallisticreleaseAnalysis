classdef SessionScanFT
    %FTSEPERATEDAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        force_origin
        force
        torque_origin
        FTrot_M = ...   % global: x-right, y-front, z-up, FT_base x-backup, y-frontup, z-left
            [0          0           cosd(180)
            -sind(45)   cosd(45)    0
            cosd(45)    sind(45)    0];
        RDT
        FT
        elapse
    end
    
    methods
        function obj = SessionScanFT(inputArg1)
            %FTSEPERATEDAT Construct an instance of this class
            %   Detailed explanation goes here
            fdir = '/Users/cleave/Documents/projPitt/Ballistic_release_data/FT.data';
            fname = 'KingKong.FT.01857.csv';
            Data = readtable([fdir '/' fname]);
            obj.force_origin = [Data.Fx' - Data.Fx0'
                                Data.Fy' - Data.Fy0'
                                Data.Fz' - Data.Fz0'];
            obj.torque_origin = [Data.Tx' - Data.Tx0'
                                Data.Ty'  - Data.Ty0'
                                Data.Tz'  - Data.Tz0'];
            obj.RDT = [Data.RDT];           % read-time sequence
            obj.FT = [Data.FT];             % FT sequence
            obj.elapse = [Data.elapse];     % read time elapse
            obj = forceFTconvert(obj);
            plotForce(obj)

        end
        
        function plotForceOrigin(obj)
           figure();
           SAMPLE_R = 1; % plot from every 100 data points
           plot(obj.force_origin(:,1:SAMPLE_R:end)');
           legend('x','y','z');
           xlabel('read_timepints');
           ylabel('force / N'); 
           title('original force');
        end
        function plotForce(obj, sample_r)
           figure();
           if (nargin < 2)
            sample_r = 5; % plot from every 100 data points
           end
           plot(obj.force(:,1:sample_r:end)');
           legend('x','y','z');
           xlabel('read timepints');
           ylabel('force / N'); 
           title('original force');
        end
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
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

