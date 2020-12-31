classdef FTseperateDat
    %FTSEPERATEDAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        force_origin
        force
        torque_origin
        FTrot_M = ...   % global: x-right, y-front, z-up, FT_base x-backup, y-frontup, z-left
            [0          0           cosd(180)
            cosd(135)   cosd(45)    0
            cosd(45)    cosd(45)    0];
        seq1            % the sending sequence
        seq2            % the sending sequence2
    end
    
    methods
        function obj = FTseperateDat(inputArg1)
            %FTSEPERATEDAT Construct an instance of this class
            %   Detailed explanation goes here
            fdir = '/Users/cleave/Documents/projPitt/Ballistic_release_data/FT.data';
            fname = '20201229192124.dat';
            Data = readtable([fdir '/' fname]);
            obj.force_origin = [Data.Fx' - Data.Fx0'
                                Data.Fy' - Data.Fy0'
                                Data.Fz' - Data.Fz0'];
            obj.torque_origin = [Data.Tx' - Data.Tx0'
                                Data.Ty'  - Data.Ty0'
                                Data.Tz'  - Data.Tz0'];
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
        function plotForce(obj)
           figure();
           SAMPLE_R = 5; % plot from every 100 data points
           plot(obj.force(:,1:SAMPLE_R:end)');
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
    end
end

