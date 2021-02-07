classdef SessionScanWam
    %SESSIONSCANWAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % DT sequence:
        %   time, jpOutput, jvOutput, toolPositionOutput,
        %   toolVelocityOutput, wamJTOutput, wamRDTOutput
        DOF = 4 % WAM4
        Data
        time    % time from starting? the wam
        jp      % joint position
        jv      % joint velocity 
        tp      % tool position
        tv      % tool velocity
        jt      % joint torque
        rdt     % read-time sequence
        it      % iteration_number
        cf      % perturbation force
    end
    
    methods
        function obj = SessionScanWam(ss_num)
            %SESSIONSCANWAM scan from the file and get a variable
            %   read from file, the file do not have form head. The file
            %   header read from %DT sequence% described above.
            %FTSEPERATEDAT Construct an instance of this class
            %   read file according to the data sequence
            fdir = '/Users/cleave/Documents/projPitt/Ballistic_release_data/WAM.data';
            %fname = '20210127aft00.csv';
            fname = sprintf('KingKongWAM%05d.csv', ss_num);
            filename = [fdir '/' fname];
            Data = readmatrix(filename);
            obj.Data = Data;    % consider remove this due to memory taking?
            if (size(Data,2) ~= 24)
                msg = 'Data dimension is not 24, check!';
                error(msg);
            end
            DOF = obj.DOF;
            %%%% after ss1898
            idx_time = 1;
            idx_jp   = 2         :   2-1+DOF;        % 2 : 5
            idx_jv   = 2+DOF     :   2-1+2*DOF;      % 6 : 9
            idx_tp   = 2+2*DOF   :   2+2*DOF+2;      % 10: 12
            idx_tv   = 2+2*DOF+3 :   2+2*DOF+5;      % 13: 15
            idx_jt   = 2+2*DOF+6 :   2+3*DOF+6-1;    % 16: 19
           % idx_rdt  = 2+3*DOF+6;                    % 20
            idx_cf   = 2+3*DOF+6 :   2+3*DOF+9-1;    % 20: 22
            idx_it   = 2+3*DOF+9;                    % 23
            idx_rdt  = 2+3*DOF+10;                   % 24
            %%% before ss1898
            
            
            
            try
            obj.time = Data(:,idx_time);
            obj.jp   = Data(:,idx_jp);   
            obj.jv   = Data(:,idx_jv);    
            obj.tp   = Data(:,idx_tp);     
            obj.tv   = Data(:,idx_tv);     
            obj.jt   = Data(:,idx_jt);    
            obj.rdt  = Data(:,idx_rdt);   
            obj.it   = Data(:,idx_it);
            obj.cf   = Data(:,idx_cf);
            catch
                
            end
            obj = convert0tonan_RDT(obj);
        end
        
        function obj = convert0tonan_RDT(obj)
            rdt = double(obj.rdt);
            rdt_idx = find([rdt==0]);
            rdt(rdt_idx) = NaN;
            obj.rdt = rdt;
        end
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

