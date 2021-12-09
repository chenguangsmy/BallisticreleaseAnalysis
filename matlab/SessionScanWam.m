classdef SessionScanWam
    %SESSIONSCANWAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % DT sequence:
        %   time, jpOutput, jvOutput, toolPositionOutput,
        %   toolVelocityOutput, wamJTOutput, wamRDTOutput
        DOF = 4         % WAM4
        ss_num
        Data            % raw data reading % Removed to save space
        time            % time from starting? the wam
        jp              % joint position
        jv              % joint velocity 
        tp              % tool position
        tv              % tool velocity
        jt              % joint torque
        rdt             % read-time sequence
        it              % iteration_number
        cf              % perturbation force
        state           % robot state
        Data_pert       % perturbation dataset, have:
        Data_pert_ensemble
                        %       datMat  % perturbation data Frame
                        %       FT      % force Threshold
                        %       x0      % endpoint
    end
    
    methods
        function obj = SessionScanWam(ss_num)
            %SESSIONSCANWAM scan from the file and get a variable
            %   read from file, the file do not have form head. The file
            %   header read from %DT sequence% described above.
            %FTSEPERATEDAT Construct an instance of this class
            %   read file according to the data sequence
            obj.ss_num = ss_num;
            fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
            %fdir = ['data/'];
            %fname = '20210127aft00.csv';
            fname = sprintf('KingKongWAM%05d.csv', ss_num);
            filename = [fdir '/' fname];
            if ~exist(filename, 'file')
                % try to convert from binary file...
                fnametmp = sprintf('KingKongWAM%05d.bin', ss_num);
                filenametmp = [fdir '/' fnametmp];
                exportWambin2csv(filenametmp,filename)
            end
            Data = readmatrix(filename);
%             obj.Data = Data;    % consider remove this due to memory taking?
            %if (size(Data,2) ~= 24)
            %    msg = 'Data dimension is not 24, check!';
            %    error(msg);
            %end
            DOF = obj.DOF;
            if (ss_num>=1898)
            %%%% after ss1898
                idx_time = 1;
                idx_jp   = 2         :   2-1+DOF;        % 2 : 5
                idx_jv   = 2+DOF     :   2-1+2*DOF;      % 6 : 9
                idx_tp   = 2+2*DOF   :   2+2*DOF+2;      % 10: 12
                idx_tv   = 2+2*DOF+3 :   2+2*DOF+5;      % 13: 15
                idx_jt   = 2+2*DOF+6 :   2+3*DOF+6-1;    % 16: 19
                idx_cf   = 2+3*DOF+6 :   2+3*DOF+9-1;    % 20: 22
                idx_it   = 2+3*DOF+9;                    % 23
                idx_rdt  = 2+3*DOF+10;                   % 24
                if(size(Data,2) >= 25)
                    idx_state  = 2+3*DOF+11;             % 25 Pull out state
                end
            %%% before ss1898
            else
                idx_time = 1;
                idx_jp   = 2         :   2-1+DOF;        % 2 : 5
                idx_jv   = 2+DOF     :   2-1+2*DOF;      % 6 : 9
                idx_tp   = 2+2*DOF   :   2+2*DOF+2;      % 10: 12
                idx_tv   = 2+2*DOF+3 :   2+2*DOF+5;      % 13: 15
                idx_jt   = 2+2*DOF+6 :   2+3*DOF+6-1;    % 16: 19
                idx_rdt  = 2+3*DOF+6;                    % 20
            end
            
            
            try
            obj.time = Data(:,idx_time)';
            obj.jp   = Data(:,idx_jp)';   
            obj.jv   = Data(:,idx_jv)';    
            obj.tp   = Data(:,idx_tp)';     
            obj.tv   = Data(:,idx_tv)';     
            obj.jt   = Data(:,idx_jt)';    
            obj.rdt  = Data(:,idx_rdt)';   
            obj.it   = Data(:,idx_it)';
            obj.cf   = Data(:,idx_cf)';
            if(exist('idx_state','var'))
                obj.state = Data(:,idx_state)';
            end
            catch
                
            end
        end
        
        function obj = concatinateTrials2File(obj, tarL_list, fTh_list, rdt_ranges_all)
            % Define which pararamiters are at which index
%             figure;
                        
            % tarL_list, fTh_list, rdt_ranges, have same length
            saveSepFile_flag = 1;   % if a seperate file needed
            for conditioni = 1:length(tarL_list)
                rdt_ranges = rdt_ranges_all{conditioni};
                % get rdt_idx from rdt_ranges
                rdt_idx = [];
                rdt_ranges = rdt_ranges_all{conditioni};
                for pair_i = 1:size(rdt_ranges,2)
                    rdt_idx = [rdt_ranges(1,pair_i)+2*500: ...
                               rdt_ranges(2,pair_i)];
                    
                    Data_pert(conditioni,pair_i).FT     = fTh_list(conditioni);
                    Data_pert(conditioni,pair_i).x0     = tarL_list(conditioni);
                    % Data_pert(conditioni).datMat = obj.Data(rdt_idx,:);
                    Data_pert(conditioni,pair_i).time = obj.time(rdt_idx);
                    Data_pert(conditioni,pair_i).tp = obj.tp(rdt_idx,:);
                    Data_pert(conditioni,pair_i).tv = obj.tv(rdt_idx,:);
                    Data_pert(conditioni,pair_i).cf = obj.cf(rdt_idx,:);
                    Data_pert(conditioni,pair_i).it = obj.it(rdt_idx,:);
                    Data_pert(conditioni,pair_i).rdt = obj.rdt(rdt_idx,:);
                    
                    %                 tmp = Data_pert(conditioni,pair_i).tp(:,2);
                    %                 plot(tmp);hold on;
                    %                 plot([1,1].*length(tmp)-500*5,[0, 1],'-k','linewidth',2.5); xlim([0 length(tmp)]); ylim([0.49 0.5]);hold off;
                    %                 pause();

                end
            end

            obj.Data_pert = Data_pert;
            
        end
        
        function obj = concatinateTrials2File_ensemble(obj, tarL_list, fTh_list, rdt_ranges_all, rdt_mov_all)
            % Define which pararamiters are at which index
%             figure;
                        
            % tarL_list, fTh_list, rdt_ranges, have same length
            saveSepFile_flag = 1;   % if a seperate file needed
            for conditioni = 1:length(tarL_list)
                rdt_ranges = rdt_ranges_all{conditioni};
                % get rdt_idx from rdt_ranges
                rdt_idx = [];
                rdt_ranges = rdt_ranges_all{conditioni};
                rdt_mov = rdt_mov_all{conditioni};
                for pair_i = 1:size(rdt_ranges,2)
                    rdt_idx = [rdt_ranges(1,pair_i):... % Added 5 to take every 5th sample
                                        rdt_ranges(2,pair_i)];
                
%                 Data_pert_ensemble(conditioni,pair_i).FT     = fTh_list(conditioni);
%                 Data_pert_ensemble(conditioni,pair_i).x0     = tarL_list(conditioni);
% %                 Data_pert(conditioni).datMat = obj.Data(rdt_idx,:);
%                 Data_pert_ensemble(conditioni,pair_i).time = obj.time(rdt_idx) - obj.time(rdt_mov(pair_i));
%                 Data_pert_ensemble(conditioni,pair_i).tp = obj.tp(rdt_idx,:);
%                 Data_pert_ensemble(conditioni,pair_i).tv = obj.tv(rdt_idx,:);
%                 Data_pert_ensemble(conditioni,pair_i).cf = obj.cf(rdt_idx,:);
%                 Data_pert_ensemble(conditioni,pair_i).it = obj.it(rdt_idx,:);
%                 Data_pert_ensemble(conditioni,pair_i).rdt = obj.rdt(rdt_idx,:);
%                 Data_pert_ensemble(conditioni,pair_i).rdt = obj.state(rdt_idx);

                Data_pert_ensemble(conditioni).FT     = fTh_list(conditioni);
                Data_pert_ensemble(conditioni).x0     = tarL_list(conditioni);
                
                Data_pert_ensemble(conditioni).time_r(pair_i,:) = obj.time(rdt_idx) - obj.time(rdt_mov(pair_i));
                
                Data_pert_ensemble(conditioni).z_r_1(pair_i,:) = obj.tp(rdt_idx,1);
                Data_pert_ensemble(conditioni).z_r_2(pair_i,:) = obj.tp(rdt_idx,2);
                
                Data_pert_ensemble(conditioni).u_r_1(pair_i,:) = obj.cf(rdt_idx,1);
                Data_pert_ensemble(conditioni).u_r_2(pair_i,:) = obj.cf(rdt_idx,2);
                
                Data_pert_ensemble(conditioni).v_r_1(pair_i,:) = obj.tv(rdt_idx,1);
                Data_pert_ensemble(conditioni).v_r_2(pair_i,:) = obj.tv(rdt_idx,2);
                                
                Data_pert_ensemble(conditioni).it_r(pair_i,:) = obj.it(rdt_idx,:);
                
                Data_pert_ensemble(conditioni).rdt_r(pair_i,:) = obj.rdt(rdt_idx,:);
                Data_pert_ensemble(conditioni).state_r(pair_i,:) = obj.state(rdt_idx);
                
                
%                 tmp = Data_pert(conditioni,pair_i).tp(:,2);
%                 plot(tmp);hold on;
%                 plot([1,1].*length(tmp)-500*5,[0, 1],'-k','linewidth',2.5); xlim([0 length(tmp)]); ylim([0.49 0.5]);hold off;
%                 pause();

                end
            end

            obj.Data_pert_ensemble = Data_pert_ensemble;
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

