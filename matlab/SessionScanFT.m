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
        FTrot_M1 = ...  % another position"raise-hand" 
                   ...  % global: x-right, y-front, z-up, FT_base x-rightback, y-rightfront, z-up
            [cosd(315) -sind(315)    0
             sind(315)  cosd(315)    0
                    0          0     1];
        FTrot_M2 = ... % another position "joystick", "raise-hand horizontal mount"
                    ... % global: x-right, y-front, z-up, FT_base x-rightfront, y-leftfront, z-up
            [cosd(45) -sind(45)    0
             sind(45)  cosd(45)    0
                    0          0     1];
        RDT
        FT
        elapse
        brtime % corresponding to black rock time
        force0

        cali_list = [4315 4312 4324 4327 4335 4338 4349 4354]; 
        % move this to calibration file later, for MVF which need original displacement
    end
    
    methods
        function obj = SessionScanFT(ss_num, fname_other)
            %FTSEPERATEDAT Construct an instance of this class
            %   Detailed explanation goes here
            fdir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
            %fdir = ['data/'];
            %fname = 'KingKongFT01865.csv';
            fname = sprintf('KingKongFT%05d.csv', ss_num);
            Data = readtable([fdir '/' fname]);
            if (exist('fname_other', 'var'))
                Data = readtable([fdir '/' fname_other]);
            end
            % deal with exception that end line is invalid 
            if (sum(isnan(Data{end,:})))
                Data = Data(1:end-1,:);
            end
            Data = dealRDTError(Data);
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
            obj.RDT = [Data.RDT]';           % read-time sequence
            obj.FT = [Data.FT]';             % FT sequence
            obj.elapse = [Data.elapse]';     % read time elapse, wrongly elapse... Change code! 

            if (sum(ss_num == obj.cali_list)) 
                obj.force_net = obj.force_origin;
            end
            obj = forceFTconvert(obj);
            %plotForce(obj)
            obj = obj.dealwithExceptions(ss_num);
            
            obj = intropTime(obj);
            
            ifplot = 0;
            if (ifplot)
                clf;
                axh(1) = subplot(2,1,1); plot(obj.RDT, obj.force_net(2,:));
                axh(2) = subplot(2,1,2); plot(obj.RDT, obj.torque_origin(2,:));
            end
            
            % this part hard code to deal with exception 
            
        end
        
        function obj = dealwithExceptions(obj, ss_num)
            % deal with exceptions 
            switch ss_num
             case  3775
                 obj.elapse(502620) = nanmean(obj.elapse(502619:502621));
             case  3778
                 obj.elapse(199100) = nanmean(obj.elapse(199099:199101));
             case  3803
                 obj.elapse(125225) = nanmean(obj.elapse(125224:125226));
                 obj.elapse(214050) = nanmean(obj.elapse(214049:214051));
                 obj.elapse(535351) = nanmean(obj.elapse(535350:535352));
                 obj.elapse(1354865) = nanmean(obj.elapse(1354864:1354866));
             case  3820
                 obj.elapse(83250) = nanmean(obj.elapse(83249:83251));
             case  3851
                 obj.elapse(656798) = nanmean(obj.elapse(656797:656799));
             case  3852
                 obj.elapse(1140790) = nanmean(obj.elapse(1140789:1140791));
                 obj.elapse(1352024) = nanmean(obj.elapse(1352023:1352025));
             case  3905
                 obj.elapse(613100) = nanmean(obj.elapse(613099:613101));
             case  3910
                 obj.elapse(368420) = nanmean(obj.elapse(368419:368421));
             case  3914
                 obj.elapse(191400) = nanmean(obj.elapse(191399:191401));
             case  3917
                 obj.elapse(602778) = nanmean(obj.elapse(602777:602779));
                 obj.elapse(680951) = nanmean(obj.elapse(680950:680952));
             case  3920
                 obj.elapse(92050) = nanmean(obj.elapse(92049:92051));
                 obj.elapse(746502)= nanmean(obj.elapse(746501:746503));
             case  3938
                 obj.elapse(650240) = nanmean(obj.elapse(650239:650241));
             case 4062
                 obj.elapse(90800) = nanmean(obj.elapse(90799:90801));
            end
        end
        
        function obj = intropTime(obj)
            % obj = intropTime(obj)
            % introp the time of the data from RDT. 
            % elapse was pfem time, and this time would be the same for
            % every batch (n==10) time. 

            time = obj.elapse; 
            rdt  = obj.RDT;
            if(~isempty(setdiff(unique(diff(rdt)),1)))
                %                 rdt1 = rdt(1):rdt(end);
                %                 time1= interp1(rdt(10:10:end),time(10:10:end), rdt1, 'linear', 'extrap');
                timediffidx = find(diff([time(1) time]) > 0.003);
                timediffdx = (timediffidx-1);            % select the last timepoint
                timediffidx = timediffdx(timediffdx>1);
                time1 = interp1(rdt(timediffidx),time(timediffidx), rdt, 'linear', 'extrap');
                ifplot = 1;
                if (ifplot)
                    clf; hold on;
                    plot(rdt, time, 'b*');
                    plot(rdt(timediffidx), time(timediffidx), 'mo');
                    plot(rdt, time1, 'r.');
                    legend('pfem time', 'difftime', 'reconstructed time' );
                end
            else
                time1= interp1(rdt(10:10:end),time(10:10:end), rdt, 'linear', 'extrap'); 
                                % has to be 10 because the last datapoint
                                % in the batch is the 'realtime'. The first
                                % should have time before that. 
                ifplot = 1;
                if (ifplot)
                    clf; hold on;
                    plot(rdt, time, 'b*');
                    plot(rdt, time1, 'r.');
                    legend('pfem time', 'difftime', 'reconstructed time' );
                end
            end
            

            
            obj.elapse = time1;
        end
        
        function plotForceOrigin(obj)
           figure();
           SAMPLE_R = 1; % plot from every 700 data points
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
            %obj.force = obj.FTrot_M * obj.force_origin;
%             obj.force = obj.FTrot_M * obj.force_net;
%             obj.force = obj.FTrot_M1 * obj.force_net; % raise-up configuration
            obj.force = obj.FTrot_M2 * obj.force_net; % joystick configuration

            %obj.force_net = obj.FTrot_M * obj.force_net;
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
        % plots here
        function plotSmoothedForce(obj)
            % smooth data
            force_smt = smoothdata(obj.force,2, 'movmean',20); % window 20
            %force_smt(1,:) = smoothdata(obj.force(1,:),'movmean',20);
            % plot
            figure(); hold on;
            plot(force_smt(1,:));
            plot(force_smt(2,:));
            plot(force_smt(3,:));
        end
        
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
        function plotForcexy_ss2271(obj)
            force0 = mean(obj.force_origin(:,1:50), 2);
            forceNet = obj.force_origin(:,:) - force0;
            force_xy = zeros(size(forceNet(1:2,:),2), 1);
            for time_i = 1:size(obj.force_origin, 2)
                force_xy(time_i) = norm(forceNet(1:2,time_i), 2);
            end
            force_z  = forceNet(3,:)';
            plot([force_xy, force_z]);
        end
    end
end

function Data = dealRDTError(Data)
    ifplot = 0;
    % sometime the FT has non-uniqe values, this is an error
    FT = [Data.FT]';
    [FTunq, idx_raw, idx_clean] = unique(FT);
    if (length(FTunq) < length(FT))
        fprintf('WARNING: FT time skew detected, abort %d data point!', ...
            length(FT) - length(FTunq));
    end
    idx_valid = idx_raw;
    Data = Data(idx_valid, :);
    
    % RDT non-unique is also an error
    RDT = [Data.RDT]';
    [RDTunq, idx_raw, idx_clean] = unique(RDT);
    if (length(RDTunq) < length(RDT))
        fprintf('WARNING: RDT time skew detected, abort %d data point!', ...
            length(RDT) - length(RDTunq));
    end
    idx_valid = idx_raw;
    Data = Data(idx_valid, :);
    
    % some times the elapse(time) has non-monotonic increasing, this is the
    % computer error.
    elapse = [Data.elapse]';
    elapse_diff = [0 diff(elapse) ];
    decrease_idx = elapse_diff<0;
    if (sum(elapse_diff<0)~=0)
        fprintf('WARNING: elapse time skew detected, abort %d data point!', ...
            sum(elapse_diff<0));
    end
    idx_valid = 0;    
    if (ifplot)
        axh(1) = subplot(2,1,1);
        plot([0, elapse_diff], '.');
        axh(2) = subplot(2,1,2);
        plot(elapse, '.');
        linkaxes(axh, 'x');
    end
    idx_valid = setdiff(1:size(elapse,2), find(decrease_idx)+1);
    Data = Data(idx_valid, :);
    % check again!
    elapse = [Data.elapse]';
    elapse_diff = [0 diff(elapse) ];
    decrease_idx = elapse_diff<0;
    if (sum(elapse_diff<0)~=0)
        fprintf('WARNING: elapse time skew detected, abort %d data point!', ...
            sum(elapse_diff<0));
    end
    if (ifplot)
        axh(1) = subplot(2,1,1);
        plot([0, elapse_diff], 'r.');
        axh(2) = subplot(2,1,2);
        plot(elapse, 'r.');
        linkaxes(axh, 'x');
    end
end
