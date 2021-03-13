classdef TrialScan
    %TRIALSCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        %%% trial variables
         %  general for all tasks
        tNo     % trial number
        bNo     % block number
        bgn     % bgn_idx in sessionScanObj 
        edn     % end_idx, avlid confliction
        bgn_t   % time for high_sample
        edn_t 
        outcome
        comboNo % get from intermediate data
        comboTT % combo of task targets
        states
        tarR    % target-rotation
        tarL    % target-length
        tarP    % target-position 
         % state indexes
        idx_bgn
        idx_prt
        idx_fcr
        idx_mov
        idx_hld
        idx_end
        idx_rst
        time_orn        % time origin
        time            % time after aligned 
        xyi             % xy index
        xyn             % xy name
         % specific ballistic-realease
        fTh             % force-threshold
        force           % 
        force_h         % from NetFT, 3-by-n
        force_t
        position
        position_h      % from WAM
        position_t
        velocity
        velocity_h
        velocity_t
        
        % predicted varialbes:
          % predicted linear variables just as Scott Did in thesis::pg44
        pred_x0
        pred_K
        pred_D
        pred_A
        pred_J
        pred_S %snap
        
        % perturbation related variables
        pert_t_bgn                      % perturbation time start
        pert_t_edn                      % perturbation time end
        pert_rdt_bgn                    % perturbation ReadTime bgn, RDT is SessionScanWam::rdt 
        pert_rdt_edn                    % perturbation ReadTime edn
        
    end
    
    methods
        %%% process
        function obj = TrialScan(sessionScanObj,trialNo)
            %TRIALSCAN Construct an instance of this class
            %   Construct a trial instance from the SessionScan class
            %   object, and trialNo
            %   asign all general variables, and specific ballistic-release
            %   variables from the trial
                    %%% task states, is it useful?
            ST_BGN = 1; 
            ST_PRT = 2;
            ST_FCR = 3; % force ramp
            ST_MOV = 4;
            ST_HLD = 5;
            ST_END = 6;
            ST_RST = 7;

            obj.tNo = trialNo;
            idx = find(sessionScanObj.Data.TrialNo == trialNo);
            obj.bgn = idx(1);
            obj.edn = idx(end);                             % end_idx, avlid confliction
            obj.bgn_t = sessionScanObj.Data.Time(obj.bgn);  % time for high_sample
            obj.edn_t = sessionScanObj.Data.Time(obj.edn);
            obj.outcome = unique(sessionScanObj.Data.OutcomeMasks.Success(obj.bgn:obj.edn));
            obj.comboNo = sessionScanObj.Data.ComboNo(obj.edn);
            obj.states  = unique(sessionScanObj.Data.TaskStateCodes.Values(obj.bgn:obj.edn));
            maskMov     = sessionScanObj.Data.TaskStateMasks.Move;
            maskTrial   = false(size(sessionScanObj.Data.TaskJudging.Target(5, :)));
            maskTrial(obj.bgn:obj.edn) = 1;   
            obj.tarR    = unique(sessionScanObj.Data.TaskJudging.Target(5, maskMov & maskTrial));          % target-rotation
            %obj.tarR    = unique(sessionScanObj.Data.TaskJudging.Target(5, obj.bgn:obj.edn));          % target-rotation
            obj.tarL    = unique(sessionScanObj.Data.TaskJudging.Target(6, maskMov & maskTrial));    % target-length
            %obj.tarL    = sort(unique(nonzeros(sessionScanObj.Data.TaskJudging.Target(6, obj.bgn:obj.edn))));    % target-length
            if isempty(obj.tarL)
                obj.tarL = nan;
            end
            [x, y]      = pol2cart(obj.tarR, obj.tarL);                       % TODO: consider the tarR convert to degree 
            obj.tarP    = [x, y];
             % state indexes, if multile there, select the first one
            tSCV = sessionScanObj.Data.TaskStateCodes.Values(obj.bgn:obj.edn); % taskStateCodeValues
            obj.idx_bgn = find_first_safe(tSCV, ST_BGN); % the first value
            obj.idx_prt = find_first_safe(tSCV, ST_PRT);
            obj.idx_fcr = find_first_safe(tSCV, ST_FCR); 
            obj.idx_mov = find_first_safe(tSCV, ST_MOV); 
            obj.idx_hld = find_first_safe(tSCV, ST_HLD); 
            obj.idx_end = find_first_safe(tSCV, ST_END); 
            obj.idx_rst = find_first_safe(tSCV, ST_RST); 
            obj.time_orn= sessionScanObj.time(obj.bgn:obj.edn);
            obj.time    = sessionScanObj.time(obj.bgn:obj.edn) - sessionScanObj.time(obj.bgn);       % time after aligned 
             % specific ballistic-realease
            obj.fTh = unique(nonzeros(sessionScanObj.Data.TaskJudging.Target(4, obj.bgn:obj.edn)));  % problematic here, if fTh==0; return 0 here         % force-threshold
            if isempty(obj.fTh)
                obj.fTh = nan;
            end
            obj.comboTT = getComboTT(obj,sessionScanObj);
            obj.force    = sessionScanObj.force(:,obj.bgn:obj.edn);
            obj.position = sessionScanObj.Data.Position.Actual(obj.bgn:obj.edn,:);
            if (~isempty(sessionScanObj.force_h))
                forceh_idx   = sessionScanObj.force_t >= obj.bgn_t & sessionScanObj.force_t <= obj.edn_t;
                obj.force_h  = sessionScanObj.force_h(:,forceh_idx);      % from NetFT
                obj.force_t  = sessionScanObj.force_t(forceh_idx) - sessionScanObj.time(obj.bgn);      % time aligned with trial 
            end
            if (~isempty(sessionScanObj.wamp_h))
                positionh_idx   = sessionScanObj.wam_t >= obj.bgn_t & sessionScanObj.wam_t <= obj.edn_t;
                obj.position_h  = sessionScanObj.wamp_h(positionh_idx,:)';     % from WAM
                obj.position_t  = sessionScanObj.wam_t(positionh_idx) - sessionScanObj.time(obj.bgn);     % time aligned with trial
            end
            if (~isempty(sessionScanObj.wamp_h))
                velocityh_idx   = sessionScanObj.wam_t >= obj.bgn_t & sessionScanObj.wam_t <= obj.edn_t;
                obj.velocity_h  = sessionScanObj.wamv_h(velocityh_idx,:)';     % from WAM
                obj.velocity_t  = sessionScanObj.wam_t(velocityh_idx) - sessionScanObj.time(obj.bgn);     % time aligned with trial
            end
            if isempty(setdiff(obj.tarR, [0,4])) %only y direction
                obj.xyi = 1;
            elseif isempty(setdiff(obj.tarR, [2, 6]))
                obj.xyi = 2;
            end
            xy_char = 'xy';
            obj.xyn = xy_char(obj.xyi);
            % other process
                if obj.tNo == 1
                    obj.outcome = 0;
                end
            % find perturbation time
            obj = findPerterbTime(obj, sessionScanObj);
        end
        function obj = findPerterbTime(obj, sessionScanObj)
            % find perturbation based on we only perturb on ForceRamp
             idx = obj.idx_fcr:obj.idx_mov;
             pert_t = obj.time_orn(idx);
            % find time  
            try
             obj.pert_t_bgn = pert_t(1);
             obj.pert_t_edn = pert_t(end);
            catch
                return
            end
              %if not bigger than 6s, abort.
             if length(pert_t) < 50*6 % 6s
                 display(['Not enough long in trial' num2str(obj.tNo)]);
                 return
             end
             % find RDT?
             pert_idx = sessionScanObj.time >= obj.pert_t_bgn & ...
                        sessionScanObj.time <= obj.pert_t_edn;
             pert_rdt = sessionScanObj.Data.Position.RDT(pert_idx);
             obj.pert_rdt_bgn = pert_rdt(1);
             obj.pert_rdt_edn = pert_rdt(end);
        end
        function obj = alignMOV(obj)
            %alignMOV align all trials at ST_MOV
            %   Just do linear shift, do NOT skew time
            time = obj.time;
            time_offset = time(obj.idx_mov);
            % if no ST_MOV, abort align
            if (~isempty(time_offset))
                obj.time = time - time_offset;
                
                if(~isempty(obj.position_t))
                    obj.position_t = obj.position_t - time_offset;
                end
                if(~isempty(obj.force_t))
                    obj.force_t = obj.force_t - time_offset;
                end
            end
        end
        
        function obj = cleanData(obj)
            % TODO: clean the trials only in specific part, avoid
            % un-related information.
        end
        
        function obj = predictImpedanceLinDev(obj)
            % obj = predictImpedanceLinDev(obj)
            % predict 1D dynamic parameters for simplicity
            % predict the Impedance using a linear derivative module, just
            % as Scott did in thesis pg27. Where:
            % F(t) = K(x0-x(t)) - Dx`(t) - Ax``(t) - Jx```(t)
            % solution:
            % Using matrix least squre regression:
            %  F = X*b
            %  b = inv(X'X)*X'*F
            %  where:
            %  b = [Kx0, -K, -D, -A, -J]';         % 5-by-1
            %  X = [1 x(1) x`(1) x``(1) x```(1) ...
            %       1 x(2) x`(2) x``(2) x```(2) ...
            %       ...
            %       1 x(n) x`(n) x``(n) x```(n)];   % n-by-5
            %  F = [F(1) F(2) ... F(n)]';           % n-by-1
            
            % check if the time aligned, if not, re-align
            t_bgn = 0;
            t_edn = 0.4;
            freq = 500;     % 500Hz, the same as WAM
            pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
            fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
            pos_t = obj.position_t(pos_t_idx);
            fce_t = obj.force_t(fce_t_idx);
            pos = obj.position_h(obj.xyi,pos_t_idx);
            fce = obj.force_h(obj.xyi,fce_t_idx);
            try
                ifsame = min(pos_t == fce_t); % one 0, all 0
            catch 
                ifsame = 0;
            end
            
            if (~ifsame)
                t_all = t_bgn:(1/freq):t_edn;
                % resample at time t
                pos_ = interp1(pos_t, pos, t_all); % .... process here
                fce_ = interp1(fce_t, fce, t_all); 
                % chances the last pos_ and fce_ is nan
                pos_ = pos_(1:end-1);
                fce_ = fce_(1:end-1);
                t_all= t_all(1:end-1);
                if sum(isnan(pos_) | isnan(fce_)) % if still have force
                    display('pos_ and fce_ have nan values, abort!');
                    return
                end
                
                if(0) % see the interp result
                    fh = figure();
                    subplot(2,1,1);
                    plot(pos_t, pos, 'o', t_all, pos_, ':.');
                    ylabel('position interp');
                    subplot(2,1,2);
                    plot(fce_t, fce, 'o', t_all, fce_, ':.');
                    ylabel('foece interp');
                end
                
            end
            x      = pos_ - pos(1);
            dx     = diff(x,1,2)/(1/freq);
            ddx    = diff(x,2,2)/((1/freq).^2);
            dddx   = diff(x,3,2)/((1/freq).^3);
            
            n = length(dddx);
            F = reshape(fce_(1:n),n,1);
            
            X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
                reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
            b = [];
            % calculate 
            b = (pinv(X'*X)*X'*F);
            % asign values
            obj.pred_K = -b(2);
            obj.pred_D = -b(3);
            obj.pred_A = -b(4);
            obj.pred_J = -b(5);
            obj.pred_x0 = b(1)/obj.pred_K;
            if(0)
                fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, J=%.3f, x0=%.3f', ...
                    obj.tNo, obj.pred_K, obj.pred_D, obj.pred_A, obj.pred_J, obj.pred_x0 );
            end
        end
        function obj = predictImpedanceLinDevS(obj)
            % obj = predictImpedanceLinDevS(obj)
            % predict 1D dynamic parameters for simplicity
            % predict the Impedance using a linear derivative module, just
            % expand the folulation Scott did in thesis pg27. Where:
            % F(t) = K(x0-x(t)) - Dx`(t) - Ax``(t) - Jx```(t) - Sx````(t)
            % solution:
            % Using matrix least squre regression:
            %  F = X*b
            %  b = inv(X'X)*X'*F
            %  where:
            %  b = [Kx0, -K, -D, -A, -J, -S]';         % 6-by-1
            %  X = [1 x(1) x`(1) x``(1) x```(1) ...
            %       1 x(2) x`(2) x``(2) x```(2) ...
            %       ...
            %       1 x(n) x`(n) x``(n) x```(n)];   % n-by-6
            %  F = [F(1) F(2) ... F(n)]';           % n-by-1
            
            % check if the time aligned, if not, re-align
            t_bgn = 0;
            t_edn = 0.4;
            freq = 500;     % 500Hz, the same as WAM
            pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
            fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
            pos_t = obj.position_t(pos_t_idx);
            fce_t = obj.force_t(fce_t_idx);
            pos = obj.position_h(obj.xyi,pos_t_idx);
            fce = obj.force_h(obj.xyi,fce_t_idx);
            try
                ifsame = min(pos_t == fce_t); % one 0, all 0
            catch 
                ifsame = 0;
            end
            
            if (~ifsame)
                t_all = t_bgn:(1/freq):t_edn;
                % resample at time t
                pos_ = interp1(pos_t, pos, t_all); % .... process here
                fce_ = interp1(fce_t, fce, t_all); 
                % chances the last pos_ and fce_ is nan
                pos_ = pos_(1:end-1);
                fce_ = fce_(1:end-1);
                t_all= t_all(1:end-1);
                if sum(isnan(pos_) | isnan(fce_)) % if still have force
                    display('pos_ and fce_ have nan values, abort!');
                    return
                end
                
                if(0) % see the interp result
                    fh = figure();
                    subplot(2,1,1);
                    plot(pos_t, pos, 'o', t_all, pos_, ':.');
                    ylabel('position interp');
                    subplot(2,1,2);
                    plot(fce_t, fce, 'o', t_all, fce_, ':.');
                    ylabel('foece interp');
                end
                
            end
            x      = pos_ - pos(1);
            dx     = diff(x,1,2)/(1/freq);
            ddx    = diff(x,2,2)/((1/freq).^2);
            dddx   = diff(x,3,2)/((1/freq).^3);
            ddddx  = diff(x,4,2)/((1/freq).^4);
            
            n = length(ddddx);
            F = reshape(fce_(1:n),n,1);
            X = [ones(n,1), ...
                reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
                reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1), ...
                reshape(ddddx(1:n),n,1)];
            b = [];
            % calculate 
            b = (pinv(X'*X)*X'*F);
            % asign values
            obj.pred_K = -b(2);
            obj.pred_D = -b(3);
            obj.pred_A = -b(4);
            obj.pred_J = -b(5);
            obj.pred_S = -b(6);
            obj.pred_x0 = b(1)/obj.pred_K;
            if(0)
                fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, J=%.3f, S=%.3f, x0=%.3f', ...
                    obj.tNo, obj.pred_K, obj.pred_D, obj.pred_A, obj.pred_J, obj.pred_S, obj.pred_x0 );
            end
        end
        
        function comboTT = getComboTT(obj,sessionScanObj)
        	% defined: targets = [obj.tarR, obj.tarL, obj.fTh];
            if length(unique(obj.tarR)) > 1
                display(['trial' num2str(obj.tNo) ' have more than 1 tarR']);
            end
            tarR = max(obj.tarR);
            tarL = obj.tarL;
            fTh  = obj.fTh;
            tarR_all = sessionScanObj.tarRs; 
            tarL_all = sessionScanObj.tarLs; 
            fTh_all  = sessionScanObj.fThs; 
            try % sometrials do not have tarR
                tarR_idx = find(tarR == tarR_all);
                tarL_idx = find(tarL == tarL_all);
                fTh_idx  = find(fTh  == fTh_all );
                comboTT = (tarR_idx-1) * length(tarL_all) * length(fTh_all) + ...
                          (tarL_idx-1) * length(fTh_all) + ...
                          fTh_idx;
            catch
                comboTT = nan;
            end
            if isempty(comboTT)
                comboTT = nan;
            end
        end
        
        %%% plot
        function axh = plotPredictedForceOnPosition(obj)
            % use regression terms to get the predicted force
            t_bgn = 0;
            t_edn = 0.4;
            freq = 500;     % 500Hz, the same as WAM
            pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
            fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
            pos_t = obj.position_t(pos_t_idx);
            fce_t = obj.force_t(fce_t_idx);
            pos = obj.position_h(obj.xyi,pos_t_idx);
            pos = pos - pos(1); % remove the offset
            fce = obj.force_h(obj.xyi,fce_t_idx);
            try
                ifsame = min(pos_t == fce_t); % one 0, all 0
            catch 
                ifsame = 0;
            end
            
            
            if (~ifsame)
                t_all = t_bgn:(1/freq):t_edn;
                % resample at time t
                pos_ = interp1(pos_t, pos, t_all); % .... process here
                fce_ = interp1(fce_t, fce, t_all); 
                % chances the last pos_ and fce_ is nan
                pos_ = pos_(1:end-1);
                fce_ = fce_(1:end-1);
                t_all= t_all(1:end-1);
                if sum(isnan(pos_) | isnan(fce_)) % if still have force
                    err('pos_ and fce_ have nan values, abort!');
                end
            end
            x      = pos_;
            dx     = diff(x,1,2)/(1/freq);
            ddx    = diff(x,2,2)/((1/freq).^2);
            dddx   = diff(x,3,2)/((1/freq).^3);
            
            n = length(dddx);
            X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
                reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
            b = [obj.pred_x0*obj.pred_K;
                -obj.pred_K;
                -obj.pred_D;
                -obj.pred_A;
                -obj.pred_J];
            F = X*b;
            t_all_ = t_all(1:length(F));
            axh = figure('Visible', 'on');
            hold on;
            plot(t_all, fce_, 'b', 'LineWidth', 3); 
            plot(t_all_, F, 'r--', 'LineWidth', 3); 
            xlabel('time (s)');
            ylabel('force (N)');
            legend('origin force', 'regressed force');
            title(['origin and regress force trial' num2str(obj.tNo)]);
            
        end
        function axh = plotPredictedForceOnPositionS(obj)
            % use regression terms to get the predicted force
            % with snap term
            t_bgn = 0;
            t_edn = 0.4;
            freq = 500;     % 500Hz, the same as WAM
            pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
            fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
            pos_t = obj.position_t(pos_t_idx);
            fce_t = obj.force_t(fce_t_idx);
            pos = obj.position_h(obj.xyi,pos_t_idx);
            pos = pos - pos(1); % remove the offset
            fce = obj.force_h(obj.xyi,fce_t_idx);
            try
                ifsame = min(pos_t == fce_t); % one 0, all 0
            catch 
                ifsame = 0;
            end
            
            
            if (~ifsame)
                t_all = t_bgn:(1/freq):t_edn;
                % resample at time t
                pos_ = interp1(pos_t, pos, t_all); % .... process here
                fce_ = interp1(fce_t, fce, t_all); 
                % chances the last pos_ and fce_ is nan
                pos_ = pos_(1:end-1);
                fce_ = fce_(1:end-1);
                t_all= t_all(1:end-1);
                if sum(isnan(pos_) | isnan(fce_)) % if still have force
                    err('pos_ and fce_ have nan values, abort!');
                end
            end
            x      = pos_;
            dx     = diff(x,1,2)/(1/freq);
            ddx    = diff(x,2,2)/((1/freq).^2);
            dddx   = diff(x,3,2)/((1/freq).^3);
            ddddx  = diff(x,4,2)/((1/freq).^4);
            n = length(ddddx);
            X = [ones(n,1), ...
                reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
                reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)...
                reshape(ddddx(1:n),n,1)];
            b = [obj.pred_x0*obj.pred_K;
                -obj.pred_K;
                -obj.pred_D;
                -obj.pred_A;
                -obj.pred_J;
                -obj.pred_S];
            F = X*b;
            t_all_ = t_all(1:length(F));
            axh = figure('Visible', 'on');
            hold on;
            plot(t_all, fce_, 'b', 'LineWidth', 3); 
            plot(t_all_, F, 'r--', 'LineWidth', 3); 
            xlabel('time (s)');
            ylabel('force (N)');
            legend('origin force', 'regressed force');
            title(['origin and regress force trial' num2str(obj.tNo)]);
            
        end
        function axh = plotPredictedForceOnPosition_lowsample(obj, sampleRate)
            % only plot, but not generate data here
            % predict and plot the Force and Position on a low sample rate
            % for example, FT runs at 500Hz whereas RTMA runs 50Hz if we
            % low sample at the rate of 500/50, it would be runs well 
            if (nargin<2)
                sampleRate = 10;
            end
            % %% part1: generate the low_sampled data. 
            t_bgn = 0;
            t_edn = 0.4;
            freq = 500;     % 500Hz, the same as WAM
            pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
            fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
            pos_t = obj.position_t(pos_t_idx);
            fce_t = obj.force_t(fce_t_idx);
            pos = obj.position_h(obj.xyi,pos_t_idx);
            fce = obj.force_h(obj.xyi,fce_t_idx);
            try
                ifsame = min(pos_t == fce_t); % one 0, all 0
            catch 
                ifsame = 0;
            end
            
            if (~ifsame)
                t_all = t_bgn:(1/freq):t_edn;
                % resample at time t
                pos_ = interp1(pos_t, pos, t_all); % .... process here
                fce_ = interp1(fce_t, fce, t_all); 
                % chances the last pos_ and fce_ is nan
                pos_ = pos_(1:end-1);
                fce_ = fce_(1:end-1);
                t_all= t_all(1:end-1);
                if sum(isnan(pos_) | isnan(fce_)) % if still have force
                    display('pos_ and fce_ have nan values, abort!');
                    return
                end
            end
            pos    = pos_ - pos(1);
            fce    = fce_;
            posR   = pos(1:sampleRate:end);
            fceR   = fce(1:sampleRate:end);
            resample_freq = 500/sampleRate;
            
            % %% part2: estimate through lines
            
            x      = posR;
            dx     = diff(x,1,2)/((1/resample_freq).^1);
            ddx    = diff(x,2,2)/((1/resample_freq).^2);
            dddx   = diff(x,3,2)/((1/resample_freq).^3);
            
            n = length(dddx);
            F = reshape(fceR(1:n),n,1);
            X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
                reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
            b = [];
            % calculate 
            b = (pinv(X'*X)*X'*F);
            % asign values
            pred_K = -b(2);
            pred_D = -b(3);
            pred_A = -b(4);
            pred_J = -b(5);
            pred_x0 = b(1)/obj.pred_K;
            if(0)
                fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, J=%.3f, x0=%.3f', ...
                    obj.tNo, obj.pred_K, obj.pred_D, obj.pred_A, obj.pred_J, obj.pred_x0 );
            end
            % %% calculating the predicted value
            X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
                reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
            b = [pred_x0*pred_K;
                -pred_K;
                -pred_D;
                -pred_A;
                -pred_J];
            F = X*b;
            t_all_ = t_all(1:length(F));
            axh = figure('Visible', 'on');
            hold on;
            % %% part3: figure the low_sampled data
            plot(t_all, fce_, ':.');   % origin
            t_rsp = t_bgn:(1/resample_freq):t_edn;             % t_resample
            t_prd = t_rsp;                                     % for prediction
            if length(t_rsp) > length(fceR)
                t_rsp = t_rsp(1:length(fceR));
            end
            plot(t_rsp, fceR, 'k', 'LineWidth', 3);                   % low-sample
            if length(t_prd) > length(F)
                t_prd = t_prd(1:length(F));
            end
            plot(t_prd, F, 'b', 'LineWidth', 3);
            xlabel('time (s)');
            ylabel('force (N)');
            legend('origin force', 'low sample force', 'pred force');
            title(['origin and regress force trial' num2str(obj.tNo)]);
        end
        function axh = plotPredictedForceODE(obj)
            % solve the force using differential equation
        end
    end
end

function rlt = find_first_safe(dat, mark)
    if (isempty(dat))
        rlt = [];
    end
    find_all = find(dat == mark);
    if (isempty(find_all))
        rlt = [];
    else
        rlt = find_all(1);
    end
    
end