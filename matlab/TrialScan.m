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
        outcomeo  % outcome from offline judge
%         comboNo % get from intermediate data
%         comboTT % combo of task targets
        states
        states_arr
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
        opt
        opth
        fTh             % force-threshold
        position_offset % steady position before release, as wam uses impedance control
        data
        ifpert
        pert_f
        wamKp
        wamBp
        % predicted varialbes:
          % predicted linear variables just as Scott Did in thesis::pg44
          % ==> Now use cleave updated Scott's method?
          % F(t) + Ms*x``(t) = Kx0 - Kx(t) - Dx`(t)
          % Only focus on mass on robot side, but can write as this form
          % F(t): censored force, need 0-aligned;
          % Ms:   subject side mass, get from Immediate release mass
          % x``(t): mass side acceleration;
          
        pred_x0
        pred_K
        pred_D
        pred_A
        pred_J
        pred_S %snap
        % perturbation related variables
        pert_iter
        perturbation_length = 1000
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
            ST_FHD = 4; % force hold
            ST_MOV = 5;
            ST_HLD = 6;
            ST_END = 7;
            ST_RST = 8;
            
            obj.tNo = trialNo;
            idx = find(sessionScanObj.Data.TrialNo == trialNo);
            obj.bgn = idx(1);
            obj.edn = idx(end);                             % end_idx, avlid confliction
            obj.bgn_t = sessionScanObj.time(obj.bgn);  % time for high_sample
            obj.edn_t = sessionScanObj.time(obj.edn);
            obj.outcome = unique(sessionScanObj.Data.OutcomeMasks.Success(obj.bgn:obj.edn));
%             obj.comboNo = sessionScanObj.Data.ComboNo(obj.edn);
            obj.states  = unique(sessionScanObj.Data.TaskStateCodes.Values(obj.bgn:obj.edn));
            obj.states_arr= sessionScanObj.Data.TaskStateCodes.Values(obj.bgn:obj.edn);
            
            %%%%%%%%%%%%%%%% block of getting data %%%%%%%%%%%%%%%%%%%%%%%%
            data_idx = sessionScanObj.data.t >= obj.bgn_t & sessionScanObj.data.t <= obj.edn_t;
            obj.data.t  = sessionScanObj.data.t(data_idx);
            obj.data.x  = sessionScanObj.data.x(:,data_idx);
            obj.data.v  = sessionScanObj.data.v(:,data_idx);
            obj.data.Fp = sessionScanObj.data.Fp(:,data_idx);
            obj.data.tq = sessionScanObj.data.tq(:,data_idx);
            obj.data.ts = sessionScanObj.data.ts(:,data_idx);
            obj.data.f  = sessionScanObj.data.f(:,data_idx);
            obj.data.ftq  = sessionScanObj.data.ftq(:,data_idx);
            if (~isempty(sessionScanObj.opt))   % if use optotrak
                data.optx=sessionScanObj.data.optx(:,data_idx);
                data.opty=sessionScanObj.data.opty(:,data_idx);
                data.optz=sessionScanObj.data.optz(:,data_idx);
                if_splineInterp = 1; 
                ifplot = 0;
                if (if_splineInterp) 
                    for marker_i = 1:size(data.optx,1)
                        % do the splineInterp here, for x, y, and z
                        val_idx = (~isnan(data.optx(marker_i,:))) ... 
                                 & (~isnan(data.opty(marker_i,:))) ...
                                 & (~isnan(data.optz(marker_i,:))); 
                        if (sum(val_idx) ~= 0)
                            intpx = spline(obj.data.t(val_idx), ...
                                data.optx(marker_i,val_idx),...
                                obj.data.t);
                            intpy = spline(obj.data.t(val_idx), ...
                                data.opty(marker_i,val_idx),...
                                obj.data.t);
                            intpz = spline(obj.data.t(val_idx), ...
                                data.optz(marker_i,val_idx),...
                                obj.data.t);
                        else
                            intpx = nan(size(obj.data.t));
                            intpy = nan(size(obj.data.t));
                            intpz = nan(size(obj.data.t));
                        end
                        if (ifplot)
                            clf; 
                            subplot(3,1,1); title('x'); hold on;
                            plot(obj.data.t, data.optx(marker_i,:), 'r.', 'MarkerSize', 10);
                            plot(obj.data.t, intpx, 'b.');
                            
                            subplot(3,1,2); title('y'); hold on;
                            plot(obj.data.t, data.opty(marker_i,:), 'r.', 'MarkerSize', 10);
                            plot(obj.data.t, intpy, 'b.');
                            
                            subplot(3,1,3); title('z'); hold on;
                            plot(obj.data.t, data.optz(marker_i,:), 'r.', 'MarkerSize', 10);
                            plot(obj.data.t, intpz, 'b.');
                        end
                        data.optx(marker_i,:) = intpx;
                        data.opty(marker_i,:) = intpy;
                        data.optz(marker_i,:) = intpz;
                        if(marker_i == 1)
                            obj.data.ox(1:3,:) = [  data.optx(marker_i,:); 
                                                    data.opty(marker_i,:); 
                                                    data.optz(marker_i,:); ];
                        end
                        
                    end
                end
            end 
            
            %%% optional: when emg data exists 
            if (isfield(sessionScanObj.data, 'emg'))
                obj.data.emg = sessionScanObj.data.emg(:,data_idx);
            end
            
            %%%%%%%%%%%%%%%% deal with perturbation, etc %%%%%%%%%%%%%%%%%%
            if isfield(sessionScanObj.Data.TaskJudging, 'ifpert')
                obj.ifpert  = ...
                double(unique(sessionScanObj.Data.TaskJudging.ifpert(obj.bgn:obj.edn)));
                if ~isempty(setdiff(obj.ifpert,0))
                    obj.ifpert = setdiff(obj.ifpert, 0);
                end
            else 
                % look for perturbation message in wam.cf
                obj.ifpert = ...
                    obj.findPerturbationinWAMcf(sessionScanObj);
            end
          
            obj.ifpert = obj.ifpert * obj.findPerturbationinWAMcf(sessionScanObj); % keep 0 when failed to pert
%             try 
%                 obj.ifpert  = ...
%                 double(unique(sessionScanObj.Data.TaskJudging.ifpert(obj.bgn:obj.edn)));
%             catch
%                 obj.ifpert = [];
%             end
            maskMov     = sessionScanObj.Data.TaskStateMasks.Move;
            maskBegin   = sessionScanObj.Data.TaskStateMasks.Begin;
            maskHold    = sessionScanObj.Data.TaskStateMasks.Hold;
            maskTrial   = false(size(sessionScanObj.Data.TaskJudging.Target(5, :)));
            maskTrial(obj.bgn:obj.edn) = 1;   
            obj.tarR    = unique(sessionScanObj.Data.TaskJudging.Target(5, maskMov & maskTrial));          % target-rotation
            obj.tarL    = unique(sessionScanObj.Data.TaskJudging.Target(6, maskMov & maskTrial));    % target-length % some old trials
            if isempty(obj.tarL)
                obj.tarL = nan;
            end
            [x, y]      = pol2cart(obj.tarR, obj.tarL);                       % TODO: consider the tarR convert to degree 
            obj.tarP    = [x, y];
             % state indexes, if multile there, select the first one
            ts = sessionScanObj.Data.TaskStateCodes.Values(obj.bgn:obj.edn); % taskStateCodeValues
            obj.idx_bgn = find_first_safe(ts, ST_BGN); % the first value
            obj.idx_prt = find_first_safe(ts, ST_PRT);
            obj.idx_fcr = find_first_safe(ts, ST_FCR); 
            obj.idx_mov = find_first_safe(ts, ST_MOV); 
            obj.idx_hld = find_first_safe(ts, ST_HLD); 
            obj.idx_end = find_first_safe(ts, ST_END); 
            obj.idx_rst = find_first_safe(ts, ST_RST); 
            obj.time_orn= sessionScanObj.time(obj.bgn:obj.edn);
            obj.time    = sessionScanObj.time(obj.bgn:obj.edn) - sessionScanObj.time(obj.bgn);       % time after aligned 
             % specific ballistic-realease
            obj.fTh = unique(nonzeros(sessionScanObj.Data.TaskJudging.Target(4, obj.bgn:obj.edn)));  % problematic here, if fTh==0; return 0 here         % force-threshold
            if isempty(obj.fTh)
                obj.fTh = nan;
            end
            
            if length(sessionScanObj.pertCond.wamKp) == 1
                obj.wamKp    = sessionScanObj.pertCond.wamKp;
            else
                obj.wamKp    = unique(sessionScanObj.Data.TaskJudging.wamKp(maskTrial & maskHold));
            end
            
            if length(sessionScanObj.pertCond.wamBp) == 1
                obj.wamBp    = sessionScanObj.pertCond.wamBp;
            else
                 obj.wamBp    = unique(sessionScanObj.Data.TaskJudging.wamBp(maskTrial & maskHold));
            end

            try
                obj.pert_f   = unique(sessionScanObj.Data.TaskJudging.pertdf_mag(maskTrial & maskHold));
            catch
                 obj.pert_f  = [];
            end
            
%             obj.comboTT = getComboTT(obj,sessionScanObj);
            
            
            if isempty(setdiff(obj.tarR, [0,4])) %only y direction
                obj.xyi = 1;
            elseif isempty(setdiff(obj.tarR, [2, 6]))
                obj.xyi = 2;
            end
            
            ifplot = 0;
            if (ifplot)
                axh(1) = subplot(2,1,1);
                plot(obj.data.t, obj.data.x);
                axh(2) = subplot(2,1,2);
                plot(obj.data.t, obj.data.f);
                linkaxes(axh, 'x');
            end
            
%             if (~isempty(sessionScanObj.opt))
%                 opt_idx   = sessionScanObj.opt.data.t >= obj.bgn_t & sessionScanObj.opt.data.t <= obj.edn_t;
%                 opt_idxh  = sessionScanObj.opt.datah.t >= obj.bgn_t & sessionScanObj.opt.datah.t <= obj.edn_t;
%                 obj.opt.data.t = sessionScanObj.opt.data.t(opt_idx) - sessionScanObj.time(obj.bgn);
%                 obj.opt.data.rdt = sessionScanObj.opt.data.rdt(opt_idx);
%                 obj.opt.data.x = sessionScanObj.opt.data.x(opt_idx);
%                 obj.opt.data.y = sessionScanObj.opt.data.y(opt_idx);
%                 obj.opt.data.z = sessionScanObj.opt.data.z(opt_idx);
%                 
%                 obj.opt.datah.t = sessionScanObj.opt.datah.t(opt_idxh) - sessionScanObj.time(obj.bgn);
%                 obj.opt.datah.rdt = sessionScanObj.opt.datah.rdt(opt_idxh);
%                 obj.opt.datah.x = sessionScanObj.opt.datah.x(opt_idxh);
%                 obj.opt.datah.y = sessionScanObj.opt.datah.y(opt_idxh);
%                 obj.opt.datah.z = sessionScanObj.opt.datah.z(opt_idxh);
%             end
            
            xy_char = 'xy';
            obj.xyn = xy_char(obj.xyi);
            % other process
                if obj.tNo == 1
                    obj.outcome = 0;
                end
            % find (double check) perturbation time
%             if (obj.ifpert==1 || obj.ifpert==0)
%                 obj = findStepPerterbTime(obj, sessionScanObj);
%             else
%                 obj = findStocPerterbTime(obj, sessionScanObj);
%             end
            %obj = findStepx0PerterbTime(obj, sessionScanObj);
            %obj.ifpert = obj.ifpert * obj.findPerturbationinWAMcf(sessionScanObj); % check it again for sanity!
        end
        function ifpert = findPerturbationinWAMcf(obj, sessionScanObj)
            % find if being perturbed via looking at wam.cf data 
            % applicable for some sessions without ifpert variable
            wam_t = obj.data.t;
            wam_cf = obj.data.Fp;
            wam_ts = obj.data.ts;
            position_h = obj.data.x;
            ifpert = ~sum(sum(abs(wam_cf)))==0;
            
            if (~isempty(wam_ts) && ~isempty(wam_cf))
                pertsig = setdiff([(unique(wam_ts(wam_cf(1,:)~=0)))... % pert at x
                                   (unique(wam_ts(wam_cf(2,:)~=0)))... % pert at y
                                   (unique(wam_ts(wam_cf(3,:)~=0)))]... % pert at z
                                   , 3); 
                if (isempty(pertsig))
                    obj.ifpert = 0;
                    ifpert = 0;
                end
                ifplot = 0;
                if(ifplot)
                    axhl(1) = subplot(3,1,1);
                    plot(wam_t, position_h(2,:));
                    axhl(2) = subplot(3,1,2);
                    plot(wam_t, wam_cf(2,:));
                    axhl(3) = subplot(3,1,3); hold on;
                    plot(wam_t, wam_ts, 'r');
                    linkaxes(axhl, 'x');
                end

            end
            
        end
%         function pt = getPerturbationPattern(obj)
%             % get the perturbation pattern
%             % 1. Stochastic force
%             % 2. step force
%             pert_fce = unique(obj.pertfce_h); 
%             pert_fce(pert_fce<0.1) = 0; % basiclly will not exert force level < 0.1; 
%             if ((isempty(obj.pert_dx0) || isempty(setdiff(obj.pert_dx0, 0))) && ~isempty(setdiff(pert_fce,0)) )
%                 pt = 2;
%                 return 
%             end
%             % 3. step x0
%             pt = 0;
%             if (obj.ifpert && isempty(setdiff(pert_fce ,0)))
%                 pt = 3;
%                 return
%             end
%         end
%         function obj = findStocPerterbTime(obj, sessionScanObj)
%             % find perturbation based on we only perturb on ForceRamp
%             % This is only for the elongated stochastic perturbation.
%              %idx = obj.idx_fcr:obj.idx_mov;
%              idx = obj.idx_fcr:obj.idx_end;
%              pert_t = obj.time_orn(idx);
%             % find time  
%             try
%              obj.pert_t_bgn = pert_t(1);
%              obj.pert_t_edn = pert_t(end);
%             catch
%                 return
%             end
%               %if not bigger than 6s, abort.
%              if length(pert_t) < 50*6 % 6s
%                  display(['Not enough long in trial' num2str(obj.tNo)]);
%                  return
%              end
%              % find RDT?
%              pert_idx = sessionScanObj.time >= obj.pert_t_bgn & ...
%                         sessionScanObj.time <= obj.pert_t_edn;
%              pert_rdt = sessionScanObj.Data.Position.RDT(pert_idx);
%              obj.pert_rdt_bgn = pert_rdt(1);    % force increasing time (pert start)
%              obj.pert_rdt_edn = pert_rdt(end);  % release time (pert finished)
%         end
%         function obj = findStepPerterbTime(obj, sessionScanObj)
%             % find perturbation based on we only perturb on ForceRamp
%             % This is only for the STEP PERTURBATION.
%             wam_pert_signal = obj.data.Fp(2,:);
%             wam_pert_init_idx = find(wam_pert_signal == 0 & [diff(wam_pert_signal) 0]~=0); 
%             wam_pert_edn_idx  = find(wam_pert_signal ~= 0 & [diff(wam_pert_signal) 0]~=0); 
%             if obj.ifpert == 0 || isempty(wam_pert_init_idx)
% %                 obj.pert_rdt_bgn = [];
% %                 obj.pert_rdt_edn = [];
%                 % The comming line, only works at the real experiment,
%                 % if were conducting spring test, can comment it out.
%                 %obj.ifpert = 0; % re-write (some trial should be perturbed, but did not wait until it)
%                 return
%             end
%             fin_pert_init_idx = wam_pert_init_idx(end); % ONLY use the final perturb 
%             % dangerous in the task analysis here, as subject may change
%             % strategy to fulfill the perturbation requirements. 
%             fin_pert_edn_idx  = wam_pert_edn_idx(end);
%             if fin_pert_edn_idx < fin_pert_init_idx
%                 display(['trial' obj.tNo 'has unfinished perturbation']);
%                 obj.pert_rdt_bgn = [];
%                 obj.pert_rdt_edn = [];
%             end
%             idxh = fin_pert_init_idx:fin_pert_edn_idx;
%             pert_t = obj.data.t(idxh);
%             % find time  
%             try
%              obj.pert_t_bgn = pert_t(1);
%              obj.pert_t_edn = pert_t(end);
%             catch
%                 return
%             end
% %              % find RDT
% %              obj.pert_rdt_bgn = obj.wamrdt(fin_pert_init_idx);    % force increasing time (pert start)
% %              obj.pert_rdt_edn = obj.wamrdt(fin_pert_edn_idx);  % release time (pert finished)
% %              
% %              % position offset
% %              pos_bef_pert = obj.position_h(2, fin_pert_init_idx-100:fin_pert_init_idx-1); 
% %              obj.position_offset = mean(pos_bef_pert);
%         end
%         function obj = findStepx0PerterbTime(obj, sessionScanObj)
%             % find perturbation based on we only perturb on x0
%             % As we did not record the x0 (at code version Sep7th, 2021),
%             % we use the torque (jt) as a marker to recognize the
%             % perturbation
%             % This is only for the STEP PERTURBATION.
%             wam_pert_signal = obj.pert_iter;
%             % find index
%             wam_pert_signald= [0 diff(wam_pert_signal)'];
%             wam_pert_signaldd= [0 diff(wam_pert_signald)];
%             [idx] = find(wam_pert_signaldd == 1); % hope normal activity will not cross 5*std
%             if (length(idx)<1) % did not detected the positive and negative edge
%                 %                 obj.pert_rdt_bgn = [];
%                 obj.pert_rdt_edn = [];
%                 % The comming line, only works at the real experiment,
%                 % if were conducting spring test, can comment it out.
%                 %obj.ifpert = 0; % re-write (some trial should be perturbed, but did not wait until it)
%                 return
%             end
%             wam_pert_init_idx = idx; 
%             wam_pert_edn_idx  = wam_pert_init_idx + 1000 - 1; 
%             fin_pert_init_idx = wam_pert_init_idx(end); % ONLY use the final perturb 
%             % dangerous in the task analysis here, as subject may change
%             % strategy to fulfill the perturbation requirements. 
%             fin_pert_edn_idx  = wam_pert_edn_idx(end);
%             if fin_pert_edn_idx < fin_pert_init_idx
%                 display(['trial' obj.tNo 'has unfinished perturbation']);
%                 obj.pert_rdt_bgn = [];
%                 obj.pert_rdt_edn = [];
%             end
%             idxh = fin_pert_init_idx:fin_pert_edn_idx;
%             pert_t = obj.position_t(idxh);
%             % find time  
%             try
%              obj.pert_t_bgn = pert_t(1);
%              obj.pert_t_edn = pert_t(end);
%             catch
%                 return
%             end
%              % find RDT
%              obj.pert_rdt_bgn = obj.wamrdt(fin_pert_init_idx);    % force increasing time (pert start)
%              obj.pert_rdt_edn = obj.wamrdt(fin_pert_edn_idx);  % release time (pert finished)
%              
%              % position offset
%              pos_bef_pert = obj.position_h(2, fin_pert_init_idx-100:fin_pert_init_idx-1); 
%              obj.position_offset = mean(pos_bef_pert);
%         end
        function obj = alignMOV(obj)
            %alignMOV align all trials at ST_MOV
            %   Just do linear shift, do NOT skew time
            time = obj.time;
            idx = find(obj.data.ts == 5);
            if isempty(idx)
                time_offset = obj.data.t(end); % not at all moved
            else
                time_offset = obj.data.t(idx(1));
            end
            % if no ST_MOV, abort align
            if (~isempty(time_offset))
                obj.time = time - time_offset;
                
                if(~isempty(obj.data.t))
                    obj.data.t_shift = obj.data.t - time_offset;
                end
                if(~isempty(obj.opt))
                    obj.opt.data.t = obj.opt.data.t - time_offset;
                    obj.opt.datah.t = obj.opt.datah.t - time_offset;
                end
            end
        end
        function obj = alignPertInit(obj, sessionScanObj)
            %alignPertInit align all trials at the perturbation start
            % specific for the STEP PERTURBATION
            % Just do linear shift, do not skew time
            time = obj.time;
            % should use the # to find the exact perturb time
            if obj.ifpert == 0 || isempty(obj.pert_rdt_bgn) % there is no perturbation in the current trial
                return
            end 
            idx_t = find(obj.pert_rdt_bgn == obj.wamrdt);
            time_offset = obj.position_t(idx_t); 
            
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
%         function obj = cleanData(obj)
%             % TODO: clean the trials only in specific part, avoid
%             % un-related information.
%         end
        
%         function obj = predictImpedanceLinDev(obj)
%             % obj = predictImpedanceLinDev(obj)
%             % predict 1D dynamic parameters for simplicity
%             % predict the Impedance using a linear derivative module, just
%             % as Scott did in thesis pg27. Where:
%             % F(t) = K(x0-x(t)) - Dx`(t) - Ax``(t) - Jx```(t)
%             % solution:
%             % Using matrix least squre regression:
%             %  F = X*b
%             %  b = inv(X'X)*X'*F
%             %  where:
%             %  b = [Kx0, -K, -D, -A, -J]';         % 5-by-1
%             %  X = [1 x(1) x`(1) x``(1) x```(1) ...
%             %       1 x(2) x`(2) x``(2) x```(2) ...
%             %       ...
%             %       1 x(n) x`(n) x``(n) x```(n)];   % n-by-5
%             %  F = [F(1) F(2) ... F(n)]';           % n-by-1
%             
%             % check if the time aligned, if not, re-align
%             t_bgn = 0;
%             t_edn = 1.4; % 0.4
%             freq = 500;     % 500Hz, the same as WAM
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             fce_t = obj.force_t(fce_t_idx);
%             try % for the normal trials
%                 pos = obj.position_h(obj.xyi,pos_t_idx);
%                 fce = [sind(-45), cosd(45)] * obj.force_h(1:2,fce_t_idx);
%             catch % for the simulated trials
%                 pos = obj.position_h(:,pos_t_idx);
%                 fce = obj.force_h(:,fce_t_idx);
%             end
%             try
%                 ifsame = min(pos_t == fce_t); % one 0, all 0
%             catch 
%                 ifsame = 0;
%             end
%             
%             if (~ifsame)
%                 t_all = t_bgn:(1/freq):t_edn;
%                 % resample at time t
%                 pos_ = interp1(pos_t, pos, t_all); % .... process here
%                 fce_ = interp1(fce_t, fce, t_all); 
%                 % chances the last pos_ and fce_ is nan
%                 pos_ = pos_(1:end-1);
%                 fce_ = fce_(1:end-1);
%                 t_all= t_all(1:end-1);
%                 if sum(isnan(pos_) | isnan(fce_)) % if still have force
%                     display('pos_ and fce_ have nan values, abort!');
%                     return
%                 end
%                 
%                 if(0) % see the interp result
%                     fh = figure();
%                     subplot(2,1,1);
%                     plot(pos_t, pos, 'o', t_all, pos_, ':.');
%                     ylabel('position interp');
%                     subplot(2,1,2);
%                     plot(fce_t, fce, 'o', t_all, fce_, ':.');
%                     ylabel('foece interp');
%                 end
%             else 
%                 pos_ = pos;
%                 fce_ = fce;
%             end
%             x      = pos_ - pos(1);                 % unit: m
%             dx     = diff(x,1,2)/(1/freq);          % unit: m/s
%             ddx    = diff(x,2,2)/((1/freq).^2);     % unit: m/s^2
%             dddx   = diff(x,3,2)/((1/freq).^3);     % unit: m/s^3
%             islowpass = 1;
%             if (islowpass)
%                 x_lowpass = lowpass(x, 4, 500);
%                 x_lowpass = x_lowpass(1:500);
%                 dx   = diff(x_lowpass,1,2)/(1/freq);
%                 ddx  = diff(x_lowpass,1,2)/((1/freq).^2);
%                 dddx = diff(x_lowpass,1,2)/((1/freq).^3);
%                 f_lowpass = lowpass(fce_, 4, 500);
%                 fce_ = f_lowpass(1:500);
%                 fce_0 = mean(fce_(end-99:end));
%                 fce_ = fce_ - fce_0;
%             end
%             n = length(dddx);
%             F = reshape(fce_(1:n),n,1);
%             
%             X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
%                 reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
%             b = [];
%             % calculate 
%             b = (inv(X'*X)*X'*F); % looks un-safe, is other ways to do least square?
%             % asign values
%             obj.pred_K = -b(2);
%             obj.pred_D = -b(3);
%             obj.pred_A = -b(4);
%             obj.pred_J = -b(5);
%             obj.pred_x0 = b(1)/obj.pred_K;
%             if(1)
%                 fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, J=%.3f, x0=%.3f', ...
%                     obj.tNo, obj.pred_K, obj.pred_D, obj.pred_A, obj.pred_J, obj.pred_x0 );
%             end
%         end
%         function obj = predictImpedanceLinDev2ndOrder(obj)
%             % obj = predictImpedanceLinDev(obj)
%             % predict 1D dynamic parameters for simplicity
%             % predict the Impedance using a linear derivative module, just
%             % as Scott did in thesis pg27. Where:
%             % F(t) = K(x0-x(t)) - Dx`(t) - Ax``(t)
%             % solution:
%             % Using matrix least squre regression:
%             %  F = X*b
%             %  b = inv(X'X)*X'*F
%             %  where:
%             %  b = [Kx0, -K, -D, -A, -J]';         % 5-by-1
%             %  X = [1 x(1) x`(1) x``(1) x```(1) ...
%             %       1 x(2) x`(2) x``(2) x```(2) ...
%             %       ...
%             %       1 x(n) x`(n) x``(n) x```(n)];   % n-by-4
%             %  F = [F(1) F(2) ... F(n)]';           % n-by-1
%             
%             % check if the time aligned, if not, re-align
%             t_bgn = 0;
%             t_edn = 1.4; % 0.4
%             freq = 500;     % 500Hz, the same as WAM
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             fce_t = obj.force_t(fce_t_idx);
%             try % for the normal trials
%                 pos = obj.position_h(obj.xyi,pos_t_idx);
%                 fce = [sind(-45), cosd(45)] * obj.force_h(1:2,fce_t_idx);
%             catch % for the simulated trials
%                 pos = obj.position_h(:,pos_t_idx);
%                 fce = obj.force_h(:,fce_t_idx);
%             end
%             try
%                 ifsame = min(pos_t == fce_t); % one 0, all 0
%             catch 
%                 ifsame = 0;
%             end
%             
%             if (~ifsame)
%                 t_all = t_bgn:(1/freq):t_edn;
%                 % resample at time t
%                 pos_ = interp1(pos_t, pos, t_all); % .... process here
%                 fce_ = interp1(fce_t, fce, t_all); 
%                 % chances the last pos_ and fce_ is nan
%                 pos_ = pos_(1:end-1);
%                 fce_ = fce_(1:end-1);
%                 t_all= t_all(1:end-1);
%                 if sum(isnan(pos_) | isnan(fce_)) % if still have force
%                     display('pos_ and fce_ have nan values, abort!');
%                     return
%                 end
%                 
%                 if(0) % see the interp result
%                     fh = figure();
%                     subplot(2,1,1);
%                     plot(pos_t, pos, 'o', t_all, pos_, ':.');
%                     ylabel('position interp');
%                     subplot(2,1,2);
%                     plot(fce_t, fce, 'o', t_all, fce_, ':.');
%                     ylabel('foece interp');
%                 end
%             else 
%                 pos_ = pos;
%                 fce_ = fce;
%             end
%             x      = pos_ - pos(1);                 % unit: m
%             dx     = diff(x,1,2)/(1/freq);          % unit: m/s
%             ddx    = diff(x,2,2)/((1/freq).^2);     % unit: m/s^2
%             islowpass = 1;
%             if (islowpass)
%                 %x_lowpass = lowpass(x, 4, 500);
%                 x_lowpass = x; 
%                 x_lowpass = x_lowpass(1:500);
%                 dx   = diff(x_lowpass,1,2)/(1/freq);
%                 ddx  = diff(dx,1,2)/((1/freq));
%                 f_lowpass = lowpass(fce_, 4, 500);
%                 fce_ = f_lowpass(1:500);
%                 fce_0 = mean(fce_(end-99:end));
%                 fce_ = fce_ - fce_0;
%             end
%             n = length(ddx);
%             F = reshape(fce_(1:n),n,1);
%             
%             X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
%                 reshape(ddx(1:n),n,1)];
%             b = [];
%             % calculate 
%             b = (inv(X'*X)*X'*F); % looks un-safe, is other ways to do least square?
%             % asign values
%             obj.pred_K = -b(2);
%             obj.pred_D = -b(3);
%             obj.pred_A = -b(4);
%             obj.pred_x0 = b(1)/obj.pred_K;
%             if(1)
%                 fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, x0=%.3f', ...
%                     obj.tNo, obj.pred_K, obj.pred_D, obj.pred_A, obj.pred_x0 );
%             end
%         end
%         function obj = predictImpedanceLinDev2ndOrderFixM(obj)
%             % obj = predictImpedanceLinDev(obj)
%             % predict 1D dynamic parameters for simplicity
%             % predict the Impedance using a linear derivative module, just
%             % as Scott did in thesis pg27. Where:
%             % F(t) = K(x0-x(t)) - Dx`(t) - Ax``(t)
%             % Make a alteration so that: 
%             % F(t)/M = K/M(x0-x(t)) - D/Mx`(t) - x``(t)
%             % solution:
%             % Using matrix least squre regression:
%             %  F/m + x``(t) = X*b
%             %  b = inv(X'X)*X'*F
%             %  where:
%             %  b = [Kx0, -K, -D]';          % 3-by-1
%             %  X = [1 x(1) x`(1)  ...
%             %       1 x(2) x`(2)  ...
%             %       ...
%             %       1 x(n) x`(n)];          % n-by-4
%             %  F = [F(1) F(2) ... F(n)]';   % n-by-1
%             
%             % check if the time aligned, if not, re-align
%             t_bgn = -1;
%             t_edn = 2.4; % 0.4
%             freq = 500;     % 500Hz, the same as WAM
%             M = 3; % assume 3 kg
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             fce_t = obj.force_t(fce_t_idx);
%             try % for the normal trials
%                 pos = obj.position_h(obj.xyi,pos_t_idx);
%                 fce = [sind(-45), cosd(45)] * obj.force_h(1:2,fce_t_idx);
%             catch % for the simulated trials
%                 pos = obj.position_h(:,pos_t_idx);
%                 fce = obj.force_h(:,fce_t_idx);
%             end
%             try
%                 ifsame = min(pos_t == fce_t); % one 0, all 0
%             catch 
%                 ifsame = 0;
%             end
%             
%             if (~ifsame)
%                 t_all = t_bgn:(1/freq):t_edn;
%                 % resample at time t
%                 pos_ = interp1(pos_t, pos, t_all); % .... process here
%                 fce_ = interp1(fce_t, fce, t_all); 
%                 % chances the last pos_ and fce_ is nan
%                 pos_ = pos_(1:end-1);
%                 fce_ = fce_(1:end-1);
%                 t_all= t_all(1:end-1);
%                 if sum(isnan(pos_) | isnan(fce_)) % if still have force
%                     display('pos_ and fce_ have nan values, abort!');
%                     return
%                 end
%                 
%                 if(0) % see the interp result
%                     fh = figure();
%                     subplot(2,1,1);
%                     plot(pos_t, pos, 'o', t_all, pos_, ':.');
%                     ylabel('position interp');
%                     subplot(2,1,2);
%                     plot(fce_t, fce, 'o', t_all, fce_, ':.');
%                     ylabel('foece interp');
%                 end
%             else 
%                 pos_ = pos;
%                 fce_ = fce;
%             end
%             x      = pos_ - pos(1);                 % unit: m
%             dx     = diff(x,1,2)/(1/freq);          % unit: m/s
%             ddx    = diff(x,2,2)/((1/freq).^2);     % unit: m/s^2 
%             % try forward and backward diffrentiation ==> seems useless.
%             x       = pos_ - pos(1); 
%             dx_f    = [diff(x,1,2)/(1/freq) 0];         % unit: m/s
%             dx_b    = [0 diff(x,1,2)/(1/freq)];  % backward differntiation
%             dx_     = (dx_f + dx_b)/2;   
%             ddx_f   = [diff(dx_,1,2)/(1/freq) 0];
%             ddx_b   = [0 diff(dx_,1,2)/(1/freq)];
%             ddx_    = (ddx_f + ddx_b)/2;   
%             islowpass = 0;
%             if (islowpass)
%                 x_lowpass = lowpass(x, 4, 500);
%                 x_lowpass = x_lowpass(1:500);
%                 f_lowpass = lowpass(fce_, 4, 500);
%             else
%                 x_lowpass = x; %x_lowpass = x(1:500);
%                 f_lowpass = fce_;
%             end
%             ddx_length = length(ddx);
%             %fce_ = f_lowpass(1:ddx_length);
%             %fce_0 = mean(fce_(end-99:end));
%             %fce_ = fce_ - fce_0;
%                 
%             
%             M = fce_(1:end-2) ./ ddx;
%             M = mean(M(~isnan(M)));
%             f_m_a = ddx;
%             n = length(ddx(501:end));
%             F = reshape(f_m_a(end-n+1:end),n,1);
%             
%             X = [ones(n,1), reshape(x(end-n+1:end),n,1), reshape(dx(end-n+1:end),n,1)];
%             %X = [ones(n,1), reshape(x(end-n+1-2:end-2),n,1), reshape(dx(end-n+1-1:end-1),n,1)];
%             %X = [ones(n,1), reshape(x(end-n+1:end),n,1), reshape(dx(end-n+1-1:end-1),n,1)];
%             b = [];
%             % calculate 
%             b = (inv(X'*X)*X'*F); % looks un-safe, is other ways to do least square?
%             b = ((X'*X)\(X'*F));
%             % asign values
%             obj.pred_K = -b(2) * M;
%             obj.pred_D = -b(3) * M;
%             obj.pred_x0 = b(1)/obj.pred_K * M;
%             obj.pred_A = M;
%             if(1)
%                 fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, x0=%.3f', ...
%                     obj.tNo, obj.pred_K, obj.pred_D, M, obj.pred_x0 );
%             end
%             
%             % get back to original points and get r2
%             F_pred = X*b;
%             mean(F - F_pred)
%         end
%         function obj = predictImpedanceLinDevS(obj)
%             % obj = predictImpedanceLinDevS(obj)
%             % predict 1D dynamic parameters for simplicity
%             % predict the Impedance using a linear derivative module, just
%             % expand the folulation Scott did in thesis pg27. Where:
%             % F(t) = K(x0-x(t)) - Dx`(t) - Ax``(t) - Jx```(t) - Sx````(t)
%             % solution:
%             % Using matrix least squre regression:
%             %  F = X*b
%             %  b = inv(X'X)*X'*F
%             %  where:
%             %  b = [Kx0, -K, -D, -A, -J, -S]';         % 6-by-1
%             %  X = [1 x(1) x`(1) x``(1) x```(1) ...
%             %       1 x(2) x`(2) x``(2) x```(2) ...
%             %       ...
%             %       1 x(n) x`(n) x``(n) x```(n)];   % n-by-6
%             %  F = [F(1) F(2) ... F(n)]';           % n-by-1
%             
%             % check if the time aligned, if not, re-align
%             t_bgn = 0;
%             t_edn = 0.4;
%             freq = 500;     % 500Hz, the same as WAM
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             fce_t = obj.force_t(fce_t_idx);
%             pos = obj.position_h(obj.xyi,pos_t_idx);
%             fce = obj.force_h(obj.xyi,fce_t_idx);
%             try
%                 ifsame = min(pos_t == fce_t); % one 0, all 0
%             catch 
%                 ifsame = 0;
%             end
%             
%             if (~ifsame)
%                 t_all = t_bgn:(1/freq):t_edn;
%                 % resample at time t
%                 pos_ = interp1(pos_t, pos, t_all); % .... process here
%                 fce_ = interp1(fce_t, fce, t_all); 
%                 % chances the last pos_ and fce_ is nan
%                 pos_ = pos_(1:end-1);
%                 fce_ = fce_(1:end-1);
%                 t_all= t_all(1:end-1);
%                 if sum(isnan(pos_) | isnan(fce_)) % if still have force
%                     display('pos_ and fce_ have nan values, abort!');
%                     return
%                 end
%                 
%                 if(0) % see the interp result
%                     fh = figure();
%                     subplot(2,1,1);
%                     plot(pos_t, pos, 'o', t_all, pos_, ':.');
%                     ylabel('position interp');
%                     subplot(2,1,2);
%                     plot(fce_t, fce, 'o', t_all, fce_, ':.');
%                     ylabel('foece interp');
%                 end
%                 
%             end
%             x      = pos_ - pos(1);
%             dx     = diff(x,1,2)/(1/freq);
%             ddx    = diff(x,2,2)/((1/freq).^2);
%             dddx   = diff(x,3,2)/((1/freq).^3);
%             ddddx  = diff(x,4,2)/((1/freq).^4);
%             
%             n = length(ddddx);
%             F = reshape(fce_(1:n),n,1);
%             X = [ones(n,1), ...
%                 reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
%                 reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1), ...
%                 reshape(ddddx(1:n),n,1)];
%             b = [];
%             % calculate 
%             b = (pinv(X'*X)*X'*F);
%             % asign values
%             obj.pred_K = -b(2);
%             obj.pred_D = -b(3);
%             obj.pred_A = -b(4);
%             obj.pred_J = -b(5);
%             obj.pred_S = -b(6);
%             obj.pred_x0 = b(1)/obj.pred_K;
%             if(0)
%                 fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, J=%.3f, S=%.3f, x0=%.3f', ...
%                     obj.tNo, obj.pred_K, obj.pred_D, obj.pred_A, obj.pred_J, obj.pred_S, obj.pred_x0 );
%             end
%         end
        function obj = judgeOffline(obj)
            % OBJ = JUDGEOFFLINE(OBJ)
            % The robot and judge model's lagging makes the judgement does
            % not work well. In this function, I'll offline judge the
            % position and the velocity to re-define the obj.outcome 
            
            % 1. offline judge parameters define
            time_scale = [0.4 0.6];
            pos_dist = 0.01;        % ±1cm
            vel_dist = 0.05;       % 5mm/s ??? is it 5cm or 5mm?
            TS_MOV = 5;
            pos_offset = 0.48;
            idx = obj.data.t_shift > time_scale(1) & obj.data.t_shift < time_scale(2);
            
            % 2. judge begin
                % kick the unfinished trials
            if isempty(obj.data.ts==TS_MOV)
                obj.outcomeo = 0; 
                return;
            end
                % kick the not enough long trials
            if sum(idx) < 100
                obj.outcomeo = 0; 
                return;
            end
            
            if unique(obj.data.ts(idx)==8)
                obj.outcomeo = 0; 
                return;
            end
                
            
                % begin judge
            
            jd(1) = prod(abs(obj.data.x(2,idx) - pos_offset -obj.tarL) < pos_dist);
            jd(2) = prod(abs(obj.data.v(2,idx)- 0 ) < vel_dist);
            obj.outcomeo = jd(1) * jd(2); % should I use &&?
            
            ifplot = 0;
            if(ifplot)
                clf;
                axh(1) = subplot(2,1,1); hold on;
                plot(obj.data.t_shift, obj.data.x(2,:));
                line(time_scale, (pos_offset+obj.tarL+pos_dist)*[1 1], 'color', 'r');
                line(time_scale, (pos_offset+obj.tarL-pos_dist)*[1 1], 'color', 'r');
                title(['trial ' num2str(obj.tNo) 'outcome ' num2str(obj.outcome)]);
                ylabel('position (m)');
                grid on;
                
                axh(2) = subplot(2,1,2); hold on;
                plot(obj.data.t_shift, obj.data.v(2,:));
                line(time_scale, (0+0+vel_dist)*[1 1], 'color', 'r');
                line(time_scale, (0+0-vel_dist)*[1 1], 'color', 'r');
                ylabel('velocity (m/s)');
                grid on;
                
                linkaxes(axh, 'x');
                xlim([0 0.7]);
            end
        end
        function obj = predictMass(obj, ifplot)
            % This function use the force immediate before and after the
            % release to estimate the mass of the subject side. The mass is
            % with the FT subject side plate.  
            % obj = predictMass(obj, ifplot)
            % The raw force data is needed. Current version use 1000Hz force data
            % 1. get raw data
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            % other defined values:
            Mr = 1/(1.413-0.3549); % mass on the robot side % value may inaccurate
            Ms0= 0.3549*Mr;
            y_raw = obj.force_h(2,:);
            yt_raw = obj.force_t;
            y = y_raw; % only on y direction
            y_t = yt_raw;
            % 2. time truncate
            t_min = 0; t_max = 1.5;
            y_truncate = y(y_t>=t_min & y_t<=t_max);
            y_t_truncate = y_t(y_t>=t_min & y_t<=t_max);
            % 3. time shift according to the sudden change peak
            y_d = diff(y_truncate);
            [y_dmin, y_didx] = min(y_d);
            duration = floor(1.2/mean(diff(y_t_truncate)));
            try
                y_shift = y_truncate(y_didx:y_didx+duration-1);
                y_t_shift = y_t_truncate(y_didx:y_didx+duration-1) - y_t_truncate(y_didx); % new 0-nize
            catch 
                display('too short data after release'); 
                return
            end
            % 4. regression model
            % y = e^(\ksai t)*sin(\omega t + \theta)
            [pks,locs] = findpeaks(lowpass(y_shift, pi/100));
            pks_num = length(pks);
            try 
                B0 = -0.69/y_t_shift(locs(2)-locs(1));  % about 1 sec decrease half;
            catch 
                B0 = -0.69/1;
            end
            B3 = max(y_shift);
            B1 = 2*pi/(range(y_t_shift)/pks_num*2); % The period is around 1
            B2 = pi/2; % looks like its dropping value
            X = y_t_shift;
            % 5. regression
            myFit = NonLinearModel.fit(X,y_shift, 'y_shift ~ exp(b0*x1)*sin(b1*x1 + b2)*b3', [B0, B1, B2, B3]);
            y_fit = myFit.Fitted;
            if (ifplot)
%                 subplot(4,1,1);
%                 plot(yt_raw, y_raw);
%                 subplot(4,1,2);
%                 plot(y_t_truncate(1:end-1), diff(y_truncate));
%                 subplot(4,1,3);
%                 plot(y_t_shift, y_shift);
%                 subplot(4,1,4);
%                 plot(y_t_shift, y_shift, 'b');
%                 hold on
%                 plot(y_t_shift', y_fit, 'r');
%                 
                figure(); hold on;
                plot(yt_raw - y_t_truncate(y_didx), y_raw, 'b');
                plot(y_t_shift', myFit.Fitted,'r');
                xlim([-0.1 1]);
                legend('raw data', 'damped sinusoid fit')
                title('force immediatly after release');
                xlabel('time at release (s)'); ylabel('Censored force (N)');
            end
            % 6. Get the Force0+
            Force_0p = y_fit(1);
            % 7. Compare it with Force0-
            force_0n_arr = mean(y_raw(yt_raw>-0.1 & yt_raw<0)); 
            Force_0n = mean(force_0n_arr)*ones(size(Force_0p));
            % Force_0n/Force0p = 1+Ms/Mr;
            Ms = (Force_0n/Force_0p - 1)*Mr; % only known Mr could we ger Ms
            obj.pred_A = Ms;
        end
        function obj = predictStiffDampx0(obj, ifplot)
            % Use regression method (Scott's proach) to estimate stiffness,
            % damping and equilibrium position
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            % 1. get raw data
            f_raw = obj.force_h(2,:);
            ft_raw= obj.force_t;
            x_raw = obj.position_h(2,:);
            xt_raw= obj.position_t;
            x_offset = mean(x_raw(xt_raw>-0.1 & xt_raw<0));
            if (isempty(obj.pred_A))
                try
                    obj = obj.predictMass();
                catch 
                    disp('mass fit failure, return blank');
                    return
                end
            end
            % 2. time truncate
            t_min = 0; t_max = 2; 
            f_truncate = f_raw(ft_raw>=t_min & ft_raw<=t_max);
            f_t_truncate = ft_raw(ft_raw>=t_min & ft_raw<=t_max);
            x_truncate = x_raw(xt_raw>=t_min & xt_raw<=t_max)-x_offset;
            x_t_truncate = xt_raw(xt_raw>=t_min & xt_raw<=t_max);
            % 3. time shift according to the sudden change peak
                % the FT
            f_d = diff(f_truncate);
            [~, y_didx] = min(f_d);
            duration = 1.2; % seconds
            duration_ct = floor(duration/mean(diff(f_t_truncate)));
            f_shift = f_truncate(y_didx:y_didx+duration_ct-1);
            f_t_shift = f_t_truncate(y_didx:y_didx+duration_ct-1) - f_t_truncate(y_didx); % new 0-nize
                % the robot
            x_d = diff(diff(x_truncate(x_t_truncate < 0.03)));
            [~, x_didx] = max(x_d);
            duration_ct = floor(duration/mean(diff(x_t_truncate)));
            x_shift = x_truncate(x_didx:x_didx+duration_ct-1);
            x_t_shift = x_t_truncate(x_didx:x_didx+duration_ct-1) - x_t_truncate(x_didx); % new 0-nize
            
            % 4a. First choice, try resample and do the regression
            iffilterF = 1; %smooth F for better recognize
            if (iffilterF)
            % 4b. If not, try use filtered force and do regression
            % raw force?
            % y = e^(\ksai t)*sin(\omega t + \theta)
            % 4a-1 fit the force
            [pks,locs] = findpeaks(lowpass(x_shift, pi/100), 'MinPeakDistance',100);
            pks_num = length(pks);
            try 
                B0 = -0.69/f_t_shift(locs(2)-locs(1));  % about 1 sec decrease half;
            catch 
                B0 = -0.69/1;
            end
            B3 = max(f_shift);
            B1 = 2*pi/(range(f_t_shift)/pks_num);%*2); % The period is around 1
            B2 = pi/2; % looks like its dropping value
            X1 = f_t_shift; 
            X2 = min(X1):mean(diff(X1)):max(X1);
            % 5. regression
            %myFit1 = NonLinearModel.fit(X1,f_shift, 'f_t_shift ~ exp(b0*x1)*sin(b1*x1 + b2)*b3', [B0, B1, B2, B3]);
            myFit1 = NonLinearModel.fit(X2,f_shift, 'f_t_shift ~ exp(b0*x1)*sin(b1*x1 + b2)*b3', [B0, B1, B2, B3]);
            f_fit = myFit1.Fitted;
            % 4a-2 fit the displacement
            % y = e^(\ksai t)*sin(\omega t + \theta) + a
            B3 = x_shift(end);
            X2 = x_t_shift;
            %   PROBLEMATIC EQUATION HERE! 
            myFit2 = NonLinearModel.fit(X2,x_shift, 'f_t_shift ~ exp(b0*x1)*sin(b1*x1 + b2)*b3+b3', [B0, B1, B2, B3]);
            xfit = myFit2.Fitted;
            end
            
            % resample 
            time = 0:1/500:1.5;
            %tsinf = timeseries(f_shift',f_t_shift);
            %tsinf = timeseries(f_fit,f_t_shift);
            tsinf = myFit1.predict(time');
            tsinx = myFit2.predict(time');
            %tsinx = timeseries(x_shift',x_t_shift);
            %tsinx = timeseries(xfit,x_t_shift);
            %tsoutf = resample(tsinf,time);
            %tsoutx = resample(tsinx,time);
            % regression 
            % Ft+Mx``(t) = Kx0 - Kx(t) - Dx`(t);
            x = tsinx;%tsoutx.Data;
            dx = diff(x)./[diff(time)]';
            ddx = diff(dx)./[diff(time(1:end-1))]';
            x = x(1:length(ddx));
            dx = dx(1:length(ddx));
            f = tsinf(1:length(ddx));%tsoutf.Data(1:length(ddx));
            m = obj.pred_A;
            try 
                y = (f + m*ddx)';
            catch
                display(['Error: m is:' num2str(m)]);
                return
            end
            X = [ones(size(x))'; x'; dx'];
            % A: Kx0, -K, -D
            A = y*X'/(X*X');
            
            obj.pred_K = -A(2);
            obj.pred_D = -A(3);
            obj.pred_x0= A(1)/obj.pred_K;
           
            if (ifplot) % need to plot the raw data, interpretation, and the result
                 subplot(2,1,1); hold on
                 plot(x_t_shift, x_shift, 'b');
                 plot(time, tsinx, 'r');
                 xlim([0,0.5]); 
                 %legend('raw data', 'damped sinusoid fit')
                 title('raw data and fitted data');
                 ylabel('displacement (m)')
                 subplot(2,1,2); hold on;
                 plot(f_t_shift, f_shift, 'b');
                 plot(time, tsinf, 'r');
                 xlim([0,0.5]);
                 legend('raw data', 'damped sinusoid fit')
                 %title('force');
                 ylabel('censored force (N)')
            end
        end
%         function obj = aligntime_wamft(obj)
%             % obj = aligntime_wamft(obj)
%             % align the time of wam and force transducer. 
%             % 1. Find the part they are the same
%             disp('Align the time of wam and ft now');
%             wamt = obj.position_t; 
%             ftt  = obj.force_t;
%             
%             % 2. check whether the time are overlepped, if not, rely on wam
%             % time base. 
%             
%             % 3. save the trancated FT and wam_related variables in obj
%             
%         end

%         function comboTT = getComboTT(obj,sessionScanObj)
%         	% defined: targets = [obj.tarR, obj.tarL, obj.fTh];
%             if length(unique(obj.tarR)) > 1
%                 display(['trial' num2str(obj.tNo) ' have more than 1 tarR']);
%             end
%             tarR = max(obj.tarR);
%             tarL = obj.tarL;
%             fTh  = obj.fTh;
%             tarR_all = sessionScanObj.tarRs; 
%             tarL_all = sessionScanObj.tarLs; 
%             fTh_all  = sessionScanObj.fThs; 
%             try % sometrials do not have tarR
%                 tarR_idx = find(tarR == tarR_all);
%                 tarL_idx = find(tarL == tarL_all);
%                 fTh_idx  = find(fTh  == fTh_all );
%                 comboTT = (tarR_idx-1) * length(tarL_all) * length(fTh_all) + ...
%                           (tarL_idx-1) * length(fTh_all) + ...
%                           fTh_idx;
%             catch
%                 comboTT = nan;
%             end
%             if isempty(comboTT)
%                 comboTT = nan;
%             end
%         end
        
%         function mov_time = getMoveTime(obj)
%             if obj.outcome == 1
%                 time_bgn = obj.time(obj.idx_mov); %...
%                 time_edn = obj.time(obj.idx_hld); %...
%                 mov_time = time_edn - time_bgn;
%             else
%                 mov_time = -1;
%             end
%         end
        
%         function mov_time = getMoveTimeArr(obj)
%             trials_num = length(obj);
%             mov_time = zeros(1,trials_num);
%             for trial_i = 1:length(obj)
%                 trial_tmp = obj(trial_i);
%                 if trial_tmp.outcome == 1
%                     time_bgn = trial_tmp.time(trial_tmp.idx_mov); %...
%                     time_edn = trial_tmp.time(trial_tmp.idx_hld); %...
%                     mov_time(trial_i) = time_edn - time_bgn;
%                 else
%                     mov_time(trial_tmp) = -1;
%                 end
%             end
%         end
%         function frc = getforceVecBeforeRelease(obj)
%             % return 100ms force average before release
%             frc = zeros(3,1);
%             time_bgn = obj.time(obj.idx_mov) - 0.2; %...
%             time_edn = obj.time(obj.idx_mov); %...
%             time_idx = obj.time > time_bgn & obj.time < time_edn;
%             frc = mean(obj.force(:,time_idx), 2); % row avg
%         end
%         function obj = simuTrialusingODE(obj, k, b, m)
%             % obj = simuTrialusingODE(obj, k, b, m)
%             % generate a trial using a specified stiffness, damping and
%             % mass. 
%             
%             %F_list = [3 5 10 13 15 20];
%             %figure(1); hold on;
%             %figure(2); hold on;
%             %k = 226.156; 
%             % b = 10.171;
%             if ~exist('m', 'var')
%                 m = 3.011;
%             end
%             x0 = 0.05; % 5cm;
%             F  = k * 0.05;
%             %m  = 1;     % kg
%             ks = k;   % N/m  % assume subject can alter his/her stiffness
%             kr = b;    % N/(m*s)
%             % specify ODE
%             % m*x'' = ks(x0 - x) - kr*x'
%             % x'(0) = 0;
%             % x(0) = 0;
%             syms x(t)
%             Dx = diff(x);
%             
%             % ode = diff(x, t, 2) == 1/m * ks * (x0 - x) - kr * Dx;
%             ode = diff(x, t, 2) == eval('1/m') * eval('ks')*(eval('x0') - x) - eval('kr') * Dx;
%             cond1 = x(0) == 0;
%             cond2 = Dx(0)== 0;
%             
%             conds = [cond1 cond2];
%             % get solutions
%             xSol(t) = dsolve(ode, conds);
%             xSol = simplify(xSol);
%             
%             t = -2:0.002:3;
%             x_result = xSol(t);
%             x_result_numerical= double(x_result);
%             x_resultd_numerical= diff(x_result_numerical);
%             tdiff = t(1:end-1);
%             x_result_numerical(t<0) = 0;
%             x_resultd_numerical(tdiff<0) = 0;
%             figure();
%             subplot(3,1,1);
%             plot(t, x_result_numerical);
%             title('position');
%             xlim([-0.2, 0.8]);
%             ylabel('m');
%             subplot(3,1,2);
%             plot(t(1:end-1), diff(x_result_numerical)./diff(t));
%             title('velocity');
%             xlim([-0.2, 0.8]);
%             ylabel('m/s');
%             velocity = diff(x_result_numerical)./diff(t);
%             position = x_result_numerical;  
%             t_vel = t(1:end-1);
%             force = diff(velocity)./diff(t_vel) * m;
%             subplot(3,1,3);
%             plot(t_vel(1:end-1), force);
%             title('force');
%             ylabel('N');
%             xlim([-0.2, 0.8]); 
%             
%             % saving datas in obj
%             obj.position_t = t;
%             obj.force_t = t_vel(1:end-1);
%             obj.position_h = position; 
%             obj.force_h = force; 
%             %obj.force_h(obj.force_t<0) = F;
%         end
%         
        function [dat] = export_as_formatted(obj, ifplot)
            % export as a t(trials_num)-by-p(perturbation options) cell mat
            % for each cell, the data format are each trial, which contains:
            %   x: 3-by-N matrix, robot endpoint
            %   v: 3-by-N matrix, robot velocity
            %   f: 3-by-N matrix, force transducer force
            %   Fpert: 1-by-N matrix, perturbation force, 2-by-N in stoc(xy)
            %   ts: 1-by-N matrix, task states
            %   time: 1-by-N matrix, time 
            %   movement onset: the mask that robot start move
            %  -[ ] emg: 8-by-N matrix, emg data
        % consider tiey up a single trial
%         timepoints = obj.data.t;
%         ts_time = timepoints(find(diff(obj.data.ts)));  % time when ts change
        if (~exist('ifplot', 'var'))
            ifplot = 0;
        end
        
        ts_valid = [1:7];
        idx = ismember(obj.data.ts, ts_valid);
%         
%         wam_t = timepoints(ismember(obj.data.ts, ts_valid));
%         wam_t = reshape(wam_t, 1, length(wam_t));
%         time_diff = diff(wam_t);
%         time_diff = [time_diff(1) time_diff];
%         % exclude the time skewed too much
%         tolerance = 0.5;
%         tolerance_range = [1-tolerance, 1+tolerance]*mode(time_diff);
%         vidx = time_diff>tolerance_range(1) & time_diff<tolerance_range(2); % valid idx
%         ctn_flag = diff(find(vidx)); % if has value >1, discontinue
%         if (sum(ctn_flag>1))
%             disp(['IS discontinuety exist in wam time trial' num2str(obj.tNo) ', ABORT!']);
%         end

        % %%% package the data into format
        
        %dat.t = obj.data.t(idx);
        dat.t = obj.data.t(idx);
        %wamt_org= dat.t;
        %t_const = min(dat.t):2e-3:max(dat.t);
        %dat.t = t_const;
        %dat.x = obj.position_h(:,vidx);
        %dat.x = interp1(obj.position_t', obj.position_h', dat.t)';
        dat.x = obj.data.x(:,idx);
        dat.v = obj.data.v(:,idx);
        dat.f = obj.data.f(:,idx);
        dat.ftq=obj.data.ftq(:,idx);
        dat.Fp= obj.data.Fp(:,idx);
        dat.ts= obj.data.ts(:,idx);
        dat.tq= obj.data.tq(:,idx);
        if isfield(obj.data, 'emg')
            dat.emg=obj.data.emg(:,idx);
        end
        if isfield(obj.data, 'ox')
            dat.ox = obj.data.ox(:,idx);
            dat.ov = [ diff(obj.data.ox(1,idx)) ./ diff(obj.data.t(idx));
                       diff(obj.data.ox(2,idx)) ./ diff(obj.data.t(idx));
                       diff(obj.data.ox(3,idx)) ./ diff(obj.data.t(idx));];
            dat.ov = dat.ov(:, [1,1:end]); % make same length
        end
        %dat.x = interp1(wamt_org, dat.x', dat.t, 'linear', 'extrap');
        %dat.v = obj.velocity_h(:,vidx);
        %dat.v = interp1(obj.position_t', obj.velocity_h', dat.t)';
%         try
%             dat.f = (interp1(obj.force_t', obj.force_h', dat.t, 'linear', 'extrap'))'; %... need intropolate
%         catch
%             [C,IA,IC] = unique(obj.force_t);
%             dat.f = (interp1(obj.force_t(IA)', obj.force_h(:,IA)', dat.t, 'linear', 'extrap'))'; %... need intropolate
%             disp(['force_t have' num2str(sum(diff(obj.force_t)==0)) ' overlapping values, ! Need to check!']);
%         end
%         
%         dat.Fp= obj.pertfce_h(:,vidx);
%         dat.Fp= interp1(obj.position_t', obj.pertfce_h', dat.t)';
        %dat.FP = interp1(wamt_org, dat.FP, dat.t, 'linear', 'extrap');
%         ts = reshape(obj.data.ts, 1, length(obj.data.ts));
%         dat.ts = interp1(obj.position_t',ts,dat.t, 'linear', 'extrap');
        %dat.ts = round(interp1(wamt_org, dat.ts, dat.t, 'linear', 'extrap'));
%         dat.tq= obj.wamtqe_h(:,vidx);
%         dat.tq= interp1(obj.position_t', obj.wamtqe_h', dat.t)';
        %dat.tq= interp1(wamt_org, dat.tq, dat.t, 'linear', 'extrap');
        if (obj.ifpert==0 || obj.ifpert==1)
            dat.mvst= (dat.ts==5 | dat.ts==6); % moveent start
        elseif(obj.ifpert==2) % stochastic pert
            dat.mvst= (dat.ts==5 | dat.ts==6);
        end
        % the optotrak data
%         if ~isempty(obj.opt)
%             ox = interp1(obj.opt.datah.t, obj.opt.datah.x, dat.t, 'linear', 'extrap');
%             oy = interp1(obj.opt.datah.t, obj.opt.datah.y, dat.t, 'linear', 'extrap');
%             oz = interp1(obj.opt.datah.t, obj.opt.datah.z, dat.t, 'linear', 'extrap');
%             dat.ox = [ox; oy; oz];
%         end
        

        % deal with errors (that back to ts3)
        if (~isempty(find(dat.ts==3 & [diff([1 dat.ts]) == -1], 1 )))
            ts = dat.ts;
            %idx_forcefail = [find(ts==3 & [diff([1 ts]) == -1] )];  % index that back to ts3
            %idx_forcefail = idx_forcefail(end);
            idx_forceadv = find(ts==4 & [diff([1 ts]) == 1]);      % going into force hold zone
            idx_forceadv = idx_forceadv(end);
            idx_forcein = find(ts==4);                             % in force zone (4)
            idx_forcein  = idx_forcein(1);
            %dat.Fp(2,idx_forcein:idx_forcefail) = 0;
            dat.Fp(1:2,idx_forcein:idx_forceadv) = 0;

        end
        
%         % check the plot
%         ifplot = 0;
%         if (ifplot)
%             clf;
%             axhl(1) = subplot(3,1,1); hold on;
%             plot(obj.time, obj.states_arr, 'r');
%             plot(dat.t, dat.ts, 'b');
%             legend('RTMA-TS', 'WAM-TS');
%             title('task states');
%             axhl(2) = subplot(3,1,2); hold on;
%             %plot(dat.t, dat.Fp, 'r');
%             plot(dat.t, dat.f, 'r');
%             %plot(obj.force_t, obj.force_h(:,:), 'b');
%             title('perturbation force');
%             axhl(3) = subplot(3,1,3); hold on;
%             plot(dat.t, dat.x(2,:), 'r');
%             %plot(obj.position_t, obj.position_h(2,:), 'b');
%             title('endpoint position');
%             %linkaxes(axhl, 'x');
%             suptitle(['trial ' num2str(obj.tNo)]);
%         end
        
        % if pert happen eailier
        %if ~isempty(setdiff(unique((dat.ts(dat.Fp(2,:)~=0))),3)) && obj.ifpert==1
        %if ~isempty(setdiff(unique((dat.ts(dat.Fp(2,:)~=0))),4)) && obj.ifpert==1 % perturb at ts4
        %if ~isempty(setdiff(unique((dat.ts(dat.Fp(2,:)~=0))),[3,4])) && obj.ifpert==1 % perturb at either ts3 or ts4
        %if ~isempty(setdiff(unique((dat.ts(dat.Fp(2,:)~=0))),[3,4,5])) && obj.ifpert==1 % perturb at either ts3 or ts4
        if ~isempty(setdiff(unique((dat.ts(dat.Fp(2,:)~=0))),[3,4,5,6,7])) && obj.ifpert==1 % perturb at either ts3 or ts4, 5, 6, 7
            dat.Fp = zeros(size(dat.Fp));
        end
        
        % another way for avoid too eairly pert
        %if abs(nanmean(dat.f(2,dat.Fp(2,:)~=0) ))<5 && obj.ifpert==1
        %    dat.Fp = zeros(size(dat.Fp));
        %end
        
        ifplot = 0;
        outcome_name = 'sf';
        if (ifplot)
%             subplot(2,1,1); 
%             plot(obj.position_t, obj.position_h(2,:));
%             subplot(2,1,2);
%             plot(obj.position_t, obj.pertfce_h(2,:));
%             plot(obj.force_t, obj.force_h(2,:));
%             
%             subplot(2,1,1);
%             plot(dat.t, dat.x);
%             subplot(2,1,2);
%             plot(dat.t, dat.f);
%             
%             subplot(2,1,1);
%             plot(obj.position_t, obj.position_h(2,:));
%             subplot(2,1,2);
%             plot(dat.t, dat.x(2,:));
              clf;
%               time = dat.t;
              t = obj.data.t_shift(idx);
              axh(1) = subplot(4,1,1);
              plot(t, dat.Fp);
%               plot(t, dat.mvst);
              title(['trial' num2str(obj.tNo) ' :' outcome_name(2-obj.outcome)]);
              grid on;
              ylabel('Fp (N)' );
              axh(2) = subplot(4,1,2); hold on;
              plot(t, dat.x(2,:));
              line([0.4 1.0], (0.48+obj.tarL+0.01)*[1 1], 'color', 'r');
              line([0.5 1.0], (0.48+obj.tarL-0.01)*[1 1], 'color', 'r');
              line([0.4 1.0], (0.48+obj.tarL+0.005)*[1 1], 'color', 'g');
              line([0.5 1.0], (0.48+obj.tarL-0.005)*[1 1], 'color', 'g');
              ylabel('position (m)');
              grid on;
              axh(3) = subplot(4,1,3); hold on;
              plot(t, dat.v(2,:));
              line([0.5 1.0], [0.05 0.05], 'color', 'r');
              line([0.4 1.0], [-0.05 -0.05], 'color', 'r');
              ylabel('velocity (m/s)');
              grid on;
              axh(4) = subplot(4,1,4);  hold on;
              plot(t, dat.f(2,:), 'Marker', '.');
              grid on;
              ylabel('Force (N)')
              
              linkaxes(axh, 'x');
              % xlim for better read
%               xlim([[-0.01 0.02]]);
%               xlim([-0.2 1]);
              xlim([-0.2 2]);

%             subplot(2,1,1);
%             plot(obj.force_t', obj.force_h');
%             subplot(2,1,2);
%             plot(dat.t, dat.f);
%%%%%%%%%%%%%% condition: abs(t(dat.Fp(2,:)==-12) - 0.2) < 0.02
        end
        ifplot = 0; 
        if (isfield(dat, 'ox'))
        if (~isempty(dat.ox))
            
             if (ifplot)
              clf;
%               time = dat.t;
              t = obj.data.t_shift(idx);
              axh(1) = subplot(4,1,1);
              plot(t, dat.Fp);
%               plot(t, dat.mvst);
              title(['trial' num2str(obj.tNo) ' :' outcome_name(2-obj.outcome)]);
              grid on;
              ylabel('Fp (N)' );
              axh(2) = subplot(4,1,2); hold on;
              plot(t, dat.x(2,:), 'b.');    % robot
              plot(t, dat.x(2,:), 'r.');    % optotrak
              line([0.4 1.0], (0.48+obj.tarL+0.01)*[1 1], 'color', 'r');
              line([0.5 1.0], (0.48+obj.tarL-0.01)*[1 1], 'color', 'r');
              line([0.4 1.0], (0.48+obj.tarL+0.005)*[1 1], 'color', 'g');
              line([0.5 1.0], (0.48+obj.tarL-0.005)*[1 1], 'color', 'g');
              ylabel('position (m)');
              grid on;
              axh(3) = subplot(4,1,3); hold on;
              lnh(1) = plot(t, dat.v(2,:), 'b.');    % robot
              lnh(2) = plot(t, dat.ov(2,:), 'r.');    % optotrak
              line([0.5 1.0], [0.05 0.05], 'color', 'r');
              line([0.4 1.0], [-0.05 -0.05], 'color', 'r');
              legend(lnh, 'robot', 'optotrak');
              ylabel('velocity (m/s)');
              grid on;
              axh(4) = subplot(4,1,4);  hold on;
              plot(t, dat.f(2,:), 'Marker', '.');
              grid on;
              ylabel('Force (N)')
              
              linkaxes(axh, 'x');
              % xlim for better read
%               xlim([[-0.01 0.02]]);
%               xlim([-0.2 1]);
              xlim([-0.2 2]);

%             subplot(2,1,1);
%             plot(obj.force_t', obj.force_h');
%             subplot(2,1,2);
%             plot(dat.t, dat.f);
%%%%%%%%%%%%%% condition: abs(t(dat.Fp(2,:)==-12) - 0.2) < 0.02
            end
        end % end of if
        end
        end
        
        function params = export_trial_params(obj)
            % return the parameters in this trial 
            params.wamKp = obj.wamKp;
            params.wamBp = obj.wamBp;
            params.pertf = obj.pert_f;
            
        end
        
        
        function obj = dealForceException(obj)
            % obj = dealForceException(obj)
            % for sessions that not collect reasonable force. 
            % minus the force bias from the ts7
            idx_ = obj.data.ts == 7 & ~isnan(obj.data.x(1,:));
            
            ifplot = 0;
            if(ifplot)
                plot(obj.data.t_shift, idx_);
            end
            
            force_offset = mean(obj.data.f(:,idx_),2);
            obj.data.f = obj.data.f - force_offset;
            
            if (ifplot)
                plot(obj.data.t_shift, obj.data.f(2,:));
            end
        end
        %%% plot
%         function axh = visualizeFrcVelDelay(obj)
%             % axh = visualizeFrcVelDelay(obj)
%             % Visualize force and velocity as a overlap lines, and by this
%             % determine whether the alignment is good. 
%             % In the current situation, force and neural activity only has
%             % alignment according to the RTMA system, which could be
%             % delayed in network. 
%             
%             % pre-process data
%             % same magnitude
%             t_bgn = 0;
%             t_edn = 0.4;
%             freq = 500;     % 500Hz, the same as WAM
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             
%             fce_t = obj.force_t(fce_t_idx);
%             pos = obj.position_h(obj.xyi,pos_t_idx);
%             fce = obj.force_h(obj.xyi,fce_t_idx);
%             vel = obj.velocity_h(obj.xyi,pos_t_idx);
%             
%             % here only for y direction
%             pos_reMag = 1/range(pos) * pos; 
%             vel_reMag = 1/range(vel) * vel;
%             frc_reMag = 1/range(fce) * fce; 
%             pos_reNorm = pos_reMag - mean(pos_reMag);
%             frc_reNorm = frc_reMag - mean(frc_reMag);
%             
%             axh = figure(); hold on;
%             plot(pos_t, pos_reNorm(:)); 
%             plot(fce_t, frc_reMag(:));
%             plot(pos_t, vel_reMag(:));
%             
%             legend('position\_renorm', 'force\_remag', 'velocity\_remag');
%             title('timeAlign comparation in trial');
%             xlabel('time after release');
%         end
%         function axh = plotPredictedForceOnPosition(obj)
%             % use regression terms to get the predicted force
%             t_bgn = 0;
%             t_edn = 0.4;
%             freq = 500;     % 500Hz, the same as WAM
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             fce_t = obj.force_t(fce_t_idx);
%             pos = obj.position_h(obj.xyi,pos_t_idx);
%             pos = pos - pos(1); % remove the offset
%             fce = obj.force_h(obj.xyi,fce_t_idx);
%             try
%                 ifsame = min(pos_t == fce_t); % one 0, all 0
%             catch 
%                 ifsame = 0;
%             end
%             
%             
%             if (~ifsame)
%                 t_all = t_bgn:(1/freq):t_edn;
%                 % resample at time t
%                 pos_ = interp1(pos_t, pos, t_all); % .... process here
%                 fce_ = interp1(fce_t, fce, t_all); 
%                 % chances the last pos_ and fce_ is nan
%                 pos_ = pos_(1:end-1);
%                 fce_ = fce_(1:end-1);
%                 t_all= t_all(1:end-1);
%                 if sum(isnan(pos_) | isnan(fce_)) % if still have force
%                     err('pos_ and fce_ have nan values, abort!');
%                 end
%             end
%             x      = pos_;
%             dx     = diff(x,1,2)/(1/freq);
%             ddx    = diff(x,2,2)/((1/freq).^2);
%             dddx   = diff(x,3,2)/((1/freq).^3);
%             
%             n = length(dddx);
%             X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
%                 reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
%             b = [obj.pred_x0*obj.pred_K;
%                 -obj.pred_K;
%                 -obj.pred_D;
%                 -obj.pred_A;
%                 -obj.pred_J];
%             F = X*b;
%             t_all_ = t_all(1:length(F));
%             axh = figure('Visible', 'on');
%             hold on;
%             plot(t_all, fce_, 'b', 'LineWidth', 3); 
%             plot(t_all_, F, 'r--', 'LineWidth', 3); 
%             xlabel('time (s)');
%             ylabel('force (N)');
%             legend('origin force', 'regressed force');
%             title(['origin and regress force trial' num2str(obj.tNo)]);
%             
%         end
%         function axh = plotPredictedForceOnPositionS(obj)
%             % use regression terms to get the predicted force
%             % with snap term
%             t_bgn = 0;
%             t_edn = 0.4;
%             freq = 500;     % 500Hz, the same as WAM
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             fce_t = obj.force_t(fce_t_idx);
%             pos = obj.position_h(obj.xyi,pos_t_idx);
%             pos = pos - pos(1); % remove the offset
%             fce = obj.force_h(obj.xyi,fce_t_idx);
%             try
%                 ifsame = min(pos_t == fce_t); % one 0, all 0
%             catch 
%                 ifsame = 0;
%             end
%             
%             
%             if (~ifsame)
%                 t_all = t_bgn:(1/freq):t_edn;
%                 % resample at time t
%                 pos_ = interp1(pos_t, pos, t_all); % .... process here
%                 fce_ = interp1(fce_t, fce, t_all); 
%                 % chances the last pos_ and fce_ is nan
%                 pos_ = pos_(1:end-1);
%                 fce_ = fce_(1:end-1);
%                 t_all= t_all(1:end-1);
%                 if sum(isnan(pos_) | isnan(fce_)) % if still have force
%                     err('pos_ and fce_ have nan values, abort!');
%                 end
%             end
%             x      = pos_;
%             dx     = diff(x,1,2)/(1/freq);
%             ddx    = diff(x,2,2)/((1/freq).^2);
%             dddx   = diff(x,3,2)/((1/freq).^3);
%             ddddx  = diff(x,4,2)/((1/freq).^4);
%             n = length(ddddx);
%             X = [ones(n,1), ...
%                 reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
%                 reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)...
%                 reshape(ddddx(1:n),n,1)];
%             b = [obj.pred_x0*obj.pred_K;
%                 -obj.pred_K;
%                 -obj.pred_D;
%                 -obj.pred_A;
%                 -obj.pred_J;
%                 -obj.pred_S];
%             F = X*b;
%             t_all_ = t_all(1:length(F));
%             axh = figure('Visible', 'on');
%             hold on;
%             plot(t_all, fce_, 'b', 'LineWidth', 3); 
%             plot(t_all_, F, 'r--', 'LineWidth', 3); 
%             xlabel('time (s)');
%             ylabel('force (N)');
%             legend('origin force', 'regressed force');
%             title(['origin and regress force trial' num2str(obj.tNo)]);
%             
%         end
%         function axh = plotPredictedForceOnPosition_lowsample(obj, sampleRate)
%             % only plot, but not generate data here
%             % predict and plot the Force and Position on a low sample rate
%             % for example, FT runs at 500Hz whereas RTMA runs 50Hz if we
%             % low sample at the rate of 500/50, it would be runs well 
%             if (nargin<2)
%                 sampleRate = 10;
%             end
%             % %% part1: generate the low_sampled data. 
%             t_bgn = 0;
%             t_edn = 0.4;
%             freq = 500;     % 500Hz, the same as WAM
%             pos_t_idx = obj.position_t>=t_bgn & obj.position_t<=t_edn;
%             fce_t_idx = obj.force_t>=t_bgn & obj.force_t<=t_edn;
%             pos_t = obj.position_t(pos_t_idx);
%             fce_t = obj.force_t(fce_t_idx);
%             pos = obj.position_h(obj.xyi,pos_t_idx);
%             fce = obj.force_h(obj.xyi,fce_t_idx);
%             try
%                 ifsame = min(pos_t == fce_t); % one 0, all 0
%             catch 
%                 ifsame = 0;
%             end
%             
%             if (~ifsame)
%                 t_all = t_bgn:(1/freq):t_edn;
%                 % resample at time t
%                 pos_ = interp1(pos_t, pos, t_all); % .... process here
%                 fce_ = interp1(fce_t, fce, t_all); 
%                 % chances the last pos_ and fce_ is nan
%                 pos_ = pos_(1:end-1);
%                 fce_ = fce_(1:end-1);
%                 t_all= t_all(1:end-1);
%                 if sum(isnan(pos_) | isnan(fce_)) % if still have force
%                     display('pos_ and fce_ have nan values, abort!');
%                     return
%                 end
%             end
%             pos    = pos_ - pos(1);
%             fce    = fce_;
%             posR   = pos(1:sampleRate:end);
%             fceR   = fce(1:sampleRate:end);
%             resample_freq = 500/sampleRate;
%             
%             % %% part2: estimate through lines
%             
%             x      = posR;
%             dx     = diff(x,1,2)/((1/resample_freq).^1);
%             ddx    = diff(x,2,2)/((1/resample_freq).^2);
%             dddx   = diff(x,3,2)/((1/resample_freq).^3);
%             
%             n = length(dddx);
%             F = reshape(fceR(1:n),n,1);
%             X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
%                 reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
%             b = [];
%             % calculate 
%             b = (pinv(X'*X)*X'*F);
%             % asign values
%             pred_K = -b(2);
%             pred_D = -b(3);
%             pred_A = -b(4);
%             pred_J = -b(5);
%             pred_x0 = b(1)/obj.pred_K;
%             if(0)
%                 fprintf('trial%03d, K=%.3f, B=%.3f, M=%.3f, J=%.3f, x0=%.3f', ...
%                     obj.tNo, obj.pred_K, obj.pred_D, obj.pred_A, obj.pred_J, obj.pred_x0 );
%             end
%             % %% calculating the predicted value
%             X = [ones(n,1), reshape(x(1:n),n,1), reshape(dx(1:n),n,1), ...
%                 reshape(ddx(1:n),n,1), reshape(dddx(1:n),n,1)];
%             b = [pred_x0*pred_K;
%                 -pred_K;
%                 -pred_D;
%                 -pred_A;
%                 -pred_J];
%             F = X*b;
%             t_all_ = t_all(1:length(F));
%             axh = figure('Visible', 'on');
%             hold on;
%             % %% part3: figure the low_sampled data
%             plot(t_all, fce_, ':.');   % origin
%             t_rsp = t_bgn:(1/resample_freq):t_edn;             % t_resample
%             t_prd = t_rsp;                                     % for prediction
%             if length(t_rsp) > length(fceR)
%                 t_rsp = t_rsp(1:length(fceR));
%             end
%             plot(t_rsp, fceR, 'k', 'LineWidth', 3);                   % low-sample
%             if length(t_prd) > length(F)
%                 t_prd = t_prd(1:length(F));
%             end
%             plot(t_prd, F, 'b', 'LineWidth', 3);
%             xlabel('time (s)');
%             ylabel('force (N)');
%             legend('origin force', 'low sample force', 'pred force');
%             title(['origin and regress force trial' num2str(obj.tNo)]);
%         end
%         function axh = plotPredictedForceODE(obj)
%             % solve the force using differential equation
%             % wait to be written
%             obj;
%             axh = -1;
%         end
        function [axh, lnh] = plotRobotEndpointTraj(obj, axh, clr)
            %printf('begin the plot');
            if nargin == 1
                axh = figure();
                clr = 'r';
            end
            center = [-0.513, 0.483];
            figure(axh);
            % only plot the movement idx
            idx = find(obj.data.ts == 5 | obj.data.ts == 6);
            idxstt = idx(1); 
            idxedn = idx(end); 
%             timestt = obj.data.t_shift(idxstt); 
%             timeedn = obj.data.t_shift(idxedn);
%             [~, idxstt_] = min(abs(timestt - obj.position_t)); 
%             [~, idxedn_] = min(abs(timeedn - obj.position_t)); 
            
            % plot(obj.position_h(1,:) - center(1), obj.position_h(2,:) - center(2), '*', 'color', clr); % all the points
            hold on;
            lnh = plot(obj.data.x(1,idxstt:idxedn) - center(1), ...
                obj.data.x(2,idxstt:idxedn) - center(2), ...
                'color', clr); % all the points
            xlim([-0.12, 0.12]); ylim([-0.12, 0.12]);
            title('');
            xlabel('x (m)');
            ylabel('y (m)');
            
        end
%         function axh = plotRobotEndpointTrajRot(obj, axh, rot, cl)
%             % obj.plotRobotEndpointTrajRot(axh, rot, cl)
%             %printf('begin the plot');
%             
%             if nargin == 1
%                 axh = figure();
%             elseif nargin == 2
%                 rot = 0;
%             end
%             % data
%             if exist('cl', 'var')
%                 line_col = cl;
%             else
%                 if obj.outcome == 1
%                     line_col = 'b';
%                     if (rot ~= 0)
%                         line_col = 'g';
%                     end
%                 end
%             end
%             if isempty(obj.idx_mov) || isempty(obj.idx_end)
%                 return
%             end
%             time_bgn = obj.time(obj.idx_mov); %...
%             time_edn = obj.time(obj.idx_end); %...
%             time_idx = (obj.position_t > time_bgn) & (obj.position_t < time_edn); 
%             center = [-0.513, 0.483];
%             x = obj.position_h(1,time_idx) - center(1); 
%             y = obj.position_h(2,time_idx) - center(2); 
%             p_rot = [cos(rot), -sin(rot); sin(rot), cos(rot)] * [x; y];
%             
%             % plot
%             figure(axh);
%             plot(p_rot(1,:), p_rot(2,:), line_col);
%             xlim([-0.12, 0.12]); ylim([-0.12, 0.12]);
%             title('');
%             xlabel('x (m)');
%             ylabel('y (m)');
%         
%         end
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