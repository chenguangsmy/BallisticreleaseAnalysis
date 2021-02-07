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
        comboNo
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
        
    end
    
    methods
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
            obj.tarR    = unique(sessionScanObj.Data.TaskJudging.Target(5, obj.bgn:obj.edn));          % target-rotation
            obj.tarL    = sort(unique(nonzeros(sessionScanObj.Data.TaskJudging.Target(6, obj.bgn:obj.edn))));    % target-length
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