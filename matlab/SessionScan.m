classdef (HandleCompatible)SessionScan < handle
    %VARSCAN scanning some variables in the formatted data
    %   Detailed explanation goes here

    properties
        %%% marking vars
        ssnum
        %%% session stat
        Data
        time
        trials_num
        duration
        duration_avg
        %%% task targetst
        tarRs
        tarLs
        fThs
        %%% task variables 
        hand_pos        % position read from WAM endpoint
        hand_pos_offset % the center_pos for WAM endpoint
        force           % force in the force transducer
        FTrot_M = ...   % global: x-right, y-front, z-up, FT_base x-backup, y-frontup, z-left
            [0          0           cosd(180)
            -sind(45)   cosd(45)    0
            cosd(45)    sind(45)    0];
        taskState
        trials TrialScan% member function
        %%% perturbation variables
        pert_state = 3  % only perturb at force ramp pert_state == 3
        pert_time       % perturbation start and end time
        pert_rdt        % perturbation start read time, (same with wam)
        pertCond
        %%% other modules
        ft              % object of force
        wam             % object of wam
        emg             % object of emg
        opt             % object of opt (OPTOTRAK)
        %force_h         % time for wam seperate data
        force_t         % time for ft seperate data
        %wamp_h          % highly sampled from wam
        %wamv_h
        %wamt_h
        wam_t
        %wam_ts
        emg_h
        emg_t
        data            % all the data including force and wam information (aligned)
        %%% other variables
        endpoint0 = [-0.517 0.481 0.0] 
        col_vec = colormap('lines');
        badTrials = [1];       % bad trial, cull in data
%         badTrials = [1, 278];       % FOR SS1898
        stateNames = {'Begin', 'Present', 'FrcRamp', 'FrcHold', 'Move', 'Hold', 'End', 'Reset'};
        
    end
    
    methods
        %%% process
        function obj = SessionScan(ss_num, badTrials) %(inputArg1,inputArg2)
            %VARSCAN Construct an instance of this class
            obj.ssnum = ss_num;
            %   Detailed explanation goes here
            % load the INTERMEDIATE file for the time 
            file_name = ['KingKong.0' num2str(obj.ssnum) '.mat']; % an examplary trial
            file_dir_formmed = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/';
            file_dir_int = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate/';
            data1 = load([file_dir_int  file_name], 'Data');
            wam_sendt = data1.Data.QL.Headers.BURT_STATUS.send_time;
            wam_recvt = data1.Data.QL.Headers.BURT_STATUS.recv_time;
            wam_rdt   = double(data1.Data.QL.Data.BURT_STATUS.RDT);
            wam_intm = [wam_sendt; wam_recvt; wam_rdt];
            ft_sendt = data1.Data.QL.Headers.FORCE_SENSOR_DATA.send_time;
            ft_recvt = data1.Data.QL.Headers.FORCE_SENSOR_DATA.recv_time;
            % load synchrony signals from the intermediate message
            try  
                % from message
                MID_NETBOX = 67;
                msgidx = data1.Data.QL.Headers.TIME_SYNC.src_mod_id == MID_NETBOX;
                ft_synctime1 = data1.Data.QL.Data.TIME_SYNC.tleading(msgidx); % also, choose MID to do this
                ft_synctime2 = data1.Data.QL.Data.TIME_SYNC.tlasting(msgidx);
                ft_synctime = mean([ft_synctime1 ft_synctime2]);
                ifplot = 1;
                if (ifplot)
                    clf;
                    subplot(2,1,1); 
                    plot(ft_synctime1, ft_synctime2, '*');
                    title('tBeforeSend vs tAfterSend');
                    subplot(2,1,2); 
                    plot(1:length(ft_synctime1), ft_synctime2-ft_synctime1);
                    title('duration in two times');
                end
            catch 
                disp('Hardware synchrony Message error, please check');
            end
            
            try 
                % from hardware
                
            end
            
            ft_rdt   = double(data1.Data.QL.Data.FORCE_SENSOR_DATA.rdt_sequence);
            ft_intm  = [ft_sendt; ft_recvt; ft_rdt];
            %clear file_name  file_dir data1
            
            
            fname0 = ([file_dir_formmed '/' file_name]);
            fname1 = ([file_dir_int '/' file_name]);
            flag_progress = 1;      % progress display
            try 
                load(fname, 'Data');
            catch 
                load(fname0, 'Data');
            end
            try 
                Data1 = load(fname1, 'Data');
                Data1 = Data1.Data; % load the Intermediate data, not sure useful or not.
            catch
                disp('no intermediate data');
            end
            % objects
            flag_noFW_data = 0;
            try
                if (flag_progress)
                    disp('Loading raw FT and WAM...');
                end
                obj.ft = SessionScanFT(obj.ssnum);
                obj.wam = SessionScanWam(obj.ssnum);
            catch 
                disp('no ft and wam data here! ');
                flag_noFW_data = 1;
            end
            % cg: currently we are not focusing on the EMG and OPT -Nov2021
            try 
                if (flag_progress)
                    disp('Loading (intermed/raw) EMG...');
                end
                obj.emg = SessionScanEMG(obj.ssnum);
            catch
%                 disp('no EMG data here! ')
            end
            
            try 
                if (flag_progress)
                    disp('Loading (intermed/raw) OPTOTRAK...');
                end
                obj.opt = SessionScanOPT(obj.ssnum);
            catch 
%                 disp('no OPT data here! '); 
            end
                
            % other data

            obj.Data = Data;
            % obj.time = Data.Time; % how to make this time same with the bk time?
            timetmp = Data.SpikeTimestamp;  % blackrock time
            timetmp1 = interp1(find(~isnan(timetmp)),timetmp(~isnan(timetmp)), 1:length(timetmp), 'linear', 'extrap');
            ifplot = 0;
            if(1)
                clf;
                hold on;
                plot(1:length(timetmp), timetmp, '*');
                plot(1:length(timetmp), timetmp1, '.');
                legend('original', 'processed');
                title('original and extrapolated spike times');
            end
            obj.time = timetmp1; % have problems, don't do that
            
            TrialNo = Data.TrialNo;
            obj.trials_num = max(TrialNo);
            try
                obj.hand_pos_offset = Data.Position.Center(:,~isnan(Data.Position.Center(1,:)));
                obj.hand_pos_offset = obj.hand_pos_offset(:,1);
            catch
                disp('no position info!');
            end
            obj.taskState.Values = Data.TaskStateCodes.Values;              % cg: what does it mean 99?
            obj.taskState.trialNo = Data.TrialNo;
            obj.force = rotateAxisForce(obj);                                % cg: .... change it into aligned one
            obj.duration = max(obj.time);
            obj.tarRs = unique(obj.Data.TaskJudging.Target(5,obj.Data.TaskStateMasks.Move)); %have 0
            obj.tarLs = unique(obj.Data.TaskJudging.Target(6,obj.Data.TaskStateMasks.Move)); % deviate 0
            obj.fThs  = unique(obj.Data.TaskJudging.Target(4,obj.Data.TaskStateMasks.Move)); % deviate 0
            try 
                obj.pertCond.pertdx0_mag = unique(Data.TaskJudging.pertdx0_mag(obj.Data.TaskStateMasks.FrcHold));
                obj.pertCond.wamKp = unique(Data.TaskJudging.wamKp(obj.Data.TaskStateMasks.FrcHold));
                obj.pertCond.wamBp = unique(Data.TaskJudging.wamBp(obj.Data.TaskStateMasks.FrcHold));
            catch 
                obj.pertCond.pertdx0_mag = [];
                obj.pertCond.wamKp = [];
                obj.pertCond.wamBp = [];
            end
            try 
                obj.pertCond.pertdf= Data.TaskJudging.pertdf_mag;
            catch 
                obj.pertCond.pertdf= [];
            end
                

            % processing 
            trials_all = setdiff(unique(TrialNo), 0);
            [obj, syncflag] = updateTimeBlackRock(obj);
%             if (~isempty(obj.ft))
%                 if (flag_progress)
%                     disp('FT High Sample...');
%                 end
%             %    obj = highSampleForce(obj, obj.ft, ft_intm);
%             end
%             if (~isempty(obj.wam))
%                 if (flag_progress)
%                     disp('WAM High Sample...');
%                 end
%             %    obj = highSampleWam(obj, obj.wam, wam_intm);
%             end
%             if (~isempty(obj.emg))
%                 if (flag_progress)
%                     disp('EMG High Sample...');
%                 end
%                 obj = highSampleEmg(obj, obj.emg);
%             end
            
            % should interpolate all the data here !!!!!
            obj = interpData(obj);

            
            if (flag_progress)
                    disp('Trialfy...');
            end
            
            percent_prev = 0;
            for trial_i = 1:length(trials_all)
                trial_percent = floor(trial_i/length(trials_all) * 20)*5;
                if ~(trial_percent == percent_prev) % avoid showing a lot.
                    fprintf('  %02d%%...', trial_percent);
                    percent_prev = trial_percent;
                end
                obj.trials(trial_i) = TrialScan(obj, trial_i);
                % align to mov
                obj.trials(trial_i) = alignMOV(obj.trials(trial_i));
                %obj.trials(trial_i) = alignPertInit(obj.trials(trial_i));
                if (trial_i == length(trials_all)) % last trial
                    fprintf('  100%%  FINISHED!\n');
                end
            end
            if (nargin>1) %specify bad trials
                for trial_i = badTrials
                    obj.trials(trial_i).outcome = 0;
                end
            else
                obj.trials(1).outcome = 0; % 1st trials are bad for force align.
            end
            %
            % plots:
            %   whole session plot
%             axh = taskForceData(obj);
            % axh = taskForceDatah(obj, axh);
            
%             axh = taskEndpointPosition(obj);
            %axh = taskEndpointPositionh(obj, axh);
            % taskEndpointPosition_relative(obj);
%             taskStateMuskFig(obj);
            % taskJointPosition_relateve(obj);
            %   trialfy plot
            % axh = plotTrialfyPositionh_xy(obj);
%             axh = plotTrialfyForce_xy(obj);
%             axh = plotTrialfyForceh_xy(obj);
            %[axhF, axhP] = plotSameTrial(obj);
            
            % remove error task conditions to avoid code error 
            obj = obj.dealingSessionsExceptions(); % nothing here 

        end
        function obj = dealingSessionsExceptions(obj)
            % solve some data-code inconsistant problem, specify for each
            % sessions

        end
        function [sT, tT, sR] = getConditionalSucessTrials(obj) 
            % tobe re-written for the multi sessions processing
% %             % getConditionalSucessTrials(obj) 
% %             % This function is for display the sucess rate of the trials . 
% %             % [sucessTrials, totalTrials, sucessRate] = getConditionalSucessTrials(obj) 
% %             % sT: sucessTrials, 4(directions)-by-n(targetnum)-by-m(forcenum)
% %             % tT: totalTrials, 4(directions)-by-n(targetnum)-by-m(forcenum)
% %             % sR: sucessRate, 4(directions)-by-n(targetnum)-by-m(forcenum)
% %             
% %             tard = setdiff([obj.tarRs], 0); % potential BUG here if both 0 and others are targets direction!
% %             if isempty(tard)
% %                 tard = 0;
% %             end
% %             tar_num = length(tard);
% %             tarl = obj.tarLs;
% %             tarl_num = length(tarl);
% %             fThs = obj.fThs;
% %             fThs_num = length(fThs);
% %             sT = zeros(tar_num, tarl_num, fThs_num);    % sucessful trials
% %             tT = zeros(tar_num, tarl_num, fThs_num);    % total trials
% %             sR = zeros(tar_num, tarl_num, fThs_num); 
% %             % copy all the tarR, tarL, fTh from all trials first
% %             tarR = zeros(1, obj.trials_num);
% %             tarL = zeros(1, obj.trials_num);
% %             fTh  = zeros(1, obj.trials_num);
% %             for trial_i = 1:obj.trials_num
% %                 if isempty(obj.trials(trial_i).tarR)
% %                     tarR(trial_i) = -1;
% %                 else
% %                     tarR(trial_i) = obj.trials(trial_i).tarR;
% %                 end
% %                 if isempty(obj.trials(trial_i).tarL)
% %                     tarL(trial_i) = -1;
% %                 else
% %                     tarL(trial_i) = obj.trials(trial_i).tarL;
% %                 end
% %                 if isempty(obj.trials(trial_i).fTh) 
% %                     fTh(trial_i) = -1;
% %                 else
% %                     fTh(trial_i) = obj.trials(trial_i).fTh;
% %                 end
% %             end
% %             for tard_i = 1:tar_num
% %                 for tarl_i = 1:tarl_num
% %                     for fThi = 1:fThs_num
% %                         tT(tard_i, tarl_i, fThi) = ...
% %                             sum(tarR == tard(tard_i) &...
% %                             tarL == tarl(tarl_i) & ...
% %                             fTh  == fThs(fThi));
% %                         sT(tard_i, tarl_i, fThi) = ...
% %                             sum(tarR == tard(tard_i) &...
% %                             tarL == tarl(tarl_i) & ...
% %                             fTh  == fThs(fThi) & ...
% %                             [obj.trials.outcome] == 1);
% %                         sR(tard_i, tarl_i, fThi) = sT(tard_i, tarl_i, fThi)/tT(tard_i, tarl_i, fThi);
% %                     end
% %                 end
% %             end
% %             sR_2d = reshape(sR(1,:,:), size(sR, 2), size(sR, 3));
% %             sR_table = [[obj.fThs]', sR_2d'];
% %             % display rate using table
% %             display(['For session' num2str(obj.ssnum)]);
% %             VarNames = {'Force (N)', 'tar 2.5 (cm)', 'tar 5.0 (cm)', 'tar 7.5 (cm)', 'tar 10.0 (cm)'}; % could be different when task diff
% %             T = table(sR_table(:,1), sR_table(:,2), sR_table(:,3), sR_table(:,4), sR_table(:,5), 'VariableNames', VarNames)
        end
        function [time_mean] = getConditionaltime(obj) 
            % todo... form to multi-conditional function
% %             % [trialTime] = getConditionaltime(obj) 
% %             tard = setdiff([obj.tarRs], 0); % potential BUG here if both 0 and others are targets direction!
% %             if isempty(tard)
% %                 tard = 0;
% %             end
% %             tar_num = length(tard);
% %             tarl = obj.tarLs;
% %             tarl_num = length(tarl);
% %             fThs = obj.fThs;
% %             fThs_num = length(fThs);
% %             trialTime = zeros(tar_num, tarl_num, fThs_num);
% %             % copy all the tarR, tarL, fTh from all trials first
% %             tarR = zeros(1, obj.trials_num);
% %             tarL = zeros(1, obj.trials_num);
% %             fTh  = zeros(1, obj.trials_num);
% %             for trial_i = 1:obj.trials_num
% %                 if isempty(obj.trials(trial_i).tarR)
% %                     tarR(trial_i) = -1;
% %                 else
% %                     tarR(trial_i) = obj.trials(trial_i).tarR;
% %                 end
% %                 if isempty(obj.trials(trial_i).tarL)
% %                     tarL(trial_i) = -1;
% %                 else
% %                     tarL(trial_i) = obj.trials(trial_i).tarL;
% %                 end
% %                 if isempty(obj.trials(trial_i).fTh) 
% %                     fTh(trial_i) = -1;
% %                 else
% %                     fTh(trial_i) = obj.trials(trial_i).fTh;
% %                 end
% %             end
% %             for tard_i = 1:tar_num
% %                 for tarl_i = 1:tarl_num
% %                     for fThi = 1:fThs_num
% %                         trialid = ...
% %                             (tarR == tard(tard_i) &...
% %                             tarL == tarl(tarl_i) & ...
% %                             fTh  == fThs(fThi));
% %                         time_all = [obj.trials(trialid).edn_t] - [obj.trials(trialid).bgn_t];
% %                         time_mean(tard_i, tarl_i, fThi) = mean(time_all);
% %                     end
% %                 end
% %             end
% %             time_2d = reshape(time_mean(1,:,:), size(time_mean, 2), size(time_mean, 3));
% %             time_table = [[obj.fThs]', time_2d'];
% %             % display rate using table
% %             display(['For session' num2str(obj.ssnum)]);
% %             VarNames = {'Force (N)', 'tar 2.5 (cm)', 'tar 5.0 (cm)', 'tar 7.5 (cm)', 'tar 10.0 (cm)'}; % could be different when task diff
% %             T = table(time_table(:,1), time_table(:,2), time_table(:,3), time_table(:,4), time_table(:,5), 'VariableNames', VarNames)
        end
%         function obj = highSampleForce(obj, force_obj, ft_intm)
%             % highSampleForce(OBJ,FORCE_OBJ,FT_INTM)
%             % FILL THE FORCE DATA INTO OBJ
% 
%              obj.force_h = force_obj.force;         % rotated
%              if isempty(obj.force_t)
%                 disp('ERROR: no detected aligned force_time. SessionScan::forceHignSample');
%              end
%         end
%         function obj = highSampleWam(obj, wam_obj, wam_intm)
%             % highSampleWam(OBJ,WAM_OBJ,WAM_INTM) 
%             % FILL THE WAM DATA INTO OBJ
%              
%             obj.wamp_h = wam_obj.tp;
%             obj.wamv_h = wam_obj.tv;
%             obj.wamt_h = wam_obj.jt;
%             obj.wam_ts = wam_obj.state;
%             if (isempty(obj.wam_t))
%                 disp('ERROR: no detected aligned wam_t. SessionScan::highSampleWam');
%             end
%         end
%         function obj = highSampleEmg(obj, emg_obj)
%             % obj = highSampleEmg(obj, emg_obj)
%             % align EMG to task time. EMG data comes from a  seperate 
%             % emg.obj file. the seperate emg.obj file should  extract from 
%             % SessionScanEMG 
%             % the after sampled data will still be the same sampling rate
%             % as the EMG file, but the time has being aligned with the task
%             % time.
%             brtime_receieved = obj.Data.SpikeTimestamp;
%             disp([num2str(sum(isnan(brtime_receieved))) ' timepoints were nan']);
%             brtime_receieved(1) = brtime_receieved(2) - mean(diff(brtime_receieved), 'omitnan');
%             brrawtime_aligned = interp1(brtime_receieved, obj.Data.Time, emg_obj.brtime, 'linear', 'extrap');
%             obj.emg_h = emg_obj.dat;
%             obj.emg_t = brrawtime_aligned;
%             
%         end
        function obj = interpData(obj)
            % OBJ = INTERPDATA(OBJ)
            % INTERPOLATE DATA TO ALIGNED DATA IN FORMAT 
            %   CONDITION: wam_t and force_t was aligned to the blackrock
            %   time
            %   TODO: align force_t into wam_t time, interpolate and save
            %   data in the obj.data
%             obj.force_t
%             obj.wam_t
            force_h = interp1(obj.force_t', obj.ft.force', obj.wam_t', 'linear', 'extrap')'; 
            ifplot = 1;
            if (ifplot)
                clf;
                axh(1) = subplot(3,1,1);  hold on;
                plot(obj.force_t, obj.ft.force(1,:), 'r.');
                plot(obj.wam_t, force_h(1,:), 'b.');
                axh(2) = subplot(3,1,2);  hold on;
                plot(obj.force_t, obj.ft.force(2,:), 'r.');
                plot(obj.wam_t, force_h(2,:), 'b.');
                axh(3) = subplot(3,1,3);  hold on;
                plot(obj.force_t, obj.ft.force(3,:), 'r.');
                plot(obj.wam_t, force_h(3,:), 'b.');
                linkaxes(axh, 'x');
            end

            obj.data.t = obj.wam_t;
            obj.data.x = obj.wam.tp;
            obj.data.v = obj.wam.tv;
            obj.data.Fp = obj.wam.cf;
            obj.data.tq = obj.wam.jt;
            obj.data.jp = obj.wam.jp;
            obj.data.ts = obj.wam.state;
            obj.data.f = force_h;
            
            if (ifplot)
                clf;
                axh(1) = subplot(5,1,1);  hold on;
                plot(obj.data.t, obj.data.x(2,:), 'r.'); 
                ylabel('position (m)');
                axh(2) = subplot(5,1,2);  hold on;
                plot(obj.data.t, obj.data.v(2,:), 'b.');
                ylabel('velocity (m/s)');
                axh(3) = subplot(5,1,3);  hold on;
                plot(obj.data.t, obj.data.f(2,:), 'r.');
                ylabel('force (N)');
                axh(4) = subplot(5,1,4);  hold on;
                plot(obj.data.t, obj.data.Fp(2,:), 'r.');
                ylabel('PertForce (N)');
                axh(5) = subplot(5,1,5);  hold on;
                plot(obj.data.t, obj.data.tq(4,:), 'b.');
                ylabel('torque (Nm)');
                linkaxes(axh, 'x');
            end
        end
        function force = rotateAxisForce(obj) % convert from select into world axis
            % force = rotateAxisForce(obj) 
            % This is because the force here is from RTMA message
            force = obj.FTrot_M * obj.Data.Force.Sensor(1:3,:);
        end
%         function [resample_t, resample_f] = trialDataResampleFT(obj, trial_idx, timezone)
%             % [resample_t, resample_f] = trialDataResampleFT(obj, trial_idx)
%             % for all trials resample the original data and time
%             % origin data and time are Nonuniformly sampled
%             % due to the FT is non-normally distributed (code bug of
%             % FTmodule -cg), I used resample to get code data
%             if nargin<2
%                 trial_idx = [obj.trials.outcome]==1;
%             end
%             trial_idx_num = find(trial_idx);
%             trials = (obj.trials(trial_idx));
%             %display(['Enter function Resample;']);
%             if exist('timezone', 'var')
%                 tz_bgn = timezone(1);
%                 tz_edn = timezone(2);
%             else
%                 tz_bgn = -0.5;
%                 tz_edn =  1.0;
%             end
%             resample_freq = 500;   % 500Hz
%             resample_t = [tz_bgn: (1/resample_freq): tz_edn];
%             resample_t = resample_t(2:end); % looks like this one is loger 1 element than Ty? how to deal withit?
%             resample_f = zeros(length(trials), length(resample_t), 3); % resample_value, respectively, x, y, z
% 
%             for trial_i = 1:length(trials) % for all trials
%                 % select specific timezone
%                 idx_t = trials(trial_i).force_t > tz_bgn & trials(trial_i).force_t <= tz_edn;
%                 irregTx = trials(trial_i).force_t(idx_t); % problem here, not wanted force_t
%                 
%                 for dim_i = 1:3 % x, y, z seperately
%                     x = trials(trial_i).force_h(dim_i,idx_t);
%                     try
%                     [y, Ty] = resample(x,irregTx,resample_freq);            % the non-uniform resample
%                     %[y, Ty] = resample(y,irregTx,resample_freq);
%                     catch
%                         y = [];
%                         Ty = [];
%                         display(['Unable to resample trial' num2str(trial_i) ' dim' num2str(dim_i)]);
%                     end
%                     
%                     if (0) % visualize the resample result
%                         figure(); 
%                         hold on;
%                         plot(irregTx,x,'.-', Ty,y,'o-')
%                         legend('Original','Resampled')
%                         plot(Ty);
%                     end
%                     try
%                         resample_f(trial_i,:,dim_i) = y;
%                     catch
%                     %    display(['trial:' num2str(trial_idx_num(trial_i))]);
%                         y = resample(y,size(resample_f,2),length(y));         % may cause time skew here!
%                         resample_f(trial_i,:,dim_i) = y;
%                     end
%                 end
%                 %resample_t = Ty;  % still have 1ms variance between different trials, why?
%             end
%         end
function [t, f_mat] = alignTrialDataF(obj, trial_idx, timezone)
    %  [t, f_mat] = alignTrialDataF(obj, trial_idx, timezone)
            % for all trials resample the original data and time
            % origin data and time are Nonuniformly sampled
            if nargin<2
                trial_idx = find([obj.trials.outcome]==1);
            end
            trial_idx_num = trial_idx;
            trialstmp = (obj.trials(trial_idx_num));
            if exist('timezone', 'var')
                tz_bgn = timezone(1);
                tz_edn = timezone(2);
            else
                tz_bgn = -0.5;
                tz_edn =  0.7;
            end

            for trial_i = 1:length(trialstmp) % for all trials
                % select specific timezone
                idx_t = trialstmp(trial_i).data.t_shift > tz_bgn & trialstmp(trial_i).data.t_shift <= tz_edn;
                if (trial_i) == 1
                    t = trialstmp(trial_i).data.t_shift(idx_t);
                end

                for dim_i = 1:3 % x, y, z seperately
                    f_mat(trial_i,:,dim_i) = trialstmp(trial_i).data.f(dim_i,idx_t);
                end
            end
end

function [t, x_mat, v_mat] = alignTrialDataX(obj, trial_idx, timezone)
        %  [t, x_mat, v_mat] = alignTrialDataX(obj, trial_idx, timezone)
            % for all trials resample the original data and time
            % origin data and time are Nonuniformly sampled
            if nargin<2
                trial_idx = find([obj.trials.outcome]==1);
            end
            trial_idx_num = trial_idx;
            trialstmp = (obj.trials(trial_idx_num));
            if exist('timezone', 'var')
                tz_bgn = timezone(1);
                tz_edn = timezone(2);
            else
                tz_bgn = -0.5;
                tz_edn =  0.7;
            end

            for trial_i = 1:length(trialstmp) % for all trials
                % select specific timezone
                idx_t = trialstmp(trial_i).data.t_shift > tz_bgn & trialstmp(trial_i).data.t_shift <= tz_edn;
                if (trial_i) == 1
                    t = trialstmp(trial_i).data.t_shift(idx_t);
                end

                for dim_i = 1:3 % x, y, z seperately
                    x_mat(trial_i,:,dim_i) = trialstmp(trial_i).data.x(dim_i,idx_t);
                    v_mat(trial_i,:,dim_i) = trialstmp(trial_i).data.v(dim_i,idx_t);
                end
            end
end
%         function [resample_t, resample_p, resample_v] = trialDataAlignWAM(obj, trial_idx, timezone)
%             % for all trials align the original data and time into one
%             % matrix
%             % Assuming they are uniformlly sampled within one trial
%             % pick correct ones and save them in the matrix
%             if nargin<2
%                 trial_idx = [obj.trials.outcome]==1;
%             end
%             if exist('timezone', 'var')
%                 tz_bgn = timezone(1);
%                 tz_edn = timezone(2);
%             else
%                 tz_bgn = -0.5; %time-zone
%                 tz_edn =  0.8;
%             end
%             trial_idx_num = find(trial_idx);
%             trials = (obj.trials(trial_idx));
%             display(['Enter function Resample;']);
% 
%             resample_freq = 500;   % 500Hz, same with WAM
%             resample_t = [tz_bgn: (1/resample_freq): tz_edn];
%             [~,idx_tmp] = sort(abs(resample_t - 0),2,'ascend'); 
%             resample_p = nan(length(trials), length(resample_t), 3); % x, y z
%             resample_v = nan(length(trials), length(resample_t), 3); % x, y z
%             resample_left_num = sum(resample_t<=0);         % left contain zero
%             %resample_left_idx = -resample_left_num+1:0;
%             resample_right_num = sum(resample_t>0);
%             %resample_right_idx = 1:resample_right_num;
%             % resample for all the positions
% 
%             for trial_i = 1:length(trials) % for all trials
%                 wam_t = trials(trial_i).position_t;
%                 % find the cloest to zero 1
%                 [~,idx_Tmp] = sort(abs(wam_t - 0),2,'ascend');
%                 resample_Left_num = idx_Tmp(1); 
%                 resample_Right_num = length(wam_t) - idx_Tmp(1);
%                 % pair index
%                 left_num = min(resample_left_num, resample_Left_num);
%                 right_num = min(resample_right_num, resample_Right_num);
%                 resample_left_idx = -left_num + idx_tmp(1)+1;     % avoid 0
%                 resample_right_idx= right_num + idx_tmp(1);
%                 resample_Left_idx = -left_num + idx_Tmp(1)+1;
%                 resample_Right_idx= right_num + idx_Tmp(1);
%                 % put index in the resampe_p and resample_v
%                 resample_p(trial_i,resample_left_idx:resample_right_idx,1) = ...
%                     trials(trial_i).position_h(1,resample_Left_idx:resample_Right_idx);
%                 resample_p(trial_i,resample_left_idx:resample_right_idx,2) = ...
%                     trials(trial_i).position_h(2,resample_Left_idx:resample_Right_idx);
%                 resample_p(trial_i,resample_left_idx:resample_right_idx,3) = ...
%                     trials(trial_i).position_h(3,resample_Left_idx:resample_Right_idx);
%                 resample_v(trial_i,resample_left_idx:resample_right_idx,1) = ...
%                     trials(trial_i).velocity_h(1,resample_Left_idx:resample_Right_idx);
%                 resample_v(trial_i,resample_left_idx:resample_right_idx,2) = ...
%                     trials(trial_i).velocity_h(2,resample_Left_idx:resample_Right_idx);
%                 resample_v(trial_i,resample_left_idx:resample_right_idx,3) = ...
%                     trials(trial_i).velocity_h(3,resample_Left_idx:resample_Right_idx);
%             end
% 
%         end
%         function [obj, axh_list] = batchPredImpedanceLinDev_4th(obj)
%             axh_list = zeros(size(obj.trials));
%             for trial_i = 1:length(obj.trials)
%                 if obj.trials(trial_i).outcome==1 
%                     try
%                         obj.trials(trial_i) = obj.trials(trial_i).predictImpedanceLinDev();
%                         %axh_list(trial_i) = obj.trials(trial_i).plotPredictedForceOnPosition();
%                     catch
%                         display(['unable to calculate in trial' num2str(obj.trials(trial_i).tNo)]);
%                     end
%                 end
%             end
%         end
%         function [obj, axh_list] = batchPredImpedanceLinDev_5th(obj)
%             axh_list = zeros(size(obj.trials));
%             for trial_i = 1:length(obj.trials)
%                 if obj.trials(trial_i).outcome==1 
%                     try
%                         obj.trials(trial_i) = obj.trials(trial_i).predictImpedanceLinDevS();
%                         %axh_list(trial_i) = obj.trials(trial_i).plotPredictedForceOnPosition();
%                     catch
%                         display(['unable to calculate in trial' num2str(obj.trials(trial_i).tNo)]);
%                     end
%                 end
%             end
%         end
        function pert_idx = getPerturbedTrialIdx(obj, outcome)
            % PERT_IDX = GETPERTURBEDTRIALIDX(OBJ,OUTCOME)
            % return the trials index that has perturbations 
            % default, outcome = 1 (suceed trials);
            % when specify outcome = 0, return unsucessful trials; 
            % when specify outcome = -1, return all perturbed trials 
            if ~exist('outcome', 'var')
                outcome = 1;
            end
                pert_idx = find([obj.trials.ifpert]);     % only depend on ifpert, is it safe?
                pert_idx = setdiff(pert_idx, 1);          % first trial always have bad output, remvoe it
                trials_fin  = find([obj.trials(:).outcome] == outcome);
                if (outcome ~= -1) 
                    pert_idx = intersect(pert_idx, trials_fin);
                end
                % ! no use todo the following part.
%                 empty_idx = zeros(size(obj.trials)); 
%                 for trial_i = 1:length(obj.trials)
%                     empty_idx(trial_i) = isempty(obj.trials(trial_i).pert_t_bgn);
%                 end
%                 pert_idx = setdiff(pert_idx, find(empty_idx));
%                 
            return
        end
%         function yield_positions = getPertPos(obj, filtHz)
%             ifLowpass = 1; % read out the `$low_pass_freq` low pass force
%             if ~exist('filtHz', 'var')
%                 low_pass_freq = 10; % Hz 
%             else
%                 low_pass_freq = filtHz;
%             end
%             % return each trial Position according to the task condition
%             % (if more than 1 task conditions, may cause error).
%             % The position was calculated the average value from epecified
%             % time zone, of the last time during perturbation. (on the last
%             % datapoint in the pulse).
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot position
%             
%             % plot color
%             if ~exist('col_i', 'var')
%                 col_i = 1;
%             end
%             % align for the perturbation time
%             for trial_i = 1:length(obj.trials)
%                 if (obj.trials(trial_i).ifpert)
%                     obj.trials(trial_i) = alignPertInit(obj.trials(trial_i), obj);
%                 end
%             end
%             yield_positions = [];
%             % get the mean
%             % col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%             hold on;
%             % get the perturbation time
%                 % go with the first peturbed trial
%             obj = obj.updatePertEachTrial;
%             trials_idx = obj.getPerturbedTrialIdx();
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             pert_tz = [obj.trials(trials_idx(1)).pert_t_bgn obj.trials(trials_idx(1)).pert_t_edn]; % timezone 
%                                                                             %... This line always get wrong for no reason, but not when execute seperately
%             clearance = 0.2;
%             tz_interest = [pert_tz(1)-clearance, pert_tz(2)+clearance];
%             [resample_t, resample_p, ~] = trialDataAlignWAM(obj, trials_idx,tz_interest);
%             freq = 1/mean(diff(resample_t));
%             respp = zeros(size(resample_p));
%             for axi = 1:size(resample_p,3) % xyz
%                 for ti = 1:size(resample_p,1)
%                     if (ifLowpass)
%                         respp(ti,:,axi) = smooth(resample_p(ti,:,axi), ceil(freq/low_pass_freq));
%                     end
%                 end
%             end
%             if_substeady = 2;
%             if (if_substeady==1)    % relatively error
%                 respp( :,:, 1) = respp(:,:,1) - mean(respp(:,1:50,1), 2, 'omitnan');
%                 respp( :,:, 2) = respp(:,:,2) - mean(respp(:,1:50,2), 2, 'omitnan');
%                 respp( :,:, 3) = respp(:,:,3) - mean(respp(:,1:50,3), 2, 'omitnan');
%             elseif (if_substeady == 2)  % absolute error
%                 respp( :,:, 1) = respp(:,:,1) - 0.482;
%                 respp( :,:, 2) = respp(:,:,2) - 0.482;
%                 respp( :,:, 3) = respp(:,:,3) - 0.482;
%             elseif (if_substeady == 3)  % absolute value
%                 respp( :,:, 1) = respp(:,:,1);
%                 respp( :,:, 2) = respp(:,:,2);
%                 respp( :,:, 3) = respp(:,:,3);
%             end
%             % find stady values of resample_p
%             steadyVal = findSteadyValue(respp(:,:,2), 50, 1);
%             yield_positions = steadyVal;
%         end
%         function yield_force = getPertFce(obj, filtHz)
%             % return each trial censored Force according to the task condition
%             % (if more than 1 task conditions, may cause error).
%             % The position was calculated the average value from epecified
%             % time zone, of the last time during perturbation. (on the last
%             % datapoint in the pulse).
%             ifLowpass = 1; % read out the `$low_pass_freq` low pass force
%             if ~exist('filtHz', 'var')
%                 low_pass_freq = 10; % Hz 
%             else
%                 low_pass_freq = filtHz;
%             end
%           
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot position
%             
%             % plot color
%             if ~exist('col_i', 'var')
%                 col_i = 1;
%             end
%             % align for the perturbation time
%             for trial_i = 1:length(obj.trials)
%                 if (obj.trials(trial_i).ifpert)
%                     obj.trials(trial_i) = alignPertInit(obj.trials(trial_i), obj);
%                 end
%             end
%             yield_force = [];
%             % get the mean
%             % col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%             trials_idx = obj.getPerturbedTrialIdx();
%             % get the perturbation time
%                 % go with the first peturbed trial
%             obj = obj.updatePertEachTrial;
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             pert_tz = [obj.trials(trials_idx(1)).pert_t_bgn obj.trials(trials_idx(1)).pert_t_edn]; % timezone
%             clearance = 0.2;
%             tz_interest = [pert_tz(1)-clearance, pert_tz(2)+clearance];
%             [resample_t, resample_f] = trialDataResampleFT(obj, trials_idx, tz_interest);
%             freq = 1/mean(diff(resample_t));
%             if (ifLowpass)
%                 resample_fF = zeros(size(resample_f));
%                 for axi = 1:size(resample_f, 3)     % xyz 
%                     for ti = 1:size(resample_f,1) % trial_i
%                         resample_fF(ti,:,axi) = smooth(resample_f(ti,:,axi), ceil(freq/low_pass_freq));
%                     end
%                 end
% 
%             end
%             if_substeady = 1;
%             if (if_substeady==1)    % relatively error
%                 resample_fF( :,:, 1) = resample_fF(:,:,1) - mean(resample_fF(:,1:50,1), 2, 'omitnan');
%                 resample_fF( :,:, 2) = resample_fF(:,:,2) - mean(resample_fF(:,1:50,2), 2, 'omitnan');
%                 resample_fF( :,:, 3) = resample_fF(:,:,3) - mean(resample_fF(:,1:50,3), 2, 'omitnan');
%             elseif (if_substeady == 2)  % absolute error
%                 resample_fF( :,:, 1) = resample_fF(:,:,1) - 0.482;
%                 resample_fF( :,:, 2) = resample_fF(:,:,2) - 0.482;
%                 resample_fF( :,:, 3) = resample_fF(:,:,3) - 0.482;
%             end
%             % low-pass filter the fce
%             % find stady values of resample_p
%             steadyVal = findSteadyValue(resample_fF(:,:,2), 200, 1);
%             yield_force = steadyVal;
%         end
%         function command_force = getPertFce_cmd(obj)
%             % extract command force directly from the session data
%             for trial_i = 1:length(obj.trials)
%                 if (obj.trials(trial_i).ifpert)
%                     obj.trials(trial_i) = alignPertInit(obj.trials(trial_i), obj);
%                 end
%             end
%             
%             trials_idx = obj.getPerturbedTrialIdx();
%             command_force = (zeros(size(trials_idx)))';
%             % get the perturbation time
%                 % go with the first peturbed trial
%             pert_time_idx = find(obj.trials(trials_idx(1)).pertfce_h ~= 0);
%             for trial_i = 1:length(trials_idx)
%                 try
%                     command_force(trial_i) = setdiff(unique(obj.trials(trials_idx(trial_i)).pertfce_h(pert_time_idx)), 0);
%                 catch 
%                     disp(['no pertfce on trial' num2str(trials_idx(trial_i)) ', give a 0!']);
%                     command_force(trial_i) = 0;
%                 end
%             end
%         end
        %%% other process
        function obj_new = ConcatTrials(obj1, obj2, trial_idx1, trial_idx2)
            trials = [obj1.trials(trial_idx1) obj2.trials(trial_idx2)];
            obj_new = obj1;
            obj_new.trials = trials;
        end
%         function obj = processSession_n(obj)
%             if obj.ssnum == 1965
%                 % all trial 
%                 for trial_i = 1:length(obj.trials)
%                     % set fTh = 15;
%                     obj.trials(trial_i).fTh = 15;
%                     % all trial set force_y+5.4;
%                     obj.trials(trial_i).force(2,:) = obj.trials(trial_i).force(2,:)+5.4;
%                     obj.trials(trial_i).force_h(2,:) = obj.trials(trial_i).force_h(2,:)+5.4;
%                 end
%             end
%         end
%         function obj = updatePertEachTrial(obj)
%             % obj = updatePertEachTrial()
%             % Updating if_pert variable depend on whether perturbed or not,
%             % because some time a trial should be perturbed, however as the
%             % force threshold never achieved enough long, the perturbed was
%             % not achieved. 
%             % Use function TrialScan.findStepPerterbTime()
%             
%             %trial_num = length(obj.trials);
%             
%             trial_idx = obj.getPerturbedTrialIdx();
%             for trial_i = trial_idx
%                 if obj.trials(trial_i).ifpert && (obj.trials(trial_i).outcome==1)
%                     switch obj.trials(trial_i).getPerturbationPattern
%                         case 2
%                             obj.trials(trial_i) = obj.trials(trial_i).findStepPerterbTime();
%                         case 3
%                         obj.trials(trial_i) = obj.trials(trial_i).findStepx0PerterbTime();
%                     end
%                 else
%                     obj.trials(trial_i).pert_t_bgn = nan; % will this cause error?
%                 end
%             end
%         end
        %%% communicate 
%         function obj = generateWamPertData(obj) % old function for ensemble analysis
%             % send data into wam function to help SessionScanWam generate perturbation-only data
%             % perturbation rdt already saved in TrialScan
%             tarL_list = [];
%             fTh_list = [];
%             rdt_ranges_all = {};
%             i = 0;
%             for tarLi = obj.tarLs
%                 for fThi = obj.fThs
%                     i = i+1;
%                     trials_idx = [obj.trials.tarL] == tarLi &...
%                                     [obj.trials.fTh] == fThi &...
%                                     [obj.trials.outcome] == 1;      % only sucessful trials
%                     rdt_ranges = [obj.trials(trials_idx).pert_rdt_bgn;...  
%                                     obj.trials(trials_idx).pert_rdt_edn];
%                     tarL_list = [tarL_list tarLi];
%                     fTh_list  = [fTh_list fThi];
%                     rdt_ranges_all{i} = rdt_ranges;
%                 end
%             end
%             obj.wam = obj.wam.concatinateTrials2File(tarL_list, fTh_list, rdt_ranges_all);
%             % for each trial condition, concatinate a structure
%             
%              % call generateWamPertData()
%         end
        
        function [obj] = add_to_PertData(obj, new_Data_pert)
            n_test1 = length(obj.wam.Data_pert);
            n_test2 = length(new_Data_pert);

            for i = 1:n_test2
                obj.wam.Data_pert(i+n_test1) = new_Data_pert(i);
            end
            
        end
        
%         function obj = generateWamPertData_ensemble(obj)
%             % send data into wam function to help SessionScanWam generate perturbation-only data
%             % perturbation rdt already saved in TrialScan
%             
%             % Find and talk out start up time in zero state
%             dexStart = max(find(obj.wam.state == 0))+1;
% %             dexEnd = max(find(abs(diff(obj.wam.state(1:499999))) == 6)); % Temporarly cut end due to change in sampling rate
%             dexEnd = max(find(abs(diff(obj.wam.state)) == 6)); % Keep all
% 
%             state = obj.wam.state(dexStart:dexEnd);
%             
%             % Figure out state
%             success = zeros(size(state));
%             last_dexForceRamp = 0;
%             for i = 1:length(state)-1
%                 % Set last to 1 if its the start of the force hold
%                 if(state(i) == 2 && state(i+1) == 3)
%                     last_dexForceRamp = i;
%                 end
%                 
%                 % Set state last to zero if trial progreses to fast
%                 if((state(i+1)-state(i)~=0))
%                     if((state(i+1)-state(i))~=1)
%                         % leave success at zero
%                         last_dexForceRamp = 0;
%                     end
%                 end
%                 
%                 % If both checks are passed and 6 is reached consider it a
%                 % success
%                 if(last_dexForceRamp ~= 0 && state(i+1) == 6)
%                     success(last_dexForceRamp:i) = 1;
%                 end
%             end
%             
%             % Find trial onset
%             dexForceRampStart = find(state.*success == 3 & [0;abs(diff(state))>=1])+dexStart;
%             dexMovOnset = find(state.*success == 4 & [0;abs(diff(state))>=1])+dexStart-1;
%             dexHoldEnd = find(state.*success == 6 & [0;abs(diff(state))>=1])+dexStart-2;
%             
%             state_advance = zeros(size(state));
%             state_advance(dexForceRampStart) = 1;
%             state_advance(dexMovOnset) = 2;
%             state_advance(dexHoldEnd) = 3;
%             
%                             
%             % Chop first trial
%             dexForceRampStart = dexForceRampStart(2:end-3);
%             dexMovOnset = dexMovOnset(2:end-3);
%             dexHoldEnd = dexHoldEnd(2:end-3);
%             
% %             % Check if the size of the dexs are the same. If not check
% %             % order
% %             if(0~=sum(diff([length(dexForceRampStart),length(dexMovOnset),length(dexHoldEnd)])))
% %                 if(dexMovOnset(1) > dexForceRampStart(1)) % remove extra trial
% %                     dexMovOnset = dexMovOnset(2:end);
% %                     dexHoldEnd = dexHoldEnd(2:end);
% %                 end
% %                 
% %                 if(length(dexHoldEnd) > length(dexForceRampStart)+100) % If weird problem
% %                     dexHoldEnd = dexHoldEnd(1:length(dexMovOnset));
% %                     dexForceRampStart = dexForceRampStart(1:length(dexMovOnset));
% %                 end
% %             end
%             
% %             % Check for failed trials that did not go through all states
% %             % If length is all the same dont bother checking
% %             if ~(length(dexHoldEnd) == length(dexMovOnset))
% %                 dexSkip = find(dexHoldEnd(1:length(dexMovOnset)) - dexMovOnset<0,1); 
% %                 dexForceRampStart = dexForceRampStart([1:dexSkip-1,dexSkip+1:length(dexForceRampStart)]);
% %                 dexHoldEnd = dexHoldEnd([1:dexSkip-1,dexSkip+1:length(dexHoldEnd)]);
% %             end
% 
%             
%             % Check the ranges were chosen correctly
%             figure; 
%             ax1 = subplot(3,1,1);
%             plot(obj.wam.time,obj.wam.tp(:,2),'linewidth',2.5); hold on;
%             plot(obj.wam.time(dexForceRampStart),obj.wam.tp(dexForceRampStart(1), 2),'o','linewidth',2.5); hold on;
%             plot(obj.wam.time(dexMovOnset),obj.wam.tp(dexMovOnset, 2),'*','linewidth',2.5); grid on;
%             plot(obj.wam.time(dexHoldEnd),obj.wam.tp(dexHoldEnd, 2),'+','linewidth',2.5); grid on;
%             ylabel('Postion');
%             
%             ax2 = subplot(3,1,2);
%             plot(obj.wam.time,obj.wam.state,'linewidth',2.5); hold on; grid on; ylim([0 8]);
% %             plot(obj.wam.time(dexStart:dexEnd),success);
%             plot(obj.wam.time(dexForceRampStart),obj.wam.state(dexForceRampStart),'o','linewidth',2.5); hold on;
%             plot(obj.wam.time(dexMovOnset),obj.wam.state(dexMovOnset),'*','linewidth',2.5); grid on;
%             plot(obj.wam.time(dexHoldEnd),obj.wam.state(dexHoldEnd),'+','linewidth',2.5); grid on;
%             ylabel('State');
%             
%             ax3 = subplot(3,1,3);
%             plot(obj.wam.time(1:end-1),1./diff(obj.wam.time),'.','linewidth',2.5);
%             ylabel('Sampling rate (samples/s)'); xlabel('Time (s)');
%             
%             linkaxes([ax1, ax2, ax3],'x');
%         
%             tarL_list = [];
%             fTh_list = [];
%             rdt_ranges_all = {};
%             i = 0;
%             for tarLi = obj.tarLs
%                 for fThi = obj.fThs
%                     i = i+1;
%                     % Remove trials_idx for now add to select successful
%                     % and right contion in the onset time (e.g. obj.trials(trials_idx).idx_fcr)
% %                     trials_idx = [obj.trials.tarL] == tarLi &...
% %                                     [obj.trials.fTh] == fThi;%; &...
% %                                     [obj.trials.outcome] == 1;      % only sucessful trials
% %                                 
% %                     % Find movement onsent time and align them            
% %                     rdt_tmp = [obj.trials.idx_fcr;...
% %                                     obj.trials.idx_mov;...
% %                                     obj.trials.idx_end];
% %                     rdt_ranges = [obj.trials(trials_idx).bgn] + ...
% %                         [[obj.trials(trials_idx).idx_mov] -  min(rangeDiff(1,:));...
% %                         [obj.trials(trials_idx).idx_mov] +  min(rangeDiff(2,:))];
% %                     rdt_mov = [obj.trials(trials_idx).bgn] + [obj.trials(trials_idx).idx_mov]; % Movement onset use to alter time vector
% 
%                     rangeDiff = diff([dexForceRampStart';...
%                                dexMovOnset';...
%                                dexHoldEnd']);
%                              
%                     ranges = [[dexMovOnset'] -  min(rangeDiff(1,:));...
%                               [dexMovOnset'] +  min(rangeDiff(2,:))];
%                     tarL_list = [tarL_list tarLi];
%                     fTh_list  = [fTh_list fThi];
%                     rdt_ranges_all{i} = ranges;
%                     dexMovOnset_all{i} = dexMovOnset;
%                 end
%             end
%             obj.wam = obj.wam.concatinateTrials2File_ensemble(tarL_list, fTh_list, rdt_ranges_all, dexMovOnset_all);
%             % for each trial condition, concatinate a structure
%             % call generateWamPertData_ensemble()
%             
%         end
        
%         function [obj] = add_to_ensemble(obj, new_Data_pert_ensemble)
%             
%             for i = 1:length(new_Data_pert_ensemble)
%                 n = length(obj.wam.Data_pert_ensemble(i).time_r) - length(new_Data_pert_ensemble(i).time_r);
%                 
%                 obj.wam.Data_pert_ensemble(i).time_r = [obj.wam.Data_pert_ensemble(i).time_r;...
%                                          obj.get_chopOrPadd(new_Data_pert_ensemble(i).time_r,n)];
%                                      
%                 obj.wam.Data_pert_ensemble(i).z_r_1 = [obj.wam.Data_pert_ensemble(i).z_r_1;...
%                                         obj.get_chopOrPadd(new_Data_pert_ensemble(i).z_r_1,n)];
%                 obj.wam.Data_pert_ensemble(i).z_r_2 = [obj.wam.Data_pert_ensemble(i).z_r_2;...
%                                         obj.get_chopOrPadd(new_Data_pert_ensemble(i).z_r_2,n)];  
%                                     
%                 obj.wam.Data_pert_ensemble(i).u_r_1 = [obj.wam.Data_pert_ensemble(i).u_r_1;...
%                                         obj.get_chopOrPadd(new_Data_pert_ensemble(i).u_r_1,n)];
%                 obj.wam.Data_pert_ensemble(i).u_r_2 = [obj.wam.Data_pert_ensemble(i).u_r_2;...
%                                         obj.get_chopOrPadd(new_Data_pert_ensemble(i).u_r_2,n)];  
%                                     
%                 obj.wam.Data_pert_ensemble(i).v_r_1 = [obj.wam.Data_pert_ensemble(i).v_r_1;...
%                                         obj.get_chopOrPadd(new_Data_pert_ensemble(i).v_r_1,n)];          
%                 obj.wam.Data_pert_ensemble(i).v_r_2 = [obj.wam.Data_pert_ensemble(i).v_r_2;...
%                                         obj.get_chopOrPadd(new_Data_pert_ensemble(i).v_r_2,n)];  
%                                     
%                 obj.wam.Data_pert_ensemble(i).it_r = [obj.wam.Data_pert_ensemble(i).it_r;...
%                                        obj.get_chopOrPadd(new_Data_pert_ensemble(i).it_r,n)]; 
%                                    
%                 obj.wam.Data_pert_ensemble(i).rdt_r = [obj.wam.Data_pert_ensemble(i).rdt_r;...
%                                         obj.get_chopOrPadd(new_Data_pert_ensemble(i).rdt_r,n)]; 
%                                     
%                 obj.wam.Data_pert_ensemble(i).state_r = [obj.wam.Data_pert_ensemble(i).state_r;...
%                                           obj.get_chopOrPadd(new_Data_pert_ensemble(i).state_r,n)]; 
%             end
% 
%         % Old
%             % R1 = size(ss1972.wam.Data_pert_ensemble,2);
%             % R2 = size(ss1973.wam.Data_pert_ensemble,2);
%             % n =  size(ss1973.wam.Data_pert_ensemble(1).tp,1) - size(ss1972.wam.Data_pert_ensemble(54).tp,1);
%             %
%             % for i = 1:R2
%             %     ss1972.wam.Data_pert_ensemble(R1+i).FT = 10;
%             %     ss1972.wam.Data_pert_ensemble(R1+i).x0 = 0.05;
%             %
%             %     ss1972.wam.Data_pert_ensemble(R1+i).time = [ss1973.wam.Data_pert_ensemble(i).time; zeros(n,1)];
%             %     ss1972.wam.Data_pert_ensemble(R1+i).tp = [ss1973.wam.Data_pert_ensemble(i).tp; zeros(n,3)];
%             %     ss1972.wam.Data_pert_ensemble(R1+i).tv = [ss1973.wam.Data_pert_ensemble(i).tv; zeros(n,3)];
%             %     ss1972.wam.Data_pert_ensemble(R1+i).cf = [ss1973.wam.Data_pert_ensemble(i).cf; zeros(n,3)];
%             %     ss1972.wam.Data_pert_ensemble(R1+i).it = [ss1973.wam.Data_pert_ensemble(i).it; zeros(n,1)];
%             %     ss1972.wam.Data_pert_ensemble(R1+i).rdt = [ss1973.wam.Data_pert_ensemble(i).rdt; zeros(n,1)];
%             %
%             % end
%             
%         end
        
%         function [mat] = get_chopOrPadd(obj,mat,n)
%             if(n==0)
%                 mat = mat;
%             elseif(n > 0)
%                 mat = [zeros(size(mat,1),n),mat];
%             elseif(n < 0) 
%                 mat = mat(:,1:end+n);
%             end  
%         end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EXPORT CODE FOR FURTHER ANALYSIS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [cellsmat] = export_as_formatted(obj, ifplot)
            % ...TODO: RE-EDIT IT TO GET BETTER OUTPUT
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            % export as a t(trials_num)-by-p(perturbation options) cell mat
            % for each cell, the data format are each trial, which contains:
            %   x: 3-by-N matrix, robot endpoint
            %   v: 3-by-N matrix, robot velocity
            %   f: 3-by-N matrix, force transducer force
            %   Fpert: 1-by-N matrix, perturbation force, in stoc 2-by-N(xy)
            %   ts: 1-by-N matrix, task states
            %   time: 1-by-N matrix, time 
            %   movement onset: the mask that robot start move
            %  -[ ] emg: 8-by-N matrix, emg data
            pert_trials = [obj.trials.ifpert];
            t_idx = cell(1,3);
            for p_i = 1:3
                 t_idx{p_i} = find(pert_trials==p_i-1); % 0,nopert; 1, pulse; 2, stoc
            end
            

            if ~isempty(t_idx{3}) % stoc-perturbed trials. 
                % Assume only stocpert do not in the same session with step ones
                cellsmat = cell(length(obj.tarLs), max(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),15),3);
                for tl_i = 1:length(obj.tarLs)
                    %trial_list = setdiff(find([obj.trials.tarL] == obj.tarLs(tl_i)),1);
                    trial_list = find([obj.trials.tarL] == obj.tarLs(tl_i) & [obj.trials.outcome] == 1);
                    trial_list = setdiff(trial_list,1);
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{tl_i,t_i,3} = t_tmp.export_as_formatted;
                    end
                end
            else
                cellsmat = cell(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),3);
                for p_i = 1:2
                    %trial_list = setdiff(t_idx{p_i},1);
                    trial_list = find([obj.trials.outcome] == 1); % no failed trials 
                    %trial_list = find([obj.trials.outcome] ~= -1); % spring test parameter selection
                    trial_list = intersect(trial_list, t_idx{p_i});
                    %%%%%%%%%%%%%%%%%% add some exceptions here %%%%%%%%%%%%%%%%%%
                    switch obj.ssnum
                        case 3402
                            trial_list = setdiff(trial_list,18);
                        case 3385
                            trial_list = setdiff(trial_list,28);
                        case 3387
                            trial_list = setdiff(trial_list,9);
                        case 3361
                            trial_list = setdiff(trial_list,9);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %trial_list = setdiff(trial_list,[1 2 51:60]);
                    trial_list = setdiff(trial_list,[1]);
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{t_i,p_i} = t_tmp.export_as_formatted;  % each trial
                        xlim([-5 -4])
                        ifplot = true;
                        if (ifplot) 
                        subplot(2,1,1);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.x(2,:));
                        subplot(2,1,2);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.f(2,:));
                        end 
                    end
                end
            end
            
            % plot out the time skew
            ifplot = 0; %-test
            
            if (ifplot) 
                plt_offset = 2e-3/10;
                for p_i = 1:3
                    subplot(1,3,p_i); hold on;
                    switch p_i
                        case 1
                            title('no pert');
                        case 2
                            title('pulse pert');
                        case 3
                            title('stoc pert');
                    end
                    if (p_i ~= 3)
                        for t_i = 1:length(t_idx{p_i})
                            if(isempty(cellsmat{t_i, p_i}))
                                continue;
                            end
                            dt = diff(cellsmat{t_i,p_i}.t);
                            dt = [dt(1) dt];
                            plot(t_i*plt_offset + dt);
                        end
                    else
                        if isempty(t_idx{p_i}) 
                            continue
                        end
                        tiofst = 0; % plot offset
                        for tl_i = 1:size(cellsmat, 1)
                            for t_i = 1:length(cellsmat(tl_i,:,p_i))
                                if ~isempty(cellsmat{tl_i,t_i,p_i})
                                dt = diff(cellsmat{tl_i,t_i,p_i}.t);
                                dt = [dt(1) dt];
                                tiofst = tiofst+1;
                                plot(tiofst*plt_offset + dt);
                                end
                            end
                        end
                    end
                end
            end
        end
        function [cellsmat] = export_as_formatted_4(obj, ifplot)
            % THIS FUNCTION IS DESIGNED FOR EXPORTING TRIALS WITH step perturbation
            % THE EACH EXPORTED DATA COLUMN IS: 
            %   1. NO PULSE;
            %   2. PULSE;
            %   3. STOC;
            %   4. STEP;
            
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            % export as a t(trials_num)-by-p(perturbation options) cell mat
            % for each cell, the data format are each trial, which contains:
            %   x: 3-by-N matrix, robot endpoint
            %   v: 3-by-N matrix, robot velocity
            %   f: 3-by-N matrix, force transducer force
            %   Fpert: 1-by-N matrix, perturbation force, in stoc 2-by-N(xy)
            %   ts: 1-by-N matrix, task states
            %   time: 1-by-N matrix, time 
            %   movement onset: the mask that robot start move
            %  -[ ] emg: 8-by-N matrix, emg data
            pert_trials = [obj.trials.ifpert];
            t_idx = cell(1,4);
            ifpert_keys = [0 1 2 4];
            for p_i = 1:4
                 t_idx{p_i} = find(pert_trials==ifpert_keys(p_i)); % 0,nopert; 1, pulse; 2, stoc; 3, step; 4. square
            end
            

            if ~isempty(t_idx{3}) % stoc-perturbed trials. 
                % Assume only stocpert do not in the same session with step ones
                cellsmat = cell(length(obj.tarLs), max(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),15),3);
                for tl_i = 1:length(obj.tarLs)
                    %trial_list = setdiff(find([obj.trials.tarL] == obj.tarLs(tl_i)),1);
                    trial_list = find([obj.trials.tarL] == obj.tarLs(tl_i) & [obj.trials.outcome] == 1);
                    trial_list = setdiff(trial_list,1);
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{tl_i,t_i,3} = t_tmp.export_as_formatted;
                    end
                end
            else
                cellsmat = cell(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3}), length(t_idx{4})]),3);
                for p_i = [1:2 4]
                    %trial_list = setdiff(t_idx{p_i},1);
                    trial_list = find([obj.trials.outcome] == 1); % no failed trials 
                    %trial_list = find([obj.trials.outcome] ~= -1); % spring test parameter selection
                    trial_list = intersect(trial_list, t_idx{p_i});
                    %%%%%%%%%%%%%%%%%% add some exceptions here %%%%%%%%%%%%%%%%%%
                    switch obj.ssnum
                        case 3402
                            trial_list = setdiff(trial_list,18);
                        case 3385
                            trial_list = setdiff(trial_list,28);
                        case 3387
                            trial_list = setdiff(trial_list,9);
                        case 3361
                            trial_list = setdiff(trial_list,9);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %trial_list = setdiff(trial_list,[1 2 51:60]);
                    trial_list = setdiff(trial_list,[1]);
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{t_i,p_i} = t_tmp.export_as_formatted;  % each trial
                        xlim([-5 -4])
                        ifplot = true;
                        if (ifplot) 
                        subplot(2,1,1);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.x(2,:));
                        subplot(2,1,2);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.f(2,:));
                        end 
                    end
                end
            end
            
            % plot out the time skew
            ifplot = 0; %-test
            
            if (ifplot) 
                plt_offset = 2e-3/10;
                for p_i = 1:3
                    subplot(1,3,p_i); hold on;
                    switch p_i
                        case 1
                            title('no pert');
                        case 2
                            title('pulse pert');
                        case 3
                            title('stoc pert');
                    end
                    if (p_i ~= 3)
                        for t_i = 1:length(t_idx{p_i})
                            if(isempty(cellsmat{t_i, p_i}))
                                continue;
                            end
                            dt = diff(cellsmat{t_i,p_i}.t);
                            dt = [dt(1) dt];
                            plot(t_i*plt_offset + dt);
                        end
                    else
                        if isempty(t_idx{p_i}) 
                            continue
                        end
                        tiofst = 0; % plot offset
                        for tl_i = 1:size(cellsmat, 1)
                            for t_i = 1:length(cellsmat(tl_i,:,p_i))
                                if ~isempty(cellsmat{tl_i,t_i,p_i})
                                dt = diff(cellsmat{tl_i,t_i,p_i}.t);
                                dt = [dt(1) dt];
                                tiofst = tiofst+1;
                                plot(tiofst*plt_offset + dt);
                                end
                            end
                        end
                    end
                end
            end
        end
        function [cellsmat] = export_as_formatted_5(obj, ifplot)
            % THIS FUNCTION IS DESIGNED FOR EXPORTING TRIALS WITH NO
            % RELEASE
            % THE EACH EXPORTED DATA COLUMN IS: 
            %   1. NO PULSE;
            %   2. PULSE;
            %   3. STOC;
            %   4. NO RELEASE;
            
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            % export as a t(trials_num)-by-p(perturbation options) cell mat
            % for each cell, the data format are each trial, which contains:
            %   x: 3-by-N matrix, robot endpoint
            %   v: 3-by-N matrix, robot velocity
            %   f: 3-by-N matrix, force transducer force
            %   Fpert: 1-by-N matrix, perturbation force, in stoc 2-by-N(xy)
            %   ts: 1-by-N matrix, task states
            %   time: 1-by-N matrix, time 
            %   movement onset: the mask that robot start move
            %  -[ ] emg: 8-by-N matrix, emg data
            pert_trials = [obj.trials.ifpert];
            t_idx = cell(1,4);
            ifpert_keys = [0 1 2 5];
            for p_i = 1:4
                 t_idx{p_i} = find(pert_trials==ifpert_keys(p_i)); % 0,nopert; 1, pulse; 2, stoc; 3, step; 5. no release
            end

            if ~isempty(t_idx{3}) % stoc-perturbed trials. 
                % Assume only stocpert do not in the same session with step ones
                cellsmat = cell(length(obj.tarLs), max(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),15),3);
                for tl_i = 1:length(obj.tarLs)
                    %trial_list = setdiff(find([obj.trials.tarL] == obj.tarLs(tl_i)),1);
                    trial_list = find([obj.trials.tarL] == obj.tarLs(tl_i) & [obj.trials.outcome] == 1);
                    trial_list = setdiff(trial_list,1);
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{tl_i,t_i,3} = t_tmp.export_as_formatted;
                    end
                end
            else
                cellsmat = cell(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3}), length(t_idx{4})]),3);
                for p_i = [1:2 4]
                    %trial_list = setdiff(t_idx{p_i},1);
                    trial_list = find([obj.trials.outcome] == 1); % only select sucessful trials 
                    % the session-specific cases here
                    if (sum(obj.ssnum == [3602 3603]))
                        trial_list = find([obj.trials.outcome] ~= -1);
                    end
                    %trial_list = find([obj.trials.outcome] ~= -1); % spring test parameter selection
                    trial_list = intersect(trial_list, t_idx{p_i});
                    %%%%%%%%%%%%%%%%%% add some exceptions here %%%%%%%%%%%%%%%%%%
                    switch obj.ssnum
                        case 3402
                            trial_list = setdiff(trial_list,18);
                        case 3385
                            trial_list = setdiff(trial_list,28);
                        case 3387
                            trial_list = setdiff(trial_list,9);
                        case 3361
                            trial_list = setdiff(trial_list,9);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %trial_list = setdiff(trial_list,[1 2 51:60]);
                    if (sum(obj.ssnum ~= [3603 3602]))
                        trial_list = setdiff(trial_list,[1]);
                    end
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{t_i,p_i} = t_tmp.export_as_formatted;  % each trial
                        xlim([-5 -4])
                        ifplot = true;
                        if (ifplot) 
                        subplot(2,1,1);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.x(2,:));
                        subplot(2,1,2);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.f(2,:));
                        end 
                    end
                end
            end
            
            % plot out the time skew
            ifplot = 0; %-test
            
            if (ifplot) 
                plt_offset = 2e-3/10;
                for p_i = 1:3
                    subplot(1,3,p_i); hold on;
                    switch p_i
                        case 1
                            title('no pert');
                        case 2
                            title('pulse pert');
                        case 3
                            title('stoc pert');
                    end
                    if (p_i ~= 3)
                        for t_i = 1:length(t_idx{p_i})
                            if(isempty(cellsmat{t_i, p_i}))
                                continue;
                            end
                            dt = diff(cellsmat{t_i,p_i}.t);
                            dt = [dt(1) dt];
                            plot(t_i*plt_offset + dt);
                        end
                    else
                        if isempty(t_idx{p_i}) 
                            continue
                        end
                        tiofst = 0; % plot offset
                        for tl_i = 1:size(cellsmat, 1)
                            for t_i = 1:length(cellsmat(tl_i,:,p_i))
                                if ~isempty(cellsmat{tl_i,t_i,p_i})
                                dt = diff(cellsmat{tl_i,t_i,p_i}.t);
                                dt = [dt(1) dt];
                                tiofst = tiofst+1;
                                plot(tiofst*plt_offset + dt);
                                end
                            end
                        end
                    end
                end
            end
        end
        function [cellsmat] = export_as_formatted_5_failedTrials(obj, ifplot)
            % This is work only for the failed trials
            % If the failed trial is less than 15, repeat it to be 15... 
            % THIS FUNCTION IS DESIGNED FOR EXPORTING TRIALS WITH NO
            % RELEASE
            % THE EACH EXPORTED DATA COLUMN IS: 
            %   1. NO PULSE;
            %   2. PULSE;
            %   3. STOC;
            %   4. NO RELEASE;
            
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            % export as a t(trials_num)-by-p(perturbation options) cell mat
            % for each cell, the data format are each trial, which contains:
            %   x: 3-by-N matrix, robot endpoint
            %   v: 3-by-N matrix, robot velocity
            %   f: 3-by-N matrix, force transducer force
            %   Fpert: 1-by-N matrix, perturbation force, in stoc 2-by-N(xy)
            %   ts: 1-by-N matrix, task states
            %   time: 1-by-N matrix, time 
            %   movement onset: the mask that robot start move
            %  -[ ] emg: 8-by-N matrix, emg data
            pert_trials = [obj.trials.ifpert];
            t_idx = cell(1,4);
            ifpert_keys = [0 1 2 5];
            for p_i = 1:4
                 t_idx{p_i} = find(pert_trials==ifpert_keys(p_i)); % 0,nopert; 1, pulse; 2, stoc; 3, step; 5. no release
            end

            if ~isempty(t_idx{3}) % stoc-perturbed trials. 
                % Assume only stocpert do not in the same session with step ones
                cellsmat = cell(length(obj.tarLs), max(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),15),3);
                for tl_i = 1:length(obj.tarLs)
                    %trial_list = setdiff(find([obj.trials.tarL] == obj.tarLs(tl_i)),1);
                    trial_list = find([obj.trials.tarL] == obj.tarLs(tl_i) & [obj.trials.outcome] == 1);
                    trial_list = setdiff(trial_list,1);
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{tl_i,t_i,3} = t_tmp.export_as_formatted;
                    end
                end
            else
                cellsmat = cell(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3}), length(t_idx{4})]),3);
                for p_i = [1:2 4]
                    %trial_list = setdiff(t_idx{p_i},1);
                    %trial_list = find([obj.trials.outcome] == 1); % only select sucessful trials 
                    trial_list = find([obj.trials.outcome] == 0); % only select failed trials 
                    % the session-specific cases here
                    if (sum(obj.ssnum == [3602 3603]))
                        trial_list = find([obj.trials.outcome] ~= -1);
                    end
                    %trial_list = find([obj.trials.outcome] ~= -1); % spring test parameter selection
                    trial_list = intersect(trial_list, t_idx{p_i});
                    %%%%%%%%%%%%%%%%%% add some exceptions here %%%%%%%%%%%%%%%%%%
                    switch obj.ssnum
                        case 3402
                            trial_list = setdiff(trial_list,18);
                        case 3385
                            trial_list = setdiff(trial_list,28);
                        case 3387
                            trial_list = setdiff(trial_list,9);
                        case 3361
                            trial_list = setdiff(trial_list,9);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %trial_list = setdiff(trial_list,[1 2 51:60]);
                    if (sum(obj.ssnum ~= [3603 3602]))
                        trial_list = setdiff(trial_list,[1]);
                    end
                    for t_i = 1:length(trial_list)
                        t_tmp = obj.trials(trial_list(t_i));
                        cellsmat{t_i,p_i} = t_tmp.export_as_formatted;  % each trial
                        xlim([-5 -4])
                        ifplot = true;
                        if (ifplot) 
                        subplot(2,1,1);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.x(2,:));
                        subplot(2,1,2);
                        plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.f(2,:));
                        end 
                    end
                end
                
                %%%% repeat the last in conidtion (failed trials)
                for trial_i = 1:15
                    if isempty(cellsmat{trial_i,2})
                        cellsmat(trial_i,2) = cellsmat(trial_i-1,2);
                    end
                end
            end
            
            % plot out the time skew
            ifplot = 0; %-test
            
            if (ifplot) 
                plt_offset = 2e-3/10;
                for p_i = 1:3
                    subplot(1,3,p_i); hold on;
                    switch p_i
                        case 1
                            title('no pert');
                        case 2
                            title('pulse pert');
                        case 3
                            title('stoc pert');
                    end
                    if (p_i ~= 3)
                        for t_i = 1:length(t_idx{p_i})
                            if(isempty(cellsmat{t_i, p_i}))
                                continue;
                            end
                            dt = diff(cellsmat{t_i,p_i}.t);
                            dt = [dt(1) dt];
                            plot(t_i*plt_offset + dt);
                        end
                    else
                        if isempty(t_idx{p_i}) 
                            continue
                        end
                        tiofst = 0; % plot offset
                        for tl_i = 1:size(cellsmat, 1)
                            for t_i = 1:length(cellsmat(tl_i,:,p_i))
                                if ~isempty(cellsmat{tl_i,t_i,p_i})
                                dt = diff(cellsmat{tl_i,t_i,p_i}.t);
                                dt = [dt(1) dt];
                                tiofst = tiofst+1;
                                plot(tiofst*plt_offset + dt);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function [cellsmat] = export_as_formatted_hybridss(obj, ifplot)
            % hybrids can have pulse, stoc and no pulse pert
            if (~exist('ifplot', 'var'))
                ifplot = 0;
            end
            % export as a t(trials_num)-by-p(perturbation options) cell mat
            % for each cell, the data format are each trial, which contains:
            %   x: 3-by-N matrix, robot endpoint
            %   v: 3-by-N matrix, robot velocity
            %   f: 3-by-N matrix, force transducer force
            %   Fpert: 1-by-N matrix, perturbation force, in stoc 2-by-N(xy)
            %   ts: 1-by-N matrix, task states
            %   time: 1-by-N matrix, time 
            %   movement onset: the mask that robot start move
            %  -[ ] emg: 8-by-N matrix, emg data
            pert_trials = [obj.trials.ifpert];
            t_idx = cell(1,3);
            for p_i = 1:3
                 t_idx{p_i} = find(pert_trials==p_i-1); % 0,nopert; 1, pulse; 2, stoc
            end
            
            t_idx{3} = t_idx{3}(2:end); % works for the hybrid pert, to avoid error
            
            cellsmat = cell(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),3);
            if ~isempty(t_idx{3}) % stoc-perturbed trials. 
                % stoc trials are in the same sessions for step trials (for
                % spring testing)
                for tl_i = 1:length(obj.tarLs)
                    %trial_list = setdiff(find([obj.trials.tarL] == obj.tarLs(tl_i)),1);
                    trial_list = find([obj.trials.tarL] == obj.tarLs(tl_i) & [obj.trials.outcome] == 1);
                    trial_list = setdiff(trial_list,1);
                    for t_i = 1:length(t_idx{3})
                        t_tmp = obj.trials(trial_list==t_idx{3}(t_i));
                        %cellsmat{tl_i,t_i,3} = t_tmp.export_as_formatted;
                        cellsmat{t_i,3} = t_tmp.export_as_formatted;
                    end
                end
            end
            for p_i = 1:2
                %trial_list = setdiff(t_idx{p_i},1);
                trial_list = find([obj.trials.outcome] == 1); % no failed trials 
                %trial_list = find([obj.trials.outcome] ~= -1);
                trial_list = intersect(trial_list, t_idx{p_i});
                %trial_list = setdiff(trial_list,[1 2 51:60]);
                trial_list = setdiff(trial_list,[1]);
                for t_i = 1:length(trial_list)
                    t_tmp = obj.trials(trial_list(t_i));
                    cellsmat{t_i,p_i} = t_tmp.export_as_formatted;  % each trial
                    xlim([-5 -4])
                    ifplot = true;
                    if (ifplot) 
                    subplot(2,1,1);
                    plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.x(2,:));
                    subplot(2,1,2);
                    plot(cellsmat{t_i, p_i}.t, cellsmat{t_i, p_i}.f(2,:));
                    end 
                end
            end

            
            % plot out the time skew
            ifplot = 0; %-test
            
            if (ifplot) 
                plt_offset = 2e-3/10;
                for p_i = 1:3
                    subplot(1,3,p_i); hold on;
                    switch p_i
                        case 1
                            title('no pert');
                        case 2
                            title('pulse pert');
                        case 3
                            title('stoc pert');
                    end
                    if (p_i ~= 3)
                        for t_i = 1:length(t_idx{p_i})
                            if(isempty(cellsmat{t_i, p_i}))
                                continue;
                            end
                            dt = diff(cellsmat{t_i,p_i}.t);
                            dt = [dt(1) dt];
                            plot(t_i*plt_offset + dt);
                        end
                    else
                        if isempty(t_idx{p_i}) 
                            continue
                        end
                        tiofst = 0; % plot offset
                        for tl_i = 1:size(cellsmat, 1)
                            for t_i = 1:length(cellsmat(tl_i,:,p_i))
                                if ~isempty(cellsmat{tl_i,t_i,p_i})
                                dt = diff(cellsmat{tl_i,t_i,p_i}.t);
                                dt = [dt(1) dt];
                                tiofst = tiofst+1;
                                plot(tiofst*plt_offset + dt);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        %%% with perturbations
%         function pert_ct = countPerturbation(obj)
%             % pert_ct = countPerturbation(); % return the trial # being
%             % perturbed
%             pert_ct = 0;
%             for trial_i = 1:length(obj.trials)
%                 ifpert = obj.trials(trial_i).ifpert;
%                 pert_ct = pert_ct + double(ifpert);
%             end
%         end
        %%% plot all session as a line
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CHECK THE TASK PERFORMANCE DATA  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dispTaskScanTrials(obj)
            % I need true combo here! see the ProcessRawData
            all_fTH = unique([obj.trials.fTh]);
            all_fTH = all_fTH(~isnan(all_fTH));
            all_tarL = unique([obj.trials.tarL]);
            all_tarL = all_tarL(~isnan(all_tarL));
            % for each task condition, print the fin_trials/all_trials
            for fTH_i = 1:length(all_fTH)
                for tarL_i = 1:length(all_tarL)
                    trial_all_idx = ([obj.trials.fTh] == all_fTH(fTH_i) ...
                        & [obj.trials.tarL] == all_tarL(tarL_i));
                    trial_fin_idx = ([obj.trials.fTh] == all_fTH(fTH_i) ...
                        & [obj.trials.tarL] == all_tarL(tarL_i) ...
                        & [obj.trials.outcome] == 1);
                    fprintf('FT: %02.0fN, Dist: %.02fm, trials: %02d / %02d\n', ...
                        all_fTH(fTH_i), all_tarL(tarL_i), ...
                        sum(trial_fin_idx), sum(trial_all_idx));
                end
            end
        end
        function axh = plotTaskTrialFinish(obj)
            combo_array = [obj.trials.comboTT];
            combo_all = unique(combo_array(~isnan([combo_array])));
            combo_min = min(combo_all);
            trial_sucess_idx = ([obj.trials.outcome]==1);
            trial_failure_idx= ([obj.trials.outcome]==0);
            axh = figure();
            hold on;
            plot(find(trial_sucess_idx), ...
                [obj.trials(trial_sucess_idx).comboTT]-combo_min+1, 'go'); %suceed trials
            plot(find(trial_failure_idx), ...
                [obj.trials(trial_failure_idx).comboTT]-combo_min+1, 'r*'); %failure trials
            xlabel('trial number');
            ylabel('target type number');
            title('Task trials finish condition');
        end
        function plotTaskStateMuskFig(obj)
            % TODO: ... THINK HOW TO MAKE THIS EASIER TO UNDERSTAND!!! 
            % display task masks in y-axis, blue: suceed, red: failure
            fields_t = fieldnames(obj.Data.TaskStateMasks); 
            fields_num = length(fields_t);
            % when they are the same length
            var_length = length(eval(['obj.Data.TaskStateMasks.' fields_t{1}]));
            var_ST = nan(1, var_length);
            for ii = 1:fields_num
                idx = [eval(['obj.Data.TaskStateMasks.' fields_t{ii}]) == 1];
                var_ST(1,idx) = ii;
            end
            
            % suceed color code 
            fields_s = fieldnames(obj.Data.TaskStateOutcomeMasks);
            fields_num = length(fields_s);
            % when they are the same length
            var_length = length(eval(['obj.Data.TaskStateOutcomeMasks.' fields_s{1}]));
            var_SF = nan(1, var_length);
            for ii = 1:fields_num
                idx = [eval(['obj.Data.TaskStateOutcomeMasks.' fields_s{ii}]) == 1];
                var_SF(1,idx) = ii;
            end
            col = ['bkrk']; 
            % plot var_all
            figure(); hold on;
            for ii = unique(var_SF(~isnan(var_SF)))
                idx = find(var_SF == ii);
                plot(idx, var_ST(idx), [col(ii) '*']);
                plot(idx, var_ST(idx), [col(ii) ]);
            end
            % notations
            yticks([1 2 3 4 5 6 7 8]);
            yticklabels(fields_t');
            xlabel('time pts');
            title('states though time');
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CHECK MESSAGE DATA AND HIGH SAMPLED DATA 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fh = plotTaskJointPosition(obj)  
            % plot the joint 4 position throughout this session
            fh = figure(); hold on;
            plot(obj.time, obj.Data.Position.JointPosition(4,:), 'b-o');
            plot(obj.data.t, obj.data.jp(4,:), '.');
            legend('msg', 'hi-sp');
            xlabel('time (s)');
            ylabel('pos (rad)');
            title('Joint 4 position')
        end
        function fh = plotTaskJointPosition_all(obj) 
            % plot 4 joints position
            fh = figure();
            for ji = 1:4
                axh(ji) = subplot(4,1,ji); hold on; grid on;
                plot(obj.time, obj.Data.Position.JointPosition(ji,:), 'b-o');
                plot(obj.data.t, obj.data.jp(ji,:), '.');
                if ji == 1
                    legend('msg', 'hi-sp');
                end
                xlabel('time (s)');
                ylabel('pos (rad)');
                title(['Joint ' num2str(ji) ' position'])
            end
            linkaxes(axh, 'x');
        end
%         function axh = plotTaskJointPositionh(obj, axh)
%             if nargin < 2
%                 axh = figure();
%             else 
%                 figure(axh); hold on;
%             end
%             position = obj.wam.jp'; 
%             %position = obj.wam.jt'; 
%             for sbpi = 1:4
%                 subplot(4,1,sbpi);
%                 plot(obj.wam_t, (position(sbpi,:))', '.');
%                 title(['J' num2str(sbpi)]);
%             end
%             ylabel('endpoint positions');
%             xlabel('time points');
%             %legend('x', 'y', 'z');
%             title('relative endpoint positions');
%         end 
%         function axh = taskJointPosition_relateve(obj) % Plot relative position. 
%             % plot 
%             figure();
%             axh(1) = subplot(2,1,1);
%             % plot([diff(obj.Data.Position.JointPosition,1,2)-obj.Data.Position.JointVelocity(:,2:end)]');
%             position_offset = obj.Data.Position.JointPosition(:,~isnan(obj.Data.Position.JointPosition(1,:)));
%             position_offset = repmat(position_offset(:,1),1,size(obj.Data.Position.JointPosition,2));
%             plot(obj.time, (obj.Data.Position.JointPosition - position_offset)');
%             legend('J1', 'J2', 'J3', 'J4');
%             axh(2) = subplot(2,1,2);
%             % plot([diff(obj.Data.Position.JointVelocity,1,2)-obj.Data.Position.JointTorque(:,2:end)]');
%             plot(obj.Data.Position.JointVelocity');
%             legend('J1', 'J2', 'J3', 'J4');
%             % notation
%             %set(axish(1), 'Ylim', [-0.02 0.02]);
%             ylabel(axh(1), 'joints positions');
%             %set(axish(2), 'Ylim', [-0.3 0.3]);
%             ylabel(axh(2), 'joints velocities');
%             xlabel(axh(2), 'time points');
%             
%         end
        function fh = plotTaskEndpointPosition(obj)
             % plot the y-axis position throughout this session
            fh = figure(); hold on;
            plot(obj.time, obj.Data.Position.Actual(:,2), 'b-o');
            plot(obj.data.t, obj.data.x(2,:), '.');
            legend('msg', 'hi-sp');
            xlabel('time (s)');
            ylabel('pos (m)');
            title('y axis position')
        end 
        function fh = plotTaskEndpointPosition_all(obj)
             % plot the xyz-axis position throughout this session
             fh = figure();
             axis_name = 'xyz';
             for ai = 1:3 % xyz
                 axh(ai)=subplot(3,1,ai);hold on;grid on;
                 plot(obj.time, obj.Data.Position.Actual(:,ai), 'b-o');
                 plot(obj.data.t, obj.data.x(ai,:), '.');
                 if ai == 1
                     legend('msg', 'hi-sp');
                 end
                 xlabel('time (s)');
                 ylabel('pos (m)');
                 title([ axis_name(ai) ' axis position']);
             end
             linkaxes(axh, 'x');
        end 
        function fh = plotTaskEndpointVelocity(obj)
             % plot the y-axis velocity throughout this session
            fh = figure(); hold on;
            plot(obj.time, obj.Data.Velocity.Actual(2,:), 'b-o');
            plot(obj.data.t, obj.data.v(2,:), '.');
            legend('msg', 'hi-sp');
            xlabel('time (s)');
            ylabel('vel (m/s)');
            title('y axis velocity')
        end 
        function fh = plotTaskEndpointVelocity_all(obj)
             % plot the xyz-axis velocity throughout this session
             fh = figure();
             axis_name = 'xyz';
             for ai = 1:3 % xyz
                 axh(ai)=subplot(3,1,ai);hold on;grid on;
                 plot(obj.time, obj.Data.Velocity.Actual(ai,:), 'b-o');
                 plot(obj.data.t, obj.data.v(ai,:), '.');
                 if ai == 1
                     legend('msg', 'hi-sp');
                 end
                 xlabel('time (s)');
                 ylabel('vel (m)');
                 title([ axis_name(ai) ' axis velocity']);
             end
             linkaxes(axh, 'x');
        end 
%         function axh = plotTaskEndpointPositionh(obj, axh)
%             if nargin < 2
%                 axh = figure();
%             else 
%                 figure(axh); hold on;
%             end
%             position = obj.wamp_h'; 
%             plot(obj.wam_t, (position)', '.');  
%             ylabel('endpoint positions');
%             xlabel('time points');
%             legend('x', 'y', 'z');
%             title('relative endpoint positions');
%         end 
        function fh = plotTaskEndpointForce(obj)
             % plot the y-axis force throughout this session
            fh = figure(); 
            
            axh(1) = subplot(2,1,1); hold on;
            plot(obj.time, obj.force(2,:), 'b-o');
            plot(obj.data.t, obj.data.f(2,:), '.');
            legend('msg', 'hi-sp');
            xlabel('time (s)');
            ylabel('force (N)');
            title('y axis force')
            
            axh(2) = subplot(2,1,2); hold on;
            plot(obj.time, sqrt(obj.force(1,:).^2+...
                                obj.force(2,:).^2+...
                                obj.force(3,:).^2), 'b-o');
            plot(obj.data.t, sqrt(obj.data.f(1,:).^2+...
                                obj.data.f(2,:).^2+...
                                obj.data.f(3,:).^2), '.');
            xlabel('time (s)');
            ylabel('force (N)');
            title('normalized force')
            
            linkaxes(axh, 'x');
                            
        end 
        function fh = plotTaskEndpointForce_all(obj)
             % plot the xyz-axis force throughout this session
             fh = figure();
             axis_name = 'xyz';
             for ai = 1:3 % xyz
                 axh(ai)=subplot(3,1,ai);hold on;grid on;
                 plot(obj.time, obj.force(ai,:), 'b-o');
                 plot(obj.data.t, obj.data.f(ai,:), '.');
                 if ai == 1
                     legend('msg', 'hi-sp');
                 end
                 xlabel('time (s)');
                  ylabel('force (N)');
                 title([ axis_name(ai) ' axis force']);
             end
             linkaxes(axh, 'x');
        end 
%         function axh = plotTaskForceData(obj, axh)
%             if nargin<2
%                 axh = figure();
%             else
%                 figure(axh); hold on;
%             end
%             % force = obj.Data.Force.Sensor(1:3,:); 
%             force = obj.force;
%             plot(obj.time, force');
%             ylabel('force (N)');
%             xlabel('time');
%             legend('x', 'y', 'z'); % remember to alter the axis 
%             title('Force data');
%         end
%         function axh = plotTaskForceDatah(obj, axh)
%             if nargin<2
%                 axh = figure();
%             else
%                 figure(axh); hold on;
%             end
%             % force = obj.Data.Force.Sensor(1:3,:); 
%             force = obj.data.f;
%             plot(obj.data.t, force', '.');
%             ylabel('force (N)');
%             xlabel('time');
%             legend('x', 'y', 'z'); % remember to alter the axis 
%             title('Force data high sample');
%         end
%         function axh = taskForceDataMag(obj, axh)
%             if nargin<2
%                 axh = figure();
%             else
%                 figure(axh); hold on;
%             end
%             % force = obj.Data.Force.Sensor(1:3,:); 
%             force = obj.force;
%             forceMag = sqrt(force(1,:).^2 + force(2,:).^2);
%             plot(obj.time, forceMag');
%             ylabel('force (N)');
%             xlabel('time');
%             legend('x', 'y', 'z'); % remember to alter the axis 
%             title('Force data');
%         end
%         function taskEndpointPosition_relative(obj)
%             figure();
%             position = obj.Data.Position.Actual'; 
%             % Use first element as offset
%             % position_offset = position(:,~isnan(position(1,:)));
%             position_offset = obj.hand_pos_offset;
%             position_offset = repmat(position_offset(:,1),1,size(position,2));
%             plot(obj.time, (position - position_offset)');  
%             ylabel('relative endpoint positions');
%             xlabel('time points');
%             legend('x', 'y', 'z');
%             title('relative endpoint positions');
%         end
%         function axh = plotTaskEndpointForceh(obj, axh)
%             if nargin < 2
%                 axh = figure();
%             else 
%                 figure(axh); hold on;
%             end
%             %position = obj.wamp_h'; 
% 
%             plot(obj.data.t, obj.data.f(2,:));  
%             ylabel('endpoint positions');
%             xlabel('time points');
%             legend('x', 'y', 'z');
%             title('relative endpoint positions');
%         end 
        function axh = plotTaskjointTorqeh(obj, axh)
            if nargin < 2
                axh = figure();
            else 
                figure(axh); hold on;
            end
            %position = obj.wamp_h'; 
            axh1 = subplot(4,1,1);
            plot(obj.wam.time, smooth((obj.data.tq(1,:)),1));  
            ylabel('J1'); grid on;
            axh2 = subplot(4,1,2);
            plot(obj.wam.time, smooth((obj.data.tq(2,:)),2)); 
            ylabel('J2'); grid on;
            axh3 = subplot(4,1,3);
            plot(obj.wam.time, smooth((obj.data.tq(3,:)),3));  
            ylabel('J3'); grid on;
            axh4 = subplot(4,1,4);
            plot(obj.wam.time, smooth((obj.data.tq(4,:)),4));  
            ylabel('J4'); grid on;
            linkaxes([axh1 axh2 axh3 axh4], 'xy');
            xlabel('time points');
            %legend('x', 'y', 'z');
            title('joint 4 torque positions');
        end 
        function axh = plotAddTrialMark(obj, axh)
            % read already exist figure, add trial mark on it
            if (~exist('axh', 'var'))
                disp('no figure was readed, ABORT!');
            end
            if isa(axh, 'matlab.ui.Figure')
                axh = figure(axh); hold on; % stack
            elseif isa(axh, 'matlab.graphics.axis.Axes')
                subplot(axh); hold on;
            end
            for trial_i = 1:length(obj.trials)
                t = obj.trials(trial_i).time_orn(1);
                line([t t], [-1 1]);
                text(t, 7, ['trial' num2str(trial_i)]);
            end
            
            legend off;
        end
        function axh = taskEPP_FToverlap_ns(obj) % overlapping endpoint position and FT in one axis, non-scale
            figure(); hold on;
            position = obj.Data.Position.Actual'; 
            % Use first element as offset
            % position_offset = position(:,~isnan(position(1,:)));
            position_offset = obj.hand_pos_offset;
            position_offset = repmat(position_offset(:,1),1,size(position,2));
            % normalize the (position - position_offset)/range
            position_centered = position - position_offset;
            position_nan_idx = isnan(position_centered(1,:));
            position_range = range(position_centered(:,~position_nan_idx)');
            position_norm = position_centered./repmat(position_range,size(position_centered,2),1)';
            % normalize the force data
            force = obj.force;
            force_nan_idx = isnan(force(1,:));
            force_range = range(force(:,~force_nan_idx)');
            force_norm = force./repmat(force_range,size(force,2),1)';
            % convert the force data into nan when position is nan. 
            force_norm(:, position_nan_idx) = nan; % convert same size as position
            % plot(obj.time, ([position_norm; force_norm])');  
            ylabel_str = 'xy';
            for ii = 1:2 % x- and y- axis
                axh(ii) = subplot(2,1,ii); hold on;
                plot(obj.time, position_centered(ii,:)' * 10); % in-acurate maxium values. 
                plot(obj.time, force_norm(ii,:)', '--'); 
                ylabel(ylabel_str(ii)); 
                legend(['P' ylabel_str(ii)], ['F' ylabel_str(ii)]);
            end
            xlabel('time points');
            
            title('norm Positions and Forces');
            
        end
%         function axh = addmark_STMOV(obj, axh) % add lines showing mov state. 
%             % not good for now, as I do not align it good with time.
%             % how to addline without add the legend???
%             mov_mask = obj.Data.TaskStateMasks.Move;
%             mov_diff = [0 diff(mov_mask)]; 
%             mov_idx = find((mov_diff == 1) & (mov_mask == 1));
%             % plot bars in the axh
%             for axh_i = 1:length(axh)
%                 ylim_range = get(axh(axh_i), 'ylim');
%                 v= ver('MATLAB'); 
%                 for ii = 1:length(mov_idx)
%                     if str2double(v.Version) <= 9.0 % less than matlab 2018
%                         line(axh(axh_i), [mov_idx(ii) mov_idx(ii)], [ylim_range(1) ylim_range(2)]);
%                     else
%                         xline(axh(axh_i), mov_idx(ii));
%                     end
%                 end
%             end
%         end
        
        %%% plot overlapped release curves
%         function axh = plotTrialfyPosition(obj, axh)
%             if nargin < 2
%                 axh = figure();
%             else
%                 figure(axh);
%             end
%             hold on;
%             trials = obj.trials;
%             for trial_i = 1:length(trials)
%                 plot(trials(trial_i).time, trials(trial_i).position);
%             end
%             xlabel('time');
%             ylabel('position');
%             title('all trials position');
%         end
%         function axh = plotTrialfyForce(obj, axh)
%             if nargin < 2
%                 axh = figure();
%             else
%                 figure(axh);
%             end
%             hold on;
%             trials = obj.trials;
%             for trial_i = 1:length(trials)
%                 plot(trials(trial_i).time, trials(trial_i).force);
%             end
%             xlabel('time');
%             ylabel('force');
%             title('all trials force');
%         end
%         function axh = plotTrialfyForce_xy(obj, axh)
%             if nargin < 2
%                 axh = figure();
%             else
%                 figure(axh);
%             end
%             
%             hold on;
%             trials = obj.trials;
%             subplot(2,1,1); hold on;
%             for trial_i = 1:length(trials)
%                 plot(trials(trial_i).time, trials(trial_i).force(1,:));
%             end
%             xlim([-1 1]);
%             xlabel('time');
%             ylabel('force x');
%             title('force xy (m)');
%             subplot(2,1,2); hold on;
%             for trial_i = 1:length(trials)
%                 plot(trials(trial_i).time, trials(trial_i).force(2,:));
%             end
%             xlim([-1 1]);
%             title('force y');
%             xlabel('time');
%             ylabel('position y (m)');
%         end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOT TRIALFY DATA ON SPECIFY ALIGNED EVENTS 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function axh = plotTrialfyPositionh(obj, axh)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).data.t_shift, trials(trial_i).data.x(2,:));
            end
            xlabel('time');
            ylabel('position');
            title('all trials position');
        end
        function axh = plotTrialfyPositionh_all(obj, axh)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            axh_arr = 'xyz';
            trials = obj.trials;
            for axi = 1:3
                subplot(3,1,axi); hold on;
                for trial_i = 1:length(trials)
                    plot(trials(trial_i).data.t_shift, trials(trial_i).data.x(axi,:));
                end
                xlim([-1 1]);
                xlabel('time');
                ylabel('position (m)');
                title(['position ' axh_arr(axi)]);
                %             subplot(2,1,2); hold on;
                %             for trial_i = 1:length(trials)
                %                 plot(trials(trial_i).data.t_shift, trials(trial_i).data.x(2,:));
                %             end
                %             xlim([-1 1]);
                %             title('position y');
                %             xlabel('time');
                %             ylabel('position y (m)');
            end
            % title('all trials position');
        end
        function axh = plotTrialfyForceh(obj, axh)
            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).data.t_shift, trials(trial_i).data.f(2,:));
            end
            xlim([-1 1]);
            xlabel('time (s)');
            ylabel('force (N)');
            title('all trials force raw');
        end
        function axh = plotTrialfyForceh_all(obj, axh)
            % ...TODO: should change into xyz
            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            axh_names = 'xyz';
            for axi = 1:3
                axh(axi) = subplot(3,1,axi); hold on;
                for trial_i = 1:length(trials)
                    plot(trials(trial_i).data.t_shift, trials(trial_i).data.f(axi,:));
                end
                xlabel('time');
                ylabel('force (N)');
                title(['force ' axh_names(axi)]);
                
            end
%             xlim([-1 1]);
%             subplot(2,1,2); hold on;
%             for trial_i = 1:length(trials)
%                 plot(trials(trial_i).data.t_shift, trials(trial_i).data.f(2,:));
%             end
%             xlabel('time');
%             ylabel('force y (N)');
%             title('force y ');
            linkaxes(axh, 'x');
            xlim([-1 1]);
        end
%         function axh = plotTrialfyForceh_xy(obj, axh)
%             % ...TODO: should change into xyz
%             if nargin < 2
%                 axh = figure();
%             else
%                 figure(axh);
%             end
%             hold on;
%             trials = obj.trials;
%             subplot(2,1,1); hold on;
%             for trial_i = 1:length(trials)
%                 plot(trials(trial_i).data.t_shift, trials(trial_i).data.f(1,:));
%             end
%             xlabel('time');
%             ylabel('force x (N)');
%             title('force x ');
%             xlim([-1 1]);
%             subplot(2,1,2); hold on;
%             for trial_i = 1:length(trials)
%                 plot(trials(trial_i).data.t_shift, trials(trial_i).data.f(2,:));
%             end
%             xlabel('time');
%             ylabel('force y (N)');
%             title('force y ');
%             xlim([-1 1]);
%         end
%         
%         function axh = plotMeantrialForce(obj, axh, cor_i)
%             % plot the meaned trial Force according to the task condition
%              % plot the meaned trial Velocity according to the task condition
%             % if did not calculate the mean and var, calculate
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(cor_i)
%                 cor_i = 0;
%             end
%             if length(cor_i) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % conditions:
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             if obj.ssnum == 2459 % 
%                 all_tarL = setdiff(all_tarL, [0.01]);
%             end
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             %axhf = figure();
%             l_h = [];
%             labels = {};
%             label_i = 0;
%             for fTH_i = 1:length(all_fTH)
%                 for tarL_i = 1:length(all_tarL)
%                     %col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
% 
%                     %trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) ...
%                         & [obj.trials.tarL]==all_tarL(tarL_i)...
%                         & [obj.trials.outcome]==1; 
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTH(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [t, f_mat] = alignTrialDataF(obj, trials_idx);
%                     force_mean = mean(f_mat(:,:,xyi)); %only y direction
%                     force_std = std(f_mat(:,:,xyi));  
%                     % mean line
%                     l_h = [l_h plot(t, force_mean, 'LineWidth', 3, 'Color', obj.col_vec(col_i,:))];
%                     % 1std shade
%                     force_up = force_mean + force_std;
%                     force_dn = force_mean - force_std;
%                     [axh, msg] = jbfill(t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
%                 end
%             end
%             xlabel('time aligned at MOV signal');
%             ylabel([xy_char(xyi) ' dir force (N)']);
%             xlim([-0.5 0.4]);
%             title('Force signal');
%             %legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
%             legend(l_h, labels);
%         end
%         function axh = plotMeantrialForcePert(obj, axh, col_i)
%             % plot the meaned trial Force according to the task condition
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             if ~exist('axhp', 'var')
%                 axh = figure();
%             else
%                 % for figure
%                 if strcmp(axh.Type, 'figure')
%                     axh = figure(axh);
%                     hold on;
%                 elseif strcmp(axh.Type, 'axes')
%                 % for axis
%                     axh = subplot(axh);
%                     hold on;
%                 else
%                     disp('Uknown handle type, return!');
%                     return;
%                 end
%             end
%             % plot color
%             if ~exist('col_i', 'var')
%                 col_i = 1;
%             end
%             
%             xy_char = 'xy';
%             %axhf = figure();
%             l_h = [];
%             labels = {};
%             label_i = 0;
%             for trial_i = 1:length(obj.trials)
%                 if (obj.trials(trial_i).ifpert)
%                     obj.trials(trial_i) = alignPertInit(obj.trials(trial_i), obj);
%                 end
%             end
%             for fTH_i = 1:length(all_fTH)
%                 for tarL_i = 1:length(all_tarL)
%                     %col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
% 
%                     %trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) ...
%                         & [obj.trials.tarL]==all_tarL(tarL_i)...
%                         & [obj.trials.outcome]==1 ...
%                         & [obj.trials.ifpert]==1; 
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTH(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, resample_f] = trialDataResampleFT(obj, trials_idx);
%                     % substract the <0 average value
%                     time_idx = resample_t < 0;
%                     for trial_i = 1:size(resample_f, 1)
%                         resample_f(trial_i,:,xyi) = resample_f(trial_i,:,xyi); %- mean(resample_f(trial_i,time_idx,xyi));
%                     end
%                     force_mean = mean(resample_f(:,:,xyi)); %only y direction
%                     force_std = std(resample_f(:,:,xyi));  
%                     % mean line
%                     l_h = [l_h plot(resample_t, force_mean, 'LineWidth', 3, 'Color', obj.col_vec(col_i,:))];
%                     % 1std shade
%                     force_up = force_mean + force_std;
%                     force_dn = force_mean - force_std;
%                     [axhf, msg] = jbfill(resample_t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
%                 end
%             end
%             xlabel('time aligned at MOV signal');
%             ylabel([xy_char(xyi) ' dir force (N)']);
%             xlim([-0.05 0.95]);
%             title('Force signal');
%             %legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
%             legend(l_h, labels);
%         end
%         function axh = plotMeantrialPos(obj, axh, col_i)
%             % plot the meaned trial Position according to the task condition
%             % if did not calculate the mean and var, calculate
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(col_i)
%                 col_i = 1;
%             end
%             if length(col_i) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % plot color
%             if ~exist('col_i', 'var')
%                 col_i = 1;
%             end
%             
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             if obj.ssnum == 2459
%                 all_tarL = setdiff(all_tarL, [0.01]);
%             end
%             
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot position
%             % axhp = figure(); % when plot in a single figure
%             l_h = [];
%             labels = {};
%             label_i = 0;
%             for fTH_i = 1:length(all_fTH)
%                 for tarL_i = 1:length(all_tarL)
%                     %col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
%                     trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTH(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, resample_p, ~] = alignTrialDataX(obj, trials_idx);
%                     pos_mean = mean(resample_p(:,:,xyi), 1, 'omitnan') - obj.endpoint0(xyi); %only y direction
%                     pos_std = std(resample_p(:,:,xyi), 1, 'omitnan');  
%                     % mean line
%                     l_h = [l_h plot(resample_t, pos_mean, 'LineWidth', 3, 'Color', obj.col_vec(col_i,:))];
%                     % 1std shade
%                     pos_up = pos_mean + pos_std/2;
%                     pos_dn = pos_mean - pos_std/2;
%                     [axh, msg] = jbfill(resample_t, pos_up, pos_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
%                 end
%             end
%             xlabel('time aligned at MOV signal');
%             ylabel([xy_char(xyi) 'dir position (m)']);
%             xlim([-0.2 0.5]);
%             title('Position signal');
%             %legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
%             legend(l_h, labels);
%         end
%         function axh = plotMeantrialPosPert(obj, axhp, col_i)
%             % plot the meaned trial Position according to the task condition
%             % if did not calculate the mean and var, calculate
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             if obj.ssnum == 2459 % 
%                 all_tarL = setdiff(all_tarL, [0.01]);
%             end
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot position
%             if ~exist('axhp', 'var')
%                 axhp = figure();
%             else
%                 % for figure
%                 if strcmp(axhp.Type, 'figure')
%                     axhp = figure(axhp);
%                     hold on;
%                 elseif strcmp(axhp.Type, 'axes')
%                 % for axis
%                     axhp = subplot(axhp);
%                     hold on;
%                 else
%                     disp('Uknown handle type, return!');
%                     return;
%                 end
%             end
%             % plot color
%             if ~exist('col_i', 'var')
%                 col_i = 1;
%             end
%             l_h = [];
%             labels = {};
%             label_i = 0;  
%             % align for the perturbation time
%             for trial_i = 1:length(obj.trials)
%                 if (obj.trials(trial_i).ifpert)
%                     obj.trials(trial_i) = alignPertInit(obj.trials(trial_i), obj);
%                 end
%             end
%             % get the mean
%             for fTH_i = 1:length(all_fTH)
%                 for tarL_i = 1:length(all_tarL)
%                     %col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
%                     trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i) & [obj.trials.ifpert];
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTH(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, resample_p, ~] = trialDataAlignWAM(obj, trials_idx);
%                     % substract average before pert
%                     resample_t_0idx = find(resample_t == 0);
%                     if resample_t_0idx >= 50 
%                         resample_t_idx = (resample_t_0idx - 49):resample_t_0idx;
%                     else 
%                         resample_t_idx = 1:resample_t_0idx;
%                     end
%                     % process each trial to subtract the average
%                     for trial_i = 1:size(resample_p(:,:,xyi),1)
%                         data = resample_p(trial_i,resample_t_idx,xyi);
%                         resample_p(trial_i, :, xyi) = resample_p(trial_i, :, xyi) - mean(data);
%                     end
%                     % pos_mean = mean(resample_p(:,:,xyi), 'omitnan') - obj.endpoint0(xyi); %only y direction
%                     pos_mean = mean(resample_p(:,:,xyi), 'omitnan'); %only y direction
%                     pos_std = std(resample_p(:,:,xyi), 'omitnan');  
%                     % mean line
%                     l_h = [l_h plot(resample_t, pos_mean, 'LineWidth', 3, 'Color', obj.col_vec(col_i,:))];
%                     % 1std shade
%                     force_up = pos_mean + pos_std;
%                     force_dn = pos_mean - pos_std;
%                     [axh, msg] = jbfill(resample_t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
%                 end
%             end
%             xlabel('time aligned at MOV signal');
%             ylabel([xy_char(xyi) 'dir position (m)']);
%             xlim([-0.1 0.7]);
%             title('Position signal');
%             %legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
%             legend(l_h, labels);
%         end
%         function axh = plotMeantrialVel(obj, axh, col_i)
%             % plot the meaned trial Velocity according to the task condition
%             % if did not calculate the mean and var, calculate
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(col_i)
%                 col_i = 1;
%             end
%             if length(col_i) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % plot color
%             if ~exist('col_i', 'var')
%                 col_i = 1;
%             end
%             
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot velocity
%             %axhv = figure();
%             l_h = [];
%             labels = {};
%             label_i = 0;
%             for fTH_i = 1:length(all_fTH)
%                 for tarL_i = 1:length(all_tarL)
%                     %col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
%                     trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTH(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, ~, resample_v] = alignTrialDataX(obj, trials_idx);
%                     vel_mean = mean(resample_v(:,:,xyi), 'omitnan'); %only y direction
%                     vel_std = std(resample_v(:,:,xyi), 'omitnan');  
%                     % mean line
%                     try
%                         col_tmp = obj.col_vec(col_i,:);
%                     catch
%                         col_tmp = obj.col_vec(mod(col_i, size(obj.col_vec, 1)),:);
%                     end
%                     l_h = [l_h plot(resample_t, vel_mean, 'LineWidth', 3, 'Color', col_tmp)];
%                     % 1std shade
%                     vel_up = vel_mean + vel_std;
%                     vel_dn = vel_mean - vel_std;
%                     [~, ~] = jbfill(resample_t, vel_up, vel_dn, col_tmp, col_tmp, 1, 0.3);
%                 end
%             end
%             xlabel('time aligned at MOV signal');
%             ylabel([xy_char(xyi) 'dir velocity (m/s)']);
%             xlim([-0.2 0.5]);
%             title('Velocity signal');
%             %legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
%             legend(l_h, labels);
%         end
%         
%         %%%%%%%%%%%%%%%%%%%% useless code block here!!! %%%%%%%%%%%%%%%%%%%
%         function axhf = plotMeantrialForce_sameCond(obj)
%             % plot the meaned trial Force according to the task condition
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % axhf = figure();
%             for tarL_i = 1:length(all_tarL)
%                 axhf(tarL_i) = figure(tarL_i); hold on;
%                 l_h = [];
%                 labels = {};
%                 label_i = 0;
%                 trials_all = obj.trials([obj.trials.tarL] == all_tarL(tarL_i));
%                 all_fTh = unique([trials_all.fTh]); % condition specific all targets.
%                 for fTH_i = 1:length(all_fTh)
%                     col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
% 
%                     %trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     trials_idx = [obj.trials.fTh]==all_fTh(fTH_i) ...
%                         & [obj.trials.tarL]==all_tarL(tarL_i)...
%                         & [obj.trials.outcome]==1; 
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTh(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, resample_f] = trialDataResampleFT(obj, trials_idx);
%                     force_mean = mean(resample_f(:,:,xyi)); %only y direction
%                     force_std = std(resample_f(:,:,xyi));  
%                     % mean line
%                     try
%                         col_tmp = obj.col_vec(col_i,:);
%                     catch
%                         col_tmp = obj.col_vec(mod(col_i-1, size(obj.col_vec, 1))+1,:);
%                     end
%                     l_h = [l_h plot(resample_t, force_mean, 'LineWidth', 3, 'Color', col_tmp)];
%                     % 1std shade
%                     %force_up = force_mean + force_std/2;
%                     %force_dn = force_mean - force_std/2;
%                     %[axhf, msg] = jbfill(resample_t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
%                 end
%                 legend(l_h, labels);
%                 xlabel('time aligned at MOV signal');
%                 ylabel([xy_char(xyi) ' dir force (N)']);
%                 xlim([-0.5 0.4]);
%                 title('Force signal');
%             end
% 
%         end
%         function axhp = plotMeantrialPos_sameCond(obj)
%             % plot the meaned trial Position according to the task condition
%             % if did not calculate the mean and var, calculate
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot position
%             for tarL_i = 1:length(all_tarL)
%                 axhp(tarL_i) = figure();
%                 l_h = [];
%                 labels = {};
%                 label_i = 0;
%                 trials_all = obj.trials([obj.trials.tarL] == all_tarL(tarL_i));
%                 all_fTh = unique([trials_all.fTh]); % condition specific all targets.
%                 for fTH_i = 1:length(all_fTh)
%                     col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
%                     trials_idx = [obj.trials.fTh]==all_fTh(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTh(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, resample_p, ~] = trialDataAlignWAM(obj, trials_idx);
%                     force_mean = mean(resample_p(:,:,xyi), 'omitnan') - obj.endpoint0(xyi); %only y direction
%                     force_std = std(resample_p(:,:,xyi), 'omitnan');  
%                     % mean line
%                     try
%                         col_tmp = obj.col_vec(col_i,:);
%                     catch
%                         col_tmp = obj.col_vec(mod(col_i-1, size(obj.col_vec, 1))+1,:);
%                     end
%                     l_h = [l_h plot(resample_t, force_mean, 'LineWidth', 3, 'Color', col_tmp)];
%                     % 1std shade
%                     %force_up = force_mean + force_std/2;
%                     %force_dn = force_mean - force_std/2;
%                     %[axh, msg] = jbfill(resample_t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
%                 end
%                 legend(l_h, labels);
%                 xlabel('time aligned at MOV signal');
%                 ylabel([xy_char(xyi) 'dir position (m)']);
%                 xlim([-0.2 0.5]);
%                 title('Position signal');
%             end
%             
%         end
%         function axhv = plotMeantrialVel_sameCond(obj, invert)
%             % plot the meaned trial Velocity according to the task condition
%             % if did not calculate the mean and var, calculate
%             % 'invert == -1' will flirt over the trial curve
%             if ~exist('invert', 'var') 
%                 invert = 1;
%             end
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot velocity
%             
%             for tarL_i = 1:length(all_tarL)
%                 axhv(tarL_i) = figure();
%                 l_h = [];
%                 labels = {};
%                 label_i = 0;
%                 trials_all = obj.trials([obj.trials.tarL] == all_tarL(tarL_i));
%                 all_fTh = unique([trials_all.fTh]);
%                 for fTH_i = 1:length(all_fTh)
%                     col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
%                     trials_idx = [obj.trials.fTh]==all_fTh(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTh(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, ~, resample_v] = trialDataAlignWAM(obj, trials_idx);
%                     vel_mean = mean(resample_v(:,:,xyi), 'omitnan'); %only y direction
%                     force_std = std(resample_v(:,:,xyi), 'omitnan');  
%                     % mean line
%                     try
%                         col_tmp = obj.col_vec(col_i,:);
%                     catch
%                         col_tmp = obj.col_vec(mod(col_i-1, size(obj.col_vec, 1))+1,:);
%                     end
%                     l_h = [l_h plot(resample_t, invert*vel_mean, 'LineWidth', 3, 'Color', col_tmp)];
%                     % 1std shade
%                     %force_up = force_mean + force_std/2;
%                     %force_dn = force_mean - force_std/2;
%                     %[axh, msg] = jbfill(resample_t, force_up, force_dn, col_tmp, col_tmp, 1, 0.3);
%                 end
%             xlabel('time aligned at MOV signal');
%             ylabel([xy_char(xyi) 'dir velocity (m/s)']);
%             xlim([-0.2 0.5]);
%             title('Velocity signal');
%             %legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
%             legend(l_h, labels);
%             end
%         end
%         function axhv = plotMeantrialVel_sameCond_overlap(obj, invert, axhv, col_i)
%             % plot the meaned trial Velocity according to the task condition
%             % if did not calculate the mean and var, calculate
%             % 'invert == -1' will flirt over the trial curve
%             % Required each session tobe only one directional target.
%             if ~exist('invert', 'var') 
%                 invert = 1;
%             end
%             if ~exist('axhv', 'var')
%                 axhv_flag = -1;
%             elseif ~isstruct(axhv)
%                 axhv_flag = -1;
%             else
%                 axhv_flag = 1;
%             end
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             all_tarR = unique([obj.trials.tarR]);
%             all_tarR = all_tarR(~isnan(all_tarR));
%             % assume this session only have x- or y- trials
%             if isempty(setdiff(all_tarR, [0,4])) %only y direction
%                 xyi = 1;
%             elseif isempty(setdiff(all_tarR, [2, 6]))
%                 xyi = 2;
%             end
%             xy_char = 'xy';
%             % plot velocity
%             % check if axhv is the same length with all_tarL. 
%             %   if so, plot on the axhv. If not, plot on new figures 
%             if axhv_flag == -1 %| length(axhv.axih)~=length(all_tarL)
%                 pltopt = 'new'; 
%                 axhv = [];
%                 axhv.figh = figure();
%             else
%                 pltopt = 'ovl'; % overlap
%                 axhv.figh = figure(axhv.figh);
%             end
%             for tarL_i = 1:length(all_tarL)
%                 if strcmp(pltopt, 'new')
%                     axhv.axih(tarL_i) = subplot(length(all_tarL), 1, tarL_i);
%                 elseif strcmp(pltopt, 'ovl')
%                     axhv.axih(tarL_i) = subplot(axhv.axih(tarL_i));
%                 end
%                 l_h = [];
%                 labels = {};
%                 label_i = 0;
%                 trials_all = obj.trials([obj.trials.tarL] == all_tarL(tarL_i));
%                 all_fTh = unique([trials_all.fTh]);
%                 for fTH_i = 1:length(all_fTh)
%                     % col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
%                     hold on;
%                     trials_idx = [obj.trials.fTh]==all_fTh(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
%                     label_i = label_i + 1;
%                     labels{label_i} = [num2str(all_fTh(fTH_i)), 'N, ', num2str(all_tarL(tarL_i)*100) 'cm'];
%                     [resample_t, ~, resample_v] = trialDataAlignWAM(obj, trials_idx);
%                     vel_mean = mean(resample_v(:,:,xyi), 'omitnan'); %only y direction
%                     force_std = std(resample_v(:,:,xyi), 'omitnan');  
%                     % mean line
%                     try
%                         col_tmp = obj.col_vec(col_i,:);
%                     catch
%                         col_tmp = obj.col_vec(mod(col_i-1, size(obj.col_vec, 1))+1,:);
%                     end
%                     l_h = [l_h plot(resample_t, invert*vel_mean, 'LineWidth', 3, 'Color', col_tmp)];
%                     % 1std shade
%                     %force_up = force_mean + force_std/2;
%                     %force_dn = force_mean - force_std/2;
%                     %[axh, msg] = jbfill(resample_t, force_up, force_dn, col_tmp, col_tmp, 1, 0.3);
%                 end
%             xlabel('time aligned at MOV signal');
%             ylabel([xy_char(xyi) 'dir velocity (m/s)']);
%             xlim([-0.2 0.5]);
%             title(['3, 12, 21N at ' num2str(all_tarL(tarL_i)*100) 'cm']);
%             %legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
%             %legend(l_h, labels);
%             end
%         end
%         
%         function [axhf, axhp, axhv] = plotSameTrial(obj)
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             outcome = [obj.trials.outcome];
%             % fTH now only have 2 vals
%             % tarL now only have 2 vals
%             axhf = figure(); % force figure
%             set(axhf, 'Visible', 'off');
%             axhp = figure(); % position figure
%             set(axhp, 'Visible', 'off');
%             axhv = figure(); % velocity figure
%             set(axhv, 'Visible', 'off');
%             for fTH_i = 1:length(all_fTH)
%                 fprintf('\n plotting fTH:%2.1f  ', all_fTH(fTH_i));
%                 for tarL_i = 1:length(all_tarL)
%                     fprintf('target: %0.2f', all_tarL(tarL_i));
%                     trials_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)));
%                     trial_num = sum(trials_idx);
%                     trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) & outcome==1); %succeed
%                     %trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) ); %all trials
%                     col_i = (fTH_i-1)*length(all_fTH) + tarL_i;
%                     if col_i>=7
%                         col_i = mod(col_i-1,8)+1;
%                     end
%                     for trial_i = find(trialsS_idx)
%                         % plot force
%                         % data
%                         try
%                             force = obj.trials(trial_i).force_h;
%                             time = obj.trials(trial_i).force_t;
%                         catch
%                             force = obj.trials(trial_i).force;
%                             time = obj.trials(trial_i).time;
%                         end
%                         % figure
%                         set(0, 'CurrentFigure', axhf);
%                         % x, y seperately
%                         subplot(2,1,1); hold on; plot(time, force(1,:), 'Color', obj.col_vec(col_i,:));
%                         subplot(2,1,2); hold on; plot(time, force(2,:), 'COlor', obj.col_vec(col_i,:));
%                         % plot position
%                         % data
%                         try
%                             position = obj.trials(trial_i).position_h;
%                             time = obj.trials(trial_i).position_t;
%                         catch
%                             position = obj.trials(trial_i).position;
%                             time = obj.trials(trial_i).time;
%                         end
%                         set(0, 'CurrentFigure', axhp);
%                         % x, y seperately
%                         subplot(2,1,1); hold on; plot(time, position(1,:), 'Color', obj.col_vec(col_i,:));
%                         subplot(2,1,2); hold on; plot(time, position(2,:), 'Color', obj.col_vec(col_i,:));
%                         try
%                             velocity = obj.trials(trial_i).velocity_h;
%                             time = obj.trials(trial_i).position_t;
%                         catch
%                             velocity = obj.trials(trial_i).velocity;
%                             time = obj.trials(trial_i).time;
%                         end
%                         set(0, 'CurrentFigure', axhv);
%                         subplot(2,1,1); hold on; plot(time, velocity(1,:), 'Color', obj.col_vec(col_i,:));
%                         subplot(2,1,2); hold on; plot(time, velocity(2,:), 'Color', obj.col_vec(col_i,:));
%                     end
%                 end
%             end
%             
%             xrangeF = [-0.5, 0.5];
%             figure(axhf); 
%             subplot(2,1,1);
%             title('force data x');
%             xlim(xrangeF);
%             subplot(2,1,2);
%             title('force data y');
%             xlim(xrangeF);
%             
%             xrangeP = [-0.2, 0.8];
%             figure(axhp);
%             subplot(2,1,1);
%             title('robot position x');
%             xlim(xrangeP);
%             subplot(2,1,2);
%             title('robot position y');
%             xlim(xrangeP);
%             
%             xrangeV = [-0.2, 0.8];
%             figure(axhv);
%             subplot(2,1,1);
%             title('robot velocity x');
%             xlim(xrangeV);
%             subplot(2,1,2);
%             title('robot velocity y');
%             xlim(xrangeV);
%             
%             set(axhf, 'Visible', 'on');
%             set(axhp, 'Visible', 'on');
%             set(axhv, 'Visible', 'on');
%         end
%         function [axhf, axhp, axhv] = plotSameTrial_y(obj)
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             outcome = [obj.trials.outcome];
%             % fTH now only have 2 vals
%             % tarL now only have 2 vals
%             axhf = figure(); % force figure
%             set(axhf, 'Visible', 'off');
%             axhp = figure(); % position figure
%             set(axhp, 'Visible', 'off');
%             axhv = figure(); % velocity figure
%             set(axhv, 'Visible', 'off');
%             for fTH_i = 1:length(all_fTH)
%                 fprintf('\n plotting fTH:%2.1f  ', all_fTH(fTH_i));
%                 for tarL_i = 1:length(all_tarL)
%                     fprintf('target: %0.2f', all_tarL(tarL_i));
%                     trials_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)));
%                     trial_num = sum(trials_idx);
%                     trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) & outcome==1); %succeed
%                     %trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) ); %all trials
%                     col_i = (fTH_i-1)*length(all_fTH) + tarL_i;
%                     if col_i>=7
%                         col_i = mod(col_i,8);
%                     end
%                     for trial_i = find(trialsS_idx)
%                         % plot force
%                         % data
%                         try
%                             force = obj.trials(trial_i).force_h;
%                             time = obj.trials(trial_i).force_t;
%                         catch
%                             force = obj.trials(trial_i).force;
%                             time = obj.trials(trial_i).time;
%                         end
%                         % figure
%                         set(0, 'CurrentFigure', axhf);
%                         % x, y seperately
%                         hold on; plot(time, force(2,:), '.', 'COlor', obj.col_vec(col_i,:));
%                         % plot position
%                         % data
%                         try
%                             position = obj.trials(trial_i).position_h;
%                             time = obj.trials(trial_i).position_t;
%                         catch
%                             position = obj.trials(trial_i).position;
%                             time = obj.trials(trial_i).time;
%                         end
%                         set(0, 'CurrentFigure', axhp);
%                         % x, y seperately
%                         hold on; plot(time, position(2,:), '.', 'Color', obj.col_vec(col_i,:));
%                         try
%                             velocity = obj.trials(trial_i).velocity_h;
%                             time = obj.trials(trial_i).position_t;
%                         catch
%                             velocity = obj.trials(trial_i).velocity;
%                             time = obj.trials(trial_i).time;
%                         end
%                         set(0, 'CurrentFigure', axhv);
%                         hold on; plot(time, velocity(2,:), '.', 'Color', obj.col_vec(col_i,:));
%                     end
%                     figure(axhf); title('Force');
%                     figure(axhp); title('Position');
%                     figure(axhv); title('Velocity');
%                 end
%             end
%             
%             xrangeF = [-0.2, 0.4];
%             figure(axhf); 
%             title('force data y');
%             xlim(xrangeF);
%             
%             xrangeP = [-0.2, 0.4];
%             figure(axhp);
%             title('robot position y');
%             xlim(xrangeP);
%             
%             xrangeV = [-0.2, 0.4];
%             figure(axhv);
%             title('robot velocity y');
%             xlim(xrangeV);
%             
%             set(axhf, 'Visible', 'on');
%             set(axhp, 'Visible', 'on');
%             set(axhv, 'Visible', 'on');
%         end
%         function [axhf, axhp, axhv] = plotSameTrial_failure(obj)
%             all_fTH = unique([obj.trials.fTh]);
%             all_fTH = all_fTH(~isnan(all_fTH));
%             all_tarL = unique([obj.trials.tarL]);
%             all_tarL = all_tarL(~isnan(all_tarL));
%             outcome = [obj.trials.outcome];
%             % fTH now only have 2 vals
%             % tarL now only have 2 vals
%             axhf = figure(); % force figure
%             set(axhf, 'Visible', 'off');
%             axhp = figure(); % position figure
%             set(axhp, 'Visible', 'off');
%             axhv = figure(); % velocity figure
%             set(axhv, 'Visible', 'off');
%             for fTH_i = 1:length(all_fTH)
%                 fprintf('\n plotting fTH:%2.1f  ', all_fTH(fTH_i));
%                 for tarL_i = 1:length(all_tarL)
%                     fprintf('target: %0.2f', all_tarL(tarL_i));
%                     trials_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)));
%                     trial_num = sum(trials_idx);
%                     trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) & outcome==0); %failed
%                     %trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) ); %all trials
%                     col_i = (fTH_i-1)*length(all_fTH) + tarL_i;
%                     if col_i>=7
%                         col_i = mod(col_i,8);
%                     end
%                     for trial_i = find(trialsS_idx)
%                         % plot force
%                         % data
%                         try
%                             force = obj.trials(trial_i).force_h;
%                             time = obj.trials(trial_i).force_t;
%                         catch
%                             force = obj.trials(trial_i).force;
%                             time = obj.trials(trial_i).time;
%                         end
%                         % figure
%                         set(0, 'CurrentFigure', axhf);
%                         % x, y seperately
%                         subplot(2,1,1); hold on; plot(time, force(1,:), 'Color', obj.col_vec(col_i,:));
%                         subplot(2,1,2); hold on; plot(time, force(2,:), 'COlor', obj.col_vec(col_i,:));
%                         % plot position
%                         % data
%                         try
%                             position = obj.trials(trial_i).position_h;
%                             time = obj.trials(trial_i).position_t;
%                         catch
%                             position = obj.trials(trial_i).position;
%                             time = obj.trials(trial_i).time;
%                         end
%                         set(0, 'CurrentFigure', axhp);
%                         % x, y seperately
%                         subplot(2,1,1); hold on; plot(time, position(1,:), 'Color', obj.col_vec(col_i,:));
%                         subplot(2,1,2); hold on; plot(time, position(2,:), 'Color', obj.col_vec(col_i,:));
%                         try
%                             velocity = obj.trials(trial_i).velocity_h;
%                             time = obj.trials(trial_i).position_t;
%                         catch
%                             velocity = obj.trials(trial_i).velocity;
%                             time = obj.trials(trial_i).time;
%                         end
%                         set(0, 'CurrentFigure', axhv);
%                         subplot(2,1,1); hold on; plot(time, velocity(1,:), 'Color', obj.col_vec(col_i,:));
%                         subplot(2,1,2); hold on; plot(time, velocity(2,:), 'Color', obj.col_vec(col_i,:));
%                     end
%                 end
%             end
%             
%             xrangeF = [-0.5, 0.5];
%             figure(axhf); 
%             subplot(2,1,1);
%             title('force data x');
%             xlim(xrangeF);
%             subplot(2,1,2);
%             title('force data y');
%             xlim(xrangeF);
%             
%             xrangeP = [-0.2, 0.8];
%             figure(axhp);
%             subplot(2,1,1);
%             title('robot position x');
%             xlim(xrangeP);
%             subplot(2,1,2);
%             title('robot position y');
%             xlim(xrangeP);
%             
%             xrangeV = [-0.2, 0.8];
%             figure(axhv);
%             subplot(2,1,1);
%             title('robot velocity x');
%             xlim(xrangeV);
%             subplot(2,1,2);
%             title('robot velocity y');
%             xlim(xrangeV);
%             
%             set(axhf, 'Visible', 'on');
%             set(axhp, 'Visible', 'on');
%             set(axhv, 'Visible', 'on');
%         end
%         function axh = plotEndPointvsX0(obj)
%             combos = unique([obj.trials.comboTT]);
%             combos = combos(~isnan(combos));
%             axh = figure(); hold on;
%             hist_h = zeros(1, length(combos));
%             for i = 1:length(combos)
%                 combo_i = combos(i);
%                 trials_idx = ([obj.trials.comboTT] == combo_i);
%                 x0 = [obj.trials(trials_idx).pred_x0];
%                 hist_h(i) = histogram(x0);
%             end
%             tarL_all = [];
%             fTH_all  = [];
%             for combo_i = combos
%                 trials_idx = ([obj.trials.comboTT] == combo_i);
%                 tarL = unique([obj.trials(trials_idx).tarL]);
%                 tarL_all = [tarL_all tarL];
%                 fTh = unique([obj.trials(trials_idx).fTh]);
%                 fTH_all = [fTH_all fTh];
%             end
%             legend(hist_h, [num2str(tarL_all(1)) 'm, ' num2str(fTH_all(1)) 'N'], ...
%                 [num2str(tarL_all(2)) 'm, ' num2str(fTH_all(2)) 'N'], ...
%                 [num2str(tarL_all(3)) 'm, ' num2str(fTH_all(3)) 'N'], ...
%                 [num2str(tarL_all(4)) 'm, ' num2str(fTH_all(4)) 'N']);
%             title(['Session' num2str(obj.ssnum) 'x0 prediction with different task setting']);
%             ylabel('count');
%             xlabel('x0 prediction');
%         end
%         function axh = plotPredStiffness_LinDev_hist(obj)
%             combos = unique([obj.trials.comboTT]);
%             combos = combos(~isnan(combos));
%             axh = figure(); hold on;
%             hist_h = zeros(1, length(combos));
%             for i = 1:length(combos)
%                 combo_i = combos(i);
%                 trials_idx = ([obj.trials.comboTT] == combo_i);
%                 pred_K = [obj.trials(trials_idx).pred_K];
%                 hist_h(i) = histogram(pred_K);
%             end
%             tarL_all = [];
%             fTH_all  = [];
%             for combo_i = combos
%                 trials_idx = ([obj.trials.comboTT] == combo_i &...
%                     [obj.trials.outcome] == 1);
%                 tarL = unique([obj.trials(trials_idx).tarL]);
%                 tarL_all = [tarL_all tarL];
%                 fTh = unique([obj.trials(trials_idx).fTh]);
%                 fTH_all = [fTH_all fTh];
%             end
%             legend(hist_h, [num2str(tarL_all(1)) 'm, ' num2str(fTH_all(1)) 'N'], ...
%                 [num2str(tarL_all(2)) 'm, ' num2str(fTH_all(2)) 'N'], ...
%                 [num2str(tarL_all(3)) 'm, ' num2str(fTH_all(3)) 'N'], ...
%                 [num2str(tarL_all(4)) 'm, ' num2str(fTH_all(4)) 'N']);
%             title(['Session' num2str(obj.ssnum) ' K prediction with different task setting']);
%             ylabel('count');
%             xlabel('K prediction');
%         end
%         function axh = plotPredStiffness_LinDev_box(obj)
%             combos = unique([obj.trials.comboTT]);
%             combos = combos(~isnan(combos));
%             axh = figure(); hold on;
%             combo_all = [];
%             predK_all = [];
%             for i = 1:length(combos)
%                 combo_i = combos(i);
%                 trials_idx = ([obj.trials.comboTT] == combo_i);
%                 pred_K = [obj.trials(trials_idx).pred_K];
%                 combo_all = [combo_all repmat(combo_i, 1, length(pred_K))];
%                 predK_all = [predK_all pred_K];
%             end
%             tarL_all = [];
%             fTH_all  = [];
%             boxchart(combo_all, predK_all);
%             for combo_i = combos
%                 trials_idx = ([obj.trials.comboTT] == combo_i &...
%                     [obj.trials.outcome] == 1);
%                 tarL = unique([obj.trials(trials_idx).tarL]);
%                 tarL_all = [tarL_all tarL];
%                 fTh = unique([obj.trials(trials_idx).fTh]);
%                 fTH_all = [fTH_all fTh];
%             end
%             xticks([combos]);
%             xticklabels({[num2str(tarL_all(1)) 'm' num2str(fTH_all(1)) 'N'],...
%                 [num2str(tarL_all(2)) 'm' num2str(fTH_all(2)) 'N'],...
%                 [num2str(tarL_all(3)) 'm' num2str(fTH_all(3)) 'N'],...
%                 [num2str(tarL_all(4)) 'm' num2str(fTH_all(4)) 'N']});
%             title(['Session' num2str(obj.ssnum) ' K prediction with different task setting']);
%             ylabel('Stiffness N/m');
%             xlabel('task settings');
%         end
%         function axh = plotPredX0_LinDev_box(obj)
%             combos = unique([obj.trials.comboTT]);
%             combos = combos(~isnan(combos));
%             axh = figure(); hold on;
%             combo_all = [];
%             predX0_all = [];
%             for i = 1:length(combos)
%                 combo_i = combos(i);
%                 trials_idx = ([obj.trials.comboTT] == combo_i);
%                 pred_X0 = [obj.trials(trials_idx).pred_x0];
%                 combo_all = [combo_all repmat(combo_i, 1, length(pred_X0))];
%                 predX0_all = [predX0_all pred_X0];
%             end
%             tarL_all = [];
%             fTH_all  = [];
%             boxchart(combo_all, predX0_all);
%             for combo_i = combos
%                 trials_idx = ([obj.trials.comboTT] == combo_i &...
%                     [obj.trials.outcome] == 1);
%                 tarL = unique([obj.trials(trials_idx).tarL]);
%                 tarL_all = [tarL_all tarL];
%                 fTh = unique([obj.trials(trials_idx).fTh]);
%                 fTH_all = [fTH_all fTh];
%             end
%             ylim([0 0.2]);
%             xticks([combos]);
%             xticklabels({[num2str(tarL_all(1)) 'm' num2str(fTH_all(1)) 'N'],...
%                 [num2str(tarL_all(2)) 'm' num2str(fTH_all(2)) 'N'],...
%                 [num2str(tarL_all(3)) 'm' num2str(fTH_all(3)) 'N'],...
%                 [num2str(tarL_all(4)) 'm' num2str(fTH_all(4)) 'N']});
%             title(['Session' num2str(obj.ssnum) ' x0 prediction with different task setting']);
%             ylabel('Eq m');
%             xlabel('task settings');
%         end
%         function axh = plotForceVecBeforeRelease(obj, axh, rot)
%             if nargin == 1
%                 axh = figure();
%             elseif nargin == 2
%                 rot = 0; % rotation from the robot to the subject 
%             end
%             % data
%             sucessful_idx = find([obj.trials.outcome]);
%             forceVectors = zeros(3, length(sucessful_idx));
%             for trial_i = 1:length(sucessful_idx)
%                 forceVectors(:,trial_i) = ...
%                     obj.trials(sucessful_idx(trial_i)).getforceVecBeforeRelease();
%             end
%             forceVectors = [[cos(rot), -sin(rot); sin(rot), cos(rot)] * ...
%                 forceVectors(1:2,:); ...
%                 forceVectors(3,:)]; % rotation;
%             
%             % plot
%             figure(axh);
%             for trial_i = 1:length(sucessful_idx)
%                 line([0 forceVectors(1, trial_i)], [0 forceVectors(2, trial_i)], ...
%                     'lineWidth', 3, 'color', 'b');
%             end
%             xlabel('x component'); 
%             ylabel('y component');
%             title('force vector before release');
%         end
%         %%% plot overlapped pert responses
%         function axh = plotStepPertResponse_raw(obj, axh, color_arr)
%             % axh = plotStepPertResponse_raw % plot the raw response of
%             % step perturbation, of all trials in this session 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%              
%             % check input 
%             if ~exist('axh', 'var')
%                 axh = 0;
%             end
%             if ~exist('color_arr', 'var')
%                 color_arr = [0 0 0];
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             obj = updatePertEachTrial(obj);
%             %trials_pert = find([obj.trials(:).ifpert]);
%             %trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             trials_pert = obj.getPerturbedTrialIdx();
%             % plot
%             
%             % figure properties
%             
%             % for each trial
%             for trial_i = trials_pert
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).position_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 %if (isempty(obj.trials(trial_i).pert_t_bgn))
%                 %    time0 = pert_t_last;
%                 %end
%                 resp_p = obj.trials(trial_i).position_h(2,:);
%                 resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 position_offset = 0.482;
%                 %position_offset = mean(obj.trials(trial_i).position_h(2,...
%                 %    find(time>time0-0.1 & time<time0)));
%                 resp_p_net= resp_p - position_offset; %obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time-time0, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                 else
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                 end
%                 pert_t_last = time0;
%             end
%             if flag_stak == 0
%                 xlim([-0.2, 0.8]);
%                 ylabel('endpoint position (m)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.2, 0.6]);
%                 %ylim([-0.015, 0.015]);
%                 ylabel('endpoint position (m)');
%                 xlabel('time'); 
%                 title(['session' num2str(obj.ssnum) ' perturbation']);
%             end
%             
%         end
%         function axh = plotStepPertResponse_rawV(obj, axh, color_arr)
%             % axh = plotStepPertResponse_rawV % plot the raw velocity response of
%             % step perturbation, of all trials in this session 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%              
%             % check input 
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(color_arr)
%                 color_arr = 0;
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             obj = updatePertEachTrial(obj);
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             
%             % plot
%             
%             % figure properties
%             
%             % for each trial
%             for trial_i = trials_pert
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).position_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 if (isempty(obj.trials(trial_i).pert_t_bgn))
%                     time0 = pert_t_last;
%                 end
%                 resp_p = obj.trials(trial_i).position_h(2,:);
%                 resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 position_offset = mean(obj.trials(trial_i).position_h(2,...
%                     find(time>time0-0.1 & time<time0)));
%                 resp_p_net= resp_p;% - position_offset; %obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time-time0, resp_v, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                 else
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time-time0, resp_v);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                 end
%                 pert_t_last = time0;
%             end
%             if flag_stak == 0
%                 xlim([-0.2, 0.8]);
%                 ylabel('endpoint velocity (m/s)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.2, 0.6]);
%                 %ylim([-0.015, 0.015]);
%                 ylabel('endpoint velocity (m/s)');
%                 xlabel('time'); 
%                 title(['session' num2str(obj.ssnum) ' perturbation']);
%             end
%             
%         end
%         function axh = plotStepPertResponse_rawF(obj, axh, color_arr)
%             % axh = plotStepPertResponse_rawF % plot the raw response of
%             % step perturbation, of all trials in this session 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%              
%             % check input 
%             if ~exist('axh', 'var')
%                 axh = 0;
%             end
%             if ~exist('color_arr', 'var')
%                 color_arr = [0 0 0];
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             obj = updatePertEachTrial(obj);
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             % plot
%             
%             % figure properties
%             
%             % for each trial
%             for trial_i = trials_pert
%                 if trial_i == 1
%                     continue;
%                 end
%                 % make the time aligned for the perturbation
%                 %time = obj.trials(trial_i).position_t;
%                 time = obj.trials(trial_i).force_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 %resp_p = obj.trials(trial_i).position_h(2,:);
%                 %resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 %resp_f = obj.trials(trial_i).force_h(2,:);
%                 resp_f = obj.trials(trial_i).force_h(2,:);
%                 %resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time-time0, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                 else
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time-time0, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                 end
%             end
%             if flag_stak == 0
%                 xlim([-0.2, 0.8]);
%                 ylabel('force (N)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.2, 0.6]);
%                 %ylim([-0.015, 0.015]);
%                 ylabel('force (N)');
%                 xlabel('time');
%                 title(['session' num2str(obj.ssnum) ' perturbation']);
%             end
%             
%         end       
%         function axh = plotReleaseResponse_rawF(obj, axh, color_arr)
%             % axh = plotReleaseResponse_rawF % plot the raw response of
%             % ballistic release, in terms of force. of all trials in this session 
%             % Assuming all trials are at the same task condition
%             % It requires the same stiffness so that response magnitudes are
%             % the same.
%              
%             % check input 
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(color_arr)
%                 color_arr = 0;
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             %obj = updatePertEachTrial(obj);
%             %trials_pert = find([obj.trials(:).ifpert]);
%             trials_all = find([obj.trials(:).outcome]); % outcome not 0
%             % plot
%             
%             % figure properties
%             
%             % for each trial
%             for trial_i = trials_all
%                 if trial_i == 1
%                     continue;
%                 end
%                 % make the time aligned for the perturbation
%                 %time = obj.trials(trial_i).position_t;
%                 time = obj.trials(trial_i).force_t;
%                 %time0 = obj.trials(trial_i).pert_t_bgn;
%                 %resp_p = obj.trials(trial_i).position_h(2,:);
%                 %resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 %resp_f = obj.trials(trial_i).force_h(2,:);
%                 resp_f = obj.trials(trial_i).force_h(2,:);
%                 %resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                 else
%                     %plot(time-time0, resp_p, 'color', color_arr);
%                     plot(time, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                 end
%             end
%             if flag_stak == 0
%                 xlim([-0.1, 1]);
%                 ylabel('force (N)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.1, 1]);
%                 %ylim([-0.015, 0.015]);
%                 ylabel('force (N)');
%                 xlabel('time');
%                 title(['session' num2str(obj.ssnum) ' perturbation']);
%             end
%             
%         end
%         function axh = plotReleaseResponse_rawP(obj, axh, color_arr)
%             % axh = plotReleaseResponse_rawP % plot the raw response of
%             % ballistic release, in terms of position of all trials in this session 
%             % Assuming all trials are at the same task condition
%             % It requires the same stiffness so that response magnitudes are
%             % the same.
%              
%             % check input 
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(color_arr)
%                 color_arr = 0;
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             %obj = updatePertEachTrial(obj);
%             %trials_pert = find([obj.trials(:).ifpert]);
%             trials_all = find([obj.trials(:).outcome]); % outcome not 0
%             % plot
%             
%             % figure properties
%             
%             % for each trial
%             for trial_i = trials_all
%                 if trial_i == 1
%                     continue;
%                 end
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).position_t;
%                 %time = obj.trials(trial_i).force_t;
%                 %time0 = obj.trials(trial_i).pert_t_bgn;
%                 resp_p = obj.trials(trial_i).position_h(2,:);
%                 %resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 %resp_f = obj.trials(trial_i).force_h(2,:);
%                 %resp_f = obj.trials(trial_i).force_h(2,:);
%                 position_offset = 0.482;
%                 resp_p_net= resp_p - position_offset; %obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     %plot(time, resp_p, 'color', color_arr);
%                     %plot(time, resp_f, 'color', color_arr);
%                     plot(time, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                 else
%                     %\plot(time, resp_p, 'color', color_arr);
%                     %plot(time, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                 end
%             end
%             if flag_stak == 0
%                 xlim([-0.1, 1]);
%                 ylabel('position (m)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.1, 1]);
%                 %ylim([-0.015, 0.015]);
%                 ylabel('position (m)');
%                 xlabel('time');
%                 title(['session' num2str(obj.ssnum) ' perturbation']);
%             end
%             
%         end
%         function axh = plotReleaseResponse_rawV(obj, axh, color_arr)
%             % axh = plotReleaseResponse_rawV % plot the raw response of
%             % ballistic release, in terms of velocity of all trials in this session 
%             % Assuming all trials are at the same task condition
%             % It requires the same stiffness so that response magnitudes are
%             % the same.
%              
%             % check input 
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(color_arr)
%                 color_arr = 0;
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             %obj = updatePertEachTrial(obj);
%             %trials_pert = find([obj.trials(:).ifpert]);
%             trials_all = find([obj.trials(:).outcome]); % outcome not 0
%             % plot
%             
%             % figure properties
%             
%             % for each trial
%             for trial_i = trials_all
%                 if trial_i == 1
%                     continue;
%                 end
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).position_t;
%                 %time = obj.trials(trial_i).force_t;
%                 %time0 = obj.trials(trial_i).pert_t_bgn;
%                 %resp_p = obj.trials(trial_i).position_h(2,:);
%                 resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 %resp_f = obj.trials(trial_i).force_h(2,:);
%                 %resp_f = obj.trials(trial_i).force_h(2,:);
%                 %resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     %plot(time, resp_p, 'color', color_arr);
%                     %plot(time, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net, 'color', color_arr);
%                     plot(time, resp_v, 'color', color_arr); %velocity
%                 else
%                     %\plot(time, resp_p, 'color', color_arr);
%                     %plot(time, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net);
%                     plot(time, resp_v, 'color', color_arr);
%                 end
%             end
%             if flag_stak == 0
%                 xlim([-0.1, 1]);
%                 ylabel('position (m)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.1, 1]);
%                 %ylim([-0.015, 0.015]);
%                 ylabel('position (m)');
%                 xlabel('time');
%                 title(['session' num2str(obj.ssnum) ' perturbation']);
%             end
%             
%         end
%         function [axh, val, lnh] = plotStepPertResponse_raw_subavg(obj, axh, color_arr)
%             % [axh, val, lnh] = plotStepPertResponse_raw % plot the raw response of
%             % step perturbation, of all trials in this session 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%             %   axh: the axis handle, for future plot. 
%             %   val: the value (avg) of peak, val{1} positive, val{2}
%             %           negative
%             %   lnh: the line handle, for the futrue legend on
%              
%             % check input 
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(color_arr)
%                 color_arr = 0;
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             % obj = updatePertEachTrial(obj);
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             trials_fin  = find([obj.trials(:).outcome] == 1);
%             pert_pltavg_pos = []; % put values here, pltavg: pleateau avg
%             pert_pltavg_neg = []; 
%             % for each trial
%             trials_list = intersect(trials_pert, trials_fin);
%             for trial_i = trials_list %trials_pert
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).position_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 resp_p = obj.trials(trial_i).position_h(2,:); % resting position
%                 if min(time-time0) < -0.1
%                     time_idx = (time-time0) > -0.1 & (time-time0) < 0;
%                 else
%                     [~, time_idx] = min(abs(time - time0));
%                 end
%                 resp_p0= mean(obj.trials(trial_i).position_h(2,time_idx));
%                 resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 
%                 time_pleateu_idx = time-time0>0.2 & time-time0 < 0.3;
%                 pert_pleateu_avg = mean(obj.trials(trial_i).position_h(2,time_pleateu_idx)) - resp_p0;
%                 if pert_pleateu_avg > 0
%                     pert_pltavg_pos = [pert_pltavg_pos, pert_pleateu_avg];
%                 else
%                     pert_pltavg_neg = [pert_pltavg_neg, pert_pleateu_avg];
%                 end
%                 %resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     plot(time-time0, resp_p-resp_p0, 'color', color_arr);
%                     %plot(time-time0, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_p-resp_p0, 'color', color_arr);
%                     end
%                 else
%                     plot(time-time0, resp_p-resp_p0, 'color', color_arr);
%                     %plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_p-resp_p0, 'color', color_arr);
%                     end
%                 end
%             end
%             if flag_stak == 0
%                 xlim([-0.2, 0.8]);
%                 ylabel('endpoint position (m)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.2, 0.6]);
%                 %ylim([-0.015, 0.015]);
%             end
%             
%             % send out the peak values
%             
%             val{1} = pert_pltavg_pos;
%             val{2} = pert_pltavg_neg;
%         end
%         function [axh, val, lnh] = plotStepPertResponse_raw_pertfce(obj, axh, color_arr)
%             % [axh, val, lnh] = plotStepPertResponse_raw % plot the raw response of
%             % step perturbation, of all trials in this session 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%             %   axh: the axis handle, for future plot. 
%             %   val: the value (avg) of peak, val{1} positive, val{2}
%             %           negative
%             %   lnh: the line handle, for the futrue legend on
%              
%             % check input 
%             if ~exist('axh','var')
%                 axh = 0;
%             end
%             if ~exist('color_arr', 'var')
%                 color_arr = [0 0 0];
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             % obj = updatePertEachTrial(obj);
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             trials_fin  = find([obj.trials(:).outcome] == 1);
%             pert_pltavg_pos = []; % put values here, pltavg: pleateau avg
%             pert_pltavg_neg = []; 
%             % for each trial
%             trials_list = intersect(trials_pert, trials_fin);
%             for trial_i = trials_list %trials_pert
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).position_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 if (isempty(time0))
%                     time0 = pert_t_last;
%                 end
%                 resp_fp = obj.trials(trial_i).pertfce_h(2,:); % resting position
%                 if min(time-time0) < -0.1
%                     time_idx = (time-time0) > -0.1 & (time-time0) < 0;
%                 else
%                     [~, time_idx] = min(abs(time - time0));
%                 end
%                 resp_p0= mean(obj.trials(trial_i).position_h(2,time_idx));
%                 resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 
%                 time_pleateu_idx = time-time0>0.2 & time-time0 < 0.3;
%                 pert_pleateu_avg = mean(obj.trials(trial_i).position_h(2,time_pleateu_idx)) - resp_p0;
%                 if pert_pleateu_avg > 0
%                     pert_pltavg_pos = [pert_pltavg_pos, pert_pleateu_avg];
%                 else
%                     pert_pltavg_neg = [pert_pltavg_neg, pert_pleateu_avg];
%                 end
%                 %resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     plot(time-time0, resp_fp, 'Color', color_arr);
%                     %plot(time-time0, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_fp, 'color', color_arr);
%                     end
%                 else
%                     plot(time-time0, resp_fp, 'color', color_arr);
%                     %plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_fp, 'color', color_arr);
%                     end
%                 end
%                 pert_t_last = time0;
%             end
%             if flag_stak == 0
%                 xlim([-0.2, 0.8]);
%                 ylabel('perturbation force (N)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.2, 0.6]);
%                 %ylim([-0.015, 0.015]);
%             end
%             
%             % send out the peak values
%             
%             val{1} = pert_pltavg_pos;
%             val{2} = pert_pltavg_neg;
%         end
%         function [axh, val, lnh] = plotStepPertResponse_rawFce(obj, axh, color_arr, low_pass_freq)
%             % [axh, val, lnh] = plotStepPertResponse_raw % plot the raw response of
%             % step perturbation, of all trials in this session 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%             %   axh: the axis handle, for future plot. 
%             %   val: the value (avg) of peak, val{1} positive, val{2}
%             %           negative
%             %   lnh: the line handle, for the futrue legend on
%              
%             % check input 
%             if isempty(axh)
%                 axh = 0;
%             end
%             if isempty(color_arr)
%                 color_arr = 0;
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = 0;
%             else
%                 ifcolor = 1;
%             end
%             if ~exist('low_pass_freq', 'var')
%                 ifLowpass = false;
%             else
%                 ifLowpass = true;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             % obj = updatePertEachTrial(obj);
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             trials_fin  = find([obj.trials(:).outcome] == 1);
%             pert_pltavg_pos = []; % put values here, pltavg: pleateau avg
%             pert_pltavg_neg = []; 
%             % for each trial
%             trials_list = intersect(trials_pert, trials_fin);
%             for trial_i = trials_list %trials_pert
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).force_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 if (isempty(obj.trials(trial_i).pert_t_bgn))
%                     continue;
%                 end
%                 resp_f = obj.trials(trial_i).force_h(2,:); % resting position
%                 time_idx = (time-time0) > -0.1 & (time-time0) < 0;
%                 resp_f0= mean(obj.trials(trial_i).force_h(2,time_idx)); % may need change
%                 resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 freq = 1/mean(diff(time));
%                 if (ifLowpass)
%                     resp_f = smooth(resp_f, floor(freq/low_pass_freq));
%                 end
%                 
%                 time_pleateu_idx = time-time0>0.2 & time-time0 < 0.3;
%                 pert_pleateu_avg = mean(obj.trials(trial_i).force_h(2,time_pleateu_idx)) - resp_f0;
%                 if pert_pleateu_avg > 0
%                     pert_pltavg_pos = [pert_pltavg_pos, pert_pleateu_avg];
%                 else
%                     pert_pltavg_neg = [pert_pltavg_neg, pert_pleateu_avg];
%                 end
%                 %resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     plot(time-time0, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_f, 'color', color_arr);
%                     end
%                 else
%                     plot(time-time0, resp_f, 'color', color_arr);
%                     %plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_f, 'color', color_arr);
%                     end
%                 end
%             end
%             if flag_stak == 0
%                 xlim([-0.2, 0.8]);
%                 ylabel('endpoint position (m)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.2, 0.6]);
%                 %ylim([-0.015, 0.015]);
%             end
%             
%             % send out the peak values
%             
%             val{1} = pert_pltavg_pos;
%             val{2} = pert_pltavg_neg;
%         end
%         function [axh, val, lnh] = plotStepPertResponse_rawFce_subavg(obj, axh, color_arr, low_pass_freq)
%             % [axh, val, lnh] = plotStepPertResponse_raw % plot the raw response of
%             % step perturbation, of all trials in this session 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%             %   axh: the axis handle, for future plot. 
%             %   val: the value (avg) of peak, val{1} positive, val{2}
%             %           negative
%             %   lnh: the line handle, for the futrue legend on
%              
%             % check input 
%             if ~exist('axh', 'var')
%                 axh = 0;
%             end
%             if ~exist('color_arr','var')
%                 color_arr = 0;
%             end
%             if length(color_arr) ~= 3
%                 ifcolor = [0 0 0];
%             else
%                 ifcolor = 1;
%             end
%             if ~exist('low_pass_freq', 'var')
%                 ifLowpass = false;
%             else
%                 ifLowpass = true;
%             end
%             if isa(axh, 'matlab.ui.Figure')
%                 axh = figure(axh); hold on; % stack
%                 flag_stak = 1; 
%             elseif isa(axh, 'matlab.graphics.axis.Axes')
%                 subplot(axh); hold on;
%                 flag_stak = 1;
%             else
%                 axh = figure(); hold on;
%                 flag_stak = 0;
%             end
%             % list all trials that being perturbed
%             % obj = updatePertEachTrial(obj);
%             trials_pert = find([obj.trials(:).ifpert]);
%             trials_pert = setdiff(trials_pert, 1); % remove first trial as unstable
%             trials_fin  = find([obj.trials(:).outcome] == 1);
%             pert_pltavg_pos = []; % put values here, pltavg: pleateau avg
%             pert_pltavg_neg = []; 
%             % for each trial
%             trials_list = intersect(trials_pert, trials_fin);
%             for trial_i = trials_list %trials_pert
%                 % make the time aligned for the perturbation
%                 time = obj.trials(trial_i).force_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 if (isempty(obj.trials(trial_i).pert_t_bgn))
%                     time0 = pert_t_last;
%                 end
%                 resp_f = obj.trials(trial_i).force_h(2,:); % resting position
%                 time_idx = (time-time0) > -0.1 & (time-time0) < 0;
%                 resp_f0= mean(obj.trials(trial_i).force_h(2,time_idx)); % may need change
%                 resp_v = obj.trials(trial_i).velocity_h(2,:);
%                 freq = 1/mean(diff(time));
%                 if (ifLowpass)
%                     resp_f = smooth(resp_f, floor(freq/low_pass_freq));
%                 end
%                 
%                 time_pleateu_idx = time-time0>0.2 & time-time0 < 0.3;
%                 pert_pleateu_avg = mean(obj.trials(trial_i).force_h(2,time_pleateu_idx)) - resp_f0;
%                 if pert_pleateu_avg > 0
%                     pert_pltavg_pos = [pert_pltavg_pos, pert_pleateu_avg];
%                 else
%                     pert_pltavg_neg = [pert_pltavg_neg, pert_pleateu_avg];
%                 end
%                 %resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 if ifcolor == 1
%                     plot(time-time0, resp_f-resp_f0, 'color', color_arr);
%                     %plot(time-time0, resp_p_net, 'color', color_arr);
%                     %plot(time-time0, resp_v, 'color', color_arr); %velocity
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_f-resp_f0, 'color', color_arr);
%                     end
%                 else
%                     plot(time-time0, resp_f-resp_f0);
%                     %plot(time-time0, resp_p_net);
%                     %plot(time-time0, resp_v, 'color', color_arr);
%                     if trial_i == trials_list(1)
%                         lnh = plot(time-time0, resp_f-resp_f0);
%                     end
%                 end
%                 pert_t_last = time0;
%             end
%             if flag_stak == 0
%                 xlim([-0.2, 0.8]);
%                 ylabel('sensored force (N)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             else
%                 xlim([-0.2, 0.6]);
%                 %ylim([-0.015, 0.015]);
%             end
%             
%             % send out the peak values
%             
%             val{1} = pert_pltavg_pos;
%             val{2} = pert_pltavg_neg;
%         end
%         function axh = plotStepPertResponse_raw_inv(obj)
%             % axh = plotStepPertResponse_raw_inv % plot the raw response of
%             % step perturbation, of all trials in this session 
%             % Invert the negative perturbation. 
%             % Assuming all trials are at the same task condition (e.g.
%             % requires the same stiffness so that response magnitudes are
%             % the same.
%             
%             % check input 
%             % ...
%             % list all trials that being perturbed
%             % obj = updatePertEachTrial(obj);
%             trials_pert = obj.getPerturbedTrialIdx();
%             % plot
%             axh = figure(); hold on;
%             % figure properties
%             %color = ['rg']; %r- front; g-back
%             % for each trial
%             for trial_i = trials_pert
%                 % make the time aligned for the perturbation
%                 time_tmp = obj.trials(trial_i).position_t;
%                 time0 = obj.trials(trial_i).pert_t_bgn;
%                 pertidx = find(time_tmp == time0);
%                 resp_p = obj.trials(trial_i).position_h(2,:);
%                 resp_p_net= resp_p - obj.trials(trial_i).position_offset;
%                 pert_sign = sign(obj.trials(trial_i).pertfce_h(pertidx+1));
%                 if pert_sign == 1
%                     color_i = 'r';
%                 else
%                     color_i = 'g';
%                 end
%                 %plot each trial's perturbation response
%                 %plot(time-time0, resp_p);
%                 
%                 plot(time_tmp-time0, pert_sign*resp_p_net, 'color', color_i);
%             end
% 
%                 xlim([-0.2, 0.6]);
%                 ylabel('endpoint position (m)');
%                 xlabel('time');
%                 title(['step pert response for session' num2str(obj.ssnum)]);
%             
%         end
        
        % exception figures for specific sessions:
        function axh = plotRecordedEndPointPosition(obj)
            % for testing if the recorded endpoint position is the actual
            % endpoint position.  
            if obj.ssnum == 1931
            % Moving endpoint +y direction seperately +2cm and +4cm in a
            % series of stiffness values, and plot the endpoint position. 
            % see RSHJournal in 2020-02-19
            
            % see if right session
            
            % the timepoints
            idx_y2cm = ...
                [17518	30069	41191	51895	61208	
                76916	97865	108529	121239	132379	
                149379	160977	171482	182723	195245
                209273	220054	230728	244754	258236	
                274799	286833	299005	315583	329559	
                345472	359380	372882	386543	400317	];
            idx_y2cm0= ...
                [15169	27988	37924	49741	58777	
                71535	94548	106463	118717	130513	
                146341	158098	169219	180817	193081
                207714	216974	228342	241758	255598	
                271838	284246	296775	312843	327006	
                342348	356617	370249	383523	398407	];
            idx_y4cm = ...
                [24318  35789   46173   55960    67382 
                86651	103401	114704  126840  138965 
                154638	166043	176801	188376	200769 
                214114	225145	237655	251155	265044 
                280841  292825	306851  322512	336332	
                353173	366509	379671	394934	409571	];
            idx_y4cm0 = ...
                [19720  32640   42790   53780   62890 
                80827	101112  111708  122994  135267 
                152425	163047	174486	185108	197908 
                212059	222592	233587	247945	262185 
                277513	289391	302883	318693	333144 
                349605	363426	376401	390174	406454 ];
            Kx0 = [0, 500, 1000, 1500, 2000, 2500];
            Kx0_mat = repmat(Kx0, 5, 1);
            % the y position
            y2cm = (obj.wamp_h(idx_y2cm',2) - obj.wamp_h(idx_y2cm0',2))/0.01;
            y4cm = (obj.wamp_h(idx_y4cm',2) - obj.wamp_h(idx_y4cm0',2))/0.01;
            y2cm_= (2 - y2cm);
            y4cm_= (4 - y4cm);
            % plot the point 
            axh = figure();
            hold on;
            dth1 = plot(Kx0_mat(:), y2cm_, '.', 'MarkerSize', 10);
            refline;
            dth2 = plot(Kx0_mat(:), y4cm_, '.', 'MarkerSize', 10);
            refline; 
            ax = gca;
            ax.XGrid = 'off';
            ax.YGrid = 'on';
            legend([dth2, dth1], 'x=4cm', 'x=2cm');
            xlim([-100, 2600]);
            %ylim([0, 0.05]);
            ylim([-0.5, 1.5]);
            xlabel('Kx N/m');
            ylabel('error cm');
            title('Measurement error');
            end
            
            if obj.ssnum == 1934
            % Moving endpoint +y direction seperately +0.5cm and +1cm in a
            % series of stiffness values, and plot the endpoint position. 
            % see RSHJournal in 2020-02-23
            
            % see if right session
            
            % the timepoints
            idx_y_5cm0= ...
                [[43166,53354,63125,72204,81101]	
                [89780,106318,116129,124344,132576]	
                [143097,161871,171277,179980,190026]
                [200306,217437,225979,234183,242343]	
                [250829,264622,272934,281923,292109]	
                [304025,314258,322328,331534,340066]	];
            idx_y_5cm = ...
                [[46052,55908,65660,74400,83165]	
                [91784,107985,117586,125554,134182]	
                [144755,163484,172414,181250,191748]
                [202181,218802,227254,234940,243731]	
                [253131,265462,274172,282870,293272]	
                [305771,315204,323240,332306,340835]	];
            idx_y1cm0 = ...
                [[49192,58610,67509,76634,85295] 
                [100851,111147,119954,128350,138854] 
                [156048,166584,175455,184630,194119] 
                [212717,221631,230274,238428,246602] 
                [259946,268730,277192,287910,296635] 
                [309366,318166,326597,335150,343811] ];
            idx_y1cm = ...
                [[51273,61318,69998,78522,87425] 
                [103086,113438,121879,130060,140598] 
                [158274,167181,176712,185685,195731] 
                [214212,222956,231080,239580,247831] 
                [261365,269693,277854,288737,297834]	
                [310857,319299,328166,336458,345126]	];
            Kx0 = [0, 500, 1000, 1500, 2000, 2500];
            Kx0_mat = repmat(Kx0, 5, 1);
            % the y position
            y2cm = (obj.wamp_h(idx_y_5cm',2) - obj.wamp_h(idx_y_5cm0',2))/0.01;
            y4cm = (obj.wamp_h(idx_y1cm',2) - obj.wamp_h(idx_y1cm0',2))/0.01;
            y2cm_= (0.5 - y2cm);
            y4cm_= (1 - y4cm);
            % plot the point 
            axh = figure();
            hold on;
            dth1 = plot(Kx0_mat(:), y2cm_, '.', 'MarkerSize', 10);
            refline;
            dth2 = plot(Kx0_mat(:), y4cm_, '.', 'MarkerSize', 10);
            refline; 
            ax = gca;
            ax.XGrid = 'off';
            ax.YGrid = 'on';
            legend([dth2, dth1], 'x=1cm', 'x=0.5cm');
            xlim([-100, 2600]);
            %ylim([0, 0.05]);
            ylim([-0.2, 0.5]);
            xlabel('Kx N/m');
            ylabel('error cm');
            title('Measurement error');
            end
            
            if obj.ssnum == 1935
            % Moving endpoint +y direction using force measurement in a
            % series of stiffness values, and plot the endpoint position. 
            % see RSHJournal in 2020-02-23
            
            % see if right session
            
            % the timepoints
            idx_y_5cm0= ...
                [[43166,53354,63125,72204,81101]	
                [89780,106318,116129,124344,132576]	
                [143097,161871,171277,179980,190026]
                [200306,217437,225979,234183,242343]	
                [250829,264622,272934,281923,292109]	
                [304025,314258,322328,331534,340066]	];
            idx_y_5cm = ...
                [[46052,55908,65660,74400,83165]	
                [91784,107985,117586,125554,134182]	
                [144755,163484,172414,181250,191748]
                [202181,218802,227254,234940,243731]	
                [253131,265462,274172,282870,293272]	
                [305771,315204,323240,332306,340835]	];
            idx_y1cm0 = ...
                [[49192,58610,67509,76634,85295] 
                [100851,111147,119954,128350,138854] 
                [156048,166584,175455,184630,194119] 
                [212717,221631,230274,238428,246602] 
                [259946,268730,277192,287910,296635] 
                [309366,318166,326597,335150,343811] ];
            idx_y1cm = ...
                [[51273,61318,69998,78522,87425] 
                [103086,113438,121879,130060,140598] 
                [158274,167181,176712,185685,195731] 
                [214212,222956,231080,239580,247831] 
                [261365,269693,277854,288737,297834]	
                [310857,319299,328166,336458,345126]	];
            Kx0 = [0, 500, 1000, 1500, 2000, 2500];
            Kx0_mat = repmat(Kx0, 5, 1);
            % the y position
            y2cm = (obj.wamp_h(idx_y_5cm',2) - obj.wamp_h(idx_y_5cm0',2))/0.01;
            y4cm = (obj.wamp_h(idx_y1cm',2) - obj.wamp_h(idx_y1cm0',2))/0.01;
            y2cm_= (0.5 - y2cm);
            y4cm_= (1 - y4cm);
            % plot the point 
            axh = figure();
            hold on;
            dth1 = plot(Kx0_mat(:), y2cm_, '.', 'MarkerSize', 10);
            refline;
            dth2 = plot(Kx0_mat(:), y4cm_, '.', 'MarkerSize', 10);
            refline; 
            ax = gca;
            ax.XGrid = 'off';
            ax.YGrid = 'on';
            legend([dth2, dth1], 'x=1cm', 'x=0.5cm');
            xlim([-100, 2600]);
            %ylim([0, 0.05]);
            ylim([-0.2, 0.5]);
            xlabel('Kx N/m');
            ylabel('error cm');
            title('Measurement error');
            end
            
            if obj.ssnum == 1937 % ........ remember to do it later today-cg, tell the difference of force 
            % Moving endpoint +y direction using force measurement in a
            % series of stiffness values, and plot the endpoint position. 
            % see RSHJournal in 2020-02-24
            
            % see if right session
            
            % the timepoints
            mark_points_pidx = [[98711,106238,115289,121453,128060,134232,140822,146860,154532,164123;170039,176956,183303,190006,195974,202767,207357,213104,218773,225426;229908,234749,240283,245042,251434,256594,262023,267078,273125,278855;283931,288132,294089,298720,303960,308201,314518,320204,326198,332360;336335,340619,346080,350330,355241,360779,365989,370401,375789,380547]];
            mark_points_fidx = [[131649,141811,153296,161952,171277,181191,190031,199230,208382,221735;230949,239667,249602,257913,267645,277198,283658,292322,298482,307948;314907,320869,328917,335411,344123,352020,360516,366502,376120,381122;389864,396727,404140,411256,418624,424603,432600,440484,448752,457617;463822,470138,475846,483106,489594,496613,504303,510249,517776,524549]];
            mark_point_p = obj.wamp_h(mark_points_pidx,2)-0.481;
            mark_point_f = obj.force_h(2,mark_points_fidx)+3.2; % add 3 newton bias
            mark_point_f_theoretic = ones(5,10)*16;
            figure();
            hold on;
            plot(mark_point_f(:), mark_point_p(:), '*');
            plot(mark_point_f_theoretic(:), mark_point_p(:), '*');
            title('force vs position');
            figure();
            hold on;
            title('position vs stiffness');
            plot(mark_point_p(:),mark_point_f(:)./mark_point_p(:), '*');
            plot( mark_point_p(:), mark_point_f_theoretic(:)./mark_point_p(:), '*');
            legend('measured', 'theoretical');
            
            ax = gca;
            ax.XGrid = 'off';
            ax.YGrid = 'on';
            xlabel('position cm');
            ylabel('Stiffness');
            title('Measurement error');
            end
        end
        
        %% update time from the blackrock hardware recording 
        function [obj flag] = updateTimeBlackRock(obj)
            % THIS FUNCTION IS FOR ALIGNED THE MESSAGE TIME 
            % For the WAM, replace the BURT_STATUS message time to the blackrock time
            % For the FT, replace the PFEM_TIME to the blackrock time 
            %clear; close all; clc;
            ss_num = obj.ssnum;
            fname_sync = sprintf('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/KingKongTSync.%05d.mat', ss_num);
            fname_intm = sprintf('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate/KingKong.%05d.mat', ss_num);
            if ~exist(fname_sync, 'file') 
                disp('ERROR: no TSync file saved, please check!');
                flag = -1; 
                return
            end
                
            if ~exist(fname_intm, 'file')
                disp('ERROR: no intermediate file saved, please check!');
                flag = -2;
                return
            end
            
            d = load(fname_sync);

            try
                %%% Use the saved synchronize data;
                dataTs = d.data;
            catch
                %%%% Data save failure in the TSync, use the saved temp data;
                i = 1;
                d.data.eventsL = [];
                d.data.eventsT = [];
                d.data.eventsTrials = [];
                
                while (isfield(d, ['datatmp' num2str(i)]))
                    datatmp_sync = eval(['d.datatmp' num2str(i)]);
                    d.data.eventsL = [d.data.eventsL datatmp_sync.eventsL];
                    d.data.eventsT = [d.data.eventsT datatmp_sync.eventsT];
                    d.data.eventsTrials = [d.data.eventsTrials datatmp_sync.eventsTrials];
                    i = i+1;
                end
                dataTs = d.data;
            end
            
                
            d = load(fname_intm);
            dataMsT = d.Data.QL.Data.TIME_SYNC;
            dataMsTh= d.Data.QL.Headers.TIME_SYNC;
            dataTConfig = d.Data.QL.Data.TRIAL_CONFIG;
            dataTConfigh= d.Data.QL.Headers.TRIAL_CONFIG;
            dataMsB = d.Data.QL.Data.BURT_STATUS;
            clear d
            
            % defines num
            MID_FT  = 67;   % NETBOX
            MID_WAM = 62;   % ROBOT
            mid_type = [MID_FT, MID_WAM];
            
            %%%%%%%%%%%%%%%%% DATA READING PART %%%%%%%%%%%%%%%%%%%%
            %%%%% 1. read times from blackrock and check it value (ifplot)
            %%%%%%%%%%%%%%%%%%%%
            eventsL = dataTs.eventsL;   % events_label, 
            eventsT = dataTs.eventsT;   % events_time;
            if (isfield(dataTs, 'eventsTrials')) % trial list
                eventTrials= dataTs.eventsTrials;
            end
            
            % pulse from the blackrock
            events_type = unique(eventsL);
            for etype = events_type
                bk_time{etype} = eventsT(eventsL == etype);
                bk_trials{etype} = eventTrials(eventsL == etype); % after ss3090
            end
            
            %   etype:  1-FT, 
            %   etype:  2-ROBOT
            ifplot = 0;
            if (ifplot)
                clf;
                axh(1) = subplot(2,1,1); hold on;
                plot(bk_time{1}, 'r*'); plot(bk_time{2}, 'b*');
                legend('FT', 'WAM')
                axh(2) = subplot(2,1,2);
                stem(bk_time{1} - bk_time{2});
                linkaxes(axh, 'x');
                title(axh(1),'sync signal times to BlackRock');
                title(axh(2), 'syn signal difference to Blackrock');
                
                clf;
                axh(1) = subplot(2,1,1); hold on;
                plot(diff(bk_time{1}), 'r*');
                axh(2) = subplot(2,1,2);
                plot(diff(bk_time{2}), 'b*');
                linkaxes(axh, 'x');
                title(axh(1),'time diff between ft on BlackRock');
                title(axh(2), 'time pulse between wam on Blackrock');
            end
            
            %%%%% 2. read times message and check it value (in ifplot)
            % message time for each computer
            time_leading = dataMsT.tleading;
            time_lasting = dataMsT.tlasting;
            % correspond the pulse time to each trial number
            if (isfield(dataMsT, 'trial_no'))
                times_trialno= dataMsT.trial_no; % after ss3110
            else 
                % need to use the other information to decide trials 
                msg_recvt = dataMsTh.recv_time;
                trial_no = dataTConfig.trial_no;
                trial_no_t = dataTConfigh.recv_time;
                trial_no_t_extrap = interp1(1:length(trial_no_t), trial_no_t, 1:(length(trial_no_t)+1), 'linear', 'extrap');
                % figure(); hold on;
                % plot(msg_recvt, 'r*');
                % plot(trial_no_t, 'b*');
                
                [n, tidx ] = histc(msg_recvt, trial_no_t_extrap);
                times_trialno = trial_no(tidx);
            end
            msg_mid = dataMsTh.src_mod_id;
            %mid_type = unique(msg_mid);
            mid_type = [MID_FT, MID_WAM];
            for mtype_idx = 1:length(mid_type)
                tleading{mtype_idx} = time_leading(msg_mid == mid_type(mtype_idx));
                tlasting{mtype_idx} = time_lasting(msg_mid == mid_type(mtype_idx));
                t_means{mtype_idx}  = (tleading{mtype_idx}+tlasting{mtype_idx})/2;
                t_error{mtype_idx}  =-(tleading{mtype_idx}-tlasting{mtype_idx});
                if (exist('times_trialno', 'var'))
                    t_msg_trialidx{mtype_idx} = times_trialno(msg_mid == mid_type(mtype_idx)); 
                end
            end
            
            ifplot = 0;
            if (ifplot)
                clf; % pfem time precision
                axh(1) = subplot(2,1,1); hold on;
                plot(tleading{1}, 'r*'); plot(tlasting{1}, 'b*');
                legend('tleading', 'tlasting')
                axh(2) = subplot(2,1,2);
                stem(1e6*(tleading{1} - tlasting{1}));
                ylim([-1 1]*1e1);
                linkaxes(axh, 'x');
                title(axh(1), 'tleading and tlasting (s)');
                title(axh(2), 'time elapse in pfem (us)');
                
                clf; % ernie time precision
                axh(1) = subplot(2,1,1); hold on;
                tleading{2} - tleading{2}(1)
                plot(tleading{2}, 'r*'); plot(tlasting{2}, 'b*');
                legend('tleading', 'tlasting')
                axh(2) = subplot(2,1,2);
                stem(1e6*(tleading{2} - tlasting{2}));
                ylim([-1 1]*1e1);
                linkaxes(axh, 'x');
                title(axh(1), 'tleading and tlasting (s)');
                title(axh(2), 'time elapse in ernie (us)');
                
                clf;
                axh(1) = subplot(2,1,1);
                plot(diff(tlasting{1}), '*');
                ylabel('sec');
                title('intertrial t\_lasting duration, FT');
                axh(2) = subplot(2,1,2);
                plot(diff(tlasting{2}), '*');
                title('intertrial t\_lasting duration, WAM');
                ylabel('sec');
                linkaxes(axh, 'x');
                
            end
            
            % select the time will be align with the blackrock time, here
            % use the average between 't_leading' and 't_lasting'
            t_interest = t_means;
            
            % - Think this condition as the last one did not record in BK.
            if (exist('t_msg_trialidx', 'var') && exist('bk_trials', 'var'))
                % 1. find the intersect of trials
                [trial_its, idx_msg1, idx_bk1] = intersect(t_msg_trialidx{1}, bk_trials{1}); % FT
                [trial_its, idx_msg2, idx_bk2] = intersect(t_msg_trialidx{2}, bk_trials{2}); % FT
                % assuem every FT sync signal has a WAM sync signal
                
                % 2. change to-aligned data into certain trials
                t_interest{1} = t_interest{1}(idx_msg1);
                t_interest{2} = t_interest{2}(idx_msg2);
                bk_time{1} = bk_time{1}(idx_bk1);
                bk_time{2} = bk_time{2}(idx_bk2);
                
                
            else % old way to deal with the two message do not have the same length problem
                if (length(t_interest{1}) == 1 + length(bk_time{1}))
                    disp('ERROR(MAYBE): Message and pulse size inconsistance, in SessionScan::updateTimeBlackRock()');
                    
                    t_interest{1} = t_interest{1}(1:end-1);
                    t_interest{2} = t_interest{2}(1:end-1);
                end
             end
            
            ifplot = 1;
            if (ifplot)
                clf;
                hold on;
                plot(bk_time{1}, t_interest{1}, 'r*'); 
                plot(bk_time{2}, t_interest{2}, 'b*');
                legend('FT', 'WAM');
                xlabel('BK time');
                ylabel('each computer time');
            end
            % intropolate each data time to the bk_time;
            
            %%%%%%%%%%%%%% TIME EXPORT PART %%%%%%%%%%%%%%%%
            %%% 1. the FT time
            try
            obj.force_t = interp1(t_interest{1}, bk_time{1}, obj.ft.elapse, 'linear', 'extrap'); 
            catch 
                if (ismember(obj.ssnum, [3197 3198]))
                    obj.force_t = interp1(t_interest{2}, bk_time{2}, obj.ft.elapse, 'linear', 'extrap'); 
                    % this trial do not have the FT pulse, use WAM pulse to
                    % debug...
                end
            end
            
            ifplot = 1;
            if (ifplot)
                clf; 
                hold on;
                plot(t_interest{1}, bk_time{1}, 'b*');
                plot(obj.ft.elapse, obj.force_t, 'r.');
            end
            
            %%% 2. the wam time
            % 1. plot the wam time with the ernie computer time
            wam_time = dataMsB.wamt;
            ernie_time=dataMsB.erniet;
            
            if (length(wam_time) > length(unique(wam_time)))
                [C, IA, ~] = unique(wam_time);
                wam_time = wam_time(IA);
                ernie_time = ernie_time(IA);
            end
            clf; 
            ernie_timeh = interp1(wam_time, ernie_time, obj.wam.time, 'linear', 'extrap');
            obj.wam_t  = interp1(t_interest{2}, bk_time{2}, ernie_timeh, 'linear', 'extrap'); 
            if (size(obj.wam_t,1) > size(obj.wam_t,2)) % a thin matrix
                obj.wam_t = obj.wam_t';
            end
            ifplot = 1;
            if(ifplot)
                clf;
                axh(1) = subplot(2,2,1);
                hold on;
                plot(ernie_time, wam_time, 'b*');
                plot(ernie_timeh, obj.wam.time, 'r.');
                xlabel('ernie time'); ylabel('wam time');
                axh(2) = subplot(2,2,3);
                hold on;
                plot(t_interest{2}, bk_time{2}, 'b*');
                plot(ernie_timeh, obj.wam_t, 'r.');
                xlabel('ernie time'); ylabel('BK time');
                axh(3) = subplot(2,2,4); 
                plot(obj.wam.time, obj.wam_t, 'r.'); 
                xlabel('wam time'); ylabel('BK time');
                linkaxes([axh(1), axh(2)], 'x');
                linkaxes([axh(2), axh(3)], 'y');
                
            end
            flag = 1;
            return;
        end
    end

end

function [steadyValue, duration] = findSteadyValue(intMat, durat, ifplot)
    % steadyValue = findSteadyValue(interestMat)
    % for a interestMatrix (intMat) with m-by-n values (which m serves as
    % series and n is time series), find if the n is steady in the last
    % some duration. 
    % Will find the duration (data points) and report it if not specify a
    % duration,
    % minimum duration: 50, for 0.1s when the sampling rate is 500Hz
    if ~exist('durat', 'var')
        durat = -1;
    end
    if ~exist('ifplot', 'var')
        ifplot = 0;
    end
    min_duration = 50;
    [r,c] = size(intMat);
    if c<min_duration
        disp(['ERROR in findSteadyValue: not enough length of data']);
        return;
    end
    if r == 0
        disp('ERROR in findSteadyValue: not enough data points');
        return;
    end
    % define a threshold th, smaller than which will be regarded as steady
    % th = k*std(dat_value); k=0.05;
% % %     k = 1e-10;
% % %     for try_i = 1:100
% % %         ths = zeros(r,1);                   % thresholds for selection
% % %         row_idx = ones(size(intMat,1), 1);  % index, that which data have big undulate
% % %         thresholds = sort(range(intMat,2), 'ascend');
% % %         threshold  = ((thresholds(end)-thresholds(1))*1/3+thresholds(1)) * k;
% % %         for r_i = 1:r
% % %             row_idx(r_i) = range(intMat(r_i,:))>threshold;
% % %         end
% % %         for r_i = 1:r
% % %             if (row_idx(r_i))
% % %                 ths(r_i) = std(intMat(r_i,:), 'omitnan');
% % %             else
% % %                 ths(r_i) = inf;
% % %             end
% % %         end
% % %         % sort and choose the second leatest to avoid when th==0
% % %         th_tmp = sort(unique(ths), 'ascend');
% % %         try
% % %             th = th_tmp(2);
% % %         catch % smaller than 2
% % %             th = th_tmp(1);
% % %         end
% % %         % choose index on the ones have significant change.
% % %         
% % %         % find the steady state index from end of the matrix
% % %         idx = prod([zeros(r,1) diff(intMat, 1, 2)] < th, 1);
% % %         
% % %         %    k = k+0.01;
% % %         %end
% % %         % find the largest consequtive
% % %         i = find(diff(idx));
% % %         n = [i numel(idx)] - [0 i];
% % %         c = arrayfun(@(X) X-1:-1:0, n , 'un',0);
% % %         y = cat(2,c{:});
% % %         if(sum(idx)) >= 1200
% % %             k = k*2;
% % %         end
% % %         if(sum(idx)) > 100 && (sum(idx)) <1200
% % %             break
% % %         end
% % %         %         k = k*2;
% % %     end
% % %     [~, idx_stt] = max(y.*idx);
% % %     idx_edn = idx_stt + y(idx_stt);
% % %     % change idx_stt according to durat specify
% % %     if (durat == -1) % unspecified length, use all steady state
% % %         duration = idx_edn - idx_stt + 1;
% % %     else
% % %         idx_edn = idx_edn - 10; % random offset, make sure right value
% % %         idx_stt = idx_edn - durat + 1;
% % %         duration = durat;
% % %     end
    % get data by summation
    data = sum(abs(intMat'), 2);
    datad = diff(data');
    % detect edge 
    edge_arr = [datad(1) datad];
    edge_abs = abs(edge_arr);
    % 1. find boundaries of edge (two maximum change) 
    [edge_sort, maxidx] = sort(edge_abs, 'descend');
    cond = true; 
    i = 2;
    while(cond)
        edge1idx = maxidx(1);
        edge2idx = maxidx(i);
        if (abs(edge1idx-edge2idx) > c/2) % two edges are faraway
            cond = false;
        else 
            i = i+1; % refind edge2
        end
    end
    % 2. find data steady
    [n, x] = hist(edge_arr); 
    x_boundDist = range(x)/length(x);
    diff_tolerance = x(n == max(n)) + [-1 1]*x_boundDist;
    arr_idx1 = [edge_arr > diff_tolerance(1)] &...
        [edge_arr<diff_tolerance(2)] & ...
        [1:c] > min([edge1idx edge2idx]) &...
        [1:c] < max([edge1idx edge2idx]);
    
    % find the largest consequtive
    i = [find(diff(arr_idx1))];
    n = [i numel(arr_idx1)] - [0 i];
    cels = arrayfun(@(X) X-1:-1:0, n , 'un',0);
    y = cat(2,cels{:});

    [~, idx_stt] = max(y.*arr_idx1);
    idx_edn = idx_stt + y(idx_stt);
    % change idx_stt according to durat specify
    if (durat == -1) % unspecified length, use all steady state
        duration = idx_edn - idx_stt + 1;
    else
        idx_edn = idx_edn - 10; % random offset, make sure right value
        idx_stt = idx_edn - durat + 1;
        duration = durat;
    end

    if idx_edn - idx_stt < 50
        disp('CONDITION: not enough long steady state, use default 50');
        duration = 50;
        idx_stt = idx_edn - duration + 1;
    end
    steadyValue = mean(intMat(:, idx_edn-duration+1:idx_edn),2);
    % HOW TO FIND THE LAST CONSECUTIVE 1s? 
    if (ifplot)
        figure();
        hold on; 
        plot(intMat'); 
        line([idx_stt, idx_stt], [min(intMat(:)), max(intMat(:))]);
        line([idx_edn, idx_edn], [min(intMat(:)), max(intMat(:))]);
        title('data selection demo');
    end
    return
end
