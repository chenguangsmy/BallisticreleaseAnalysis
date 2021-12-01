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
        force_t         % time for ft seperate data
        wam_t
        emg_h
        emg_t
        data            % all the data including force and wam information (aligned)
        %%% other variables
        endpoint0 = [-0.517 0.481 0.0] 
        col_vec = colormap('lines');
        badTrials = [1];       % bad trial, cull in data
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

            flag_progress = 0;      % progress display

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
                if (flag_progress)
                    if ~(trial_percent == percent_prev) % avoid showing a lot.
                        fprintf('  %02d%%...', trial_percent);
                        percent_prev = trial_percent;
                    end
                end
                obj.trials(trial_i) = TrialScan(obj, trial_i);
                % align to mov
                obj.trials(trial_i) = alignMOV(obj.trials(trial_i));
                %obj.trials(trial_i) = alignPertInit(obj.trials(trial_i));
                if (flag_progress)
                    if (trial_i == length(trials_all)) % last trial
                        fprintf('  100%%  FINISHED!\n');
                    end
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
        function displayBlockCondition(obj)
            % displayBlockCondition(obj)
            all_trials = length(obj.trials);
            fin_trials = sum([obj.trials.outcome]==1);
            rate = fin_trials/all_trials;
            if (length(obj.tarLs) == 1)
                fprintf("tar: %.1f(cm), F: %d(N): %d/%d, rate: %f \n" ,...
                    obj.tarLs(1)*100, obj.fThs(1), fin_trials, all_trials, rate);
            else 
                fprintf("Stoc, F: %d(N): %d/%d, rate: %f \n" ,...
                 obj.fThs(1), fin_trials, all_trials, rate);
            end
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
            return
        end

        %%% other process
        function obj_new = ConcatTrials(obj1, obj2, trial_idx1, trial_idx2)
            trials = [obj1.trials(trial_idx1) obj2.trials(trial_idx2)];
            obj_new = obj1;
            obj_new.trials = trials;
        end
        function [obj] = add_to_PertData(obj, new_Data_pert)
            n_test1 = length(obj.wam.Data_pert);
            n_test2 = length(new_Data_pert);

            for i = 1:n_test2
                obj.wam.Data_pert(i+n_test1) = new_Data_pert(i);
            end
            
        end
        

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
            %   Fp: 1-by-N matrix, perturbation force, in stoc 2-by-N(xy)
            %   ts: 1-by-N matrix, task states
            %   t: 1-by-N matrix, time 
            %   mvst: the mask that robot can freely move
            trial_perturbs = [obj.trials.ifpert];
            pert_max = max(3, (max(trial_perturbs)+1)); % 0,nopert; 1, pulse; 2, stoc
            t_idx = cell(2,pert_max);                   % succ/failure * pert_types
            for sf = 1:2
                for p_i = 1:pert_max
                     t_idx{sf,p_i} = intersect(find(trial_perturbs==p_i-1), ...
                                    find([obj.trials.outcome]==2-sf));  % 1 or 0
                     t_idx{sf,p_i} = setdiff(t_idx{sf,p_i}, 1);
                end
            end
            
            % TODO: consider remove correction trials here? 
            
            % export part ...
            
            if ~isempty([t_idx{1:2,3}]) % stoc-perturbed trials. 
                % Stocpert do not in the same session with step ones in current experiments
                cellsmat = cell(2, length(obj.tarLs), 15, 3); % think!!! 
                for sf = 1:2 % suceed and failed trials
                    for tl_i = 1:length(obj.tarLs)
                        %trial_list = setdiff(find([obj.trials.tarL] == obj.tarLs(tl_i)),1);
                        trial_list = intersect(find([obj.trials.tarL] == obj.tarLs(tl_i)) , t_idx{sf,3});
%                         trial_list = setdiff(trial_list,1);
                        for t_i = 1:length(trial_list)
                            t_tmp = obj.trials(trial_list(t_i));
                            cellsmat{sf,tl_i,t_i,3} = t_tmp.export_as_formatted;
                        end
                    end
                end
            else
                trial_max = max([length(t_idx{1,1}), length(t_idx{1,2}),...
                                length(t_idx{2,1}), length(t_idx{2,2})]);
                cellsmat = cell(2,length(obj.tarLs),trial_max,3); % for block design, only fill idx_tarLs==1 
                for sf =  1:2
                    for p_i = 1:2
                        %trial_list = setdiff(t_idx{p_i},1);
%                         trial_list = find([obj.trials.outcome] == 1); % no failed trials
                        trial_list = t_idx{sf,p_i};
                        %trial_list = find([obj.trials.outcome] ~= -1); % spring test parameter selection
%                         trial_list = intersect(trial_list, t_idx{p_i});
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
%                         trial_list = setdiff(trial_list,[1]);
                        for t_i = 1:length(trial_list)
                            t_tmp = obj.trials(trial_list(t_i));
                            cellsmat{sf,1,t_i,p_i} = t_tmp.export_as_formatted;  % each trial
                            xlim([-5 -4])
                            ifplot = true;
                            if (ifplot)
                                subplot(2,1,1);
                                plot(cellsmat{sf,1,t_i, p_i}.t, cellsmat{sf,1,t_i, p_i}.x(2,:));
                                subplot(2,1,2);
                                plot(cellsmat{sf,1,t_i, p_i}.t, cellsmat{sf,1,t_i, p_i}.f(2,:));
                            end
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
                            if(isempty(cellsmat{1,t_i, p_i}))
                                continue;
                            end
                            dt = diff(cellsmat{1,t_i,p_i}.t);
                            dt = [dt(1) dt];
                            plot(t_i*plt_offset + dt);
                        end
                    else
                        if isempty(t_idx{p_i}) 
                            continue
                        end
                        tiofst = 0; % plot offset
                        for tl_i = 1:size(cellsmat, 2)
                            for t_i = 1:length(cellsmat(1,tl_i,:,p_i))
                                if ~isempty(cellsmat{1,tl_i,t_i,p_i})
                                dt = diff(cellsmat{1,tl_i,t_i,p_i}.t);
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
