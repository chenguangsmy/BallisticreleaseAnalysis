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
        tarFs
        %%% task variables 
        hand_pos        % position read from WAM endpoint
        hand_pos_offset % the center_pos for WAM endpoint
        force           % force in the force transducer
%         FTrot_M = ...   % global: x-right, y-front, z-up, FT_base x-backup, y-frontup, z-left
%            [0          0           -cosd(180)
%             +sind(45)  -cosd(45)    0
%             -cosd(45)   -sind(45)    0];
        FTrot_M = ...
            [cosd(45) -sind(45)    0
             sind(45)  cosd(45)    0
                    0          0     1];
        taskState
        trials TrialScan% member function
        %%% perturbation variables
%       pert_state = 3  % only perturb at force ramp pert_state == 3
%         pert_time       % perturbation start and end time
%         pert_rdt        % perturbation start read time, (same with wam)
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
        emg_evl
        opt_t           % optotrak time
        data            % all the data including force and wam information (aligned)
        %%% other variables
        endpoint0 = [-0.517 0.481 0.0] 
        col_vec = colormap('lines');
        badTrials = [1];       % bad trial, cull in data
        stateNames = {'Begin', 'Present', 'FrcRamp', 'FrcHold', 'Move', 'Hold', 'End', 'Reset'};
        Fs = 2000; % sample frequency
    end
    
    methods
        %%% process
        function obj = SessionScan(ss_num, if_reload) %(inputArg1,inputArg2)
            %VARSCAN Construct an instance of this class
            %   Detailed explanation goes here
            % load the INTERMEDIATE file for the time 
            % if_reload=0, set to 0 if I need to reload the data
            if (~exist('if_reload', 'var'))
                if_reload = 0;
            end
            obj.ssnum = ss_num;
            file_name = ['KingKong.0' num2str(obj.ssnum) '.mat']; % the previous 'formatted' version
            file_dir_formmed = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data';
            file_dir_int = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/Intermediate';
            filename_fmt = ['KingKong.0' num2str(obj.ssnum) '.format.mat']; % my formated version
            if exist(fullfile(file_dir_formmed, filename_fmt), 'file') && ~if_reload% ifloadagain
                obj = obj.loadData();
            else
            data1 = load([file_dir_int  '/' file_name], 'Data');
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
                ifplot = 0;
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
            % cg: we use EMG and OPT as a reference in the movement Aug10th, f2022
            ifemg = 1;
            if(ifemg) % the key to load EMG data 
            try 
                if (flag_progress)
                    disp('Loading (intermed/raw) EMG...');
                end
                obj.emg = SessionScanEMG(obj.ssnum, 0); % normal blocks
%                 obj.emg = SessionScanEMG(obj.ssnum, 1); % reload EMG
                obj.emg_t = obj.emg.data.t;
                obj.emg_h = obj.emg.data.emg;
                obj.emg_evl=obj.emg.data.emgevl; % envolope
            catch

%                 disp('no EMG data here! ')
            end
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
            obj.tarFs  = unique(obj.Data.TaskJudging.Target(4,obj.Data.TaskStateMasks.Move)); % deviate 0
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
            % syncflag: 
            %           -1: no TSync file 
            %           -2: no intermediate file 
            %           -3: no FT pulse was detected 
            if syncflag == -3
                % under this case, map the FT time into the BK time 
                
                % rely on the message to get the data aligned
                obj = obj.forceHighSample(obj.ft);%, ft_intm);
            end
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

            obj = dealingOPTMissingShoulderMarker(obj); 
            ifdotrialfy = 1;
            if (ifdotrialfy)
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
%                 obj.trials(trial_i) = judgeOffline(obj.trials(trial_i));
                %obj.trials(trial_i) = alignPertInit(obj.trials(trial_i));
                
                if (flag_progress)
                    if (trial_i == length(trials_all)) % last trial
                        fprintf('  100%%  FINISHED!\n');
                    end
                end
            end
            
            % deal with the trial sucess/failure
            [s,f] = readManualSetsf(obj.ssnum);
            
            if (~isempty(s)) %specify good trials
                disp('Manual set trials suceed: ');
                for trial_i = s
                    fprintf('%d ', trial_i);
                    obj.trials(trial_i).outcome = 1;
                end
            end
            if (~isempty(f)) %specify bad trials
                disp('Manual set trials failed: ');
                for trial_i = f
                    fprintf('%d ', trial_i);
                    obj.trials(trial_i).outcome = 0;
                end
            end
            
            
            obj.trials(1).outcome = 0; % 1st trials are bad for force align.
            
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
            obj = obj.dealingEMGOpenLoopNoise();   % bad EMG to 0

            end % from 'ifloadagain'
            obj.saveFile(); 
            if_calibSS = findCalibSS(obj.ssnum);
            if (~if_calibSS)
                if (size(obj.data.ox,2)~=0)
                   obj = obj.dealingOPTMarkersErr();      % bad OPT markers
                end
            end
            end
        end
        function obj = dealingSessionsExceptions(obj)
            % solve some data-code inconsistant problem, specify for each
            % sessions
            
            % deal with force exceptions
            force_exception_sessions = [3683:3691, ...
                                        3722:3725, ...
                                        3727:3728, ...
                                        3737, 3740, ...
                                        3766, 3767, 3768:3778, ...
                                        3793:3795, ...
                                        3803:3812, 3856:3860, ...
                                        3873:3884, ...
                                        3906:3912, 3925:3937, ...
                                        3987:4002, 4046:4053];
            % in these sessions, I wrongly calibrate the force, that the
            % collected force is biased for certain value. To deal with
            % this exception, the only way is to add the force value of ts7
            % to avoid the force calibration problem, do the force
            % exception for every single trial
            if (sum(obj.ssnum == force_exception_sessions) ~= 0 )
                for trial_i = 1:obj.trials_num
                    obj.trials(trial_i) = obj.trials(trial_i).dealForceException();
                end
            end

        end
        function obj = dealingOPTMissingShoulderMarker(obj)
            % when the case the shoulder marker is not shown 
            % shoulder marker index-3
            % another shoulder marker index-4
            % Assume a constant displacemnt, convert another marker's
            % position to the current marker

            % 1. Check the displacement between markers 
            ifplot = 1; 
            if (ifplot)
              clf;
              figure();
              axh(1) = subplot(4,1,1); hold on; 
              plot(obj.data.t, obj.data.ox(1,:,3), 'r.');
              plot(obj.data.t, obj.data.ox(1,:,4), 'g.');
              legend('marker3', 'marker4');
              title('position-x (m)');
              ylabel('x (m)');
              xlabel('time (s)');
              axh(2) = subplot(4,1,2); hold on; 
              plot(obj.data.t, obj.data.ox(2,:,3), 'r.');
              plot(obj.data.t, obj.data.ox(2,:,4), 'g.');
%                 plot(obj.data.t, obj.data.x(1,:));
%                 plot(obj.data.t, obj.data.tNo);
              title('position-y (m)');
              ylabel('y (m)');
              xlabel('time (s)');
              axh(3) = subplot(4,1,3); hold on;
              plot(obj.data.t, obj.data.ox(3,:,3), 'r.');
              plot(obj.data.t, obj.data.ox(3,:,4), 'g.');
              title('position-z (m)');
              ylabel('z (m)');
              xlabel('time (s)');
              axh(4) = subplot(4,1,4); hold on;
              plot(obj.data.t, obj.data.ox(1,:,3) - obj.data.ox(1,:,4), 'c.');
              plot(obj.data.t, obj.data.ox(2,:,3) - obj.data.ox(2,:,4), 'm.');
              plot(obj.data.t, obj.data.ox(3,:,3) - obj.data.ox(3,:,4), 'y.');
              legend('dx', 'dy', 'dz');
              title('difference');
              ylabel('\Delta x (m)');
              xlabel('time (s)');
              linkaxes(axh, 'x');
            end
            
             % 2. Add the shoulder marker information by another marker

            delta_pos = (obj.data.ox(:,:,3) - obj.data.ox(:,:,4)); % stop here and debug...!!!
%             delta_pos = reshape(delta_pos, size(delta_pos,[2,3]));
            delta_pos_avg = mean(delta_pos, 2, 'omitnan');

            idx_n3v4 = isnan(obj.data.ox(1,:,3)) & ~isnan(obj.data.ox(1,:,4));
            idx_n4v3 = isnan(obj.data.ox(1,:,4)) & ~isnan(obj.data.ox(1,:,3));
            
            % length 
            if (sum(idx_n3v4)>0)
                Delta_pos_avg1 = repmat(delta_pos_avg,1,[sum(idx_n3v4)]);
                Delta_pos_avg2 = zeros(3,sum(idx_n3v4),1);
                Delta_pos_avg2(:,:,1) = Delta_pos_avg1;
                obj.data.ox(:,idx_n3v4,3) = obj.data.ox(:,idx_n3v4,4) + Delta_pos_avg2;
            end
            if (sum(idx_n4v3)>0)
                Delta_pos_avg1 = repmat(delta_pos_avg,1,[sum(idx_n4v3)]);
                Delta_pos_avg2 = zeros(3,sum(idx_n4v3),1);
                Delta_pos_avg2(:,:,1) = Delta_pos_avg1;
                obj.data.ox(:,idx_n4v3,4) = obj.data.ox(:,idx_n4v3,3) - Delta_pos_avg2;
            end

            % see it again
            if (ifplot)
%               clf;
              figure();
              axh(1) = subplot(4,1,1); hold on; 
              plot(obj.data.t, obj.data.ox(1,:,3), 'r.');
              plot(obj.data.t, obj.data.ox(1,:,4), 'g.');
              legend('marker3', 'marker4');
              title('position-x (m)');
              ylabel('x (m)');
              xlabel('time (s)');
              axh(2) = subplot(4,1,2); hold on; 
              plot(obj.data.t, obj.data.ox(2,:,3), 'r.');
              plot(obj.data.t, obj.data.ox(2,:,4), 'g.');
              title('position-y (m)');
              ylabel('y (m)');
              xlabel('time (s)');
              axh(3) = subplot(4,1,3); hold on;
              plot(obj.data.t, obj.data.ox(3,:,3), 'r.');
              plot(obj.data.t, obj.data.ox(3,:,4), 'g.');
              title('position-z (m)');
              ylabel('z (m)');
              xlabel('time (s)');
              axh(4) = subplot(4,1,4); hold on;
              plot(obj.data.t, obj.data.ox(1,:,3) - obj.data.ox(1,:,4), 'c.');
              plot(obj.data.t, obj.data.ox(2,:,3) - obj.data.ox(2,:,4), 'm.');
              plot(obj.data.t, obj.data.ox(3,:,3) - obj.data.ox(3,:,4), 'y.');
              legend('dx', 'dy', 'dz');
              title('difference');
              ylabel('\Delta x (m)');
              xlabel('time (s)');
              linkaxes(axh, 'x');
            end
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
% %             tarFs = obj.tarFs;
% %             tarFs_num = length(tarFs);
% %             sT = zeros(tar_num, tarl_num, tarFs_num);    % sucessful trials
% %             tT = zeros(tar_num, tarl_num, tarFs_num);    % total trials
% %             sR = zeros(tar_num, tarl_num, tarFs_num); 
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
% %                     for fThi = 1:tarFs_num
% %                         tT(tard_i, tarl_i, fThi) = ...
% %                             sum(tarR == tard(tard_i) &...
% %                             tarL == tarl(tarl_i) & ...
% %                             fTh  == tarFs(fThi));
% %                         sT(tard_i, tarl_i, fThi) = ...
% %                             sum(tarR == tard(tard_i) &...
% %                             tarL == tarl(tarl_i) & ...
% %                             fTh  == tarFs(fThi) & ...
% %                             [obj.trials.outcome] == 1);
% %                         sR(tard_i, tarl_i, fThi) = sT(tard_i, tarl_i, fThi)/tT(tard_i, tarl_i, fThi);
% %                     end
% %                 end
% %             end
% %             sR_2d = reshape(sR(1,:,:), size(sR, 2), size(sR, 3));
% %             sR_table = [[obj.tarFs]', sR_2d'];
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
                    obj.tarLs(1)*100, obj.tarFs(1), fin_trials, all_trials, rate);
            else 
                fprintf("Stoc, F: %d(N): %d/%d, rate: %f \n" ,...
                 obj.tarFs(1), fin_trials, all_trials, rate);
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
% %             tarFs = obj.tarFs;
% %             tarFs_num = length(tarFs);
% %             trialTime = zeros(tar_num, tarl_num, tarFs_num);
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
% %                     for fThi = 1:tarFs_num
% %                         trialid = ...
% %                             (tarR == tard(tard_i) &...
% %                             tarL == tarl(tarl_i) & ...
% %                             fTh  == tarFs(fThi));
% %                         time_all = [obj.trials(trialid).edn_t] - [obj.trials(trialid).bgn_t];
% %                         time_mean(tard_i, tarl_i, fThi) = mean(time_all);
% %                     end
% %                 end
% %             end
% %             time_2d = reshape(time_mean(1,:,:), size(time_mean, 2), size(time_mean, 3));
% %             time_table = [[obj.tarFs]', time_2d'];
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
            % resample all the data in a same frequency: 2000Hz
            % 
            t_min = min(obj.wam_t); 
            t_max = max(obj.wam_t); 
            Fs = 2000; 
            pts = floor((t_max - t_min))/(1/Fs);
            t_max_ = t_min + (pts-1)*(1/Fs);
            sample_t = linspace(t_min, t_max_, pts);
            force_h = interp1(obj.force_t', obj.ft.force', sample_t', 'linear', nan)'; 
            torque_h = interp1(obj.force_t', obj.ft.torque_origin', sample_t', 'linear', nan)';
            
            %%% deal with the task state masks
            ts_masks = obj.Data.TaskStateMasks;
            ts_masks_mat = zeros(8,length(obj.Data.TaskStateMasks.Begin));
            ts_masks_list={'Begin','Present','FrcRamp','FrcHold','Move','Hold','End','Reset'};
            for mask_i = 1:8
                ts_masks_mat(mask_i,:) = mask_i * eval(['obj.Data.TaskStateMasks.' ts_masks_list{mask_i}]);
            end
            ts_ = sum(ts_masks_mat);
            idx_nan = isnan(obj.Data.SpikeTimestamp);
            ts_h = interp1(obj.Data.SpikeTimestamp(~idx_nan),...
                ts_(~idx_nan), sample_t, 'next');
            %%% deal with the messaged position and velocity
            x_msg_h = interp1(obj.Data.SpikeTimestamp(~idx_nan)', ...
                obj.Data.Position.Actual(~idx_nan,:), sample_t', 'previous');
            x_msg_h = x_msg_h';
            v_msg_h = interp1(obj.Data.SpikeTimestamp(~idx_nan)', ...
                obj.Data.Velocity.Actual(:,~idx_nan)', sample_t', 'previous');
            v_msg_h = v_msg_h';
            f_rotate = obj.rotateAxisForce();
            f_msg_h = interp1(obj.Data.SpikeTimestamp(~idx_nan)', ...
                f_rotate(1:3,~idx_nan)', sample_t', 'previous');
            f_msg_h = f_msg_h';

            ifplot = 0;
            if (ifplot)
                clf;
                axh(1) = subplot(3,1,1);  hold on;
                plot(obj.force_t, obj.ft.force(1,:), 'r.');
                plot(sample_t, force_h(1,:), 'b.');
                axh(2) = subplot(3,1,2);  hold on;
                plot(obj.force_t, obj.ft.force(2,:), 'r.');
                plot(sample_t, force_h(2,:), 'b.');
                axh(3) = subplot(3,1,3);  hold on;
                plot(obj.force_t, obj.ft.force(3,:), 'r.');
                plot(sample_t, force_h(3,:), 'b.');
                linkaxes(axh, 'x');
            end
            
            ifplot = 1;     % torque
            if (ifplot)
                clf;
                axh(1) = subplot(3,1,1);  hold on;
                plot(obj.force_t, obj.ft.torque_origin(1,:), 'r.');
                plot(sample_t, torque_h(1,:), 'b.');
                axh(2) = subplot(3,1,2);  hold on;
                plot(obj.force_t, obj.ft.torque_origin(2,:), 'r.');
                plot(sample_t, torque_h(2,:), 'b.');
                axh(3) = subplot(3,1,3);  hold on;
                plot(obj.force_t, obj.ft.torque_origin(3,:), 'r.');
                plot(sample_t, torque_h(3,:), 'b.');
                linkaxes(axh, 'x');
            end

            obj.data.t = sample_t;
            obj.data.x = interp1(obj.wam_t', obj.wam.tp', sample_t')';
            obj.data.v = interp1(obj.wam_t', obj.wam.tv', sample_t')';
            obj.data.Fp = interp1(obj.wam_t', obj.wam.cf', sample_t')';
            obj.data.tq = interp1(obj.wam_t', obj.wam.jt', sample_t')';
            obj.data.jp = interp1(obj.wam_t', obj.wam.jp', sample_t')';
            obj.data.ts = interp1(obj.wam_t', obj.wam.state', sample_t', 'previous')';
            obj.data.tsf = ts_h; % get from the formatted data
            obj.data.x_msg = x_msg_h;
            obj.data.v_msg = v_msg_h;
            obj.data.f = force_h;
            obj.data.ftq= torque_h; % force transducer cencored torque
            obj.data.f_msg = f_msg_h;
            idx_spiket_nnan = ~isnan(obj.Data.SpikeTimestamp);
            obj.data.tNo = int32(interp1(obj.Data.SpikeTimestamp(idx_spiket_nnan)', double(obj.Data.TrialNo(idx_spiket_nnan))', sample_t', 'previous')');
            obj.data.sNo = int32(ones(size(sample_t))*obj.ssnum);
            
            % deal with exception of session 4339
            ifplot = 1; % if check
            if unique(obj.data.sNo) == 4339 % I missed saving the wam.bin
                % detect when the release happens 
                fce_thr = 5*std(diff(obj.data.f'));
                obj.data.ts = obj.data.tsf;
                fcediff = [0 diff(obj.data.f(1,:))];
                for trial_i = unique(obj.data.tNo)
                    trial_idx = obj.data.tNo == trial_i;
                    fcechgidx = find(fcediff< -fce_thr(1) & obj.data.tsf == 5 & obj.data.tNo == trial_i);
                    statechgidx= find(obj.data.tsf == 5 & obj.data.tNo == trial_i);
                    if (~isempty(statechgidx))
                    obj.data.ts(statechgidx(1):fcechgidx(1)) = 4;
                    if (ifplot)
                        clf;
                        axh(1) = subplot(3,1,1); 
                        plot(obj.data.t(trial_idx), obj.data.ts(trial_idx)); 
                        axh(2) = subplot(3,1,2); 
                        plot(obj.data.t(trial_idx), obj.data.f(1,trial_idx));
                        axh(3) = subplot(3,1,3); 
                        plot(obj.data.t(trial_idx), fcediff(trial_idx));
                        linkaxes(axh, 'x');
                    end
                    end
                end

            end

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
            
            if (~isempty(obj.opt))
                opt_t1 = interp1(1:length(obj.opt_t(~isnan(obj.opt_t))), obj.opt_t(~isnan(obj.opt_t)), 1:length(obj.opt_t), 'linear', 'extrap'); % not guareentee!  
                obj.opt_t(1) = opt_t1(1); %??? possible???
%                 data.optx = interp1(obj.opt_t, obj.opt.datah.x', obj.wam_t', 'linear', 'extrap')'; 
%                 data.opty = interp1(obj.opt_t, obj.opt.datah.y', obj.wam_t', 'linear', 'extrap')'; 
%                 data.optz = interp1(obj.opt_t, obj.opt.datah.z', obj.wam_t', 'linear', 'extrap')'; 
                data.optx = interp1(obj.opt_t, obj.opt.datah.x', sample_t', 'linear')'; 
                data.opty = interp1(obj.opt_t, obj.opt.datah.y', sample_t', 'linear')'; 
                data.optz = interp1(obj.opt_t, obj.opt.datah.z', sample_t', 'linear')'; 
                % in case of obj.opt_t has nan values
%                 opt_t1 = interp1(1:length(obj.opt_t(~isnan(obj.opt_t))), obj.opt_t(~isnan(obj.opt_t)), 1:length(obj.opt_t), 'linear', 'extrap'); % not guareentee!  
%                 obj.opt_t(1) = opt_t1(1); %??? possible???
%                 data.optx = interp1(obj.opt_t, obj.opt.datah.x', obj.wam_t', 'spline', 'extrap')'; 
%                 data.opty = interp1(obj.opt_t, obj.opt.datah.y', obj.wam_t', 'spline', 'extrap')'; 
%                 data.optz = interp1(obj.opt_t, obj.opt.datah.z', obj.wam_t', 'spline', 'extrap')'; 
                
                nmarkers = size(data.optx,1);
                if (ifplot)
                    clf; 
                    for i = 1:nmarkers
                        axh(i) = subplot(nmarkers,1,i); hold on; grid on; 
                        %plot(obj.opt.datah.bkt,obj.opt.datah.x(i,:),'marker', '.', 'Color', 'r');
                        plot(sample_t, data.optx(i,:), '.', 'Color', 'b');
                        
                        %plot(obj.opt.datah.bkt,obj.opt.datah.y(i,:),'marker', '.', 'Color', 'r');
                        plot(sample_t, data.opty(i,:), '.', 'Color', 'b');
                        
                        %plot(obj.opt.datah.bkt,obj.opt.datah.z(i,:),'marker', '.', 'Color', 'r');
                        plot(sample_t, data.optz(i,:), '.', 'Color', 'b');
                    end
                    
                end
                obj.data.optx = data.optx;
                obj.data.opty = data.opty;
                obj.data.optz = data.optz;
                obj.data.ox = zeros([3, size(data.optx, [2,1])]);
                obj.data.ox(1,:,:) = data.optx';
                obj.data.ox(2,:,:) = data.opty';
                obj.data.ox(3,:,:) = data.optz';
            else 
                obj.data.ox = zeros(1,0,3);
            end
            
            if (~isempty(obj.emg_t) && ~isempty(obj.emg_h))
                obj.data.emg = interp1(obj.emg_t', obj.emg_h', sample_t', 'linear', nan)'; % for data not collected, set to nan directly
%                 obj.data.emgevl=interp1(obj.emg_t', obj.emg_evl', sample_t', 'linear', 'extrap')';
                
                ifplot = 0;
                if (ifplot)
                    clf;
                    for i = 1:8
                        axh(i) = subplot(8,1,i); hold on; grid on;
                        plot(obj.emg_t, obj.emg_h(i,:), 'marker', '.', 'Color', 'r');
                        plot(obj.data.t, obj.data.emg(i,:), 'marker', '.', 'Color', 'b');
%                         ylim([-5000 5000]);
                        ylim([-1 1]*5);
                    end
                    linkaxes(axh, 'xy');
                end
            end

            % get from trial specific...

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
        function delay_idx = getDelayedTrialIdx(obj)
            % delay_idx = getDelayTrialIdx(obj) 
            % Return the delayed trials indexes.
            % The delay was detected in some trials, at the release point, 
            % it did not release immediately, it will delay the data a bit. 
            delay_idx = [];
            
            trial = obj.trials;
            for trial_i = 1:length(trial)
                
                ts_valid = [1:7];
                idx = ismember(trial(trial_i).data.ts, ts_valid);
        
                ts_interest = 5; % the first 5 is the release point 
                release_idx = find(trial(trial_i).data.ts(idx)==ts_interest); 
                if (isempty(release_idx)) % no ts5
                    continue;
                end
                release_idx = release_idx(1);
                t_release = trial(trial_i).data.t(idx);
                t_release = t_release(release_idx);
                t_shift = trial(trial_i).data.t(idx) - t_release; 
                Fp = trial(trial_i).data.Fp(:,idx); 
                x = trial(trial_i).data.x(:,idx); 
                v = trial(trial_i).data.v(:,idx); 
                f = trial(trial_i).data.f(:,idx); 
                
                v_threshold = 5e-4; 
                val = sum(v(:,release_idx + 2).^2);
                ifplot = 0;
                if val < v_threshold
                    delay_idx = [delay_idx, trial_i];
                    ifplot = 1;
                end
                
                if (ifplot)
                    figure(); 
                    axh(1) = subplot(4,1,1);     
                    grid on;
                    hold on;
                    plot(t_shift, Fp, 'Marker', '.'); 
                    axh(2) = subplot(4,1,2);     
                    grid on;
                    hold on;
                    plot(t_shift, x, 'Marker', '.');
                    axh(3) = subplot(4,1,3);     
                    grid on;
                    hold on;
                    plot(t_shift, v, 'Marker', '.');
                    plot(t_shift(release_idx), v(:,release_idx), 'Marker', 'o', 'MarkerSize', 5);
                    plot(t_shift(release_idx+2), v(:,release_idx+2), 'Marker', 'o', 'MarkerSize', 5);
                    axh(4) = subplot(4,1,4);     
                    grid on;
                    hold on;
                    plot(t_shift, f, 'Marker', '.');
                    linkaxes(axh, 'x');
                    xlim([-0.02 0.04]);
                    sgtitle(['trial' num2str(trial_i)]);
                end
            end
            
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
        function recognizeBadEmgTrials(obj, vThreshold)
            % recognize which trial have abnormal EMG, and write that trial
            % is not able to use.  
            % Write the trial list in the config file:
            %  manualSetEMGTrialf.conf

            % 1. Specify the threshold
            if (~exist('vThreshold', 'var'))
                vThreshold = 0.1; % mV
            end
            % 2. Get the EMG bad time series 
            t_list = obj.emg.scanThresholdCrossing(vThreshold); 

            % 3. Convert the EMG bad time into the trial 
            trials_list = cell(8,1);
            for ch_i = 1:8
                if isempty(t_list{ch_i})
                    continue
                end
                % 2.2 for each trial
                trials_list_tmp = [];
                for trial_i = 1:length(obj.trials)  
                    t_stt = obj.trials(trial_i).bgn_t; % ...? start time
                    t_edn = obj.trials(trial_i).edn_t; % ...? end time 
                    % if there is time in the zone 
                    if (sum(t_list{ch_i}>t_stt & t_list{ch_i}<t_edn))
                        trials_list_tmp = [trials_list_tmp trial_i];
                    end
                end
                trials_list{ch_i} = trials_list_tmp;
            end
                % concatinate all the trials in this session
                trials_list_all = unique([trials_list{:}]);
                if (isempty(trials_list_all))
                    return; % no trials should be marked as 'wrong EMG'
                end

            % 4. Read the *.conf, check the *.conf, write the *.conf
                % read the *.conf
                filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetEMGTrialsf.conf';
                fid = fopen(filename);
                C = textscan(fid, '%s\n','CommentStyle','#');
                fclose(fid);
                for li = 1:size(C{1}, 1)
                    str = C{1}{li};
                    freadtmp = textscan(str,'%d,%s');
                    
                    ss_num = freadtmp{1}(1);
                    if ss_num ~= obj.ssnum
                        continue;
                    else
                        % this part is not really useful as it read trials
                        % (unnessesary for the current part).
                        ch_trials_str = freadtmp{2}{1};
                        trials_lists = textscan(ch_trials_str,'%s', 'Delimiter',';');
                        for ch_i = 1:8
                            ch_trials{ch_i} = [];
                            if (~isempty(trials_lists{1}{ch_i}))
                                 strtmp = textscan(trials_lists{1}{ch_i},'%d', 'Delimiter',',');
                                 ch_trials{ch_i} = strtmp{1}';
                            end
                        end
                    end
                    trial_str = [ch_trials{:}];
                end
                if exist('trial_str', 'var')    % defined by the config file
                    trials_idx = double(trial_str)';
                    clear trial_str;
                else % default value
                    trials_idx = [];    % The right order
                end

                if ~isempty(trials_idx)
                    fprintf('SessionScan::RecognizeBadEmgTrials: already marked tirals!\n');
                    return
                end
                % write the *.conf
                str_towrite = [];
                for ch_i = 1:length(trials_list)
                    for t_i = 1:length(trials_list{ch_i})
  %                     str_towrite = [str_towrite ',' num2str(trials_list_all(t_i))];
                        str_towrite = [str_towrite num2str(trials_list{ch_i}(t_i)) ','];
                    end
                    if (~isempty(trials_list{ch_i})) 
                        str_towrite = [str_towrite(1:end-1) ';']; % end-1 for remove end ','
                    else
                        str_towrite = [str_towrite ';'];
                    end
                end
                fid = fopen(filename, 'a');
                fprintf(fid,'%s',[num2str(obj.ssnum), ',', str_towrite]);
                fclose(fid);
        end
        function trials_list = readBadEmgTrials(obj)
            % read which trial have abnormal EMG, from config file:
            %  manualSetEMGTrialf.conf

            % 0. default value, blank
            trials_list = cell(8,1);
            % 1. Read the *.conf, check the *.conf, write the *.conf
                % read the *.conf
                filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetEMGTrialsf.conf';
                fid = fopen(filename);
                C = textscan(fid, '%s\n','CommentStyle','#');
                fclose(fid);


                for li = 1:size(C{1}, 1)
                    str = C{1}{li};
                    freadtmp = textscan(str,'%d,%s');
                    
                    ss_num = freadtmp{1}(1);
                    if ss_num ~= obj.ssnum
                        continue;
                    else
                        % this part is not really useful as it read trials
                        % (unnessesary for the current part).
                        ch_trials_str = freadtmp{2}{1};
                        trials_lists = textscan(ch_trials_str,'%s', 'Delimiter',';');
                        for ch_i = 1:8
                            if (~isempty(trials_lists{1}{ch_i}))
                                 strtmp = textscan(trials_lists{1}{ch_i},'%d', 'Delimiter',',');
                            else 
                                strtmp{1} = [];
                            end
                            ch_trials{ch_i} = strtmp{1}';
                        end
                        trials_list = ch_trials;
                    end
                end
            % have bug if there are repeat ss_num in the file
        end
        function obj = dealingEMGOpenLoopNoise(obj)
            % deal with EMG open-loop noise by set specific trials data to Nan
            % read the config from config file 
            % the trials will be set according to the pre-setted config
            trials_list = readBadEmgTrials(obj);
            
            for ch_i = 1:8
                for trial_i = 1:length(trials_list{ch_i})
                    % for specific trial, nan the emg data
                    trial_idx = trials_list{ch_i}(trial_i);
                    obj.trials(trial_idx).data.emg(ch_i,:) = nan;
                    try 
                        obj.trials(trial_idx).data.emgevl(ch_i,:) = nan;
                    end
                end
            end
        end
        function obj = dealingOPTMarkersErr(obj)
            rot_sessions = [4218 4219 4220 ...
                            4225 4226 ...
                            4237 4238 4239 ...
                            4303 4304 ... 
                            4313 4314 ...
                            4328 4329 ...
                            4339 4340 4341 ...
                            4354 4355 4356 ...
                            ];
            % replace the recording error trial with the closest other
            % trial data 
            % especially useful for marker2(elbow) and marker3(shoulder)

            % 1. read bad trials
            marker_val = [obj.trials.opt_v]; % valid
            markers_list = [];
            trials_list = [];
            marker_max = 3; 
            if (sum(rot_sessions==obj.ssnum))
                marker_max = 2; % 3rd marker cannot be recorded.
            end
            for marker_i = 1:marker_max
                trials_list_tmp = find(marker_val(marker_i,:) == 0);
                trials_nstt_tmp = find(isnan([obj.trials.tarF]));
                trials_list_tmp = setdiff(trials_list_tmp, trials_nstt_tmp);
                markers_list_tmp = ones(size(trials_list_tmp))*marker_i;
                markers_list = [markers_list, markers_list_tmp];
                trials_list = [trials_list, trials_list_tmp];
            end
            % the return should be an n-by-1 values (marker), and an n-by-1
            % values (trials_list)
            % 2. for current trial, find the similarest trial. 
            trials_all = find([obj.trials.outcome] == 1);
            
            for trial_i = 1:length(trials_list) % not consider 1st trial
                clear data_org
                trial_i_tmp = trials_list(trial_i); 
                marker_i_tmp = markers_list(trial_i);
                if trial_i_tmp == 1
                    continue
                end
                ifplot = 1;
                if (ifplot)     % see if this trial has a lost data
                    clf; 
                    subplot(3,1,1); hold on;
                    plot(obj.trials(trial_i_tmp).data.t_shift,...
                        reshape(obj.trials(trial_i_tmp).data.ox(1,:,:), size(obj.trials(trial_i_tmp).data.ox(1,:,:), [2,3])), '.');
                    legend('marker1', 'marker2', 'marker3');
                    ylabel('x'); 
                    subplot(3,1,2); hold on;
                    plot(obj.trials(trial_i_tmp).data.t_shift,...
                        reshape(obj.trials(trial_i_tmp).data.ox(2,:,:), size(obj.trials(trial_i_tmp).data.ox(2,:,:), [2,3])), '.');
                    ylabel('y'); 
                    subplot(3,1,3); hold on;
                    plot(obj.trials(trial_i_tmp).data.t_shift,...
                        reshape(obj.trials(trial_i_tmp).data.ox(3,:,:), size(obj.trials(trial_i_tmp).data.ox(3,:,:), [2,3])), '.');
                    ylabel('z'); 
                    xlabel('time (s)');
                end
                marker_lost = markers_list(trial_i); 
                trial_lost = trial_i_tmp; % the index of trials_lost                
                trials_rest= setdiff(trials_all, unique(trials_list)); %? how to get a array with certain trials that 'able to use'?
                %find the possible trials for the current trial & marker
                trials_rest_cond = ...
                    [obj.trials(trials_rest).tarF] == obj.trials(trial_i_tmp).tarF & ...
                    [obj.trials(trials_rest).tarL] == obj.trials(trial_i_tmp).tarL; 
                
                % if there is the same condition 
                if sum(trials_rest_cond)>1
                    trials_rest_cond_idx = trials_rest(trials_rest_cond);
                else
                    if markers_list(trial_i) == 1
                        disp('DealingWithOPTRecordingError: cannot find suitable replacement for Marker1!')
                        trials_rest_cond_idx = []; % should I do this...? cg20221012
                    else % marker 2 or 3 
                        trials_rest_cond = ...
                            [obj.trials(trials_rest).tarL] == obj.trials(trial_i_tmp).tarL; 
                            trials_rest_cond_idx = trials_rest(trials_rest_cond);
                    end

                    if sum(trials_rest_cond) == 0
                        % rare case no data in one condition
                        trials_rest_cond = true(size(trials_rest));
                        trials_rest_cond_idx = trials_rest(trials_rest_cond);
                    end
                end
                t_idx = obj.trials(trial_lost).data.t_shift>-0.5 & obj.trials(trial_lost).data.t_shift<0;
                ox_contnan = isnan(obj.trials(trial_lost).data.ox(marker_lost,t_idx,1)); % has nan value

                if sum(ox_contnan) == 0
                    % if the position reading is valid before release
                    % find depend on the position of the 0.5s before release
                    x_bef_rel = mean(obj.trials(trial_lost).data.ox(:,t_idx,marker_lost),2);
                    x_bef_rel_others = zeros(3,sum(trials_rest_cond));
                    for trial_j = 1:length(trials_rest_cond_idx)
                        trial_subs = trials_rest_cond_idx(trial_j);     % trial substitute
                        t_idx_j = obj.trials(trial_subs).data.t_shift>-0.5 & obj.trials(trial_subs).data.t_shift<0;
                        x_bef_rel_others(:,trial_j) = mean(obj.trials(trial_subs).data.ox(:,t_idx_j,marker_lost),2);
                    end
                    % find the cloest one
                    x_pos_dist = x_bef_rel_others - repmat(x_bef_rel,1,length(trials_rest_cond_idx));
                    x_pos_dist_card = vecnorm(x_pos_dist);
                    [minval, minidx] = min(x_pos_dist_card);
                    %
                    trial_idx_sup = trials_rest_cond_idx(minidx);
                    % interp the old data into the new
                    data_org(:,:,marker_i_tmp) = obj.trials(trial_lost).data.ox(:,:,marker_i_tmp);
                    trial_lost;
                    if (~isempty(trial_idx_sup))
                        obj.trials(trial_lost).data.ox(:,:,marker_i_tmp) = (interp1(...
                            obj.trials(trial_idx_sup).data.t_shift', obj.trials(trial_idx_sup).data.ox(:,:,marker_i_tmp)',...
                            obj.trials(trial_lost).data.t_shift', 'linear', 'extrap'))';
                    end

                else
                    % else, if the position reading is nan before release
                    % use trial 1 to suppliment
                    %trials_rest_cond_idx
                    trial_idx_sup = trials_rest_cond_idx(1);
                    data_org(:,:,marker_i_tmp) = obj.trials(trial_lost).data.ox(:,:,marker_i_tmp);
                    obj.trials(trial_lost).data.ox(:,:,marker_i_tmp) = (interp1(...
                        obj.trials(trial_idx_sup).data.t_shift', obj.trials(trial_idx_sup).data.ox(:,:,marker_i_tmp)',...
                        obj.trials(trial_lost).data.t_shift', 'linear', 'extrap'))';
                end
%                 ifplot = 0;
                if(ifplot)
                    clear axhtmp lnhtmp;
                    axhtmp(1) = subplot(3,2,1);
                    plot(obj.trials(trial_idx_sup).data.t_shift,obj.trials(trial_idx_sup).data.ox(1,:,marker_i_tmp), 'b.');
                    title('suppliment trial');
                    ylabel('x');
                    axhtmp(3) = subplot(3,2,3);
                    plot(obj.trials(trial_idx_sup).data.t_shift,obj.trials(trial_idx_sup).data.ox(2,:,marker_i_tmp), 'b.');
                    ylabel('y');
                    axhtmp(5) = subplot(3,2,5);
                    plot(obj.trials(trial_idx_sup).data.t_shift,obj.trials(trial_idx_sup).data.ox(3,:,marker_i_tmp), 'b.');
                    ylabel('z');
                    xlabel('t (s)');
                    
                    axhtmp(2) = subplot(3,2,2); hold on;
                    lnhtmp{1} = plot(obj.trials(trial_lost).data.t_shift,data_org(1,:,marker_i_tmp), 'g.');
                    lnhtmp{2} = plot(obj.trials(trial_lost).data.t_shift,obj.trials(trial_lost).data.ox(1,:,marker_i_tmp), 'r.');
                    legend([lnhtmp{1}(1), lnhtmp{2}(1)], 'raw', 'replaced');
                    title('lost trial reconstruct');ylabel('x');
                    axhtmp(4) = subplot(3,2,4); hold on;
                    lnhtmp{1} = plot(obj.trials(trial_lost).data.t_shift,data_org(2,:,marker_i_tmp), 'g.');
                    lnhtmp{2} = plot(obj.trials(trial_lost).data.t_shift,obj.trials(trial_lost).data.ox(2,:,marker_i_tmp), 'r.');
                    ylabel('y');
                    axhtmp(6) = subplot(3,2,6); hold on;
                    lnhtmp{1} = plot(obj.trials(trial_lost).data.t_shift,data_org(3,:,marker_i_tmp), 'g.');
                    lnhtmp{2} = plot(obj.trials(trial_lost).data.t_shift,obj.trials(trial_lost).data.ox(3,:,marker_i_tmp), 'r.');
                    ylabel('z');
                    xlabel('t (s)');
                    linkaxes(axhtmp, 'x');
                    xlim([-0.5 2]);
                    sgtitle({'reconstructing data on other trial', ...
                        ['ss' num2str(obj.ssnum)  ...
                        ' tr' num2str(trial_i_tmp) ...
                        ' mkr' num2str(marker_i_tmp)]});
                end



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
            % This line is for the multiple perturbation 
            trial_perturbs(trial_perturbs~=0 & trial_perturbs~=2) = 1;
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
                            cellsmat{sf,tl_i,t_i,3} = t_tmp.export_as_formatted(ifplot);
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
                            cellsmat{sf,1,t_i,p_i} = t_tmp.export_as_formatted(ifplot);  % each trial
%                             xlim([-5 -4])
                            %ifplot = true;
%                             if (ifplot)
%                                 subplot(2,1,1);
%                                 plot(cellsmat{sf,1,t_i, p_i}.t, cellsmat{sf,1,t_i, p_i}.x(2,:));
%                                 subplot(2,1,2);
%                                 plot(cellsmat{sf,1,t_i, p_i}.t, cellsmat{sf,1,t_i, p_i}.f(2,:));
%                             end
                        end
                    end
                end
            end
            
            % plot out the time skew
%             ifplot = 0; %-test
            
%             if (ifplot) 
%                 plt_offset = 2e-3/10;
%                 for p_i = 1:3
%                     subplot(1,3,p_i); hold on;
%                     switch p_i
%                         case 1
%                             title('no pert');
%                         case 2
%                             title('pulse pert');
%                         case 3
%                             title('stoc pert');
%                     end
%                     if (p_i ~= 3)
%                         for t_i = 1:length(t_idx{p_i})
%                             if(isempty(cellsmat{1,t_i, p_i}))
%                                 continue;
%                             end
%                             dt = diff(cellsmat{1,t_i,p_i}.t);
%                             dt = [dt(1) dt];
%                             plot(t_i*plt_offset + dt);
%                         end
%                     else
%                         if isempty(t_idx{p_i}) 
%                             continue
%                         end
%                         tiofst = 0; % plot offset
%                         for tl_i = 1:size(cellsmat, 2)
%                             for t_i = 1:length(cellsmat(1,tl_i,:,p_i))
%                                 if ~isempty(cellsmat{1,tl_i,t_i,p_i})
%                                 dt = diff(cellsmat{1,tl_i,t_i,p_i}.t);
%                                 dt = [dt(1) dt];
%                                 tiofst = tiofst+1;
%                                 plot(tiofst*plt_offset + dt);
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
        end
        function [cellsmat, paramsmat] = export_as_formatted_hybridss(obj, ifplot)
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
            %   emg: 8-by-N matrix, emg data
            pert_trials = [obj.trials.ifpert];
            t_idx = cell(1,3);
            for p_i = 1:3
                 t_idx{p_i} = find(pert_trials==p_i-1); % 0,nopert; 1, pulse; 2, stoc
            end
            
%             t_idx{3} = t_idx{3}(2:end); % works for the hybrid pert, to avoid error
            % if trials are not enough, commit the upper line
            t_idx{3} = setdiff(t_idx{3}, find([obj.trials.outcome]==0));
            cellsmat = cell(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),3);
            paramsmat =cell(max([length(t_idx{1}), length(t_idx{2}), length(t_idx{3})]),3);
            if ~isempty(t_idx{3}) % stoc-perturbed trials. 
                % stoc trials are in the same sessions for step trials (for
                % spring testing)
                for tl_i = 1:length(obj.tarLs)
                    %trial_list = setdiff(find([obj.trials.tarL] == obj.tarLs(tl_i)),1);
                    % exception dealing 
                    if (obj.ssnum >= 3803 && obj.ssnum <= 3830 && obj.tarLs(tl_i) == 0.25) % wrong trial
                        continue
                    end
                   % trial_list = find([obj.trials.tarL] == obj.tarLs(tl_i)& [obj.trials.outcome] == 1);% for subject test 
                   trial_list = find([obj.trials.outcome] == 1);% for spring test 
                    trial_list = setdiff(trial_list,1);
                    for t_i = 1:length(t_idx{3})
                        trial_idx = trial_list(trial_list==t_idx{3}(t_i));
                        t_tmp = obj.trials(trial_idx);
                        %cellsmat{tl_i,t_i,3} = t_tmp.export_as_formatted;
                        cellsmat{t_i,3} = t_tmp.export_as_formatted(ifplot);  
                        paramsmat{t_i,3} = t_tmp.export_trial_params();
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
                    cellsmat{t_i,p_i} = t_tmp.export_as_formatted(ifplot);  % each trial
                    paramsmat{t_i,p_i} = t_tmp.export_trial_params();
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
%             plot(obj.time, obj.Data.Position.Actual(:,1), 'b-o');
            plot(obj.data.t, obj.data.x(1,:), '.');
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
        function fh = plotTaskEndpointPosition_all_opt(obj)
            % plot the xyz-axis position throughout this session
            fh = figure();
            axis_name = 'xyz';
            for ai = 1:3 % xyz
                axh(ai)=subplot(3,1,ai);hold on;grid on;
%                 plot(obj.time, obj.Data.Position.Actual(:,ai), 'b-o');
                switch ai
                    case 1
                        plot(obj.data.t, obj.data.optx(1,:), '.');
                    case 2
                        plot(obj.data.t, obj.data.opty(1,:), '.');
                    case 3
                        plot(obj.data.t, obj.data.optz(1,:), '.');
                end
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
            plot(obj.data.t, obj.data.v(1,:), '.');
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
%             plot(obj.time, obj.force(1,:), 'b-o');
            plot(obj.data.t, obj.data.f(1,:), '.');
            legend('msg', 'hi-sp');
            xlabel('time (s)');
            ylabel('force (N)');
            title('y axis force')
            
            axh(2) = subplot(2,1,2); hold on;
%             plot(obj.time, sqrt(obj.force(1,:).^2+...
%                                 obj.force(2,:).^2+...
%                                 obj.force(3,:).^2), 'b-o');
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

        function fh = plotForceEMGEVL(obj, fh)
            % plot the force as well as the EMG
            % aim to show the 'MVF'
            if (~exist('fh', 'var'))
                fh = figure();
            else
                fh = figure(fh);
            end
            axh(1) = subplot(5,1,1);
            plot(obj.data.t(1:end-1) - obj.data.t(1), obj.data.f(1,1:end-1), 'linewidth', 2);
            grid on;
            xlabel('time (s)');
            ylabel('force (N)');

            axh(2) = subplot(5,1,2); hold on;
            plot(obj.data.t(1:end-1) - obj.data.t(1), obj.data.emgevl(1,1:end-1), 'r');
            plot(obj.data.t(1:end-1) - obj.data.t(1),-obj.data.emgevl(2,1:end-1), 'b');
            grid on;
            xlabel('time (s)');
            ylabel('EMG (mV)');

            axh(3) = subplot(5,1,3); hold on;
            plot(obj.data.t(1:end-1) - obj.data.t(1), obj.data.emgevl(3,1:end-1), 'r');
            plot(obj.data.t(1:end-1) - obj.data.t(1),-obj.data.emgevl(4,1:end-1), 'b');
            grid on;
            xlabel('time (s)');
            ylabel('EMG (mV)');

            axh(4) = subplot(5,1,4); hold on;
            plot(obj.data.t(1:end-1) - obj.data.t(1), obj.data.emgevl(5,1:end-1), 'r');
            plot(obj.data.t(1:end-1) - obj.data.t(1),-obj.data.emgevl(6,1:end-1), 'b');
            grid on;
            xlabel('time (s)');
            ylabel('EMG (mV)');

            axh(5) = subplot(5,1,5); hold on;
            plot(obj.data.t(1:end-1) - obj.data.t(1), obj.data.emgevl(7,1:end-1), 'r');
            plot(obj.data.t(1:end-1) - obj.data.t(1),-obj.data.emgevl(8,1:end-1), 'b');
            grid on;
            xlabel('time (s)');
            ylabel('EMG (mV)');

            linkaxes(axh, 'x');
        end
        function fh = plotForceEMG(obj, fh)
            % plot the force as well as the EMG 
            % aim to show the 'MVF'  
            if (~exist('fh', 'var'))
                fh = figure();
            else
                fh = figure(fh); 
            end
            axh(1) = subplot(5,1,1);
            plot(obj.data.t(1:end-1), obj.data.f(1,1:end-1));
            grid on;


            axh(2) = subplot(5,1,2); hold on;
            plot(obj.data.t(1:end-1), obj.data.emg(1,1:end-1), 'r');
            plot(obj.data.t(1:end-1),-obj.data.emg(2,1:end-1), 'b');


            axh(3) = subplot(5,1,3); hold on;
            plot(obj.data.t(1:end-1), obj.data.emg(3,1:end-1), 'r');
            plot(obj.data.t(1:end-1),-obj.data.emg(4,1:end-1), 'b');


            axh(4) = subplot(5,1,4); hold on;
            plot(obj.data.t(1:end-1), obj.data.emg(5,1:end-1), 'r');
            plot(obj.data.t(1:end-1),-obj.data.emg(6,1:end-1), 'b');


            axh(5) = subplot(5,1,5); hold on;
            plot(obj.data.t(1:end-1), obj.data.emg(7,1:end-1), 'r');
            plot(obj.data.t(1:end-1),-obj.data.emg(8,1:end-1), 'b');


            linkaxes(axh, 'x');
        end
        function fh = plotForceEMGEVL_overlap(obj, fh) 
            % plot the force as well as the EMG 
            % aim to show the 'MVF'  
            if (~exist('fh', 'var'))
                fh = figure();
            else
                fh = figure(fh); 
            end
            axh(1) = subplot(9,1,1);
            plot(obj.data.t(1:end-1), obj.data.f(1,1:end-1));
%             ylim([-30 30]);
            ylim([-200 200]);


            axh(2) = subplot(9,1,2); hold on;
            plot(obj.data.t(1:end-1), obj.data.emg(1,1:end-1), 'r');
            plot(obj.data.t(1:end-1), obj.data.emgevl(1,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            axh(3) = subplot(9,1,3); hold on;
            plot(obj.data.t(1:end-1),obj.data.emg(2,1:end-1), 'r');
            plot(obj.data.t(1:end-1),obj.data.emgevl(2,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            axh(4) = subplot(9,1,4); hold on;
            plot(obj.data.t(1:end-1),obj.data.emg(3,1:end-1), 'r');
            plot(obj.data.t(1:end-1),obj.data.emgevl(3,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            axh(5) = subplot(9,1,5); hold on;
            plot(obj.data.t(1:end-1),obj.data.emg(4,1:end-1), 'r');
            plot(obj.data.t(1:end-1),obj.data.emgevl(4,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            axh(6) = subplot(9,1,6); hold on;
            plot(obj.data.t(1:end-1), obj.data.emg(5,1:end-1), 'r');
            plot(obj.data.t(1:end-1), obj.data.emgevl(5,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            axh(7) = subplot(9,1,7); hold on;
            plot(obj.data.t(1:end-1),obj.data.emg(6,1:end-1), 'r');
            plot(obj.data.t(1:end-1),obj.data.emgevl(6,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            axh(8) = subplot(9,1,8); hold on;
            plot(obj.data.t(1:end-1), obj.data.emg(7,1:end-1), 'r');
            plot(obj.data.t(1:end-1), obj.data.emgevl(7,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            axh(9) = subplot(9,1,9); hold on;
            plot(obj.data.t(1:end-1),obj.data.emg(8,1:end-1), 'r');
            plot(obj.data.t(1:end-1),obj.data.emgevl(8,1:end-1), 'b', 'linewidth', 2);
            ylim([-1 1]);

            linkaxes(axh, 'x');
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
        function axh = plotTrialfyPFh(obj, axh)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).data.t_shift, trials(trial_i).data.Fp(2,:));
            end
            xlabel('time');
            ylabel('Pert Force');
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
        
        function axh = plotTrialfyPositionh_all_opt(obj, axh)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            axh_arr = 'xyz';
            trials = obj.trials;
            for axi = 1:3
                axh(axi) = subplot(3,1,axi); hold on;
                for trial_i = 1:length(trials)
                    plot(trials(trial_i).data.t_shift, trials(trial_i).data.ox(axi,:,1));
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
            linkaxes(axh(:), 'x');
        end
        
        
        function axh = plotTrialfyVelocityh_all(obj, axh)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            axh_arr = 'xyz';
            trials = obj.trials;
            for axi = 1:3
                axh(axi) = subplot(3,1,axi); hold on;
                for trial_i = 1:length(trials)
                    plot(trials(trial_i).data.t_shift, trials(trial_i).data.v(axi,:));
                end
                xlim([-1 1]);
                xlabel('time');
                ylabel('velocity (m/s)');
                title(['vel ' axh_arr(axi)]);
            end
            linkaxes(axh, 'xy');
        end
        function axh = plotTrialfyVelocityh(obj, axh, tidx)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
             if ~exist('tidx','var')
                tidx = 1:length(obj.trials);
            end
            hold on;
            trials = obj.trials(tidx);
            for trial_i = 1:length(trials)
                plot(trials(trial_i).data.t_shift, trials(trial_i).data.v(2,:)); % only y position here
            end
            xlabel('time');
            ylabel('velocity');
            title('all trials velocity');
        end
        function axh = plotTrialfyForceh(obj, axh, tidx)
            if nargin < 2
                axh = figure(); hold on;
            else
                figure(axh); hold on;
            end
            if ~exist('tidx','var')
                tidx = 1:length(obj.trials);
            end
            hold on;
            trials = obj.trials(tidx);
            for trial_i = 1:length(trials)
                plot(trials(trial_i).data.t_shift, trials(trial_i).data.f(1,:));
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
        

        function axh = plotEachTrialTsForcePosition(obj) 
            trials_sel = [obj.trials.tarF] == 15 & [obj.trials.tarL] == 0.075;
            trials_list = find(trials_sel);
            for trial_i = 1:length(trials_list)
                trial_idx = trials_list(trial_i);
                fh = figure(); 
                axh(1) = subplot(4,1,1); 
                plot(obj.trials(trial_idx).data.t_shift, obj.trials(trial_idx).data.ts);
                ylabel('State'); ylim([0 8])

                axh(2) = subplot(4,1,2); 
                plot(obj.trials(trial_idx).data.t_shift, obj.trials(trial_idx).data.f(1,:));
                ylabel('F(dir-x) (N)');


                axh(3) = subplot(4,1,3); 
                plot(obj.trials(trial_idx).data.t_shift, obj.trials(trial_idx).data.ox(1,:,1));
                ylabel('x(OPT) (m)');

                axh(4) = subplot(4,1,4); 
                plot(obj.trials(trial_idx).data.t_shift, obj.trials(trial_idx).data.x(1,:));
                ylabel('x(WAM) (m)');

                linkaxes(axh, 'x');
                sgtitle(fh, ['ss' num2str(obj.ssnum) ' tr' num2str(obj.trials(trial_idx).tNo)]);
            end


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
            MID_FT  = 67;       % NETBOX
            MID_WAM = 62;       % ROBOT
            MID_OPTOTRAK = 74;  % OPTOTRAK
            mid_type = [MID_FT, MID_WAM, MID_OPTOTRAK];
            
            %%%%%%%%%%%%%%%%% DATA READING PART %%%%%%%%%%%%%%%%%%%%
            %%%%% 1. read times from blackrock and check it value (ifplot)
            %%%%%%%%%%%%%%%%%%%%
            eventsL = dataTs.eventsL;   % events_label, 
            eventsT = dataTs.eventsT;   % events_time;
            if (isfield(dataTs, 'eventsTrials')) % trial list
                eventTrials= dataTs.eventsTrials;
            end
            
            
            % pulse from the blackrock
            if max(eventsL) == 16 % did not record the right pulse, find the falling edge
                % seperte the data into 1:16 
                eventsLidx = find(eventsL == 16); 
                eventsT_bin = unique(eventsT); % should be same length as eventsL_bin
                eventsL_bin = ones(length(eventsLidx), 16);
                for i = 1:length(eventsLidx)
                    if i == 1 
                        eventsLidx_tmp = 1:eventsLidx(1);
                    else
                        eventsLidx_tmp = eventsLidx(i-1)+1:eventsLidx(i);
                    end
                    binstr = eventsL(eventsLidx_tmp); % a string with some subset of 1:16
                    binstr_0 = setdiff(1:16, binstr); 
                    for j = 1:length(binstr_0)
                        % binstr_0(j) % 1, 2, or 3
                        eventsL_bin(i, binstr_0(j)) = 0;
                    end
                end
                eventsL_bin = [ones(1,16); eventsL_bin];
                eventsL_bin_diff = diff(eventsL_bin);
                eventsL_bin_diff = eventsL_bin_diff(1:end-1,:); % which may contains '-1'; 
                eventsL_bin_ridx = find(sum(eventsL_bin_diff == -1,2))
                % find the falling edge from the eventsL_bin 
                eventsL_falling = [];
                eventsT_falling = [];
                for i = 1:length(eventsL_bin_ridx)
                    eventsL_tmp = find(eventsL_bin_diff(i,:) == -1)
                    eventsT_tmp = eventsT_bin(i);
                    eventsL_falling = [eventsL_falling eventsL_tmp];
                    eventsT_falling = [eventsT_falling eventsT_tmp];
                end
                eventsL = eventsL_falling;
                eventsT = eventsT_falling;
            end
            events_type = unique(eventsL);
            for etypei = 1:length(events_type)
                bk_time{etypei} = eventsT(eventsL == events_type(etypei));
                bk_trials{etypei} = eventTrials(eventsL == events_type(etypei)); % after ss3090
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
            mid_type = [MID_FT, MID_WAM, MID_OPTOTRAK];
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
                try
                [trial_its, idx_msg2, idx_bk2] = intersect(t_msg_trialidx{2}, bk_trials{2}); % FT
                catch 
                    disp('force time sync ERROR! use erroneous time');
                    idx_msg2 = []; 
                    idx_bk2 = []; 
                end
%                 [trial_its, idx_msg3, idx_bk3] = intersect(bk_trials{3}, bk_trials{3}); % OPTOTRAK
                % assuem every FT sync signal has a WAM sync signal
                
                % 2. change to-aligned data into certain trials
                t_interest{1} = t_interest{1}(idx_msg1);
                t_interest{2} = t_interest{2}(idx_msg2);
                
                bk_time{1} = bk_time{1}(idx_bk1);
                if ~isempty(idx_bk2)
                    bk_time{2} = bk_time{2}(idx_bk2);
                else
                    bk_time{2} = [];
                end
                
                % no need this part again as the optotrak sync message has
                % been updated.
                if_OPT = 0;
                if (length(bk_trials)>=3)
                    if obj.ssnum < 3957 
                        if_OPT = 0;
                    else
                        if_OPT = 1;
                    end
%                     [trial_its, idx_msg3, idx_bk3] = intersect(t_msg_trialidx{3}, bk_trials{3}); % OPTOTRAK, buttom off
%                     bk_time{3} = bk_time{3}(idx_bk3); %buttom off
                    [trial_its, idx_msg3, idx_bk3] = intersect(t_msg_trialidx{3}, bk_trials{4}); % OPTOTRAK, buttom on
                    bk_time{3} = bk_time{4}(idx_bk3); % buttom on
                    t_interest{3} = t_interest{3}(idx_msg3);
                end
                
                
            else % old way to deal with the two message do not have the same length problem
                if (length(t_interest{1}) == 1 + length(bk_time{1}))
                    disp('ERROR(MAYBE): Message and pulse size inconsistance, in SessionScan::updateTimeBlackRock()');
                    
                    t_interest{1} = t_interest{1}(1:end-1);
                    t_interest{2} = t_interest{2}(1:end-1);
                    t_interest{3} = t_interest{3}(1:end-1);
                end
             end
            
            ifplot = 0;
            if (ifplot)
                clf;
                hold on;
                plot(bk_time{1}, t_interest{1} - t_interest{1}(1), 'r*'); 
                plot(bk_time{2}, t_interest{2} - t_interest{2}(1), 'b*');
                if (if_OPT)
                    plot(bk_time{3}, t_interest{3} - t_interest{3}(1), 'g*');
                    legend('FT', 'WAM', 'OPTOTRAK');
                else 
                    legend('FT', 'WAM');
                end
                
                xlabel('BK time');
                ylabel('each computer time (shifted)');
            end
            % intropolate each data time to the bk_time;
            
            %%%%%%%%%%%%%% TIME EXPORT PART %%%%%%%%%%%%%%%%
            %%% 1. the FT time
            try
            obj.force_t = interp1(t_interest{1}, bk_time{1}, obj.ft.elapse, 'linear', 'extrap'); 
            obj.ft.brtime = obj.force_t;
            catch 
                if (ismember(obj.ssnum, [3197 3198]))
                    obj.force_t = interp1(t_interest{2}, bk_time{2}, obj.ft.elapse, 'linear', 'extrap'); 
                    obj.ft.brtime = obj.force_t;
                    % this trial do not have the FT pulse, use WAM pulse to
                    % debug...
                end
            end
            
            ifplot = 0;
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
            if (isempty(obj.wam.time)) % When did not save wam file 
                obj.wam.time = wam_time(1):1/500:wam_time(end);
                ernie_timeh = interp1(wam_time, ernie_time, obj.wam.time, 'linear', 'extrap');
                obj.wam_t = interp1(t_interest{2}, bk_time{2}, ernie_timeh, 'linear', 'extrap');
                obj.wam.tp = nan(3,length(obj.wam.time));
                obj.wam.jp = nan(4,length(obj.wam.time));
                obj.wam.tv = nan(3,length(obj.wam.time));
                obj.wam.jt= nan(4,length(obj.wam.time));
                obj.wam.cf = zeros(3,length(obj.wam.time));
                obj.wam.state = nan(1,length(obj.wam.time));
            end
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
            
            %%% 1. the OPT time
%             try
%                 obj.opt_t = interp1(t_interest{3}, bk_time{3}, obj.opt.datah.t, 'linear', 'extrap'); 
%             catch 
%                 display('wrong in opt time!');
%             end
            if isempty(obj.opt) 
                obj.opt_t = []; 
                disp('no OPT data was recorded');
                flag = -4; % ?write flag?
                return
            end
            obj.opt_t = obj.opt.datah.t;
            
            ifplot = 1;
            if (ifplot && if_OPT)
                clf; 
                hold on;
                plot(t_interest{3}, bk_time{3}, 'b*');
                plot(obj.opt.datah.t, obj.opt_t, 'r.');
            end
            
            
            flag = 1;
            
            if (isempty(setdiff(eventsL, 2)))   % error message: no FT pulse recorded! 
                disp('no FT pulse was recorded');
                flag = -3; 
                % not quit, but finish the rest 
            end
            
            return;
        end
        %% update time from the message 
        function obj = forceHighSample(obj, force_obj)%, ft_intm)
            % align force to higher resolution according to a seperate
            % ft_obj file. the seperate ft_obj file should extract from
            % object of SessionScanFT
            % check variable 
            align_frdt = obj.Data.Force.RDTSeq;
            align_Frdt = reshape(force_obj.RDT, 1, length(force_obj.RDT));
            align_time = obj.time;
            %plot(align_frdt, align_time, '*', align_Frdt, align_Time, 'o');
            % aim: find all non-NaN value of frdts, fill the corresponding
            %      time to align_Time, and linearly fill each interval 
            %%% use each interval seperately
             if (length(align_frdt)<2)
                 msg = 'Data error: not enough sample, aborted!';
                 error(msg);
            end
            %align_time_= align_time(~isnan(align_frdt));
%             [~,idxFrdt] = intersect(align_Frdt,align_frdt,'stable'); %??? check this line
%             [~,idxfrdt]          = intersect(align_frdt,align_Frdt,'stable'); %??? check this line
%             align_Time = nan(1, length(align_Frdt));           % ??? can this work
%             align_Time(idxFrdt) = align_time(idxfrdt);
             align_frdt = double(align_frdt);
             align_frdt(align_frdt==0) = nan;
             align_frdt_idx = [min(find(~isnan(align_frdt))) max(find(~isnan(align_frdt)))];
             align_Time = interp1(align_frdt(align_frdt_idx), obj.time(align_frdt_idx), force_obj.RDT, 'linear', 'extrap');
             obj.force_t = align_Time;
            % loop each interval
%             for i = 1:length(idxFrdt)-1
%                 idx_l = idxFrdt(i); 
%                 idx_r = idxFrdt(i+1);
%                 aligned_tmp = interp1q([idx_l, idx_r]', [align_Time(idx_l),align_Time(idx_r)]', (idx_l:idx_r)');
%                 %plot(1:(idx_r-idx_l+1), align_Time(idx_l:idx_r), 'o', 1:(idx_r-idx_l+1), aligned_tmp', '*');
%                 %title(['i = ' num2str(i) ' of ' num2str(length(idxFrdt))]);
%                 align_Time(idx_l:idx_r) = aligned_tmp';
%             end
            ifplot = 0;
            if(ifplot) % validation figure
                %figure();
                clf;
                plot(align_frdt, align_time, 'o', align_Frdt, align_Time, '.');
                legend('reconstructed (fake) time', 'introplated time'); 
                xlabel('RDT seq num'); ylabel('time'); 
                title('validation interp');
                
                clf;
                title('lo- hi- sampled data');
                plot(obj.time, obj.force(2,:), '*'); hold on;
                plot(align_Time, force_obj.force(2,:), '.');
                legend('lo-sample, RTMA', 'hi-sample, reconstructed time');
                
            end
             %obj.force_h = force_obj.force_net;    % not rotated
%              obj.force_h = force_obj.force;         % rotated
%              if isempty(obj.force_t)
%                 disp('Force did not aligned with blackrock signal, align using messages');
%                 obj.force_t = align_Time; % not work after ss 3046
%                 obj.force_t = reshape(obj.force_t, 1, length(obj.force_t));
%              end
        end
    
        %% save and load the data directly
        function saveFile(obj)
        % save the data into the .mat file 
            % ssprop: sessoin properties;
            % data: time-sequence data; 
            % trials: 
            ssprop.ss_num = obj.ssnum;
            ssprop.stateNames = obj.stateNames;
            ssprop.Fs = obj.Fs; 
            data = obj.data; 
            trials = obj.trials;
            filename = sprintf('KingKong.%05d.format.mat', obj.ssnum);
            save(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/', ...
                filename], 'ssprop', 'data', 'trials', '-v7.3');
        end

        function obj = loadData(obj)
        % load the data from the .mat file 
            % ssprop: sessoin properties;
            % data: time-sequence data; 
            % trials: 

            filename = sprintf('KingKong.%05d.format.mat', obj.ssnum);
            load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/', ...
                filename], 'ssprop', 'data', 'trials');

            obj.ssnum = ssprop.ss_num;
            obj.stateNames = ssprop.stateNames;
            obj.Fs = ssprop.Fs; 
            obj.data = data; 
            obj.trials = trials;
        end
    end

end
function [s,f] = readManualSetsf(ssnum)

    %  for manual sucessful trials; 
    filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetTrials.conf';
    fid = fopen(filename);
    %C = textscan(fid, '%s', 'delimiter',sprintf('\n')); 
    C = textscan(fid, '%s\n','CommentStyle','#'); 
    fclose(fid);
    for li = 1:size(C{1}, 1)
        str = C{1}{li};
        freadtmp = textscan(str,'%d,');
        ss_num = freadtmp{1}(1);
        if ss_num ~= ssnum
            continue;
        end
        trials_num = freadtmp{1}(2:end);
    end
    if exist('trials_num', 'var')
        s = double(trials_num)';
        clear trials_num;
    else 
        s = [];
    end
    
    %  for manual failure trials;
    filename = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/config/manualSetTrialf.conf';
    fid = fopen(filename);
    C = textscan(fid, '%s\n','CommentStyle','#'); 
    fclose(fid);
    for li = 1:size(C{1}, 1)
        str = C{1}{li};
        freadtmp = textscan(str,'%d,');
        ss_num = freadtmp{1}(1);
        if ss_num ~= ssnum
            continue;
        end
        trials_num = freadtmp{1}(2:end);
    end
    if exist('trials_num', 'var')
        f = double(trials_num)';
    else 
        f = [];
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
function if_calibSS = findCalibSS(ss_num)
    % remove this into some .conf file in the future
    ss_list = [4315, 4312 ...
        4324, 4327 ...
        4335, 4338 ...
        4349, 4354 ...
        ];
    if sum(ss_num==ss_list)
        if_calibSS = true;
    else
        if_calibSS = false;
    end
end