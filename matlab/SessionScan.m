classdef (HandleCompatible)SessionScan < handle
    %VARSCAN scanning some variables in the formatted data
    %   Detailed explanation goes here

    properties
        %%% session stat
        Data
        time
        trials_num
        duration
        duration_avg
        %%% task targets
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
        %%% other modules
        ft              % object of force
        wam             % object of wam
        force_h         % time for wam seperate data
        force_t         % time for ft seperate data
        wamp_h          % highly sampled from wam
        wamv_h
        wamt_h
        wam_t
        
        %%% other variables
        endpoint0 = [-0.517 0.483 0.001]
        col_vec = [1 0 0
                    0 1 0
                    0 0 1 
                    0 1 1
                    1 0 1
                    1 1 0
                    1 1 1]   % color for plot, rgb cmyk
        %badTrials = [1];       % bad trial, cull in data
        badTrials = [1, 278];       % FOR SS1898
        stateNames = {'Begin', 'Present', 'FrcRamp', 'Move', 'Hold', 'End', 'Reset'};
        
    end
    
    methods
        %%% process
        function obj = SessionScan(ss_num, badTrials) %(inputArg1,inputArg2)
            %VARSCAN Construct an instance of this class
            %   Detailed explanation goes here
            file_name = ['KingKong.0' num2str(ss_num) '.mat']; % an examplary trial
            file_dir = '/Users/cleave/Documents/projPitt/Ballistic_release_data/Formatted';
            fname0 = ([file_dir '/' file_name]);
            flag_progress = 1;      % show something to make me less anxious
            try 
                load(fname);
            catch 
                load(fname0);
            end
            % objects
            noFW_data = 0;
            try
                if (flag_progress)
                    display('Loading raw FT and WAM...');
                end
                obj.ft = SessionScanFT(ss_num);
                obj.wam = SessionScanWam(ss_num);
            catch 
                display('no ft and wam data here! ');
                noFW_data = 1;
            end
            % other data

            obj.Data = Data;
            obj.time = Data.Time;
            TrialNo = Data.TrialNo;
            obj.trials_num = max(TrialNo);
            obj.hand_pos_offset = Data.Position.Center(:,~isnan(Data.Position.Center(1,:)));
            obj.hand_pos_offset = obj.hand_pos_offset(:,1); 
            obj.taskState.Values = Data.TaskStateCodes.Values;
            obj.taskState.Outcome = Data.OutcomeMasks;
            obj.taskState.trialNo = Data.TrialNo;
            obj.force = forceFTconvert(obj);
            obj.duration = max(obj.time);
            targets_idx = ~isnan(obj.Data.TaskJudging.Target(1,:));
            obj.tarRs = unique(obj.Data.TaskJudging.Target(5,targets_idx)); %have 0
            obj.tarLs = setdiff(unique(obj.Data.TaskJudging.Target(6,targets_idx)),0); % deviate 0
            obj.fThs  = setdiff(unique(obj.Data.TaskJudging.Target(4,targets_idx)),0); % deviate 0

            
            % execution functions 
            % trialTimeAverage(obj); % how to use class function?

            % processing 
            obj = convert0toNan(obj);
            trials_all = setdiff(unique(TrialNo), 0);
            if (~isempty(obj.ft))
                if (flag_progress)
                    display('FT High Sample...');
                end
                obj = forceHighSample(obj, obj.ft);
            end
            if (~isempty(obj.wam))
                if (flag_progress)
                    display('WAM High Sample...');
                end
                obj = wamHighSample(obj, obj.wam);
            end
            if (flag_progress)
                    display('Trialfy...');
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
            
        end
        function obj = convert0toNan(obj) % dealing with some Nan-int confliction
            rdt = double(obj.Data.Position.RDT);
            rdt_idx = find([rdt==0]);
            rdt(rdt_idx) = NaN;
            obj.Data.Position.RDT = rdt;
            
            rdt = double(obj.Data.Force.RDTSeq);
            rdt_idx = find([rdt==0]);
            rdt(rdt_idx) = NaN;
            obj.Data.Force.RDTSeq = rdt;
        end
        function trialTimeAverage(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.duration_avg = obj.duration/double(obj.trials_num);
        end
        function obj = forceHighSample(obj, force_obj)
            % align force to higher resolution according to a seperate
            % ft_obj file. the seperate ft_obj file should extract from
            % object of SessionScanFT
            % check variable 
            align_frdt = obj.Data.Force.RDTSeq;
            align_Frdt = force_obj.RDT;
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
            [~,idxFrdt] = intersect(align_Frdt,align_frdt,'stable'); %??? check this line
            [~,idxfrdt]          = intersect(align_frdt,align_Frdt,'stable'); %??? check this line
            align_Time = nan(1, length(align_Frdt));           % ??? can this work
            align_Time(idxFrdt) = align_time(idxfrdt);
            % loop each interval
            for i = 1:length(idxFrdt)-1
                idx_l = idxFrdt(i); 
                idx_r = idxFrdt(i+1);
                aligned_tmp = interp1q([idx_l, idx_r]', [align_Time(idx_l),align_Time(idx_r)]', (idx_l:idx_r)');
                %plot(1:(idx_r-idx_l+1), align_Time(idx_l:idx_r), 'o', 1:(idx_r-idx_l+1), aligned_tmp', '*');
                %title(['i = ' num2str(i) ' of ' num2str(length(idxFrdt))]);
                align_Time(idx_l:idx_r) = aligned_tmp';
            end
            if(0) % validation figure
                figure();
                plot(align_frdt, align_time, 'o', align_Frdt, align_Time, '*');
                legend('RTMA time', 'forceT time'); 
                xlabel('RDT seq num'); ylabel('time'); 
                title('validation interp');
            end
             obj.force_h = force_obj.force;
             obj.force_t = align_Time;
        end
        function obj = wamHighSample(obj, wam_obj)
            % align robot movement to higher resolution according to a
            % seperate wam.obj file. the seperate wam.obj file should
            % extract from SessionScanWam
            
            align_wrdt = obj.Data.Position.RDT;
            align_Wrdt = wam_obj.rdt;
            align_time = obj.time;

             wam_obj = convert0tonan_RDT(wam_obj);
             align_wrdt = obj.Data.Position.RDT;
  
             align_Wrdt = wam_obj.rdt;
             align_time = obj.time;
            % aim: find the max and min non-NaN value of wrdt, and time, 
            %       Apply interval to all align_Wrdt

            if (length(align_wrdt)<2)
                 msg = 'Data error: not enough sample, aborted!';
                 error(msg);
            end
            [~,idxWrdt] = intersect(align_Wrdt,align_wrdt,'stable'); %??? check this line
            [~,idxwrdt]          = intersect(align_wrdt,align_Wrdt,'stable'); %??? check this line
            align_Time = nan(1, length(align_Wrdt));           % ??? can this work
            align_Time(idxWrdt) = align_time(idxwrdt);
            % loop each interval
            for i = 1:length(idxWrdt)-1
                idx_l = idxWrdt(i); 
                idx_r = idxWrdt(i+1);
                aligned_tmp = interp1q([idx_l, idx_r]', [align_Time(idx_l),align_Time(idx_r)]', (idx_l:idx_r)');
                %plot(1:(idx_r-idx_l+1), align_Time(idx_l:idx_r), 'o', 1:(idx_r-idx_l+1), aligned_tmp', '*');
                %title(['i = ' num2str(i) ' of ' num2str(length(idxFrdt))]);
                align_Time(idx_l:idx_r) = aligned_tmp';
            end
            if(0) % validation figure
                figure();
                plot(align_wrdt, align_time, 'o', align_Wrdt, align_Time, '*');
                legend('RTMA time', 'WAM time'); 
                xlabel('RDT seq num'); ylabel('time'); 
                title('validation interp');
            end
            %size(align_Time)
            obj.wamp_h = wam_obj.tp;
            obj.wamv_h = wam_obj.tv;
            obj.wamt_h = wam_obj.jt;
            obj.wam_t = align_Time;
        end
        function force = forceFTconvert(obj) % convert from select into world axis
            force = obj.FTrot_M * obj.Data.Force.Sensor(1:3,:);
        end
        function [resample_t, resample_f] = trialDataResampleFT(obj, trial_idx)
            % [resample_t, resample_f] = trialDataResampleFT(obj, trial_idx)
            % for all trials resample the original data and time
            % origin data and time are Nonuniformly sampled
            % due to the FT is non-normally distributed (code bug of
            % FTmodule -cg), I used resample to get code data
            if nargin<2
                trial_idx = [obj.trials.outcome]==1;
            end
            trial_idx_num = find(trial_idx);
            trials = (obj.trials(trial_idx));
            %display(['Enter function Resample;']);
            tz_bgn = -1.0;
            tz_edn =  0.5;
            resample_freq = 500;   % 500Hz
            resample_t = [tz_bgn: (1/resample_freq): tz_edn];
            resample_t = resample_t(2:end); % looks like this one is loger 1 element than Ty? how to deal withit?
            resample_f = zeros(length(trials), length(resample_t), 3); % resample_value, respectively, x, y, z

            for trial_i = 1:length(trials) % for all trials
                % select specific timezone
                idx_t = trials(trial_i).force_t > tz_bgn & trials(trial_i).force_t <= tz_edn;
                irregTx = trials(trial_i).force_t(idx_t); % problem here, not wanted force_t
                
                for dim_i = 1:3 % x, y, z seperately
                    x = trials(trial_i).force_h(dim_i,idx_t);
                    try
                    [y, Ty] = resample(x,irregTx,resample_freq);            % the non-uniform resample
                    %[y, Ty] = resample(y,irregTx,resample_freq);
                    catch
                        y = [];
                        Ty = [];
                        display(['Unable to resample trial' num2str(trial_i) ' dim' num2str(dim_i)]);
                    end
                    
                    if (0) % visualize the resample result
                        figure(); 
                        hold on;
                        plot(irregTx,x,'.-', Ty,y,'o-')
                        legend('Original','Resampled')
                        plot(Ty);
                    end
                    try
                        resample_f(trial_i,:,dim_i) = y;
                    catch
                    %    display(['trial:' num2str(trial_idx_num(trial_i))]);
                        y = resample(y,size(resample_f,2),length(y));         % may cause time skew here!
                        resample_f(trial_i,:,dim_i) = y;
                    end
                end
                %resample_t = Ty;  % still have 1ms variance between different trials, why?
            end
        end
        function [resample_t, resample_p, resample_v] = trialDataAlignWAM(obj, trial_idx)
            % for all trials align the original data and time into one
            % matrix
            % Assuming they are uniformlly sampled within one trial
            % pick correct ones and save them in the matrix
            if nargin<2
                trial_idx = [obj.trials.outcome]==1;
            end
            trial_idx_num = find(trial_idx);
            trials = (obj.trials(trial_idx));
            display(['Enter function Resample;']);
            tz_bgn = -0.5; %time-zone
            tz_edn =  0.8;
            resample_freq = 500;   % 500Hz, same with WAM
            resample_t = [tz_bgn: (1/resample_freq): tz_edn];
            [~,idx_tmp] = sort(abs(resample_t - 0),2,'ascend'); 
            resample_p = nan(length(trials), length(resample_t), 3); % x, y z
            resample_v = nan(length(trials), length(resample_t), 3); % x, y z
            resample_left_num = sum(resample_t<=0);         % left contain zero
            %resample_left_idx = -resample_left_num+1:0;
            resample_right_num = sum(resample_t>0);
            %resample_right_idx = 1:resample_right_num;
            % resample for all the positions

            for trial_i = 1:length(trials) % for all trials
                wam_t = trials(trial_i).position_t;
                % find the cloest to zero 1
                [~,idx_Tmp] = sort(abs(wam_t - 0),2,'ascend');
                resample_Left_num = idx_Tmp(1); 
                resample_Right_num = length(wam_t) - idx_Tmp(1);
                % pair index
                left_num = min(resample_left_num, resample_Left_num);
                right_num = min(resample_right_num, resample_Right_num);
                resample_left_idx = -left_num + idx_tmp(1)+1;     % avoid 0
                resample_right_idx= right_num + idx_tmp(1);
                resample_Left_idx = -left_num + idx_Tmp(1)+1;
                resample_Right_idx= right_num + idx_Tmp(1);
                % put index in the resampe_p and resample_v
                resample_p(trial_i,resample_left_idx:resample_right_idx,1) = ...
                    trials(trial_i).position_h(1,resample_Left_idx:resample_Right_idx);
                resample_p(trial_i,resample_left_idx:resample_right_idx,2) = ...
                    trials(trial_i).position_h(2,resample_Left_idx:resample_Right_idx);
                resample_p(trial_i,resample_left_idx:resample_right_idx,3) = ...
                    trials(trial_i).position_h(3,resample_Left_idx:resample_Right_idx);
                resample_v(trial_i,resample_left_idx:resample_right_idx,1) = ...
                    trials(trial_i).velocity_h(1,resample_Left_idx:resample_Right_idx);
                resample_v(trial_i,resample_left_idx:resample_right_idx,2) = ...
                    trials(trial_i).velocity_h(2,resample_Left_idx:resample_Right_idx);
                resample_v(trial_i,resample_left_idx:resample_right_idx,3) = ...
                    trials(trial_i).velocity_h(3,resample_Left_idx:resample_Right_idx);
            end

        end
        %%% communicate 
        function obj = generateWamPertData(obj)
            % send data into wam function to help SessionScanWam generate perturbation-only data
            % perturbation rdt already saved in TrialScan
            tarL_list = [];
            fTh_list = [];
            rdt_ranges_all = {};
            i = 0;
            for tarLi = obj.tarLs
                for fThi = obj.fThs
                    i = i+1;
                    trials_idx = [obj.trials.tarL] == tarLi &...
                                    [obj.trials.fTh] == fThi &...
                                    [obj.trials.outcome] == 1;      % only sucessful trials
                    rdt_ranges = [obj.trials(trials_idx).pert_rdt_bgn;...
                                    obj.trials(trials_idx).pert_rdt_edn];
                    tarL_list = [tarL_list tarLi];
                    fTh_list  = [fTh_list fThi];
                    rdt_ranges_all{i} = rdt_ranges;
                end
            end
            obj.wam = obj.wam.concatinateTrials2File(tarL_list, fTh_list, rdt_ranges);
            % for each trial condition, concatinate a structure
            
             % call generateWamPertData()
        end
        %%% plot
        function taskScanTrials(obj)
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
        function taskStateMuskFig(obj)
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
            end
            % notations
            yticks([1 2 3 4 5 6 7]);
            yticklabels(fields_t');
            xlabel('time pts');
            title('states though time');
        end
        function taskJointPosition(obj) % need further changing. 
            % plot 
            figure();
            axish(1) = subplot(2,1,1);
            % plot([diff(obj.Data.Position.JointPosition,1,2)-obj.Data.Position.JointVelocity(:,2:end)]');
            plot(obj.time, obj.Data.Position.JointPosition');
            legend('J1', 'J2', 'J3', 'J4');
            axish(2) = subplot(2,1,2);
            % plot([diff(obj.Data.Position.JointVelocity,1,2)-obj.Data.Position.JointTorque(:,2:end)]');
            plot(obj.Data.Position.JointVelocity');
            legend('J1', 'J2', 'J3', 'J4');
            % notation
            %set(axish(1), 'Ylim', [-0.02 0.02]);
            ylabel(axish(1), 'joints positions');
            %set(axish(2), 'Ylim', [-0.3 0.3]);
            ylabel(axish(2), 'joints velocities');
            xlabel(axish(2), 'time points');
            
        end
        function axh = taskJointPosition_relateve(obj) % Plot relative position. 
            % plot 
            figure();
            axh(1) = subplot(2,1,1);
            % plot([diff(obj.Data.Position.JointPosition,1,2)-obj.Data.Position.JointVelocity(:,2:end)]');
            position_offset = obj.Data.Position.JointPosition(:,~isnan(obj.Data.Position.JointPosition(1,:)));
            position_offset = repmat(position_offset(:,1),1,size(obj.Data.Position.JointPosition,2));
            plot(obj.time, (obj.Data.Position.JointPosition - position_offset)');
            legend('J1', 'J2', 'J3', 'J4');
            axh(2) = subplot(2,1,2);
            % plot([diff(obj.Data.Position.JointVelocity,1,2)-obj.Data.Position.JointTorque(:,2:end)]');
            plot(obj.Data.Position.JointVelocity');
            legend('J1', 'J2', 'J3', 'J4');
            % notation
            %set(axish(1), 'Ylim', [-0.02 0.02]);
            ylabel(axh(1), 'joints positions');
            %set(axish(2), 'Ylim', [-0.3 0.3]);
            ylabel(axh(2), 'joints velocities');
            xlabel(axh(2), 'time points');
            
        end
        function axh = taskForceData(obj, axh)
            if nargin<2
                axh = figure();
            else
                figure(axh); hold on;
            end
            % force = obj.Data.Force.Sensor(1:3,:); 
            force = obj.force;
            plot(obj.time, force');
            ylabel('force (N)');
            xlabel('time');
            legend('x', 'y', 'z'); % remember to alter the axis 
            title('Force data');
        end
        function axh = taskForceDatah(obj, axh)
            if nargin<2
                axh = figure();
            else
                figure(axh); hold on;
            end
            % force = obj.Data.Force.Sensor(1:3,:); 
            force = obj.force_h;
            plot(obj.force_t, force', '.');
            ylabel('force (N)');
            xlabel('time');
            legend('x', 'y', 'z'); % remember to alter the axis 
            title('Force data high sample');
        end
        function taskEndpointPosition_relative(obj)
            figure();
            position = obj.Data.Position.Actual'; 
            % Use first element as offset
            % position_offset = position(:,~isnan(position(1,:)));
            position_offset = obj.hand_pos_offset;
            position_offset = repmat(position_offset(:,1),1,size(position,2));
            plot(obj.time, (position - position_offset)');  
            ylabel('relative endpoint positions');
            xlabel('time points');
            legend('x', 'y', 'z');
            title('relative endpoint positions');
        end
        function axh = taskEndpointPosition(obj)
            axh = figure();
            position = obj.Data.Position.Actual'; 
            plot(obj.time, (position)');  
            ylabel('endpoint positions');
            xlabel('time points');
            legend('x', 'y', 'z');
            title('relative endpoint positions');
        end 
        function axh = taskEndpointPositionh(obj, axh)
            if nargin < 2
                axh = figure();
            else 
                figure(axh); hold on;
            end
            position = obj.wamp_h'; 
            plot(obj.wam_t, (position)', '.');  
            ylabel('endpoint positions');
            xlabel('time points');
            legend('x', 'y', 'z');
            title('relative endpoint positions');
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
        function axh = addmark_STMOV(obj, axh) % add lines showing mov state. 
            % how to addline without add the legend???
            mov_mask = obj.Data.TaskStateMasks.Move;
            mov_diff = [0 diff(mov_mask)]; 
            mov_idx = find((mov_diff == 1) & (mov_mask == 1));
            % plot bars in the axh
            for axh_i = 1:length(axh)
                ylim_range = get(axh(axh_i), 'ylim');
                v= ver('MATLAB'); 
                for ii = 1:length(mov_idx)
                    if str2double(v.Version) <= 9.0 % less than matlab 2018
                        line(axh(axh_i), [mov_idx(ii) mov_idx(ii)], [ylim_range(1) ylim_range(2)]);
                    else
                        xline(axh(axh_i), mov_idx(ii));
                    end
                end
            end
        end
        function axh = plotTrialfyPosition(obj, axh)
            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).time, trials(trial_i).position);
            end
            xlabel('time');
            ylabel('position');
            title('all trials position');
        end
        function axh = plotTrialfyForce(obj, axh)
            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).time, trials(trial_i).force);
            end
            xlabel('time');
            ylabel('force');
            title('all trials force');
        end
        function axh = plotTrialfyForce_xy(obj, axh)
            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            
            hold on;
            trials = obj.trials;
            subplot(2,1,1); hold on;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).time, trials(trial_i).force(1,:));
            end
            xlim([-1 1]);
            xlabel('time');
            ylabel('force x');
            title('force xy (m)');
            subplot(2,1,2); hold on;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).time, trials(trial_i).force(2,:));
            end
            xlim([-1 1]);
            title('force y');
            xlabel('time');
            ylabel('position y (m)');
        end
        function axh = plotTrialfyPositionh(obj, axh)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).position_t, trials(trial_i).position_h);
            end
            xlabel('time');
            ylabel('position');
            title('all trials position');
        end
        function axh = plotTrialfyPositionh_xy(obj, axh)

            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            trials = obj.trials;

            subplot(2,1,1); hold on;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).position_t, trials(trial_i).position_h(1,:));
            end
            xlim([-1 1]);
            xlabel('time');
            ylabel('position');
            title('position x (m)');
            subplot(2,1,2); hold on;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).position_t, trials(trial_i).position_h(2,:));
            end
            xlim([-1 1]);
            title('position y');
            xlabel('time');
            ylabel('position y (m)');
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
                plot(trials(trial_i).force_t, trials(trial_i).force_h);
            end
            xlabel('time');
            ylabel('force');
            title('all trials force raw');
        end
        function axh = plotTrialfyForceh_xy(obj, axh)
            if nargin < 2
                axh = figure();
            else
                figure(axh);
            end
            hold on;
            trials = obj.trials;
            subplot(2,1,1); hold on;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).force_t, trials(trial_i).force_h(1,:));
            end
            xlabel('time');
            ylabel('force x (N)');
            title('force x ');
            xlim([-1 1]);
            subplot(2,1,2); hold on;
            for trial_i = 1:length(trials)
                plot(trials(trial_i).force_t, trials(trial_i).force_h(2,:));
            end
            xlabel('time');
            ylabel('force y (N)');
            title('force y ');
            xlim([-1 1]);
        end
        function axhf = plotMeantrialForce(obj)
            % plot the meaned trial Force according to the task condition
            all_fTH = unique([obj.trials.fTh]);
            all_fTH = all_fTH(~isnan(all_fTH));
            all_tarL = unique([obj.trials.tarL]);
            all_tarL = all_tarL(~isnan(all_tarL));
            all_tarR = unique([obj.trials.tarR]);
            all_tarR = all_tarR(~isnan(all_tarR));
            % assume this session only have x- or y- trials
            if isempty(setdiff(all_tarR, [0,4])) %only y direction
                xyi = 1;
            elseif isempty(setdiff(all_tarR, [2, 6]))
                xyi = 2;
            end
            xy_char = 'xy';
            axhf = figure();
            l_h = [];
            for fTH_i = 1:length(all_fTH)
                for tarL_i = 1:length(all_tarL)
                    col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
                    hold on;
                    %trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
                    trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) ...
                        & [obj.trials.tarL]==all_tarL(tarL_i)...
                        & [obj.trials.outcome]==1; 
                    [resample_t, resample_f] = trialDataResampleFT(obj, trials_idx);
                    force_mean = mean(resample_f(:,:,xyi)); %only y direction
                    force_std = std(resample_f(:,:,xyi));  
                    % mean line
                    l_h = [l_h plot(resample_t, force_mean, 'LineWidth', 3, 'Color', obj.col_vec(col_i,:))];
                    % 1std shade
                    force_up = force_mean + force_std/2;
                    force_dn = force_mean - force_std/2;
                    [axhf, msg] = jbfill(resample_t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
                end
            end
            xlabel('time aligned at MOV signal');
            ylabel([xy_char(xyi) ' dir force (N)']);
            xlim([-0.5 0.4]);
            title('Force signal');
            legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
        end
        function axhp = plotMeantrialPos(obj)
            % plot the meaned trial Position according to the task condition
            % if did not calculate the mean and var, calculate
            all_fTH = unique([obj.trials.fTh]);
            all_fTH = all_fTH(~isnan(all_fTH));
            all_tarL = unique([obj.trials.tarL]);
            all_tarL = all_tarL(~isnan(all_tarL));
            all_tarR = unique([obj.trials.tarR]);
            all_tarR = all_tarR(~isnan(all_tarR));
            % assume this session only have x- or y- trials
            if isempty(setdiff(all_tarR, [0,4])) %only y direction
                xyi = 1;
            elseif isempty(setdiff(all_tarR, [2, 6]))
                xyi = 2;
            end
            xy_char = 'xy';
            % plot position
            axhp = figure();
            l_h = [];
            for fTH_i = 1:length(all_fTH)
                for tarL_i = 1:length(all_tarL)
                    col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
                    hold on;
                    trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
                    [resample_t, resample_p, ~] = trialDataAlignWAM(obj, trials_idx);
                    force_mean = mean(resample_p(:,:,xyi), 'omitnan') - obj.endpoint0(xyi); %only y direction
                    force_std = std(resample_p(:,:,xyi), 'omitnan');  
                    % mean line
                    l_h = [l_h plot(resample_t, force_mean, 'LineWidth', 3, 'Color', obj.col_vec(col_i,:))];
                    % 1std shade
                    force_up = force_mean + force_std/2;
                    force_dn = force_mean - force_std/2;
                    [axh, msg] = jbfill(resample_t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
                end
            end
            xlabel('time aligned at MOV signal');
            ylabel([xy_char(xyi) 'dir position (m)']);
            xlim([-0.2 0.5]);
            title('Position signal');
            legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
        end
        function axhv = plotMeantrialVel(obj)
            % plot the meaned trial Velocity according to the task condition
            % if did not calculate the mean and var, calculate
            all_fTH = unique([obj.trials.fTh]);
            all_fTH = all_fTH(~isnan(all_fTH));
            all_tarL = unique([obj.trials.tarL]);
            all_tarL = all_tarL(~isnan(all_tarL));
            all_tarR = unique([obj.trials.tarR]);
            all_tarR = all_tarR(~isnan(all_tarR));
            % assume this session only have x- or y- trials
            if isempty(setdiff(all_tarR, [0,4])) %only y direction
                xyi = 1;
            elseif isempty(setdiff(all_tarR, [2, 6]))
                xyi = 2;
            end
            xy_char = 'xy';
            % plot velocity
            axhv = figure();
            l_h = [];
            for fTH_i = 1:length(all_fTH)
                for tarL_i = 1:length(all_tarL)
                    col_i = (fTH_i-1)*length(all_tarL) + tarL_i;
                    hold on;
                    trials_idx = [obj.trials.fTh]==all_fTH(fTH_i) & [obj.trials.tarL]==all_tarL(tarL_i);
                    [resample_t, ~, resample_v] = trialDataAlignWAM(obj, trials_idx);
                    force_mean = mean(resample_v(:,:,xyi), 'omitnan'); %only y direction
                    force_std = std(resample_v(:,:,xyi), 'omitnan');  
                    % mean line
                    l_h = [l_h plot(resample_t, force_mean, 'LineWidth', 3, 'Color', obj.col_vec(col_i,:))];
                    % 1std shade
                    force_up = force_mean + force_std/2;
                    force_dn = force_mean - force_std/2;
                    [axh, msg] = jbfill(resample_t, force_up, force_dn, obj.col_vec(col_i,:), obj.col_vec(col_i,:), 1, 0.3);
                end
            end
            xlabel('time aligned at MOV signal');
            ylabel([xy_char(xyi) 'dir velocity (m/s)']);
            xlim([-0.2 0.5]);
            title('Velocity signal');
            legend(l_h, {'10N5cm', '10N10cm', '20N5cm', '20N10cm'});
        end
        function [axhf, axhp, axhv] = plotSameTrial(obj)
            all_fTH = unique([obj.trials.fTh]);
            all_fTH = all_fTH(~isnan(all_fTH));
            all_tarL = unique([obj.trials.tarL]);
            all_tarL = all_tarL(~isnan(all_tarL));
            outcome = [obj.trials.outcome];
            % fTH now only have 2 vals
            % tarL now only have 2 vals
            axhf = figure(); % force figure
            set(axhf, 'Visible', 'off');
            axhp = figure(); % position figure
            set(axhp, 'Visible', 'off');
            axhv = figure(); % velocity figure
            set(axhv, 'Visible', 'off');
            for fTH_i = 1:length(all_fTH)
                fprintf('\n plotting fTH:%2.1f  ', all_fTH(fTH_i));
                for tarL_i = 1:length(all_tarL)
                    fprintf('target: %0.2f', all_tarL(tarL_i));
                    trials_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)));
                    trial_num = sum(trials_idx);
                    trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) & outcome==1); %succeed
                    %trialsS_idx = (([obj.trials.fTh]==all_fTH(fTH_i)) & ([obj.trials.tarL]==all_tarL(tarL_i)) ); %all trials
                    col_i = (fTH_i-1)*length(all_fTH) + tarL_i;
                    if col_i>=7
                        col_i = mod(col_i,8);
                    end
                    for trial_i = find(trialsS_idx)
                        % plot force
                        % data
                        try
                            force = obj.trials(trial_i).force_h;
                            time = obj.trials(trial_i).force_t;
                        catch
                            force = obj.trials(trial_i).force;
                            time = obj.trials(trial_i).time;
                        end
                        % figure
                        set(0, 'CurrentFigure', axhf);
                        % x, y seperately
                        subplot(2,1,1); hold on; plot(time, force(1,:), 'Color', obj.col_vec(col_i,:));
                        subplot(2,1,2); hold on; plot(time, force(2,:), 'COlor', obj.col_vec(col_i,:));
                        % plot position
                        % data
                        try
                            position = obj.trials(trial_i).position_h;
                            time = obj.trials(trial_i).position_t;
                        catch
                            position = obj.trials(trial_i).position;
                            time = obj.trials(trial_i).time;
                        end
                        set(0, 'CurrentFigure', axhp);
                        % x, y seperately
                        subplot(2,1,1); hold on; plot(time, position(1,:), 'Color', obj.col_vec(col_i,:));
                        subplot(2,1,2); hold on; plot(time, position(2,:), 'Color', obj.col_vec(col_i,:));
                        try
                            velocity = obj.trials(trial_i).velocity_h;
                            time = obj.trials(trial_i).position_t;
                        catch
                            velocity = obj.trials(trial_i).velocity;
                            time = obj.trials(trial_i).time;
                        end
                        set(0, 'CurrentFigure', axhv);
                        subplot(2,1,1); hold on; plot(time, velocity(1,:), 'Color', obj.col_vec(col_i,:));
                        subplot(2,1,2); hold on; plot(time, velocity(2,:), 'Color', obj.col_vec(col_i,:));
                    end
                end
            end
            
            xrangeF = [-0.5, 0.5];
            figure(axhf); 
            subplot(2,1,1);
            title('force data x');
            xlim(xrangeF);
            subplot(2,1,2);
            title('force data y');
            xlim(xrangeF);
            
            xrangeP = [-0.2, 0.8];
            figure(axhp);
            subplot(2,1,1);
            title('robot position x');
            xlim(xrangeP);
            subplot(2,1,2);
            title('robot position y');
            xlim(xrangeP);
            
            xrangeV = [-0.2, 0.8];
            figure(axhv);
            subplot(2,1,1);
            title('robot velocity x');
            xlim(xrangeV);
            subplot(2,1,2);
            title('robot velocity y');
            xlim(xrangeV);
            
            set(axhf, 'Visible', 'on');
            set(axhp, 'Visible', 'on');
            set(axhv, 'Visible', 'on');
        end
    end

end

