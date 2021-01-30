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
        %%% other modules
        ft              % object of force
        wam             % object of wam
        force_h         % time for wam seperate data
        force_t         % time for ft seperate data
        wamp_h          % highly sampled from wam
        wamv_h
        wamt_h
        wam_t
        
    end
    
    methods
        %%% process
        function obj = SessionScan(ss_num) %(inputArg1,inputArg2)
            %VARSCAN Construct an instance of this class
            %   Detailed explanation goes here
            file_name = ['KingKong.0' num2str(ss_num) '.mat']; % an examplary trial
            file_dir = '/Users/cleave/Documents/projPitt/Ballistic_release_data/Formatted';
            fname0 = ([file_dir '/' file_name]);
            try 
                load(fname);
            catch 
                load(fname0);
            end
            % objects
            obj.ft = SessionScanFT(ss_num);
            obj.wam = SessionScanWam(ss_num);
            % other data
            obj.time = Data.Time;
            TrialNo = Data.TrialNo;
            obj.duration = max(obj.time);
            obj.trials_num = max(TrialNo);
            obj.Data = Data;
            obj.hand_pos_offset = Data.Position.Center(:,~isnan(Data.Position.Center(1,:)));
            obj.hand_pos_offset = obj.hand_pos_offset(:,1); 
            obj.taskState.Values = obj.Data.TaskStateCodes.Values;
            obj.taskState.Outcome = obj.Data.OutcomeMasks;
            obj.taskState.trialNo = obj.Data.TrialNo;
            % execution functions 
            % trialTimeAverage(obj); % how to use class function?

            % processing 
            obj = convert0toNan(obj);
            forceFTconvert(obj);
            obj = forceHighSample(obj, obj.ft);
            obj = wamHighSample(obj, obj.wam);
            trials_all = setdiff(unique(TrialNo), 0);
            for trial_i = 1:length(trials_all)
                obj.trials(trial_i) = TrialScan(obj, trial_i);
                % align to mov
                obj.trials(trial_i) = alignMOV(obj.trials(trial_i));
            end
            %
            % plots
            axh = taskForceData(obj);
            axh = taskForceDatah(obj, axh);
            
            axh = taskEndpointPosition(obj);
            axh = taskEndpointPositionh(obj, axh);
            % taskEndpointPosition_relative(obj);
            % taskStateMuskFig(obj);
            % taskJointPosition_relateve(obj);
            

            
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
        %%% plot
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
            plot(obj.force_t, force');
            ylabel('force (N)');
            xlabel('time');
            legend('x', 'y', 'z'); % remember to alter the axis 
            title('Force data high sample');
        end
        function forceFTconvert(obj) % convert from select into world axis
            obj.force = obj.FTrot_M * obj.Data.Force.Sensor(1:3,:);
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
            plot(obj.wam_t, (position)');  
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
        function plotMeantrial(obj)
            % plot the meaned trial according to the task condition
        end
    end

end

