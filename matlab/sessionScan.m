classdef sessionScan < handle
    %VARSCAN scanning some variables in the formatted data
    %   Detailed explanation goes here

    properties
        Data
        trials_num
        duration
        duration_avg
        hand_pos        % position read from WAM endpoint
        hand_pos_offset % the center_pos for WAM endpoint
        force           % force in the force transducer
        FTrot_M = ...   % global: x-right, y-front, z-up, FT_base x-backup, y-frontup, z-left
            [0          0           cosd(180)
            cosd(75)    -sind(75)    0
            sind(75)    cosd(75)    0];
    end
    
    methods
        
        function obj = sessionScan(ss_num) %(inputArg1,inputArg2)
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
            
            Time = Data.Time;
            TrialNo = Data.TrialNo;
            obj.duration = max(Time);
            obj.trials_num = max(TrialNo);
            obj.Data = Data;
            obj.hand_pos_offset = Data.Position.Center(:,~isnan(Data.Position.Center(1,:)));
            obj.hand_pos_offset = obj.hand_pos_offset(:,1); 
            % execution functions 
            % trialTimeAverage(obj); % how to use class function?

            % processing 
            forceFTconvert(obj);
            % plots
            taskForceData(obj);
            %taskEndpointPosition(obj);
            %taskStateMuskFig(obj);
            taskJointPosition(obj);
        end
        
        function trialTimeAverage(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.duration_avg = obj.duration/double(obj.trials_num);
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
            plot(obj.Data.Position.JointPosition');
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
        function taskJointPosition_relateve(obj) % Plot relative position. 
            % plot 
            figure();
            axish(1) = subplot(2,1,1);
            % plot([diff(obj.Data.Position.JointPosition,1,2)-obj.Data.Position.JointVelocity(:,2:end)]');
            position_offset = obj.Data.Position.JointPosition(:,~isnan(obj.Data.Position.JointPosition(1,:)));
            position_offset = repmat(position_offset(:,1),1,size(obj.Data.Position.JointPosition,2));
            plot((obj.Data.Position.JointPosition - position_offset)');
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
        function taskForceData(obj)
            figure();
            % force = obj.Data.Force.Sensor(1:3,:); 
            force = obj.force;
            plot(force');
            ylabel('force (N)');
            xlabel('time points');
            legend('x', 'y', 'z'); % remember to alter the axis 
            title('Force data before convert axis');
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
            plot((position - position_offset)');  
            ylabel('relative endpoint positions');
            xlabel('time points');
            legend('x', 'y', 'z');
            title('relative endpoint positions');
        end
        function taskEndpointPosition(obj)
            figure();
            position = obj.Data.Position.Actual'; 
            plot((position)');  
            ylabel('endpoint positions');
            xlabel('time points');
            legend('x', 'y', 'z');
            title('relative endpoint positions');
        end
    end

end

