classdef SessionsScan
    %SESSIONSSCAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cond
        % cond.subject
        % cond.direction
        % cond.fce
        % cond.disp
        % cond.tNo
        % cond.ssNo
        % cond.sf
        trials TrialScan
        export_cond
        filename
    end
    
    methods
        function obj = SessionsScan() %inputArg1,inputArg2
            %SESSIONSSCAN Construct an instance of this class
            %   read several sessions, and put all the trials together with
            %   conditions (the trial descriptions) aligned. 
            %
            % 
            % 1. Have the subjects - directions - sessions trials
            % subject 1, CZ
%             ss_num{1} = [4216 4217 4218 4219 4220]; % CZ            
%             ss_num{1} = [ 4217];
%             ss_num{1} = [4300 4301 4303 4304]; %CZ testing
%             ss_num{2} = [4310 4311 4313 4314]; % HA testing, 4315 & 4312 as MVF
            ss_num = {...[4300 4301 4303 4304] ...  %CZ testing
                      [4310 4311 4313 4314]...      %HA testing
                      [4324 4325 4327 4328]...      %NN testing 
                      [4336 4337 4339 4340 4341]};  %HM testing 
                      
            % subject 2, HA
%             ss_num{2} = [4204 4205 4206 4208]
%             ss_num{2} = [4222 4223 4224 4225 4226];   % HA 
%             ss_num{2} = [4235 4236 4237 4238 4239];   % HA Update condition
            % subject 3, CZ2
%             ss_num{3} = [4204 4205 4206 4208] % 4204:4211

            % the 1 direction with 2 pulses (pulse before release + pulse
            % at position hold 
%             ss_num{1} = [4253 4256 4259 4262 4257 4260 4263 4258 4261];
%             ss_num{2} = [4274 4265 4266 4267 4269 4270 4271 4272 4273];
%             ss_num{1} = [4274 4265 4266 4267 4269 4270 4271 4272 4273];
%             ss_num{1} = [4274 4265 4266 4263 4258];
            % export setup
            obj.export_cond.subject = [1:size(ss_num)];
            obj.export_cond.direction = [0 2 4 6];
            obj.export_cond.fce = [15 20 25];
            obj.export_cond.disp = [0.025 0.05 0.075];
            obj.export_cond.pert = 1;
            % 2. Read the sessions 
            % count how many sessions in total 
            ss_num_total = 0;
            for subj_i = 1:length(ss_num)
                ss_num_total = ss_num_total + length(ss_num{subj_i});
            end

            % loading one by one 
            sessions_load_i = 0;
            for subj_i = 1:length(ss_num) 
                for ss_i = 1:length(ss_num{subj_i})
                    sessions_load_i = sessions_load_i + 1;
                    % display
                    disp(['Loading Sessions ' num2str(sessions_load_i) '/' ...
                        num2str(ss_num_total)]);
                    % loading
                    ss_tmp{subj_i}{ss_i} = SessionScan(ss_num{subj_i}(ss_i));
                end
            end
            
            % 3. Concatenate trials and with their subject/direction
            % information attached

            obj.cond.subject = [];
            obj.cond.direction = [];
            obj.cond.fce = [];
            obj.cond.disp = [];
            obj.cond.tNo = [];
            obj.cond.ssNo = [];
            obj.cond.sf = [];
            obj.cond.pert = [];
            
            for subj_i = 1:length(ss_num)
                for ss_i = 1:length(ss_num{subj_i})
                    obj.trials = [obj.trials ss_tmp{subj_i}{ss_i}.trials];
                    ss_cond.subject = [subj_i]*ones(size(ss_tmp{subj_i}{ss_i}.trials));
                    ss_cond.direction = [ss_tmp{subj_i}{ss_i}.trials.tarR];
                    ss_cond.fce = [ss_tmp{subj_i}{ss_i}.trials.tarF]; %? why 2 values?
                    ss_cond.disp = [ss_tmp{subj_i}{ss_i}.trials.tarL];
                    ss_cond.tNo = [ss_tmp{subj_i}{ss_i}.trials.tNo];
                    ss_cond.ssNo = [ss_num{subj_i}(ss_i)]*ones(size(ss_tmp{subj_i}{ss_i}.trials));
                    ss_tmp{subj_i}{ss_i}.trials(2).outcome = 0;             % manual fail the trial avoid the nan ox data
                    ss_cond.sf = [ss_tmp{subj_i}{ss_i}.trials.outcome];   % relying on online judge
%                     ss_cond.sf = [[ss_tmp{subj_i}{ss_i}.trials.outcomeo] | [ss_tmp{subj_i}{ss_i}.trials.outcome]];   % relying on offline | online judge
%                     ss_cond.pert = [1]*ones(size(ss_tmp{subj_i}{ss_i}.trials)); % perterturbation type... edt
                    ss_cond.pert = [ss_tmp{subj_i}{ss_i}.trials.ifpert]; % perterturbation type... edt
                    
                    % concatinate into conditions 
                    obj.cond.subject = [obj.cond.subject ss_cond.subject];
                    obj.cond.direction = [obj.cond.direction ss_cond.direction];
                    obj.cond.fce = [obj.cond.fce ss_cond.fce];
                    obj.cond.disp = [obj.cond.disp ss_cond.disp];
                    obj.cond.tNo = [obj.cond.tNo ss_cond.tNo];
                    obj.cond.ssNo = [obj.cond.ssNo ss_cond.ssNo];
                    obj.cond.sf = [obj.cond.sf ss_cond.sf];
                    obj.cond.pert = [obj.cond.pert ss_cond.pert];
                end
            end

            session_min = min(obj.cond.ssNo);
            session_max = max(obj.cond.ssNo);
            obj.filename = sprintf('ss%4d_%4d', session_min, session_max);

            % set every first trial in the session 0 


            obj = obj.SessionsSpecifyRotation();

            obj.SessionsExport();
        end
        
        function data = SessionsExport(obj)
            %SessionsExport Export to formatted data matrix
            %   the formatted data matrix is a 6-D matrix, which is: 
            %   subj - direction - force - distance - trials - pert
            % 
            % In experiment with one direction multiple pulse, perturbation
            % could be 4: 
            %       [nopert, pertPreMotion, pertInMotion, pertPostMotion]
            % 
            % In experiment with 4 directions no pulse, perturbation is 1
            data = cell(length(unique(obj.export_cond.subject)), ... % subjects
                4, ...                                       % directions
                3, ...                                       % fce
                3, ...                                       % dist
                9, ...15, ...                                      % trials
                4);                                          % pert
%                 1);                                          % pert
            % one direction
            pert_export_code = [0 1 6]; % each pulse
%             trials_req =       [20 10 10];      % each perturb
             trials_req =       [9 nan nan];      % each perturb
            % four directions
%             pert_export_code = [0 nan nan];   
%             trials_req =       [15 nan nan];
            for subj_i = 1:length(obj.export_cond.subject)
                for direction_i = 1:4
                    for fce_i = 1:3
                        for disp_i = 1:3
                            for pert_i = 1:3
                                trialMask = [obj.cond.subject == subj_i &...
                                             obj.cond.direction == obj.export_cond.direction(direction_i) & ...
                                             obj.cond.fce == obj.export_cond.fce(fce_i) & ...
                                             obj.cond.disp == obj.export_cond.disp(disp_i) & ... 
                                             obj.cond.pert == pert_export_code(pert_i) & ...
                                             obj.cond.sf == 1
                                             ];


                                if (sum(trialMask)) == 0 % no trials
                                    continue;
                                end
                                
                                trial_idx = find(trialMask);
                                % if trial is enough, get trials, 
                                if (sum(trialMask))>=trials_req(pert_i)
                                    trial_idx = find(trialMask);
                                    trial_idx = trial_idx(1:trials_req(pert_i));
                                % if trial is not enough, get more trials
                                % from repeating
                                else
                                    trials_qualify_num = sum(trialMask);
                                    trials_lack = trials_req(pert_i) - trials_qualify_num;

                                    for trials_lacki = 1:trials_lack
                                        trial_idx(trials_qualify_num+trials_lacki) = ...
                                            trial_idx(trials_lacki);
                                        disp(['put trial' num2str(trials_lacki) 'in slot' num2str(trials_qualify_num+trials_lacki)]);
                                    end
                                end

                                for trial_idx_dest = 1:length(trial_idx)
                                    trial_idx_from = trial_idx(trial_idx_dest);
                                    data{subj_i,direction_i,fce_i,disp_i,trial_idx_dest,pert_i} = ...
                                        obj.trials(trial_idx_from).export_as_formatted(); % need edition.
%                                         obj.trials(trial_idx_from).export_as_formatted(1); % need edition.
                                end
                            end
                        end
                    end
                end
            end
            
            save(['data/processedData/' obj.filename '.mat'], 'data', '-v7.3');
        end
        
        function data = SessionsExportf(obj)
            %SessionsExport Export to formatted data matrix, only failed
            %trials
            %   the formatted data matrix is a 6-D matrix, which is: 
            %   subj - direction - force - distance - trials - pert
            % 
            % In experiment with one direction multiple pulse, perturbation
            % could be 4: 
            %       [nopert, pertPreMotion, pertInMotion, pertPostMotion]
            % 
            % In experiment with 4 directions no pulse, perturbation is 1
            data = cell(length(unique(obj.export_cond.subject)), ... % subjects
                4, ...                                       % directions
                3, ...                                       % fce
                3, ...                                       % dist
                9, ...                                      % trial
                4);                                          % pert
%               1);                                          % pert
            TRIALS_REQ = 9; % 15

            for subj_i = 1:length(obj.export_cond.subject)
                for direction_i = 1:4
                    for fce_i = 1:3
                        for disp_i = 1:3
                            for pert_i = 1
                                trialMask = [obj.cond.subject == subj_i &...
                                             obj.cond.direction == obj.export_cond.direction(direction_i) & ...
                                             obj.cond.fce == obj.export_cond.fce(fce_i) & ...
                                             obj.cond.disp == obj.export_cond.disp(disp_i) & ... 
                                             obj.cond.pert == pert_i & ...
                                             obj.cond.sf == 0
                                             ];


                                if (sum(trialMask)) == 0 % no trials
                                    continue;
                                end
                                
                                trial_idx = find(trialMask);
                                % if trial is enough, get trials, 
                                if (sum(trialMask))>=TRIALS_REQ
                                    trial_idx = find(trialMask);
                                    trial_idx = trial_idx(1:TRIALS_REQ);
                                % if trial is not enough, get more trials
                                % from repeating
                                else
                                    trials_qualify_num = sum(trialMask);
                                    trials_lack = TRIALS_REQ - trials_qualify_num;

                                    for trials_lacki = 1:trials_lack
                                        trial_idx(trials_qualify_num+trials_lacki) = ...
                                            trial_idx(trials_lacki);
                                        disp(['put trial' num2str(trials_lacki) 'in slot' num2str(trials_qualify_num+trials_lacki)]);
                                    end
                                end

                                for trial_idx_dest = 1:length(trial_idx)
                                    trial_idx_from = trial_idx(trial_idx_dest);
                                    data{subj_i,direction_i,fce_i,disp_i,trial_idx_dest,pert_i} = ...
                                        obj.trials(trial_idx_from).export_as_formatted(); % need edition.
%                                         obj.trials(trial_idx_from).export_as_formatted(1); % need edition.
                                end
                            end
                        end
                    end
                end
            end
            
            save(['data/processedData/' obj.filename 'f.mat'], 'data', '-v7.3');
        end
        

        function obj = SessionsSpecifyRotation(obj)
            % Specify the roated sessions in this function
            % Where the 'target rotation' array will be changed. 
            % 
            % This function explains which sessions the subect is sitting
            % on the left of the robot, doing subjects' front and back
            % using robot's right and front. 
            % It will change all the target rotation in this session (0, 4)
            % to (2, 6) 

            sessions_rot = [4218 4219 4220 ...
                            4225 4226 ...
                            4237 4238 4239 ...
                            4303 4304 ... 
                            4313 4314 ...
                            4327 4328 ...
                            4339 4340 4341]; 
                % first assume these sessions. These can be read from .conf
                % in the future. 
                
            % look for sessions that intersect
            sessions_in = unique(obj.cond.ssNo);
            sessions_int= intersect(sessions_rot, sessions_in);
            for ss_i = sessions_int
                % do the 0 --> 2
                idx = obj.cond.ssNo == ss_i & obj.cond.direction == 0;
                obj.cond.direction(idx) = 2;

                % do the 4 --> 6
                idx = obj.cond.ssNo == ss_i & obj.cond.direction == 4;
                obj.cond.direction(idx) = 6;
            end
        end
    end
end

