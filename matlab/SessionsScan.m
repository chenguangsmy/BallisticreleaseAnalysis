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
        sessions
        subjprop % subject properties
        trials TrialScan
        export_cond
        filename
    end
    
    methods
        function obj = SessionsScan(ss_input, ifsave) %inputArg1,inputArg2
            %SESSIONSSCAN Construct an instance of this class
            %   read several sessions, and put all the trials together with
            %   conditions (the trial descriptions) aligned.
            %
            %
            % 1. Have the subjects - directions - sessions trials
            if (exist('ss_input', 'var'))
                ss_num = ss_input;
            else

                ss_num = {  ...                               % avoid initials in final code
                    % [4401 4402 4404 4405 4426 4427] ...     % BH conducting
                    % [4408 4409 4412 4414] ...               % MR conducting
                    % [4418 4419 4422 4421] ...               % NN conducting
                    % [4432 4434 4437 4438] ...               % HM conducting
                    % [4379 4380 4382 4383] ...               % CZ testing
                    % [4385 4386 4387 4388] ...               % HA testing
                    % [4446 4448 4450 4451]  ...              % FM conducting
                    % [4455 4456 4458 4459]  ...              % FT conduting
                    % [4463 4464 4466 4467]  ...              % QX conducting
                    % [4472 4473 4476 4477] ...               % VC
                    % [4481 4482 4483 4485 4486] ...          % DS
                    % [4491 4492 4494 4495] ...               % BW
                    % [4500 4501 4503 4504] ...               % AS
                    % [4512 4513 4515 4516] ...               % XZ
                    % [4520 4521 4523 4524] ...               % ZC
                    % [4530 4531 4533 4534] ...               % KO
                    % [4542 4543 4545 4546] ...               % SL
                    % [4558 4560 4562 4563] ...               % AK
                    % [4568 4569] ... % unfinished            % AR
                    [4573 4574 4576 4577] ...                 % RL
                    [4583 4584 4586 4587] ...                 % HD
                    };
            end

            if (~exist('ifsave', 'var'))
                ifsave = 0;
            end

            allsubj.height = [170.2000  175.0000  182.9000  160.0000  183.0000  172.0000  180.0000  183.0000  182.0000  163.0000  180.0000  182.0000  188.0000  163.0000  180.3000  162.6000  160.0000  157.5000    168.0000  162.0000];
            allsubj.weight = [63.2000   65.7000   85.8000   50.2000   82.0000   77.8000   90.0000   78.5000   81.3000   60.0000   75.0000  113.0000  101.5000   55.2000   65.0000   62.1000    58.0000   58.5000     57.3000   46.4000];
            allsubj.gender = 'MMMMMMMMMFMMMFMFFFFF';
            subj_sel = [19 20];

            obj.subjprop.height = allsubj.height(subj_sel);
            obj.subjprop.weight = allsubj.height(subj_sel);
            obj.export_cond.subject = 1:size(ss_num,2);
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
                    disp(['Loading Sessions ' num2str(sessions_load_i) '/' ...
                        num2str(ss_num_total)]);
                    % ss_tmp{subj_i}{ss_i} = SessionScan(ss_num{subj_i}(ss_i),1); % if trying re-read
                    ss_tmp{subj_i}{ss_i} = SessionScan(ss_num{subj_i}(ss_i));
                    obj.sessions{subj_i}{ss_i} = ss_tmp{subj_i}{ss_i};
                end
            end

            % 3. Concatenate trials and with their subject/direction
            % information attached

            obj.cond.subject    = [];
            obj.cond.direction  = [];
            obj.cond.fce        = [];
            obj.cond.disp       = [];
            obj.cond.tNo        = [];
            obj.cond.ssNo       = [];
            obj.cond.sf         = [];
            obj.cond.pert       = [];

            for subj_i = 1:length(ss_num)
                for ss_i            = 1:length(ss_num{subj_i})
                    obj.trials      = [obj.trials ss_tmp{subj_i}{ss_i}.trials];
                    ss_cond.subject = subj_i*ones(size(ss_tmp{subj_i}{ss_i}.trials));
                    ss_cond.direct  = [ss_tmp{subj_i}{ss_i}.trials.tarR];
                    ss_cond.fce     = [ss_tmp{subj_i}{ss_i}.trials.tarF]; 
                    ss_cond.disp    = [ss_tmp{subj_i}{ss_i}.trials.tarL];
                    ss_cond.tNo     = [ss_tmp{subj_i}{ss_i}.trials.tNo];
                    ss_cond.ssNo    = [ss_num{subj_i}(ss_i)]*ones(size(ss_tmp{subj_i}{ss_i}.trials));
                    ss_cond.sf      = [ss_tmp{subj_i}{ss_i}.trials.outcome];    % relying on online judge
                    ss_tmp{subj_i}{ss_i}.trials(2).outcome = 0;                 % manual fail the trial avoid the nan ox data
                    
                    
                    ss_cond.pert = [ss_tmp{subj_i}{ss_i}.trials.ifpert]; % perterturbation type... edt

                    % concatinate into conditions
                    obj.cond.subject    = [obj.cond.subject     ss_cond.subject];
                    obj.cond.direction  = [obj.cond.direction   ss_cond.direct];
                    obj.cond.fce        = [obj.cond.fce         ss_cond.fce];
                    obj.cond.disp       = [obj.cond.disp        ss_cond.disp];
                    obj.cond.tNo        = [obj.cond.tNo         ss_cond.tNo];
                    obj.cond.ssNo       = [obj.cond.ssNo        ss_cond.ssNo];
                    obj.cond.sf         = [obj.cond.sf          ss_cond.sf];
                    obj.cond.pert       = [obj.cond.pert        ss_cond.pert];
                end
            end

            session_min = min(obj.cond.ssNo);
            session_max = max(obj.cond.ssNo);
            obj.filename = sprintf('ss%4d_%4d', session_min, session_max);

            % set every first trial in the session 0

            obj = obj.SessionsSpecifyRotation();

            if (ifsave)
                obj.SessionsExport();
            end
        end
        
        function [performc] = getTrialPerformance(obj) 
            % tobe re-written for the multi sessions processing
            performc = cell(max(obj.export_cond.subject), 1);
            for subj_i = 1:max(obj.export_cond.subject)
                performc{subj_i}.time_all   = zeros(1,1);
                performc{subj_i}.time_eachd = zeros(1,4);
                performc{subj_i}.rate_all   = zeros(1,1);
                performc{subj_i}.rate_eachd = zeros(1,4);
                performc{subj_i}.rate_eachc = zeros(1,4,3,3);
                performc{subj_i}.trial_all  = zeros(1,1);
                
                % time-all
                trial_idx = obj.cond.subject == subj_i;
                performc{subj_i}.time_all = max([obj.trials(trial_idx).edn_t]) - ...
                                            min([obj.trials(trial_idx).bgn_t]);
                trialLength = [obj.trials.edn_t] - [obj.trials.bgn_t];
                
                % rate-all
                performc{subj_i}.rate_all = sum([obj.cond.sf(trial_idx)] & ~isnan([obj.cond.direction(trial_idx)]))/sum(~isnan(obj.cond.direction(trial_idx)));
                performc{subj_i}.trial_all = sum(~isnan(obj.cond.direction(trial_idx))); % trial started (excluded from the rest trials

                for d_i = 1:length(obj.export_cond.direction)
                    trial_idx = obj.cond.subject == subj_i & ...
                                obj.cond.direction == obj.export_cond.direction(d_i);
                    performc{subj_i}.time_eachd(d_i) = sum(trialLength(trial_idx));
                    performc{subj_i}.rate_eachd(d_i) = sum([obj.cond.sf(trial_idx)]) / ...
                        sum(trial_idx);

                    for f_i = 1:length(obj.export_cond.fce)
                        for t_i = 1:length(obj.export_cond.disp)
                            trial_idx = obj.cond.subject == subj_i & ...
                                      obj.cond.direction == obj.export_cond.direction(d_i) & ...
                                            obj.cond.fce == obj.export_cond.fce(f_i) & ...
                                           obj.cond.disp == obj.export_cond.disp(t_i);  

                            performc{subj_i}.rate_eachc(1,d_i,f_i,t_i) = ...
                                sum([obj.cond.sf(trial_idx)])/sum(trial_idx);
                        end
                    end
                end
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to export data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function data = SessionsExport(obj)
            % SessionsExport Export to formatted data matrix
            %   Formatted data matrix is a 6-D matrix: 
            %   subj - direction - force - distance - trials - pert
            % 
            % In experiment with one direction multiple pulse, perturbation
            % could be 4: 
            %       [nopert, pertPreMotion, pertInMotion, pertPostMotion]
            % 
            % In experiment with 4 directions no pulse, perturbation is 1
            data = cell(...
                        length(unique(obj.export_cond.subject)), ... % subjects
                        4, ...                                       % directions
                        3, ...                                       % fce
                        3, ...                                       % dist
                        9, ...                                       % trials
                        4);                                          % pert
            data_index_ss = zeros(size(data));
            data_index_tr = zeros(size(data));
            % one direction
            pert_export_code = [0 1 6]; % each pulse
            trials_req =       [9 nan nan];      % each perturb
            for subj_i = 1:length(obj.export_cond.subject)
                for direction_i = 1:4
                    for fce_i = 1:3
                        for disp_i = 1:3
                            for pert_i = 1:3
                                trialMask = [obj.cond.subject   == subj_i &...
                                             obj.cond.direction == obj.export_cond.direction(direction_i) & ...
                                             obj.cond.fce       == obj.export_cond.fce(fce_i) & ...
                                             obj.cond.disp      == obj.export_cond.disp(disp_i) & ... 
                                             obj.cond.pert      == pert_export_code(pert_i) & ...
                                             obj.cond.sf        == 1
                                             ];


                                if (sum(trialMask)) == 0 % no trials
                                    continue;
                                end
                                
                                trial_idx = find(trialMask);
                                % if trial is enough, get trials, 
                                if (sum(trialMask))>=trials_req(pert_i)
                                    trial_idx = find(trialMask);
                                    trial_idx = trial_idx(end-trials_req(pert_i)+1:end);
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
                                    data_index_ss(subj_i,direction_i,fce_i,disp_i,trial_idx_dest,pert_i) = obj.trials(trial_idx_from).ssnum;
                                    data_index_tr(subj_i,direction_i,fce_i,disp_i,trial_idx_dest,pert_i) = obj.trials(trial_idx_from).tNo;
                                end
                            end
                        end
                    end
                end
            end
            subjprop = obj.subjprop;
            save(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/' obj.filename '.mat'], 'data', 'data_index_ss', 'data_index_tr', '-v7.3');
            save(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/' obj.filename '.mat'], 'subjprop', '-append');
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
            TRIALS_REQ = 9;

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
            % Mark the session which subject are doing front/back movements
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
                            4328 4329 ...
                            4339 4340 4341 ...
                            4354 4355 4356 ...
                            4382 4383 ...
                            4387 4388 ...
                            4404 4405 4427 ...
                            4412 4414 ...
                            4421 4422 ...
                            4437 4438 ...
                            4450 4451 ...
                            4458 4459 ...
                            4466 4467 ...
                            4476 4477 ...
                            4485 4486 ...
                            4494 4495 ...
                            4503 4504 ...
                            4515 4516 ...
                            4523 4524 ...
                            4533 4534 ...
                            4545 4546 ...
                            4562 4563 ...
                            4576 4577 ...
                            4586 4587 ...
                            ]; 
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

