%% A list of data that worth for behavioral analysis 
clear; close all; clc
% 
subject_dex = [2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 20 21]';  % subj11 has 5 sessions in 4 directions 
session_dex= [  [4408 4409 4412 4414]; ...               % MR conducting
                [4418 4419 4422 4421]; ...               % NN conducting
                [4432 4434 4437 4438]; ...               % HM conducting
                [4379 4380 4382 4383]; ...               % CZ testing
                [4385 4386 4387 4388]; ...               % HA testing
                [4446 4448 4450 4451]; ...               % FM conducting
                [4455 4456 4458 4459]; ...               % direction +X, -X, +Y, -Y
                [4463 4464 4466 4467]; ...
                [4472 4473 4476 4477]; ...
                [4491 4492 4494 4495]; ... 
                [4500 4501 4503 4504]; ...               % AS
                [4512 4513 4515 4516]; ...               % XZ
                [4520 4521 4523 4524]; ...               % ZC
                [4530 4531 4533 4534]; ...               % KO
                [4542 4543 4545 4546]; ...               % SL
                [4558 4560 4562 4563]; ...               % AK
                [4573 4574 4576 4577]; ...                 % RL
                [4583 4584 4586 4587];
                ];

% seems got error in session 4563! and later!

%% grab the data from a session that has subject 8 direction 1 and all the
% trials (sucessful and un-sucessful) 

for subj_i = [18]
    for dir_i = 2:4


        session_i = session_dex(subj_i, dir_i);
        ss = SessionScan(session_i, 1);
        
        % save to...
        headers.trialHeader = [];
        recordedData = [];
        trial_num = length(ss.trials);

        for trial_dex = 2:trial_num
            headers.trialHeader(trial_dex).tNo = ss.trials(trial_dex).tNo;
            %     headers.trialHeader(trial_dex).bNo = ss.trials(trial_dex).bNo;
            headers.trialHeader(trial_dex).ssnum = ss.trials(trial_dex).ssnum;
            %     headers.trialHeader(trial_dex).bgn = ss.trials(trial_dex).bgn;
            %     headers.trialHeader(trial_dex).edn = ss.trials(trial_dex).edn;
            headers.trialHeader(trial_dex).bgn_t = ss.trials(trial_dex).bgn_t;
            headers.trialHeader(trial_dex).edn_t = ss.trials(trial_dex).edn_t;
            headers.trialHeader(trial_dex).outcome = ss.trials(trial_dex).outcome;
            %     headers.trialHeader(trial_dex).outcomeo = ss.trials(trial_dex).outcomeo;
            headers.trialHeader(trial_dex).states = ss.trials(trial_dex).states;
            %     headers.trialHeader(trial_dex).states_arr = ss.trials(trial_dex).states_arr;
            %     headers.trialHeader(trial_dex).tarR = ss.trials(trial_dex).tarR;
            headers.trialHeader(trial_dex).tarL = ss.trials(trial_dex).tarL;
            %     headers.trialHeader(trial_dex).tarP = ss.trials(trial_dex).tarP;
            headers.trialHeader(trial_dex).tarF = ss.trials(trial_dex).tarF;
            %     headers.trialHeader(trial_dex).idx_bgn = ss.trials(trial_dex).idx_bgn;
            %     headers.trialHeader(trial_dex).idx_prt = ss.trials(trial_dex).idx_prt;
            %     headers.trialHeader(trial_dex).tarHD = ss.trials(trial_dex).tarHD;
            %     headers.trialHeader(trial_dex).idx_fcr = ss.trials(trial_dex).idx_fcr;
            %     headers.trialHeader(trial_dex).idx_mov = ss.trials(trial_dex).idx_mov;
            %     headers.trialHeader(trial_dex).idx_hld = ss.trials(trial_dex).idx_hld;
            %     headers.trialHeader(trial_dex).idx_end = ss.trials(trial_dex).idx_end;
            %     headers.trialHeader(trial_dex).idx_rst = ss.trials(trial_dex).idx_rst;
            %     headers.trialHeader(trial_dex).time_orn = ss.trials(trial_dex).time_orn;
            %     headers.trialHeader(trial_dex).time = ss.trials(trial_dex).time;
            %     headers.trialHeader(trial_dex).xyi = ss.trials(trial_dex).xyi;
            %     headers.trialHeader(trial_dex).xyn = ss.trials(trial_dex).xyn;
            %     headers.trialHeader(trial_dex).opt_v = ss.trials(trial_dex).opt_v;
            %     headers.trialHeader(trial_dex).opt = ss.trials(trial_dex).opt;
            %     headers.trialHeader(trial_dex).opth = ss.trials(trial_dex).opth;
            %     headers.trialHeader(trial_dex).position_offset = ss.trials(trial_dex).position_offset;
            recordedData(trial_dex).data = ss.trials(trial_dex).data;
            %     headers.trialHeader(trial_dex).ifpert = ss.trials(trial_dex).ifpert;
            %     headers.trialHeader(trial_dex).pert_f = ss.trials(trial_dex).pert_f;
            %     headers.trialHeader(trial_dex).wamKp = ss.trials(trial_dex).wamKp;
            %     headers.trialHeader(trial_dex).wamBp = ss.trials(trial_dex).wamBp;
            %     headers.trialHeader(trial_dex).pred_x0 = ss.trials(trial_dex).pred_x0;
            %     headers.trialHeader(trial_dex).pred_K = ss.trials(trial_dex).pred_K;
            %     headers.trialHeader(trial_dex).pred_D = ss.trials(trial_dex).pred_D;
            %     headers.trialHeader(trial_dex).pred_A = ss.trials(trial_dex).pred_A;
            %     headers.trialHeader(trial_dex).pred_J = ss.trials(trial_dex).pred_J;
            %     headers.trialHeader(trial_dex).pred_S = ss.trials(trial_dex).pred_S;
            %     headers.trialHeader(trial_dex).pert_iter = ss.trials(trial_dex).pert_iter;
            %     headers.trialHeader(trial_dex).perturbation_length = ss.trials(trial_dex).perturbation_length;
        end

        headers.trialHeader = headers.trialHeader(2:end);
        recordedData = recordedData(2:end); 


        % back to block num here
        switch (dir_i) 
            case 1
                F_tars = [20        25      20      25      15      25      15      20      15       20];
                x_tars = [0.05  0.075    0.025   0.025      0.05    0.05    0.025   0.075   0.075   0.05];
            case 2 %...
                F_tars = [20        15      25      15      20      20      25      25      15       20];
                x_tars = [0.05  0.025    0.025   0.075      0.025   0.075   0.05   0.075    0.05   0.05];
            case 3%...
                F_tars = [20        25      20      25      15      25      15      20      15       20];
                x_tars = [0.05  0.075    0.025   0.025      0.05    0.05    0.025   0.075   0.075   0.05]

            case 4 %...
                F_tars = [20        15      25      15      20      20      25      25      15       20];
                x_tars = [0.05  0.025    0.025   0.075      0.025   0.075   0.05   0.075    0.05   0.05];

        end
        trial_dex = 1;
        for b_i = 1:10
            while (1)
                if headers.trialHeader(trial_dex).tarF == F_tars(b_i) && headers.trialHeader(trial_dex).tarL == x_tars(b_i)
                    headers.trialHeader(trial_dex).bNo = b_i;
                    trial_dex = trial_dex + 1;
                else
                    break
                end

                if trial_dex > length(headers.trialHeader)
                    break
                end
            end
        end

        % remove trials that don't release at all;
        trial_num = length(headers.trialHeader);
        trials_sel = ones(1,trial_num);
        for trial_dex = 1:trial_num
            headers.trialHeader(trial_dex).failReason = '';
            if isempty(intersect(headers.trialHeader(trial_dex).states, [5 6]))
                trials_sel(trial_dex) = 0;
                headers.trialHeader(trial_dex).failReason = 'FH';
            end

            if isempty(intersect(headers.trialHeader(trial_dex).states, 6))
                headers.trialHeader(trial_dex).failReason = 'TO';
            end

            if isempty(setdiff([1:8 99], headers.trialHeader(trial_dex).states))
                headers.trialHeader(trial_dex).failReason = 'CR';
            end
        end

        headers.trialHeader = headers.trialHeader(logical(trials_sel));
        recordedData        = recordedData(logical(trials_sel));

        fdir = 'test_data';
        fname = ['exampleData_subj' num2str(subject_dex(subj_i)) 'dir' num2str(dir_i) '.mat'];
        save([fdir '/' fname], 'headers', 'recordedData', '-v7.3');

    end
end
