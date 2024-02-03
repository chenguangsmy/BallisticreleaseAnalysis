% 3. do impedance estimation on each trial, check K, B, M and fitting

%% laod data 
clear all; close all; clc;
% [13, 2], [15, 2], [17, 2], [20,2] [13, 4], [15, 4], [17, 4], [20,4]
subject_dex = [2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 20 21];  % subj11 has 5 sessions in 4 directions 
dir_dex     = [1 2 3 4];
session_dex= [  [4408 4409 4412 4414]; ...               % MR conducting
                [4418 4419 4422 4421]; ...               % NN conducting
                [4432 4434 4437 4438]; ...               % HM conducting
                [4379 4380 4382 4383]; ...               % CZ testing
                [4385 4386 4387 4388]; ...               % HA testing
                [4446 4448 4450 4451]; ...               % FM conducting
                [4455 4456 4458 4459]; ...               % FT direction +X, -X, +Y, -Y
                [4463 4464 4466 4467]; ...               % QX
                [4472 4473 4476 4477]; ...               % VC
                [4491 4492 4494 4495]; ...               % DS  
                [4500 4501 4503 4504]; ...               % AS
                [4512 4513 4515 4516]; ...               % XZ
                [4520 4521 4523 4524]; ...               % ZC
                [4530 4531 4533 4534]; ...               % KO
                [4542 4543 4545 4546]; ...               % SL
                [4558 4560 4562 4563]; ...               % AK
                [4573 4574 4576 4577]; ...                 % RL
                [4583 4584 4586 4587];
                ];

subject_arr = repmat(subject_dex, 4, 1);
subject_arr = subject_arr(:);

dir_arr = repmat(dir_dex,1,18);
dir_arr = dir_arr'; 


session_arr = session_dex; 
session_arr = session_dex(:);


% 1-6, 11-18 don't have estimations 
for session_i = [41:72] % wrong when session_i == 17
    fname = ['exampleData_subj' num2str(subject_arr(session_i)) 'dir' num2str(dir_arr(session_i)), '.mat']
    fdir  = 'test_data';
    load([fdir '/' fname], 'headers', 'recordedData');

    if (isfield(headers.trialHeader, 'estK'))
        continue;
    end

    % make plot 
    fh = figure('unit', 'inch', 'position', [0 0 4 6]);
    axh(1) = subplot(2,1,1);
    hold on;
    set(axh(1), 'linewidth', 1);
    set(axh(1), 'fontsize', 15);
    ylabel('x (m)');
    axh(2) = subplot(2,1,2);
    hold on;
    set(axh(2), 'linewidth', 1);
    set(axh(2), 'fontsize', 15);
    ylabel('F (N)');
    xlabel('t (s)');

    t_range = [-0.05 0.5];
    for trial_dex = 1:length(headers.trialHeader)

        t_orig = recordedData(trial_dex).data.t_shift;
        t_dex  = t_orig > t_range(1) & t_orig < t_range(2);
        t = t_orig(t_dex);
        x = recordedData(trial_dex).data.x(1,t_dex);
        x_hold = mean(x(t>t_range(1) & t<0));
        x = x - x_hold;

        F = recordedData(trial_dex).data.f(1,t_dex);
        F_hold = mean(F(t>t_range(1) & t<0));
        F_iv = F_hold - F;

        %     saveas(fh, ['trials_fig/trial' num2str(trial_dex) '.png'])
        data_est_UPs = iddata(x',F_iv',1/2000);
        sysUP_s = tfest(data_est_UPs,2,0);
        opt = predictOptions('InitialCondition','z');
        [yp,~,~] = predict(sysUP_s,data_est_UPs,0,opt);
        %         preddisp_interp_t(trial_i,:) = yp.OutputData;
        [NUM_UPs,DEN_UPs] = tfdata(sysUP_s);
        headers.trialHeader(trial_dex).estK = DEN_UPs{1}(3)/NUM_UPs{1}(3);
        headers.trialHeader(trial_dex).estB = DEN_UPs{1}(2)/NUM_UPs{1}(3);
        headers.trialHeader(trial_dex).estM = DEN_UPs{1}(1)/NUM_UPs{1}(3);
        headers.trialHeader(trial_dex).estFIT = sysUP_s.Report.Fit.FitPercent;

        cla(axh(1));
        cla(axh(2));

        plot(axh(1), t, x, 'b', 'linewidth', 1);
        plot(axh(1), t, yp.OutputData, '--', 'linewidth', 1)
        plot(axh(2), t, F_iv, 'b', 'linewidth', 1);

        if (headers.trialHeader(trial_dex).outcome == 1)
            plot(axh(1), t, x, 'b', 'linewidth', 1);
            plot(axh(1), t, yp.OutputData, '--', 'linewidth', 1)
            plot(axh(2), t, F_iv, 'b', 'linewidth', 1);
            sgtitle(['trial ' num2str(headers.trialHeader(trial_dex).tNo)]);
        else % fail
            plot(axh(1), t, x, 'r', 'linewidth', 1);
            plot(axh(1), t, yp.OutputData, '--', 'linewidth', 1)
            plot(axh(2), t, F_iv, 'r', 'linewidth', 1);
            sgtitle({ ['trial ' num2str(headers.trialHeader(trial_dex).tNo)]; headers.trialHeader(trial_dex).failReason});
        end
        
        title(axh(1), {['F ' num2str(headers.trialHeader(trial_dex).tarF)]; ['x ' num2str(headers.trialHeader(trial_dex).tarL)]});
        title(axh(2), {['K ' num2str(headers.trialHeader(trial_dex).estK)]; ['M' num2str(headers.trialHeader(trial_dex).estM)]; ['FIT' num2str(headers.trialHeader(trial_dex).estFIT)]});

        picname = ['dir' num2str(dir_arr(session_i)) ...
            'cond' num2str(headers.trialHeader(trial_dex).bNo) ... 
            'subj' num2str(dir_arr(session_i)) ...
            'trial' num2str(headers.trialHeader(trial_dex).tNo) ...
            '.png'];
        saveas(fh, ['trials_fig/' picname])
    end

    save([fdir '/' fname], 'headers', '-append', '-v7.3');

end