% export data as Federico's request 

%% 1. all the data with previous setting, set up as the conventional form
% but raw EMG
% sss1 = SessionsScan();

% check the data here...
% clear; 
% load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' ...
%     'processedData/ss4310_4356.mat']);

for subj_i = 1:4
    for dir_i = 1:4
        for fce_i = 1:3
            for disp_i = 1:3
                for trial_i = 1:9
                    dattmp = data{subj_i,dir_i,fce_i,disp_i,trial_i,1};
                    clf; clear axh;
                    axh(1) = subplot(3,1,1); 
                    plot(dattmp.t, dattmp.ox(1,:,1));
                    grid on;
                    axh(2) = subplot(3,1,2); 
                    plot(dattmp.t, dattmp.f(1,:));
                    grid on;
                    axh(3) = subplot(3,1,3); 
                    plot(dattmp.t, abs(dattmp.emg));
                    linkaxes(axh, 'x');
                    xlim([-0.5 0.5]);
                    grid on;
                end
            end
        end
    end
end
%% 2. all raw EMG in session-based

field_keep = {'t', 'x', 'v', 'ts', 'f', 'tNo', 'sNo', 'ox', 'emg'};
ss_list = { 4310 4313 4311 4314;                  %HA testing % MVF 4315 4312
            4325 4328 4326 4329;                  %NN testing % MVF 4324 4327
            4336 4337 [4339 4340]  4341;             %HM testing % MVF 4335 4338
            4351 [4352 4353] 4355 4356};

raw_data = cell(4,4); 
for subj_i = 1:4
    for dir_i = 1:4
        for session_i = 1:length(ss_list{subj_i,dir_i})
            clear data datatmp; 
            sstmp = SessionScan(ss_list{subj_i,dir_i}(session_i),1);
            for field_i = 1:length(field_keep)
                datatmp.(field_keep{field_i}) = sstmp.data.(field_keep{field_i});
            end
            data{session_i} = datatmp;
        end
        raw_data{subj_i,dir_i} = data;
    end
end
data_all = raw_data;
save(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' ...
    'processedData/ss4310_4356_raw.mat'], 'data_all', '-v7.3');

clear data raw_data data_all;
ss_mvf_list = {4315 4312;
               4324 4327;
               4335 4338;
               4349 4354};
mvf_data = cell(4,2);
for subj_i = 1:4
    for dir_i = 1:2
        clear data datatmp; 
        session_i = ss_mvf_list{subj_i,dir_i};
        sstmp = SessionScan(session_i);
        for field_i = 1:length(field_keep)
                datatmp.(field_keep{field_i}) = sstmp.data.(field_keep{field_i});
        end
        mvf_data{subj_i,dir_i} = datatmp;
    end
end
save(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' ...
    'processedData/ss4310_4356_raw.mat'], 'mvf_data', '-v7.3', '-append');

%% check the data here
% load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' ...
%     'processedData/ss4310_4356_raw.mat']);
for subj_i = 1:4
    for dir_i = 1:2
%         dattmp = data_all{subj_i,dir_i}{1};
        dattmp = mvf_data{subj_i,dir_i};
        clf; clear axh;
        axh(1) = subplot(3,1,1);
        plot(dattmp.t, dattmp.ox(1,:,1));
        grid on;
        axh(2) = subplot(3,1,2);
        plot(dattmp.t, dattmp.f(1,:));
        grid on;
        axh(3) = subplot(3,1,3);
%         plot(dattmp.t, abs(dattmp.emg));
        plot(dattmp.t, dattmp.emg);
        linkaxes(axh, 'x');
        grid on;
    end
end