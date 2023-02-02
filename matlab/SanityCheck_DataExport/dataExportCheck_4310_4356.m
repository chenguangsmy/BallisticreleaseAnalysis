% in data ss4310_4356, Federico saids some trial has no EMG data. 
% let me see

blp = ballisticReleaseTaksPlots

% run this cause error because of :'lack of valid data'
blp.plotEMG_release_1 

% Check them by each of the trials 
obj = blp; 
pert_i = 1
for subj_i = 4
    for dir_i = 1
        for fce_i = 1
            for dist_i = 1
                figure(); hold on;
                for trial_i = 1:9
                    
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,:,1))
plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,:,1))
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(1:8,:))
                end
            end
        end
    end
end


%% The file should be: 
data_dir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
data_name = 'ss4310_4356_emg';
load([data_dir '/' data_name] );

obj.data = data; 
col_type = colormap('lines');
%%
close all;
pert_i = 1
for subj_i = 3%1:4
    figure(); hold on;
    for dir_i = 2%1:4
        for fce_i = 1
            for dist_i = 1
                
                for trial_i = 1:9
                    
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,:,1))
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ox(1,:,1), 'color', col_type(dir_i,:))
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,:), 'color', col_type(dir_i,:))
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.f(1,:), 'color', col_type(dir_i,:))
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(:,:));
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.x(1,:,1))
% plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(1:8,:))
plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.ts, 'color', col_type(dir_i,:))
%                     for ch_i = 1:8
%                         plot(obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.t,ch_i + obj.data{subj_i,dir_i,fce_i,dist_i,trial_i,pert_i}.emg(ch_i,:));
%                     end
                end
            end
        end
    end
    grid on;
    title(['subject', num2str(subj_i)]);
end
% conclusion, subject3 direction 2 has problem with shift, also has no ts
% here 
% problem corresponding session: 4434?
%[4432 4434 4437 4438] ...              %HM conducting

% subject 2 data are small. 

%% Do the specific session- 4434
sstmp = SessionScan(4438)
sstmp.plotTrialfyPositionh