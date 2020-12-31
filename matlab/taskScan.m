% Scan specific data
file_dir = '/Users/cleave/Documents/projPitt/Ballistic_release_data/Formatted';
fname = 'KingKong.01545.mat';
dst_file = [file_dir '/' fname(1:end-4) '_tf.mat'];
load([file_dir '/' fname]);

rot_M = ... % global: x-right, y-front, z-up, FT_base x-backup, y-frontup, z-left
    [0          0           cosd(180)
    cosd(135)   cosd(45)    0
    cosd(45)    cosd(45)    0];

Data.Force.Sensor(7:9,:) = rot_M*Data.Force.Sensor(1:3,:);
save([file_dir '/' fname, 'r.mat'],'Data'); 
% TODO: convert following function into a class member function -cg
ballisticRelease_trialfy([file_dir '/' fname 'r.mat'], dst_file);

load(dst_file);
%% summary
trial_num = length(trial);
fth_all = sort(unique([trial(~isnan([trial.fth])).fth]));
tar_all = sort(unique([trial.target]));
tarl_all = sort(unique(nonzeros([trial.targetl]')),1,'descend'); % longer target easier
fth_num = length(fth_all);
tar_num = length(tar_all);
tarl_num = length(tarl_all);
SAMPLE_T = 0.02;
disp(['fth_all:     ' num2str(fth_all)]);
disp(['tar_all:     ' num2str(tar_all)]);
disp(['tarl_all:    ' num2str(tarl_all')]);


%% show sucess rate on each condition
success_table = zeros(fth_num,tar_num);
success_table_idx = cell(fth_num,tar_num);
for ii = 1:fth_num % force threshold
    for jj = 1:tar_num % target
        all_trials_num = sum([trial.fth] == fth_all(ii) & [trial.target] == tar_all(jj));
        fin_trials_num = sum([trial.fth] == fth_all(ii) & [trial.target] == tar_all(jj) & [trial.outcome] == 1);
        success_table(ii,jj) = fin_trials_num;
        success_table_idx{ii, jj} = find([trial.fth] == fth_all(ii) & [trial.target] == tar_all(jj) & [trial.outcome] == 1);
    end
end

%% align according to ST_MOV
trial = ballisticRelease_aligntrials(trial, ST_MOV, 'min');

%% show force on each condition
f_h1 = figure;
for ii = 1:fth_num % force threshold
   for jj = 1:tar_num % target
        % plot force here
        ii = 1; jj = 1;
        
        %subplot(fth_num,tar_num,(ii-1)*tar_num+jj);
        qtrials = trial([trial.fth] == fth_all(ii) & [trial.target] == tar_all(jj) & [trial.targetl] == tarl_all(1)); % qualified trials
        hold on;
        if length(qtrials) == 0, continue; end % no trials in this condition, skip
        % title and xylabel
        xlabel('time/s');
        ylabel('force/N');
        title_str = ['fth' num2str(fth_all(ii)) ', tar' num2str(tar_all(jj))];
        title(title_str);
        % x and y range; x(time/s)
        %xlim([0, length(qtrials(1).time)]);
        x_ct = length(qtrials(1).time)/2 + 1;
        xlim(([0,length(qtrials(1).time)]-0.5*length(qtrials(1).time))/50);
        %x_tick = [-50 0 50]' + x_ct;
        %x_ticklabel = cell(size(x_tick)); 
        x_tick = [-1 0 1]';
            for tt = 1:length(x_tick) 
                % x_ticklabel{tt,1} = num2str((x_tick(tt)-x_ct)*SAMPLE_T); 
            end
        xticks(x_tick); 
        %xticklabels(x_ticklabel);
        if (1)% plot each line
            for kk = 1:length(qtrials)
                x = qtrials(kk).time(1,:)-qtrials(kk).time(1,x_ct);
                plot(x, qtrials(kk).force(7,:),'r');
                plot(x, qtrials(kk).force(8,:),'g');
                plot(x, qtrials(kk).force(9,:),'b');
                plot(x, sqrt(qtrials(kk).force(7,:).^2+qtrials(kk).force(7,:).^2),'m');
            end
        else % plot mean
            f_mat = [qtrials.force]; 
            f_matr = reshape(f_mat', size(qtrials(1).force,2), length(qtrials), size(qtrials(1).force,1));
            plot(mean(f_matr(:,:,7)'),'r');
            plot(mean(f_matr(:,:,8)'),'g');
            plot(mean(f_matr(:,:,9)'),'b');
        end
    end
end
suptitle('force with task');



%% show trajectory on each condition
f_h2 = figure;
for ii = 1:fth_num % force threshold
    for jj = 1:tar_num % target
        % plot force here
        subplot(fth_num,tar_num,(ii-1)*tar_num+jj);
        hold on;
        qtrials = trial([trial.fth] == fth_all(ii) & [trial.target] == tar_all(jj) & [trial.targetl] == tarl_all(1)); % qualified trials
        if length(qtrials) == 0, continue; end % no trials in this condition, skip
        % title and xylabel
        xlabel('time/s');
        ylabel('pos/?m');
        title_str = ['fth' num2str(fth_all(ii)) ', tar' num2str(tar_all(jj))];
        title(title_str);
        % x and y range; x(time/s)
        xlim([0, length(qtrials(1).time)]);
        x_ct = length(qtrials(1).time)/2 + 1;
        x_tick = [-50 0 50]' + x_ct;
        x_ticklabel = cell(size(x_tick)); 
            for tt = 1:length(x_tick) 
                x_ticklabel{tt,1} = num2str((x_tick(tt)-x_ct)*SAMPLE_T); 
            end
        xticks(x_tick); xticklabels(x_ticklabel);
        if (1) % plot each line
            for kk = 1:length(qtrials)
                plot(qtrials(kk).positionAct(1,:),'r');
                plot(qtrials(kk).positionAct(2,:),'g');
                plot(qtrials(kk).positionAct(3,:),'b');
            end
        else % plot mean
            p_mat = [qtrials.positionAct]; 
            p_matr = reshape(p_mat', size(qtrials(1).positionAct,2), length(qtrials), size(qtrials(1).positionAct,1));
            plot(mean(p_matr(:,:,1)'),'r');
            plot(mean(p_matr(:,:,2)'),'g');
            plot(mean(p_matr(:,:,3)'),'b');
        end
    end
end
suptitle('position with task');

%% show velocity on each condition
f_h3 = figure;
for ii = 1:fth_num % force threshold
    for jj = 1:tar_num % target
        % plot force here
        subplot(fth_num,tar_num,(ii-1)*tar_num+jj);
        hold on;
        qtrials = trial([trial.fth] == fth_all(ii) & [trial.target] == tar_all(jj) & [trial.targetl] == tarl_all(1)); % qualified trials
        if length(qtrials) == 0, continue; end % no trials in this condition, skip

        % x and y range; x(time/s)
        xlim([0, length(qtrials(1).time)]);
        x_ct = length(qtrials(1).time)/2 + 1;
        x_tick = [-50 0 50]' + x_ct;
        x_ticklabel = cell(size(x_tick)); 
            for tt = 1:length(x_tick) 
                x_ticklabel{tt,1} = num2str((x_tick(tt)-x_ct)*SAMPLE_T); 
            end
        xticks(x_tick); xticklabels(x_ticklabel);
        if (1)% plot each line
            for kk = 1:length(qtrials)
                plot(qtrials(kk).velocity(1,:),'r');
                plot(qtrials(kk).velocity(2,:),'g');
                % plot(qtrials(kk).velocity(3,:),'b');
            end
        else % plot mean
            v_mat = [qtrials.velocity]; 
            v_matr = reshape(v_mat', size(qtrials(1).velocity,2), length(qtrials), size(qtrials(1).velocity,1));
            plot(mean(v_matr(:,:,1)'),'r');
            plot(mean(v_matr(:,:,2)'),'g');
            plot(mean(v_matr(:,:,3)'),'b');
        end
    end
end
suptitle('velocity with task');

% shade them

%%% still: need function to plot meaned;