%function figure_compare_spring_subject()
% plot panels with spring teset as well as the subject being perturbed 
% 1. load data from processed .mat files
%       if failed, process from the raw data
% 2. plot using corresponding functions in SessionScan
%       and tidy up the plots

%% plot spring test when the xs0 and sr0 equilibrium at k0
addpath('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab');
interDataDir = '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData';
% force levels         | 3N |     9N |    15N |   21N |
force_list =            [3          9       15      21]; % before-release force level
sessions_mat_spring = [ 2611    2610    2612    2613;           %Ks = 160N/m
                        2614    2615    2617    2619];          %Ks = 320N/m
sessions_mat_human = [  2468    2469    2470    2471            %tar = 2.5cm
                        2457    2460    2456    2464            %tar = 5.0cm
                        2458    2461    2463    2466            %tar = 7.5cm
                        2459    2462    2455    2467];          %tar =10.0cm

% load files for spring test
sessions_all = sessions_mat_spring(:);
if ~exist('ss2611', 'var')
    try 
        disp('spring test files not in workspace, try read from disk...');
        %load('../data/processedData/ss2611Springtests.mat'); % problematic in loading
        load([interDataDir '/' 'ss2611Springtests.mat']); % problematic in loading
        no_matfile = 0;
    catch 
        no_matfile = 1;
    end
else 
    no_matfile = 0;
end
if no_matfile % read files fron raw data
    disp('subject test files not in disk, try load from raw files...');
    for session_i = 1:length(sessions_all)
        eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
            num2str(sessions_all(session_i)) ');']);
    end
    save([interDataDir '/' 'ss2611Springtests.mat'], 'ss26*');
end
% load files for human experiments 
sessions_all = sessions_mat_human(:);
if ~exist('ss2468', 'var')
    try 
        disp('subject test files not in workspace, try read from disk...');
        %load('../data/processedData/ss2455lowImpefront.mat'); % problematic in loading
        load([interDataDir '/' 'ss2455lowImpefront.mat']); % problematic in loading
        no_matfile = 0;
    catch 
        no_matfile = 1;
    end
else 
    no_matfile = 0;
end
if no_matfile
    disp('subject test files not in disk, try load from raw files...');
    for session_i = 1:length(sessions_all)
         eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
             num2str(sessions_all(session_i)) ');']);
    end
    % now delete bulk variables in the class, ssScan.ft, ssScan.wam
    for session_i = 1:length(sessions_all)
        eval(['ss' num2str(sessions_all(session_i)) '.ft=[];']);
        eval(['ss' num2str(sessions_all(session_i)) '.wam=[];']);
    end
    save([interDataDir '/' 'ss2455lowImpefront.mat'], 'ss24*');
end

% prepare to plot 
figure();
color_arr = colormap('lines');
% % 1st line, spring Ks = 160N/m
for force_i = 1:size(force_list,2)
    axh = subplot(3,4,force_i); hold on;
    eval(['sstmp = ss' num2str(sessions_mat_spring(1, force_i)) ';']);
    sstmp.plotStepPertResponse_raw(axh,color_arr(end,:)); 
    title(['force ' num2str(force_list(force_i)) 'N']);
    set(gca, 'Ygrid', 'on');
    ylabel('');
    if force_i == 1
        ylabel('K=160 position(m)');
    else 
        yticklabels({});
    end
    yticks([-0.05:0.01:0]);
    xlim([-0.1 1]);
    ylim([-0.05 0.01]); % posision
    xlabel('');
end
% % 2nd line, spring Ks = 320N/m
for force_i = 1:size(force_list,2)
    axh = subplot(3,4,4+force_i); hold on;
    eval(['sstmp = ss' num2str(sessions_mat_spring(2, force_i)) ';']);
    sstmp.plotStepPertResponse_raw(axh,color_arr(end,:)); % how to subtract the mean?
    title('');
    set(gca, 'Ygrid', 'on');
    ylabel('');
    if force_i == 1
        ylabel('K=320 position(m)');
    else
        yticklabels({});
    end
    yticks([-0.05:0.01:0]);
    xlim([-0.1 1]);
    ylim([-0.05 0.01]); % posision
    xlabel('');
end

% % 3nd line, the real experiment data
force_all = 3:6:21; 
disp_pert_levels = zeros(4,4);
disp_stady_idx = 399; % a random from 390~400

for fce_i = 1:size(sessions_mat_human, 2)
    axh = subplot(3,4, fce_i+8); hold on; 
    for dist_i = 1:size(sessions_mat_human, 1)
        sstmp = eval(['ss' num2str(sessions_mat_human(dist_i, fce_i))]);
        axh = plotMeantrialPosPert(sstmp, axh, dist_i); % perturbation
    end
    ylim([-0.05 0.01]);
    xlim([-0.1 1]);
    xlabel('time'); 
    title('');
    ylabel('');
    if (fce_i == 1)
        ylabel('subject position(m)');
    else
        yticklabels({});
    end
    yticks([-0.05:0.01:0]);
    set(gca, 'Ygrid', 'on');
    legend off;
    %title([num2str(force_all(fce_i)) 'N']);
    if (fce_i == 4)
        label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
        label_names = {'2.5cm', '5cm', '7.5cm', '10cm'};
        legend(label_handles, label_names);
    end
    
    % assume stiffness is the before-release (force/disp-stiff_robot)
    label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
    lines_x = [label_handles.XData];
    lines_y = [label_handles.YData];
    data_length = length(label_handles(1).XData);
    lines_x = reshape(lines_x, data_length, 4);
    lines_y = reshape(lines_y, data_length, 4);
    disp_pert_levels(:,fce_i) = lines_y(disp_stady_idx,:)';  %%% do it later
end

stiffness_pert_sanity = -15./disp_pert_levels - 300*ones(4,4);

% sanity checks
varname = who;
%plot_force_sanitycheck(varname);
%stiffness_sanityCheck = stiffnessSanityCheck();

%end
%% a Force sanity check, whether the before-perturbation exerted force is with task condition 
%function plot_force_sanitycheck(varname)

%for i = 1:size(varname,1)
%    if iscellstr(varname)
%      evalin('base',[varname{i} ' = ' mat2str(i) ';'])
%    elseif ischar(varname)
%      evalin('base',[deblank(varname(i,:)) ' = ' mat2str(i) ';'])
%    end
%end

figure();
% first row
for force_i = 1:size(force_list,2)
    axh = subplot(3,4,force_i); hold on;
    eval(['sstmp = ss' num2str(sessions_mat_spring(1, force_i)) ';']);
    sstmp.plotMeantrialForcePert(axh,14);
    title(['force ' num2str(force_list(force_i)) 'N']);
    set(gca, 'Ygrid', 'on');
    ylabel('');
    if force_i == 1
        ylabel('K=160 force(N)');
    else
        yticklabels({});
    end
    xlim([-0.1 1]);
    %ylim([-0.05 0.01]); % posision
    ylim([-20 40]); % posision
    xlabel('');
end
% 2nd row
for force_i = 1:size(force_list,2)
    axh = subplot(3,4,4+force_i); hold on;
    eval(['sstmp = ss' num2str(sessions_mat_spring(2, force_i)) ';']);
    sstmp.plotMeantrialForcePert(axh,14); % how to subtract the mean?
    title('');
    set(gca, 'Ygrid', 'on');
    ylabel('');
    if force_i == 1
        ylabel('K=320 force(N)');
    else
        yticklabels({});
    end
    xlim([-0.1 1]);
    %ylim([-0.05 0.01]); % posision
    ylim([-20 40]); % posision
    xlabel('');
end

for fce_i = 1:size(sessions_mat_human, 1)
    axh = subplot(3,4, fce_i+8); hold on;
    for dist_i = 1:size(sessions_mat_human, 1)
        sstmp = eval(['ss' num2str(sessions_mat_human(fce_i, dist_i))]);
        axh = plotMeantrialForcePert(sstmp, axh, dist_i); % perturbation
    end
    ylim([-20 40]);
    xlim([-0.1 1]);
    xlabel('time');
    title('');
    ylabel('');
    if (fce_i == 1)
        ylabel('subject exerting force(N)');
    else
        yticklabels({});
    end
    set(gca, 'Ygrid', 'on');
    legend off;
    %title([num2str(force_all(fce_i)) 'N']);
    if (fce_i == 4)
        label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
        label_names = {'2.5cm', '5cm', '7.5cm', '10cm'};
        legend(label_handles, label_names);
    end
    
end
%end

%% ballistic release sanity check, whether the force/position variate with the stiffness requirements
%function stiffness_sanityCheck = stiffnessSanityCheck()
displacement_levels = zeros(4,4);
force_levels = zeros(4,4);
sessions_mat_human = [[2468, 2457, 2458, 2459]; [2469, 2460, 2461, 2462]; [2470, 2456, 2463, 2455]; [2471, 2464, 2466, 2467]];
sessions_all = sessions_mat_human(:);
for session_i = 1:length(sessions_all)
    eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
        num2str(sessions_all(session_i)) ');']);
end
for fce_i = 1:size(sessions_mat_human, 1)
    axh = subplot(2,4, fce_i); hold on; 
    for dist_i = 1:size(sessions_mat_human, 1)
        sstmp = eval(['ss' num2str(sessions_mat_human(fce_i, dist_i))]);
        plotMeantrialPos(sstmp, axh, dist_i); % perturbation
    end
    ylim([-0.01 0.15]);
    xlim([-0.1 1]);
    xlabel('time'); 
    title('');
    ylabel('');
    if (fce_i == 1)
        ylabel('subject exerting force(N)');
    else
        yticklabels({});
    end
    set(gca, 'Ygrid', 'on');
    legend off;
    %title([num2str(force_all(fce_i)) 'N']);
    pos_stady_idx = 600;
    if (fce_i == 4)
        label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
        label_names = {'2.5cm', '5cm', '7.5cm', '10cm'};
        legend(label_handles, label_names);
    end
    % get data from the average, to sanity check the perturbation 
    axh = subplot(2,4, fce_i);
    label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
    lines_x = [label_handles.XData];
    lines_y = [label_handles.YData];
    data_length = length(label_handles(1).XData);
    lines_x = reshape(lines_x, data_length, 4);
    lines_y = reshape(lines_y, data_length, 4);
    displacement_levels(:,fce_i) = lines_y(pos_stady_idx,:)';  %%% do it later
end

for fce_i = 1:size(sessions_mat_human, 1)
    axh = subplot(2,4, fce_i+4); hold on; 
    for dist_i = 1:size(sessions_mat_human, 1)
        sstmp = eval(['ss' num2str(sessions_mat_human(fce_i, dist_i))]);
        plotMeantrialForce(sstmp, axh, dist_i); % perturbation
    end
    ylim([0 25]);
    xlim([-0.1 1]);
    xlabel('time'); 
    title('');
    ylabel('');
    if (fce_i == 1)
        ylabel('subject exerting force(N)');
    else
        yticklabels({});
    end
    set(gca, 'Ygrid', 'on');
    legend off;
    %title([num2str(force_all(fce_i)) 'N']);
    fce_stady_idx = 20;
    if (fce_i == 4)
        label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
        label_names = {'2.5cm', '5cm', '7.5cm', '10cm'};
        legend(label_handles, label_names);
    end
    axh = subplot(2,4, fce_i+4);
    label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
    lines_x = [label_handles.XData];
    lines_y = [label_handles.YData];
    data_length = length(label_handles(1).XData);
    lines_x = reshape(lines_x, data_length, 4);
    lines_y = reshape(lines_y, data_length, 4);
    force_levels(:,fce_i) = lines_y(fce_stady_idx,:)';  %%% do it later
        
end

stiffness_sanityCheck = force_levels./displacement_levels;

figure(); hold on;
plot(stiffness_pert_sanity(:), stiffness_sanityCheck(:), '*');
plot(1:900, 1:900);
title('sanity check stiffness');
xlabel('overall - robot (N/m)');
ylabel('force\_threshold/ stop\_osition (N/m)');

%end