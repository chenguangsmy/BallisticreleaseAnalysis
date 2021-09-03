%% look the difference between low-impedance robot setting

% see 3N series
ss2468 = SessionScan(2468);
ss2457 = SessionScan(2457);
ss2458 = SessionScan(2458);
ss2459 = SessionScan(2459);
sessions_all = [2457, 2458, 2459];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);

[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2457, axh2, col_array(1,:));
[axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2458, axh, col_array(2,:));
[axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2459, axh, col_array(3,:));
legend([lnh1, lnh2, lnh3], {'3N 5cm', '3N 7.5cm', '3N 10cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');

% 9N series
ss2469 = SessionScan(2469);
ss2460 = SessionScan(2460);
ss2461 = SessionScan(2461);
ss2462 = SessionScan(2462);
sessions_all = [2460, 2461, 2462];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);
[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2460, axh2, col_array(1,:));
[axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2461, axh, col_array(2,:));
[axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2462, axh, col_array(3,:));
legend([lnh1, lnh2, lnh3], {'9N 5cm', '9N 7.5cm', '9N 10 cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');

% 15N series
% ss2470 = SessionScan(2470);
% ss2456 = SessionScan(2456);
% ss2463 = SessionScan(2463);
% ss2455 = SessionScan(2455);
ss2479 = SessionScan(2479);
ss2480 = SessionScan(2480);
sessions_all = [2456, 2463, 2455];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);
% [axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2456, axh2, col_array(1,:));
% [axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2463, axh, col_array(2,:));
% [axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2455, axh, col_array(3,:));
[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2479, axh2, col_array(1,:));
[axh, val3, lnh2] = plotStepPertResponse_raw_subavg(ss2480, axh, col_array(3,:));
%legend([lnh1, lnh2, lnh3], {'15N 5cm', '15N 7.5cm', '15N 10 cm'});
legend([lnh1, lnh2], {'15N 5cm', '15N 10 cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');
pos_mean = mean([val2{1} val1{1}]); 
neg_mean = mean([val2{2} val1{2}]);
pos_1 = mean([val1{1}]); 
neg_1 = mean([val1{2}]);
val12_var = std(val1{2});
display('');
mean(val1{2}) - mean(val2{2})

% 21N series
ss2471 = SessionScan(2471);
ss2464 = SessionScan(2464);
ss2466 = SessionScan(2466);
ss2467 = SessionScan(2467);
sessions_all = [2464, 2466, 2467];

axh2 = figure();
tar_arr = [5 7.5 10];
col_array = colormap(lines);
[axh, val1, lnh1] = plotStepPertResponse_raw_subavg(ss2464, axh2, col_array(1,:));
[axh, val2, lnh2] = plotStepPertResponse_raw_subavg(ss2466, axh, col_array(2,:));
[axh, val3, lnh3] = plotStepPertResponse_raw_subavg(ss2467, axh, col_array(3,:));
legend([lnh1, lnh2, lnh3], {'21N 5cm', '21N 7.5cm', '21N 10 cm'});
xlabel('time after perturbation (s)');
ylabel('displacement');
title('raw data on step-perturbation');

%% figure for position against time categorized by target as subplots 
sessions_all = [[2468, 2457, 2458, 2459]; [2469, 2460, 2461, 2462]; [2470, 2456, 2463, 2455]; [2471, 2464, 2466, 2467]];
%sessions_all = [[2481, 2482, 2483]; [2484, 2486, 2487]; [2479, 2488, 2480]; [2489, 2490, 2492]]; % backward movement, forward perturbation
sessions_all_arr = sessions_all(:);
for session_i = 1:length(sessions_all_arr)
    sstmp = sessions_all_arr(session_i);
    %['ss' num2str(sstmp) ' = SessionScan(' num2str(sstmp) ');']
    eval(['ss' num2str(sstmp) ' = SessionScan(' num2str(sstmp) ');']);
end
%% plot...
force_all = 3:6:21; 
figure();
sessions_all = [[2468, 2457, 2458, 2459]; [2469, 2460, 2461, 2462]; [2470, 2456, 2463, 2455]; [2471, 2464, 2466, 2467]];
%sessions_all = [[2481, 2482, 2483]; [2484, 2486, 2487]; [2479, 2488, 2480]; [2489, 2490, 2492]]';
for fce_i = 1:size(sessions_all, 1)
    axh = subplot(1, size(sessions_all, 1), fce_i); hold on; 
    for dist_i = 1:size(sessions_all, 1)
        sstmp = eval(['ss' num2str(sessions_all(fce_i, dist_i))]);
        %axh = plotMeantrialPosPert(sstmp, axh, dist_i); % perturbation
        %axh = plotMeantrialForcePert(sstmp, axh, dist_i); % perturbation
        %axh = plotMeantrialVel(sstmp, axh, dist_i); % perturbation
        %plotMeantrialVel(sstmp, axh, dist_i); % perturbation
        plotMeantrialForce(sstmp, axh, dist_i);
        %axh = plotMeantrialPos(sstmp, axh, dist_i); % release
    end
    %ylim([-0.045, 0.005]); % pos, -
    %xlim([-0.1, 0.8]);
    %ylim([-0.005, 0.035]);
    xlim([-0.1, 0.7]); % perturb
    %ylim([-0.1, 0.5]); % velocity
    ylim([-4, 25]);    % force, -
    xlabel('time'); 
    if (fce_i == 1)
        %ylabel('perturbation response position');
        ylabel(' ');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(force_all(fce_i)) 'N']);
end

% figure for position against time categorized by target as subplots 
sessions_all = [[2468, 2457, 2458, 2459]; [2469, 2460, 2461, 2462]; [2470, 2456, 2463, 2455]; [2471, 2464, 2466, 2467]]';
%sessions_all = [[2481, 2482, 2483]; [2484, 2486, 2487]; [2479, 2488, 2480]; [2489, 2490, 2492]]';
dist_all = [0.025 0.05, 0.075, 0.10]; 
%dist_all = [0.05, 0.075, 0.10];
figure();
for dist_i = 1:size(sessions_all, 2)
    axh = subplot(1, size(sessions_all, 2), dist_i); hold on; 
    for fce_i = 1:size(sessions_all, 2)
        sstmp = eval(['ss' num2str(sessions_all(dist_i, fce_i))]);
        %axh = plotMeantrialPosPert(sstmp, axh, fce_i);
        %axh = plotMeantrialVel(sstmp, axh, fce_i); % perturbation
        plotMeantrialForce(sstmp, axh, fce_i);
        %axh = plotMeantrialForcePert(sstmp, axh, fce_i); % perturbation
    end
    %ylim([-0.045, 0.005]); % pos -
    %ylim([-0.005, 0.045]); % pos +
    xlim([-0.1, 0.7]); % perturb
    %ylim([-0.1, 0.5]); % velocity
    ylim([-4, 25]); % force - 
    xlabel('time'); 
    if (dist_i == 1)
        ylabel('perturbation response position');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(dist_all(dist_i)) 'm']);
end

%%
figure();
for dist_i = 1:size(sessions_all, 1)
    axh = subplot(1, size(sessions_all, 1), dist_i); hold on; 
    for fce_i = 1:size(sessions_all, 2)
        sstmp = eval(['ss' num2str(sessions_all(dist_i, fce_i))]);
        axh = plotMeantrialForcePert(sstmp, axh, fce_i);
    end
    %ylim([-0.045, 0.005]);
    ylim([-8, 3]);
    xlabel('time'); 
    if (dist_i == 1)
        ylabel('perturbation response position');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(dist_all(dist_i)) 'm']);
end
figure();
for fce_i = 1:size(sessions_all, 2)
    axh = subplot(1, size(sessions_all, 2), fce_i); hold on; 
    for dist_i = 1:size(sessions_all, 1)
        sstmp = eval(['ss' num2str(sessions_all(dist_i, fce_i))]);
        axh = plotMeantrialForcePert(sstmp, axh, dist_i);
    end
    %ylim([-0.045, 0.005]);
    ylim([-8, 3]);
    xlabel('time'); 
    if (dist_i == 1)
        ylabel('perturbation response position');
    end
    set(gca, 'Ygrid', 'on');
    title([num2str(force_all(fce_i)) 'N']);
end

%% %% Look through the trajectory of Burdet's experiment (using WAM) 
ss2455 = SessionScan(2455);
ss2456 = SessionScan(2456);
%% plot all the trajectory out
axh = figure();
lnh_lgd = zeros(2,1); 
arr_lgd = {'150/m','300/m'};
for trial_i = 1:length(ss2455.trials)
    [axh, lnh] = ss2455.trials(trial_i).plotRobotEndpointTraj(axh, 'b');
    if trial_i == 1
        lnh_lgd(1) = lnh;
    end
end
for trial_i = 1:length(ss2456.trials)
    [axh, lnh] = ss2456.trials(trial_i).plotRobotEndpointTraj(axh, 'r');
    if trial_i == 1
        lnh_lgd(2) = lnh;
    end
end
legend(lnh_lgd, arr_lgd);
title('trajectories of "insteability dyamics"');

%% plot SPRING TEST with perturbations, non-0 offset
% 2524, 2527, 2529, 2530
% color 
color_arr = colormap('lines');
% read data
ss2524 = SessionScan(2524); ss2527 = SessionScan(2527); 
ss2529 = SessionScan(2529); ss2530 = SessionScan(2530);
% plot out
figure();
ss2524.plotStepPertResponse_raw(subplot(1,2,1), color_arr(1,:));
ss2524.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(1,:));
figure();
ss2527.plotStepPertResponse_raw(subplot(1,2,1), color_arr(2,:));
ss2527.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(2,:));
figure();
ss2529.plotStepPertResponse_raw(subplot(1,2,1), color_arr(3,:));
ss2529.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(3,:));
figure();
ss2530.plotStepPertResponse_raw(subplot(1,2,1), color_arr(4,:));
ss2530.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(4,:));
% mean subtract overlap
% ??? to be continued...

% ss2536; ss2537; ss2538; ss2539; ss2541; ss2540
ss2536 = SessionScan(2536); ss2537 = SessionScan(2537); ss2538 = SessionScan(2538); 
ss2539 = SessionScan(2539); ss2541 = SessionScan(2541); ss2540 = SessionScan(2540); 
figure();
ss2536.plotStepPertResponse_raw(subplot(1,2,1), color_arr(1,:));
ss2536.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(1,:));
figure();
ss2537.plotStepPertResponse_raw(subplot(1,2,1), color_arr(2,:));
ss2537.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(2,:));
figure();
ss2538.plotStepPertResponse_raw(subplot(1,2,1), color_arr(3,:));
ss2538.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(3,:));
figure();
ss2539.plotStepPertResponse_raw(subplot(1,2,1), color_arr(4,:));
ss2539.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(4,:));
figure();
ss2541.plotStepPertResponse_raw(subplot(1,2,1), color_arr(5,:));
ss2541.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(5,:));
figure();
ss2540.plotStepPertResponse_raw(subplot(1,2,1), color_arr(6,:));
ss2540.plotStepPertResponse_rawF(subplot(1,2,2),color_arr(6,:));

%% spring perturbation test, with x0 start at 0 positions
 sessions_mat = [2558 2560 2562 2563 2564;   % 160N/m
                 2547 2546 2605 2545 2602;   % 320N/m, external spring
                 2548 2549 2606 2551 2552;   % 640N/m
                 2557 2556 2555 2554 2553;]; % 960N/m
%sessions_mat = [2611 2610 2612 2613;        % 160N/m
%                2614 2615 2617 2619];       % 320N/m
                
sessions_all = sessions_mat(:);
for session_i = 1:length(sessions_all)
     eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
         num2str(sessions_all(session_i)) ');']);
end
%%
color_arr = colormap('lines');
figure();
for stiffness_i = 1:size(sessions_mat,1)
    axh = subplot(1,4,stiffness_i); hold on;
    for ss_i = 1:size(sessions_mat,2)
        if sessions_mat(stiffness_i,ss_i) == 2550 || ...
           sessions_mat(stiffness_i,ss_i) == 2543 || ...
           sessions_mat(stiffness_i,ss_i) == 2544 % bad session
            continue
        end
        %num2str(sessions_mat(stiffness_i, ss_i)); % display which session
        eval(['sstmp = ss' num2str(sessions_mat(stiffness_i, ss_i)) ';']);
        %sstmp.plotStepPertResponse_rawF(axh,color_arr(ss_i,:));
        sstmp.plotStepPertResponse_raw(axh,color_arr(ss_i,:));
    end
    title(['stiffness ' num2str(SpringStiff_list(stiffness_i)) 'N/m']);
    set(gca, 'Ygrid', 'on');
    xlim([-0.1 1]);
    ylim([-0.07 0.01]); % posision
    % ylim([-5 30]); % force
end

%% spring ballistic release test, with x0 start at 0 positions
sessions_mat = [2584 2585 2586 2587;   % 320    N/m, external spring
                2588 2589 2590 2591;    % 160   N/m
                2592 2593 2594 2596;    % 106.7 N/m
                -1   2597 2598 2599];    % 89    N/m
% read data                
sessions_all = sessions_mat(:);
for session_i = 1:length(sessions_all)
 try
    eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
        num2str(sessions_all(session_i)) ');']);
 catch 
     display('unable to load this session!');
 end
end
%% plots 
SpringStiff_list = [320, 160, 107, 89]; % spring test with fixed stiffness
targets_list = [10, 7.5, 5.0, 2.5]/100; %targets position (m)

color_arr = colormap('lines');
figure();
%for stiffness_i = 1:size(sessions_mat,1)
for target_i = 1:length(targets_list)
    target_tmp = targets_list(target_i);
    axh = subplot(1,4,target_i); hold on;
    %for ss_i = 1:size(sessions_mat,2)
    for stiffness_i = 1:size(sessions_mat,1)
        if sessions_mat(stiffness_i,target_i) == -1 
            continue
        end
        %num2str(sessions_mat(stiffness_i, ss_i)); % display which session
        eval(['sstmp = ss' num2str(sessions_mat(stiffness_i, target_i)) ';']);
        %sstmp.plotSameTrial();
        %axh=sstmp.plotReleaseResponse_rawF(axh, color_arr(stiffness_i,:));
        %axh=sstmp.plotReleaseResponse_rawP(axh, color_arr(stiffness_i,:));
        axh=sstmp.plotReleaseResponse_rawV(axh, color_arr(stiffness_i,:));
    end
    title(['target ' num2str(targets_list(target_i)) 'm']);
    set(gca, 'Ygrid', 'on');
    xlim([-0.1 1.2]);
    %ylim([-0.01 0.16]);
    ylim([-1.1 1.2]);
end

%% Perturb position summary 
% Only focus on the front and back movement for now.
% Sessions choose from: 9N/21N * 0.05m & 0.10m. 
% |sessions: | 9N*10cm  | 9N*5cm   | 21N*10cm | 21N*5cm  |
% | -------- | -------- | -------- | -------- | -------- |
% |front     | 2462     | 2460     | 2467     | 2464     | 
% |back      | 2487     | 2484     | 2492     | 2489     |
% |left      | 2657     | 2655     | 2658     | 2661     |
% |right     | 2674     | 2672     | 2675     | 2676     |

sessions_mat = [2462 2460 2467 2464;    % front
                2487 2484 2492 2489;    % back
                2657 2655 2658 2661;    % left
                2674 2672 2675 2676;];  % right
% read data
sessions_all = sessions_mat(:);
for session_i = 1:length(sessions_all)
    try
        eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
            num2str(sessions_all(session_i)) ');']);
    catch
        display('unable to load this session!');
    end
end
% only look through the front now
pert_avg_mat = nan(size(sessions_mat));
pert_std_mat = nan(size(sessions_mat));
for direction_i = 1:size(sessions_mat, 1)
    for session_i = 1:size(sessions_mat, 2) 
        sstmp = eval(['ss' num2str(sessions_mat(direction_i,session_i))]);
        positions = sstmp.getPosPert();
        positions_nnan = positions(~isnan(positions));
        if positions_nnan < 0
            positions_nnan = -positions_nnan;
        end
        pert_avg_mat(direction_i, session_i) = mean(positions_nnan, 'omitnan');
        pert_std_mat(direction_i, session_i) = std(positions_nnan, 'omitnan');
    end
end
axh = figure(); hold on;
sessions_num = size(sessions_mat, 2);
% bar(1:sessions_num, pert_avg_mat(1,:)); 
% errorbar(1:sessions_num, pert_avg_mat(1,:), pert_std_mat(1,:));
% bar(5:4+sessions_num, pert_avg_mat(2,:));
% errorbar(5:4+sessions_num, pert_avg_mat(2,:), pert_std_mat(2,:));
bar(1:sessions_num, pert_avg_mat(:,:));% errorbar(pert_std_mat(1:2,:));
title('perturbed average position');


%% left movement subject 
ss2655 = SessionScan(2655); ss2657 = SessionScan(2657); 
ss2658 = SessionScan(2658); ss2661 = SessionScan(2661); 
sessions_all = [2655        2657        2658        2661];
color_arr = colormap('lines');
%for session_i = 1:size(sessions_mat,1)
    %axh = subplot(1,4,session_i); hold on;
    axh = figure();
    for ss_i = 1:size(sessions_all,2)
        %num2str(sessions_mat(stiffness_i, ss_i)); % display which session
        eval(['sstmp = ss' num2str(sessions_all(ss_i)) ';']);
        %sstmp.plotStepPertResponse_rawF(axh,color_arr(ss_i,:));
        sstmp.plotStepPertResponse_raw(axh,color_arr(ss_i,:));
    end
    %title(['stiffness ' num2str(SpringStiff_list(session_i)) 'N/m']);
    set(gca, 'Ygrid', 'on');
    xlim([-0.1 1]);
    ylim([-0.07 0.01]); % posision
    % ylim([-5 30]); % force
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mass test served for Scottt's method, predict that higher mass on subject side, lower initial force.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss2697 = SessionScan(2697); % no mass on
ss2700 = SessionScan(2700); % 277g
ss2698 = SessionScan(2698); % 522g
ss2699 = SessionScan(2699); % 799g
% in a single plot
axh1 = figure(1); 
color_arr = colormap('lines');
ss2697.plotReleaseResponse_rawF(axh1,[0.5 0.5 0.5])
ss2700.plotReleaseResponse_rawF(axh1,color_arr(1,:));
ss2698.plotReleaseResponse_rawF(axh1,color_arr(2,:));
ss2699.plotReleaseResponse_rawF(axh1,color_arr(3,:));
xlim([-0.1 2]); title('Force response with different subject mass');

axh2 = figure(2); 
axh21 = subplot(1,4,1); axh22 = subplot(1,4,2); axh23 = subplot(1,4,3); axh24 = subplot(1,4,4); 
ss2697.plotReleaseResponse_rawF(axh21,[0.5 0.5 0.5]); xlim([-0.05 0.1]);  ylim([0, 22]);
title('+0g'); yticks([0:5:20]); axh21.YGrid = 'on';
ss2700.plotReleaseResponse_rawF(axh22,color_arr(1,:)); xlim([-0.05 0.1]); ylim([0, 22]);
title('+277g'); ylabel(''); yticks([0:5:20]);axh22.YTickLabel=''; axh22.YGrid = 'on';
ss2698.plotReleaseResponse_rawF(axh23,color_arr(2,:)); xlim([-0.05 0.1]); ylim([0, 22]);
title('+522g'); ylabel(''); yticks([0:5:20]);axh23.YTickLabel=''; axh23.YGrid = 'on';
ss2699.plotReleaseResponse_rawF(axh24,color_arr(3,:)); xlim([-0.05 0.1]); ylim([0, 22]);
title('+799g'); ylabel(''); yticks([0:5:20]);axh24.YTickLabel=''; axh24.YGrid = 'on';

%% continue mass test, Using regression to diagonize the sudden force change. 
% 1. get raw data
sstmp = ss2699; %97, 70 98 99
y_raw = sstmp.trials(3).force_h(2,:);
yt_raw = sstmp.trials(3).force_t;
y = sstmp.trials(3).force_h(2,:);
y_t = sstmp.trials(3).force_t;
% 2. select time
t_min = 0; t_max = 2.5; 
y = y(y_t>=t_min & y_t<=t_max);
y_t = y_t(y_t>=t_min & y_t<=t_max);
subplot(4,1,1); 
plot(y_t, y);
subplot(4,1,2);
plot(y_t(1:end-1), diff(y));
% 3. find sudden change peak, select again
y_d = diff(y); 
[y_dmin, y_didx] = min(y_d);
duration = 2000;
y = y(y_didx:y_didx+2000-1);
y_t = y_t(y_didx:y_didx+2000-1) - y_t(y_didx); % new 0-nize
subplot(4,1,3);
plot(y_t, y);
% 4. regression model
% y = e^(\ksai t)*sin(\omega t + \theta)
B0 = -0.69/1;  % about 1 sec decrease half;
B3 = 15;
B1 = 1*2*pi; % The period is around 1
B2 = pi/2; % looks like its dropping value
X = y_t;
% 5. regression
myFit = NonLinearModel.fit(X,y, 'y ~ exp(b0*x1)*sin(b1*x1 + b2)*b3', [B0, B1, B2, B3]);
subplot(4,1,4);
plot(y_t, y, 'b');
hold on
plot(y_t', myFit.Fitted, 'r');
% 6. Get the answer 
myFit.Fitted(1)
B3
% plot (overlap) the raw force and the net force
figure(); hold on;
plot(yt_raw - y_t(y_didx), y_raw, 'b');
% get the before release value
y_raw_0 = mean(y_raw(yt_raw<0 & yt_raw>-0.1));
y_t_temp = yt_raw(yt_raw<0 & yt_raw>-0.1);
plot(yt_raw(yt_raw<0 & yt_raw>-0.1)', y_raw_0*ones(size(y_t_temp)),'r');
plot(y_t', myFit.Fitted,'r');
plot([y_t(1) y_t(1)], [y_raw_0 myFit.Fitted(1)],'r')
xlim([-0.1 2]);
legend('raw data', 'damped sinusoid fit')
title('force immediatly after release');
xlabel('time at release (s)'); ylabel('Censored force (N)');
%% continue mass test, collect forces, and plot out a figure 
%Force_0p = [];
%Force_0p = [Force_0p myFit.Fitted(1)];
Force_0n = 20*ones(size(Force_0p));
% Force_0n/Force0p = 1+Ms/Mr;
% Mr = Ms/(F(0-) - F(0+))/F(0+)
% (Force_0n - Force_0p)/Force_0p = m/Mr + m_add/Mr; % 1/Mr slope, m/Mr intercept
y_mass = ((Force_0n - Force_0p) ./ Force_0p);
x_mass = [0, 0.277, 0.522, 0.799]; % m_add
figure();
plot(y_mass, x_mass, '*'); 
xlim([0, 1.2]); ylim([-0.3, 0.9]);
lsline();
% m/Mr = 0.3549; 1/Mr = 1.413-0.3549;
Mr = 1/(1.413-0.3549); m = 0.3549*Mr;
ylabel('added mass (kg)'); 
xlabel('(F_{0-} - F_{0+}) /F_{0+}');
title('mass prediction from force')

%% Scott's method with update, Using the immediate release to know the mass

sstmp = ss2627;
for trial_i = 2:length(sstmp.trials)-1
    sstmp.trials(trial_i) = sstmp.trials(trial_i).predictStiffDampx0();
end
stiffness = [sstmp.trials.pred_K];
damping = [sstmp.trials.pred_D];
x0 = [sstmp.trials.pred_x0];
m = [sstmp.trials.pred_A];

%% stiffness varying 
ss_num = [2585 2589 2593 2597];
%ss_num =[2593 2597];
ref = [89 107 160 320];
color_arr = colormap('lines');
for session_i = 1:length(ss_num)
     eval(['sstmp = SessionScan(' num2str(ss_num(session_i)) ');']);
    for trial_i = 2:length(sstmp.trials)-1
        sstmp.trials(trial_i) = sstmp.trials(trial_i).predictStiffDampx0();
    end
    stiffness_all{session_i} = [sstmp.trials(2:end).pred_K];
    damping_all{session_i} = [sstmp.trials(2:end).pred_D];
    x0_all{session_i} = [sstmp.trials(2:end).pred_x0];
    m_all{session_i} = [sstmp.trials(2:end).pred_A];
    val_mean(session_i) = mean(stiffness_all{session_i});
    val_std(session_i) = std(stiffness_all{session_i});
end
axh=subplot(1,3,1); hold on;
bar(1:4, val_mean, 'FaceColor', color_arr(1,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks(ref);
axh.YGrid='on';
ylabel('stiffness (N/m)');
xticks('');

% x0 varying 
ss_num = [2584 2585 2586 2587];
ref = [2.5 5 7.5 10]/100;
for session_i = 1:length(ss_num)
     eval(['sstmp = SessionScan(' num2str(ss_num(session_i)) ');']);
    for trial_i = 2:length(sstmp.trials)-1
        sstmp.trials(trial_i) = sstmp.trials(trial_i).predictStiffDampx0();
    end
    stiffness_all{session_i} = [sstmp.trials(2:end).pred_K];
    damping_all{session_i} = [sstmp.trials(2:end).pred_D];
    x0_all{session_i} = [sstmp.trials(2:end).pred_x0];
    m_all{session_i} = [sstmp.trials(2:end).pred_A];
    val_mean(session_i) = mean(x0_all{session_i});
    val_std(session_i) = std(x0_all{session_i});
end
axh=subplot(1,3,2); hold on;
bar(1:4, val_mean, 'FaceColor', color_arr(2,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks(ref);
axh.YGrid='on';
ylabel('x_0 position (m)');
title('various prediction on regression method');
xticks('');

% mass varying
ss_num = [2697 2700 2698 2699];
ref = ([0 277 522 799]+300) / 1000;
for session_i = 1:length(ss_num)
     eval(['sstmp = SessionScan(' num2str(ss_num(session_i)) ');']);
    for trial_i = 2:length(sstmp.trials)-1
        sstmp.trials(trial_i) = sstmp.trials(trial_i).predictStiffDampx0();
    end
    stiffness_all{session_i} = [sstmp.trials(2:end).pred_K];
    damping_all{session_i} = [sstmp.trials(2:end).pred_D];
    x0_all{session_i} = [sstmp.trials(2:end).pred_x0];
    m_all{session_i} = [sstmp.trials(2:end).pred_A];
    val_mean(session_i) = mean(m_all{session_i});
    val_std(session_i) = std(m_all{session_i});
end
axh = subplot(1,3,3); hold on;
bar(1:4, val_mean, 'FaceColor', color_arr(3,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks(ref);
axh.YGrid='on';
ylabel('M_s (kg)');
xticks('');

%% 
% sanity check to make sure the above is right
% stiffness
ss_num = [2585 2589 2593 2597];
%ss_num =[2593 2597];
ref = [89 107 160 320];
color_arr = colormap('lines');
for session_i = 1:length(ss_num)
     eval(['ss' num2str(ss_num(session_i)) ' = SessionScan(' num2str(ss_num(session_i)) ');']);
end
for session_i = 1:length(ss_num)
     eval(['sstmp = ss' num2str(ss_num(session_i)) ';']);
    for trial_i = 2:length(sstmp.trials)-1
        sstmp.trials(trial_i) = sstmp.trials(trial_i).predictStiffDampx0();
    end
    stiffness_all{session_i} = [sstmp.trials(2:end).pred_K];
    damping_all{session_i} = [sstmp.trials(2:end).pred_D];
    x0_all{session_i} = [sstmp.trials(2:end).pred_x0];
    m_all{session_i} = [sstmp.trials(2:end).pred_A];
end
figure();
axh=subplot(1,4,1); hold on; % stiffness
for session_i = 1:4
    val_mean(session_i) = mean(stiffness_all{session_i});
    val_std(session_i) = std(stiffness_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(1,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks(ref);
ylim([0 350]);
axh.YGrid='on';
ylabel('stiffness (N/m)');
xticks('');

axh=subplot(1,4,2); hold on; % damping
for session_i = 1:4
    val_mean(session_i) = mean(damping_all{session_i});
    val_std(session_i) = std(damping_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(1,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([-5 0 5]);
ylim([-5 5]);
axh.YGrid='on';
ylabel('damping (Ns/m)');
xticks('');

axh=subplot(1,4,3); hold on; % mass
for session_i = 1:4
    val_mean(session_i) = mean(m_all{session_i});
    val_std(session_i) = std(m_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(1,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([0 0.3]);
ylim([0 2]);
axh.YGrid='on';
ylabel('mass (kg)');
xticks('');

axh=subplot(1,4,4); hold on; % x0 position
for session_i = 1:4
    val_mean(session_i) = mean(x0_all{session_i});
    val_std(session_i) = std(x0_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(1,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([0 0.075]);
ylim([0 0.1]);
axh.YGrid='on';
ylabel('x_0 (m)');
xticks('');

%%
% x0 varying 
ss_num = [2584 2585 2586 2587];
ref = [2.5 5 7.5 10]/100;
for session_i = 1:length(ss_num)
     eval(['ss' num2str(ss_num(session_i)) ' = SessionScan(' num2str(ss_num(session_i)) ');']);
end
for session_i = 1:length(ss_num)
     eval(['sstmp = ss' num2str(ss_num(session_i)) ';']);
    for trial_i = 2:length(sstmp.trials)-1
        sstmp.trials(trial_i) = sstmp.trials(trial_i).predictStiffDampx0(1);
    end
    stiffness_all{session_i} = [sstmp.trials(2:end).pred_K];
    damping_all{session_i} = [sstmp.trials(2:end).pred_D];
    x0_all{session_i} = [sstmp.trials(2:end).pred_x0];
    m_all{session_i} = [sstmp.trials(2:end).pred_A];
    val_mean(session_i) = mean(x0_all{session_i});
    val_std(session_i) = std(x0_all{session_i});
end
% axh=subplot(1,3,2); hold on;
% bar(1:4, val_mean, 'FaceColor', color_arr(2,:)); 
% for i = 1:4
%     errorbar(i, val_mean(i), val_std(i),'k');
% end
% yticks(ref);
% axh.YGrid='on';
% ylim([0 0.1]);
% ylabel('x_0 position (m)');
% title('various prediction on regression method');
% xticks('');
figure();
axh=subplot(1,4,1); hold on; % stiffness
for session_i = 1:4
    val_mean(session_i) = mean(stiffness_all{session_i});
    val_std(session_i) = std(stiffness_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(2,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([0 320]);
ylim([0 350]);
axh.YGrid='on';
ylabel('stiffness (N/m)');
xticks('');

axh=subplot(1,4,2); hold on; % damping
for session_i = 1:4
    val_mean(session_i) = mean(damping_all{session_i});
    val_std(session_i) = std(damping_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(2,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([-5 0 5]);
ylim([-5 5]);
axh.YGrid='on';
ylabel('damping (Ns/m)');
xticks('');

axh=subplot(1,4,3); hold on; % mass
for session_i = 1:4
    val_mean(session_i) = mean(m_all{session_i});
    val_std(session_i) = std(m_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(2,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([0 0.3]);
ylim([0 2]);
axh.YGrid='on';
ylabel('mass (kg)');
xticks('');

axh=subplot(1,4,4); hold on; % x0 position
for session_i = 1:4
    val_mean(session_i) = mean(x0_all{session_i});
    val_std(session_i) = std(x0_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(2,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks(ref);
ylim([0 0.1]);
axh.YGrid='on';
ylabel('x_0 (m)');
xticks('');

%% 
% mass varying
ss_num = [2697 2700 2698 2699];
ref = ([0 277 522 799]+300) / 1000;
for session_i = 1:length(ss_num)
     eval(['ss' num2str(ss_num(session_i)) ' = SessionScan(' num2str(ss_num(session_i)) ');']);
end
for session_i = 1:length(ss_num)
     eval(['sstmp = ss' num2str(ss_num(session_i)) ';']);
    for trial_i = 2:length(sstmp.trials)-1
        sstmp.trials(trial_i) = sstmp.trials(trial_i).predictStiffDampx0(1);
    end
    stiffness_all{session_i} = [sstmp.trials(2:end).pred_K];
    damping_all{session_i} = [sstmp.trials(2:end).pred_D];
    x0_all{session_i} = [sstmp.trials(2:end).pred_x0];
    m_all{session_i} = [sstmp.trials(2:end).pred_A];
    val_mean(session_i) = mean(m_all{session_i});
    val_std(session_i) = std(m_all{session_i});
end
% axh = subplot(1,3,3); hold on;
% bar(1:4, val_mean, 'FaceColor', color_arr(3,:)); 
% for i = 1:4
%     errorbar(i, val_mean(i), val_std(i),'k');
% end
% yticks(ref);
% axh.YGrid='on';
% ylabel('M_s (kg)');
% xticks('');
figure();
axh=subplot(1,4,1); hold on; % stiffness
for session_i = 1:4
    val_mean(session_i) = mean(stiffness_all{session_i});
    val_std(session_i) = std(stiffness_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(3,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([0 160]);
ylim([0 350]);
axh.YGrid='on';
ylabel('stiffness (N/m)');
xticks('');

axh=subplot(1,4,2); hold on; % damping
for session_i = 1:4
    val_mean(session_i) = mean(damping_all{session_i});
    val_std(session_i) = std(damping_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(3,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([-5 0 5]);
ylim([-5 5]);
axh.YGrid='on';
ylabel('damping (Ns/m)');
xticks('');

axh=subplot(1,4,3); hold on; % mass
for session_i = 1:4
    val_mean(session_i) = mean(m_all{session_i});
    val_std(session_i) = std(m_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(3,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks(ref);
ylim([0 2]);
axh.YGrid='on';
ylabel('mass (kg)');
xticks('');

axh=subplot(1,4,4); hold on; % x0 position
for session_i = 1:4
    val_mean(session_i) = mean(x0_all{session_i});
    val_std(session_i) = std(x0_all{session_i});
end
bar(1:4, val_mean, 'FaceColor', color_arr(3,:)); 
for i = 1:4
    errorbar(i, val_mean(i), val_std(i),'k');
end
yticks([0 0.12]);
ylim([0 0.15]);
axh.YGrid='on';
ylabel('x_0 (m)');
xticks('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable test of which damping parameter is suitable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 300ms perturbations with various damping of 0Ns/m to 20Ns/m
ss_num = [2707 2703 2704 2705 2706];

color_arr = colormap('lines');
for session_i = 1:length(ss_num)
     eval(['ss' num2str(ss_num(session_i)) ' = SessionScan(' num2str(ss_num(session_i)) ');']);
end
axh = figure();
xlim([-0.1 0.7]);
for session_i = 1:length(ss_num)
    eval(['sstmp = ss' num2str(ss_num(session_i)) ';']);
    axh_lines = sstmp.plotStepPertResponse_rawF(axh, color_arr(session_i,:));
end
xlim([-0.1 0.7]); title('robot damping from 0 to 20');
grid on;
axh = figure();
for session_i = 1:length(ss_num)
     eval(['sstmp = ss' num2str(ss_num(session_i)) ';']);    axh_lines = sstmp.plotStepPertResponse_rawV(axh, color_arr(session_i,:));
end
xlim([-0.1 0.7]); title('robot damping from 0 to 20');

%% 300ms perturbations with 10Ns/m damping and various pert_mag of [15N, 10N, 5N, 3N]
ss_num = [2704 2708 2709 2710]% 2711 2712];
axh = figure();
color_arr = colormap('lines');
for session_i = 1:length(ss_num)
     eval(['ss' num2str(ss_num(session_i)) ' = SessionScan(' num2str(ss_num(session_i)) ');']);
end
for session_i = 1:length(ss_num)
     eval(['sstmp = ss' num2str(ss_num(session_i)) ';']);    axh_lines = sstmp.plotStepPertResponse_raw_subavg(axh, color_arr(session_i,:));
end
xlim([-0.1 0.7]); title('robot pertF: 15, 10, 5, 3');
axh = figure();
for session_i = 1:length(ss_num)
     eval(['sstmp = ss' num2str(ss_num(session_i)) ';']);    axh_lines = sstmp.plotStepPertResponse_raw(axh, color_arr(session_i,:));
    axh_lines = sstmp.plotStepPertResponse_rawF(axh, color_arr(session_i,:));
end
xlim([-0.1 0.7]); title('robot pertF: 15, 10, 5, 3');
%% the pert magnitude vs response
%ss2721 = SessionScan(2721);
%ss2722 = SessionScan(2722);
fig1 = ss2722.plotStepPertResponse_raw_subavg(1, [0.5 0.5 0.5]);
axObjs = fig1.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
time2 = 0.7;
positions_y1 = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    [~, x_idx2] = min(abs(xdat-time2));
    positions_y1(data_i) = mean(ydat(x_idx1:x_idx2));
end
%% % 5N
%ss2723 = SessionScan(2723); %within 5N
fig2 = ss2723.plotStepPertResponse_rawFce_subavg(1, [0.5 0.0 0.5], 10); % 10Hz filter
fig3 = ss2723.plotStepPertResponse_raw_pertfce(1, [0.8 0.2 0.8]);
axObjs = fig2.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
time2 = 0.7;
force_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    [~, x_idx2] = min(abs(xdat-time2));
    force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
end
axObjs = fig3.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
force_command = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    force_command(data_i) = mean(ydat(x_idx1));
end

figure(); plot(force_command, -force_censor, '*'); hold on;
plot(force_command, force_command, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('force censored (N)'); 
legend('force censored', 'command reference');
title('force check at robot with no stiffness'); 

%
%%ss2724 = SessionScan(2724); %within 5N
%% ss2727 = SessionScan(2727); %within 5N
fig2 = ss2724.plotStepPertResponse_rawFce_subavg(1, [0.5 0.0 0.5], 10); % 10Hz filter
fig2 = ss2727.plotStepPertResponse_rawFce_subavg(fig2, [0.5 0.0 0.5], 10); % 10Hz filter
fig3 = ss2724.plotStepPertResponse_raw_pertfce(1, [0.8 0.2 0.8]);
fig3 = ss2727.plotStepPertResponse_raw_pertfce(fig3, [0.8 0.2 0.8]);
axObjs = fig2.Children;
dataObjs = axObjs.Children;
time1 = 0.5;
time2 = 0.6;
force_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    [~, x_idx2] = min(abs(xdat-time2));
    force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
end
axObjs = fig3.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
force_command = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    force_command(data_i) = mean(ydat(x_idx1));
end

figure(); plot(force_command, -force_censor, '*'); hold on;
plot(force_command, force_command, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('force censored (N)'); 
ylim([-0.3 0.3]);
legend('force censored', 'command reference');
title('force check at robot with no stiffness'); 

%%
% ss2725 = SessionScan(2725); %within 5N
% ss2726 = SessionScan(2726);
fig2 = ss2725.plotStepPertResponse_rawFce_subavg(1, [0.5 0.0 0.5], 1); % 10Hz filter
fig2 = ss2726.plotStepPertResponse_rawFce_subavg(fig2, [0.5 0.0 0.5], 1); % 10Hz filter
fig3 = ss2725.plotStepPertResponse_raw_pertfce(1, [0.8 0.2 0.8]);
fig3 = ss2726.plotStepPertResponse_raw_pertfce(fig3, [0.8 0.8 0.2]);
axObjs = fig2.Children;
dataObjs = axObjs.Children;
time1 = 1.2;
time2 = 1.4;
force_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    [~, x_idx2] = min(abs(xdat-time2));
    force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
end
axObjs = fig3.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
force_command = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    force_command(data_i) = mean(ydat(x_idx1));
end

figure(); plot(force_command, -force_censor, '*'); hold on;
plot(force_command, force_command, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('force censored (N)'); 
ylim([-0.3 0.3]);
legend('force censored', 'command reference');
title('force check at robot with no stiffness'); 

%% 
% with force offset 
% ss2732 = SessionScan(2732); %within 5N
% ss2726 = SessionScan(2726);
fig2 = ss2732.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fig3 = ss2732.plotStepPertResponse_raw_pertfce(1, [0.8 0.2 0.8]);

axObjs = fig2.Children;
dataObjs = axObjs.Children;
time1 = 0.68;
time2 = 0.7;
force_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    [~, x_idx2] = min(abs(xdat-time2));
    force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
end
axObjs = fig3.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
force_command = zeros(size(dataObjs));
f_offset = 10;
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    force_command(data_i) = mean(ydat(x_idx1));
end

figure(); plot(force_command+f_offset, -force_censor, '*'); hold on;
plot(force_command+f_offset, force_command+f_offset, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('force censored (N)'); 
%ylim([-0.3 0.3]);
legend('force censored', 'command reference');
title('force check at robot with no stiffness'); 

% ss2733 = SessionScan(2733);
fig2 = ss2733.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fig3 = ss2733.plotStepPertResponse_raw_pertfce(1, [0.8 0.2 0.8]);

axObjs = fig2.Children;
dataObjs = axObjs.Children;
time1 = 0.68;
time2 = 0.7;
force_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    [~, x_idx2] = min(abs(xdat-time2));
    force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
end
axObjs = fig3.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
force_command = zeros(size(dataObjs));
f_offset = -10;
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    force_command(data_i) = mean(ydat(x_idx1));
end

figure(); plot(force_command+f_offset, -force_censor, '*'); hold on;
plot(force_command+f_offset, force_command+f_offset, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('force censored (N)'); 
%ylim([-0.3 0.3]);
legend('force censored', 'command reference');
title('force check at robot with no stiffness'); 

%%
% ss2734 = SessionScan(2734);
fig2 = ss2734.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
xlim([-0.1 3]);
fig3 = ss2734.plotStepPertResponse_raw_pertfce(1, [0.8 0.2 0.8]);

axObjs = fig2.Children;
dataObjs = axObjs.Children;
time1 = 1.6;
time2 = 1.7;
force_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    [~, x_idx2] = min(abs(xdat-time2));
    force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
end
axObjs = fig3.Children;
dataObjs = axObjs.Children;
time1 = 0.6;
force_command = zeros(size(dataObjs));
f_offset = 0;
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    force_command(data_i) = mean(ydat(x_idx1));
end

figure(); plot(force_command+f_offset, -force_censor, '*'); hold on;
plot(force_command+f_offset, force_command+f_offset, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('force censored (N)'); 
%ylim([-0.3 0.3]);
legend('force censored', 'command reference');
title('force check at robot with no stiffness'); 