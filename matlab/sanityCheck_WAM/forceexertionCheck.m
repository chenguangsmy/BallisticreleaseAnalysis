sd% force exertion check
% Generating force using WAM and check if WAM is generating the wanted force 

% sessions and its description:
% large scales:
%   | ss2722 | -19:2:19N            | count up   |
%   | ss2729 | -19:2:19N            | count down |
%   | ss2734 | [-40:-21, 21:40]N    | count down |

% small scales:
%   | ss2723 | -5:0.25:5N           | count up   |
%   | ss2724 | -1.25:0.0626:1.25N   | count up   |
%   | ss2725 | -0.4:0.02:0.4 N      | count up   |
%   | ss2728 | -5:0.25:5N           | count down |
%   | ss2727 | -1.25:0.0626:1.25N   | count down |
%   | ss2726 | -0.4:0.02:0.4 N      | count down |

% overlap force: (already exerting force and generate another)
%   | ss2732 | 10+(-20:1:20)N       | count down |
%   | ss2733 | -10+(-20:1:20)N      | count down |

%% 1. plot large scales 
ss_list = [2740]; % may inaccurate as the force at 0.
% 2722 2729 2734
for session_i = 1:length(ss_list)
     eval(['ss' num2str(ss_list(session_i)) ' = SessionScan(' num2str(ss_list(session_i)) ');']);
end

%% plot raw data
% fig01 = ss2722.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
% fig01 = ss2734.plotStepPertResponse_rawFce(fig01, [0.5 0.0 0.5], 10); % 10Hz filter
% fig01 = ss2729.plotStepPertResponse_rawFce(fig01, [0.5 0.0 0.5], 10); % 10Hz filter
% fig02 = ss2722.plotStepPertResponse_raw_pertfce(1, [0.8 0.2 0.8]);
% fig02 = ss2734.plotStepPertResponse_raw_pertfce(fig02, [0.8 0.2 0.8]);
% fig02 = ss2729.plotStepPertResponse_raw_pertfce(fig02, [0.8 0.2 0.8]);
close all; clc;
figure();
fig01 = subplot(3,1,1);
fig01 = ss2739.plotStepPertResponse_rawFce(fig01, [0.5 0.0 0.5], 10); % 10Hz filter
%fig01 = ss2737.plotStepPertResponse_rawFce(fig01, [0.5 0.0 0.5], 10);
xlim([-3 3]);
fig02 = subplot(3,1,2);
fig02 = ss2739.plotStepPertResponse_raw_pertfce(fig02, [0.8 0.2 0.8]);
%fig02 = ss2737.plotStepPertResponse_raw_pertfce(fig02, [0.8 0.2 0.8]);
xlim([-3 3]);
fig03 = subplot(3,1,3);
fig03 = ss2739.plotStepPertResponse_raw(fig03, [0.5 0.0 0.5]); % 10Hz filter
%fig03 = ss2737.plotStepPertResponse_raw(fig03, [0.5 0.0 0.5]);
xlim([-3 3]);
% function, extract values from plot...
% extract values
%axObjs = fig01.Children;
dataObjs = fig01.Children;
time1 = [0.6 0.7];
time2 = [1.7 1.7];
force_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time2(1)));
    [~, x_idx2] = min(abs(xdat-time2(2))); 
    force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
    if force_censor(data_i) < 5 % small forces use different time
        [~, x_idx1] = min(abs(xdat-time1(1)));
        [~, x_idx2] = min(abs(xdat-time1(2))); 
        force_censor(data_i) = mean(ydat(x_idx1:x_idx2));
    end
end
%axObjs = fig03.Children;
dataObjs = fig03.Children;
time1 = [0.6 0.7];
time2 = [1.7 1.7];
position_censor = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time2(1)));
    [~, x_idx2] = min(abs(xdat-time2(2))); 
    [~, x_idx0] = min(abs(xdat - 0)); %before release 
    position_censor(data_i) = mean(ydat(x_idx1:x_idx2));
    if abs(position_censor(data_i) - 0.48) < 0.1 % small forces use different time
        [~, x_idx1] = min(abs(xdat-time1(1)));
        [~, x_idx2] = min(abs(xdat-time1(2))); 
        position_censor(data_i) = mean(ydat(x_idx1:x_idx2)) - ydat(x_idx0);
    end
end
%axObjs = fig02.Children;
dataObjs = fig02.Children;
f_offset = -20;
time1 = 0.6;
force_command = zeros(size(dataObjs));
for data_i = 1:length(dataObjs)
    xdat = dataObjs(data_i).XData;
    ydat = dataObjs(data_i).YData;
    [~, x_idx1] = min(abs(xdat-time1));
    force_command(data_i) = mean(ydat(x_idx1)) + f_offset;
end
% summary data and plot 
figure(); 
subplot(3,1,1); hold on;
plot(force_command, -force_censor, '*'); 
plot(force_command, force_command, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('force censored (N)'); 
legend('force censored', 'command reference');
title('force check at robot with no stiffness'); 
subplot(3,1,2); hold on;
plot(force_command, (-force_censor - force_command), '*'); 
title('force generation error');
xlabel('force command (N)'); ylabel('force error (N)'); 
grid on;
subplot(3,1,3); hold on;
plot(force_command, position_censor, '*'); 
plot(force_command, force_command/960, 'Color', [0.5 0.5 0.5]); 
legend('sensored position', 'predict by spring and FT');
%plot(force_command, force_command, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command (N)'); ylabel('position censored (m)'); 
title('force check at robot with no stiffness'); 
%% look at the position
% ss2736 = SessionScan(2736);
fce_cmd = ss2736.getPertFce_cmd(); 
fce_csd = -ss2736.getPertFce(10);
pos_csd = ss2736.getPertPos();
figure();
subplot(1,2,1); 
plot(fce_cmd, fce_csd, '*'); hold on;
plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
xlabel('force command'); ylabel('force measured');
subplot(1,2,2);
plot(fce_csd, pos_csd, '*'); hold on;
hline = refline();
stiff_posind = (hline.XData(2) - hline.XData(1))/(hline.YData(2) - hline.YData(1));
text(0,0.46,['slope=' num2str(stiff_posind)]);
%plot(fce_cmd, fce_csd/960, 'Color', [0.5 0.5 0.5]); hold on;
xlabel('force command'); ylabel('position measured');

figure();

posidx = fce_cmd>0; negidx = fce_cmd<0;
subplot(2,2,1)
plot(pos_csd(posidx), fce_cmd(posidx),'r*');
subplot(2,2,2)
plot(pos_csd(negidx),fce_cmd(negidx), 'r*');
subplot(2,2,3)
plot( pos_csd(posidx), fce_csd(posidx),'b*');
subplot(2,2,4)
plot(pos_csd(negidx), fce_csd(negidx), 'b*');
stiff_posind = (hline.XData(2) - hline.XData(1))/(hline.YData(2) - hline.YData(1));

plot(fce_csd, pos_csd, 'b*'); hold on;
legend('fce_cmd', 'fce_csd');

%% ss2739
ss2739 = SessionScan(2739);
fce_cmd = ss2739.getPertFce_cmd() +20; 
fce_csd = -ss2739.getPertFce(5);
pos_csd = ss2739.getPertPos();
figure();
subplot(1,2,1); 
plot(fce_cmd, fce_csd, '*'); hold on;
plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
xlabel('force command'); ylabel('force measured');
subplot(1,2,2);
plot(fce_csd, pos_csd, '*'); hold on;
hline = refline();
stiff_posind = (hline.XData(2) - hline.XData(1))/(hline.YData(2) - hline.YData(1));
text(20,0.50,['slope=' num2str(stiff_posind)]);
%plot(fce_cmd, fce_csd/960, 'Color', [0.5 0.5 0.5]); hold on;
xlabel('force command'); ylabel('position measured');

%% ss2738 --- not working
ss2738 = SessionScan(2738);
fce_cmd = ss2738.getPertFce_cmd() -20; 
fce_csd = -ss2738.getPertFce(5);
pos_csd = ss2738.getPertPos();
figure();
subplot(1,2,1); 
plot(fce_cmd, fce_csd, '*'); hold on;
plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
xlabel('force command'); ylabel('force measured');
title('command force vs censored');
subplot(1,2,2);
plot(fce_csd, pos_csd, '*'); hold on;
hline = refline();
stiff_posind = (hline.XData(2) - hline.XData(1))/(hline.YData(2) - hline.YData(1));
text(-20,0.46,['slope=' num2str(stiff_posind)]);
title('censored force vs censored position');
%plot(fce_cmd, fce_csd/960, 'Color', [0.5 0.5 0.5]); hold on;
xlabel('force command'); ylabel('position measured');

%% ss2741: clamping the robot, ss2743, better clamp
ss2743 = SessionScan(2743);
ss2743.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fce_cmd = ss2743.getPertFce_cmd(); 
fce_csd = -ss2743.getPertFce(100);
pos_csd = ss2743.getPertPos(5);
figure();
plot(fce_cmd, fce_csd, '*', 'MarkerSize', 5); hold on;
plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('force measured');
title('command force vs censored');
figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');

%% ss2744: clamping the robot, with block force transducer
ss2749 = SessionScan(2749);
ss2749.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fce_cmd = ss2749.getPertFce_cmd(); 
fce_csd = -ss2749.getPertFce(100);
pos_csd = ss2749.getPertPos(5);
figure();
plot(fce_cmd, fce_csd, '*', 'MarkerSize', 5); hold on;
plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
fce_cmd_pos = fce_cmd(fce_cmd>0);
fce_cmd_neg = fce_cmd(fce_cmd<0);
fce_csd_pos = fce_csd(fce_cmd>0);
fce_csd_neg = fce_csd(fce_cmd<0);
P1 = polyfit(fce_cmd_pos, fce_csd_pos,1);
P2 = polyfit(fce_cmd_neg, fce_csd_neg,1);
fcecsd_fit1 = P1(1)*fce_cmd_pos + P1(2);
fcecsd_fit2 = P2(1)*fce_cmd_neg + P2(2);
plot(fce_cmd_pos, fcecsd_fit1); 
plot(fce_cmd_neg, fcecsd_fit2); 
legend('datapoints', 'commands', 'pos points fit', 'neg points fit');
grid on;
text(fce_cmd_pos(2), fce_csd_pos(2), ['slope:' num2str(P1)]);
text(fce_cmd_neg(2), fce_csd_neg(2), ['slope:' num2str(P2)]);
xlabel('force command'); ylabel('force measured');
title('command force vs censored');

figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');

%% ss2750: clamping the robot, with change the toolplate length 50 41cm, turned off the jaccobian update
%ss2750 = SessionScan(2750);
ss2750.plotStepPertResponse_rawFce(1, [0.5 0.5 0.5]); % 10Hz filter
fce_cmd = ss2750.getPertFce_cmd(); 
fce_csd = -ss2750.getPertFce(200);
pos_csd = ss2750.getPertPos(500);
figure(); subplot(1,2,1);
plot(fce_cmd, fce_csd, 'r*', 'MarkerSize', 5); hold on;
plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
fce_cmd_pos = fce_cmd(fce_cmd>0);
fce_cmd_neg = fce_cmd(fce_cmd<0);
fce_csd_pos = fce_csd(fce_cmd>0);
fce_csd_neg = fce_csd(fce_cmd<0);
P1 = polyfit(fce_cmd_pos, fce_csd_pos,1);
P2 = polyfit(fce_cmd_neg, fce_csd_neg,1);
fcecsd_fit1 = P1(1)*fce_cmd_pos + P1(2);
fcecsd_fit2 = P2(1)*fce_cmd_neg + P2(2);
%plot(fce_cmd_pos, fcecsd_fit1); 
%plot(fce_cmd_neg, fcecsd_fit2); 
%legend('datapoints', 'commands');%, 'pos points fit', 'neg points fit');
legend('F_{csd}', 'F_{csd} == F_{pert}');
grid on;
%text(fce_cmd_pos(2), fce_csd_pos(2), ['slope:' num2str(P1)]);
%text(fce_cmd_neg(2), fce_csd_neg(2), ['slope:' num2str(P2)]);
%xlabel('force command'); ylabel('force measured');
%title('command force vs censored');
xlabel('F_{pert} (N)'); ylabel('F_{csd} (N)');
title('perturbed vs censored force');
subplot(1,2,2)
stem(fce_cmd, fce_csd-fce_cmd); grid on;
%xlabel('force command'); ylabel('force difference');
%title('force error');
xlabel('F_{pert} (N)'); ylabel('F_{csd} - F_{pert} (N)');
title('force error');
figure();
%plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
stem(fce_cmd, pos_csd , 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
%xlabel('force command'); ylabel('position measured');
xlabel('F_{pert} (N)'); ylabel('x (m)');
title('command force vs position'); 
ss2750.plotTaskEndpointPositionh(); grid on;
ylim([0.478 0.496]);
title('measured endpoint position while arm clamped')
ylabel('x (m)'); xlabel('time (s)');
%% ss2751: clamping the robot, with change the toolplate length 50 41cm, turned on the update jaccobian
ss2751 = SessionScan(2751);
ss2751.plotStepPertResponse_rawFce(1, [0 0 0]); % 10Hz filter
grid on; ylim([0 80]); xlim([-0.1 0.6]);
fce_cmd = ss2751.getPertFce_cmd(); 
fce_csd = -ss2751.getPertFce(100);
pos_csd = ss2751.getPertPos(5);
figure();
plot(fce_cmd, fce_csd, 'b*', 'MarkerSize', 5); hold on;
plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
fce_cmd_pos = fce_cmd(fce_cmd>0);
fce_cmd_neg = fce_cmd(fce_cmd<0);
fce_csd_pos = fce_csd(fce_cmd>0);
fce_csd_neg = fce_csd(fce_cmd<0);
P1 = polyfit(fce_cmd_pos, fce_csd_pos,1);
P2 = polyfit(fce_cmd_neg, fce_csd_neg,1);
fcecsd_fit1 = P1(1)*fce_cmd_pos + P1(2);
fcecsd_fit2 = P2(1)*fce_cmd_neg + P2(2);
%plot(fce_cmd_pos, fcecsd_fit1); 
%plot(fce_cmd_neg, fcecsd_fit2); 
legend('datapoints', 'commands');%, 'pos points fit', 'neg points fit');
grid on;
%text(fce_cmd_pos(2), fce_csd_pos(2), ['slope:' num2str(P1)]);
%text(fce_cmd_neg(2), fce_csd_neg(2), ['slope:' num2str(P2)]);
xlabel('force command'); ylabel('force measured');
title('command force vs censored');

figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');

%% ss2754: clamping the robot sleeve, testing position accuracy
ss2754 = SessionScan(2754);
ss2754.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fce_cmd = ss2754.getPertFce_cmd(); 
fce_csd = -ss2754.getPertFce(100);
pos_csd = ss2754.getPertPos(5);
figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');
figure();
ss2754.plotTaskEndpointPositionh(); grid on;
ylim([0.474 0.492]);
xlim([0 350]);

%% ss2755: clamping the robot sleeve, testing position accuracy
ss2755 = SessionScan(2755);
ss2755.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fce_cmd = ss2755.getPertFce_cmd(); 
fce_csd = -ss2755.getPertFce(100);
pos_csd = ss2755.getPertPos(5);
figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');
ss2755.plotTaskEndpointPositionh(); grid on;
ylim([0.474 0.492]);
xlim([0 350]);


%% ss2757: clamping the robot sleeve, near side <-->ss2754
ss2757 = SessionScan(2757);
ss2757.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fce_cmd = ss2757.getPertFce_cmd(); 
fce_csd = -ss2757.getPertFce(100);
pos_csd = ss2757.getPertPos(5);
figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');
ss2757.plotTaskEndpointPositionh(); grid on;
ylim([0.474 0.492]);
xlim([0 350]);

%% ss2758: clamping the robot sleeve, after autotention <---> 2757
ss2758 = SessionScan(2758);
ss2758.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5], 10); % 10Hz filter
fce_cmd = ss2758.getPertFce_cmd(); 
fce_csd = -ss2758.getPertFce(100);
pos_csd = ss2758.getPertPos(5);
figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');
ss2758.plotTaskEndpointPositionh(); grid on;
ylim([0.474 0.492]);
xlim([0 350]);

%% ss2767: positions, 2771, positions without joint position torque
%  Step perturbation, various dx0: -8-8cm, Kx = 300N/m
ss2767 = SessionScan(2767);
ss2767.plotTaskEndpointPositionh()
ss2767.plotStepPertResponse_raw(1, [0 0 0]); 
xlim([-2 2]); grid on;title('ss2767');
ss2767.plotStepPertResponse_rawV(1, [0 0 0]); 
xlim([-2 2]); grid on;title('ss2767');
ss2767.plotTaskjointTorqeh()
ylim([-5 5])

%% ss2769
% Step perturbation, various Kx, B=20Ns/m, dx0 = 3cm. 
% Kx = 0:300:2700N/m
ss2769 = SessionScan(2769);
ss2769.plotTaskEndpointPositionh()
ss2769.plotStepPertResponse_raw(1, [0 0 0]); 
xlim([-2 2]); grid on;title('ss2769');
ss2769.plotStepPertResponse_rawV(1, [0 0 0]); 
xlim([-2 2]); grid on; title('ss2769');
ss2769.plotTaskjointTorqeh();
ylim([-5 5]); title('ss2769');

ss2770 = SessionScan(2770); % here no Kq
ss2770.plotTaskEndpointPositionh()
%%%%%-----
axh2770 = ss2770.plotStepPertResponse_raw(1, [0 0 0]); 
xlim([-0.1 0.6]); grid on; hold on;
axh_lines = axh2770.Children.Children;
x_ref = [-0.1:0.002:0.6]; y = zeros(size(x_ref));
y(x_ref>0) = 0.03;
axh_linesref = plot(x_ref, y, 'r', 'LineWidth', 2);
legend([axh_linesref, axh_lines(1)], '\Delta x_0', '\Delta x');
xlabel('time (s)'); ylabel('absolute position (m)');
%%%%% -----
ss2770.plotStepPertResponse_rawV(1, [0 0 0]); 
xlim([-2 2]); grid on;title('ss2770');
ss2770.plotTaskjointTorqeh();title('ss2770');
ylim([-5 5])

ss2772 = SessionScan(2772); % here no Kq
ss2772.plotStepPertResponse_rawV(1, [0 0 0]); 
xlim([-2 2]); grid on;title('ss2772');
ss2771 = SessionScan(2771); % here no Kq
ss2771.plotTaskEndpointPositionh()
ss2771.plotStepPertResponse_raw(1, [0 0 0]); 
xlim([-2 2]); grid on;
ss2771.plotStepPertResponse_rawV(1, [0 0 0]); 
xlim([-2 2]); grid on;title('ss2771');
ss2771.plotTaskjointTorqeh(); title('ss2771');
ylim([-5 5])

axh = ss2773.plotTaskjointTorqeh(); title('ss2773');
ss2773.plotAddTrialMark(axh.Children(1));
ylim([-50 50])

axh = ss2774.plotTaskjointTorqeh(); title('ss2774');
ss2774.plotAddTrialMark(axh.Children(1));
ylim([-50 50])

axh = ss2777.plotTaskjointTorqeh(); title('ss2777');
ss2777.plotAddTrialMark(axh.Children(1));
ylim([-50 50])

axh = ss2778.plotTaskjointTorqeh(); title('ss2778');
ss2778.plotAddTrialMark(axh.Children(1));
ylim([-5 5])

ss2779 = SessionScan(2779)
axh = ss2779.plotTaskjointTorqeh(); title('ss2779');
ss2779.plotAddTrialMark(axh.Children(1));
ylim([-5 5])

%% position accuracy
%ss2780 = SessionScan(2780);
ss2780.plotStepPertResponse_raw(1, [0.5 0.5 0.5])
xlim([-0.5 2.5]);
pos_csd = ss2780.getPertPos(500);
pos_cmd = [ss2780.trials.pert_dx0];
subplot(1,2,1);
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd(2:end), pos_csd', '*'); hold on;
plot(pos_cmd(2:end), pos_cmd(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2);
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd(2:end), pos_csd' - pos_cmd(2:end));
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
%xlim([-0.08 0.08]); ylim([-0.006 0.006]);

%ss2783 = SessionScan(2783);
ss2783.plotStepPertResponse_raw(1, [0.5 0.5 0.5])
xlim([-0.5 2.5]);
pos_csd = ss2783.getPertPos(500);
pos_cmd = [ss2783.trials.pert_dx0];
subplot(1,2,1);
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd(2:end), pos_csd', '*'); hold on;
plot(pos_cmd(2:end), pos_cmd(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2);
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd(2:end), pos_csd' - pos_cmd(2:end));
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
xlim([-0.08 0.08]); ylim([-0.006 0.006]);

%ss2785 = SessionScan(2785);
ss2785.plotStepPertResponse_raw(1, [0.5 0.5 0.5])
xlim([-0.5 2.5]);
pos_csd = ss2785.getPertPos(500);
pos_cmd = [ss2785.trials.pert_dx0];
subplot(1,2,1);
plot(pos_cmd(2:end), pos_csd, '*'); hold on;
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd(2:end), pos_cmd(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2);
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd(2:end), pos_csd' - pos_cmd(2:end));
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
%xlim([-0.08 0.08]); ylim([-0.006 0.006]);
%% position accuracy test, not spring test.
clf; clc; close all; clear;
ss2780 = SessionScan(2780);
ss2783 = SessionScan(2783);
ss2785 = SessionScan(2785);
% Give a plot in one canvas
pos_csd300 = ss2780.getPertPos(500);
pos_cmd300 = [ss2780.trials.pert_dx0];
pos_csd2500 = ss2783.getPertPos(500);
pos_cmd2500 = [ss2783.trials.pert_dx0];
pos_csd5000 = ss2785.getPertPos(500);
pos_cmd5000 = [ss2785.trials.pert_dx0];
%fce_csd5000 = ss2785.getPertFce(500);
subplot(1,2,1); hold on;
plot(pos_cmd300(2:end), pos_csd300, 'r*'); 
plot(pos_cmd2500(2:end), pos_csd2500, 'g*'); 
plot(pos_cmd5000(2:end), pos_csd5000, 'b*'); 
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd300(2:end), pos_cmd300(2:end), 'Color', [0.5 0.5 0.5]);
% xlabel('command position x0 change (m)')
% ylabel('measured position x0 change (m)')
% title('position command and yield position');
xlabel('\Delta x_0 (m)')
ylabel('\Delta x (m)')
title('position control when no load');
legend('K=300', 'K=2500', 'K=5000', 'refline');
grid on;
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2); hold on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd300(2:end), pos_csd300' - pos_cmd300(2:end), 'r');
stem(pos_cmd2500(2:end), pos_csd2500' - pos_cmd2500(2:end), 'g');
stem(pos_cmd5000(2:end), pos_csd5000' - pos_cmd5000(2:end), 'b');
legend('K300', 'K2500', 'K5000');
ylim([-0.006 0.006]);
% xlabel('command position x0 change (m)')
% ylabel('the position error (m)')
ylabel('\Delta x - \Delta x_0 (m)')
xlabel('\Delta x_0 (m)')
title('position error');
grid on;

figure();
axh1 = subplot(3,1,1);
axh2 = subplot(3,1,2);
axh3 = subplot(3,1,3);
ss2780.plotStepPertResponse_raw(axh1, [0.5 0.5 0.5]);
ss2783.plotStepPertResponse_raw(axh2, [0.5 0.5 0.5]);
ss2785.plotStepPertResponse_raw(axh3, [0.5 0.5 0.5]);

figure();
axh1 = subplot(3,1,1);
axh2 = subplot(3,1,2);
axh3 = subplot(3,1,3);
ss2780.plotStepPertResponse_rawV(axh1, [0.5 0.5 0.5]);
ss2783.plotStepPertResponse_rawV(axh2, [0.5 0.5 0.5]);
ss2785.plotStepPertResponse_rawV(axh3, [0.5 0.5 0.5]);

%% Spring test displacement perturbation
%ss2786 = SessionScan(2786);
ss2786.plotStepPertResponse_raw(1, [0.5 0.5 0.5])
ss2786.plotStepPertResponse_rawF(1, [0.5 0.5 0.5])
xlim([-0.5 2.5]);
pos_csd = ss2786.getPertPos(500);
fce_csd = ss2786.getPertFce(500);
pos_cmd = [ss2786.trials.pert_dx0];
subplot(1,2,1);
plot(pos_cmd(2:end), pos_csd, '*'); hold on;
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd(2:end), pos_cmd(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2);
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd(2:end), pos_csd' - pos_cmd(2:end));
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
%xlim([-0.08 0.08]); ylim([-0.006 0.006]);
figure();
subplot(1,2,1); plot(pos_cmd(2:end)', fce_csd);
subplot(1,2,2); plot(pos_csd, fce_csd);

% The setted stiffness of 2500
%%
ss2787.plotStepPertResponse_raw(1, [0.5 0.5 0.5])
ss2787.plotStepPertResponse_rawF(1, [0.5 0.5 0.5])
xlim([-0.5 2.5]);
pos_csd = ss2787.getPertPos(500);
fce_csd = ss2787.getPertFce(500);
pos_cmd = [ss2787.trials.pert_dx0];
subplot(1,2,1);
plot(pos_cmd(2:end), pos_csd, '*'); hold on;
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd(2:end), pos_cmd(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2);
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd(2:end), pos_csd' - pos_cmd(2:end));
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
%xlim([-0.08 0.08]); ylim([-0.006 0.006]);
figure();
axh1 = subplot(1,2,1); plot(pos_cmd(2:end)', fce_csd, '*');
axh2 = subplot(1,2,2); plot(pos_csd, fce_csd, '*');
linkaxes([axh1, axh2], 'xy');

%
%% ss2788 = SessionScan(2788)
ss2788.plotStepPertResponse_raw(1, [0.5 0.5 0.5])
ss2788.plotStepPertResponse_rawF(1, [0.5 0.5 0.5])
xlim([-0.5 2.5]);
pos_csd = ss2788.getPertPos(500);
fce_csd = ss2788.getPertFce(500);
pos_cmd = [ss2788.trials.pert_dx0];
subplot(1,2,1);
plot(pos_cmd(2:end), pos_csd, '*'); hold on;
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd(2:end), pos_cmd(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2);
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd(2:end), pos_csd' - pos_cmd(2:end));
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
%xlim([-0.08 0.08]); ylim([-0.006 0.006]);
figure();
axh1 = subplot(1,2,1); plot(pos_cmd(2:end)', fce_csd, '*');
axh2 = subplot(1,2,2); plot(pos_csd, fce_csd, '*');
linkaxes([axh1, axh2], 'xy');

%%  plot the three spring tests in one figure 
ss2786 = SessionScan(2786);
ss2787 = SessionScan(2787);
ss2788 = SessionScan(2788);
%
pos_csd300 = ss2786.getPertPos(500);
pos_cmd300 = [ss2786.trials.pert_dx0];
fce_csd300 = -ss2786.getPertFce(500);
pos_csd2500 = ss2787.getPertPos(500);
pos_cmd2500 = [ss2787.trials.pert_dx0];
fce_csd2500 = -ss2787.getPertFce(500);
pos_csd5000 = ss2788.getPertPos(500);
pos_cmd5000 = [ss2788.trials.pert_dx0];
fce_csd5000 = -ss2788.getPertFce(500);
subplot(1,2,1); hold on;
grid on;
%plot(pos_cmd100(2:end), pos_csd100, 'c*'); 
plot(pos_cmd300(2:end), pos_csd300, 'r*'); 
plot(pos_cmd2500(2:end), pos_csd2500, 'g*'); 
plot(pos_cmd5000(2:end), pos_csd5000, 'b*'); 
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd300(2:end), pos_cmd300(2:end), 'Color', [0.5 0.5 0.5]);
% legend('K300', 'K2500', 'K5000', 'refline');
% xlabel('command position x0 change (m)')
% ylabel('measured position x0 change (m)')
% title('position command and yield position');
legend('K=300', 'K=2500', 'K=5000', 'refline');
xlabel('\Delta x_0 (m)')
ylabel('\Delta x (m)')
title('\Delta x is a porpotion of \Delta x_0');
xlim([-0.08 0.08]); ylim([-0.04 0.04]);
subplot(1,2,2); hold on;
grid on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
%stem(pos_cmd100(2:end), pos_csd100' - pos_cmd100(2:end), 'c');
% stem(pos_cmd300(2:end), pos_csd300' - pos_cmd300(2:end), 'r');
% stem(pos_cmd2500(2:end), pos_csd2500' - pos_cmd2500(2:end), 'g');
% stem(pos_cmd5000(2:end), pos_csd5000' - pos_cmd5000(2:end), 'b');
plot(pos_cmd300(2:end), pos_csd300' - pos_cmd300(2:end), 'r*');
plot(pos_cmd2500(2:end), pos_csd2500' - pos_cmd2500(2:end), 'g*');
plot(pos_cmd5000(2:end), pos_csd5000' - pos_cmd5000(2:end), 'b*');
legend('K300', 'K2500', 'K5000');
% xlabel('command position x0 change (m)')
% ylabel('the position error (m)')
% title('position error');
xlabel('\Delta x_0 (m)')
ylabel('\Delta x - \Delta x_0 (m)')
title('the displacement taken by the robot');
%% plot the measured impedance... 
ss2789 = SessionScan(2789);
pos_csd100 = ss2789.getPertPos(500);
pos_cmd100 = [ss2789.trials.pert_dx0];
fce_csd100 = -ss2789.getPertFce(500);
% subplot(1,3,3); hold on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
figure(); hold on;
plot(pos_csd100, fce_csd100, 'c*');
plot(pos_csd300, fce_csd300, 'r*');
plot(pos_csd2500, fce_csd2500, 'g*');
plot(pos_csd5000, fce_csd5000, 'b*');
grid on;
axh = refline;
legend('K=100', 'K=300', 'K=2500', 'K=5000');
% xlabel('censored position x0 change (m)')
% ylabel('censored force (N)')
% title('stiffness estimation');
xlabel('\Delta x (m)')
ylabel('\Delta F_{csd} (N)')
title('\Delta F = K_s \Delta x');
% read k from the axh
k = zeros(size(axh)); b = zeros(size(axh));
for i = 1:4
    k(i) = diff(axh(i).YData)/diff(axh(i).XData);
    b(i) = axh(i).YData(1) - k(i)*axh(i).XData(1);
end


%% 
% use force perturbation again
ss2790 = SessionScan(2790)
ss2790.plotStepPertResponse_raw(1, [0.5 0.5 0.5])
ss2790.plotStepPertResponse_rawF(1, [0.5 0.5 0.5])
xlim([-0.5 2.5]);
pos_csd = ss2790.getPertPos(500);
fce_csd = -ss2790.getPertFce(500);
fce_cmd = [ss2790.trials.pert_f];
subplot(1,2,1);
plot(fce_cmd(2:end), fce_csd, '*'); hold on;
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(fce_cmd(2:end), fce_cmd(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2);
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(fce_cmd(2:end), fce_csd' - fce_cmd(2:end));
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
%xlim([-0.08 0.08]); ylim([-0.006 0.006]);
figure();
axh1 = subplot(1,2,1); plot(pos_csd, fce_cmd(2:end)', '*');
axh2 = subplot(1,2,2); plot(pos_csd, fce_csd, '*');
linkaxes([axh1, axh2], 'xy');

%% 
clf; clc; clear;
ss2790 = SessionScan(2790);
ss2791 = SessionScan(2791);
ss2795 = SessionScan(2795);
% plot all spring test in this figure
pos_csd300 = ss2790.getPertPos(500);
fce_cmd300 = [ss2790.trials.pert_f];
fce_csd300 = -ss2790.getPertFce(500);
pos_csd2500 = ss2791.getPertPos(500);
fce_cmd2500 = [ss2791.trials.pert_f];
fce_csd2500 = -ss2791.getPertFce(500);
pos_csd5000 = ss2795.getPertPos(500);
fce_cmd5000 = [ss2795.trials.pert_f];
fce_csd5000 = -ss2795.getPertFce(500);
subplot(1,2,1); hold on;
grid on;
%plot(pos_cmd100(2:end), pos_csd100, 'c*'); 
plot(fce_cmd300(2:end), fce_csd300, 'r*'); 
plot(fce_cmd2500(2:end), fce_csd2500, 'g*'); 
plot(fce_cmd5000(2:end), fce_csd5000, 'b*'); 
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(fce_cmd300(2:end), fce_cmd300(2:end), 'Color', [0.5 0.5 0.5]);
%legend('K300', 'K2500', 'refline');
% legend('K300', 'K2500', 'K5000', 'refline');
% xlabel('command force change (m)')
% ylabel('measured force change (m)')
% title('force command and yield force');
legend('K=300', 'K=2500', 'K=5000', 'refline');
xlabel('F_{pert} (N)')
ylabel('F_{csd} (N)')
title('F_{csd} is a porpotion of F_{pert}');
%xlim([-0.08 0.08]); ylim([-0.04 0.04]);
subplot(1,2,2); hold on;
grid on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
%stem(pos_cmd100(2:end), pos_csd100' - pos_cmd100(2:end), 'c');
% stem(fce_cmd300(2:end), fce_csd300' - fce_cmd300(2:end), 'r');
% stem(fce_cmd2500(2:end), fce_csd2500' - fce_cmd2500(2:end), 'g');
% stem(fce_cmd5000(2:end), fce_csd5000' - fce_cmd5000(2:end), 'b');
plot(fce_cmd300(2:end), fce_csd300' - fce_cmd300(2:end), 'r*');
plot(fce_cmd2500(2:end), fce_csd2500' - fce_cmd2500(2:end), 'g*');
plot(fce_cmd5000(2:end), fce_csd5000' - fce_cmd5000(2:end), 'b*');
% legend('K300', 'K2500', 'K5000');
% xlabel('command force (N)')
% ylabel('the force error (N)')
% title('force error');
legend('K=300', 'K=2500', 'K=5000');
xlabel('F_{pert} (N)')
ylabel('F_{csd} - F_{pert} (N)')
title('the force take by the robot');

%% plot the measured impedance...
ss2794 = SessionScan(2794);
pos_csd100 = ss2794.getPertPos(500);
fce_cmd100 = [ss2794.trials.pert_f];
fce_csd100 = -ss2794.getPertFce(500);
% subplot(1,3,3); hold on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
figure(); hold on;
plot(pos_csd100, fce_csd100, 'c*');
plot(pos_csd300, fce_csd300, 'r*');
plot(pos_csd2500, fce_csd2500, 'g*');
plot(pos_csd5000, fce_csd5000, 'b*');
grid on;
axh = refline;
legend('K100', 'K300', 'K2500', 'K5000');
% xlabel('censored position x0 change (m)')
% ylabel('censored force (N)')
% title('stiffness estimation');
% read k from the axh
xlabel('\Delta x (m)')
ylabel('\Delta F_{csd} (N)')
title('{\Delta F}=K_s{\Delta x}');
k = zeros(size(axh)); b = zeros(size(axh));
for i = 1:length(axh)
    k(i) = diff(axh(i).YData)/diff(axh(i).XData);
    b(i) = axh(i).YData(1) - k(i)*axh(i).XData(1);
end

%% 1. plot large scales 
ss_list = [2786 2787 2788 2789 2790 2791 2794 2795]; % may inaccurate as the force at 0.
ss_list = [ 2790 2791 2794 2795]; % may inaccurate as the force at 0.
% 2722 2729 2734
for session_i = 1:length(ss_list)
     eval(['ss' num2str(ss_list(session_i)) ' = SessionScan(' num2str(ss_list(session_i)) ');']);
end

%% Do the comparition measurement;
% 1. For the force perturbation
% ss2794 = SessionScan(2794);
% ss2790 = SessionScan(2790);
% ss2791 = SessionScan(2791);
% ss2795 = SessionScan(2795);

pos_csd100f = ss2794.getPertPos(500);
fce_csd100f = -ss2794.getPertFce(100);
fce_cmd100f = ss2794.getPertFce_cmd(); %fce_cmd100f = fce_cmd100f(2:end)';
pos_csd300f = ss2790.getPertPos(500);
fce_csd300f = -ss2790.getPertFce(500);
fce_cmd300f = ss2790.getPertFce_cmd(); %fce_cmd300f = fce_cmd300f(2:end)';
pos_csd2500f = ss2791.getPertPos(500);
fce_csd2500f = -ss2791.getPertFce(500);
fce_cmd2500f = ss2791.getPertFce_cmd(); %fce_cmd2500f = fce_cmd2500f(2:end)';
pos_csd5000f = ss2795.getPertPos(500);
fce_csd5000f = -ss2795.getPertFce(500);
fce_cmd5000f = ss2795.getPertFce_cmd(); %fce_cmd5000f = fce_cmd5000f(2:end)';

save('forcePert.mat', 'pos*', 'fce*');

pos_csd100p = ss2789.getPertPos(500);
fce_csd100p = -ss2789.getPertFce(500);
pos_csd300p = ss2786.getPertPos(500);
fce_csd300p = -ss2786.getPertFce(500);
pos_csd2500p = ss2787.getPertPos(500);
fce_csd2500p = -ss2787.getPertFce(500);
pos_csd5000p = ss2788.getPertPos(500);
fce_csd5000p = -ss2788.getPertFce(500);
pos_cmd100p = [ss2789.trials.pert_dx0];     pos_cmd100p = pos_cmd100p(2:end)';
pos_cmd300p = [ss2786.trials.pert_dx0];     pos_cmd300p = pos_cmd300p(2:end)';
pos_cmd2500p = [ss2787.trials.pert_dx0];    pos_cmd2500p = pos_cmd2500p(2:end)';
pos_cmd5000p = [ss2788.trials.pert_dx0];    pos_cmd5000p = pos_cmd5000p(2:end)';

% pert_conditions: K_r, F_pert, \Delta F_csd, \Delta x..... K_s_est 
pertF_conditions = [100*ones(size(fce_cmd100f)) fce_cmd100f fce_csd100f pos_csd100f;
    300*ones(size(fce_cmd300f)) fce_cmd300f fce_csd300f pos_csd300f;
    2500*ones(size(fce_cmd2500f)) fce_cmd2500f fce_csd2500f pos_csd2500f;
    5000*ones(size(fce_cmd5000f)) fce_cmd5000f fce_csd5000f pos_csd5000f;];

pertF_conditions(:,5) = [pertF_conditions(:,3)./pertF_conditions(:,4)];
pertF_conditions(:,6) = [pertF_conditions(:,2)./pertF_conditions(:,4) - pertF_conditions(:,1)];

% pert_conditions: K_r, \Delta_x_0, \Delta F_csd, \Delta x..... K_s_est 
pertp_conditions = [100*ones(size(pos_cmd100p)) pos_cmd100p fce_csd100p pos_csd100p;
    300*ones(size(pos_cmd300p)) pos_cmd300p fce_csd300p pos_csd300p;
    2500*ones(size(pos_cmd2500p)) pos_cmd2500p fce_csd2500p pos_csd2500p;
    5000*ones(size(pos_cmd5000p)) pos_cmd5000p fce_csd5000p pos_csd5000p;];
pertp_conditions(:,5) = [pertp_conditions(:,3)./pertp_conditions(:,4)];



fce_lim = [-8 8];
figure()
idx = pertF_conditions(:,1) == 100 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2));
axh(1) = subplot(4,2,1); plot(pertF_conditions(idx,4), pertF_conditions(idx,3),'c*');
title('force control');
idx = pertF_conditions(:,1) == 300 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2));
axh(2) = subplot(4,2,3); plot(pertF_conditions(idx,4), pertF_conditions(idx,3),'r*');
idx = pertF_conditions(:,1) ==2500 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2));
axh(3) = subplot(4,2,5); plot(pertF_conditions(idx,4), pertF_conditions(idx,3),'g*');
idx = pertF_conditions(:,1) ==5000 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2));
axh(4) = subplot(4,2,7); plot(pertF_conditions(idx,4), pertF_conditions(idx,3),'b*');
linkaxes(axh, 'xy');
idx = pertp_conditions(:,1) == 100 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2));
axh(1) = subplot(4,2,2); plot(pertp_conditions(idx,4), pertp_conditions(idx,3),'c*');
title('position control');
idx = pertp_conditions(:,1) == 300 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2));
axh(2) = subplot(4,2,4); plot(pertp_conditions(idx,4), pertp_conditions(idx,3),'r*');
idx = pertp_conditions(:,1) ==2500 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2));
axh(3) = subplot(4,2,6); plot(pertp_conditions(idx,4), pertp_conditions(idx,3),'g*');
idx = pertp_conditions(:,1) ==5000 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2));
axh(4) = subplot(4,2,8); plot(pertp_conditions(idx,4), pertp_conditions(idx,3),'b*');
linkaxes(axh, 'xy');
xlabel('\Delta x'); ylabel('\Delta F_{csd}');

% plot out of errorbar FORCE PERT
%%
fce_lim = [-8 8];
Ks_mean = 157.6*2; Ks_std = 8.8*2;
Ks_lo = Ks_mean - Ks_std;
Ks_hi = Ks_mean + Ks_std;

stiffness100F = pertF_conditions(pertF_conditions(:,1) == 100 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2)), 5);
stiffness300F = pertF_conditions(pertF_conditions(:,1) == 300 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2)), 5);
stiffness2500F = pertF_conditions(pertF_conditions(:,1) ==2500 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2)), 5);
stiffness5000F = pertF_conditions(pertF_conditions(:,1) ==5000 & (pertF_conditions(:,3) >fce_lim(1)) & (pertF_conditions(:,3) < fce_lim(2)), 5);

stiffness100P = pertp_conditions(pertp_conditions(:,1) == 100 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2)), 5);
stiffness300P = pertp_conditions(pertp_conditions(:,1) == 300 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2)), 5);
stiffness2500P = pertp_conditions(pertp_conditions(:,1) == 2500 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2)), 5);
stiffness5000P = pertp_conditions(pertp_conditions(:,1) == 5000 & (pertp_conditions(:,3) >fce_lim(1)) & (pertp_conditions(:,3) < fce_lim(2)), 5);

%% plot the estimation bar
axh = figure();
subplot(1,2,1); hold on;
series = [1 2 3 4];
X = [ones(size(stiffness100F))*1; ...
    ones(size(stiffness300F))*2; ...
    ones(size(stiffness2500F))*3; ...
    ones(size(stiffness5000F))*4]';
y = [stiffness100F; stiffness300F; stiffness2500F; stiffness5000F]';

stiffmean = [mean(stiffness100F), mean(stiffness300F), mean(stiffness2500F), mean(stiffness5000F)];
errbar_half = [std(stiffness100F), std(stiffness300F), std(stiffness2500F), std(stiffness5000F)]; % 1std
bar(series,stiffmean)                
er = errorbar(series,stiffmean,errbar_half,errbar_half);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
plot(X(X==1)+0.1, y(X==1), 'c*');
plot(X(X==2)+0.1, y(X==2), 'r*');
plot(X(X==3)+0.1, y(X==3), 'g*');
plot(X(X==4)+0.1, y(X==4), 'b*');
plot([0,5], [Ks_mean Ks_mean], 'Color', [0.5 0.5 0.5]);
plot([0,5], [Ks_lo Ks_lo], '--', 'Color', [0.5 0.5 0.5]);
plot([0,5], [Ks_hi Ks_hi], '--', 'Color', [0.5 0.5 0.5]);

ylim([200 400]); grid on;
title('force perturbation');
ylabel('\DeltaF/\Deltax');
ylabel('\DeltaF/\Deltax (N/m)');
xlabel('robot stiffness (N/m)');
xticks([1 2 3 4]); xticklabels({'100', '300', '2500', '5000'});
figure(axh);
subplot(1,2,2); hold on;
% plot out of errorbar POS PERT
series = [1 2 3 4];
X = [ones(size(stiffness100P))*1; ...
    ones(size(stiffness300P))*2; ...
    ones(size(stiffness2500P))*3; ...
    ones(size(stiffness5000P))*4]';
y = [stiffness100P; stiffness300P; stiffness2500P; stiffness5000P]';
stiffmean = [mean(stiffness100P), mean(stiffness300P), mean(stiffness2500P), mean(stiffness5000P)];
errbar_half = [std(stiffness100P), std(stiffness300P), std(stiffness2500P), std(stiffness5000P)]; % 1std
bar(series,stiffmean)                
er = errorbar(series,stiffmean,errbar_half,errbar_half);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
plot(X(X==1)+0.1, y(X==1), 'c*');
plot(X(X==2)+0.1, y(X==2), 'r*');
plot(X(X==3)+0.1, y(X==3), 'g*');
plot(X(X==4)+0.1, y(X==4), 'b*');
plot([0,5], [Ks_mean Ks_mean], 'Color', [0.5 0.5 0.5]);
plot([0,5], [Ks_lo Ks_lo], '--', 'Color', [0.5 0.5 0.5]);
plot([0,5], [Ks_hi Ks_hi], '--', 'Color', [0.5 0.5 0.5]);
ylim([200 400]); grid on;
title('position perturbation');
ylabel('\DeltaF/\Deltax (N/m)');
xlabel('robot stiffness (N/m)');
xticks([1 2 3 4]); xticklabels({'100', '300', '2500', '5000'});