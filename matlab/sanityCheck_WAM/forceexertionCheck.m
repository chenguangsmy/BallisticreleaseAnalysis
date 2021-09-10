% force exertion check
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
% ss2736
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
ss2750 = SessionScan(2750);
ss2750.plotStepPertResponse_rawFce(1, [0.5 0.5 0.5]); % 10Hz filter
fce_cmd = ss2750.getPertFce_cmd(); 
fce_csd = -ss2750.getPertFce(200);
pos_csd = ss2750.getPertPos(5);
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
legend('datapoints', 'commands');%, 'pos points fit', 'neg points fit');
grid on;
%text(fce_cmd_pos(2), fce_csd_pos(2), ['slope:' num2str(P1)]);
%text(fce_cmd_neg(2), fce_csd_neg(2), ['slope:' num2str(P2)]);
xlabel('force command'); ylabel('force measured');
title('command force vs censored');
subplot(1,2,2)
stem(fce_cmd, fce_csd-fce_cmd); grid on;
xlabel('force command'); ylabel('force difference');
title('force error');
figure();
plot(fce_cmd, pos_csd, '*', 'MarkerSize', 5); hold on;
%plot(fce_cmd, fce_cmd, 'Color', [0.5 0.5 0.5]);
grid on;
xlabel('force command'); ylabel('position measured');
title('command force vs position');

%% ss2751: clamping the robot, with change the toolplate length 50 41cm, turned on the update jaccobian
ss2751 = SessionScan(2751);
ss2751.plotStepPertResponse_rawFce(1, [0.5 0.0 0.5]); % 10Hz filter
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
ss2770.plotStepPertResponse_raw(1, [0 0 0]); 
xlim([-2 2]); grid on;
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
%%
% Give a plot in one canvas
pos_csd300 = ss2780.getPertPos(500);
pos_cmd300 = [ss2780.trials.pert_dx0];
pos_csd2500 = ss2783.getPertPos(500);
pos_cmd2500 = [ss2783.trials.pert_dx0];
pos_csd5000 = ss2785.getPertPos(500);
pos_cmd5000 = [ss2785.trials.pert_dx0];
subplot(1,2,1); hold on;
plot(pos_cmd300(2:end), pos_csd300, 'r*'); 
plot(pos_cmd2500(2:end), pos_csd2500, 'g*'); 
plot(pos_cmd5000(2:end), pos_csd5000, 'b*'); 
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd300(2:end), pos_cmd300(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,2,2); hold on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd300(2:end), pos_csd300' - pos_cmd300(2:end), 'r');
stem(pos_cmd2500(2:end), pos_csd2500' - pos_cmd2500(2:end), 'g');
stem(pos_cmd5000(2:end), pos_csd5000' - pos_cmd5000(2:end), 'b');
legend('K300', 'K2500', 'K5000');
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');

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
pos_csd100 = ss2789.getPertPos(500);
pos_cmd100 = [ss2789.trials.pert_dx0];
fce_csd100 = ss2789.getPertFce(500);
pos_csd300 = ss2786.getPertPos(500);
pos_cmd300 = [ss2786.trials.pert_dx0];
fce_csd300 = ss2786.getPertFce(500);
pos_csd2500 = ss2787.getPertPos(500);
pos_cmd2500 = [ss2787.trials.pert_dx0];
fce_csd2500 = ss2787.getPertFce(500);
pos_csd5000 = ss2788.getPertPos(500);
pos_cmd5000 = [ss2788.trials.pert_dx0];
fce_csd5000 = ss2788.getPertFce(500);
subplot(1,3,1); hold on;
plot(pos_cmd100(2:end), pos_csd100, 'c*'); 
plot(pos_cmd300(2:end), pos_csd300, 'r*'); 
plot(pos_cmd2500(2:end), pos_csd2500, 'g*'); 
plot(pos_cmd5000(2:end), pos_csd5000, 'b*'); 
%plot(pos_cmd(2:end), pos_csd'-0.482, '*'); hold on;
plot(pos_cmd300(2:end), pos_cmd300(2:end), 'Color', [0.5 0.5 0.5]);
xlabel('command position x0 change (m)')
ylabel('measured position x0 change (m)')
title('position command and yield position');
%xlim([-0.08 0.08]); ylim([-0.08 0.08]);
subplot(1,3,2); hold on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
stem(pos_cmd100(2:end), pos_csd100' - pos_cmd100(2:end), 'c');
stem(pos_cmd300(2:end), pos_csd300' - pos_cmd300(2:end), 'r');
stem(pos_cmd2500(2:end), pos_csd2500' - pos_cmd2500(2:end), 'g');
stem(pos_cmd5000(2:end), pos_csd5000' - pos_cmd5000(2:end), 'b');
legend('K300', 'K2500', 'K5000');
xlabel('command position x0 change (m)')
ylabel('the position error (m)')
title('position error');
subplot(1,3,3); hold on;
%stem(pos_cmd(2:end), pos_csd'-0.482 - pos_cmd(2:end));
plot(pos_csd100, fce_csd100, 'c*');
plot(pos_csd300, fce_csd300, 'r*');
plot(pos_csd2500, fce_csd2500, 'g*');
plot(pos_csd5000, fce_csd5000, 'b*');
legend('K300', 'K2500', 'K5000');
xlabel('censored position x0 change (m)')
ylabel('censored force (N)')
title('position error');
