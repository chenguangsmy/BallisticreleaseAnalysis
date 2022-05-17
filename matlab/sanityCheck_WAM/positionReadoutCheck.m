% In this scripts, check whether the position readout is accurate. 

% Position histerisis
% I clamped the robot, and do a positive-negative-positive command force,
% I hope the force-position relationship will shows a histerisis curve. 
sstmp = SessionScan(4005);

lw = 1;
t_range = [3767.5, 3771.5];
%plot(sstmp.data.t, sstmp.data.x(2,:));
t_idx = sstmp.data.t>t_range(1) & sstmp.data.t<t_range(2);

fh1 = figure('unit', 'inch', 'position', [0 0 6 3]);  % plot the command force, the censored force and the censored positions
axh(1) = subplot(2,1,1);  hold on;
plot(sstmp.data.t(t_idx), sstmp.data.Fp(2,t_idx), 'linewidth', lw);
plot(sstmp.data.t(t_idx), -sstmp.data.f(2,t_idx), 'linewidth', lw);
grid on;
xlabel('time (s)'); ylabel('Force (N)');
legend('command', 'censored');
title('command and censored force');
axh(2) = subplot(2,1,2); hold on;
plot(sstmp.data.t(t_idx), sstmp.data.x(2,t_idx), 'linewidth', lw);
plot(sstmp.data.t(t_idx), sstmp.data.opty(t_idx), 'linewidth', lw);
grid on;
legend('WAM', 'OPTOTRAK');
title('robot and motion capture position');
xlabel('time (s)'); ylabel('Position (m)');
linkaxes(axh, 'x'); 
xlim(t_range);

saveas(fh1,'/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_WAM/WAMpositionSanityCheck/forcePositionTimeSequence.png');

% fh2 = figure(2); % plot 1) the command force-censored force histerisis and 
%                  %      2) the command force-censored position histerisis
% axh(1) = subplot(1,3,1); hold on; grid on;
% plot(sstmp.data.Fp(2,t_idx), -sstmp.data.f(2,t_idx), 'Marker', '.');
% axh(2) = subplot(1,3,2); hold on; grid on;
% x_offset = sstmp.data.x(2,find(t_idx));
% x_shift = sstmp.data.x(2,t_idx)  - x_offset(1);
% plot(sstmp.data.Fp(2,t_idx), x_shift, 'Marker', '.');
% axh(3) = subplot(1,3,3); hold on; grid on;
% %lot(-sstmp.data.f(2,t_idx), x_shift, 'Marker', '.');

% figure3, plot the censored Force ~ x error historisis 
fh3 = figure('unit', 'inch', 'position', [0 0 4 4]);
axh(1) = subplot('position', [ 0.15 0.1 0.7 0.7]);
fce = sstmp.data.f(2,t_idx); 
x_e = sstmp.data.opty(t_idx) - sstmp.data.x(2,t_idx);
plot(fce, x_e, 'Marker', '.');
xlim([-30 30]); ylim([-4 4]*1e-3);
grid on; 
xlabel('force Censored (N)', 'fontSize', 12)
ylabel('x_{WAM} - x_{OPT} (m)', 'fontSize', 12)
sgtitle('hysteresis property of WAM position readout', 'fontSize', 14);
saveas(fh3,'/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_WAM/WAMpositionSanityCheck/positionHysteresis.png');