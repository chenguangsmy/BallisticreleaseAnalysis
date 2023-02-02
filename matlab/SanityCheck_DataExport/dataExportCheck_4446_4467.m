% 
% ballistic release data check 
blp = ballisticReleaseTaksPlots()

% 1. exception on subject 1 
% figure(); 
% trial_sel = 2; 
% trials_other = setdiff(1:9,2);
% axh(1) = subplot(2,2,1); hold on; 
% for trial_i = trials_other
%     plot(blp.data{3,2,1,1,trial_i,1}.t, blp.data{3,2,1,1,trial_i,1}.ox(1,:,1), 'Color',[0.3 0.3 0.3]);
% end
% plot(blp.data{3,2,1,1,trial_sel,1}.t, blp.data{3,2,1,1,trial_sel,1}.ox(1,:,1), 'Color',[1 0.0 0.0], 'LineWidth', 1);
% grid on;
% title('displacement (x)'); 
% ylabel('x (m)');
% xlabel('time (s)');
% 
% trial_sel = 2; 
% trials_other = setdiff(1:9,2);
% axh(1) = subplot(2,2,2); hold on; 
% for trial_i = trials_other
%     plot(blp.data{3,2,1,1,trial_i,1}.t, blp.data{3,2,1,1,trial_i,1}.ox(1,:,2), 'Color',[0.3 0.3 0.3]);
% end
% plot(blp.data{3,2,1,1,trial_sel,1}.t, blp.data{3,2,1,1,trial_sel,1}.ox(1,:,2), 'Color',[1 0.0 0.0], 'LineWidth', 1);
% grid on;
% title('displacement (y)'); 
% ylabel('x (m)');
% xlabel('time (s)');
% 
% 
% 
% axh(2) = subplot(2,2,3); hold on; 
% for trial_i = trials_other
%     plot(blp.data{3,2,1,1,trial_i,1}.t, blp.data{3,2,1,1,trial_i,1}.f(1,:), 'Color',[0.3 0.3 0.3]);
% end
% plot(blp.data{3,2,1,1,trial_sel,1}.t, blp.data{3,2,1,1,trial_sel,1}.f(1,:), 'Color',[1 0.0 0.0], 'LineWidth', 1);
% grid on;
% title('force (x)'); 
% ylabel('f (N)');
% xlabel('time (s)');
% 
% axh(2) = subplot(2,2,4); hold on; 
% for trial_i = trials_other
%     plot(blp.data{3,2,1,1,trial_i,1}.t, blp.data{3,2,1,1,trial_i,1}.f(2,:), 'Color',[0.3 0.3 0.3]);
% end
% plot(blp.data{3,2,1,1,trial_sel,1}.t, blp.data{3,2,1,1,trial_sel,1}.f(2,:), 'Color',[1 0.0 0.0], 'LineWidth', 1);
% grid on;
% title('force (y)'); 
% ylabel('f (N)');
% xlabel('time (s)');
% linkaxes(axh, 'x');
% xlim([-0.5 1.0]);

figure(); 
trial_sel = 2; 
trials_other = setdiff(1:9,2);
axh(1) = subplot(2,1,1); hold on; 
for trial_i = trials_other
    plot(blp.data{3,2,1,1,trial_i,1}.t, blp.data{3,2,1,1,trial_i,1}.ox(1,:,1), 'Color',[0.3 0.3 0.3]);
end
plot(blp.data{3,2,1,1,trial_sel,1}.t, blp.data{3,2,1,1,trial_sel,1}.ox(1,:,1), 'Color',[1 0.0 0.0], 'LineWidth', 1);
grid on;
title('displacement (x)'); 
ylabel('x (m)');
xlabel('time (s)');


axh(2) = subplot(2,1,2); hold on; 
for trial_i = trials_other
    plot(blp.data{3,2,1,1,trial_i,1}.t, blp.data{3,2,1,1,trial_i,1}.f(1,:), 'Color',[0.3 0.3 0.3]);
end
plot(blp.data{3,2,1,1,trial_sel,1}.t, blp.data{3,2,1,1,trial_sel,1}.f(1,:), 'Color',[1 0.0 0.0], 'LineWidth', 1);
grid on;
title('force (x)'); 
ylabel('f (N)');
xlabel('time (s)');

linkaxes(axh, 'x');
xlim([-0.5 1.0]);


%% 
figure(); 
trial_sel = 6; 
trials_other = setdiff(1:9,2);
axh(1) = subplot(2,1,1); hold on; 
for trial_i = trials_other
    plot(blp.data{1,2,2,1,trial_i,1}.t, blp.data{1,2,2,1,trial_i,1}.ox(1,:,1), 'Color',[0.3 0.3 0.3]);
end
plot(blp.data{1,2,2,1,trial_sel,1}.t, blp.data{1,2,2,1,trial_sel,1}.ox(1,:,1), 'Color',[1 0.0 0.0], 'LineWidth', 1);
grid on;
title('displacement (x)'); 
ylabel('x (m)');
xlabel('time (s)');


axh(2) = subplot(2,1,2); hold on; 
for trial_i = trials_other
    plot(blp.data{1,2,2,1,trial_i,1}.t, blp.data{1,2,2,1,trial_i,1}.f(1,:), 'Color',[0.3 0.3 0.3]);
end
plot(blp.data{1,2,2,1,trial_sel,1}.t, blp.data{1,2,2,1,trial_sel,1}.f(1,:), 'Color',[1 0.0 0.0], 'LineWidth', 1);
grid on;
title('force (x)'); 
ylabel('f (N)');
xlabel('time (s)');

linkaxes(axh, 'x');
xlim([-0.5 1.0]);
