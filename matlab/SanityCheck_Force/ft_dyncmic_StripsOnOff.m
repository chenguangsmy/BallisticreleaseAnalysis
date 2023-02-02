% Test if the strips on/off affects the measurement 

clear; clf; close all; clc; 

idx_on = [46794	59931 76168 ]; 
idx_off= [215549 232237 245899];

fttmp = SessionScanFT(4291);
dat_range = -200:1200; 
t_range = dat_range/1000; % s

figure('name', 'strip effect of force measurements');
title('strip effect of force measurements')
hold on; 
% axh(1) = subplot(2,1,1); title('on' ); hold on;
for i = 1:3
    lnh{i} = plot(t_range, fttmp.force_origin(:,idx_on(i)+dat_range), 'color', 'b');
end


% axh(2) = subplot(2,1,2); title('off' ); hold on;

for i = 1:3
    lnh{i+3} = plot(t_range, fttmp.force_origin(:,idx_off(i)+dat_range), 'color', 'r');
end

% linkaxes(axh, 'x');
xlabel('t (s)');
ylabel('F (N)');

legend([lnh{1}(1), lnh{4}(1)], {'strip on','strip off'});