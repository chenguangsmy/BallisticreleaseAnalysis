% some optotrak points cannot be recorded under tha WAM setting on the
% "joystick" situation. The following code will delineate this phenomenon. 

% The 
ss_num = 4062;
sstmp = SessionScan(ss_num); 

idx_nan = isnan(sstmp.data.optx);
idx_nnan=~isnan(sstmp.data.optx); 

%% 1. plot the relationship of end-point-position and the isnan
figure('name', 'value in position space');
subplot(1,2,1);
hold on;
plot(sstmp.data.x(1,idx_nnan), sstmp.data.x(2,idx_nnan), '.', 'color', 'g');
plot(sstmp.data.x(1,idx_nan), sstmp.data.x(2,idx_nan), '.', 'color', 'r');
xlabel('WAM -x (m)');
ylabel('WAM -y (m)');
legend('valid value', 'nan value');
title('nan value at the cartesian space');

subplot(1,2,2);
hold on;
plot(sstmp.data.jp(4,idx_nnan), sstmp.data.jp(3,idx_nnan), '.', 'color', 'g');
plot(sstmp.data.jp(4,idx_nan), sstmp.data.jp(3,idx_nan), '.', 'color', 'r');
xlabel('WAM -j4 (m)');
ylabel('WAM -j3 (m)');
legend('valid value', 'nan value');
title('nan value at the cartesian space');

%% 2. plot the relationship of end-point velocity and the isnan 
idx_inrange = sstmp.data.x(2,:) > 0 & sstmp.data.x(2,:) < 0.06;
figure('name', 'value in velocity space');
hold on; 
plot(sstmp.data.v(1,~idx_inrange), sstmp.data.v(2,~idx_inrange), '.', 'color', [0.5 0.5 0.5]);
plot(sstmp.data.v(1,idx_inrange&idx_nan), sstmp.data.v(2,idx_inrange&idx_nan), '.', 'color', 'r'); 
plot(sstmp.data.v(1,idx_inrange&idx_nnan), sstmp.data.v(2,idx_inrange&idx_nnan), '.', 'color', 'g'); 
xlabel('WAM -vx (m/s)');
ylabel('WAM -vy (m/s)');
legend('out-range', 'valid value', 'nan value');
title('nan value at the velocity space');

%% 3. See the overlapped point of wam and optotrak
figure(); 
hold on; 
plot(sstmp.data.t, sstmp.data.x);
plot(sstmp.data.t, [sstmp.data.optx(1,:);sstmp.data.opty(1,:);sstmp.data.optz(1,:)], '.');
title(['session' num2str(ss_num)]);