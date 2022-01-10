% Andy don't believe force could suddenly change, thus I do the following
% check:

% check condition: let the WAM be frictionless (K=0,D=0), and use hand hold
% the spring, and to the pull and see the force response.  

% cpDatarg2(3679); % copy the data with hand release
ss_tmp = SessionScan(3679);
ss_tmp1= SessionScan(3680);

%% 
figure(); 
axh(1) = subplot(2,1,1); 
plot(ss_tmp.data.t, ss_tmp.data.x(2,:));
title('y position');

axh(2) = subplot(2,1,2);
plot(ss_tmp.data.t, ss_tmp.data.f(2,:));
ylim([-25 5]);
title('y force');
xlabel('blackrock time(s)');

linkaxes(axh, 'x');
sgtitle('hand');

%%
%ss_tmp1= SessionScan(3680);
figure();
t_offset = 2255.4702;
%ss_tmp1.data.t = ss_tmp1.data.t - t_offset;


axh(1) = subplot(3,1,1); 
%plot(ss_tmp1.data.t, ss_tmp1.data.x(2,:));
plot(ss_tmp1.data.t, ss_tmp1.data.x(2,:), 'Marker', '.');
title('y position');
grid on;

axh(2) = subplot(3,1,2); 
%plot(ss_tmp1.data.t, ss_tmp1.data.x(2,:));
plot(ss_tmp1.data.t, ss_tmp1.data.v(2,:), 'Marker', '.');
title('y velocity');
grid on;

axh(3) = subplot(3,1,3);
plot(ss_tmp1.data.t, ss_tmp1.data.f(2,:), 'Marker', '.');
ylim([-25 5]);
xlim([-1 1]);
title('y force');
xlabel('blackrock time(s)');
grid on

linkaxes(axh, 'x');
sgtitle('scissor');