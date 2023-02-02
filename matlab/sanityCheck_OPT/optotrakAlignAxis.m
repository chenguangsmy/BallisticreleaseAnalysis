% optotrakAlignAxis 

% This function align axis of the optotrak system with the WAM 
% The idea is:  
% [X_wam - X_wam0] = A*[X_opt - X_opt0]

%%%%%%%%%%%%%%%%%%%%%%%% for sprigns that has two bolts connected 
sstmp = SessionScan(3983); 

t = sstmp.data.t;
X_wam = sstmp.data.x;
X_opt =[sstmp.data.optx(1,:); 
        sstmp.data.opty(1,:);
        sstmp.data.optz(1,:);];
    
fh = figure(); 
axh(1) = subplot(2,1,1);
plot(t,X_wam); 
axh(2) = subplot(2,1,2);
plot(axh(1), t, X_wam);
plot(axh(2), t, X_opt);
linkaxes(axh, 'x');

t0_range = [43204.488526237 43209.9018135048];
t0_idx = t > t0_range(1) & t < t0_range(2);
X_wam0 = nanmean(X_wam(:,t0_idx),2);
X_opt0 = nanmean(X_opt(:,t0_idx),2);

t_idx = t>4.3182*1e4 & ~isnan(X_opt(1,:));
W_mat = X_wam(:,t_idx) -X_wam0;
O_mat = X_opt(:,t_idx) -X_opt0;
A = (W_mat * O_mat')*inv(O_mat * O_mat');

X_opt1 = A*(X_opt - X_opt0) + X_wam0;

fh = figure(); 
plot(t, X_wam(2,:), '.'); hold on;
plot(t, X_opt1(2,:), '.');
legend('WAM', 'OPT')

%%  
%%%%%%%%%%%%%%%%%%%%%%%% for subject that has only one bolts connected at
%%%%%%%%%%%%%%%%%%%%%%%% handle
sstmp = SessionScan(4006); 

t = sstmp.data.t;
X_wam = sstmp.data.x;
X_opt =[sstmp.data.optx(1,:); 
        sstmp.data.opty(1,:);
        sstmp.data.optz(1,:);];
    
fh = figure(); 
axh(1) = subplot(2,1,1);
plot(t,X_wam); 
axh(2) = subplot(2,1,2);
plot(axh(1), t, X_wam);
plot(axh(2), t, X_opt);
linkaxes(axh, 'x');

t0_range = [4880 4890];
t0_idx = t > t0_range(1) & t < t0_range(2);
X_wam0 = nanmean(X_wam(:,t0_idx),2);
X_opt0 = nanmean(X_opt(:,t0_idx),2);

t_idx = t>4890 & ~isnan(X_opt(1,:));
W_mat = X_wam(:,t_idx) -X_wam0;
O_mat = X_opt(:,t_idx) -X_opt0;
A = (W_mat * O_mat')*inv(O_mat * O_mat');

X_opt1 = A*(X_opt - X_opt0) + X_wam0;

fh = figure(); 
plot(t, X_wam(2,:), '.'); hold on;
plot(t, X_opt1(2,:), '.');
legend('WAM', 'OPT')

%%
%%%%%%%%%%%%%%%%%%%%%%%% for robot remount at the horizontal places.
%%%%%%%%%%%%%%%%%%%%%%%% handle
sstmp = SessionScan(4058); 

t = sstmp.data.t;
X_wam = sstmp.data.x;
X_opt =[sstmp.data.optx(1,:); 
        sstmp.data.opty(1,:);
        sstmp.data.optz(1,:);];
    
fh = figure(); 
axh(1) = subplot(2,1,1);
plot(t,X_wam); 
axh(2) = subplot(2,1,2);
plot(axh(1), t, X_wam);
plot(axh(2), t, X_opt);
linkaxes(axh, 'x');

t0_range = [7490 7510];
t0_idx = t > t0_range(1) & t < t0_range(2);
X_wam0 = nanmean(X_wam(:,t0_idx),2);
X_opt0 = nanmean(X_opt(:,t0_idx),2);

t_idx = t>7510 & ~isnan(X_opt(1,:));
W_mat = X_wam(:,t_idx) -X_wam0;
O_mat = X_opt(:,t_idx) -X_opt0;
A = (W_mat * O_mat')*inv(O_mat * O_mat');

X_opt1 = A*(X_opt - X_opt0) + X_wam0;

fh = figure(); 
plot(t, X_wam(2,:), '.'); hold on;
plot(t, X_opt1(2,:), '.');
legend('WAM', 'OPT')

%% %%
%%%%%%%%%%%%%%%%%%%%%%%% for robot remount at the horizontal places, "inmotion"
%%%%%%%%%%%%%%%%%%%%%%%% handle
sstmp = SessionScan(4139); 

t = sstmp.data.t;
X_wam = sstmp.data.x;
X_opt =[sstmp.data.optx(1,:); 
        sstmp.data.opty(1,:);
        sstmp.data.optz(1,:);];
    
fh = figure(); 
axh(1) = subplot(2,1,1);
plot(t,X_wam); 
axh(2) = subplot(2,1,2);
plot(axh(1), t, X_wam);
plot(axh(2), t, X_opt);
linkaxes(axh, 'x');

t0_range = [586 590];
t0_idx = t > t0_range(1) & t < t0_range(2);
X_wam0 = nanmean(X_wam(:,t0_idx),2);
X_opt0 = nanmean(X_opt(:,t0_idx),2);

t_idx = t>638 & t<648 & ~isnan(X_opt(1,:));
W_mat = X_wam(:,t_idx) -X_wam0;
O_mat = X_opt(:,t_idx) -X_opt0;
A = (W_mat * O_mat')*inv(O_mat * O_mat');

X_opt1 = A*(X_opt - X_opt0) + X_wam0;

fh = figure(); 
plot(t, X_wam(2,:), '.'); hold on;
plot(t, X_opt1(2,:), '.');
legend('WAM', 'OPT')

%% %% 
%%%%%%%%%%%%%%%%%%%%%%%% for robot remount at the horizontal places, "inmotion"
%%%%%%%%%%%%%%%%%%%%%%%% handle, after I accidentally move the WAM
sstmp = SessionScan(4296); 

t = sstmp.data.t;
X_wam = sstmp.data.x;
X_opt =[sstmp.data.optx(1,:); 
        sstmp.data.opty(1,:);
        sstmp.data.optz(1,:);];
    
fh = figure(); 
axh(1) = subplot(2,1,1);
plot(t,X_wam); 
axh(2) = subplot(2,1,2);
plot(axh(1), t, X_wam);
plot(axh(2), t, X_opt);
linkaxes(axh, 'x');

t0_range = [2137.6 2137.7];
t0_idx = t > t0_range(1) & t < t0_range(2);
X_wam0 = mean(X_wam(:,t0_idx),2, 'omitnan');
X_opt0 = mean(X_opt(:,t0_idx),2, 'omitnan');

t_idx = t>2158 & t<2165 & ~isnan(X_opt(1,:));
W_mat = X_wam(:,t_idx) -X_wam0;
O_mat = X_opt(:,t_idx) -X_opt0;
A = (W_mat * O_mat')*inv(O_mat * O_mat');

X_opt1 = A*(X_opt - X_opt0) + X_wam0;

fh = figure(); 
plot(t, X_wam(2,:), '.'); hold on;
plot(t, X_opt1(2,:), '.');
legend('WAM', 'OPT')

%% %% 
%%%%%%%%%%%%%%%%%%%%%%%% for robot remount at the horizontal places, "inmotion"
%%%%%%%%%%%%%%%%%%%%%%%% handle, after I accidentally move the WAM
sstmp = SessionScan(4393); 

t = sstmp.data.t;
X_wam = sstmp.data.x;
X_opt =[sstmp.data.optx(1,:); 
        sstmp.data.opty(1,:);
        sstmp.data.optz(1,:);];
    
fh = figure(); 
axh(1) = subplot(2,1,1);
plot(t,X_wam); 
axh(2) = subplot(2,1,2);
plot(axh(1), t, X_wam);
plot(axh(2), t, X_opt);
linkaxes(axh, 'x');

t0_range = [2380 2400];
t0_idx = t > t0_range(1) & t < t0_range(2);
X_wam0 = mean(X_wam(:,t0_idx),2, 'omitnan');
X_opt0 = mean(X_opt(:,t0_idx),2, 'omitnan');

t_range = [2423 2430];
t_idx = t>t_range(1) & t<t_range(2) & ~isnan(X_opt(1,:));
W_mat = X_wam(:,t_idx) -X_wam0;
O_mat = X_opt(:,t_idx) -X_opt0;
A = (W_mat * O_mat')*inv(O_mat * O_mat');

X_opt1 = A*(X_opt - X_opt0) + X_wam0;

fh = figure(); 
plot(t, X_wam(2,:), '.'); hold on;
plot(t, X_opt1(2,:), '.');
legend('WAM', 'OPT')


%% %% 
%%%%%%%%%%%%%%%%%%%%%%%% for robot remount at the horizontal places, "inmotion"
%%%%%%%%%%%%%%%%%%%%%%%% Brady moves the robot somehow
sstmp = SessionScan(4415,1); 

t = sstmp.data.t;
X_wam = sstmp.data.x;
X_opt =[sstmp.data.optx(1,:); 
        sstmp.data.opty(1,:);
        sstmp.data.optz(1,:);];
    
fh = figure(); 
axh(1) = subplot(2,1,1);
plot(t,X_wam); 
axh(2) = subplot(2,1,2);
plot(axh(1), t, X_wam);
plot(axh(2), t, X_opt);
linkaxes(axh, 'x');

t0_range = [90 110];
t0_idx = t > t0_range(1) & t < t0_range(2);
X_wam0 = mean(X_wam(:,t0_idx),2, 'omitnan');
X_opt0 = mean(X_opt(:,t0_idx),2, 'omitnan');

t_range = [164 174];
t_idx = t>t_range(1) & t<t_range(2) & ~isnan(X_opt(1,:));
W_mat = X_wam(:,t_idx) -X_wam0;
O_mat = X_opt(:,t_idx) -X_opt0;
A = (W_mat * O_mat')*inv(O_mat * O_mat');

X_opt1 = A*(X_opt - X_opt0) + X_wam0;

fh = figure(); 
subplot(3,1,1);
plot(t, X_wam(1,:), '.'); hold on;
plot(t, X_opt1(1,:), '.');
legend('WAM', 'OPT')
xlabel('t (s)'); ylabel('pos x (m)');
subplot(3,1,2);
plot(t, X_wam(2,:), '.'); hold on;
plot(t, X_opt1(2,:), '.');
legend('WAM', 'OPT')
xlabel('t (s)'); ylabel('pos y (m)');
subplot(3,1,3);
plot(t, X_wam(3,:), '.'); hold on;
plot(t, X_opt1(3,:), '.');
legend('WAM', 'OPT')
xlabel('t (s)'); ylabel('pos z (m)');