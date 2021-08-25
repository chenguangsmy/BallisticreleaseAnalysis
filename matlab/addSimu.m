% addSimu
% simulation the 1D ballistic-release behavior use an 2nd order ODE
%% add simulations 
% specify parameters;
%F_list = [3 5 10 13 15 20];
F_list = [3 5 10 13 15 20];
figure(1); hold on;
figure(2); hold on;
for F_i = 1:length(F_list)
    m  = 1;     % kg
    F  = F_list(F_i);   % N
    x0 = 0.05;     % m
    ks = F/x0;   % N/m  % assume subject can alter his/her stiffness
    kr = 25;    % N/(m*s)
    % specify ODEds
    % m*x'' = ks(x0 - x) - kr*x'
    % x'(0) = 0;
    % x(0) = 0;
    syms x(t)
    Dx = diff(x);
    
    % ode = diff(x, t, 2) == 1/m * ks * (x0 - x) - kr * Dx;
    ode = diff(x, t, 2) == eval('1/m') * eval('ks')*(eval('x0') - x) - eval('kr') * Dx;
    cond1 = x(0) == 0;
    cond2 = Dx(0)== 0;
    
    conds = [cond1 cond2];
    % get solutions
    xSol(t) = dsolve(ode, conds);
    xSol = simplify(xSol);
    
    t = -2:0.002:3;
    x_result = xSol(t);
    x_result_numerical= double(x_result);
    x_resultd_numerical= diff(x_result_numerical);
    tdiff = t(1:end-1);
    x_result_numerical(t<0) = 0;
    x_resultd_numerical(tdiff<0) = 0;
    figure(1);
    plot(t, x_result_numerical);
    figure(2);
    plot(t(1:end-1), diff(x_result_numerical)./diff(t));
end

figure(1); 
xlim([-0.2, 0.5]); 
%ylim([0, 0.1]);
xlabel('time (s)');
ylabel('distance (m)');
legend('3N', '5N', '10N', '13N', '15N', '20N');
title('position on different force case, x0 = 5cm');
figure(2);
xlim([-0.2, 0.5]); 
%ylim([0, 0.1]);
xlabel('time (s)');
legend('3N', '5N', '10N', '13N', '15N', '20N');
ylabel('speed (m/s)');
title('speed on different force case, x0 = 5cm');


%% concatinate sessions 
% Three sessions: 2077, 2081, 2082. load all of these trials and save in
% another file

ss2077 = SessionScan(2077);
ss2081 = SessionScan(2081);
ss2082 = SessionScan(2082);
FThs = unique([ss2077.fThs, ss2081.fThs, ss2082.fThs]);
all_trials = [ss2077.trials, ss2081.trials, ss2082.trials];
for trial_i = 1:length(all_trials)
    if isempty(all_trials(trial_i).tarR)
        all_trials(trial_i).tarR = -1;
    end
end
all_trials_clean = all_trials([all_trials.tarR]==2);
ss2080 = ss2077; % combo trials here
ss2080.trials = all_trials_clean;
%ss2080.plotMeantrialForce_sameCond();
ss2080.plotMeantrialVel_sameCond();
ss2080.plotMeantrialPos_sameCond();


%% Use simulation to test Scott's method if right
trialtmp = ss2333.trials(20); % an example trial in the example session;  
trialtmp = strialtmp.simuTrialusingODE(226, 28, 3) % K, D, M
trialtmp_pred = trialtmp.predictImpedanceLinDev2ndOrderFixM()

% only shift M
M_list = 1:5; K = 320; D = 3;
prop_pred = zeros(length(M_list), 3);
for trial_i = 1:length(M_list)
    trialtmp = trialtmp.simuTrialusingODE(K, D, M_list(trial_i));
    trialpred= trialtmp.predictImpedanceLinDev2ndOrderFixM();
    prop_pred(trial_i,:) = [trialpred.pred_K, trialpred.pred_D, trialpred.pred_A];
end
% only shift K
M = 3; K_list = [160 240 320 400 480]; D = 3;
prop_pred = zeros(length(K_list), 3);
for trial_i = 1:length(K_list)
    trialtmp = trialtmp.simuTrialusingODE(K_list(trial_i), D, M);
    trialpred= trialtmp.predictImpedanceLinDev2ndOrderFixM();
    prop_pred(trial_i,:) = [trialpred.pred_K, trialpred.pred_D, trialpred.pred_A];
end
% only shift D
M = 3; K = 320; D_list = (1:5)*2;
prop_pred = zeros(length(D_list), 3);
for trial_i = 1:length(K_list)
    trialtmp = trialtmp.simuTrialusingODE(K, D_list(trial_i), M);
    trialpred= trialtmp.predictImpedanceLinDev2ndOrderFixM();
    prop_pred(trial_i,:) = [trialpred.pred_K, trialpred.pred_D, trialpred.pred_A];
end

M_list = 1*ones(20); K = 320; D = 3;
prop_pred = zeros(length(M_list), 3);
for trial_i = 1:length(M_list)
    trialtmp = trialtmp.simuTrialusingODE(K, D, M_list(trial_i));
    trialpred= trialtmp.predictImpedanceLinDev2ndOrderFixM();
    prop_pred(trial_i,:) = [trialpred.pred_K, trialpred.pred_D, trialpred.pred_A];
end

%% Simscope simulation, and prediction use Scott's method
x  = out.pos.Data;
dx = out.vel.Data;
%ddx= [0; diff(dx)]; % this one needs resample
ddx= [diff(dx); 0]; % this one needs resample
F  = out.force2.Data(2:end);
%F  = out.force.Data(1:end-1);
length(x)
length(dx)
length(ddx)
length(F)
% f = m*ddx + b*dx - k*x + k*x0
X = [ones(size(x(2:end))), x(2:end), dx(2:end), ddx(1:end-1)];
b = (X'*X) \ (X'*F);
% b(1) = kx0; b(2) = -k; b(3) = b; b(4) = m
m = b(4);
B = b(3);
K = -b(2);
x0 = b(1)/K;
fprintf('m: %fkg, B: %fN/(m/s)^-1, K: %fN/m, x0: %fm\n', m, B, K, x0);

%% plot simulation resuts in a different figure
pos_data = out.pos.Data;
vel_data = out.vel.Data;
%fce_data = out.force.Data;
fce_data = out.fce.Data;
%fce2_data= out.force2.Data;
time     = out.tout;
figure(); 
subplot(3,1,1);
plot(time, pos_data);
ylabel('position');
subplot(3,1,2);
plot(time, vel_data);
ylabel('velosity');
subplot(3,1,3);
plot(time, fce_data);
ylabel('force');
%subplot(4,1,4);
%plot(time, fce2_data);
%ylabel('force2');
%xlabel('time');
%legend('position', 'velocity', 'force');
subplot(3,1,1);
title('K 320, B 15, M 1');

%% execute a series of simulink 
% force_all = [5 10 15 20];
target_all = [2.5 5.0 7.5 10.0]/100;
stiffness_all = [320, 160, 107, 89];
stiffness_col = ['rgbc'];
for stiffness_i = 1:length(stiffness_all)
    % target_set = target_all(target_i);
    stiffness_set   = stiffness_all(stiffness_i);
    force_all       = target_all * stiffness_set;
    % stiffness_all   = force_all/target_set; 
    for target_i = 1:length(target_all)
        stiffness_set = stiffness_all(stiffness_i);
        simOut(target_i, stiffness_i) = sim('../ballisticReleaseSimu/SpringMass_2019_show');
    end
end
% plot in the figure
for target_i = 1:3
    figure(); hold on;
    target_set = target_all(target_i);
    stiffness_all = force_all/target_set; 
    for stiffness_i = 1:4
        pos_data = simOut(target_i, stiffness_i).pos.Data;
        vel_data = simOut(target_i, stiffness_i).vel.Data;
        fce_data = simOut(target_i, stiffness_i).fce.Data;
        time     = simOut(target_i, stiffness_i).tout;
        subplot(3,1,1); hold on;
        plot(time, pos_data, 'color', stiffness_col(stiffness_i));
        subplot(3,1,2); hold on;
        plot(time, vel_data, 'color', stiffness_col(stiffness_i));
        subplot(3,1,3); hold on;
        plot(time, fce_data, 'color', stiffness_col(stiffness_i));
    end
    subplot(3,1,1);
    title('K 320, B 15, M 1');
    ylabel('position');
    subplot(3,1,2);
    ylabel('velosity');
    subplot(3,1,3);
    ylabel('force');
    xlabel('time');
    legend({[num2str(stiffness_all(1)) 'N/m'],...
        [num2str(stiffness_all(2)) 'N/m'],...
        [num2str(stiffness_all(3)) 'N/m'],...
        [num2str(stiffness_all(4)) 'N/m']});
end

%% stiffness range simulation
% Assuming the subject-handle coupling is 2 spring coupling (steady-state), 
% asking how large the stiffness could be
x0_arr = [0.025, 0.05, 0.075, 0.10];
for x0i = 1:length(x0_arr)
    x0 = x0_arr(x0i); % m
    k0 = 2500; % N/m
    F  = 0:1:floor(x0*k0);
    ks = F./(x0 - F./k0);
    loglog(F, ks, '*');
    hold on;
end
grid on;
legend({'2.5cm', '5cm', '7.5cm', '10cm'})
xlabel('required Force (N)');
ylabel('theoretical stiffness (N/m)');
title('Theoretical stiffness with interacting with WAM');

%% x0 range simulation
% Assuming the subject-handle coupling is 2 spring coupling (steady-state),
% giving the existing stiffness range, ask how far the x0 could be 
x0_arr = [0.025, 0.05, 0.075, 0.10];
kmax = 400; kmin = 100;
figure(); hold on;
for x0i = 1:length(x0_arr)
    x0 = x0_arr(x0i); % m
    k0 = 2500; % N/m
    F  = 0:1:30;
    x0_prac = zeros(size(F));
    for Fi = 1:length(F)
        x1 = F(Fi)/k0;
        ks = F(Fi)/(x0 - F(Fi)/k0);
        x0_prac(Fi) = x0;
        if ks>kmax
            x0_prac(Fi) = F(Fi)/kmax + x1;
        elseif ks<kmin
            x0_prac(Fi) = F(Fi)/kmin + x1;
        end
    end
    plot(F, x0_prac, '*');
end  
grid on;
legend({'2.5cm', '5cm', '7.5cm', '10cm'})
xlabel('required Force (N)');
ylabel('possible Practical x0 (N/m)');
title('Theoretical x0 with interacting with WAM');  

%% overlap the perturb simulation
% according to the previous method, only looking at the perturbation in a
% relatively short zone. 
% 2.5cm: 3N : 8N
%   5cm: 5N : 17N
% 7.5cm: 8N : 21N
% 0.1cm: 7N : 21N
%Force_list = {[3, 6], [6, 9, 12, 15], [9:3:21], [9:3:21]};
Force_list = [3 9 15 21];
%dist = [2.5 5 7.5 10]/100;
spring_list = [320, 160, 106.67, 89]; % N/m
%dist = [10 7.5 5 2.5]/100;
k0 = 2500;
colors = colormap('lines');
stiffness0 = 300; % robot stiffness
%for dist_i = 1:length(dist)
%for spring_i = 1:length(spring_list)
for fce_i = 1:length(spring_list)
    fce_list = Force_list;%{dist_i};
    %dist = fce_list/spring_list(spring_i);
    %fce_list = dist(dist_i) * spring_list; 
    %for fce_i = 1:length(fce_list)
    for spring_i = 1:length(fce_list)
        %x0 = dist(spring_i);
        fce = fce_list(fce_i);
        dist = fce/spring_list(spring_i);
        x0 = dist;
        %stiffness = fce/(x0-fce/k0);
        %stiffness = fce/x0;
        stiffness = spring_list(spring_i);
        damping = 10;
        xr0 = fce/stiffness0;
        %stiffness_mat(fce/3, dist_i) = stiffness;
        %simout(dist_i, fce_i)=sim('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/ballisticReleaseSimu/ballisticRelease');
        simout(fce_i, spring_i)=sim('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/ballisticReleaseSimu/ballisticRelease_stepPert',...
            'FixedStep','0.002');
    end    
end
%% plot out 
figure(); 

for fce_i = 1:length(fce_list)
%for dist_i = 1:length(dist)
%for spring_i = 1:length(spring_list)
    fce_list = spring_list*dist;%dist(dist_i); 
    %subplot(1,length(fce_list),fce_i); hold on;
    %subplot(1,length(dist),dist_i); hold on;
    subplot(1,length(spring_list),fce_i); hold on;
    %figure(); hold on;
    %x0 = 0.05;
    %fce_list = Force_list{dist_i};
    %for dist_i = 1:length(dist)
    %for fce_i = 1:length(fce_list)
    for spring_i = 1:length(spring_list)
        pos = simout(fce_i, spring_i).pos.Data;
        postime= simout(fce_i, spring_i).tout;
        posidx = postime>0.5 & postime<0.8;
        pos0= mean(pos(posidx));
        vel = simout(fce_i, spring_i).vel.Data;
        time = simout(fce_i, spring_i).vel.Time - 1; % 1 for pert, 4 for release
        fce = simout(fce_i, spring_i).fce.Data;
        timeF= simout(fce_i, spring_i).fce.Time;
        %plot(time, pos-pos0, 'color', colors(fce_i,:));
        %plot(time, vel, 'color', colors(fce_i,:));
        plot(time, pos, 'color', colors(spring_i,:));
        %plot(time, pos, 'color', colors(dist_i,:));
        %plot(time, fce-fce(1), 'color', colors(fce_i,:));
        %plot(time, fce, 'color', colors(fce_i,:));
        %legend_arr{fce_i} = [num2str(fce_list(fce_i)) 'N'];
        %legend_arr{dist_i} = [num2str(dist(dist_i)) 'm'];
        legend_arr{fce_i} = [num2str(spring_list(fce_i)) 'N/m'];
    end
    legend(legend_arr);
    %ylim([-18, 35]); % force 
    %ylim([-0.02, 0.06]);
    %ylim([-0.01, 0.16]); % position
    %ylim([-0.01, 0.25]); % position, bigger
    ylim([-0.16, 0.01]); % -position
    %ylim([-0.15, 0.6]);
    %ylim([-0.35, 0.35]); % velocity
    xlim([-0.1, 1.2]);
    %ylim([-1.1, 1.2]); % velocity
    %if dist_i == 1
    if dist_i == 1
        %ylabel('position (m)')
        %ylabel('velocity (m/s)')
        ylabel('censored force (N)')
        xlabel('time at movement (s)');
    else
        %set(gca, 'yTickLabel', {});
    end
    set(gca, 'Ygrid', 'on');
    %title(['target ' num2str(dist(dist_i)*100) 'cm']);
    %title(['force ' num2str(fce_list(fce_i)) 'N']);
end

%% recognize damping through different simulation values 
% generate data

%dist = [2.5 5 7.5 10]/100;
colors = colormap('lines');
damping_list = 2.^(0:6); 
stiffness0 = 300; % robot stiffness
for damp_i = 1:length(damping_list)
        stiffness = stiffness0;
        damping = damping_list(damp_i);
        simout_damping(damp_i)=sim('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/ballisticReleaseSimu/ballisticRelease_stepPert',...
            'FixedStep','0.002');  
end
% report set damping 
damp_est = zeros(size(damping_list));
% calculate damping through: 
%   1. c = 2*k*T*(\delta / sqrt((\delta)^2 + (2*pi)^2 ));   ... delta, K, and T
%   2. delta = log(x1/x3);                                  ... peak
for damp_i = 1:length(damping_list)
    k = stiffness;
    T = 1; %...
    %x1 %= ...; % findpeak(...,1 );
    %x3 %= ...; % findpeak(...,2 );
    
    vel = simout_damping(damp_i).vel;
    vel_idx = vel.Time > 1 & vel.Time < 2;
    vel_select = vel.Data(vel_idx);
    vel_selectTime = vel.Time(vel_idx);
    [peaks, locs ] = findpeaks(vel_select);
    if length(peaks)<2
        fprintf("Overdamped in damping: %f",damping_list(damp_i));
        damp_est(damp_i) = -1;
    else
        T = vel_selectTime(locs(2)) - vel_selectTime(locs(1));
        delta = log(peaks(1)/peaks(2));
        c = 2*stiffness*(T/(2*pi))*(delta / sqrt((delta)^2 + (2*pi)^2 ));
        damp_est(damp_i) = c;
    end
end
% report "measured" damping 
damping_list
damp_est

%% only test using spring-mass-damper system to test damping
colors = colormap('lines');
damping_list = 2.^(0:6); 
stiffness0 = 300; % robot stiffness
for damp_i = 1:length(damping_list)
        stiffness = stiffness0;
        damping = damping_list(damp_i);
        simout_damping(damp_i)=sim('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/ballisticReleaseSimu/SpringMassDamper2019',...
            'FixedStep','0.001');  
end
% report set damping 
damp_est = zeros(size(damping_list));
% calculate damping through: 
%   1. c = 2*k*T*(\delta / sqrt((\delta)^2 + (2*pi)^2 ));   ... delta, K, and T
%   2. delta = log(x1/x3);                                  ... peak
for damp_i = 1:length(damping_list)
    k = stiffness;
    T = 1; %...
    %x1 %= ...; % findpeak(...,1 );
    %x3 %= ...; % findpeak(...,2 );
    
    vel = simout_damping(damp_i).vel;
    vel_idx = find(ones(size(vel.Data)));%vel.Time > 1 & vel.Time < 2;
    vel_select = vel.Data(vel_idx);
    vel_selectTime = vel.Time(vel_idx);
    [peaks, locs ] = findpeaks(vel_select);
    if length(peaks)<2
        fprintf("Overdamped in damping: %f",damping_list(damp_i));
        damp_est(damp_i) = -1;
    else
        T = vel_selectTime(locs(2)) - vel_selectTime(locs(1));
        delta = log(peaks(1)/peaks(2));
        c = 2*stiffness*(T/(2*pi))*(delta / sqrt((delta)^2 + (2*pi)^2 ));
        damp_est(damp_i) = c;
    end
end
damping_list
damp_est

%% to test perturbation at equilibrium position
SpringStiff_list = [160, 320, 640, 960];
PertFce_list = [5,10,15,20,25];
x0 = 0;
xr0 = 0;
colors = colormap('lines');
for stiffness_i = 1:length(SpringStiff_list)
    %fce_list = Force_list{dist_i};
    for PertFce_i = 1:length(PertFce_list)
        stiffness = SpringStiff_list(stiffness_i);
        PertFce = PertFce_list(PertFce_i);
        simout(stiffness_i, PertFce_i)=sim('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/ballisticReleaseSimu/ballisticRelease_stepPert',...
            'FixedStep','0.002');
    end    
end
%%
figure(); 
SpringStiff_list = [160, 320, 640, 960];
PertFce_list = [5,10,15,20,25];
for stiffness_i = 1:length(SpringStiff_list)
    axh = subplot(1,length(SpringStiff_list),stiffness_i); hold on;
    for PertFce_i = 1:length(PertFce_list)
        pos = simout(stiffness_i, PertFce_i).pos.Data;
        postime= simout(stiffness_i, PertFce_i).tout;
        posidx = postime>0.5 & postime<0.8;
        pos0= mean(pos(posidx));
        vel = simout(stiffness_i, PertFce_i).vel.Data;
        time = simout(stiffness_i, PertFce_i).vel.Time - 1;
        fce = simout(stiffness_i, PertFce_i).fce.Data;
        timeF= simout(stiffness_i, PertFce_i).fce.Time;
        %plot(time, pos-pos0, 'color', colors(PertFce_i,:));
        %plot(time, vel, 'color', colors(fce_i,:));
        %plot(time, pos, 'color', colors(fce_i,:));
        %plot(time, pos, 'color', colors(dist_i,:));
        plot(time, fce, 'color', colors(PertFce_i,:));
        legend_arr{PertFce_i} = [num2str(PertFce_list(PertFce_i)) 'N'];
        %legend_arr{dist_i} = [num2str(dist(dist_i)) 'm'];
    end
    legend(legend_arr);
    %ylim([-0.01, 0.04]);
    %ylim([-0.02, 0.06]);
    %ylim([-0.04, 0.14]); % position
    %ylim([-0.08, 0.03]); % position of perturbation
    %ylim([-0.15, 0.6]);
    ylim([-5, 30]); % force
    xlim([-0.08, 1.0]);
    %if dist_i == 1
    if stiffness_i == 1
        ylabel('Force (N)')
        %ylabel('endpoint position (m)')
        %ylabel('velocity (m/s)')
        xlabel('time at perturbation (s)');
    else
        %set(gca, 'yTickLabel', {});
    end
    set(gca, 'Ygrid', 'on');
    %title(['target ' num2str(dist(dist_i)*100) 'cm']);
    title(['stiffness ' num2str(SpringStiff_list(stiffness_i)) 'N/m']);
end
