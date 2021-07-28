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
    % specify ODE
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
force_all = [5 10 15 20];
target_all = [0.05 0.075 0.10];
stiffness_col = ['rgbc'];
for target_i = 1:3
    target_set = target_all(target_i);
    stiffness_all = force_all/target_set; 
    for stiffness_i = 1:4
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