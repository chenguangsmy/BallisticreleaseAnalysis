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
