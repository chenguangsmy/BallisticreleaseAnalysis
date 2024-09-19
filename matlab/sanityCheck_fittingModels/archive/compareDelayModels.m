load('data_tmp.mat');

% a delayed model with state-space
% Define the parameters
m = params.m;  % Mass (kg)
b = params.b;  % Damping coefficient (Ns/m)
k = params.k;  % Spring constant (N/m)
tau = params.tau;
tau = -0.3;    % delay must be negative 

% Define the state-space matrices without delay
A = [0 1; -k/m -b/m];
B = [0; 1/m];
C = [1 0];
D = 0;

% Create the state-space model (without delay)
sys = ss(A, B, C, D);

tau_list = [0.05 0.15 0.30];
% k_frac_list = [0 0.1 0.2 0.5 0.7 1];
k_frac_list = [-0.1 -0.07 -0.05  0.05  0.07 0.1];
output_cell = cell(3,6);


for tau_i = 1:3
    tau = tau_list(tau_i)
    for k_i = 1:length(k_frac_list)
%         k_0 = k*(1/(1+k_frac_list(k_i)));
        k_0 = k%*(1/(1+k_frac_list(k_i)));
        k_f = k*(k_frac_list(k_i)/(1+k_frac_list(k_i)))
        % Define the state-space matrices without delay
        A = [0 1; -k_0/m -b/m];
        B = [0; 1/m];
        C = [1 0];
        D = 0;

        % Create the state-space model (without delay)
        sys = ss(A, B, C, D);

        % Define the delayed state matrix A_d
        A_d = [0 0; -k_f/m 0];
        B_d = [0; 0];
        C_d = [0 0];
        D_d = 0;

        DelayT = struct('delay', tau, 'a', A_d, 'b', B_d, 'c', C_d, 'd', D_d);
        sys_total = delayss(A, B, C, D, DelayT);
        output = lsim(sys_total, force_interp_t, data_eg.t);
        output_cell{tau_i, k_i} = output;
    end

    
end


data_eg.x = data_eg.x - data_eg.x(1);
color_map = colormap('spring');
for tau_i = 1:3
    axh(tau_i) = subplot(3,1,tau_i); hold on;
    set(axh(tau_i), 'linewidth', 1,'fontsize', 12)
    tau = -tau_list(tau_i);
%     plot(data_eg.t, output0); hold on;
    plot(data_eg.t, data_eg.x, 'color', [0 0 1], 'linewidth', 2)
    for k_i = 1:length(k_frac_list)
        plot(data_eg.t, output_cell{tau_i, k_i}, 'linewidth', 1, ...
            'Color', color_map(k_i*40,:));
    end

    title(['Delay for ' num2str(tau_list(tau_i)) 's']);
    legend('data',...
           'K_f/K = 0', ...
           'K_f/K = 0.1', ...
           'K_f/K = 0.2', ...
           'K_f/K = 0.5', ...
           'K_f/K = 0.7', ...
           'K_f/K = 1.0');
end

xlabel('t (s)');
ylabel(axh(1), 'x (m)');
ylabel(axh(2), 'x (m)');
ylabel(axh(3), 'x (m)');

linkaxes(axh, 'xy')

sgtitle('Compare Delay Models')
% Combine both systems (using a linear time-delay system model)
%%
% calculate the fitting percent 
fit_percent = zeros(3,6)
for tau_i = 1:3
     for k_i = 1:6
%         fit_tmp = (1 - sum(abs(output_cell{tau_i, k_i} - data_eg.x'))/sum(abs(data_eg.x)))*100;
%         fit_tmp = (sum((output_cell{tau_i, k_i} - data_eg.x').^2)/length(output_cell{tau_i, k_i})); % mean squared error
%          fit_tmp = (sum((output_cell{tau_i, k_i} - data_eg.x').^2)/length(output_cell{tau_i, k_i})) / ...
%              (sum((data_eg.x' - mean(data_eg.x)).^2)/length(output_cell{tau_i, k_i}));
         fit_tmp = sqrt(sum((output_cell{tau_i, k_i} - data_eg.x').^2)/length(output_cell{tau_i, k_i})) / ...
             sqrt(sum((data_eg.x' - mean(data_eg.x)).^2)/length(output_cell{tau_i, k_i}));
        fit_percent(tau_i, k_i) = 100 - fit_tmp*100;
     end
end
%% Simulate the response to an input force
figure;
Ts = mean(diff(data_eg.t));
disp_interp_t = data_eg.x - data_eg.x(1);
force_interp_t = -(data_eg.F - mean(data_eg.F(1:200)));

output0= lsim(sys, force_interp_t, data_eg.t);
output = lsim(sys_total, force_interp_t, data_eg.t);

plot(data_eg.t, output0); hold on;
plot(data_eg.t, output); hold on;
title('Step Response of the System with Time Delay');
%%
grid on;

fold = 2;
force_interp_t_long = [force_interp_t repmat(force_interp_t(end),1, (fold-1)*length(force_interp_t))];
data_eg.t_long = [data_eg.t data_eg.t - data_eg.t(1) + data_eg.t(end) +  data_eg.t(2) - data_eg.t(1)]; % think how to make fold-able
output = lsim(sys_total, force_interp_t_long, data_eg.t_long);
% output = lsim(sys_delay, force_interp_t, data_eg.t);
plot(data_eg.t_long, output)
title('Step Response of the System with Time Delay');
grid on;
