% Compare measurements on different methods. For now, only compare
% stiffness of different measurement across subject (1) and subject (2) 

% subj(1): Chenguang
% subj(2): Himanshu

% run fede_sys_id to get the following results data 
% fede_sys_id.m
load("compareMeasurements_corsssubj.mat", "result_subj");
% result_subj{1} = results;
% result_subj{2} = results;
% save("compareMeasurements_corsssubj.mat", "result_subj", '-append');

%%%

f_val   = [ 15  15  7.5
            30  30  15
            45  45  22.5];
%         K900   K600    K300
d_val   = [ 1.6     2.5     2.5        
            3.3     5       5          
            5       7.5     7.5];

k_val = [900    600     300
         900    600     300
         900    600     300];

% do the figure of task setting 
figure('Units','inches','Position',[0 0 3 3]); hold on; 
plot(d_val(:), f_val(:), '.', 'MarkerSize', 10);

x_grids = [0:0.1:8]; 
y_grids = k_val(1,:)'*x_grids/100;
plot(x_grids, y_grids, 'LineStyle',':', 'Color', 'b')
ylim([0, 50]);
xlabel('distance (cm)');
ylabel('force (N)');
title('experiment task design');

% fig1, plot the release measurement with the pulse measurement

% say, just use subj1

results = result_subj{1};

K_est_rls_avg = results.K_up_avg;
K_est_rls_std = results.K_up_std;

K_est_pul = results.K_p{1}.p(:,:,21:28);
K_est_pul_avg = mean(K_est_pul,3);
K_est_pul_std = std(K_est_pul,[],3);

%%
% fig1, plot the bar 
figure();
axh(1) = subplot(1,2,1);
hold on;
bar(k_val', K_est_rls_avg');
yline([300 600 900]);
title('release estimation');
xticks([300 600 900]);
legend('low', 'middle', 'high');

axh(2) = subplot(1,2,2);
hold on;
bar(k_val', K_est_pul_avg');
yline([300 600 900]);
title('pulse estimation');
xticks([300 600 900]);

linkaxes(axh, 'y');

% fig2, compare across force 
figure();
bar(k_val(:,2:3)', K_est_pul_avg(:,2:3)');
% legend(d_val(:,3));
legend('7.5cm', '5cm', '2.5cm');
title('different K with same displacement');

% fig3, compare across disp 
figure();
bar(k_val(:,1:2)', K_est_pul_avg(:,1:2)');
% legend(d_val(:,3));
legend('15N', '30N', '45N');
title('different K with same force');