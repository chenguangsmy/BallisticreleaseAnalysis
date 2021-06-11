%% Process Data stiffness analysis

clear all
close all
clc
dbstop if error

addpath('func/');


%% Human analysis
% ss2005 = SessionScan(2005); % Human (test 1)
% ss2007 = SessionScan(2007); % Human (test 2)
% 
% % % Look at raw data
% axhf = plotMeantrialForce(ss2005);
% axhp = plotMeantrialPos(ss2005);
% axhv = plotMeantrialVel(ss2005);
% % 
% % % Combine preturbations data (spectral analysis)
% % ss2005.generateWamPertData(); 
% % ss2007.generateWamPertData(); 
% % ss2005.add_to_PertData(ss2007.wam.Data_pert); % (combine tests)
% % 
% % % Make Spectral plots
% % human1 = get_Z_spectral(ss2005.wam.Data_pert(1,:),ss2005.col_vec(1,:));
% % human2 = get_Z_spectral(ss2005.wam.Data_pert(2,:),ss2005.col_vec(2,:));
% % 
% % legendLabelss = {[num2str(ss2005.wam.Data_pert(1,1).FT),'N ',num2str(ss2005.wam.Data_pert(1,1).x0*100),'cm'],...
% %         [num2str(ss2005.wam.Data_pert(2,1).FT),'N ',num2str(ss2005.wam.Data_pert(2,1).x0*100),'cm']};
% % legend(legendLabelss);
% % 
% % 
% % dex1 = find((0.5<human1.TF_freq) & (human1.TF_freq<=2));
% % dex2 = find((0.5<human2.TF_freq) & (human2.TF_freq<=2));
% % 
% % figure; 
% % plot(human1.TF_freq(dex2),human1.hand_Z_22_mag(dex2)); hold on;
% % plot(human2.TF_freq(dex3),human2.hand_Z_22_mag(dex3)); hold on;
% % 
% % 
% % freq_mean = [mean(human1.hand_Z_22_mag(dex1)),mean(human2.hand_Z_22_mag(dex2))];
% % freq_sd = [std(human1.hand_Z_22_mag(dex1)),std(human2.hand_Z_22_mag(dex2))];
% %     
% % figure; 
% % errorbar([1,15],freq_mean,freq_sd,'lineStyle','none','linewidth',2.5); hold on;
% % ylim([0 800]); xlim([-5 20]); xticks([1 15]); xlabel('Target Force (N)'); ylabel('Stiffness (N/m)');
% set(gca,'fontsize',16); grid on;
% 
% % Ensemble analysis
% ss2005.generateWamPertData_ensemble();
% ss2007.generateWamPertData_ensemble();
% ss2005.add_to_ensemble(ss2007.wam.Data_pert_ensemble);
% ss2005.wam.Data_pert_ensemble.estimate = get_Z_ensemble(ss2005.wam.Data_pert_ensemble);
% 
% 
% %% Analysis with Springs
% ss1968 = SessionScan(1968);     % during this part already do [1]
% % connect data
% % Elliminates half the trials due to "failure" 
% ss1968.generateWamPertData();   %                          do [2]
% 
% % Make plots
% get_Z_spectral(ss1968.wam.Data_pert(1,:),ss1968.col_vec(1,:));
% 
% % K robot only 2928 N/m
% % Check estimate before and after realease tomorrow
% 
% % Import Springs
% ss1972 = SessionScan(1972); % 2 Spring (test 1)
% ss1972.generateWamPertData();
% ss1973 = SessionScan(1973); % 2 Spring (test 2)
% ss1973.generateWamPertData();
% 
% ss1972.add_to_PertData(ss1973.wam.Data_pert); % (combine tests)
% 
% ss1974 = SessionScan(1974); % 3 Spring (test 1)
% ss1974.generateWamPertData();
% ss1975 = SessionScan(1975); % 3 Spring (test 2)
% ss1975.generateWamPertData();
% 
% ss1974.add_to_PertData(ss1975.wam.Data_pert); % (combine tests)
% 
% spring2 = get_Z_spectral(ss1972.wam.Data_pert(1,:),ss1972.col_vec(1,:));
% spring3 = get_Z_spectral(ss1974.wam.Data_pert(1,:),ss1974.col_vec(2,:));
% 
% dex2 = find((0.5<spring2.TF_freq) & (spring2.TF_freq<=3));
% dex3 = find((0.5<spring3.TF_freq) & (spring3.TF_freq<=3));
% k = 157.6;
% k_spring2 = 2*k;
% k_spring3 = 3*k;
% 
% figure; 
% plot(spring2.TF_freq(dex2),spring2.hand_Z_22_mag(dex2)); hold on;
% plot(spring3.TF_freq(dex3),spring3.hand_Z_22_mag(dex3)); hold on;
% plot([spring2.TF_freq(dex2(1)),spring2.TF_freq(dex2(end))],[k_spring2,k_spring2],'--k','linewidth',2.5);
% plot([spring3.TF_freq(dex3(1)),spring3.TF_freq(dex3(end))],[k_spring3,k_spring3],'--k','linewidth',2.5);
% ylim([0 4*k]);
% 
% 
% freq_mean = [mean(spring2.hand_Z_22_mag(dex2)),mean(spring3.hand_Z_22_mag(dex3))];
% freq_sd = [std(spring2.hand_Z_22_mag(dex2)),std(spring3.hand_Z_22_mag(dex3))];
%     
% figure; 
% errorbar([2,3],freq_mean,freq_sd,'lineStyle','none','linewidth',2.5); hold on;
% plot([1.5 3.5],[k_spring2,k_spring2],'--k','linewidth',2.5);
% plot([1.5 3.5],[k_spring3,k_spring3],'--k','linewidth',2.5);
% ylim([0 4*k]); xlim([1.5 3.5]); xticks([2 3]); xlabel('Number of Springs'); ylabel('Stiffness (N/m)');
% set(gca,'fontsize',16);
% 
% 
% 
% 
% 
% %% Ensemble analysis
% ss1972.generateWamPertData_ensemble();
% ss1973 = SessionScan(1973);
% ss1973.generateWamPertData_ensemble();
% ss1972.add_to_ensemble(ss1973.wam.Data_pert_ensemble);
% ss1972.wam.Data_pert_ensemble.estimate = get_Z_ensemble(ss1972.wam.Data_pert_ensemble);
% 
% ss1974.generateWamPertData_ensemble();
% 
% ss1975 = SessionScan(1975);
% ss1975.generateWamPertData_ensemble();
% 
% ss1974.add_to_ensemble(ss1973.wam.Data_pert_ensemble);
% 
% ss1974.wam.Data_pert_ensemble.estimate = get_Z_ensemble(ss1974.wam.Data_pert_ensemble);
% 
% time2 = ss1972.wam.Data_pert_ensemble.estimate.time;
% time3 = ss1974.wam.Data_pert_ensemble.estimate.time;
% 
% k_hat_2 = squeeze(ss1972.wam.Data_pert_ensemble.estimate.K_hat(2,2,:));
% k_hat_3 = squeeze(ss1974.wam.Data_pert_ensemble.estimate.K_hat(2,2,:));
% 
% figure; 
% plot(time2,k_hat_2,'linewidth',2.5);hold on;
% plot(time3,k_hat_3,'linewidth',2.5);
% 
% k_spring = 2*157.6;
% k_robot = 2928;
% plot([time2(1),0,0,time2(end)],...
%     [k_spring+k_robot,k_spring+k_robot,k_spring,k_spring],'--k','linewidth',2.5);
% 
% k_spring = 3*157.6;
% k_robot = 2928;
% plot([time3(1),0,0,time3(end)],...
%     [k_spring+k_robot,k_spring+k_robot,k_spring,k_spring],'--k','linewidth',2.5);

%% Test sys ID on impulse preturbations 6/9/21 
%% Analysis with Springs
ss = SessionScan(2084);     % during this part already do [1]
% connect data
% Elliminates half the trials due to "failure" 
ss.generateWamPertData();   %                          do [2]

% Ensemble analysis
ss.generateWamPertData_ensemble();

get_Z_ensemble(ss.wam.Data_pert_ensemble);


