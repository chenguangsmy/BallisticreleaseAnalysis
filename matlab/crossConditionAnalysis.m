classdef crossConditionAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.
    
    properties
        
        depMeasures
        dexSubject
        dexForce
        dexDistance
        sfrq
        x_r
        f_target_vec
        distVal
        
    end
    
    methods
        function [this] = crossConditionAnalysis(data,dexSubject,dexForce,dexDistance,subjectType)
            
            % Define index to use (will replace soon as passed parameter)
            this.dexSubject = dexSubject;
            this.dexForce = dexForce;
            this.dexDistance = dexDistance;
            this.sfrq = 500;            
            this.f_target_vec = [15,20,25];
            this.distVal = [2.5 5 7.5];

            
            %Compute all stiffness estimates
            this.depMeasures = this.get_depMeasures(data,subjectType);
%             save('/Users/jhermus/Desktop/test_spring.mat');
%             load('/Users/jhermus/Desktop/test_spring.mat');

            this.get_mainPlot(subjectType,1);
            
%             this.plot_postionForce_pulse(data);
%             this.plot_positionForce_stocastic(data);
%             this.plot_postion_release(data);

            % Make ANOVA Plots
            % Run ANOVA between subjects

             disp([subjectType, ' estimate complete.']);
            
        end
        
        function [depMeasures] = get_depMeasures(this,data,subjectType)
            
            sizeData = size(data);
            
            for subj = this.dexSubject
                for force = this.dexForce
                    for dist = this.dexDistance
                        
                        %temporarly take each out seperately (FIX later)
                        clear dataCross Tmp_depMeasures
                        for trial = 1:sizeData(4)
                            for pret = 1:sizeData(5)
                                dataCross{trial,pret} = data{subj,force,dist,trial,pret};
                            end
                        end
                        
                        disp(['subj = ',int2str(subj),...
                              ', force = ',int2str(force),...
                              ', dist = ',int2str(dist)]);
                        depMeasures{subj,force,dist} = crossTrialAnalysis(dataCross,this.sfrq,this.f_target_vec(force),subjectType);
                        
                    end
                end
            end
            
        end
        
        function [] = get_mainPlot(this,subjectType,subNum)
            
            %% Human Estimates
            if(strcmp(subjectType,'human'))
            distVal = [2.5 5 7.5];
            springVal = [0, 160, 320, 640];
            forceVal = {['F =',int2str(this.f_target_vec(1)),'N'],...
                        ['F =',int2str(this.f_target_vec(2)),'N'],...
                        ['F =',int2str(this.f_target_vec(3)),'N']};
            xRange = [1.5 8.5];

            
%             plotMeasure = this.pullDataFromStruct();
                        
            % Summary across trial
            for subj = this.dexSubject
                for force = this.dexForce
                    for dist = this.dexDistance
                        structPlotTrial{subj,force,dist}.est_norm_release = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_norm.est_release);
                        structPlotTrial{subj,force,dist}.est_catch_release = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_catch.est_release);
                        structPlotTrial{subj,force,dist}.est_catch_pulse  = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_catch.est_pulse);
                        structPlotTrial{subj,force,dist}.est_stocastic_release  = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_stocastic.est_release);
                        structPlotTrial{subj,force,dist}.est_stocastic_stocastic  = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_stocastic.est_stocastic);
                    end
                end
            end
            
            %Summary Across conditions
            plotCon_norm_release = this.get_estStructCondition(structPlotTrial,'est_norm_release');
            plotCon_catch_release = this.get_estStructCondition(structPlotTrial,'est_catch_release');
            plotCon_catch_pulse = this.get_estStructCondition(structPlotTrial,'est_catch_pulse');
            plotCon_stocastic_release = this.get_estStructCondition(structPlotTrial,'est_stocastic_release');
            plotCon_stocastic_stocastic = this.get_estStructCondition(structPlotTrial,'est_stocastic_stocastic');
            
%             force = 1;
%             dist = 1;
%             plotCon_norm_release.x_hat
%             [cost,x_hat,x_dot_hat] = this.depMeasures{1,1,1}.costFunc(x,x_dot,x_ddot,x_r,this.f_traget_vec(force),X_hat);
%             
%             [cost,x_hat,x_dot_hat] = costFunc(this,x,x_dot,x_ddot,x_r,Fp,X_hat)
%             
%             Y_model = tf(1,[this.depMeasures{1,1,1}.m,B_opt,K_opt]);
%             h_model = impulse(Y_model,0:1/this.sfrq:L/this.sfrq);

            % Stiffness Plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,5,1); hold on;
            ax2 = subplot(1,5,2); hold on;
            ax3 = subplot(1,5,3); hold on;
            ax4 = subplot(1,5,4); hold on;
            ax5 = subplot(1,5,5); hold on;

            axes(ax1); title('Release (trial)');
            ylabel('Stiffness (N/m)');ylim([0 1500]);
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release (catch)');
            ylim([0 1500]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
                       
            axes(ax3); title('Pulse (catch)');
            ylim([0 1500]); %ylabel('Stiffness (N/m)');
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax4); title('Stocastic (release)');
            ylim([0 1500]); %ylabel('Stiffness (N/m)');
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax5); title('Stocastic (stocastic)');
            ylim([0 1500]); %ylabel('Stiffness (N/m)');
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;

            for i = this.dexForce
               axes(ax1); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_norm_release.k_hat_ave(subNum,i,:)),squeeze(plotCon_norm_release.k_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax2); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_catch_release.k_hat_ave(subNum,i,:)),squeeze(plotCon_catch_release.k_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax3); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_catch_pulse.k_hat_ave(subNum,i,:)),squeeze(plotCon_catch_pulse.k_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax4); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_stocastic_release.k_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_release.k_hat_std(subj,i,:)),'linewidth',2.5); hold on;
               axes(ax5); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_stocastic_stocastic.k_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_stocastic.k_hat_std(subj,i,:)),'linewidth',2.5); hold on;

            end
            axes(ax1); legend(forceVal);
            
            % Damping Plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,5,1); hold on;
            ax2 = subplot(1,5,2); hold on;
            ax3 = subplot(1,5,3); hold on;
            ax4 = subplot(1,5,4); hold on;
            ax5 = subplot(1,5,5); hold on;
            
            axes(ax1); title('Release (normal)');
            ylabel('Damping (N-s/m)');ylim([0 200]);
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release (catch)');
            ylim([0 200]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
                       
            axes(ax3); title('Pulse (catch)');
            ylim([0 200]); %ylabel('Stiffness (N/m)');
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax4); title('Stocastic (release)');
            ylabel('Damping (N-s/m)');ylim([0 200]);
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax5); title('Stocastic (stocastic)');
            ylim([0 200]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;

            for i = this.dexForce
               axes(ax1); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_norm_release.b_hat_ave(subNum,i,:)),squeeze(plotCon_norm_release.b_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax2); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_catch_release.b_hat_ave(subNum,i,:)),squeeze(plotCon_catch_release.b_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax3); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_catch_pulse.b_hat_ave(subNum,i,:)),squeeze(plotCon_catch_pulse.b_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax4); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_stocastic_release.b_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_release.b_hat_std(subj,i,:)),'linewidth',2.5); hold on;
%                axes(ax5); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_stocastic_stocastic.b_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_stocastic.b_hat_ave(subj,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax1); legend(forceVal);

            % X0 plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,5,1); hold on;
            ax2 = subplot(1,5,2); hold on;
            ax3 = subplot(1,5,3); hold on;
            ax4 = subplot(1,5,4); hold on;
            ax5 = subplot(1,5,5); hold on;
            
            axes(ax3); title('Pulse');
            ylabel('x_0 (m)');ylim([0 10]);
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
%             axes(ax2); title('Release');
%             ylim([0 10]); % ylabel('Stiffness (N/m)'); 
%             xlabel('Distance'); xticks(this.distVal); xlim(xRange);
%             set(gca,'fontsize',16); grid on;
%             
%             axes(ax3); title('Stocastic');
%             ylim([0 10]); %ylabel('Stiffness (N/m)');
%             xlabel('Distance'); xticks(this.distVal); xlim(xRange);
%             set(gca,'fontsize',16);grid on;
            

            for i = this.dexForce
               axes(ax3); errorbar(this.distVal(this.dexDistance),100*squeeze(plotCon_catch_pulse.x_h_ave(subNum,i,:)),100*squeeze(plotCon_catch_pulse.x_h_std(subNum,i,:)),'linewidth',2.5); hold on;
%                axes(ax2); errorbar(this.distVal(this.dexDistance),squeeze(k_release_mean(1,i,:)),squeeze(k_release_std(1,i,:)),'linewidth',2.5); hold on;
%                axes(ax3); errorbar(this.distVal(this.dexDistance),squeeze(k_stocastic_mean(1,i,:)),squeeze(k_stocastic_std(1,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax3); legend(forceVal);
            
            %% VAF ??
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,5,1); hold on;
            ax2 = subplot(1,5,2); hold on;
            ax3 = subplot(1,5,3); hold on;
            ax4 = subplot(1,5,4); hold on;
            ax5 = subplot(1,5,5); hold on;
            
            axes(ax1); title('Release (trial)');
            ylabel('VAF');ylim([0 100]);
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release (catch)');
            ylim([0 100]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
            
            axes(ax3); title('Pulse (catch)');
            ylim([0 100]); ylabel('Coherence');
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax4); title('Stocastic (release)');
            ylim([0 100]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
                        
            axes(ax5); title('Stocastic (stocastic)');
            ylim([0 1]); ylabel('Coherence');
            xlabel('Distance'); xticks(this.distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            for i = this.dexForce
                axes(ax1); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_norm_release.VAF_ave(subNum,i,:)),squeeze(plotCon_norm_release.VAF_std(subNum,i,:)),'linewidth',2.5); hold on;
                axes(ax2); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_catch_release.VAF_ave(subNum,i,:)),squeeze(plotCon_catch_release.VAF_std(subNum,i,:)),'linewidth',2.5); hold on;
                axes(ax3); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_catch_pulse.VAF_ave(subNum,i,:)),squeeze(plotCon_catch_pulse.VAF_std(subNum,i,:)),'linewidth',2.5); hold on;
                axes(ax4); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_stocastic_release.VAF_ave(subj,i,:)),squeeze(plotCon_stocastic_release.VAF_std(subj,i,:)),'linewidth',2.5); hold on;
                axes(ax5); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_stocastic_stocastic.OC_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_stocastic.OC_hat_std(subj,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax1); legend(forceVal,'location','south');

            
            % Plot summary figure for release fitting
            figure('Position',[300 314 929 420]);
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                ax(count) = subplot(length(this.dexForce),length(this.dexDistance),count); hold on;
                
                axes(ax(count)); title(['Force: ',num2str(this.f_target_vec(force)), ', Dist: ',num2str(distVal(dist))]);
                if(1==mod(count,length(this.dexDistance)))
                    ylabel('Velocity (m/s)'); %ylim([0 100]);
                end
                if(count>(length(this.dexForce)*length(this.dexDistance))-length(this.dexDistance))
                    xlabel('Time (s)'); %xticks(this.distVal); xlim(xRange);
                end
                set(gca,'fontsize',16);grid on;
                                
                for triall = 1:length(this.depMeasures{subNum,force,dist}.trial_norm.est_release)
                    structt = this.depMeasures{subNum,force,dist}.trial_norm.est_release{triall};
                    N = length(structt.dexRange);
                    plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structt.h_hat(structt.dexRange),'b');
                end
                
                % Plot mean model
                N = length(structPlotTrial{subj,force,dist}.est_norm_release.h_model_ave);
                plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structPlotTrial{subj,force,dist}.est_norm_release.h_model_ave,'r','linewidth',2.5);

                
                ylim([-0.02 0.03]);
%                 xlim([0 300]);
                count = count + 1;
                end
            end
%             axes(ax(1)); legend(forceVal,'location','south');

            
            % Plot summary figure for pulse fitting
            figure('Position',[300 314 929 420]);
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                ax(count) = subplot(length(this.dexForce),length(this.dexDistance),count); hold on;
                
                axes(ax(count)); title(['Force: ',num2str(this.f_target_vec(force)), ', Dist: ',num2str(distVal(dist))]);
                if(1==mod(count,length(this.dexDistance)))
                    ylabel('IRF'); %ylim([0 100]);
                end
                if(count>(length(this.dexForce)*length(this.dexDistance))-length(this.dexDistance))
                    xlabel('Time (s)'); %xticks(this.distVal); xlim(xRange);
                end
                set(gca,'fontsize',16);grid on;
                                
                for triall = 1:length(this.depMeasures{subNum,force,dist}.trial_catch.est_pulse)
                    structt = this.depMeasures{subNum,force,dist}.trial_catch.est_pulse{triall};
                    N = length(structt.x_dot_hat);
                    plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structt.x_dot_hat,'b');
                end
                
                % Plot mean model
                N = length(structPlotTrial{subj,force,dist}.est_catch_pulse.x_dot_hat_ave);
                plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structPlotTrial{subj,force,dist}.est_catch_pulse.x_dot_hat_ave,'r','linewidth',2.5);

                
                ylim([-0.15 0.07]);
                xlim([0 0.35]);
                count = count + 1;
                end
            end
            
%             axes(ax(1)); legend(forceVal,'location','south');

%             % Check learning
%             for i = this.dexForce
%                 figure;
%                 for j = this.dexDistance
%                     plot(squeeze(this.k_hat_pulse(1,i,j,:)),'-o'); hold on;
%                 end
%             end
                        
            %% Spring Estimates
            elseif(strcmp(subjectType,'spring'))
            springVal = [0, 160, 320, 640];
            forceVal = {['F =',int2str(this.f_target_vec(1)),'N'],...
                        ['F =',int2str(this.f_target_vec(2)),'N'],...
                        ['F =',int2str(this.f_target_vec(3)),'N']};
            xRange = [-200 800];
                        
            % Summary across trial
            for subj = this.dexSubject
                for force = this.dexForce
                    for dist = this.dexDistance
                        disp([int2str(force),' ',int2str(dist)]);
                        structPlotTrial{subj,force,dist}.est_norm_release = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_norm.est_release);
                        structPlotTrial{subj,force,dist}.est_catch_release = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_catch.est_release);
                        structPlotTrial{subj,force,dist}.est_catch_pulse  = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_catch.est_pulse);
                        structPlotTrial{subj,force,dist}.est_stocastic_release  = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_stocastic.est_release);
                        structPlotTrial{subj,force,dist}.est_stocastic_stocastic  = this.get_estStructCellTrial(this.depMeasures{subj,force,dist}.trial_stocastic.est_stocastic);
                    end
                end
            end
            
            %Summary Across conditions
            plotCon_norm_release = this.get_estStructCondition(structPlotTrial,'est_norm_release');
            plotCon_catch_release = this.get_estStructCondition(structPlotTrial,'est_catch_release');
            plotCon_catch_pulse = this.get_estStructCondition(structPlotTrial,'est_catch_pulse');
            plotCon_stocastic_release = this.get_estStructCondition(structPlotTrial,'est_stocastic_release');
            plotCon_stocastic_stocastic = this.get_estStructCondition(structPlotTrial,'est_stocastic_stocastic');
                
            % Stiffness Plot
            figure('Position',[300 314 929 420]);
            ax(1) = subplot(1,5,1); hold on;
            ax(2) = subplot(1,5,2); hold on;
            ax(3) = subplot(1,5,3); hold on;
            ax(4) = subplot(1,5,4); hold on;
            ax(5) = subplot(1,5,5); hold on;

            axes(ax(1)); title('Release (trial)');
            ylabel('Stiffness (N/m)');ylim([0 1000]);
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax(2)); title('Release (catch)');
            ylim([0 1000]); % ylabel('Stiffness (N/m)'); 
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
                       
            axes(ax(3)); title('Pulse (catch)');
            ylim([0 1000]); %ylabel('Stiffness (N/m)');
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax(4)); title('Release (stocastic)');
            ylim([0 1000]); %ylabel('Stiffness (N/m)');
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax(5)); title('Stocastic (stocastic)');
            ylim([0 1000]); %ylabel('Stiffness (N/m)');
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;

            for i = this.dexForce
               axes(ax(1)); errorbar(springVal,squeeze(plotCon_norm_release.k_hat_ave(subNum,i,:)),squeeze(plotCon_norm_release.k_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax(2)); errorbar(springVal,squeeze(plotCon_catch_release.k_hat_ave(subNum,i,:)),squeeze(plotCon_catch_release.k_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax(3)); errorbar(springVal,squeeze(plotCon_catch_pulse.k_hat_ave(subNum,i,:)),squeeze(plotCon_catch_pulse.k_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax(4)); errorbar(springVal,squeeze(plotCon_stocastic_release.k_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_release.k_hat_std(subj,i,:)),'linewidth',2.5); hold on;
               axes(ax(5)); errorbar(springVal,squeeze(plotCon_stocastic_stocastic.k_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_stocastic.k_hat_std(subj,i,:)),'linewidth',2.5); hold on;
            end
            for i = 1:5
                axes(ax(i)); plot(xRange,springVal.*ones(2,1),'--k','linewidth',2.5); hold on;
            end
            axes(ax1); legend(forceVal);
            
            % Damping Plot
            figure('Position',[300 314 929 420]);
            ax(1) = subplot(1,5,1); hold on;
            ax(2) = subplot(1,5,2); hold on;
            ax(3) = subplot(1,5,3); hold on;
            ax(4) = subplot(1,5,4); hold on;
            ax(5) = subplot(1,5,5); hold on;
            
            axes(ax(1)); title('Release (normal)');
            ylabel('Damping (N-s/m)');ylim([0 70]);
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax(2)); title('Release (catch)');
            ylim([0 70]); % ylabel('Stiffness (N/m)'); 
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
                       
            axes(ax(3)); title('Pulse (catch)');
            ylim([0 70]); %ylabel('Stiffness (N/m)');
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax(4)); title('Stocastic (release)');
            ylabel('Damping (N-s/m)');ylim([0 100]);
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax(5)); title('Stocastic (stocastic)');
            ylim([0 70]); % ylabel('Stiffness (N/m)'); 
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;

            for i = this.dexForce
               axes(ax(1)); errorbar(springVal,squeeze(plotCon_norm_release.b_hat_ave(subNum,i,:)),squeeze(plotCon_norm_release.b_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax(2)); errorbar(springVal,squeeze(plotCon_catch_release.b_hat_ave(subNum,i,:)),squeeze(plotCon_catch_release.b_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax(3)); errorbar(springVal,squeeze(plotCon_catch_pulse.b_hat_ave(subNum,i,:)),squeeze(plotCon_catch_pulse.b_hat_std(subNum,i,:)),'linewidth',2.5); hold on;
               axes(ax(4)); errorbar(springVal,squeeze(plotCon_stocastic_release.b_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_release.b_hat_std(subj,i,:)),'linewidth',2.5); hold on;
%                axes(ax5); errorbar(this.distVal(this.dexDistance),squeeze(plotCon_stocastic_stocastic.b_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_stocastic.b_hat_ave(subj,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax(1)); legend(forceVal);

            % X0 plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,5,1); hold on;
            ax2 = subplot(1,5,2); hold on;
            ax3 = subplot(1,5,3); hold on;
            ax4 = subplot(1,5,4); hold on;
            ax5 = subplot(1,5,5); hold on;
            
            axes(ax3); title('Pulse');
            ylabel('x_0 (m)');ylim([0 20]);
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
%             axes(ax2); title('Release');
%             ylim([0 10]); % ylabel('Stiffness (N/m)'); 
%             xlabel('Spring Stiffness (N/m)'); xticks(this.distVal); xlim(xRange);
%             set(gca,'fontsize',16); grid on;
%             
%             axes(ax3); title('Stocastic');
%             ylim([0 10]); %ylabel('Stiffness (N/m)');
%             xlabel('Spring Stiffness (N/m)'); xticks(this.distVal); xlim(xRange);
%             set(gca,'fontsize',16);grid on;
            

            for i = this.dexForce
               axes(ax3); errorbar(springVal,100*squeeze(plotCon_catch_pulse.x_h_ave(subNum,i,:)),100*squeeze(plotCon_catch_pulse.x_h_std(subNum,i,:)),'linewidth',2.5); hold on;
%                axes(ax2); errorbar(this.distVal(this.dexDistance),squeeze(k_release_mean(1,i,:)),squeeze(k_release_std(1,i,:)),'linewidth',2.5); hold on;
%                axes(ax3); errorbar(this.distVal(this.dexDistance),squeeze(k_stocastic_mean(1,i,:)),squeeze(k_stocastic_std(1,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax3); legend(forceVal);
            
            %% VAF ??
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,5,1); hold on;
            ax2 = subplot(1,5,2); hold on;
            ax3 = subplot(1,5,3); hold on;
            ax4 = subplot(1,5,4); hold on;
            ax5 = subplot(1,5,5); hold on;
            
            axes(ax1); title('Release (trial)');
            ylabel('VAF');ylim([0 100]);
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release (catch)');
            ylim([0 100]); % ylabel('Stiffness (N/m)'); 
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
            
            axes(ax3); title('Pulse (catch)');
            ylim([0 100]); ylabel('Coherence');
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax4); title('Stocastic (release)');
            ylim([0 100]); % ylabel('Stiffness (N/m)'); 
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
                        
            axes(ax5); title('Stocastic (stocastic)');
            ylim([0 1]); ylabel('Coherence');
            xlabel('Spring Stiffness (N/m)'); xticks(springVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            for i = this.dexForce
                axes(ax1); errorbar(springVal,squeeze(plotCon_norm_release.VAF_ave(subNum,i,:)),squeeze(plotCon_norm_release.VAF_std(subNum,i,:)),'linewidth',2.5); hold on;
                axes(ax2); errorbar(springVal,squeeze(plotCon_catch_release.VAF_ave(subNum,i,:)),squeeze(plotCon_catch_release.VAF_std(subNum,i,:)),'linewidth',2.5); hold on;
                axes(ax3); errorbar(springVal,squeeze(plotCon_catch_pulse.VAF_ave(subNum,i,:)),squeeze(plotCon_catch_pulse.VAF_std(subNum,i,:)),'linewidth',2.5); hold on;
                axes(ax4); errorbar(springVal,squeeze(plotCon_stocastic_release.VAF_ave(subj,i,:)),squeeze(plotCon_stocastic_release.VAF_std(subj,i,:)),'linewidth',2.5); hold on;
                axes(ax5); errorbar(springVal,squeeze(plotCon_stocastic_stocastic.OC_hat_ave(subj,i,:)),squeeze(plotCon_stocastic_stocastic.OC_hat_std(subj,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax1); legend(forceVal,'location','south');

            
            % Plot summary figure for release fitting
            figure('Position',[300 314 929 420]);
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                ax(count) = subplot(length(this.dexForce),length(this.dexDistance),count); hold on;
                
                axes(ax(count)); title(['Force: ',num2str(this.f_target_vec(force)), ', Dist: ',num2str(springVal(dist))]);
                if(1==mod(count,length(this.dexDistance)))
                    ylabel('Velocity (m/s)'); %ylim([0 100]);
                end
                if(count>(length(this.dexForce)*length(this.dexDistance))-length(this.dexDistance))
                    xlabel('Time (s)'); %xticks(this.distVal); xlim(xRange);
                end
                set(gca,'fontsize',16);grid on;
                                
                for triall = 1:length(this.depMeasures{subNum,force,dist}.trial_norm.est_release)
                    structt = this.depMeasures{subNum,force,dist}.trial_norm.est_release{triall};
                    N = length(structt.dexRange);
                    plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structt.h_hat(structt.dexRange),'b');
                end
                
                % Plot mean model
                N = length(structPlotTrial{subj,force,dist}.est_norm_release.h_model_ave);
                plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structPlotTrial{subj,force,dist}.est_norm_release.h_model_ave,'r','linewidth',2.5);

                
                ylim([-0.07 0.07]);
%                 xlim([0 300]);
                count = count + 1;
                end
            end
%             axes(ax(1)); legend(forceVal,'location','south');

            
            % Plot summary figure for pulse fitting
            figure('Position',[300 314 929 420]);
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                ax(count) = subplot(length(this.dexForce),length(this.dexDistance),count); hold on;
                
                axes(ax(count)); title(['Force: ',num2str(this.f_target_vec(force)), ', Dist: ',num2str(springVal(dist))]);
                if(1==mod(count,length(this.dexDistance)))
                    ylabel('IRF'); %ylim([0 100]);
                end
                if(count>(length(this.dexForce)*length(this.dexDistance))-length(this.dexDistance))
                    xlabel('Time (s)'); %xticks(this.distVal); xlim(xRange);
                end
                set(gca,'fontsize',16);grid on;
                                
                for triall = 1:length(this.depMeasures{subNum,force,dist}.trial_catch.est_pulse)
                    structt = this.depMeasures{subNum,force,dist}.trial_catch.est_pulse{triall};
                    N = length(structt.x_dot_hat);
                    plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structt.x_dot_hat,'b');
                end
                
                % Plot mean model
                N = length(structPlotTrial{subj,force,dist}.est_catch_pulse.x_dot_hat_ave);
                plot(0:(1/this.sfrq):(N-1)*(1/this.sfrq), structPlotTrial{subj,force,dist}.est_catch_pulse.x_dot_hat_ave,'r','linewidth',2.5);

                
                ylim([-0.3 0.3]);
                xlim([0 0.35]);
                count = count + 1;
                end
            end

%             x_s_known = 100*[0.0938, 0.0469, 0.0234,...
%                 0.1250, 0.0625, 0.0312 ,...
%                 0.1562, 0.0781, 0.0391];
%             k_plot_xs = [160,320,640,...
%                 160,320,640,...
%                 160,320,640];
            
%             plot(x_s_known,k_plot_xs,'.k','markerSize',30);
%             
%             if(strcmp(subjectType,'spring'))
%                 plot([0; 10].*ones(2,4),[160, 320, 640].*ones(2,4),'--k','linewidth',2.5);
%             end
                
            else
                error('Spesify subjectType');
            end
            
            %% Extra K vs. x0
%             colorVec = {'b','r','g'};
%             markerVec = {'o','+','*'};
%             figure;
%             subj = 1;
%             for force = this.dexForce
%                 for dist = this.dexDistance
%                     plot(1./squeeze(100*this.x0_hat_pulse(subj,force,dist,:)),...
%                         squeeze(this.k_hat_pulse(subj,force,dist,:)),...
%                         'color',colorVec{force},...
%                         'marker',markerVec{force},...
%                         'LineStyle','none','linewidth',1.5,'markersize',7); hold on;
%                 end
%             end
%             ylabel('Stiffness (N/m)'); ylim([0 1500]);
%             xlabel('1/x_0 (mm)'); %xlim([0 10]);
%             set(gca,'fontsize',16);
            
        end
        
        
        function [outputt]  = get_estStructCellTrial(this,structIn)
            
            for i = 1:length(structIn)
                outputt.k_hat(i) = structIn{i}.k_hat;
                
                if( isfield(structIn{i},'h_model') )
                    if(i == 1)
                        h_model = structIn{i}.h_model;
                    end
                    h_model = this.auto_padd(h_model,structIn{i}.h_model);
                end
                
               if( isfield(structIn{i},'x_dot_hat') )
                     if(i == 1)
                        x_dot_hat = structIn{i}.x_dot_hat';
                    end
                    x_dot_hat = this.auto_padd(x_dot_hat,structIn{i}.x_dot_hat');
                end
                
                if( isfield(structIn{i},'b_hat') )
                    outputt.b_hat(i) = structIn{i}.b_hat;
                end
                
                if( isfield(structIn{i},'VAF') )
                    outputt.VAF(i) = structIn{i}.VAF;
                end
                
                if( isfield(structIn{i},'x_h') )
                    outputt.x_h(i) = structIn{i}.x_h;
                end
                             
                if( isfield(structIn{i},'OC_hat') )
                    outputt.OC_hat(i) = structIn{i}.OC_hat;
                end
            end
            
            outputt.k_hat_ave = nanmean(outputt.k_hat);
            outputt.k_hat_std = nanstd(outputt.k_hat);
            
            if( isfield(structIn{i},'b_hat') )
                outputt.b_hat_ave = nanmean(outputt.b_hat);
                outputt.b_hat_std = nanstd(outputt.b_hat);
            end

            if( isfield(structIn{i},'VAF') )
                outputt.VAF_ave = nanmean(outputt.VAF);
                outputt.VAF_std = nanstd(outputt.VAF);
            end
            
            if( isfield(structIn{i},'x_h') )
                outputt.x_h_ave = nanmean(outputt.x_h);
                outputt.x_h_std = nanstd(outputt.x_h);
            end
            
            if( isfield(structIn{i},'OC_hat') )
                outputt.OC_hat_ave = nanmean(outputt.OC_hat);
                outputt.OC_hat_std = nanstd(outputt.OC_hat);
            end
            
            if( exist('h_model') )
                outputt.h_model_ave = nanmean(h_model,2);
            end
            
            if( exist('x_dot_hat') )
                outputt.x_dot_hat_ave = nanmean(x_dot_hat,2);
            end
                        
        end
        
        function [outputStruct] = get_estStructCondition(this,structPlotTrial,condition)

            for subj = this.dexSubject
                for force = this.dexForce
                    for dist = this.dexDistance
                        structIn = eval(['structPlotTrial{subj,force,dist}.',condition]);
                        
                        outputStruct.k_hat_ave(subj,force,dist) = structIn.k_hat_ave;
                        outputStruct.k_hat_std(subj,force,dist) = structIn.k_hat_std;
                        
                        if( isfield(structIn,'b_hat') )
                            outputStruct.b_hat_ave(subj,force,dist) = structIn.b_hat_ave;
                            outputStruct.b_hat_std(subj,force,dist) = structIn.b_hat_std;
                        end
                        
                        if( isfield(structIn,'VAF') )
                            outputStruct.VAF_ave(subj,force,dist) = structIn.VAF_ave;
                            outputStruct.VAF_std(subj,force,dist) = structIn.VAF_std;
                        end
                        
                        if( isfield(structIn,'x_h') )
                            outputStruct.x_h_ave(subj,force,dist) = structIn.x_h_ave;
                            outputStruct.x_h_std(subj,force,dist) = structIn.x_h_std;
                        end
                        
                        if( isfield(structIn,'OC_hat') )
                            outputStruct.OC_hat_ave(subj,force,dist) = structIn.OC_hat_ave;
                            outputStruct.OC_hat_std(subj,force,dist) = structIn.OC_hat_std;
                        end
                        
                    end
                end
            end
            
        end
        
        function [x_tot] = auto_padd(this,x_old, x_new)
            
            n_old = size(x_old,1);
            n_new = size(x_new,1);
            if(n_old > n_new) % padd x_new
                
                x_new = [x_new;nan(abs(n_old-n_new),1)];
                
            elseif(n_old < n_new) % padd x_old
                
                x_old = [x_old;nan(abs(n_old-n_new),size(x_old,2))];
                
            else
                % equal do nothing
            end
            
            x_tot = [x_old, x_new];
            
        end

            
        
        function [] = plot_postionForce_pulse(this,data)
            
            step_Pulse = 2;
            cf = 20; % cutoff freqnency
            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            
            colorVec = {'b','r','k'};
                        
            figure;
            sgtitle('Pulse Preturbation');
            subj = 1;
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                    
%                     subplot(length(this.dexForce),length(this.dexDistance),count);
                    for trial = 2:size(data,4)
                        clear tmp_data
                        if(~isempty(data{subj,force,dist,trial,step_Pulse}))
                            tmp_data = data{subj,force,dist,trial,step_Pulse};
                            dexPulse = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))+300];
                            
                            plot(tmp_data.x(2,dexPulse),colorVec{dist}); hold on;
                        end
                    end
                    count = count+1;
                    xlabel('Time'); ylabel('Postion (m)'); set(gca,'fontsize',16);     
                end
            end
            
            figure;
            sgtitle('Pulse Preturbation');
            subj = 1;
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                    
%                     subplot(length(this.dexForce),length(this.dexDistance),count);
                    for trial = 2:size(data,4)
                        clear tmp_data
                        if(~isempty(data{subj,force,dist,trial,step_Pulse}))
                            tmp_data = data{subj,force,dist,trial,step_Pulse};
                            dexPulse = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))+100];
                            v = filtfilt(b,a,tmp_data.v(2,dexPulse)); % apply fitler
                            plot(v,colorVec{dist}); hold on;
                        end
                        ylim([-0.2 0.2]);
                    end
                    count = count+1;
                    xlabel('Time'); ylabel('Velocity (m/s)'); set(gca,'fontsize',16);     
                end
            end
            
            %% Make Mean Plot (ADD!!!)
            
            figure;
            sgtitle('Pulse Preturbation');
            subj = 1;
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                    
                    subplot(length(this.dexForce),length(this.dexDistance),count);
                    for trial = 2:size(data,4)
                        clear tmp_data
                        if(~isempty(data{subj,force,dist,trial,step_Pulse}))
                        tmp_data = data{subj,force,dist,trial,step_Pulse};
                        dexPulse = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))];
                        f = filtfilt(b,a,tmp_data.f(2,dexPulse)); % apply fitler
                        plot(f); hold on;
                        end
                    end
                    count = count+1;
                    xlabel('Time'); ylabel('Force (N)'); set(gca,'fontsize',16);
                end
            end
            

            
        end
        
        function [] = plot_positionForce_stocastic(this,data)
            
            figure;
            sgtitle('Stocastic Preturbation');
            subj = 1;
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                    
                    subplot(length(this.dexForce),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,force,dist,trial,3};
                        if(~isempty(tmp_data))
                              dexStoc = [(min(find(tmp_data.ts==3)))+500 : (min(find(tmp_data.ts==4)))];
%                             dexStoc = [(min(find(tmp_data.Fp(2,:)~=0))) : (max(find(tmp_data.Fp(2,:)~=0)))];
                            plot(tmp_data.t(dexStoc),tmp_data.x(2,dexStoc)); hold on;
                        end
                    end
                    count = count+1;
                    xlabel('Time'); ylabel('Postion (m)'); set(gca,'fontsize',16);
                end
            end

            
            figure;
            sgtitle('Stocastic Preturbation');
            subj = 1;
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                    
                    subplot(length(this.dexForce),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,force,dist,trial,3};
                        if(~isempty(tmp_data))
                          dexStoc = [(min(find(tmp_data.ts==3)))+500 : (min(find(tmp_data.ts==4)))];
%                         dexStod = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))];
                        plot(tmp_data.f(2,dexStoc)); hold on;
                        end
                    end
                    count = count+1;
                    xlabel('Time'); ylabel('Force (N)'); set(gca,'fontsize',16);
                end
            end
 
        end
        
        function [] = plot_postion_release(this,data)
            
            figure;
            sgtitle('Release');
            subj = 1;
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                    
                    subplot(length(this.dexForce),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,force,dist,trial,1};
                        if(~isempty(tmp_data))
                        dexRelease = [max(find(tmp_data.ts== 4)) : max(find(tmp_data.ts == 6))];
                        plot(tmp_data.x(2,dexRelease)-0.481); hold on;
                        end
                    end
                    count = count+1;
                    ylim([0 0.1]);
                    xlabel('Time'); ylabel('Postion (m)'); set(gca,'fontsize',16);
                end
            end
            
            cf = 20; % cutoff freqnency
            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            
            figure;
            sgtitle('Release');
            subj = 1;
            count = 1;
            for force = this.dexForce
                for dist = this.dexDistance
                    
                    subplot(length(this.dexForce),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,force,dist,trial,1};
                        if(~isempty(tmp_data))
                        dexRelease = [max(find(tmp_data.ts== 4)) : max(find(tmp_data.ts == 6))];
                        v = filtfilt(b,a,tmp_data.v(2,dexRelease)); % apply fitler
                        plot(v); hold on;
                        end
                    end
                    count = count+1;
                    ylim([-0.5 0.5]);
                    xlabel('Time'); ylabel('Velocity (m/s)'); set(gca,'fontsize',16);
                end
            end

            

        end

    end
end