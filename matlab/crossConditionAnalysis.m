classdef crossConditionAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.
    
    properties
        
        dexSubject
        dexDirection
        dexDistance
        sfrq
        k_hat_pulse
        k_hat_stocastic
        k_hat_release
        x0_hat_pulse
        x_r
        k_nonNan
        k_hat_pulse_VAF
        k_hat_release_VAF
        OC_hat_stocastic
        
    end
    
    methods
        function [this] = crossConditionAnalysis(data,dexSubject,dexDirection,dexDistance,subjectType)
            
            % Define index to use (will replace soon as passed parameter)
            this.dexSubject = dexSubject;
            this.dexDirection = dexDirection;
            this.dexDistance = dexDistance;
            this.sfrq = 500;
            
            %Compute all stiffness estimates
            this.get_depMeasures(data,subjectType);
            for subNum = this.dexSubject
                this.get_mainPlot(subjectType,subNum);
            end
%             this.plot_postionForce_pulse(data);
%             this.plot_positionForce_stocastic(data);
%             this.plot_postion_release(data);

            % Make ANOVA Plots
            % Run ANOVA between subjects

             disp([subjectType, ' estimate complete.']);
            
        end
        
        function [] = get_depMeasures(this,data,subjectType)
            
            sizeData = size(data);
            f_target_vec = [15,20,25];
            
            for subj = this.dexSubject
                for dir = this.dexDirection
                    for dist = this.dexDistance
                        
                        %temporarly take each out seperately (FIX later)
                        clear dataCross Tmp_depMeasures
                        for trial = 1:sizeData(4)
                            for pret = 1:sizeData(5)
                                dataCross{trial,pret} = data{subj,dir,dist,trial,pret};
                            end
                        end
                        
                        Tmp_depMeasures = crossTrialAnalysis(dataCross,this.sfrq,f_target_vec(dir),subjectType);
                        
                        if(length(Tmp_depMeasures.k_hat_pulse)~=0)
                            this.k_hat_pulse(subj,dir,dist,:) = Tmp_depMeasures.k_hat_pulse;
                            this.x0_hat_pulse(subj,dir,dist,:) = Tmp_depMeasures.x0_hat_pulse;
                            this.k_hat_pulse_VAF(subj,dir,dist,:) = Tmp_depMeasures.k_hat_pulse_VAF;
                        end
                        
                        if(length(Tmp_depMeasures.k_hat_stocastic)~=0)
                            this.k_hat_stocastic(subj,dir,dist,:) = Tmp_depMeasures.k_hat_stocastic;
                            this.OC_hat_stocastic(subj,dir,dist,:) = Tmp_depMeasures.OC_hat_stocastic;

                        end
                        
                        if(length(Tmp_depMeasures.k_hat_release)~=0)
                            this.k_hat_release(subj,dir,dist,:) = Tmp_depMeasures.k_hat_release;
                            this.k_hat_release_VAF(subj,dir,dist,:) = Tmp_depMeasures.k_hat_release_VAF;
                        end
                        
                    end
                end
            end
            
        end
        
        function [] = get_mainPlot(this,subjectType,subNum)
            
            distVal = [2.5 5 7.5];
            %             direcVal = {'front','back','left','right'};
            direcVal = {'F = 15N','F = 20N','F = 25N'};
            
            % Look at Pulse Estimates
            k_pulse_mean = nanmean(this.k_hat_pulse(subNum,:,:,1:end),4);
            k_pulse_std = nanstd(this.k_hat_pulse(subNum,:,:,1:end),0,4);
            x0_pulse_mean = nanmean(this.x0_hat_pulse(subNum,:,:,1:end),4);
            x0_pulse_std = nanstd(this.x0_hat_pulse(subNum,:,:,1:end),0,4);
            
            k_pulse_VAF_mean = nanmean(this.k_hat_pulse_VAF(subNum,:,:,1:end),4);
            k_pulse_VAF_std = nanstd(this.k_hat_pulse_VAF(subNum,:,:,1:end),0,4);
            
            this.k_nonNan = squeeze(sum(~isnan(this.k_hat_pulse),4));
            
            k_release_mean = nanmean(this.k_hat_release(subNum,:,:,:),4);
            k_release_std = nanstd(this.k_hat_release(subNum,:,:,:),0,4);
            
            k_release_VAF_mean = nanmean(this.k_hat_release_VAF(subNum,:,:,1:end),4);
            k_release_VAF_std = nanstd(this.k_hat_release_VAF(subNum,:,:,1:end),0,4);
            
            k_stocastic_mean = nanmean(this.k_hat_stocastic(subNum,:,:,:),4);
            k_stocastic_std = nanstd(this.k_hat_stocastic(subNum,:,:,:),0,4);
            
            OC_hat_stocastic_mean = nanmean(this.OC_hat_stocastic(subNum,:,:,:),4);
            OC_hat_stocastic_std = nanstd(this.OC_hat_stocastic(subNum,:,:,:),0,4);
            
            xRange = [1.5 8.5];
                        
            %% Human Estimates
            if(strcmp(subjectType,'human'))

            % Stiffness Plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,3,1); hold on;
            ax2 = subplot(1,3,2); hold on;
            ax3 = subplot(1,3,3); hold on;
            
            axes(ax1); title('Pulse');
            ylabel('Stiffness (N/m)');ylim([0 1500]);
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release');
            ylim([0 1500]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
                       
            axes(ax3); title('Stocastic');
            ylim([0 1500]); %ylabel('Stiffness (N/m)');
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;

            for i = this.dexDirection
               axes(ax1); errorbar(distVal(this.dexDistance),squeeze(k_pulse_mean(1,i,:)),squeeze(k_pulse_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax2); errorbar(distVal(this.dexDistance),squeeze(k_release_mean(1,i,:)),squeeze(k_release_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax3); errorbar(distVal(this.dexDistance),squeeze(k_stocastic_mean(1,i,:)),squeeze(k_stocastic_std(1,i,:)),'linewidth',2.5); hold on;

            end
            axes(ax1); legend(direcVal);

            % X0 plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,3,1); hold on;
            ax2 = subplot(1,3,2); hold on;
            ax3 = subplot(1,3,3); hold on;
            
            axes(ax1); title('Pulse');
            ylabel('x_0 (m)');ylim([0 10]);
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release');
            ylim([0 10]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
            
            axes(ax3); title('Stocastic');
            ylim([0 10]); %ylabel('Stiffness (N/m)');
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            

            for i = this.dexDirection
               axes(ax1); errorbar(distVal(this.dexDistance),100*squeeze(x0_pulse_mean(1,i,:)),100*squeeze(x0_pulse_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax2); errorbar(distVal(this.dexDistance),squeeze(k_release_mean(1,i,:)),squeeze(k_release_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax3); errorbar(distVal(this.dexDistance),squeeze(k_stocastic_mean(1,i,:)),squeeze(k_stocastic_std(1,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax1); legend(direcVal);
            
            %% VAF ??
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,3,1); hold on;
            ax2 = subplot(1,3,2); hold on;
            ax3 = subplot(1,3,3); hold on;
            
            axes(ax1); title('Pulse');
            ylabel('VAF');ylim([0 100]);
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release');
            ylim([0 100]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
            
            axes(ax3); title('Stocastic');
            ylim([0 1]); ylabel('Coherence');
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            

            for i = this.dexDirection
               axes(ax1); errorbar(distVal(this.dexDistance),squeeze(k_pulse_VAF_mean(1,i,:)),squeeze(k_pulse_VAF_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax2); errorbar(distVal(this.dexDistance),squeeze(k_release_VAF_mean(1,i,:)),squeeze(k_release_VAF_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax3); errorbar(distVal(this.dexDistance),squeeze(OC_hat_stocastic_mean(1,i,:)),squeeze(OC_hat_stocastic_std(1,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax1); legend(direcVal,'location','south');


%             % Check learning
%             for i = this.dexDirection
%                 figure;
%                 for j = this.dexDistance
%                     plot(squeeze(this.k_hat_pulse(1,i,j,:)),'-o'); hold on;
%                 end
%             end
                        
            %% Spring Estimates
            elseif(strcmp(subjectType,'spring'))
                
             % Stiffness Plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,3,1); hold on;
            ax2 = subplot(1,3,2); hold on;
            ax3 = subplot(1,3,3); hold on;
            
            axes(ax1); title('Pulse');
            ylabel('Stiffness Estimate (N/m)');ylim([0 1500]);
            xlabel('Stiffness (N/m)'); xticks([160, 320, 640]); xlim([0 800]);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release');
            ylim([0 1500]); % ylabel('Stiffness (N/m)'); 
            xlabel('Stiffness (N/m)'); xticks([160, 320, 640]); xlim([0 800]);
            set(gca,'fontsize',16); grid on;
                       
            axes(ax3); title('Stocastic');
            ylim([0 1500]); %ylabel('Stiffness (N/m)');
            xlabel('Stiffness (N/m)'); xticks([160, 320, 640]); xlim([0 800]);
            set(gca,'fontsize',16);grid on;

            for i = this.dexDirection
               axes(ax1); errorbar([160, 320, 640],squeeze(k_pulse_mean(1,i,:)),squeeze(k_pulse_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax2); errorbar([160, 320, 640],squeeze(k_release_mean(1,i,:)),squeeze(k_release_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax3); errorbar([160, 320, 640],squeeze(k_stocastic_mean(1,i,:)),squeeze(k_stocastic_std(1,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax1);plot([0 800],[160, 320, 640].*ones(2,3),'--k','linewidth',2.5);
            axes(ax2);plot([0 800],[160, 320, 640].*ones(2,3),'--k','linewidth',2.5);
            axes(ax3);plot([0 800],[160, 320, 640].*ones(2,3),'--k','linewidth',2.5);

            axes(ax1); legend(direcVal);

            % X0 plot
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,3,1); hold on;
            ax2 = subplot(1,3,2); hold on;
            ax3 = subplot(1,3,3); hold on;
            
            axes(ax1); title('Pulse');
            ylabel('x_0 (mm)');ylim([0 10]);
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release');
            ylim([0 10]); % ylabel('Stiffness (N/m)'); 
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
            
            axes(ax3); title('Stocastic');
            ylim([0 10]); %ylabel('Stiffness (N/m)');
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            

            for i = this.dexDirection
               axes(ax1); errorbar(distVal(this.dexDistance),100*squeeze(x0_pulse_mean(1,i,:)),100*squeeze(x0_pulse_std(1,i,:)),'linewidth',2.5); hold on;
%                axes(ax2); errorbar(distVal(this.dexDistance),squeeze(k_release_mean(1,i,:)),squeeze(k_release_std(1,i,:)),'linewidth',2.5); hold on;
%                axes(ax3); errorbar(distVal(this.dexDistance),squeeze(k_stocastic_mean(1,i,:)),squeeze(k_stocastic_std(1,i,:)),'linewidth',2.5); hold on;
            end
            axes(ax1); legend(direcVal);
            
            %% VAF/Coherence
            figure('Position',[300 314 929 420]);
            ax1 = subplot(1,3,1); hold on;
            ax2 = subplot(1,3,2); hold on;
            ax3 = subplot(1,3,3); hold on;
            
            axes(ax1); title('Pulse');
            ylabel('VAF');ylim([0 100]);
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            
            axes(ax2); title('Release');
            ylim([0 100]); ylabel('VAF'); 
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16); grid on;
            
            axes(ax3); title('Stocastic');
            ylim([0 1]); ylabel('Coherence');
            xlabel('Distance'); xticks(distVal); xlim(xRange);
            set(gca,'fontsize',16);grid on;
            

            for i = this.dexDirection
               axes(ax1); errorbar(distVal(this.dexDistance),squeeze(k_pulse_VAF_mean(1,i,:)),squeeze(k_pulse_VAF_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax2); errorbar(distVal(this.dexDistance),squeeze(k_release_VAF_mean(1,i,:)),squeeze(k_release_VAF_std(1,i,:)),'linewidth',2.5); hold on;
               axes(ax3); errorbar(distVal(this.dexDistance),squeeze(OC_hat_stocastic_mean(1,i,:)),squeeze(OC_hat_stocastic_std(1,i,:)),'linewidth',2.5); hold on;

            end
            axes(ax1); legend(direcVal,'location','south');


                
                
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
            colorVec = {'b','r','g'};
            markerVec = {'o','+','*'};
            figure;
            subj = 1;
            for dir = this.dexDirection
                for dist = this.dexDistance
                    plot(1./squeeze(100*this.x0_hat_pulse(subj,dir,dist,:)),...
                        squeeze(this.k_hat_pulse(subj,dir,dist,:)),...
                        'color',colorVec{dir},...
                        'marker',markerVec{dist},...
                        'LineStyle','none','linewidth',1.5,'markersize',7); hold on;
                end
            end
            ylabel('Stiffness (N/m)'); ylim([0 1500]);
            xlabel('1/x_0 (mm)'); %xlim([0 10]);
            set(gca,'fontsize',16);
            
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
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
%                     subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:size(data,4)
                        clear tmp_data
                        if(~isempty(data{subj,dir,dist,trial,step_Pulse}))
                            tmp_data = data{subj,dir,dist,trial,step_Pulse};
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
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
%                     subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:size(data,4)
                        clear tmp_data
                        if(~isempty(data{subj,dir,dist,trial,step_Pulse}))
                            tmp_data = data{subj,dir,dist,trial,step_Pulse};
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
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
                    subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:size(data,4)
                        clear tmp_data
                        if(~isempty(data{subj,dir,dist,trial,step_Pulse}))
                        tmp_data = data{subj,dir,dist,trial,step_Pulse};
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
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
                    subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,dir,dist,trial,3};
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
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
                    subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,dir,dist,trial,3};
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
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
                    subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,dir,dist,trial,1};
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
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
                    subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,dir,dist,trial,1};
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