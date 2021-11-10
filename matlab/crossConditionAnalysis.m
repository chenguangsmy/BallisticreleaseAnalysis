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
        k_nonNan
        
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
            
            this.plot_postionForce_pulse(data);
%             this.plot_positionForce_stocastic(data);
%             this.plot_postion_release(data);
            
            distVal = [2.5 5 7.5];
%             direcVal = {'front','back','left','right'};
            direcVal = {'F = 15N','F = 20N','F = 25N'};

            
            %% Look at Step Estimates
            k_mean = nanmean(this.k_hat_pulse(:,:,:,1:end),4);
            k_std = nanstd(this.k_hat_pulse(:,:,:,1:end),0,4);
            
            this.k_nonNan = squeeze(sum(~isnan(this.k_hat_pulse),4));
            
            figure; 
            subplot(1,3,1);
            for i = this.dexDirection
                errorbar(distVal(this.dexDistance),squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'linewidth',2.5); hold on;
            end
            title('Pulse'); 
            ylabel('Stiffness (N/m)');ylim([0 1500]);
            xlabel('Distance'); xticks(distVal); xlim([1.5 11]);
            legend(direcVal);
            set(gca,'fontsize',16);
            if(strcmp(subjectType,'spring'))
                plot([0; 10].*ones(2,4),[160, 320, 640, 0].*ones(2,4),'--k');
            end
            
            % Tmp figure for lab presenation
%             figure; 
%             for i = this.dexDirection
%                 errorbar(distVal(this.dexDistance),squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'b','linewidth',3); hold on;
%             end
%             plot(ones(1,4).*[0;5],ones(2,1).*[160, 320, 640,0],'--k','linewidth',2.5);
%             title('Pulse'); 
%             ylabel('Stiffness (N/m)');ylim([0 1000]);
%             xlabel('Distance'); xticks(distVal); xlim([0 5]); set(gca,'Xtick',[]);
% %             legend(direcVal);
%             set(gca,'fontsize',16);
            
            
            
%             for i = this.dexDirection
%                 figure;
%                 for j = this.dexDistance
%                     plot(squeeze(this.k_hat_pulse(1,i,j,:)),'-o'); hold on;
%                 end
%             end
            
            %% Look at stocastic estimates
            k_mean = nanmean(this.k_hat_stocastic(:,:,:,:),4);
            k_std = nanstd(this.k_hat_stocastic(:,:,:,:),0,4);
            
% %             figure;
            subplot(1,3,2);
            for i = this.dexDirection
                errorbar(distVal(this.dexDistance),squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'linewidth',2.5); hold on;
            end
            title('Stocastic');
            ylabel('Stiffness (N/m)');ylim([0 1500]);
            xlabel('Distance'); xticks(distVal); xlim([1.5 11]);
%             legend(direcVal);
            set(gca,'fontsize',16);
            plot([0; 10].*ones(2,4),[160, 320, 640, 0].*ones(2,4),'--k');

            
%             for i = this.dexDirection
%                 figure;
%                 for j = this.dexDistance
%                     plot(squeeze(this.k_hat_stocastic(1,i,j,:))); ylim([0 1000]); hold on;
%                 end
%             end
            
            %% Look at release estimates
            k_mean = nanmean(this.k_hat_release(:,:,:,:),4);
            k_std = nanstd(this.k_hat_release(:,:,:,:),0,4);
            
%             figure; 
            subplot(1,3,3);
            for i = this.dexDirection
                errorbar(distVal(this.dexDistance),squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'linewidth',2.5); hold on;
            end
            title('Release');
            ylabel('Stiffness (N/m)'); ylim([0 1500]);
            xlabel('Distance'); xticks(distVal); xlim([1.5 11]);
%             legend(direcVal);
            set(gca,'fontsize',16);
            plot([0; 10].*ones(2,4),[160, 320, 640, 0].*ones(2,4),'--k');

            
%             for i = this.dexDirection
%                 figure;
%                 for j = this.dexDistance
%                     plot(squeeze(this.k_hat_release(1,i,j,:))); ylim([0 1000]); hold on;
%                 end
%             end
            
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
                        end
                        
                        if(length(Tmp_depMeasures.k_hat_stocastic)~=0)
                            this.k_hat_stocastic(subj,dir,dist,:) = Tmp_depMeasures.k_hat_stocastic;
                        end
                        
                        if(length(Tmp_depMeasures.k_hat_release)~=0)
                            this.k_hat_release(subj,dir,dist,:) = Tmp_depMeasures.k_hat_release;
                        end
                        
                    end
                end
            end
            
        end
        
        function [] = plot_postionForce_pulse(this,data)
            
            step_Pulse = 2;
            
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
                            
                            plot(tmp_data.x(2,dexPulse)); hold on;
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
                    
                    subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:size(data,4)
                        clear tmp_data
                        if(~isempty(data{subj,dir,dist,trial,step_Pulse}))
                            tmp_data = data{subj,dir,dist,trial,step_Pulse};
                            dexPulse = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))];
                            
                            plot(tmp_data.v(2,dexPulse)); hold on;
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
                        plot(tmp_data.f(2,dexPulse)); hold on;
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
                        dexRelease = [min(find(tmp_data.ts== 4)) : max(find(tmp_data.ts == 5))];
                        plot(tmp_data.x(2,dexRelease)); hold on;
                        end
                    end
                    count = count+1;
                    xlabel('Time'); ylabel('Postion (m)'); set(gca,'fontsize',16);
                end
            end

            

        end

    end
end