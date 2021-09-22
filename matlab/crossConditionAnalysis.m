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
        
    end
    
    methods
        function [this] = crossConditionAnalysis(data,dexSubject,dexDirection,dexDistance)
            
            % Define index to use (will replace soon as passed parameter)
            this.dexSubject = dexSubject;
            this.dexDirection = dexDirection;
            this.dexDistance = dexDistance;
            this.sfrq = 500;
            
            %Compute all stiffness estimates
            this.get_depMeasures(data);
            
            this.plot_postionForce_pulse(data);
            this.plot_positionForce_stocastic(data);
            this.plot_postion_release(data);
            
            
            %% Look at Step Estimates
            k_mean = nanmean(this.k_hat_pulse(:,:,:,5:end),4);
            k_std = nanstd(this.k_hat_pulse(:,:,:,5:end),0,4);
            
            figure;
            for i = this.dexDirection
                errorbar(this.dexDistance,squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'linewidth',2.5); hold on;
            end
            
            for i = this.dexDirection
                figure;
                for j = this.dexDistance
                    plot(squeeze(this.k_hat_pulse(1,i,j,:))); hold on;
                end
            end
            
            %% Look at stocastic estimates
            k_mean = nanmean(this.k_hat_stocastic(:,:,:,:),4);
            k_std = nanstd(this.k_hat_stocastic(:,:,:,:),0,4);
            
            figure;
            for i = this.dexDirection
                errorbar(this.dexDistance,squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'linewidth',2.5); hold on;
            end
            
            for i = this.dexDirection
                figure;
                for j = this.dexDistance
                    plot(squeeze(this.k_hat_stocastic(1,i,j,:))); ylim([0 1000]); hold on;
                end
            end
            
            %% Look at release estimates
            k_mean = nanmean(this.k_hat_release(:,:,:,:),4);
            k_std = nanstd(this.k_hat_release(:,:,:,:),0,4);
            
            figure;
            for i = this.dexDirection
                errorbar(this.dexDistance,squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'linewidth',2.5); hold on;
            end
            
            for i = this.dexDirection
                figure;
                for j = this.dexDistance
                    plot(squeeze(this.k_hat_release(1,i,j,:))); ylim([0 1000]); hold on;
                end
            end
            

            % Make ANOVA Plots

%             
            
            % Run ANOVA between subjects
%             this.

            disp('test');
            
        end
        
        function [] = get_depMeasures(this,data)
            
            sizeData = size(data);
            
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
                        Tmp_depMeasures = crossTrialAnalysis(dataCross,this.sfrq);
                        this.k_hat_pulse(subj,dir,dist,:) = Tmp_depMeasures.k_hat_pulse;
                        this.k_hat_stocastic(subj,dir,dist,:) = Tmp_depMeasures.k_hat_stocastic;
%                         this.k_hat_release(subj,dir,dist,:) = Tmp_depMeasures.k_hat_release;
                        
                    end
                end
            end
            
        end
        
        function [] = plot_postionForce_pulse(this,data)
            
            figure;
            sgtitle('Pulse Preturbation');
            subj = 1;
            count = 1;
            for dir = this.dexDirection
                for dist = this.dexDistance
                    
                    subplot(length(this.dexDirection),length(this.dexDistance),count);
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,dir,dist,trial,2};
                        dexPulse = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))];
                        plot(tmp_data.x(2,dexPulse)); hold on;
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
                    for trial = 2:15
                        clear tmp_data
                        tmp_data = data{subj,dir,dist,trial,2};
                        dexPulse = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))];
                        plot(tmp_data.f(2,dexPulse)); hold on;
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
                            dexStoc = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))];
                            plot(tmp_data.x(2,dexStoc)); hold on;
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
                        dexStod = [min(find(tmp_data.Fp(2,:)~=0)) : max(find(tmp_data.Fp(2,:)~=0))];
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