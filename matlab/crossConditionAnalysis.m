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
            
            %% Look at Step Estimates
            k_mean = nanmean(this.k_hat_pulse(:,:,:,5:end),4);
            k_std = nanstd(this.k_hat_pulse(:,:,:,5:end),0,4);
            
            figure;
            for i = 1:4
                errorbar([1,2,3],squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:))); hold on;
            end
            
            for i = 1:4
                figure;
                for j = 1:3
                    plot(squeeze(this.k_hat_pulse(1,i,j,:))); hold on;
                end
            end
            
            %% Loop at stocastic estimates
            k_mean = nanmean(this.k_hat_stocastic(:,:,:,:),4);
            k_std = nanstd(this.k_hat_pulse(:,:,:,:),0,4);
            
            
            figure;
            for i = 1:4
                errorbar([1,2,3],squeeze(k_mean(1,i,:)),squeeze(k_std(1,i,:)),'linewidth',2.5); hold on;
            end
            
            for i = 1:4
                figure;
                for j = 1:3
                    plot(squeeze(this.k_hat_stocastic(1,i,j,:))); ylim([0 1000]); hold on;
                end
            end
            

            % Make ANOVA Plots
%             this.plot_postion();
%             this.plot_force();
            
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
%                         k_hat_release[subj,dir,dist,:] = Tmp_depMeasure.k_hat_release;
                        
                    end
                end
            end
            
        end
        
    end
end