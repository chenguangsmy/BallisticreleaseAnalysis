classdef crossConditionAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.
    
    properties
        
        dexSubject
        dexDirection
        dexDistance
        depMeasures % Dependentmeasures (three stiffness estimates defined in crossTrialAnalysis)
        
    end
    
    methods
        function [this] = crossConditionAnalysis(data,dexSubject,dexDirection,dexDistance)
            
            % Define index to use (will replace soon as passed parameter)
            this.dexSubject = dexSubject;
            this.dexDirection = dexDirection;
            this.dexDistance = dexDistance;
            
            % Compute all stiffness estimates
            this.get_depMeasures(data);
            
            % Make ANOVA Plots
%             this.plot_postion();
%             this.plot_force();
            
            % Run ANOVA between subjects
%             this.
            
        end
        
        function [] = get_depMeasures(this,data)
            
            sizeData = size(data);
            
            for subj = this.dexSubject
                for dir = this.dexDirection
                    for dist = this.dexDistance
                        this.depMeasures{subj,dir,dist} = crossTrialAnalysis(data{subj,dir,dist,:,:});
                    end
                end
            end
            
        end
        
    end
end