% Main function to run data analysis
classdef mainAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.
    
    properties
        
        
    end
    
    methods
        function [this] = mainAnalysis()
            
            dbstop if error
            addpath('func/');

%           load('/Users/jhermus/Downloads/ss2818_2842.mat');
            % Adapted filed
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss2818_2842_adp.mat')

            % James Testing 9/22/21
%             load('/Users/jhermus/Downloads/ss2872_2876.mat');
%             dataTmp = data{1,4,2,4,2};
%             figure;
%             ax1 = subplot(3,1,1); plot(dataTmp.t, dataTmp.x(2,:));
%             ax2 = subplot(3,1,2); plot(dataTmp.t, dataTmp.f(2,:));
%             ax3 = subplot(3,1,3); plot(dataTmp.t, dataTmp.Fp(2,:));
%             linkaxes([ax1,ax2,ax3],'x');
            
            dexSubject = 1;
            dexDirection = 1:4;
            dexDistance = 1:3;
            
            % Import all human dependent measures
            depMeasures_human = crossConditionAnalysis(data, dexSubject, dexDirection, dexDistance);
            
%             % Import all spring dependent measures
%             dexSpring = 1;
%             depMeasures_spring = crossConditionAnalysis(data,dexSpring,dexDirection,dexDistance);

            % Make ANOVA Plots across human and spring trials
%             this.plot_postion();
%             this.plot_force();
%             this.plot_stiffness();
            
            % Run ANOVA between subjects
%             this.get_ANOVA_human_v_spring();
            disp('test');
            
        end
        
        
    end
end