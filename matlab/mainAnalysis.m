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


%             % Spring Test
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3334_3344.mat');
%             dexForce_spring = 1:3; % [15N, 20N, 25N]
%             data_spring = data;
%             for force = dexForce_spring
%                 for stiff = 1:3
%                     for trial = 1:15
%                         for pret = 1:3
%                             data_spring{1,force,stiff+1,trial,pret} = data{1,force,stiff,trial,pret};
%                         end
%                     end
%                 end
%             end
%             clear data
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3345.mat');
%             for force = dexForce_spring
%                 for trial = 1:15
%                     for pret = 1:3
%                     data_spring{1,force,1,trial,pret} = data{1,1,1,trial,pret};
%                     end
%                 end
%             end
%             clear data
%             dexStiff_spring = 1:4; 
%             depMeasures_spring = crossConditionAnalysis(data_spring,1,dexForce_spring,dexStiff_spring,'spring');

            % Combined Many subject Tests (12/10/2021)
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/prelimData_6subj_fine.mat');
            data_human = reshape(data(6,1,:,:,:,:),1,3,3,15,3); % Elliminate direction
            clear data
            dexSubject = 1;
            dexForce = 1:3; % [15N, 20N, 25N]
            dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
            depMeasures_human = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
            
                        
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