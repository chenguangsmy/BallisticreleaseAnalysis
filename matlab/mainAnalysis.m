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
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss2818_2842_adp.mat')
            
            % Hack this data type to be able to run analysis
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3003_3005.mat')
% 
%             for width = 1:3
%                 for height = 1:4
%                     for trial = 1:10
%                         datatmp{1,width,height,trial,2} = data{width,height}{trial};
%                     end
%                 end
%             end
%             
%             datatmp{1,width,height,trial,1}=[];
%             datatmp{1,width,height,trial,3}=[];
% 
%             clear data
%             data = datatmp;
% 
%             dexSubject = 1;
%             dexDirection = 1:3;
%             dexDistance = 1:4;
            
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/PertParamSelection_gaussian_cg.mat');
%             data = Data;
%             clear Data
%             
%             mag = 3;
%             duration = 1;
%             dist = 1:3;
%             posNeg = 1;
%             for trial = 1:15
%                 datatmp{1,1,1,trial,2} = data{mag,duration,dist,posNeg}{trial};
%             end
%             
%             datatmp{1,1,1,trial,1}=[];
%             datatmp{1,1,1,trial,3}=[];
% 
%             clear data
%             data = datatmp;
% 
%             dexSubject = 1;
%             dexDirection = 1;
%             dexDistance = 1:3;

%             % Gaussian test with springs
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/SpringGaussian3.mat');
%             data = Data;
%             clear Data
%             
%             mag = 2;
%             duration = 1;
%             for posNeg = 1:2
%                 for stiffness = 1:4
%                     for trial = 1:15
%                         datatmp{1,stiffness,posNeg,trial,2} = data{mag,duration,stiffness,posNeg}{trial};
%                     end
%                 end
%             end
%             
%             datatmp{1,1,1,trial,1}=[];
%             datatmp{1,1,1,trial,3}=[];
% 
%             clear data
%             data = datatmp;
% 
%             dexSubject = 1;
%             dexDirection = 1:4;
%             dexDistance = 1:2;

% %             James Testing 9/22/21
%             load('/Users/jhermus/Downloads/ss2872_2876.mat');
%             dataTmp = data{1,4,2,4,2};
%             figure;
%             ax1 = subplot(3,1,1); plot(dataTmp.t, dataTmp.x(2,:));
%             ax2 = subplot(3,1,2); plot(dataTmp.t, dataTmp.f(2,:));
%             ax3 = subplot(3,1,3); plot(dataTmp.t, dataTmp.Fp(2,:));
%             linkaxes([ax1,ax2,ax3],'x');
%             
%             dexSubject = 1;
%             dexDirection = 1:4;
%             dexDistance = 1:3;
            
            %% Visiting Pittsburgh Test 11/3/2021
            % Test from James Wednesday
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3307_3314.mat');
%             data_human = data;
%             clear data
%             dexSubject = 1; % [James]
%             dexForce = 1:3; % [15N, 20N, 25N]
%             dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]

            % Spring Test
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3334_3344.mat');
            data_spring = data;
            clear data
            dexSpring = 1; % [336 N/m]
            dexForce_spring = 3; % [15N, 20N, 25N]
            dexDistance_spring = 1:3; % [will vary based on spring constant]
%             Spring [4.7, 6.25, 7.81]
          
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3345.mat');
%             data_robotOnly = data;
%             clear data
%             dexRob = 1; % 
%             dexForce_Rob = 1; 
%             dexDistance_Rob = 1;
             
            % Run analysis
%             depMeasures_human = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance);
            depMeasures_spring = crossConditionAnalysis(data_spring,dexSpring,dexForce_spring,dexDistance_spring);
%             depMeasures_rob = crossConditionAnalysis(data_robotOnly,1,1,1);

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