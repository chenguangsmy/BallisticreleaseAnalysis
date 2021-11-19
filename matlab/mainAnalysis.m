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
            % Test fromt Chenguang both directions
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3353_3417.mat');
            data_human = reshape(data(1,1,:,:,:,:),1,3,3,15,3);
            clear data
            dexSubject = 1; % [Chenguang]
            dexForce = 1:3; % [15N, 20N, 25N]
            dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
            depMeasures_human2 = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
            
            
            % Test from James Wednesday
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3307_3314.mat');
            data_human = data;
            clear data
            dexSubject = 1; % [James]
            dexForce = 1:3; % [15N, 20N, 25N]
            dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
            
            depMeasures_human1 = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
             

            % Spring Test
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3334_3344.mat');
            data_spring = data;
            clear data
            dexSpring = 1; % [336 N/m]
            dexForce_spring = 1; % [15N, 20N, 25N]
            dexDistance_spring = 1:3; % [will vary based on spring constant]
%             Spring [4.7, 6.25, 7.81]
            depMeasures_spring = crossConditionAnalysis(data_spring,dexSpring,dexForce_spring,dexDistance_spring,'spring');
%             depMeasures_spring.k_nonNan
%             depMeasures_human1.k_nonNan
%             depMeasures_human2.k_nonNan
%             depMeasures_spring.k_nonNan
          
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3345.mat');
%             data_robotOnly = data;
%             clear data
%             dexRob = 1; % 
%             dexForce_Rob = 1; 
%             dexDistance_Rob = 1;
             
            % Run analysis
%             depMeasures_human = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance);
%             depMeasures_spring = crossConditionAnalysis(data_spring,dexSpring,dexForce_spring,dexDistance_spring);
%             depMeasures_rob = crossConditionAnalysis(data_robotOnly,1,1,1);

            % Spring Test - diffpreturbation param
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3318_3319.mat');
%             data_spring = data;
%             clear data
%             dexSpring = 1; % [336 N/m]
%             dexForce_spring = 1; % [15N, 20N, 25N]
%             dexDistance_spring = 1:2; % [will vary based on spring constant]
%             depMeasures_spring = crossConditionAnalysis(data_spring,dexSpring,dexForce_spring,dexDistance_spring,'spring');

            % Check out 2N stocastic
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/SpringGaussian_8N_format.mat');
            data = reshape(Data(1,1,1,:,:,:),1,1,2,15,3);
            crossConditionAnalysis(data, 1, 1, 1:2,'spring');
            
            % Check out 3N stocastic
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/SpringStoc_3N.mat');
            data = reshape(Data(1,1,1,1,:,:),1,1,1,15,3);
            crossConditionAnalysis(data, 1, 1, 1,'spring');

            
            % Test from Chenguang and subject
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3431_3440.mat');
%             clear data
            dexSubject = 1; % [first subject]
            dexForce = 1; % [15N, 20N, 25N]
            dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
            
            depMeasures_TestHuman_Overshoot1 = crossConditionAnalysis(reshape(data(1,1,1,:,:,:),1,1,3,15,3), dexSubject, dexForce, dexDistance,'human');
            depMeasures_TestHuman_Overshoot2 = crossConditionAnalysis(reshape(data(1,1,2,:,:,:),1,1,3,15,3), dexSubject, dexForce, dexDistance,'human');
%             depMeasures_TestHuman_Overshoot1.k_nonNan
%             depMeasures_TestHuman_Overshoot2.k_nonNan
%             depMeasures_TestHuman_Overshoot1 = crossConditionAnalysis(reshape(data(1,1,1,:,:,:),1,1,3,15,3), dexSubject, dexForce, dexDistance,'human');
%             depMeasures_TestHuman_Overshoot1 = crossConditionAnalysis(reshape(data(1,2,1,:,:,:),1,1,3,15,3), dexSubject, dexForce, dexDistance,'human');


%             depMeasures_Testhuman1 = crossConditionAnalysis(reshape(data(2,1,:,:,:,:),1,1,3,15,3), 1, dexForce, dexDistance,'human');
%             depMeasures_Testhuman1 = crossConditionAnalysis(reshape(data(3,1,:,:,:,:),1,1,3,15,3), 1, dexForce, dexDistance,'human');
%             depMeasures_Testhuman1 = crossConditionAnalysis(reshape(data(4,1,:,:,:,:),1,1,3,15,3), 1, dexForce, dexDistance,'human');

            % New Visual Feedback (naive subject)
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3477_3479.mat');
            dexSubject = 1; % [first subject]
            dexForce = 1; % [15N, 20N, 25N]
            dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
            
            depMeasures_TestHuman = crossConditionAnalysis(data, dexSubject, dexForce, dexDistance,'human');
            depMeasures_TestHuman.k_nonNan
            
            % New Visual/spesifict instruction Feedback (naive subject)
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/ss3491_3493.mat');
            data_human = reshape(data(1,1,1,:,:,:),1,1,3,15,3);
            dexSubject = 1; % [first subject]
            dexForce = 1; % [15N, 20N, 25N]
            dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
            
            depMeasures_TestHuman = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
            depMeasures_TestHuman.k_nonNan

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