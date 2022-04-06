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

            % Spring Test
%             load('data/ss3334_3344_alter.mat');
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
%             load('data/ss3345.mat');
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
%             load('data/prelimData_6subj_fine.mat');
%             data_human = reshape(data(:,1,:,:,:,:),6,3,3,15,3); % Elliminate direction
%             clear data
%             dexSubject = 1;%:6;
%             dexForce = 1%:3; % [15N, 20N, 25N]
%             dexDistance = 1%:3; % [2.5cm, 5cm, 7.5cm]
%             depMeasures_human = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
         
%             % Human 12N 100ms to late
%             load('data/ss3818_3828.mat');  
%             data_spring = reshape(data(1,1,:,:,:,:),1,3,3,5,13);
%             clear data
%             dexSubject = 1;
%             dexForce = 1:3;%:3; % [15N, 20N, 25N]
%             dexDistance = 1:3;%:3; % [2.5cm, 5cm, 7.5cm]
%             depMeasures_human = crossConditionAnalysis(data_spring, dexSubject, dexForce, dexDistance,'human');
%             save('/Users/jhermus/Desktop/Spring_12N_100ms_ss3818_3828_latespaceing.mat');
            
            %% Meeting Motion Pulse Analysis
            
%             % Spring 12N 100ms
%             load('data/ss3803_3812.mat');  
%             data_spring = data; % Elliminate direction
%             clear data
%             dexSubject = 1;
%             dexForce = 1:3;%:3; % [15N, 20N, 25N]
%             dexDistance = 1:3;%:3; % [2.5cm, 5cm, 7.5cm]
%             depMeasures_human = crossConditionAnalysis(data_spring, dexSubject, dexForce, dexDistance,'human');
%             save('/Users/jhermus/Desktop/Spring_12N_100ms_ss3803_3812.mat');

%             % Human 12N 100 ms 
%             load('data/ss3896_3905.mat');
%             data_human = reshape(data(1,1,:,:,:,:),1,3,3,10,8); % Elliminate direction
%             clear data
%             dexSubject = 1;
%             dexForce = 1:3 ; % [15N, 20N, 25N]
%             dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
%             depMeasures_human = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human'); 
%             save('/Users/jhermus/Desktop/Human_12N_100ms_ss3896_3905.mat');
            
%             % Springs 12N 200 ms
%             load('data/ss3873_3884.mat');
%             data_spring = data; % Elliminate direction
%             clear data
%             dexSubject = 1;
%             dexForce = 1:3;%:3; % [15N, 20N, 25N]
%             dexDistance = 1:3;%:3; % [2.5cm, 5cm, 7.5cm]
%             depMeasures_human = crossConditionAnalysis(data_spring, dexSubject, dexForce, dexDistance,'human');
% %             save('/Users/jhermus/Desktop/Spring_12N_200ms_ss3873_3884.mat');

             % Human 12N 200 ms 
%             load('data/ss3913_3921.mat');
%             data_human = reshape(data(1,1,:,:,:,:),1,3,3,10,8); % Elliminate direction
%             clear data
%             dexSubject = 1;
%             dexForce = 1:3 ; % [15N, 20N, 25N]
%             dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
%             depMeasures_human = crossConditionAnalysis(data_human, dexSubject, dexForce, dexDistance,'human');
% %             save('/Users/jhermus/Desktop/Human_12N_200ms_ss3913_3921.mat');

%             %% New data after Motion Pulse meeting (2/14/2022)
%             % Force data better aligned not pulse before release
%             % Springs 12N 200 ms
%             load('data/ss3925_3937.mat');
%             data_spring = reshape(data(1,:,:,:,1:6),1,3,3,15,6); % Elliminate direction
%             clear data
%             dexSubject = 1;
%             dexForce = 1:3; % [15N, 20N, 25N]
%             dexDistance = 1:3; % [2.5cm, 5cm, 7.5cm]
%             depMeasures_human = crossConditionAnalysis(data_spring, dexSubject, dexForce, dexDistance,'human');
% %             save('/Users/jhermus/Desktop/Spring_12N_200ms_ss3925_3937.mat');  

            % Check new postions with stocstic method
            load('data/ss4010_4013.mat');
            for force = 2
                for posture = 1:2
                    trial = 1;
                    X = data(1,1,force,posture,1,3);
                    X = X{1};
                    dex = find(X.ts == 4);
                    dex = dex(1000:end);
                    this.get_Z_spectral(X.Fp(1,dex),X.Fp(2,dex),X.x(1,dex),X.x(2,dex));
                end
            end
            %% Run from here to plot

%             save('/Users/jhermus/Desktop/test2022-2-4NoConstraint_noPrior.mat');
%             load('/Users/jhermus/Desktop/test2022-2-4.mat');

%               save('/Users/jhermus/Desktop/test2022-4-1.mat');
             
%             load('/Users/jhermus/Desktop/Spring_12N_100ms_ss3803_3812.mat');
%             load('/Users/jhermus/Desktop/Human_12N_100ms_ss3896_3905.mat');
%             load('/Users/jhermus/Desktop/Spring_12N_200ms_ss3873_3884.mat');
%             load('/Users/jhermus/Desktop/Human_12N_200ms_ss3913_3921.mat');

%             load('/Users/jhermus/Desktop/Spring_12N_100ms_ss3818_3828_latespaceing.mat');
            pathh = 'prelimSubjectsRawPlots/images/';

            % Make plots Spring
%             depMeasures_spring.get_mainPlot('spring',1,pathh);
             
            % Make plots subjects
            for i = 1
               depMeasures_human.get_mainPlot('human',i,pathh);
            end
                        
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

% sfrq = 10000;
% dt = 1/sfrq;
% fs_pulse = 10;
% bw = 0.5; % 60%
% trunk = 60; % 40 dB truncated
% tc = gauspuls('cutoff',fs_pulse,bw,[],-trunk); 
% t = -tc : dt : tc; 
% [yi,yq,ye] = gauspuls(t,fs_pulse,bw); 
% 
% figure; plot(t,yi,t,yq,[t,t],[ye,-ye],'linewidth',2.5);
% legend('Inphase','Quadrature','Envelope');