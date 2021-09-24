classdef crossTrialAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a 
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.  
   
    properties
        
        k_hat_step
        k_hat_stocastic
        k_hat_release
        
    end
    
    methods
        function [this] = crossTrialAnalysis(data)
            
            % Look at data
%             this.plot_postion();
%             this.plot_force();
            
            % Estimate stiffnesses
            this.get_k_hat_step();
%             this.get_k_hat_stocastic();
%             this.get_k_hat_release();
            
        end
        
        function [] = get_k_hat_step(this,data)
            
            test = 1;
            
        end
        
        function [] = get_k_hat_stocastic(this)
            
            
            
        end
        
        function [] = get_k_hat_release(this)
            
            
            
        end
        
    end
end


