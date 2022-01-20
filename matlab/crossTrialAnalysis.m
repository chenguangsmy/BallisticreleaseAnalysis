classdef crossTrialAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a 
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.  
   
    properties
        
        trial_norm
        trial_catch
        trial_stocastic   
        sfrq
        k_r
        x_r
        m
        x_handelStart 
        x0_hat_pulse

    end
    
    methods
        function [this] = crossTrialAnalysis(data,sfrq,f_target,distVal,subjectType)
            
            this.sfrq = sfrq;
            this.k_r = 300;
            this.x_handelStart =  0.481;
            this.x_r = this.x_handelStart-f_target/this.k_r;

            
            if(strcmp(subjectType,'human'))
                this.m = 2.15;
            elseif(strcmp(subjectType,'spring'))
                this.m = 1.15;
            else
                error('Spesify subjectType');
            end
            
            % Look at data
%             this.plot_postion();
%             this.plot_force();
            
            % Estimate stiffnesses
            this.trial_norm = this.get_trial_norm(data,f_target,distVal,subjectType);
            this.trial_catch = this.get_trial_catch(data,f_target,distVal,subjectType);
            this.trial_stocastic = this.get_trial_stocastic(data,f_target,distVal,subjectType);
            
        end
        
        %% Call general functions to create data structures

        
        function [trial_norm] = get_trial_norm(this,data,f_target,distVal,subjectType)
            
            % get_trial_norm createst a data strucutre with an estimate of pulse
            % and release
        
            trialType = 1;
            [trial_norm.est_release] = this.get_k_hat_release(data,subjectType,trialType,distVal);
            
            
        end
        
        function [trial_catch] = get_trial_catch(this,data,f_target,distVal,subjectType)
            
            % get_est_trial createst a data strucutre with an estimate of pulse
            % and release using data from the 
            
            trialType = 2;
            [trial_catch.est_release] = this.get_k_hat_release(data,subjectType,trialType,distVal);
            [trial_catch.est_pulse] = this.get_k_hat_pulse(data,f_target,subjectType,trialType);
            
        end
        
        function [trial_stocastic] = get_trial_stocastic(this,data,f_target,distVal,subjectType)
            
            % get_est_trial createst a data strucutre with an estimate of
            % from the stocastic method
            
            trialType = 3;
            [trial_stocastic.est_release] = this.get_k_hat_release(data,subjectType,trialType,distVal);
            [trial_stocastic.est_stocastic] = get_k_hat_stocastic(this,data,trialType);
            
        end
        
        function [est_pulse] = get_k_hat_pulse(this,data,f_target,subjectType,trialType)
            
            % Exclude empty cells (FIX Later)
            [N_trial] = size(data,1);
            count = 1;
            for trial = 1:N_trial-1
                if(~isempty(data{trial,trialType}))
                    if( sum(data{trial,trialType}.Fp(2,:)~=0) ~= 0 )
                        if(sum(data{trial,trialType}.Fp(2,1:30)) == 0) % find pulse to close to start
                        dexTrialNonzero(count) = trial;
                                        count = count+1;
                        end
                    end
                end
            end

            for i = 1:length(dexTrialNonzero)
                clear tmpData
                tmpData = data{dexTrialNonzero(i),trialType};
                [est_pulse{i}] = get_singleTrial_k_hat_pulse_con(this,tmpData.f(2,:),...
                                                 tmpData.Fp(2,:),...
                                                 tmpData.x(2,:),...
                                                 tmpData.ts,...
                                                 f_target,subjectType);
            end
                        
        end
        
        function [est_stocastic] = get_k_hat_stocastic(this,data,trialType)
            
            % Exclude empty cells (FIX Later)
            [N_trial] = size(data,1);
            dexTrialNonzero = [];
            count = 1;
            for trial = 1:N_trial
                if(~isempty(data{trial,trialType}))
                    dexTrialNonzero(count) = trial;
                end
                count = count+1;
            end
            
            if(~isempty(dexTrialNonzero))

                for i = 1:length(dexTrialNonzero)
                    clear tmpData
                    tmpData = data{dexTrialNonzero(i),trialType};
                    [est_stocastic{i}] = get_singleTrial_k_hat_stocastic(this,tmpData.Fp(2,:),...
                        tmpData.x(2,:),...
                        tmpData.ts);
                end
                
            end
            
        end
        
        function [est_release] = get_k_hat_release(this,data,subjectType,trialType,distVal)
               
            % Exclude empty cells (FIX Later)
            dexTrialNonzero = [];
            [N_trial] = size(data,1);
            count = 1;
            for trial = 1:N_trial % Exclude first trial
                if(~isempty(data{trial,trialType}))
                    dexTrialNonzero(count) = trial;
                end
                count = count+1;
            end
            
            if(~isempty(dexTrialNonzero))
                
                count = 1;
                for i = dexTrialNonzero
                    clear tmpData
                    tmpData = data{i,trialType};
                    est_release{count} = this.get_singleTrial_k_hat_release(tmpData.x(2,:),tmpData.ts,subjectType,distVal);
                    count = count + 1;
                end
                
            end
            
        end
        
        function [est_pulse] = get_singleTrial_k_hat_pulse(this,f,Fp,x,ts,f_target,subjectType)
            
            % Take average of last 15 measuremnts before the pulse ends
            t = 0:1/this.sfrq:(length(x)*(1/this.sfrq))-1/this.sfrq;

                % Estimate static and sliding friction value
                
                % New crop data to use only up till peak
                dexPulseStart = min(find(Fp~=0))+30;
                dexPulseEnd = max(find(Fp~=0));
                dexEndState4 = max(find(ts==4))-50;
                
                cf = 20; % cutoff freqnency
                [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
                x = filtfilt(b,a,x); % apply fitler
                                
                x_dot = this.sfrq*diff(x);
                x_dot(end+1) = x_dot(end);
                x_ddot = this.sfrq*diff(x_dot);
                        
                [pks,locs,w,p] = findpeaks(-x_ddot(dexPulseStart:dexEndState4));
                locs = locs + dexPulseStart;
                
                [B,I] = sort(p,'descend');
                locFound = sort(locs(I(1:2)));
                
                dexRange = locFound(1):locFound(2);
                
%                 figure('Position',[210 305 560 420]);
%                 ax1 = subplot(4,1,1); plot(x); hold on;
%                 plot(locFound,x(locFound),'o');hold on;
% 
%                 ax2 = subplot(4,1,2); plot(x_dot); hold on;
%                 plot(locFound,x_dot(locFound),'o');hold on;
% 
%                 ax3 = subplot(4,1,3); plot(x_ddot);hold on;
%                 plot(locFound,x_ddot(locFound),'o');hold on; 
% 
%                 ax4 = subplot(4,1,4); plot(Fp);
%                 linkaxes([ax1,ax2,ax3,ax4],'x');xlim([dexPulseStart, dexEndState4]);
                
                if(isempty(dexRange))
                    error('Cutting failed: Look for peak finding error');
                end

%                 figure; plot(x_dot(dexPulseStart:dexPulseEnd+extraDex)); hold on;
%                 plot(x_dot(dexRange),'o');

                x_old = x;
                x_dot_old = x_dot;
                x_ddot_old = x_ddot;
                Fp_old = Fp;
                x = x(dexRange);
                x_dot = x_dot(dexRange);
                x_ddot = x_ddot(dexRange);
                Fp = Fp(dexRange);
                
                % Robot info
                k_r = this.k_r;
                x_r = this.x_r; %(Important too check)
                
                % Decision Variables
                f_s = 1; % Static friction
                f_d = 1; % Dynamics friction
                k_h = 300;
                x_h = f(1)/k_h + x(1);
                B = 1;
                X0 = [k_h,x_h,B];
                
                % 300/(2.75*2*pi)^2
                
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb = [0,-2,0];
                ub = [5000,2,500];
                options = optimoptions('fmincon','Display','off');
                                
%                 options = optimoptions('fmincon','TolFun', 1e-8, 'MaxIter', 10000, ...
%                                        'MaxFunEvals', 100000, 'Display', 'off' , ...
%                                        'DiffMinChange', 0.001, 'Algorithm', 'sqp');
                [X_hat,FVAL,EXITFLAG,OUTPUT] = fmincon(@(X_hat)this.costFunc(x,x_dot,x_ddot,x_r,Fp,X_hat),X0,A,b,Aeq,beq,lb,ub,@(x)this.mycon,options);
                
                t_plot = t(dexRange)-t(dexRange(1));
                t_old = t-t(dexRange(1));

                f_threshold = (abs(x_dot)> 3e-3) | (abs(x_ddot) > 0.4);
                f_threshold_dex = find(diff(f_threshold)>0);
                                    
                if(EXITFLAG >= 1) % Local minimum successfully found
                    %                 m = X_hat(1);
                    %                 f_s = X_hat(2);
                    %                 f_d = X_hat(1);
                    k_hat = X_hat(1);
                    x_h = X_hat(2)-this.x_handelStart;
                    B = X_hat(3);
                    
                    [cost,x_hat,x_dot_hat] = this.costFunc(x,x_dot,x_ddot,x_r,Fp,X_hat);
                                        
                    [VAF] = get_VAF(this,x_dot,x_dot_hat);
        
                    colorVec = {[0.9290, 0.6940, 0.1250],...
                        [0, 0.4470, 0.7410],...
                        [0.8500, 0.3250, 0.0980]};
                    
%                     figure('Position',[771 305 560 420]);
%                     
%                     ax1 = subplot(4,1,1);
%                     plot(t_old,x_old,'color',colorVec{1},'linewidth',2);hold on;
%                     plot(t_plot,x,'linewidth',3,'color',colorVec{2}); hold on;
%                     plot(t_plot,x_hat,'.','markersize',10,'color',colorVec{3});
%                     title(['k_{hat} = ',num2str(k_hat),...
%                         ', x_{hat} = ',num2str(x_h),...
%                         ', B = ',num2str(B)]);
%                     xlabel('Time (s)'); ylabel('Postion (m)');
%                     set(gca,'fontsize',16);
%                     
%                     ax2 = subplot(4,1,2);
%                     plot(t_old,x_dot_old,'color',colorVec{1},'linewidth',2); hold on;
%                     plot(t_plot,x_dot,'linewidth',3,'color',colorVec{2}); hold on;
%                     plot(t_plot,x_dot_hat,'.','markersize',10,'color',colorVec{3});
%                     xlabel('Time (s)'); ylabel('Velocity (m/s)');
%                     title(f_target);
%                     set(gca,'fontsize',16);
% 
%                     ax3 = subplot(4,1,3); plot(t_old,x_ddot_old,'color',colorVec{1},'linewidth',2);hold on;
%                     plot(t_old(locFound),x_ddot_old(locFound),'ok');hold on;
%                     xlabel('Time (s)'); ylabel('Acc (m/s^2)'); set(gca,'fontsize',16);
% 
%                     
%                     ax4 = subplot(4,1,4); plot(t_old,Fp_old,'color',colorVec{1},'linewidth',2);
%                     xlabel('Time (s)'); ylabel('Force Pret (N)'); set(gca,'fontsize',16);
%                     linkaxes([ax1,ax2,ax3,ax4],'x'); xlim([t_old(dexPulseStart), t_old(dexEndState4-300)]);


                else % Did not converge
                     k_hat = NaN;
                     x_h = NaN;
                     B = NaN;
                     VAF = NaN;
                     x_dot_hat = NaN;
                     
%                     figure;
%                     
%                     subplot(2,1,1);
%                     plot(t_plot,x,'linewidth',3); hold on;
%                     title('Failed to converge');
%                     xlabel('Time (s)'); ylabel('Postion (m)');
%                     set(gca,'fontsize',16);
%                     
%                     subplot(2,1,2);
%                     plot(t_plot,x_dot,'linewidth',3); hold on;
%                     plot(t_plot,f_threshold*max(x_dot),':k','linewidth',3);
%                     xlabel('Time (s)'); ylabel('Velocity (m/s)');
%                     title(f_target);
%                     set(gca,'fontsize',16);
                    
                end
                
                
                
%                     figure(1);
%                     subplot(2,1,1);
%                     plot(t_plot,x,'-b'); hold on;
%                     if(EXITFLAG == 1)
%                         plot(t_plot,x_hat,':r');
%                     end
%                     xlabel('Time (s)'); ylabel('Postion (m)');
%                     set(gca,'fontsize',16);
%                     
%                     subplot(2,1,2);
%                     plot(t_plot,x_dot,'-b'); hold on;
%                     if(EXITFLAG == 1)
%                         plot(t_plot,x_dot_hat,':r');
%                         plot(t_plot(f_threshold_dex),0,'.g','markersize',10);
%                     else
%                         plot(t_plot(f_threshold_dex),0,'.k','markersize',10);
%                     end
%                     xlabel('Time (s)'); ylabel('Velocity (m/s)');
%                     set(gca,'fontsize',16); grid on;
                

                est_pulse.k_hat = k_hat;
                est_pulse.b_hat = B;
                est_pulse.x_h = x_h;
                est_pulse.x_dot = x_dot;
                est_pulse.t = t_plot;
                est_pulse.dexRange = dexRange;
                est_pulse.x_dot_hat = x_dot_hat;
                est_pulse.VAF = VAF;
                
%                 disp('test');
                
            
        end
        
        function [est_pulse] = get_singleTrial_k_hat_pulse_con(this,f,Fp,x,ts,f_target,subjectType)
            
            % Take average of last 15 measuremnts before the pulse ends
            t = 0:1/this.sfrq:(length(x)*(1/this.sfrq))-1/this.sfrq;

                % Estimate static and sliding friction value
                
                % New crop data to use only up till peak
                dexPulseStart = min(find(Fp~=0))+30;
                dexPulseEnd = max(find(Fp~=0));
                dexEndState4 = max(find(ts==4))-50;
                
                cf = 20; % cutoff freqnency
                [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
                x = filtfilt(b,a,x); % apply fitler
                                
                x_dot = this.sfrq*diff(x);
                x_dot(end+1) = x_dot(end);
                x_ddot = this.sfrq*diff(x_dot);
                        
                [pks,locs,w,p] = findpeaks(-x_ddot(dexPulseStart:dexEndState4));
                locs = locs + dexPulseStart;
                
                [B,I] = sort(p,'descend');
                locFound = sort(locs(I(1:2)));
                
                dexRange = locFound(1):locFound(2);
                
%                 figure('Position',[210 305 560 420]);
%                 ax1 = subplot(4,1,1); plot(x); hold on;
%                 plot(locFound,x(locFound),'o');hold on;
% 
%                 ax2 = subplot(4,1,2); plot(x_dot); hold on;
%                 plot(locFound,x_dot(locFound),'o');hold on;
% 
%                 ax3 = subplot(4,1,3); plot(x_ddot);hold on;
%                 plot(locFound,x_ddot(locFound),'o');hold on; 
% 
%                 ax4 = subplot(4,1,4); plot(Fp);
%                 linkaxes([ax1,ax2,ax3,ax4],'x');xlim([dexPulseStart, dexEndState4]);
                
                if(isempty(dexRange))
                    error('Cutting failed: Look for peak finding error');
                end

%                 figure; plot(x_dot(dexPulseStart:dexPulseEnd+extraDex)); hold on;
%                 plot(x_dot(dexRange),'o');

                x_old = x;
                x_dot_old = x_dot;
                x_ddot_old = x_ddot;
                Fp_old = Fp;
                x = x(dexRange);
                x_dot = x_dot(dexRange);
                x_ddot = x_ddot(dexRange);
                Fp = Fp(dexRange);
                
                % Robot info
                k_r = this.k_r;
                x_r = this.x_r; %(Important too check)
                
                % Decision Variables
                f_s = 1; % Static friction
                f_d = 1; % Dynamics friction
                k_h = 300;
                x_h = f(1)/k_h + x(1);
                B = 1;
                X0 = [k_h,x_h,B,x_dot];
                
                % 300/(2.75*2*pi)^2
                
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb = [0,-2,0,-10*ones(size(x_dot))];
                ub = [5000,2,500,10*ones(size(x_dot))];
                options = optimoptions('fmincon','Display','off');
                                
%                 options = optimoptions('fmincon','TolFun', 1e-8, 'MaxIter', 10000, ...
%                                        'MaxFunEvals', 100000, 'Display', 'off' , ...
%                                        'DiffMinChange', 0.001, 'Algorithm', 'sqp');
                [X_hat,FVAL,EXITFLAG,OUTPUT] = fmincon(@(X_hat)this.costFunc_con(x_dot,X_hat),X0,A,b,Aeq,beq,lb,ub,@(X_hat)this.nl_con(x,x_dot,x_ddot,x_r,Fp,X_hat),options);
                
                t_plot = t(dexRange)-t(dexRange(1));
                t_old = t-t(dexRange(1));

                f_threshold = (abs(x_dot)> 3e-3) | (abs(x_ddot) > 0.4);
                f_threshold_dex = find(diff(f_threshold)>0);
                                    
                if(EXITFLAG >= 1) % Local minimum successfully found
                    %                 m = X_hat(1);
                    %                 f_s = X_hat(2);
                    %                 f_d = X_hat(1);
                    k_hat = X_hat(1);
                    x_h = X_hat(2)-this.x_handelStart;
                    B = X_hat(3);
                    
                    [cost,x_hat,x_dot_hat] = this.costFunc(x,x_dot,x_ddot,x_r,Fp,X_hat);
                                        
                    [VAF] = get_VAF(this,x_dot,x_dot_hat);
        
                    colorVec = {[0.9290, 0.6940, 0.1250],...
                        [0, 0.4470, 0.7410],...
                        [0.8500, 0.3250, 0.0980]};
                    
%                     figure('Position',[771 305 560 420]);
%                     
%                     ax1 = subplot(4,1,1);
%                     plot(t_old,x_old,'color',colorVec{1},'linewidth',2);hold on;
%                     plot(t_plot,x,'linewidth',3,'color',colorVec{2}); hold on;
%                     plot(t_plot,x_hat,'.','markersize',10,'color',colorVec{3});
%                     title(['k_{hat} = ',num2str(k_hat),...
%                         ', x_{hat} = ',num2str(x_h),...
%                         ', B = ',num2str(B)]);
%                     xlabel('Time (s)'); ylabel('Postion (m)');
%                     set(gca,'fontsize',16);
%                     
%                     ax2 = subplot(4,1,2);
%                     plot(t_old,x_dot_old,'color',colorVec{1},'linewidth',2); hold on;
%                     plot(t_plot,x_dot,'linewidth',3,'color',colorVec{2}); hold on;
%                     plot(t_plot,x_dot_hat,'.','markersize',10,'color',colorVec{3});
%                     xlabel('Time (s)'); ylabel('Velocity (m/s)');
%                     title(f_target);
%                     set(gca,'fontsize',16);
% 
%                     ax3 = subplot(4,1,3); plot(t_old,x_ddot_old,'color',colorVec{1},'linewidth',2);hold on;
%                     plot(t_old(locFound),x_ddot_old(locFound),'ok');hold on;
%                     xlabel('Time (s)'); ylabel('Acc (m/s^2)'); set(gca,'fontsize',16);
% 
%                     
%                     ax4 = subplot(4,1,4); plot(t_old,Fp_old,'color',colorVec{1},'linewidth',2);
%                     xlabel('Time (s)'); ylabel('Force Pret (N)'); set(gca,'fontsize',16);
%                     linkaxes([ax1,ax2,ax3,ax4],'x'); xlim([t_old(dexPulseStart), t_old(dexEndState4-300)]);


                else % Did not converge
                     k_hat = NaN;
                     x_h = NaN;
                     B = NaN;
                     VAF = NaN;
                     x_dot_hat = NaN;
                     
%                     figure;
%                     
%                     subplot(2,1,1);
%                     plot(t_plot,x,'linewidth',3); hold on;
%                     title('Failed to converge');
%                     xlabel('Time (s)'); ylabel('Postion (m)');
%                     set(gca,'fontsize',16);
%                     
%                     subplot(2,1,2);
%                     plot(t_plot,x_dot,'linewidth',3); hold on;
%                     plot(t_plot,f_threshold*max(x_dot),':k','linewidth',3);
%                     xlabel('Time (s)'); ylabel('Velocity (m/s)');
%                     title(f_target);
%                     set(gca,'fontsize',16);
                    
                end
                
                
                
                    figure(1);
                    subplot(2,1,1);
                    plot(t_plot,x,'-b'); hold on;
                    if(EXITFLAG == 1)
                        plot(t_plot,x_hat,':r');
                    end
                    xlabel('Time (s)'); ylabel('Postion (m)');
                    set(gca,'fontsize',16);
                    
                    subplot(2,1,2);
                    plot(t_plot,x_dot,'-b'); hold on;
                    if(EXITFLAG == 1)
                        plot(t_plot,x_dot_hat,':r');
                        plot(t_plot(f_threshold_dex),0,'.g','markersize',10);
                    else
                        plot(t_plot(f_threshold_dex),0,'.k','markersize',10);
                    end
                    xlabel('Time (s)'); ylabel('Velocity (m/s)');
                    set(gca,'fontsize',16); grid on;
                

                est_pulse.k_hat = k_hat;
                est_pulse.b_hat = B;
                est_pulse.x_h = x_h;
                est_pulse.x_dot = x_dot;
                est_pulse.t = t_plot;
                est_pulse.dexRange = dexRange;
                est_pulse.x_dot_hat = x_dot_hat;
                est_pulse.VAF = VAF;
                
%                 disp('test');
                
            
        end
        
        function [cost,x_hat,x_dot_hat] = costFunc(this,x,x_dot,x_ddot,x_r,Fp,X_hat)
            
            
            % Start at the measured postion
            x_dot_hat = zeros(size(x));
            x_hat = zeros(size(x));
            x_dot_hat(1) = x_dot(1);
            x_hat(1) = x(1);
            m = this.m;
            
            for i = 1:length(x)-1
                
                % Get forces at a spesific step
                f_noFric = this.f_dynamics(x(i),Fp(i),x_r,X_hat);
                f_f = this.f_friction(x(i),Fp(i),x_r,X_hat,x_dot(i),x_ddot(i));
                x_ddot_hat = (1/m)*(f_noFric - f_f);
                
                % Dyanmic contraint
                x_dot_hat(i+1) = x_dot_hat(i) + (1/this.sfrq)*x_ddot_hat;
                x_hat(i+1) = x_hat(i) + (1/this.sfrq)*x_dot_hat(i);
                
            end
            
            % Evaluate cost function
            cost = sum(abs(x_dot_hat-x_dot).^2);

%             cost = sum(abs(x_hat-x).^2 + abs(x_dot_hat-x_dot).^2) + (1/this.sfrq)*abs(x_ddot_hat-x_ddot(end)).^2;
           
        end
        
        function [cost] = costFunc_con(this,x_dot,X_hat)
            
            x_dot_hat = X_hat(4:end);

            % Evaluate cost function
            cost = sum(abs(x_dot_hat-x_dot).^2);
           
        end
        
        function [c,ceq] = nl_con(this,x,x_dot,x_ddot,x_r,Fp,X_hat)
            
            % Start at the measured postion
            x_dot_hat = zeros(size(x));
            x_hat = zeros(size(x));
            x_dot_hat(1) = x_dot(1);
            x_hat(1) = x(1);
            m = this.m;
                        
            for i = 1:length(x)-1
                
                % Get forces at a spesific step
                f_noFric = this.f_dynamics(x(i),Fp(i),x_r,X_hat(1:3));
                f_f = this.f_friction(x(i),Fp(i),x_r,X_hat(1:3),x_dot(i),x_ddot(i));
                x_ddot_hat = (1/m)*(f_noFric - f_f);
                
                % Dyanmic contraint
                x_dot_hat(i+1) = x_dot_hat(i) + (1/this.sfrq)*x_ddot_hat;
                x_hat(i+1) = x_hat(i) + (1/this.sfrq)*x_dot_hat(i);
                
                ceq(i) = (x_dot_hat(i) + (1/this.sfrq)*x_ddot_hat) - x_dot_hat(i+1);
                
            end
            
            ceq(end+1) = x_dot(1)-x_dot_hat(1);
            
            c = [];
            
        end
        
        function [f_noFric] = f_dynamics(this,x,Fp,x_r,X_hat)
           
            k_r = this.k_r;
            k_h = X_hat(1);
            x_h = X_hat(2);
            
            f_noFric = k_h*(x_h-x) + k_r*(x_r-x) + Fp;
            
        end
        
        function [f_f] = f_friction(this,x,Fp,x_r,X_hat,x_dot,x_ddot)
            
            m = this.m;
            k_h = X_hat(1);
            x_h = X_hat(2);
            B = X_hat(3);
            
            if( (abs(x_dot)> 3e-3) || (abs(x_ddot) > 0.4) )
                f_f = (B/m)*(x_dot);%f_d;
            else
                f_f = this.f_dynamics(x,Fp,x_r,X_hat); %f_s;
            end
            
        end
        
        function [c,ceq] = mycon(this,x)
            c = 0;     % Compute nonlinear inequalities at x.
            ceq = 0;   % Compute nonlinear equalities at x.
        end 
               
        function [est_stocastic] = get_singleTrial_k_hat_stocastic(this,X1,Y1,ts)
        
%             figure; plot(X1);
%             figure; plot(Y1);
%             figure; plot(ts);
            
            % Cutting before changes during visit 11/3/21
%             dexStart = min(find(ts==3))+500;
%             dexEnd = min(find(ts==4));
            
            dexStart = max([min(find(X1~=0)),min(find(ts==4))])+500; % Preturbation does not start at state 4 always
%             dexStart = min(find(ts==4))+500;
            dexEnd = max(find(ts==4));
            X1 = X1(dexStart:dexEnd);
            Y1 = Y1(dexStart:dexEnd);
            
            %% Variables
            nfft=2000;
            nFreq = nfft/2+1;
            Hz = this.sfrq;
            
            unwrapThreshold = -50*pi/180;
            
            % Admittance
            Y11 = zeros(nFreq,1);
            
            % Impedance 
            Z11 = zeros(nFreq,1);     
            Z11_mag = zeros(nFreq,1);     
            Z11_phi = zeros(nFreq,1);
            PC11 = zeros(nFreq,1);
            
            
            % Check xcorr preturbations
%             [c,lags] = xcorr(X1,X1,220000,'coeff'); title('xcorr(X1,X1)');
%             figure; subplot(3,1,1); stem(lags,c); hold on;
%             [c,lags] = xcorr(X2,X2,220000,'coeff'); title('xcorr(X2,X2)');
%             subplot(3,1,2); stem(lags,c);
%             [c,lags] = xcorr(X1,X2,220000,'coeff'); title('xcorr(X1,X2)');
%             subplot(3,1,3); stem(lags,c);
            
            Y1 = detrend(Y1);
            
            TF_freq = Hz/2*linspace(0,1,nfft/2+1);
            
%             figure;
%             subplot(2,1,1); plot(X1);
%             subplot(2,1,2); plot(Y1);
            
            % Cross power spectral density
            [Px1x1,Fx1x1] = cpsd(X1,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Px1y1,Fx1y1] = cpsd(Y1,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            
            [Py1x1,Fy1x1] = cpsd(X1,Y1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Py1y1,Fy1y1] = cpsd(Y1,Y1,hamming(nfft/2),nfft/4,nfft,Hz);
            
            %% Ordinary Coherence
            [OCx1y1,OCx1y1_F] = mscohere(X1,Y1,hamming(nfft/2),nfft/4,nfft,Hz);
            
%             figure; semilogx(OCx1y1_F,OCx1y1); title('C: x_1 y_1'); grid on;
            
            for j=1:nFreq
                Y11(j,1) = Px1y1(j,1)/Px1x1(j,1);
            end
            
            % Invert compliance to get impedance
            Z11 = 1./Y11;
            
            %% Mag & Phase calculation for Impedance plot
            for j=1:nFreq
                Z11_mag(j,1) = abs(Z11(j,1));
                Z11_phi(j,1) = 180/pi*unwrap2(angle(Z11(j,1)),unwrapThreshold,'up');
            end
            
%             [tmp,dexK] = min(abs(TF_freq - 1));
%             [tmp,dexK] = find(0.5 < TF_freq  & TF_freq < 3);
            [tmp,dexK] = find(1 < TF_freq  & TF_freq < 1.5);

            k_hat = mean(Z11_mag(dexK)) - this.k_r;
            OC_hat = OCx1y1(dexK,1);
            
            %% Impedance plot (Diagonal)
            xLowerLim = 0.5;
            xUpperLim = 50.0;
            yLowerLim11 = 1e+2;
            yUpperLim11 = 1e+6;
            yLowerLim22 = 2e+3;
            yUpperLim22 = 3e+4;
            
%             zeta = 1;
%             b_fit = 2*zeta*this.m*sqrt(k_hat/this.m);
%             s = tf('s');
% %             G = (this.m*s^2+b_fit*s+k_hat)*(this.m*s^2+b_fit*s+k_hat)/200
% %             G = (1/10)*(s^3);
%             G = tf([this.m,b_fit,k_hat],1);
% %             opts = bodeoptions('cstprefs');
% %             opts.MagUnits = 'abs';
% %             opts.MagScale = 'log';
% %             opts.FreqUnits = 'rad/s';
% %             figure;bode(G,TF_freq,opts);
%             [MAG,PHASE] = bode(G,TF_freq*2*pi);
%             MAG = squeeze(MAG);
%             PHASE = squeeze(PHASE);
%             
%             figure; loglog(TF_freq,MAG);
%             
%             
%             figure(1); hold on;%figure('position',[440 61 560 736]); 
%             % Magnitude plot of ankle impedance
%             ax1 = subplot(3,1,1,'XScale','log','YScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_mag(:,1),'LineWidth',2); grid on; box on;
% %             plot(TF_freq,MAG); 
%             axis([xLowerLim xUpperLim yLowerLim11 yUpperLim11]);
%             plot(TF_freq(dexK),Z11_mag(dexK,1),'.','markersize',20);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
%             ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
%             title('Z11','fontWeight','bold','fontSize',16);
%             
%             ax2 = subplot(3,1,2,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_phi(:,1),'LineWidth',2); grid on; box on;
% %             plot(TF_freq,PHASE);
%             axis([xLowerLim xUpperLim 0 180]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
%             
%             % Partial Coherence Plot
%             xLowerLim = 0.5;
% %             figure; hold on;
%             set(gcf,'Color',[1,1,1]);
%             
%             ax3 = subplot(3,1,3,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(OCx1y1_F,OCx1y1,'LineWidth',2);
%             plot(OCx1y1_F(dexK),OCx1y1(dexK,1),'.','markersize',20);
% 
%             grid on;box on; ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
%             xlabel('Frequency (Hz)','fontWeight','bold','fontSize',14);
%             title('Y11 PC','fontWeight','bold','fontSize',16);
%             
%             linkaxes([ax1,ax2,ax3],'x');

            est_stocastic.k_hat = k_hat;
            est_stocastic.OC_hat = OC_hat;
            
        end
        
        function [est_release] = get_singleTrial_k_hat_release(this,x,ts,subjectType,distVal)
            
            % Take average of last 15 measuremnts before the pulse ends
            t = 0:1/this.sfrq:(length(x)*(1/this.sfrq))-1/this.sfrq;
                
            dexStart = min(find(ts== 5)) + 20;
            dexEnd = dexStart + 150;% max(find(ts== 6));
            
            cf = 20; % cutoff freqnency
            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            x = filtfilt(b,a,x); % apply fitler
            
            x_dot = this.sfrq*diff(x);
            x_dot(end+1) = x_dot(end);
            x_ddot = this.sfrq*diff(x_dot);
            x_ddot(end+1) = x_ddot(end);
            
%             [pks,locs,w,p] = findpeaks(-x_ddot(dexStart:dexEnd));
%             locs = locs + dexStart;
%             
%             [B,I] = sort(p,'descend');
%             locFound = sort(locs(I(1:2)));
    
%             dexRange = locFound(1):locFound(2);
            

                % hard code duration
              dexRange = dexStart:dexEnd;
                
%                 figure('Position',[210 305 560 420]);
%                 ax1 = subplot(3,1,1); plot(x,'linewidth',2.5); hold on;
%                 plot(dexStart,x(dexStart),'.r','markersize',25);hold on;
%                 plot(dexEnd,x(dexEnd),'.r','markersize',25);hold on;
% 
% 
%                 ax2 = subplot(3,1,2); plot(x_dot,'linewidth',2.5); hold on;
%                 plot(dexStart,x_dot(dexStart),'.r','markersize',25);hold on;
%                 plot(dexEnd,x_dot(dexEnd),'.r','markersize',25);hold on;
%                 
%                 ax3 = subplot(3,1,3); plot(x_ddot,'linewidth',2.5);hold on;
%                 plot(dexStart,x_ddot(dexStart),'.r','markersize',25);hold on;
%                 plot(dexEnd,x_ddot(dexEnd),'.r','markersize',25);hold on;
%                 
%                 linkaxes([ax1,ax2,ax3],'x');xlim([dexStart-100, dexEnd]);
                
                if(isempty(dexRange))
                    error('Cutting failed: Look for peak finding error');
                end

%                 figure; plot(x_dot(dexStart:dexEnd+extraDex)); hold on;
%                 plot(x_dot(dexRange),'o');

                x_old = x;
                x_dot_old = x_dot;
                x_ddot_old = x_ddot;
                x = x(dexRange);
                x_dot = x_dot(dexRange);
                x_ddot = x_ddot(dexRange);

                % Decision Variables
                k_h = 300;
                x_h = distVal;
                B = 1;
                X0 = [k_h,x_h,B,x_dot]; %zeros(size(x_dot))];
                
                % 300/(2.75*2*pi)^2
                
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb = [0,-2,0,-10*ones(size(x_dot))];
                ub = [5000,2,500,10*ones(size(x_dot))];
                options = optimoptions('fmincon','Display','off');
                                
%                 options = optimoptions('fmincon','TolFun', 1e-8, 'MaxIter', 10000, ...
%                                        'MaxFunEvals', 100000, 'Display', 'off' , ...
%                                        'DiffMinChange', 0.001, 'Algorithm', 'sqp');

                [X_hat,FVAL,EXITFLAG,OUTPUT] = fmincon(@(X_hat)this.costFunc_con(x_dot,X_hat),X0,A,b,Aeq,beq,lb,ub,@(X_hat)this.nl_con_release(x,x_dot,X_hat),options);
                

                t_plot = t(dexRange)-t(dexRange(1));
                t_old = t-t(dexRange(1));
                                    
                if(EXITFLAG >= 1) % Local minimum successfully found
                    %                 m = X_hat(1);
                    %                 f_s = X_hat(2);
                    %                 f_d = X_hat(1);
                    k_hat = X_hat(1);
                    x_h = X_hat(2)-this.x_handelStart;
                    B = X_hat(3);
                    
                    [cost,x_hat,x_dot_hat] = this.costFunc_release(x,x_dot,x_ddot,X_hat);
                                        
                    [VAF] = get_VAF(this,x_dot,x_dot_hat);
        
                    colorVec = {[0.9290, 0.6940, 0.1250],...
                        [0, 0.4470, 0.7410],...
                        [0.8500, 0.3250, 0.0980]};
                    
                    figure('Position',[771 305 560 420]);
                    
                    ax1 = subplot(2,1,1);
                    plot(t_old,x_old,'color',colorVec{1},'linewidth',2);hold on;
                    plot(t_plot,x,'linewidth',3,'color',colorVec{2}); hold on;
                    plot(t_plot,x_hat,'.','markersize',10,'color',colorVec{3});
                    title(['k_{hat} = ',num2str(k_hat),...
                        ', x_{hat} = ',num2str(x_h),...
                        ', B = ',num2str(B)]);
                    xlabel('Time (s)'); ylabel('Postion (m)');
                    set(gca,'fontsize',16);
                    
                    ax2 = subplot(2,1,2);
                    plot(t_old,x_dot_old,'color',colorVec{1},'linewidth',2); hold on;
                    plot(t_plot,x_dot,'linewidth',3,'color',colorVec{2}); hold on;
                    plot(t_plot,x_dot_hat,'.','markersize',10,'color',colorVec{3});
                    xlabel('Time (s)'); ylabel('Velocity (m/s)');
                    set(gca,'fontsize',16);

                    linkaxes([ax1,ax2],'x'); xlim([t_old(dexStart-100), t_old(dexEnd+100)]);


                else % Did not converge
                     k_hat = NaN;
                     x_h = NaN;
                     B = NaN;
                     VAF = NaN;
                     x_dot_hat = NaN;
                     
%                     figure;
%                     
%                     subplot(2,1,1);
%                     plot(t_plot,x,'linewidth',3); hold on;
%                     title('Failed to converge');
%                     xlabel('Time (s)'); ylabel('Postion (m)');
%                     set(gca,'fontsize',16);
%                     
%                     subplot(2,1,2);
%                     plot(t_plot,x_dot,'linewidth',3); hold on;
%                     plot(t_plot,f_threshold*max(x_dot),':k','linewidth',3);
%                     xlabel('Time (s)'); ylabel('Velocity (m/s)');
%                     title(f_target);
%                     set(gca,'fontsize',16);
                    
                end

                est_release.k_hat = k_hat;
                est_release.b_hat = B;
                est_release.x_h = x_h;
                est_release.x_dot = x_dot;
                est_release.t = t_plot;
                est_release.dexRange = dexRange;
                est_release.x_dot_hat = x_dot_hat;
                est_release.VAF = VAF;
                est_release.x_dot_tot = x_dot_old;
                
%                 disp('test');
                
            
        end
        
        function [cost,x_hat,x_dot_hat] = costFunc_release(this,x,x_dot,x_ddot,X_hat)
            
            
            % Start at the measured postion
            x_dot_hat = zeros(size(x));
            x_hat = zeros(size(x));
            x_dot_hat(1) = x_dot(1);
            x_hat(1) = x(1);
            k_h = X_hat(1);
            x_h = X_hat(2);
            B = X_hat(3);
            
            for i = 1:length(x)-1
                
                % Get forces at a spesific step
                x_ddot_hat = (1/this.m)*(k_h*(x_h-x(i)) - B*(x_dot(i)));
                
                % Dyanmic contraint
                x_dot_hat(i+1) = x_dot_hat(i) + (1/this.sfrq)*x_ddot_hat;
                x_hat(i+1) = x_hat(i) + (1/this.sfrq)*x_dot_hat(i);
                
            end
            
            % Evaluate cost function
            cost = sum(abs(x_dot_hat-x_dot).^2);
           
        end
        
        function [c,ceq] = nl_con_release(this,x,x_dot,X_hat)
            
            % Start at the measured postion
            x_dot_hat = zeros(size(x));
            x_hat = zeros(size(x));
            x_dot_hat(1) = x_dot(1);
            x_hat(1) = x(1);
            k_h = X_hat(1);
            x_h = X_hat(2);
            B = X_hat(3);
            X = X_hat(4:length(x));
            
            X_dot = X_hat(end-length(x)+1:end);

                     
            
            for i = 1:length(x)-1
                
                % Get forces at a spesific step
                  % Get forces at a spesific step
                x_ddot_hat = (1/this.m)*(k_h*(x_h-x(i)) - B*(x_dot(i)));
                
                % Dyanmic contraint
                x_dot_hat(i+1) = x_dot_hat(i) + (1/this.sfrq)*x_ddot_hat;
                x_hat(i+1) = x_hat(i) + (1/this.sfrq)*x_dot_hat(i);
                
                ceq1(i) = X - x_hat(i+1);
                ceq2(i) = X_dot(i+1) - x_dot_hat(i+1);
                
            end
            
            ceq3 = x_dot(1)-x_dot_hat(1);
            
            ceq = [ceq1,ceq2,ceq3];
            
            c = 0;
            
        end

        
        function [h_model, B_opt, K_opt, VAFirf] = fitModel(this, h_hat, L,nSampDelay)
            
            %% Part 2-2: 2nd order model approximation
            % Optimization to find the best I,B,K approximates (cost function may change..)
            
            %         M0 = 0.02;
            %         B0 = 0.2;
            %         K0 = 10;
            %
            %         LB = [0,    0.1,   0];
            %         UB = [0.15, 5.0, 100];
            
            M0 = this.m;
            B0 = 20;
            K0 = 600;

%             M0 = Z0(1);
%             B0 = Z0(2);
%             K0 = Z0(3);
            
%             LB = [-30, -400, -20000];
%             UB = [30, 400, 20000];
            LB = [-400, -20000];
            UB = [400, 20000];
            
            % Bounded Nonlinear Optimization
            options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun',0.01);
            [x fval exitflags outputs] = fminsearchbnd(@(x)this.irfFitting_knownM(x,M0,h_hat,L,nSampDelay),[B0 K0],LB,UB,options);
%             M_opt = x(1);
%             B_opt = x(2);
%             K_opt = x(3);
            B_opt = x(1);
            K_opt = x(2);
            
            Y_model = tf(1,[M0,B_opt,K_opt]);
            h_model = impulse(Y_model,0:1/this.sfrq:L/this.sfrq);
            
            N = length(h_hat);
            
%             figure;
%             plot(1:N,h_hat,'.b','markersize',10); hold on;
%             plot(nSampDelay:N,h_hat(nSampDelay:N),'ob'); hold on;
%             plot(nSampDelay:N,h_model(nSampDelay:N),'.r','markersize',10);
%             title(['k_{hat} = ',num2str(K_opt),...
%                  ', B_{hat} = ',num2str(B_opt)]);
%             legend('Total','Used in fit','h_{model}');
%             
            % VAF IRF Calculation
%             VAFirf = 100*(1-std(h_hat-h_model)^2/std(h_hat)^2);
            
            
        end
        
        function fval = irfFitting(this,x,h_hat,L,nSampDelay)
            Y_model = tf(1,[x(1),x(2),x(3)]);
            h_model = impulse(Y_model,0:1/this.sfrq:L/this.sfrq);
            
            fval=0;
            for j=nSampDelay:length(h_hat)
                        fval = fval + (h_hat(j)-h_model(j))^2;
                
%                         % Weighted cost function 2
%                         if(j<M2/2)
%                             fval = fval + (h_hat(j-M1+1,1)-h_model(j-M1+1,1))^2;
%                         else
%                             weight = 2*(1-j/M2);
%                             fval = fval + (h_hat(j-M1+1,1)-h_model(j-M1+1,1))^2*weight;
%                         end
            end
            fval = sqrt(fval/L);
        end

        function fval = irfFitting_knownM(this,x,M0,h_hat,L,nSampDelay)
            Y_model = tf(1,[M0,x(1),x(2)]);
            h_model = impulse(Y_model,0:1/this.sfrq:L/this.sfrq);
            
            fval=0;
            for j=nSampDelay:length(h_hat)
                        fval = fval + (h_hat(j)-h_model(j))^2;
            end
            fval = sqrt(fval/L);
        end
        
        function [VAFirf] = get_VAF(this,h_model,h_hat)
            
            % VAF IRF Calculation
            N = length(h_model);
            VAFirf = zeros(N,1);
            VAFirf = 100*(1-std(h_hat-h_model)^2/std(h_hat)^2);
            
        end

    end
end


