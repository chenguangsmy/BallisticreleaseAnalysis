classdef crossTrialAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a 
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.  
   
    properties
        
        k_hat_pulse
        k_hat_stocastic
        k_hat_release
        sfrq
        k_r
        m
        
    end
    
    methods
        function [this] = crossTrialAnalysis(data,sfrq,f_target)
            
            this.sfrq = sfrq;
            this.k_r = 300;
            
            % Look at data
%             this.plot_postion();
%             this.plot_force();
            
            % Estimate stiffnesses
            this.get_k_hat_pulse(data,f_target);
            this.get_k_hat_stocastic(data);
            this.get_k_hat_release(data);
            
        end
        
        function [] = get_k_hat_pulse(this,data,f_target)
            
            % Exclude empty cells (FIX Later)
            step_Pulse = 2;
            [N_trial] = size(data,1);
            count = 1;
            for trial = 1:N_trial-1
                if(~isempty(data{trial,step_Pulse}))
                    if( sum(data{trial,step_Pulse}.Fp(2,:)~=0) ~= 0 )
                        if(sum(data{trial,step_Pulse}.Fp(2,1:30)) == 0) % find pulse to close to start
                        dexTrialNonzero(count) = trial;
                                        count = count+1;
                        end
                    end
                end
            end
            
            this.k_hat_pulse = NaN*ones(1,15);
            for i = 1:length(dexTrialNonzero)
                clear tmpData
                tmpData = data{dexTrialNonzero(i),step_Pulse};
                this.k_hat_pulse(i) = get_singleTrial_k_hat_pulse(this,tmpData.f(2,:),...
                                                 tmpData.Fp(2,:),...
                                                 tmpData.x(2,:),...
                                                 tmpData.ts,...
                                                 f_target);
            end
                        
        end
        
        function [] = get_k_hat_stocastic(this,data)
            
            % Exclude empty cells (FIX Later)
            step_Stocastic = 3;
            [N_trial] = size(data,1);
            dexTrialNonzero = [];
            count = 1;
            for trial = 1:N_trial
                if(~isempty(data{trial,step_Stocastic}))
                    dexTrialNonzero(count) = trial;
                end
                count = count+1;
            end
            
            if(~isempty(dexTrialNonzero))
            
                this.k_hat_stocastic = NaN*ones(1,9);
                for i = 1:length(dexTrialNonzero)
                    clear tmpData
                    tmpData = data{dexTrialNonzero(i),step_Stocastic};
                    this.k_hat_stocastic(i) = get_singleTrial_k_hat_stocastic(this,tmpData.Fp(2,:),...
                        tmpData.x(2,:),...
                        tmpData.ts);
                end
                
            end
            
        end
        
        function [] = get_k_hat_release(this,data)
               
            % Exclude empty cells (FIX Later)
            release = 1;
            targetForce = 15;
            dexTrialNonzero = [];
            [N_trial] = size(data,1);
            count = 1;
            for trial = 2:N_trial % Exclude first trial
                if(~isempty(data{trial,release}))
                    dexTrialNonzero(count) = trial;
                end
                count = count+1;
            end
            
            if(~isempty(dexTrialNonzero))
                
                this.k_hat_release = NaN*ones(1,15);
                for i = 1:length(dexTrialNonzero)
                    clear tmpData
                    tmpData = data{dexTrialNonzero(i),release};
                    this.k_hat_release(i) = get_singleTrial_k_hat_release(this,tmpData.x(2,:),tmpData.ts,targetForce);
                end
                
            end
            
        end
        
        function [k_hat] = get_singleTrial_k_hat_pulse(this,f,Fp,x,ts,f_target)
            
            % Take average of last 15 measuremnts before the pulse ends
            t = 0:1/this.sfrq:(length(x)*(1/this.sfrq))-1/this.sfrq;
%             dexPulseStart = min(find(Fp~=0));
%             dexRange1 = [(dexPulseStart-15):dexPulseStart];
%             
%             dexPulseEnd = max(find(Fp~=0));
%             dexRange2 = [(dexPulseEnd-15):dexPulseEnd];
                        
%                 f1 = mean(f(dexRange1));
%                 f2 = mean(f(dexRange2));
%                 
%                 x1 = mean(x(dexRange1));
%                 x2 = mean(x(dexRange2));
%                 
%                 k_hat = -(f2-f1)/(x2-x1);
                
                % Look at raw data
%                             figure;
%                             ax1 = subplot(4,1,1); plot(t,f); hold on;
%                 %             plot(dexPulseStart,f1,'+',dexPulseEnd,f2,'*');
%                 
%                             ax2 = subplot(4,1,2); plot(t,x); hold on;
%                 %             plot(dexPulseStart,x1,'+',dexPulseEnd,x2,'*');
%                 
%                             ax3 = subplot(4,1,3); plot(t,Fp);
%                             ax4 = subplot(4,1,4); plot(t,ts);
%                             linkaxes([ax1,ax2,ax3,ax4],'x');

%                 figure; 
%                 ax1 = subplot(2,1,1); plot(t(1800:2050),Fp(1800:2050),'linewidth',2.5);
%                 xlabel('Time (s)'); ylabel('F_p (N)');
%                 xlim([t(1800) t(2050)]); set(gca,'fontsize',18);
%                 ax2 = subplot(2,1,2); plot(t(1800:2050),1000*(x(1800:2050)-x(1800)),'linewidth',2.5); hold on;
%                 xlabel('Time (s)'); ylabel('x (mm)');xlim([t(1800) t(2050)]);set(gca,'fontsize',18);

            % Double integrate force
%             f_tmp = f(dexPulseStart:dexPulseEnd) - mean(f(dexPulseStart-500:dexPulseStart));
%             f_int(1) = 0;
%             for i = 1:length(f_tmp)
%                 f_int(i+1) = f_int(i) + (1/this.sfrq)*(f_tmp(i));
%             end
%             
%             f_int_int(1) = 0;
%             for i = 1:length(f_tmp)
%                 f_int_int(i+1) = f_int_int(i) + (1/this.sfrq)*(f_int(i));
%             end
            
%             figure; 
%             subplot(3,1,1); plot(f);
%             subplot(3,1,2); plot(f_int);
%             subplot(3,1,3); plot(f_int_int);
                
                %% Try new idea: estimate static and sliding friction value
                
                % New crop data to use only up till peak
                dexPulseStart = min(find(Fp~=0))+30;
                dexPulseEnd = max(find(Fp~=0));
                extraDex = 100;
                
                cf = 20; % cutoff freqnency
                [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
                x = filtfilt(b,a,x); % apply fitler
                                
                x_dot = this.sfrq*diff(x);
                x_dot(end+1) = x_dot(end);
                
                % (IMPORTANT) Cut with peak for springs
                [pks,locs,w,p] = findpeaks((-sign(mean(Fp))) * x_dot(dexPulseStart:dexPulseEnd+extraDex));
                locs = locs + dexPulseStart;
                dexRange = dexPulseStart:locs(find(max(w)==w));
                
                % (IMPORTANT) Use entire pulse for humans
%                 dexRange = dexPulseStart:dexPulseEnd;
                
                if(isempty(dexRange))
                    error('Cutting failed: Look for peak finding error');
                end

%                 figure; plot(x_dot(dexPulseStart:dexPulseEnd+extraDex)); hold on;
%                 plot(x_dot(dexRange),'o');
%                 
%                 dexPulseStart = dexPulseStart + 25;
%                 dexPulseEnd = dexPulseEnd - 10;
                x = x(dexRange);
                x_dot = x_dot(dexRange);
                x_ddot = this.sfrq*diff(x_dot);
                x_ddot(end+1) = x_ddot(end);
                Fp = Fp(dexRange);
                
                % Robot info
                k_r = this.k_r;
                x_r = 0.481-f_target/k_r; %(Important too check)
                
                % Decision Variables
                f_s = 1; % Static friction
                f_d = 1; % Dynamics friction
                k_h = 300;
                x_h = f(1)/k_h + x(1);
                B = 1;
                X0 = [k_h,x_h,B];
                
                % 300/(2.75*2*pi)^2

                this.m = 1.25;
                
                A = [];
                b = [];
                Aeq = [];
                beq = [];
                lb = [0,-2,0];
                ub = [5000,2,100];
%                 options = optimoptions(fmincon,'TolFun', 0.00000001, 'MaxIter', 10000, ...
%                        'MaxFunEvals', 100000, 'Display', 'notify' , ...
%                        'DiffMinChange', 0.001, 'Algorithm', 'sqp');
                [X_hat,FVAL,EXITFLAG,OUTPUT] = fmincon(@(X_hat)this.costFunc(x,x_dot,x_ddot,x_r,Fp,X_hat),X0,A,b,Aeq,beq,lb,ub);
                
                if(EXITFLAG == 1) % Local minimum successfully found
                    %                 m = X_hat(1);
                    %                 f_s = X_hat(2);
                    %                 f_d = X_hat(1);
                    k_hat = X_hat(1);
                    x_h = X_hat(2);
                    B = X_hat(3);
                    
                    [cost,x_dot_hat] = this.costFunc(x,x_dot,x_ddot,x_r,Fp,X_hat);
                    
                    t_plot = t(dexRange)-t(dexRange(1));
                    
                    figure;
                    plot(t_plot,x_dot,'linewidth',3); hold on;
                    plot(t_plot,x_dot_hat,'.','markersize',10);
                    title(['k_{hat} = ',num2str(k_hat),...
                        ', x_{hat} = ',num2str(x_h),...
                        ', B = ',num2str(B)]);
                    xlabel('Time (s)'); ylabel('Velocity (m/s)');
                    set(gca,'fontsize',16);


                else % Did not converge
                     k_hat = NaN;
                     x_h = NaN;
                     B = NaN;
                end
                

%                    disp('test');
                
            
        end
        
        function [cost,x_dot_hat] = costFunc(this,x,x_dot,x_ddot,x_r,Fp,X_hat)
            
            
            % Start at the measured postion
            x_dot_hat = zeros(size(x));
            x_dot_hat(1) = x_dot(1);
            m = this.m;
            
            for i = 1:length(x)-1
                
                % Get forces at a spesific step
                f_noFric = this.f_dynamics(x(i),Fp(i),x_r,X_hat);
                f_f = this.f_friction(x(i),Fp(i),x_r,X_hat,x_dot(i),x_ddot(i));
                x_ddot_hat = (1/m)*(f_noFric - f_f);
                
                % Dyanmic contraint
                x_dot_hat(i+1) = x_dot_hat(i) + (1/this.sfrq)*x_ddot_hat;
                
            end
            
            % Evaluate cost function
            cost = sum(abs(x_dot_hat-x_dot).^2);% + abs(x_ddot_hat-x_ddot(end)).^2);
            
        end
        
        function [f_noFric] = f_dynamics(this,x,Fp,x_r,X_hat)
           
            k_r = this.k_r;
            k_h = X_hat(1);
            x_h = X_hat(2);
            
            f_noFric = ( k_h*(x_h-x) + k_r*(x_r-x) + Fp );
            
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
            
        function [k_hat] = get_singleTrial_k_hat_stocastic(this,X1,Y1,ts)
        
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
            
            %% Impedance plot (Diagonal)
            xLowerLim = 0.5;
            xUpperLim = 50.0;
            yLowerLim11 = 1e+2;
            yUpperLim11 = 1e+6;
            yLowerLim22 = 2e+3;
            yUpperLim22 = 3e+4;
            
%             figure(1); hold on;%figure('position',[440 61 560 736]); 
%             % Magnitude plot of ankle impedance
%             ax1 = subplot(3,1,1,'XScale','log','YScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_mag(:,1),'LineWidth',2); grid on; box on;
%             axis([xLowerLim xUpperLim yLowerLim11 yUpperLim11]);
%             plot(TF_freq(dexK),Z11_mag(dexK,1),'.','markersize',20);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
%             ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
%             title('Z11','fontWeight','bold','fontSize',16);
%             
%             ax2 = subplot(3,1,2,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_phi(:,1),'LineWidth',2); grid on; box on;
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
%             xlabel('Hz','fontWeight','bold','fontSize',14);
%             title('Y11 PC','fontWeight','bold','fontSize',16);
%             
%             linkaxes([ax1,ax2,ax3],'x');
            
        end
        
        function [k_hat] = get_singleTrial_k_hat_release(this,x,ts,targetForce)
           
            % Cutting used prior to visit 11/3/21
%             dexStart = min(find(ts== 4));
%             dexEnd = max(find(ts == 5)) - 1; % Subtract 1 for diff
            
            dexStart = min(find(ts== 5)); % (FIX) add artifical lag
            dexEnd = dexStart+300;%length(x)-1; 

            % Filter
            cf = 20; % cutoff freqnency
            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            x = filtfilt(b,a,x); % apply fitler
            
            % velocity estimate multiplied by inverse target force is an
            % estimate of the impulse response function
            h_hat = (1/targetForce)*this.sfrq*diff(x); 
%             figure; plot(h_hat);

%             figure;
%             ax1 = subplot(2,1,1); plot(h_hat(dexStart:dexEnd));
%             ax2 = subplot(2,1,2); plot(ts(dexStart:dexEnd));
%             linkaxes([ax1,ax2],'x');
            
            [h_model, M, B, K] =  this.fitModel(h_hat(dexStart:dexEnd), dexEnd-dexStart);
            
            k_hat = K;
        end
        
        function [h_model, M_opt, B_opt, K_opt, VAFirf] =  fitModel(this, h_hat, L)
            
            %% Part 2-2: 2nd order model approximation
            % Optimization to find the best I,B,K approximates (cost function may change..)
            
            %         M0 = 0.02;
            %         B0 = 0.2;
            %         K0 = 10;
            %
            %         LB = [0,    0.1,   0];
            %         UB = [0.15, 5.0, 100];
            
            M0 = 2;
            B0 = 40;
            K0 = 2500;

%             M0 = Z0(1);
%             B0 = Z0(2);
%             K0 = Z0(3);
            
            LB = [-30, -400, -20000];
            UB = [30, 400, 20000];
            
            % Bounded Nonlinear Optimization
            options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun',0.01);
            [x fval exitflags outputs] = fminsearchbnd(@(x)this.irfFitting(x,h_hat,L),[M0 B0 K0],LB,UB,options);
            M_opt = x(1);
            B_opt = x(2);
            K_opt = x(3);
            
            Y_model = tf(1,[M_opt,B_opt,K_opt]);
            h_model = impulse(Y_model,0:1/this.sfrq:L/this.sfrq);
            
%             figure;
%             plot(h_hat,'o'); hold on;
%             plot(h_model,'.');
%             
            % VAF IRF Calculation
%             VAFirf = 100*(1-std(h_hat-h_model)^2/std(h_hat)^2);
            
            
        end
        
        function fval = irfFitting(this,x,h_hat,L)
            Y_model = tf(1,[x(1),x(2),x(3)]);
            h_model = impulse(Y_model,0:1/this.sfrq:L/this.sfrq);
            
            fval=0;
            for j=1:length(h_hat)
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

    end
end


