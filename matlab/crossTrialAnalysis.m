classdef crossTrialAnalysis < handle
    %withInSubjectAnalysis This class takes in all the trials for a 
    %specific subject and a condition. It outputs three different stiffness
    %estimates (step, stochastic, and release) in one data structure.  
   
    properties
        
        k_hat_pulse
        k_hat_stocastic
        k_hat_release
        sfrq
        
    end
    
    methods
        function [this] = crossTrialAnalysis(data,sfrq)
            
            this.sfrq = sfrq;
            
            % Look at data
%             this.plot_postion();
%             this.plot_force();
            
            % Estimate stiffnesses
            this.get_k_hat_pulse(data);
            this.get_k_hat_stocastic(data);
            this.get_k_hat_release(data);
            
        end
        
        function [] = get_k_hat_pulse(this,data)
            
            % Exclude empty cells (FIX Later)
            step_Pulse = 2;
            [N_trial] = size(data,1);
            count = 1;
            for trial = 1:N_trial
                if(~isempty(data{trial,step_Pulse}))
                    dexTrialNonzero(count) = trial;
                end
                count = count+1;
            end
            
            for i = 1:length(dexTrialNonzero)  
                clear tmpData
                tmpData = data{dexTrialNonzero(i),step_Pulse};
                this.k_hat_pulse(i) = get_singleTrial_k_hat_pulse(this,tmpData.f(2,:),...
                                                 tmpData.Fp(2,:),...
                                                 tmpData.x(2,:),...
                                                 tmpData.ts);
            end
                        
        end
        
        function [] = get_k_hat_stocastic(this,data)
            
            % Exclude empty cells (FIX Later)
            step_Stocastic = 3;
            [N_trial] = size(data,1);
            count = 1;
            for trial = 1:N_trial
                if(~isempty(data{trial,step_Stocastic}))
                    dexTrialNonzero(count) = trial;
                end
                count = count+1;
            end
            
            this.k_hat_stocastic = [NaN, NaN, NaN, NaN, NaN];
            for i = 1:length(dexTrialNonzero)  
                clear tmpData
                tmpData = data{dexTrialNonzero(i),step_Stocastic};
                this.k_hat_stocastic(i) = get_singleTrial_k_hat_stocastic(this,tmpData.Fp(2,:),...
                                                 tmpData.x(2,:),...
                                                 tmpData.ts);
            end
            
        end
        
        function [] = get_k_hat_release(this,data)
               
            % Exclude empty cells (FIX Later)
            release = 1;
            [N_trial] = size(data,1);
            count = 1;
            for trial = 1:N_trial
                if(~isempty(data{trial,release}))
                    dexTrialNonzero(count) = trial;
                end
                count = count+1;
            end
            
            for i = 1:length(dexTrialNonzero)  
                clear tmpData
                tmpData = data{dexTrialNonzero(i),release};
                this.k_hat_release(i) = get_singleTrial_k_hat_release(this,tmpData.x(2,:),tmpData.ts);
            end
            
        end
        
        function [k_hat] = get_singleTrial_k_hat_pulse(this,f,Fp,x,ts)
            
            % Take average of last 15 measuremnts before the pulse ends
            dexPulseStart = min(find(Fp~=0));
            dexRange1 = [(dexPulseStart-15):dexPulseStart];
            
            dexPulseEnd = max(find(Fp~=0));
            dexRange2 = [(dexPulseEnd-15):dexPulseEnd];
            
            f1 = mean(f(dexRange1));
            f2 = mean(f(dexRange2));
            
            x1 = mean(x(dexRange1));
            x2 = mean(x(dexRange2));
            
            k_hat = -(f2-f1)/(x2-x1);
            
            % Look at raw data
%             figure;
%             ax1 = subplot(4,1,1); plot(f); hold on; 
% %             plot(dexPulseStart,f1,'+',dexPulseEnd,f2,'*'); 
%             
%             ax2 = subplot(4,1,2); plot(x); hold on; 
% %             plot(dexPulseStart,x1,'+',dexPulseEnd,x2,'*'); 
% 
%             ax3 = subplot(4,1,3); plot(Fp);
%             ax4 = subplot(4,1,4); plot(ts);
%             linkaxes([ax1,ax2,ax3,ax4],'x');
            
        end
    
        function [k_hat] = get_singleTrial_k_hat_stocastic(this,X1,Y1,ts)
        
%             figure; plot(X1);
%             figure; plot(Y1);
%             figure; plot(ts);
            
            dexStart = min(find(ts == 3)) + 500;
            dexEnd = max(find(ts == 3));
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
            
%             figure;
%             subplot(2,1,1); plot(Y1);
            
            Y1 = detrend(Y1);
            
            TF_freq = Hz/2*linspace(0,1,nfft/2+1);
            
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
            
            %% Impedance plot (Diagonal)
            xLowerLim = 0.5;
            xUpperLim = 50.0;
            yLowerLim11 = 1e+3;
            yUpperLim11 = 1e+6;
            yLowerLim22 = 2e+3;
            yUpperLim22 = 3e+4;
            
%             figure(1); hold on;
%             % Magnitude plot of ankle impedance
%             ax1 = subplot(2,1,1,'XScale','log','YScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_mag(:,1),'LineWidth',2); grid on; box on;
% %             axis([xLowerLim xUpperLim yLowerLim11 yUpperLim11]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
%             ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
%             title('Z11','fontWeight','bold','fontSize',16);
%             
%             ax2 = subplot(2,1,2,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_phi(:,1),'LineWidth',2); grid on; box on;
% %             axis([xLowerLim xUpperLim 0 180]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
%             
%             linkaxes([ax1,ax2],'x');
%             
%             %% Partial Coherence Plot
%             xLowerLim = 0.5;
%             figure(2); hold on;
%             set(gcf,'Color',[1,1,1]);
%             
%             ax1 = subplot(1,1,1,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(OCx1y1_F,OCx1y1,'LineWidth',2);
%             grid on;box on; ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
%             xlabel('Hz','fontWeight','bold','fontSize',14);
%             title('Y11 PC','fontWeight','bold','fontSize',16);
            
            [tmp,dexK] = min(abs(TF_freq - 1));
            k_hat = Z11_mag(dexK);
            
        end
        
        function [k_hat] = get_singleTrial_k_hat_release(this,x,ts)
           
            dexStart = min(find(ts== 4));
            dexEnd = max(find(ts == 5));
            [h_model, M, B, K, VAFirf] =  this.fitModel(x(dexStart:dexEnd), dexEnd-dexStart);
            
        end
        
        function [h_model, M, B, K, VAFirf] =  fitModel(this, h_hat, L)
            
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
            
            figure;
            plot(h_hat,'o'); hold on;
            plot(h_model,'.');
            
            %         figure;
            %         for i = iDex
            %             plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_hat_MA(:,i),'.r'); hold on;
            %             plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model(:,i),'-b'); hold on;
            %
            %         end
            %         %             plot3(M2,L,X(:,1)','--k','linewidth',2.5);
            %         xlabel('lag (s)'); % xlim([0 0.2]);
            %         ylabel('Time (s)');
            %         zlabel('x (m)');
            %         set(gca,'fontsize',18); grid on;
            
            %         for i = iDex
            %             for j =
            %              z_tmp() = cov(h_model(i,:),u_r());
            %             end
            %              z_r_hat(i,:) = z_r_hat(i,:) + z_tmp
            %         end
            
            % VAF IRF Calculation
            VAFirf = zeros(N,1);
            
            for i = iDex
                VAFirf(i,1) = 100*(1-std(h_hat_MA(1:M2-M1,i)-h_model(1:M2-M1,i))^2/std(h_hat_MA(1:M2-M1,i))^2);
            end
            
            %         figure; plot(this.tvec(iDex),VAFirf(iDex),'.-','linewidth',2.5,'markersize',15); grid on; xlim([0.2 this.tvec(end)]); ylim([0 100]);
            %         ylabel('%VAF_{IRF}'); xlabel('Time (s)'); set(gca,'fontsize',18);
            
            %         for r = 1:R
            %             for i = M2+1:N+M1
            %                 % Add same code for z
            %                 for j = M1:M2
            %                     y_tmp(j+1) = h_model(j+1,i)*u_r(r,i-j);
            %                 end
            %                 y(r,i) = sum(y_tmp);
            %             end
            %         end
            
            
            
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


