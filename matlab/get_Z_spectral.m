classdef get_Z_spectral < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        
    end
    
    methods
        function [this] = get_Z_spectral(Data_pert,col_vec)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % SS Dynamic Ankle Impedance Identification (Relaxed)
            % Impedance Identification in principal axis directions (DP and IE)
            
            
            % Load Anklebot Impedance Model: Measured separately without a human subject
            % load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/Stocastic/AnklebotImpedanceModel.mat');
            
            
            %% Variables
            nfft=12000;
            nFreq = nfft/2+1;
            Hz = 500;
%             runningTime = this.runTime;
            
            unwrapThreshold = -50*pi/180;
            
            % Admittance of Anklebot + Ankle
            Y11_s = zeros(nFreq,1);
            Y22_s = zeros(nFreq,1);
            Y12_s = zeros(nFreq,1);
            Y21_s = zeros(nFreq,1);
            
            % Impedance of Anklebot + Ankle
            Z11_s = zeros(nFreq,1);     Z11_s_mag = zeros(nFreq,1);     Z11_s_phi = zeros(nFreq,1);
            Z22_s = zeros(nFreq,1);     Z22_s_mag = zeros(nFreq,1);     Z12_s_phi = zeros(nFreq,1);
            Z12_s = zeros(nFreq,1);     Z12_s_mag = zeros(nFreq,1);     Z21_s_phi = zeros(nFreq,1);
            Z21_s = zeros(nFreq,1);     Z21_s_mag = zeros(nFreq,1);     Z22_s_phi = zeros(nFreq,1);
            
            % Impedance of Ankle
            Z11_l = zeros(nFreq,1);     Z11_l_mag = zeros(nFreq,1);     Z11_l_phi = zeros(nFreq,1);
            Z22_l = zeros(nFreq,1);     Z22_l_mag = zeros(nFreq,1);     Z12_l_phi = zeros(nFreq,1);
            Z12_l = zeros(nFreq,1);     Z12_l_mag = zeros(nFreq,1);     Z21_l_phi = zeros(nFreq,1);
            Z21_l = zeros(nFreq,1);     Z21_l_mag = zeros(nFreq,1);     Z22_l_phi = zeros(nFreq,1);
            
            PC11_s = zeros(nFreq,1);
            PC22_s = zeros(nFreq,1);
            PC12_s = zeros(nFreq,1);
            PC21_s = zeros(nFreq,1);            
            
            % Anklebot Data
            % data_abot = load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/Stocastic/SampleData.asc', 'r');
            % data_abot_length = length(data_abot);
            % cutOff_abot = data_abot_length-runningTime*Hz;
            
            % joint info
            % d_DP=data_abot(cutOff_abot:end-1,4);
            % d_IE=data_abot(cutOff_abot:end-1,3);
            % t_DP=data_abot(cutOff_abot+1:end,8);
            % t_IE=data_abot(cutOff_abot+1:end,7);
            
            
            % robot info
            % d_DP=data_abot(cutOff_abot:end-1,10);
            % d_IE=data_abot(cutOff_abot:end-1,13);
            % t_DP=data_abot(cutOff_abot+1:end,11);
            % t_IE=data_abot(cutOff_abot+1:end,14);
            
            % Import data and convert to X1, X2, Y1, Y2
            
            % Sanity check
%             figure;
%             plot(U_X_N); hold on;
            % Select points
%             [dex_X,dex_Y] = ginput(1);
%             dex_X = 15*Hz:(Hz*(15+runningTime));
%             dex_X = [round(dex_X):(round(dex_X)+Hz*runningTime)];
%             plot([dex_X(1), dex_X(end)],[0,0],'ok','linewidth',2.5);
            
            
%             cf = 500; % cutoff freqnency
%             [b,a] = butter(3,cf/(Hz/2)); % make filter
                  
%             X1=Data_pert.cf(:,1);     X2=Data_pert.cf(:,2);
%             Y1=Data_pert.tp(:,1);     Y2=Data_pert.tp(:,2);
            %Y1=filtfilt(b,a,X_m(dex_X));     Y2=filtfilt(b,a,Y_m(dex_X));
            
            %% Try combining all
            X1 = [];
            X2 = [];
            Y1 = [];
            Y2 = [];
            for i = 1:30
                X1 = [X1; Data_pert(i).cf(:,1)];
                X2 = [X2; Data_pert(i).cf(:,2)];
                Y1 = [Y1; detrend(Data_pert(i).tp(:,1))];
                Y2 = [Y2; detrend(Data_pert(i).tp(:,2))];
            end
            % Check xcorr preturbations
            % [c,lags] = xcorr(X1,X1,220000,'coeff'); title('xcorr(X1,X1)');
            % figure; subplot(3,1,1); stem(lags,c); hold on;
            % [c,lags] = xcorr(X2,X2,220000,'coeff'); title('xcorr(X2,X2)');
            % subplot(3,1,2); stem(lags,c);
            % [c,lags] = xcorr(X1,X2,220000,'coeff'); title('xcorr(X1,X2)');
            % subplot(3,1,3); stem(lags,c);
            
            figure; 
            subplot(2,1,1); plot(Y1);
            subplot(2,1,2); plot(Y2);
            
            Y1 = detrend(Y1);
            Y2 = detrend(Y2);
            
            TF_freq = Hz/2*linspace(0,1,nfft/2+1);
            
            % Cross power spectral density
            [Px1x1,Fx1x1] = cpsd(X1,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Px1x2,Fx1x2] = cpsd(X2,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Px1y1,Fx1y1] = cpsd(Y1,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Px1y2,Fx1y2] = cpsd(Y2,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            
            [Px2x1,Fx2x1] = cpsd(X1,X2,hamming(nfft/2),nfft/4,nfft,Hz);
            [Px2x2,Fx2x2] = cpsd(X2,X2,hamming(nfft/2),nfft/4,nfft,Hz);
            [Px2y1,Fx2y1] = cpsd(Y1,X2,hamming(nfft/2),nfft/4,nfft,Hz);
            [Px2y2,Fx2y2] = cpsd(Y2,X2,hamming(nfft/2),nfft/4,nfft,Hz);
            
            [Py1x1,Fy1x1] = cpsd(X1,Y1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Py1x2,Fy1x2] = cpsd(X2,Y1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Py1y1,Fy1y1] = cpsd(Y1,Y1,hamming(nfft/2),nfft/4,nfft,Hz);
            [Py1y2,Fy1y2] = cpsd(Y2,Y1,hamming(nfft/2),nfft/4,nfft,Hz);
            
            [Py2x1,Fy2x1] = cpsd(X1,Y2,hamming(nfft/2),nfft/4,nfft,Hz);
            [Py2x2,Fy2x2] = cpsd(X2,Y2,hamming(nfft/2),nfft/4,nfft,Hz);
            [Py2y1,Fy2y1] = cpsd(Y1,Y2,hamming(nfft/2),nfft/4,nfft,Hz);
            [Py2y2,Fy2y2] = cpsd(Y2,Y2,hamming(nfft/2),nfft/4,nfft,Hz);
            
            %% Sanity Check
            
%             % Sanity check xx
%             figure;
%             subplot(2,2,1); loglog(Fx1x1,abs(Px1x1)); title('x_1 x_1'); grid on;
%             subplot(2,2,2); loglog(Fx1x2,abs(Px1x2)); title('x_1 x_2'); grid on;
%             subplot(2,2,3); loglog(Fx2x1,abs(Px2x1)); title('x_2 x_1'); grid on;
%             subplot(2,2,4); loglog(Fx2x2,abs(Px2x2)); title('x_2 x_2'); grid on;
%             
%             % Sanity check yy
%             figure;
%             subplot(2,2,1); loglog(Fy1y1,abs(Py1y1)); title('y_1 y_1'); grid on;
%             subplot(2,2,2); loglog(Fy1y2,abs(Py1y2)); title('y_1 y_2'); grid on;
%             subplot(2,2,3); loglog(Fy2y1,abs(Py2y1)); title('y_2 y_1'); grid on;
%             subplot(2,2,4); loglog(Fy2y2,abs(Py2y2)); title('y_2 y_2'); grid on;
%             
%             % Sanity check xy
%             figure;
%             subplot(2,2,1); loglog(Fx1y1,abs(Px1y1)); title('x_1 y_1'); grid on;
%             subplot(2,2,2); loglog(Fx1y2,abs(Px1y2)); title('x_1 y_2'); grid on;
%             subplot(2,2,3); loglog(Fx2y1,abs(Px2y1)); title('x_2 y_1'); grid on;
%             subplot(2,2,4); loglog(Fx2y2,abs(Px2y2)); title('x_2 y_2'); grid on;
%             
%             % Sanity check yx
%             figure;
%             subplot(2,2,1); loglog(Fy1x1,abs(Py1x1)); title('y_1 x_1'); grid on;
%             subplot(2,2,2); loglog(Fy1x2,abs(Py1x2)); title('y_1 x_2'); grid on;
%             subplot(2,2,3); loglog(Fy2x1,abs(Py2x1)); title('y_2 x_1'); grid on;
%             subplot(2,2,4); loglog(Fy2x2,abs(Py2x2)); title('y_2 x_2'); grid on;
            
            %% Ordinary Coherence
            [OCx1x2,OCx1x2_F] = mscohere(X2,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            [OCx2x1,OCx2x1_F] = mscohere(X1,X2,hamming(nfft/2),nfft/4,nfft,Hz);
            [OCx2y1,OCx2y1_F] = mscohere(Y1,X2,hamming(nfft/2),nfft/4,nfft,Hz);
            [OCx1y2,OCx1y2_F] = mscohere(Y2,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            [OCx1y1,OCx1y1_F] = mscohere(Y1,X1,hamming(nfft/2),nfft/4,nfft,Hz);
            [OCx2y2,OCx2y2_F] = mscohere(Y2,X2,hamming(nfft/2),nfft/4,nfft,Hz);
            
            %             figure;
            %             subplot(4,4,2); semilogx(OCx1x2_F,OCx1x2); title('C: x_1 x_2');grid on;
            %             subplot(4,4,5); semilogx(OCx2x1_F,OCx2x1); title('C: x_2 x_1');grid on;
            %             subplot(4,4,7); semilogx(OCx2y1_F,OCx2y1); title('C: x_2 y_1');grid on;
            %             subplot(4,4,4); semilogx(OCx1y2_F,OCx1y2); title('C: x_1 y_2');grid on;
            %             subplot(4,4,3); semilogx(OCx1y1_F,OCx1y1); title('C: x_1 y_1');grid on;
            %             subplot(4,4,8); semilogx(OCx2y2_F,OCx2y2); title('C: x_2 y_2');grid on;

%             for j=1:nFreq
%                 Y11_s(j,1) = Px1y1(j,1)/Px1x1(j,1);
%             end
%             
%             for j=1:nFreq
%             Y11_s_mag(j,1) = 20*log10(abs(Y11_s(j,1)));
%             Y11_s_phi(j,1) = 180/pi*unwrap2(angle(Y11_s(j,1)),unwrapThreshold,'up');
%             end
%             
%             figure;plot(TF_freq,Y11_s_mag);
            
            for j=1:nFreq
                Y11_s(j,1) = (Px1y1(j,1)/Px1x1(j,1))*(1-(Px1x2(j,1)*Px2y1(j,1))/(Px2x2(j,1)*Px1y1(j,1)))/(1-OCx1x2(j,1));
                Y12_s(j,1) = (Px2y1(j,1)/Px2x2(j,1))*(1-(Px2x1(j,1)*Px1y1(j,1))/(Px1x1(j,1)*Px2y1(j,1)))/(1-OCx1x2(j,1));
                Y21_s(j,1) = (Px1y2(j,1)/Px1x1(j,1))*(1-(Px1x2(j,1)*Px2y2(j,1))/(Px2x2(j,1)*Px1y2(j,1)))/(1-OCx1x2(j,1));
                Y22_s(j,1) = (Px2y2(j,1)/Px2x2(j,1))*(1-(Px2x1(j,1)*Px1y2(j,1))/(Px1x1(j,1)*Px2y2(j,1)))/(1-OCx1x2(j,1));
                
                % Partial Coherence
                PC11_s(j,1) = abs(Px1y1(j,1)*Px2x2(j,1)-Px2y1(j,1)*Px1x2(j,1))^2/(Px2x2(j,1)*Px2x2(j,1)*Px1x1(j,1)*Py1y1(j,1))/(1-OCx2x1(j,1))/(1-OCx2y1(j,1));
                PC21_s(j,1) = abs(Px1y2(j,1)*Px2x2(j,1)-Px2y2(j,1)*Px1x2(j,1))^2/(Px2x2(j,1)*Px2x2(j,1)*Px1x1(j,1)*Py2y2(j,1))/(1-OCx2x1(j,1))/(1-OCx2y2(j,1));
                PC12_s(j,1) = abs(Px2y1(j,1)*Px1x1(j,1)-Px1y1(j,1)*Px2x1(j,1))^2/(Px1x1(j,1)*Px1x1(j,1)*Px2x2(j,1)*Py1y1(j,1))/(1-OCx1x2(j,1))/(1-OCx1y1(j,1));
                PC22_s(j,1) = abs(Px2y2(j,1)*Px1x1(j,1)-Px1y2(j,1)*Px2x1(j,1))^2/(Px1x1(j,1)*Px1x1(j,1)*Px2x2(j,1)*Py2y2(j,1))/(1-OCx1x2(j,1))/(1-OCx1y2(j,1));
                
                % Anklebot+Ankle Impedance
                K_matrix = inv([Y11_s(j,1) Y12_s(j,1); Y21_s(j,1) Y22_s(j,1)]);
                
                Z11_s(j,1) = K_matrix(1,1);
                Z12_s(j,1) = K_matrix(1,2);
                Z21_s(j,1) = K_matrix(2,1);
                Z22_s(j,1) = K_matrix(2,2);
                
                % Ankle Impedance
                Z11_l(j,1) = Z11_s(j,1);% - Z11_abot(j,1);
                Z22_l(j,1) = Z22_s(j,1);% - Z22_abot(j,1);
                Z12_l(j,1) = Z12_s(j,1);% - Z12_abot(j,1);
                Z21_l(j,1) = Z21_s(j,1);% - Z21_abot(j,1);
            end
            
            %% Mag & Phase calculation for Impedance plot
            for j=1:nFreq
                Z11_s_mag(j,1) = abs(Z11_s(j,1));
                Z22_s_mag(j,1) = abs(Z22_s(j,1));
                Z12_s_mag(j,1) = abs(Z12_s(j,1));
                Z21_s_mag(j,1) = abs(Z21_s(j,1));
                Z11_s_phi(j,1) = 180/pi*unwrap2(angle(Z11_s(j,1)),unwrapThreshold,'up');
                Z22_s_phi(j,1) = 180/pi*unwrap2(angle(Z22_s(j,1)),unwrapThreshold,'up');
                Z12_s_phi(j,1) = 180/pi*unwrap2(angle(Z12_s(j,1)),unwrapThreshold,'up');
                Z21_s_phi(j,1) = 180/pi*unwrap2(angle(Z21_s(j,1)),unwrapThreshold,'up');
                
                Z11_l_mag(j,1) = abs(Z11_l(j,1));
                Z22_l_mag(j,1) = abs(Z22_l(j,1));
                Z12_l_mag(j,1) = abs(Z12_l(j,1));
                Z21_l_mag(j,1) = abs(Z21_l(j,1));
                Z11_l_phi(j,1) = 180/pi*unwrap2(angle(Z11_l(j,1)),unwrapThreshold,'up');
                Z22_l_phi(j,1) = 180/pi*unwrap2(angle(Z22_l(j,1)),unwrapThreshold,'up');
                Z12_l_phi(j,1) = 180/pi*unwrap2(angle(Z12_l(j,1)),unwrapThreshold,'up');
                Z21_l_phi(j,1) = 180/pi*unwrap2(angle(Z21_l(j,1)),unwrapThreshold,'up');
            end
            
            %% Impedance plot (Diagonal)
            xLowerLim = 0.5;
            xUpperLim = 50.0;
            yLowerLim11 = 2e+3;
            yUpperLim11 = 1e+6;
            yLowerLim22 = 2e+3;
            yUpperLim22 = 3e+4;
            
            figure(1); hold on;
            % Magnitude plot of ankle impedance
            ax1 = subplot(2,2,1,'XScale','log','YScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z11_l_mag(:,1),'LineWidth',2,'Color', col_vec); grid on; box on;
            axis([xLowerLim xUpperLim yLowerLim11 yUpperLim11]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); 
            ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
            title('Z11','fontWeight','bold','fontSize',16);
            
            ax2 = subplot(2,2,2,'XScale','log','YScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z22_l_mag(:,1),'LineWidth',2,'Color', col_vec); grid on; box on;
            axis([xLowerLim xUpperLim yLowerLim22 yUpperLim22]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); 
            ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
            title('Z22','fontWeight','bold','fontSize',16);
            
            ax3 = subplot(2,2,3,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z11_l_phi(:,1),'LineWidth',2,'Color', col_vec); grid on; box on;
            axis([xLowerLim xUpperLim 0 180]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
            
            ax4 = subplot(2,2,4,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z22_l_phi(:,1),'LineWidth',2,'Color', col_vec); grid on; box on;
            axis([xLowerLim xUpperLim 0 180]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
            
            linkaxes([ax1,ax2,ax3,ax4],'x');
            
            %% Partial Coherence Plot
            xLowerLim = 0.5;
            figure(2); hold on;
            set(gcf,'Color',[1,1,1]);
            
            ax1 = subplot(2,2,1,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC11_s(:,1),'LineWidth',2,'Color', col_vec);
            grid on;box on; ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14); 
            title('Y11 PC','fontWeight','bold','fontSize',16);
            
            ax2 = subplot(2,2,2,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC12_s(:,1),'LineWidth',2,'Color', col_vec); hold off;
            grid on; box on; ylim([0 1]);
            axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14); 
            title('Y12 PC','fontWeight','bold','fontSize',16);
            
            ax3 = subplot(2,2,3,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC21_s(:,1),'LineWidth',2,'Color', col_vec); hold off;
            grid on;box on;ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14);
            title('Y21 PC','fontWeight','bold','fontSize',16);
            
            ax4 = subplot(2,2,4,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC22_s(:,1),'LineWidth',2,'Color', col_vec); hold off;
            grid on;box on;ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14); 
            title('Y22 PC','fontWeight','bold','fontSize',16);
            
            linkaxes([ax1,ax2,ax3,ax4],'x');
            
            %% Subtract Robot
            xLowerLim = 0.5;
            xUpperLim = 45.0;
            yLowerLim_sub = 5e+2;
            yUpperLim_sub = 1e+4;
            
            load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/model_robot.mat');
            model_mag = interp1(model_robot(:,1),model_robot(:,2),TF_freq,'linear')';
            model_phi = interp1(model_robot(:,1),model_robot(:,3),TF_freq,'linear')';

            figure(3); hold on;
            set(gcf,'Color',[1,1,1]);
                        
            ax2 = subplot(2,1,1,'XScale','log','YScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z22_l_mag(:,1)-model_mag,'LineWidth',2,'Color', col_vec); 
            grid on; box on;
            axis([xLowerLim xUpperLim yLowerLim_sub yUpperLim_sub]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); 
            ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
            title('Z22','fontWeight','bold','fontSize',16);
            
            ax4 = subplot(2,1,2,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z22_l_phi(:,1)-model_phi,'LineWidth',2,'Color', col_vec); grid on; box on;
            axis([xLowerLim xUpperLim -90 45 ]); yticks([-90: 45 :45])
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); 
            ylabel('phase (deg)','fontWeight','bold','fontSize',14);
            
            linkaxes([ax2,ax4],'x');
            
%             % Sanity check
%             opts = bodeoptions('cstprefs');
%             opts.MagUnits = 'abs';
%             opts.MagScale = 'log';
%             s = tf('s');
%             
%             figure;
%             m = 5;
%             b = 10;
%             k = 2500;
%             [MAG,PHASE,W] = bode(m*s^2 + b*s + k);
%             MAG = squeeze(MAG);
%             PHASE = squeeze(PHASE);
%             subplot(2,1,1);loglog(W,MAG,'-','linewidth',2.5); grid on; %xlim([W(2) W(end)]); ylim([10^-6 10^-1]);
%             ylabel('Amplitude (abs)'); set(gca,'fontsize',18);
% %             yticks([1e-5,1e-3,1e-1]);
%             
%             subplot(2,1,2); semilogx(W,PHASE,'linewidth',2.5); grid on;
%             ylabel('Phase (deg)'); xlabel('Frequency ()'); set(gca,'fontsize',18);
            

            
        end
    end
end


