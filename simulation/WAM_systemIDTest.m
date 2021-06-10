classdef WAM_systemIDTest < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sfrq
        dt
        tvec
        X_0
        U_tmp
    end
    
    methods
        function this = WAM_systemIDTest()
            dbstop if error
            
            disp('%%%%%%%%%%%%%Beginning of systeIDTest()%%%%%%%%%%%%%');
%             this.get_1DSim();
%             this.get_2DSim_linear();
            this.get_2DSim_WAM_SISO();      
%             this.get_2DSim_WAM_MIMO(); 
%             this.get_2DSim_circle_MIMO_spectral();
            
        end
        
        %% Get diffrent simulation types
        
        function [] = get_2DSim_linear(this)
            
            ensambleLength = 2;
            
            % Define sim parameters
            this.sfrq = 500; % Change later this will change window size
            this.dt = 1/this.sfrq;
            this.tvec = 0:this.dt:ensambleLength+this.dt;
            N = length(this.tvec); % Check N is even
            lagBuffer = 5;
            M1 = 0;
            M2 = 0.7*this.sfrq+lagBuffer;
            L = M2-M1+1;
            R = 300; % number of realizations
            
            %             % Intialize variables
            y1_r = zeros(R,N); % Postion
            y2_r = zeros(R,N); % Postion
            
            z1_r = zeros(R,N); % Postion
            z2_r = zeros(R,N); % Postion
            
            u1_r = zeros(R,N); % Preturbation
            u2_r = zeros(R,N); % Preturbation
            
            % Generate simulated data for ensamble method
            pw = PoolWaitbar(R, ['Simulating System ',int2str(R),' times.']);
            parfor i = 1:R
                [y_tmp,u_tmp] = this.get_sim_linear();
                y1_r(i,:) = y_tmp(1,:); % CHECK NORMALIZE
                y2_r(i,:) = y_tmp(2,:); % CHECK NORMALIZE
                u1_r(i,:) = u_tmp(1,:);
                u2_r(i,:) = u_tmp(2,:);
                
                increment(pw); % wait bar display
            end
%             
%             
%                         save('/Users/jhermus/Desktop/testSim.mat','y1_r','y2_r','u1_r','u2_r');
%                          load('/Users/jhermus/Desktop/testSim.mat');
            
            figure; plot(this.tvec,u1_r(1,:),'linewidth',2.5);
            ylim([-7 7]);xlim([0 0.2]);
            xlabel('Time (s)'); ylabel('U (N)');
            set(gca,'fontsize',18); grid on;
            
            % Add noise
            n1_r = 0; %1e-4*randn(size(y_r));
            n2_r = 0;
            %                 snr(y_r,n_r)
            z1_r = y1_r + n1_r;
            z2_r = y2_r + n2_r;
            
            clear y1_r y2_r
            
            %            % Filter
            %            cf = 50; % cutoff freqnency
            %            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            %            for i = 1:R
            %                z1_r_temp(i,:) = filtfilt(b,a,z1_r(i,:)); % apply fitler
            %                z2_r_temp(i,:) = filtfilt(b,a,z2_r(i,:)); % apply fitler
            %            end
            %
            %            figure;
            %            plot(this.tvec,z1_r(100,:)); hold on;
            %            plot(this.tvec,z1_r_temp(100,:));
            %
            %            figure;
            %            plot(this.tvec,z2_r(100,:)); hold on;
            %            plot(this.tvec,z2_r_temp(100,:));
            %
            %             figure;
            %             plot(this.tvec,z1_r'); hold on; xlim([0 this.tvec(end)]);
            %             plot(this.tvec,median(z1_r),'-k','linewidth',2.5);
            %             xlabel('Time (s)'); ylabel('x (m)');
            %             set(gca,'fontsize',18); grid on;
            %
            %             for i = 1:R
            %                 [w,P(i,:)] = this.get_fft(z1_r(i,:),1);
            %                 [w,P_temp(i,:)] = this.get_fft(z1_r_temp(i,:),1);
            %
            %             end
            %             figure; semilogx(w,sum(P,1),w,sum(P_temp,1));
            
            % Change sampling rate to 500 Hz
            %                 [z_r, u_r, N, L, R, M1, M2] = this.get_lowerSamplingRate(ensambleLength, z_r, u_r, N, L, R, M1, M2);
            
            %             z1_r_mean = mean(z1_r);
            %             z2_r_mean = mean(z2_r);
            %
            %             for i = 1:R
            %                 z1_r(i,:) = detrend(z1_r(i,:) - z1_r_mean);
            %                 z2_r(i,:) = detrend(z2_r(i,:) - z2_r_mean);
            %             end
            
            %             figure;
            %             plot(this.tvec,z1_r); hold on; xlim([0 this.tvec(end)]);
            %             xlabel('Time (s)'); ylabel('x (m)');
            %             set(gca,'fontsize',18); grid on;
            
%             get_Z_spectral(this,z1_r(1,:), z2_r(1,:), u1_r(1,:), u2_r(1,:));

            
            [H_hat_MA, MAwindow] = ensambleSysID_Matrix(this, z1_r, z2_r, u1_r, u2_r, N, L, R, M1, M2);
            
            %             [h_hat_MA, MAwindow] = ensambleSysID_windowed(this, z1_r, u1_r, N, L, R, M1, M2);
            
            iDex = floor(linspace(M2+1+40,N+M1-3,15));
            %             [H_model, M_model, B_model, K_model] = this.fitModel_Matrix(H_hat_MA, N, L, M1, M2, iDex);
            
            
            %             Y_model11 = tf(1,[M(1,1),B(1,1),K(1,1)]);
            %             h_model11 = impulse(Y_model11,0:1/this.sfrq:M2/this.sfrq);
            %
            %             Y_model12 = tf(1,[M(1,2),B(1,2),K(1,2)]);
            %             h_model12 = impulse(Y_model12,0:1/this.sfrq:M2/this.sfrq);
            %
            %             Y_model21 = tf(1,[M(2,1),B(2,1),K(2,1)]);
            %             h_model21 = impulse(Y_model21,0:1/this.sfrq:M2/this.sfrq);
            %
            %             Y_model22 = tf(1,[M(2,2),B(2,2),K(2,2)]);
            %             h_model22 = impulse(Y_model22,0:1/this.sfrq:M2/this.sfrq);
            
            %% Known model
            M = [1.7099, -0.2566; -0.2566, 2.1775];
            B = [5.2510, -1.0215; -1.0215, 39.0782];
            K = [105.0196, -20.4292; -20.4292, 781.5645];
            
            s = tf('s');
            %             syms a b c d s
            %             syms m11 m12 m21 m22
            %             syms b11 b12 b21 b22
            %             syms k11 k12 k21 k22
            %
            %             M = [m11 m12; m21 m22];
            %             B = [b11 b12; b21 b22];
            %             K = [k11 k12; k21 k22];
            
            Z = M*s^2+B*s+K;
            %             inv([a,b;c,d])
            
            figure; impulse(inv(Z));
            Y_model = inv(Z);
            H_true_model = impulse(Y_model,0:1/this.sfrq:M2/this.sfrq);
            
            %             [b,a] = prony(H_true_model(:,1,1),2,4);
            %             [w,P] = this.get_fft(H_true_model(:,1,1),1);
            %             figure; plot(w,P); grid on; hold on;
            %             xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
            
            %% Plot
            figure;
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3);
            ax4 = subplot(2,2,4);
            
            for i = iDex(1:end-1)
                axes(ax1);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_hat_MA(:,i,1),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_true_model(:,1,1),'-k'); hold on;
                %                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_model(:,1,1),'-b'); hold on;
                
                
                axes(ax2);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_hat_MA(:,i,2),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_true_model(:,1,2),'-k'); hold on;
                %                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_model(:,1,2),'-b'); hold on;
                
                
                axes(ax3);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_hat_MA(:,i,3),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_true_model(:,2,1),'-k'); hold on;
                %                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_model(:,2,1),'-b'); hold on;
                
                axes(ax4);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_hat_MA(:,i,4),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_true_model(:,2,2),'-k'); hold on;
                %                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_model(:,2,2),'-b'); hold on;
                
            end
            fs = 15;
            axes(ax1); view(35,25); title('11'); xlim([0 M2/this.sfrq]); zlim([-0.1 0.1]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('x (m)'); grid on; set(gca,'fontsize',fs);
            axes(ax2); view(35,25); title('12'); xlim([0 M2/this.sfrq]); zlim([-0.005 0.005]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('y (m)'); grid on; set(gca,'fontsize',fs);
            axes(ax3); view(35,25); title('21'); xlim([0 M2/this.sfrq]); zlim([-0.005 0.005]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('y (m)'); grid on; set(gca,'fontsize',fs);
            axes(ax4); view(35,25); title('22'); xlim([0 M2/this.sfrq]); zlim([-0.03 0.03]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('y (m)'); grid on; set(gca,'fontsize',fs);
            
            tmp = mean(H_hat_MA(:,365:end-5,1),2);
            tmp = tmp(1:end-5);
            
%             [b,a] = prony(H_hat_MA(:,250,1),4,4);
            [b,a] = prony(tmp(1:end),4,4);
            C_hat_d11 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c11 = d2c(C_hat_d11,'tustin');
            H_hat_c11 = impulse(C_hat_c11,0:1/this.sfrq:M2/this.sfrq);
            
            figure; plot(0:1/this.sfrq:M2/this.sfrq - 5/this.sfrq,tmp,0:1/this.sfrq:M2/this.sfrq,H_hat_c11);
            
            [b,a] = prony(H_hat_MA(:,250,2),4,4);
            C_hat_d12 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c12 = d2c(C_hat_d12,'tustin');
            H_hat_c12 = impulse(C_hat_c12,0:1/this.sfrq:M2/this.sfrq);
            
            figure; plot(0:1/this.sfrq:M2/this.sfrq,H_hat_MA(:,250,2),0:1/this.sfrq:M2/this.sfrq,H_hat_c12);
            
            [b,a] = prony(H_hat_MA(:,250,3),4,4);
            C_hat_d21 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c21 = d2c(C_hat_d21,'tustin');
            H_hat_c21 = impulse(C_hat_c21,0:1/this.sfrq:M2/this.sfrq);
            
            figure; plot(0:1/this.sfrq:M2/this.sfrq,H_hat_MA(:,250,3),0:1/this.sfrq:M2/this.sfrq,H_hat_c21);
            
            [b,a] = prony(H_hat_MA(:,250,4),4,4);
            C_hat_d22 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c22 = d2c(C_hat_d22,'tustin');
            H_hat_c22 = impulse(C_hat_c22,0:1/this.sfrq:M2/this.sfrq);
            
            figure; plot(0:1/this.sfrq:M2/this.sfrq,H_hat_MA(:,250,4),0:1/this.sfrq:M2/this.sfrq,H_hat_c22);
            
            
            Z_hat = inv([C_hat_c11,C_hat_c12;C_hat_c12,C_hat_c22]);
            
            
            figure; pzplot(Z_hat);
            figure;
            bode(Z_hat); hold on;
            bode(Z);xlim([10^-1 10^2]);
            
            %%
            syms s m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22
            
            M = [m11 m12; m21 m22];
            B = [b11 b12; b21 b22];
            K = [k11 k12; k21 k22];
            
            Z = M*s^2+B*s+K;
            C = inv(Z);
            
            [N,D] = numden(C);
            n(:,1,1) = flip(coeffs(N(1,1),s));
            d(:,1,1) = flip(coeffs(D(1,1),s));
            
            n(:,1,2) = flip(coeffs(N(1,2),s));
            d(:,1,2) = flip(coeffs(D(1,2),s));
            
            n(:,2,1) = flip(coeffs(N(2,1),s));
            d(:,2,1) = flip(coeffs(D(2,1),s));
            
            n(:,2,2) = flip(coeffs(N(2,2),s));
            d(:,2,2) = flip(coeffs(D(2,2),s));
            
            eqns = [n(:,1,1)./d(1,1,1) == C_hat_c11.numerator{1}(end-2:end)';...
                n(:,1,2)./d(1,1,1) == C_hat_c12.numerator{1}(end-2:end)';...
                n(:,2,1)./d(1,1,1) == C_hat_c21.numerator{1}(end-2:end)';...
                n(:,2,2)./d(1,1,1) == C_hat_c22.numerator{1}(end-2:end)'];
            %
            %             eqns = [d(2:end,1,1)/d(1,1,1) == C_hat_c11.Denominator{1}(2:end)';...
            %                     d(2:end,1,2)/d(1,1,2) == C_hat_c12.Denominator{1}(2:end)';...
            %                     d(2:end,1,2)/d(1,2,2) == C_hat_c22.Denominator{1}(2:end)'];
            
            vars = [m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22];
            
            [m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22] = solve(eqns,vars.');
            
            M = double([m11 m12; m21 m22]);
            B = double([b11 b12; b21 b22]);
            K = double([k11 k12; k21 k22]);

            disp('test');
            
            
        end
        
        function [] = get_2DSim_WAM_SISO(this)
            
            % To do list:
            % Convert to SISO
            % ??? PLAN
            
            movDuration = 1;
            R = 300;
            
           [y_r_1, y_r_2, u_r_1, u_r_2, f_r_1, f_r_2] = this.get_ensambles_from_simSystem2D(R,movDuration,'line_SISO');%,pretScale);            
%            save('/Users/jhermus/Desktop/TmpFigures_Ensamble/WAM_sim.mat','y_r_1','y_r_2','u_r_1','u_r_2','f_r_1','f_r_2','this');
%           load('/Users/jhermus/Desktop/TmpFigures_Ensamble/WAM_sim.mat');

            [z_r_1,z_r_2] = this.get_noise_subtract_x0(y_r_1, y_r_2, u_r_1, u_r_2, f_r_1, f_r_2);
            
            %% ID and fit tangential
            [R,N] = size(z_r_1);
            lagBuffer = 5;
            M1 = 0;
            M2 = 0.4*this.sfrq+lagBuffer;
            L = M2-M1+1;
            
            [H_hat_MA, MAwindow] = this.ensambleSysID_Matrix(z_r_1, z_r_2, u_r_1, u_r_2, N, L, R, M1, M2);
            h_hat_11 = H_hat_MA(:,:,1);
            h_hat_12 = H_hat_MA(:,:,2);
            h_hat_21 = H_hat_MA(:,:,3);
            h_hat_22 = H_hat_MA(:,:,4);
            
%             iDex = M2+1+40:N+M1-40; 
            iDex = floor(linspace(M2+1+40,N+M1-40,25));

            %% Model least squares fit
%             [K_par_hat, B_par_hat, M_par_hat] = get_parametricLeastSquares_1D(this,iDex,z_r_1, z_r_2, u_r_1, u_r_2);
                        
            %% Add filter
            cf = 15; % cutoff freqnency
            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            h_hat_11_filt = zeros(L,N);
            h_hat_12_filt = zeros(L,N);
            h_hat_21_filt = zeros(L,N);
            h_hat_22_filt = zeros(L,N);
            dexNonZero = M2-M1+1+MAwindow/2:N+M1-MAwindow/2;
            for i = dexNonZero
                h_hat_11_filt(:,i) = filtfilt(b,a,h_hat_11(:,i)); % apply fitler
                h_hat_12_filt(:,i) = filtfilt(b,a,h_hat_12(:,i)); % apply fitler
                h_hat_21_filt(:,i) = filtfilt(b,a,h_hat_21(:,i)); % apply fitler
                h_hat_22_filt(:,i) = filtfilt(b,a,h_hat_22(:,i)); % apply fitler
            end
            
            % Fit Model Prony
            M2_short = 0.4*this.sfrq+lagBuffer;
%             [K_xy_hat, h_model_11, h_model_12, h_model_21, h_model_22, VAFirf] = get_MIMO_Prony(this,iDex, h_hat_11_filt, h_hat_12_filt, h_hat_21_filt, h_hat_22_filt, M2,M2_short);
            [h_model_11, M_xy_hat, B_xy_hat, K_xy_hat,VAFirf_11] = this.fitModel(h_hat_11_filt, N, L, M1, M2_short, iDex);

            %% Estimate imporant values
%             K_ne_hat(1,1,:) =  K_ne_hat(1,1,:)-2500;
%             B_ne11 = B_ne11-40;

%             K_ne_hat = get_filtK(this,K_ne_hat,iDex,10);
                        
            test = WAM_simSystem2D(movDuration);                            
            F_contact = [mean(f_r_1); zeros(1,N)]; % Add tangential if robot dynamics are added
%             figure; plot(F_contact');                            
            s = tf('s');
            
            %% Make fit plots 
%             figure;
            % Estimate Joint impedance
            for i = iDex

                X = [mean(y_r_1(:,i));mean(y_r_2(:,i))];
                [q_1, q_2] = test.invKino(X(1),X(2));
                J = test.get_jacobian([q_1;q_2]);
                M_q = test.get_massMatrixArm([q_1;q_2]);
                
                [armPosVec_X, armPosVec_Y] = test.plotPosArm(X(1),X(2),'shoulder');
           
%                 ax = plot(armPosVec_X, armPosVec_Y, '-b',...
%                     armPosVec_X, armPosVec_Y, '.b','markersize',30,'linewidth',3);
%                 ylim([0 1.2]); xlim([-0.56 0.56]); grid on; axis equal;
%                 ylabel('Y Distance (m)','fontsize',14);
%                 xlabel('X Distance (m)','fontsize',14);
% %               legend('Arm','Zero-Force Trajectory','Crank','Constraint Path','Robot','location','northwest'); % Add ZFT
%                 ylabel('Y Distance (m)','fontsize',14);
%                 xlabel('X Distance (m)','fontsize',14);
%                 xlim([-0.6 0.6]);
                
                % Compute kinematic stiffness term
                parJ_parq1(:,:,i) = [-( test.l1a*cos(q_1) + test.l2a*cos(q_1 + q_2) ),  -test.l2a*cos(q_1 + q_2);...
                    -test.l1a*sin(q_1) - test.l2a*sin(q_1 + q_2), -test.l2a*sin(q_1 + q_2)];
                parJ_parq2(:,:,i) = [- test.l2a*cos(q_1 + q_2),  -test.l2a*cos(q_1 + q_2);...
                    - test.l2a*sin(q_1 + q_2), -test.l2a*sin(q_1 + q_2)];
                
                
                %% Expected Full Rank
                
                % Map joint stiffness to hand space use compliance mapping
                K_xy_const(:,:,i) = inv(J*inv(test.K)*J'); % - parJ_parq1(:,:,i)*F_contact(1,i) + parJ_parq2(:,:,i)*F_contact(2,i);
                B_xy_const(:,:,i) = inv(J*inv(test.B)*J');
                M_xy_const(:,:,i) = inv(J*inv(M_q)*J');
                
                if(test.t(i)<=2)
                  K_xy_const(:,:,i) = K_xy_const(:,:,i) + [2500,0;0,2500];
                  B_xy_const(:,:,i) = B_xy_const(:,:,i) + [40,0;0,40];
                end

                Z_true_model = M_xy_const(:,:,i)*s^2 + [B_xy_const(:,:,i)]*s + [K_xy_const(:,:,i)]; % Includes robot
           
%                 figure; impulse(inv(Z_true_model));
                Y_true_model = inv(Z_true_model);
                h_true_model(:,:,:,i) = impulse(Y_true_model,this.tvec(1:L));
                
                pause(0.01);
                
            end
            
            
            z_r_hat1 = zeros(size(z_r_1));
            z_r_hat2 = zeros(size(z_r_2));
            
            for r = 1:R
                for i = iDex
                    clear tmp1 tmp2
                    for j = 1:M2
                        
                        tmp1(j) = this.dt*h_hat_11_filt(j,i,:)*u_r_1(r,i-j) +...
                                  this.dt*h_hat_12_filt(j,i,:)*u_r_2(r,i-j);
                              
                        tmp2(j) = this.dt*h_hat_12_filt(j,i,:)*u_r_1(r,i-j) +...
                                  this.dt*h_hat_22_filt(j,i,:)*u_r_2(r,i-j);
                    end
                    
                    z_r_hat1(r,i) = sum(tmp1);
                    z_r_hat2(r,i) = sum(tmp2);
                    
                end
            end
            
             for r = 1:R
                VAF_output(1,r) = get_VAF(this,z_r_1(r,iDex(1):iDex(end)),z_r_hat1(r,iDex(1):iDex(end)));
                VAF_output(2,r) = get_VAF(this,z_r_2(r,iDex(1):iDex(end)),z_r_hat2(r,iDex(1):iDex(end)));
             end
                        
            for i = iDex
                VAF_true(1,i) = get_VAF(this,h_true_model(:,1,1,i),h_hat_11_filt(:,i));
                VAF_true(2,i) = get_VAF(this,h_true_model(:,1,2,i),h_hat_12_filt(:,i));
                VAF_true(3,i) = get_VAF(this,h_true_model(:,2,1,i),h_hat_21_filt(:,i));
                VAF_true(4,i) = get_VAF(this,h_true_model(:,2,2,i),h_hat_22_filt(:,i));
            end
            
            %% Check Impulse Reponse Function and their fits
            figure;
            for i = floor(linspace(iDex(1),iDex(end),15)) %iDex
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_hat_11_filt(:,i),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model_11(:,i),'-b'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,1,1,i),'-k'); hold on;
            end
            fs = 15;
            view(35,25); title(11); xlim([0 M2/this.sfrq]); %zlim([-0.011 0.011]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta x (m)'); grid on; set(gca,'fontsize',fs);

            i = iDex(10);
            figure;
            plot(this.tvec(1:L),h_hat_11_filt(:,i),'.r','markersize',20); hold on;
            plot(this.tvec(1:L),h_model_11(:,i),'-b','linewidth',2.5); hold on;
            plot(this.tvec(1:L),h_true_model(:,1,1,i),'--k','linewidth',2.5); hold on;
            xlabel('Lag (s)'); ylabel('IRF (m)'); grid on; set(gca,'fontsize',16);
            legend('Estimate','Fit Model','Known Model');

            figure('Position',[88 318 1180 420]); 
            subplot(1,3,1);
            plot(1:R,VAF_output(1,:),'b','linewidth',2.5); hold on;
            ylabel('VAF_{output}'); xlabel('Ensemble Number'); set(gca, 'fontsize', 16); ylim([0 100]);
            xlim([0 300]); 
           
            subplot(1,3,2);
            plot(this.tvec(iDex),VAFirf_11(iDex),'b','linewidth',2.5); hold on;
            ylabel('VAF_{irf}'); xlabel('Time (s)'); set(gca, 'fontsize', 16); ylim([0 100]);
            xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);

            subplot(1,3,3);
            plot(this.tvec(iDex),VAF_true(1,iDex),'b','linewidth',2.5); hold on;
            ylabel('VAF_{true}'); xlabel('Time (s)'); set(gca, 'fontsize', 16); ylim([0 100]);
            xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);
      
            
            % Only Ensemble
            figure('Position',[440 274 560 524]);
            labelss = {'11','12','21','22'};
            
            ax1 = subplot(3,1,1);
            p1=plot(this.tvec(iDex), M_xy_hat(iDex),'-b',this.tvec(iDex), squeeze(M_xy_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            set(gca,'fontsize',16); grid on; xlim([0.2 this.tvec(end)]); ylabel('Mass (kg)'); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylim([0 2.5]);
            
            ax2 = subplot(3,1,2);
            p1=plot(this.tvec(iDex), B_xy_hat(iDex),'-b',this.tvec(iDex), squeeze(B_xy_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            set(gca,'fontsize',16); grid on; xlim([0.2 this.tvec(end)]); ylabel('Damping (Ns/m)'); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); %ylim([0 100]);
            
            ax3 = subplot(3,1,3);
            p1=plot(this.tvec(iDex), squeeze(K_xy_hat(iDex)),'-b',this.tvec(iDex), squeeze(K_xy_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            grid on;  ylabel(' Stiffness (N/m)'); xlabel('Time (s)'); ylim([-1000 4000]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);
            set(gca,'fontsize',16); legend('Estimate','Known');
            
%             % With Least Squares
%             figure('Position',[440 274 560 524]);
%             labelss = {'11','12','21','22'};
%             
%             ax1 = subplot(3,1,1);
%             p1=plot(this.tvec(iDex), M_xy_hat(iDex),'-b',this.tvec(iDex), squeeze(M_xy_const(1,1,iDex)),'--b',this.tvec(iDex),M_par_hat(iDex),':b','linewidth',2.5,'markersize',15); hold on;
%             set(gca,'fontsize',16); grid on; xlim([0.2 this.tvec(end)]); ylabel('Mass (kg)'); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylim([0 2.5]);
%             
%             ax2 = subplot(3,1,2);
%             p1=plot(this.tvec(iDex), B_xy_hat(iDex),'-b',this.tvec(iDex), squeeze(B_xy_const(1,1,iDex)),'--b',this.tvec(iDex),B_par_hat(iDex),':b','linewidth',2.5,'markersize',15); hold on;
%             set(gca,'fontsize',16); grid on; xlim([0.2 this.tvec(end)]); ylabel('Damping (Ns/m)'); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); %ylim([0 100]);
%             
%             ax3 = subplot(3,1,3);
%             p1=plot(this.tvec(iDex), squeeze(K_xy_hat(iDex)),'-b',this.tvec(iDex), squeeze(K_xy_const(1,1,iDex)),'--b',this.tvec(iDex),K_par_hat(iDex),':b','linewidth',2.5,'markersize',15); hold on;
%             grid on;  ylabel(' Stiffness (N/m)'); xlabel('Time (s)'); ylim([-1000 4000]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);
%             set(gca,'fontsize',16); %legend([p1(1), p2(1), p3(1), p4(1)],{'11','12','21','22'});
%             
            disp('test');
                           
            
        end

        function [] = get_2DSim_WAM_MIMO(this)
            
            % To do list:
            % Convert to SISO
            % ??? PLAN
            
            movDuration = 1;
            R = 300;
            
           [y_r_1, y_r_2, u_r_1, u_r_2, f_r_1, f_r_2] = this.get_ensambles_from_simSystem2D(R,movDuration,'line');%,pretScale);            
%            save('/Users/jhermus/Desktop/TmpFigures_Ensamble/WAM_simMINO.mat','y_r_1','y_r_2','u_r_1','u_r_2','f_r_1','f_r_2','this');
%           load('/Users/jhermus/Desktop/TmpFigures_Ensamble/WAM_simMINO.mat');
          
            [z_r_1,z_r_2] = this.get_noise_subtract_x0(y_r_1, y_r_2, u_r_1, u_r_2, f_r_1, f_r_2);
            
            get_Z_spectral(this,z_r_ne11, z_r_ne22, u_r_ne11, u_r_ne22)
            
            %% ID and fit tangential
            [R,N] = size(z_r_1);
            lagBuffer = 5;
            M1 = 0;
            M2 = 0.5*this.sfrq+lagBuffer;
            L = M2-M1+1;
            
            [H_hat_MA, MAwindow] = this.ensambleSysID_Matrix(z_r_1, z_r_2, u_r_1, u_r_2, N, L, R, M1, M2);
            h_hat_11 = H_hat_MA(:,:,1);
            h_hat_12 = H_hat_MA(:,:,2);
            h_hat_21 = H_hat_MA(:,:,3);
            h_hat_22 = H_hat_MA(:,:,4);
            
            
            iDex = M2+1+40:N+M1-40; 
%             iDex = floor(linspace(M2+1+40,N+M1-40,25));

%             [K_hat_par, B_hat_par, M_hat_par] = this.get_parametricLeastSquares(iDex,z_r_1, z_r_2, u_r_1, u_r_2);

            
            %             [h_model_ne11, M_ne11, B_ne11, K_ne11,VAFirf_ne11] = this.fitModel(h_hat_ne11, z_r_ne11, u_r_ne11, N,L_norm, R, M1, M2_norm, iDex);
            %             [h_model_ne12, M_ne12, B_ne12, K_ne12,VAFirf_ne12] = this.fitModel(h_hat_ne12, z_r_ne12, u_r_ne12, N,L_norm, R, M1, M2_norm, iDex);
            %             [h_model_ne21, M_ne21, B_ne21, K_ne21,VAFirf_ne21] = this.fitModel(h_hat_ne21, z_r_ne21, u_r_ne21, N,L_tan, R, M1, M2_tan, iDex);
            %             [h_model_ne22, M_ne22, B_ne22, K_ne22,VAFirf_ne22] = this.fitModel(h_hat_ne22, z_r_ne22, u_r_ne22, N,L_tan, R, M1, M2_tan, iDex);
            
            %% Add filter
            cf = 15; % cutoff freqnency
            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            h_hat_11_filt = zeros(L,N);
            h_hat_12_filt = zeros(L,N);
            h_hat_21_filt = zeros(L,N);
            h_hat_22_filt = zeros(L,N);
            dexNonZero = M2-M1+1+MAwindow/2:N+M1-MAwindow/2;
            for i = dexNonZero
                h_hat_11_filt(:,i) = filtfilt(b,a,h_hat_11(:,i)); % apply fitler
                h_hat_12_filt(:,i) = filtfilt(b,a,h_hat_12(:,i)); % apply fitler
                h_hat_21_filt(:,i) = filtfilt(b,a,h_hat_21(:,i)); % apply fitler
                h_hat_22_filt(:,i) = filtfilt(b,a,h_hat_22(:,i)); % apply fitler
            end
            
            % Fit Model Prony
            M2_short = 0.4*this.sfrq+lagBuffer;
            [K_xy_hat, h_model_11, h_model_12, h_model_21, h_model_22, VAFirf] = get_MIMO_Prony(this,iDex, h_hat_11_filt, h_hat_12_filt, h_hat_21_filt, h_hat_22_filt, M2,M2_short);
            
            %% Estimate imporant values
%             K_ne_hat(1,1,:) =  K_ne_hat(1,1,:)-2500;
%             B_ne11 = B_ne11-40;

%             K_ne_hat = get_filtK(this,K_ne_hat,iDex,10);
                        
            test = WAM_simSystem2D(movDuration);                            
            F_contact = [mean(f_r_1); zeros(1,N)]; % Add tangential if robot dynamics are added
%             figure; plot(F_contact');                            
            s = tf('s');
            
            %% Model least squares fit
%             [K_par_hat, B_par_hat] = get_parametricLeastSquares(this,iDex,z_r_ne1, z_r_ne2, u_r_ne1, u_r_ne2);
            
            %% Make fit plots 
%             figure;
            % Estimate Joint impedance
            for i = iDex

                X = [mean(y_r_1(:,i));mean(y_r_2(:,i))];
                [q_1, q_2] = test.invKino(X(1),X(2));
                J = test.get_jacobian([q_1;q_2]);
                M_q = test.get_massMatrixArm([q_1;q_2]);
                
                [armPosVec_X, armPosVec_Y] = test.plotPosArm(X(1),X(2),'shoulder');
           
%                 ax = plot(armPosVec_X, armPosVec_Y, '-b',...
%                     armPosVec_X, armPosVec_Y, '.b','markersize',30,'linewidth',3);
%                 ylim([0 1.2]); xlim([-0.56 0.56]); grid on; axis equal;
%                 ylabel('Y Distance (m)','fontsize',14);
%                 xlabel('X Distance (m)','fontsize',14);
% %                 legend('Arm','Zero-Force Trajectory','Crank','Constraint Path','Robot','location','northwest'); % Add ZFT
%                 ylabel('Y Distance (m)','fontsize',14);
%                 xlabel('X Distance (m)','fontsize',14);
%                 xlim([-0.6 0.6]);
                
                % Compute kinematic stiffness term
                parJ_parq1(:,:,i) = [-( test.l1a*cos(q_1) + test.l2a*cos(q_1 + q_2) ),  -test.l2a*cos(q_1 + q_2);...
                    -test.l1a*sin(q_1) - test.l2a*sin(q_1 + q_2), -test.l2a*sin(q_1 + q_2)];
                parJ_parq2(:,:,i) = [- test.l2a*cos(q_1 + q_2),  -test.l2a*cos(q_1 + q_2);...
                    - test.l2a*sin(q_1 + q_2), -test.l2a*sin(q_1 + q_2)];
                
                % Map from hand space to joint space using impedance mapping
                if(test.t(i)<=2)
                	K_q_hat(:,:,i) = J'*squeeze(K_xy_hat(:,:,i)- [2500,0;0,2500])*J ;
                else
                    K_q_hat(:,:,i) = J'*squeeze(K_xy_hat(:,:,i))*J;
                end
%                 B_q_hat(:,:,i) = J'*squeeze(B_xy_hat(:,:,i))*J;
%                 M_q_hat(:,:,i) = J'*squeeze(M_xy_hat(:,:,i))*J;
                
                %% Expected Full Rank
                
                % Map joint stiffness to hand space use compliance mapping
                K_xy_const(:,:,i) = inv(J*inv(test.K)*J');% - parJ_parq1(:,:,i)*F_contact(1,i) + parJ_parq2(:,:,i)*F_contact(2,i);
                B_xy_const(:,:,i) = inv(J*inv(test.B)*J');
                M_xy_const(:,:,i) = inv(J*inv(M_q)*J');
                
                if(test.t(i)<=2)
                  K_xy_const(:,:,i) = K_xy_const(:,:,i) + [2500,0;0,2500];
                  B_xy_const(:,:,i) = B_xy_const(:,:,i) + [40,0;0,40];
                end
                
                K_q_const(:,:,i) = test.K;
                B_q_const(:,:,i) = test.B;
                M_q_const(:,:,i) = M_q;

                Z_true_model = M_xy_const(:,:,i)*s^2 + [B_xy_const(:,:,i)]*s + [K_xy_const(:,:,i)]; % Includes robot
           
%                 figure; impulse(inv(Z_true_model));
                Y_true_model = inv(Z_true_model);
                h_true_model(:,:,:,i) = impulse(Y_true_model,this.tvec(1:L));
                
                pause(0.01);
                
            end
            
            
            z_r_hat1 = zeros(size(z_r_1));
            z_r_hat2 = zeros(size(z_r_2));
            
            for r = 1:R
                for i = iDex
                    clear tmp1 tmp2
                    for j = 1:M2
                        
                        tmp1(j) = this.dt*h_hat_11_filt(j,i,:)*u_r_1(r,i-j) +...
                                  this.dt*h_hat_12_filt(j,i,:)*u_r_2(r,i-j);
                              
                        tmp2(j) = this.dt*h_hat_12_filt(j,i,:)*u_r_1(r,i-j) +...
                                  this.dt*h_hat_22_filt(j,i,:)*u_r_2(r,i-j);
                    end
                    
                    z_r_hat1(r,i) = sum(tmp1);
                    z_r_hat2(r,i) = sum(tmp2);
                    
                end
            end
            
             for r = 1:R
                VAF_output(1,r) = get_VAF(this,z_r_1(r,iDex(1):iDex(end)),z_r_hat1(r,iDex(1):iDex(end)));
                VAF_output(2,r) = get_VAF(this,z_r_2(r,iDex(1):iDex(end)),z_r_hat2(r,iDex(1):iDex(end)));
             end
                        
            for i = iDex
                VAF_true(1,i) = get_VAF(this,h_true_model(:,1,1,i),h_hat_11_filt(:,i));
                VAF_true(2,i) = get_VAF(this,h_true_model(:,1,2,i),h_hat_12_filt(:,i));
                VAF_true(3,i) = get_VAF(this,h_true_model(:,2,1,i),h_hat_21_filt(:,i));
                VAF_true(4,i) = get_VAF(this,h_true_model(:,2,2,i),h_hat_22_filt(:,i));
            end
            
            %% Check Impulse Reponse Function and their fits
            figure('position',[289 231 895 527]);
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3);
            ax4 = subplot(2,2,4);
            
            for i = floor(linspace(iDex(1),iDex(end),15)) %iDex
                axes(ax1);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_hat_11_filt(:,i),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model_11(:,i),'-b'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,1,1,i),'-k'); hold on;
                
                axes(ax2);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_hat_12_filt(:,i),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model_12(:,i),'-b'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,1,2,i),'-k'); hold on;

                axes(ax3);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_hat_21_filt(:,i),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model_21(:,i),'-b'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,2,1,i),'-k'); hold on;

                axes(ax4);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_hat_22_filt(:,i),'.r'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model_22(:,i),'-b'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,2,2,i),'-k'); hold on;

            end
            fs = 15;
            axes(ax1); view(35,25); title(11); xlim([0 M2/this.sfrq]); %zlim([-0.011 0.011]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta x (m)'); grid on; set(gca,'fontsize',fs);
            
            axes(ax2); view(35,25); title(12); xlim([0 M2/this.sfrq]); %zlim([-0.011 0.011]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta y (m)'); grid on; set(gca,'fontsize',fs);
            
            axes(ax3); view(35,25); title(21); xlim([0 M2/this.sfrq]); %zlim([-0.011 0.011]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta x (m)'); grid on; set(gca,'fontsize',fs);
            
            axes(ax4); view(35,25); title(22); xlim([0 M2/this.sfrq]); %zlim([-0.1 0.1]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta y (m)'); grid on; set(gca,'fontsize',fs);
            
            i = iDex(10);
            figure;
            plot(this.tvec(1:L),h_hat_21_filt(:,i),'.r','markersize',20); hold on;
            plot(this.tvec(1:L),h_model_21(:,i),'-b','linewidth',2.5); hold on;
            plot(this.tvec(1:L),h_true_model(:,2,1,i),'--k','linewidth',2.5); hold on;
            xlabel('Lag (s)'); ylabel('IRF (m)'); grid on; set(gca,'fontsize',16);
            legend('Estimate','Fit Model','Known Model');

            
            figure('Position',[88 318 1180 420]); 
            subplot(1,3,1);
            plot(1:R,VAF_output(1,:),'b','linewidth',2.5); hold on;
            plot(1:R,VAF_output(2,:),'r','linewidth',2.5);
            ylabel('VAF_{output}'); xlabel('Ensemble Number'); set(gca, 'fontsize', 16); ylim([0 100]);
            xlim([0 300]); legend('Normal','Tangential','location','southwest');
           
            subplot(1,3,2);
            plot(this.tvec(iDex),VAFirf(1,iDex),'b','linewidth',2.5); hold on;
            plot(this.tvec(iDex),VAFirf(2,iDex),'r','linewidth',2.5);
            plot(this.tvec(iDex),VAFirf(3,iDex),'g','linewidth',2.5);
            plot(this.tvec(iDex),VAFirf(4,iDex),'k','linewidth',2.5);
            legend('11','12','21','22','location','southwest');
            ylabel('VAF_{irf}'); xlabel('Time (s)'); set(gca, 'fontsize', 16); ylim([0 100]);
            xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);

            subplot(1,3,3);
            plot(this.tvec(iDex),VAF_true(1,iDex),'b','linewidth',2.5); hold on;
            plot(this.tvec(iDex),VAF_true(2,iDex),'r','linewidth',2.5);
            plot(this.tvec(iDex),VAF_true(3,iDex),'g','linewidth',2.5);
            plot(this.tvec(iDex),VAF_true(4,iDex),'k','linewidth',2.5);
            legend('11','12','21','22','location','southwest');
            ylabel('VAF_{true}'); xlabel('Time (s)'); set(gca, 'fontsize', 16); ylim([0 100]);
            xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);
      
            
            figure('Position',[440 274 560 524]);
            labelss = {'11','12','21','22'};
            
%             ax1 = subplot(3,1,1);
%             p1=plot(this.tvec(iDex), M_ne11(iDex),'.-b',this.tvec(iDex), squeeze(M_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
%             p2=plot(this.tvec(iDex), M_ne12(iDex),'.-r',this.tvec(iDex), squeeze(M_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
%             p3=plot(this.tvec(iDex), M_ne21(iDex),'.-g',this.tvec(iDex), squeeze(M_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
%             p4=plot(this.tvec(iDex), M_ne22(iDex),'.-k',this.tvec(iDex), squeeze(M_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
%             grid on;  xlim([0.2 this.tvec(end)]); ylabel('Mass (kg)'); ylim([0 5]);
%             legend([p1(1),p2(1),p3(1),p4(1)],labelss,'location','northwest');legend boxoff;
            
%             ax2 = subplot(3,1,2);
%             p1=plot(this.tvec(iDex), B_ne11(iDex),'.-b',this.tvec(iDex), squeeze(B_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
%             p2=plot(this.tvec(iDex), B_ne12(iDex),'.-r',this.tvec(iDex), squeeze(B_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
%             p3=plot(this.tvec(iDex), B_ne21(iDex),'.-g',this.tvec(iDex), squeeze(B_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
%             p4=plot(this.tvec(iDex), B_ne22(iDex),'.-k',this.tvec(iDex), squeeze(B_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
%             grid on;  xlim([0.2 this.tvec(end)]); ylabel('Damping (Ns/m)'); %ylim([0 100]);
            
%             ax3 = subplot(3,1,3);
            p1=plot(this.tvec(iDex), squeeze(K_xy_hat(1,1,iDex)),'-b',this.tvec(iDex), squeeze(K_xy_const(1,1,iDex)),'--b',this.tvec(iDex),K_hat_par(iDex,1),':b','linewidth',2.5,'markersize',15); hold on;
            p2=plot(this.tvec(iDex), squeeze(K_xy_hat(1,2,iDex)),'-r',this.tvec(iDex), squeeze(K_xy_const(1,2,iDex)),'--r',this.tvec(iDex),K_hat_par(iDex,2),':r','linewidth',2.5,'markersize',15); hold on;
            p3=plot(this.tvec(iDex), squeeze(K_xy_hat(2,1,iDex)),'-g',this.tvec(iDex), squeeze(K_xy_const(2,1,iDex)),'--g',this.tvec(iDex),K_hat_par(iDex,2),':g','linewidth',2.5,'markersize',15); hold on;
            p4=plot(this.tvec(iDex), squeeze(K_xy_hat(2,2,iDex)),'-k',this.tvec(iDex), squeeze(K_xy_const(2,2,iDex)),'--k',this.tvec(iDex),K_hat_par(iDex,3),':k','linewidth',2.5,'markersize',15);
            grid on;  ylabel(' Stiffness (N/m)'); xlabel('Time (s)'); ylim([-1000 4000]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);
            set(gca,'fontsize',16);legend([p1(1), p2(1), p3(1), p4(1)],{'11','12','21','22'});
            
%              p1=plot(this.tvec(iDex), squeeze(K_ne_hat(1,1,iDex)),'-b',this.tvec(iDex), squeeze(K_ne_const(1,1,iDex)),'--b',this.tvec(iDex), K_par_hat(iDex,1),':b','linewidth',2.5,'markersize',15); hold on;
%             p2=plot(this.tvec(iDex), squeeze(K_ne_hat(1,2,iDex)),'-r',this.tvec(iDex), squeeze(K_ne_const(1,2,iDex)),'--r',this.tvec(iDex), K_par_hat(iDex,2),':r','linewidth',2.5,'markersize',15); hold on;
%             p3=plot(this.tvec(iDex), squeeze(K_ne_hat(2,1,iDex)),'-g',this.tvec(iDex), squeeze(K_ne_const(2,1,iDex)),'--g',this.tvec(iDex), K_par_hat(iDex,3),':g','linewidth',2.5,'markersize',15); hold on;
%             p4=plot(this.tvec(iDex), squeeze(K_ne_hat(2,2,iDex)),'-k',this.tvec(iDex), squeeze(K_ne_const(2,2,iDex)),'--k',this.tvec(iDex), K_par_hat(iDex,4),':k','linewidth',2.5,'markersize',15);
%             grid on;  ylabel(' Stiffness (N/m)'); xlabel('Time (s)'); ylim([-2000 7000]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]);
%             set(gca,'fontsize',16);
            
            

%             figure;
%             plot(this.tvec(iDex), VAFirf_ne11(iDex),'.-b','linewidth',2.5,'markersize',15); hold on;
%             plot(this.tvec(iDex), VAFirf_ne12(iDex),'.-r','linewidth',2.5,'markersize',15); hold on;
%             plot(this.tvec(iDex), VAFirf_ne21(iDex),'.-g','linewidth',2.5,'markersize',15); hold on;
%             plot(this.tvec(iDex), VAFirf_ne22(iDex),'.-k','linewidth',2.5,'markersize',15);
%             grid on;  xlim([0.2 this.tvec(end)]); ylabel(' VAF_{irf}'); ylim([0 100]);
%             xlabel('Time (s)'); set(gca,'fontsize',18);
    
           figure;
%             subplot(3,1,1); plot(this.tvec(iDex), squeeze(M_q_hat(1,1,iDex)),'.-b',this.tvec(iDex),squeeze(M_q_const(1,1,iDex)),'--b',...
%                 this.tvec(iDex), squeeze(M_q_hat(1,2,iDex)),'.-r',this.tvec(iDex),squeeze(M_q_const(1,2,iDex)),'--r',...
%                 this.tvec(iDex), squeeze(M_q_hat(2,1,iDex)),'.-g',this.tvec(iDex),squeeze(M_q_const(2,1,iDex)),'--g',...
%                 this.tvec(iDex), squeeze(M_q_hat(2,2,iDex)),'.-k',this.tvec(iDex),squeeze(M_q_const(2,2,iDex)),'--k','markersize',20,'linewidth',2.5);
%             title('Joint Space'); ylim([0 1]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylabel('Inertia (Kg m^2)'); xlabel('Time (s)'); set(gca,'fontsize',18);
            
%             subplot(3,1,2); plot(this.tvec(iDex), squeeze(B_q_hat(1,1,iDex)),'.-b',this.tvec(iDex),squeeze(B_q_const(1,1,iDex)),'--b',...
%                 this.tvec(iDex), squeeze(B_q_hat(1,2,iDex)),'.-r',this.tvec(iDex),squeeze(B_q_const(1,2,iDex)),'--r',...
%                 this.tvec(iDex), squeeze(B_q_hat(2,1,iDex)),'.-g',this.tvec(iDex),squeeze(B_q_const(2,1,iDex)),'--g',...
%                 this.tvec(iDex), squeeze(B_q_hat(2,2,iDex)),'.-k',this.tvec(iDex),squeeze(B_q_const(2,2,iDex)),'--k','markersize',20,'linewidth',2.5);
%             ylim([0 10]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylabel('Damping (Nms/rad)'); xlabel('Time (s)'); set(gca,'fontsize',18);
            
%             subplot(3,1,3); 
            ax1 = plot(this.tvec(iDex), squeeze(K_q_hat(1,1,iDex)),'-b',this.tvec(iDex),squeeze(K_q_const(1,1,iDex)),'--b','markersize',20,'linewidth',2.5); hold on;
            ax2 = plot(this.tvec(iDex), squeeze(K_q_hat(1,2,iDex)),'-r',this.tvec(iDex),squeeze(K_q_const(1,2,iDex)),'--r','markersize',20,'linewidth',2.5); hold on;
            ax3 = plot(this.tvec(iDex), squeeze(K_q_hat(2,1,iDex)),'-g',this.tvec(iDex),squeeze(K_q_const(2,1,iDex)),'--g','markersize',20,'linewidth',2.5); hold on;
            ax4 = plot(this.tvec(iDex), squeeze(K_q_hat(2,2,iDex)),'-k',this.tvec(iDex),squeeze(K_q_const(2,2,iDex)),'--k','markersize',20,'linewidth',2.5);
            ylim([0 600]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylabel('Joint Stiffness (Nm/rad)'); xlabel('Time (s)'); set(gca,'fontsize',18);
            legend([ax1(1), ax2(1), ax3(1), ax4(1)], {'11','12','21','22'});
            
            % Fit optimization based
%             [opt_h_model_ne11, opt_m_ne11, opt_b_ne11, opt_k_ne11,opt_VAFirf_ne11] = this.fitModel(h_hat_ne11_filt, N, L, M1, M2, iDex);
%             
%             figure; 
%             plot(opt_k_ne11-2500); hold on;
%             plot(squeeze(K_ne_const(1,1,:)));
            
            disp('test');
                           
            
                end
   
        function [] = get_2DSim_circle_MIMO_spectral(this)
            
            turnPeriod = 13.33;
            numLoop = 45;
            
%             [z_r_ne11,u_r_ne11,y_r_ne11,z_r_ne22,u_r_ne22,y_r_ne22,t_cycle,N,R] = this.get_ensambles_from_simSystem2D(numLoop,turnPeriod,'circle');
%             % %
%             save('/Users/jhermus/Desktop/TmpFigures_Ensamble/SimSaved/circle_MIMO.mat','z_r_ne11','u_r_ne11','y_r_ne11',...
%                 'z_r_ne22','u_r_ne22','y_r_ne22','t_cycle','N','R','this');
            load('/Users/jhermus/Desktop/TmpFigures_Ensamble/SimSaved/circle_MIMO.mat');
            
            this.tvec = t_cycle;
            
            %% ID and fit tangential
            lagBuffer = 5;
            M1 = 0;
            M2_norm = 0.7*this.sfrq+lagBuffer;
            M2_tan = 0.7*this.sfrq+lagBuffer;
            M2 = M2_tan;
            L_norm = M2_norm-M1+1;
            L_tan = M2_tan-M1+1;
            L = L_tan;
            
            %% Try spectral appraoch
%           get_Z_spectral(this,z_r_ne11, z_r_ne22, u_r_ne11, u_r_ne22);
            
            binDex = floor(1:size(z_r_ne11,2)/6:size(z_r_ne11,2));
            for i = 1:length(binDex)
                [TF_freq,Z11_l_mag(:,i),Z12_l_mag(:,i),Z21_l_mag(:,i),Z22_l_mag(:,i),Z11_l_phi(:,i),Z12_l_phi(:,i),Z21_l_phi(:,i),Z22_l_phi(:,i)] = get_Z_spectral_individual(this, binDex(i), z_r_ne11, z_r_ne22, u_r_ne11, u_r_ne22);
            end
            
            
            figure('position',[289 231 895 527]);
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3);
            ax4 = subplot(2,2,4);
            
            for i = 1:length(binDex)
                axes(ax1); plot3(TF_freq,this.tvec(binDex(i))*ones(1,length(TF_freq)),Z11_l_mag(:,i),'.r'); hold on;
%                 axes(ax2); plot3(TF_freq,this.tvec(binDex(i))*ones(1,length(TF_freq)),Z11_l_mag(:,i),'.r'); hold on;
%                 axes(ax3); plot3(TF_freq,this.tvec(binDex(i))*ones(1,length(TF_freq)),Z22_l_mag(:,i),'.r'); hold on;
                axes(ax4); plot3(TF_freq,this.tvec(binDex(i))*ones(1,length(TF_freq)),Z22_l_mag(:,i),'.r'); hold on;
            end
            fs = 15;
            axes(ax1); set(ax1, 'ZScale', 'log','XScale','log'); xlim([0 M2_norm/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);
            axes(ax4); set(ax4, 'ZScale', 'log','XScale','log'); xlim([0 M2_tan/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);

            %view(35,25); title(11); xlim([0 M2_norm/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);
%             axes(ax2); view(35,25); title(12); xlim([0 M2_norm/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_t (m)'); grid on; set(gca,'fontsize',fs);
%             axes(ax3); view(35,25); title(21); xlim([0 M2_tan/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);
%             axes(ax4); view(35,25); title(22); xlim([0 M2_tan/this.sfrq]); zlim([-0.1 0.1]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_t (m)'); grid on; set(gca,'fontsize',fs);
            
            
            
            
            %% Continue ensamble method
            
            [H_hat_MA, MAwindow] = this.ensambleSysID_Matrix(z_r_ne11, z_r_ne22, u_r_ne11, u_r_ne22, N, L_tan, R, M1, M2_tan);
            h_hat_ne11 = H_hat_MA(:,:,1);
            h_hat_ne12 = H_hat_MA(:,:,2);
            h_hat_ne21 = H_hat_MA(:,:,3);
            h_hat_ne22 = H_hat_MA(:,:,4);
            
            %             [h_hat_ne11, MAwindow] = this.ensambleSysID(z_r_ne11, u_r_ne11, N, L_norm, R, M1, M2_norm);
            %             [h_hat_ne12, MAwindow] = this.ensambleSysID(z_r_ne12, u_r_ne12, N, L_norm, R, M1, M2_norm);
            %             [h_hat_ne21, MAwindow] = this.ensambleSysID(z_r_ne21, u_r_ne21, N, L_tan, R, M1, M2_tan);
            %             [h_hat_ne22, MAwindow] = this.ensambleSysID(z_r_ne22, u_r_ne22, N, L_tan, R, M1, M2_tan);
            
            M2_start = max([M2_norm, M2_tan]);
            iDex = floor(linspace(M2_start+1+40,N+M1-3,25));
            
            %             [h_model_ne11, M_ne11, B_ne11, K_ne11,VAFirf_ne11] = this.fitModel(h_hat_ne11, z_r_ne11, u_r_ne11, N,L_norm, R, M1, M2_norm, iDex);
            %             [h_model_ne12, M_ne12, B_ne12, K_ne12,VAFirf_ne12] = this.fitModel(h_hat_ne12, z_r_ne12, u_r_ne12, N,L_norm, R, M1, M2_norm, iDex);
            %             [h_model_ne21, M_ne21, B_ne21, K_ne21,VAFirf_ne21] = this.fitModel(h_hat_ne21, z_r_ne21, u_r_ne21, N,L_tan, R, M1, M2_tan, iDex);
            %             [h_model_ne22, M_ne22, B_ne22, K_ne22,VAFirf_ne22] = this.fitModel(h_hat_ne22, z_r_ne22, u_r_ne22, N,L_tan, R, M1, M2_tan, iDex);
            
            %% Add filter
            
            cf = 30; % cutoff freqnency
            [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            h_hat_ne11_filt = zeros(L,N);
            h_hat_ne12_filt = zeros(L,N);
            h_hat_ne21_filt = zeros(L,N);
            h_hat_ne22_filt = zeros(L,N);
            dexNonZero = M2-M1+1+MAwindow/2:N+M1-MAwindow/2;
            for i = dexNonZero
                h_hat_ne11_filt(:,i) = filtfilt(b,a,h_hat_ne11(:,i)); % apply fitler
                h_hat_ne12_filt(:,i) = filtfilt(b,a,h_hat_ne12(:,i)); % apply fitler
                h_hat_ne21_filt(:,i) = filtfilt(b,a,h_hat_ne21(:,i)); % apply fitler
                h_hat_ne22_filt(:,i) = filtfilt(b,a,h_hat_ne22(:,i)); % apply fitler
            end
            
            %% Check Impulse Reponse Function and their fits
            figure('position',[289 231 895 527]);
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3);
            ax4 = subplot(2,2,4);
            
            for i = iDex
                axes(ax1);
                plot3(this.tvec(1:L_norm),this.tvec(i)*ones(1,L_norm),h_hat_ne11_filt(:,i),'.r'); hold on;
                %                 plot3(this.tvec(1:L_norm),this.tvec(i)*ones(1,L_norm),h_model_ne11(:,i),'-b'); hold on;
                %
                axes(ax2);
                plot3(this.tvec(1:L_norm),this.tvec(i)*ones(1,L_norm),h_hat_ne12_filt(:,i),'.r'); hold on;
                %                 plot3(this.tvec(1:L_norm),this.tvec(i)*ones(1,L_norm),h_model_ne12(:,i),'-b'); hold on;
                
                axes(ax3);
                plot3(this.tvec(1:L_tan),this.tvec(i)*ones(1,L_tan),h_hat_ne21_filt(:,i),'.r'); hold on;
                %                 plot3(this.tvec(1:L_tan),this.tvec(i)*ones(1,L_tan),h_model_ne21(:,i),'-b'); hold on;
                
                axes(ax4);
                plot3(this.tvec(1:L_tan),this.tvec(i)*ones(1,L_tan),h_hat_ne22_filt(:,i),'.r'); hold on;
                %                 plot3(this.tvec(1:L_tan),this.tvec(i)*ones(1,L_tan),h_model_ne22(:,i),'-b'); hold on;
            end
            fs = 15;
            axes(ax1); view(35,25); title(11); xlim([0 M2_norm/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);
            axes(ax2); view(35,25); title(12); xlim([0 M2_norm/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_t (m)'); grid on; set(gca,'fontsize',fs);
            axes(ax3); view(35,25); title(21); xlim([0 M2_tan/this.sfrq]); zlim([-0.011 0.011]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);
            axes(ax4); view(35,25); title(22); xlim([0 M2_tan/this.sfrq]); zlim([-0.1 0.1]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_t (m)'); grid on; set(gca,'fontsize',fs);
            
            
            %% Prony
            figure;
            
            i = iDex(8);
            [b,a] = prony(h_hat_ne11_filt(:,i),4,4);
            C_hat_d11 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c11 = d2c(C_hat_d11,'tustin');
            H_hat_c11 = impulse(C_hat_c11,0:1/this.sfrq:M2/this.sfrq);
            
            subplot(2,2,1); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne11_filt(:,i),0:1/this.sfrq:M2/this.sfrq,H_hat_c11);
            
            [b,a] = prony(h_hat_ne12_filt(:,i),4,4);
            C_hat_d12 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c12 = d2c(C_hat_d12,'tustin');
            H_hat_c12 = impulse(C_hat_c12,0:1/this.sfrq:M2/this.sfrq);
            
            subplot(2,2,2); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne12_filt(:,i),0:1/this.sfrq:M2/this.sfrq,H_hat_c12);
            
            %             [b,a] = prony(h_hat_ne21_filt(:,i),4,4);
            %             C_hat_d21 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            %             C_hat_c21 = d2c(C_hat_d21,'tustin');
            %             H_hat_c21 = impulse(C_hat_c21,0:1/this.sfrq:M2/this.sfrq);
            %
            %            subplot(2,2,3);plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne21_filt(:,i),0:1/this.sfrq:M2/this.sfrq,H_hat_c21);
            
            [b,a] = prony(h_hat_ne22_filt(:,i),4,4);
            C_hat_d22 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c22 = d2c(C_hat_d22,'tustin');
            H_hat_c22 = impulse(C_hat_c22,0:1/this.sfrq:M2/this.sfrq);
            
            subplot(2,2,4); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne22_filt(:,i),0:1/this.sfrq:M2/this.sfrq,H_hat_c22);
            
            
            Z_hat = inv([C_hat_c11,C_hat_c12;C_hat_c12,C_hat_c22]);
            
            
            figure; pzplot(Z_hat);
            figure;
            bode(Z_hat); hold on;
            bode(Z); xlim([10^-1 10^2]);
            
            %% Make fit plots %% FIX When contraint is back on.
            K_ne11 = K_ne11-2500;
            K_ne12 = K_ne12-2500;
            B_ne11 = B_ne11-40;
            B_ne12 = B_ne12-40;
            
            iDex = iDex(1:end-1);
            
            test = simSystem2D(turnPeriod);
            thcp = mean(y_r_ne21)./test.lc - pi;
            
            % figure;
            % Estimate Joint impedance
            for i = iDex
                
                X = [test.d1 - test.lc*cos(thcp(i))'; test.d2 - test.lc*sin(thcp(i))'];
                [tmp,ee,nn] = test.get_crankPostion(X,'shoulder');
                [q_1, q_2] = test.invKino(X(1),X(2));
                J = test.get_jacobian([q_1;q_2]);
                M_q = test.get_massMatrixArm([q_1;q_2]);
                
                K_ne_hat(:,:,i) = [K_ne11(i), 0; 0, K_ne22(i)];
                B_ne_hat(:,:,i) = [B_ne11(i), 0; 0, B_ne22(i)];
                M_ne_hat(:,:,i) = [M_ne11(i), 0; 0, M_ne22(i)];
                
                % Transform from nt to xy
                K_xy_hat(:,:,i) = inv([nn,ee]')*squeeze(K_ne_hat(:,:,i))*inv([nn,ee]); % Neglect kinematic stiffness for now
                B_xy_hat(:,:,i) = inv([nn,ee]')*squeeze(B_ne_hat(:,:,i))*inv([nn,ee]);
                M_xy_hat(:,:,i) = inv([nn,ee]')*squeeze(M_ne_hat(:,:,i))*inv([nn,ee]);
                
                % Map from hand space to joint space using impedance mapping
                K_q_hat(:,:,i) = J'*squeeze(K_xy_hat(:,:,i))*J;
                B_q_hat(:,:,i) = J'*squeeze(B_xy_hat(:,:,i))*J;
                M_q_hat(:,:,i) = J'*squeeze(M_xy_hat(:,:,i))*J;
                
                %% Expected Full Rank
                
                % Map joint stiffness to hand space use compliance mapping
                K_xy_const(:,:,i) = inv(J*inv(test.K)*J');
                B_xy_const(:,:,i) = inv(J*inv(test.B)*J');
                M_xy_const(:,:,i) = inv(J*inv(M_q)*J');
                
                % Transform from xy to nt
                K_ne_const(:,:,i) = [nn,ee]'*squeeze(K_xy_const(:,:,i))*[nn,ee];
                B_ne_const(:,:,i) = [nn,ee]'*squeeze(B_xy_const(:,:,i))*[nn,ee];
                M_ne_const(:,:,i) = [nn,ee]'*squeeze(M_xy_const(:,:,i))*[nn,ee];
                
                K_q_const(:,:,i) = test.K;
                B_q_const(:,:,i) = test.B;
                M_q_const(:,:,i) = M_q;
                
                % Compute kinematic stiffness term
                %                parJ_parq1(:,:,i) = [-( test.l1a*cos(q_1) + test.l2a*cos(q_1 + q_2) ),  -test.l2a*cos(q_1 + q_2);...
                %                    -test.l1a*sin(q_1) - test.l2a*sin(q_1 + q_2), -test.l2a*sin(q_1 + q_2)];
                %                parJ_parq2(:,:,i) = [- test.l2a*cos(q_1 + q_2),  -test.l2a*cos(q_1 + q_2);...
                %                    - test.l2a*sin(q_1 + q_2), -test.l2a*sin(q_1 + q_2)];
                
                %                plot(X(1),X(2),'o'); hold on;
                %                pause();
            end
            %             axis equal; grid on;
            
            
            figure('Position',[440 274 560 524]);
            labelss = {'11','12','21','22'};
            
            ax1 = subplot(3,1,1);
            p1=plot(this.tvec(iDex), M_ne11(iDex),'.-b',this.tvec(iDex), squeeze(M_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            p2=plot(this.tvec(iDex), M_ne12(iDex),'.-r',this.tvec(iDex), squeeze(M_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
            p3=plot(this.tvec(iDex), M_ne21(iDex),'.-g',this.tvec(iDex), squeeze(M_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
            p4=plot(this.tvec(iDex), M_ne22(iDex),'.-k',this.tvec(iDex), squeeze(M_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel('Mass (kg)'); ylim([0 5]);
            legend([p1(1),p2(1),p3(1),p4(1)],labelss,'location','northwest');legend boxoff;
            
            ax2 = subplot(3,1,2);
            p1=plot(this.tvec(iDex), B_ne11(iDex),'.-b',this.tvec(iDex), squeeze(B_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            p2=plot(this.tvec(iDex), B_ne12(iDex),'.-r',this.tvec(iDex), squeeze(B_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
            p3=plot(this.tvec(iDex), B_ne21(iDex),'.-g',this.tvec(iDex), squeeze(B_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
            p4=plot(this.tvec(iDex), B_ne22(iDex),'.-k',this.tvec(iDex), squeeze(B_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel('Damping (Ns/m)'); %ylim([0 100]);
            
            ax3 = subplot(3,1,3);
            p1=plot(this.tvec(iDex), K_ne11(iDex),'.-b',this.tvec(iDex), squeeze(K_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            p2=plot(this.tvec(iDex), K_ne12(iDex),'.-r',this.tvec(iDex), squeeze(K_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
            p3=plot(this.tvec(iDex), K_ne21(iDex),'.-g',this.tvec(iDex), squeeze(K_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
            p4=plot(this.tvec(iDex), K_ne22(iDex),'.-k',this.tvec(iDex), squeeze(K_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel(' Stiffness (N/m)'); ylim([0 3000]);
            
            figure;
            plot(this.tvec(iDex), VAFirf_ne11(iDex),'.-b','linewidth',2.5,'markersize',15); hold on;
            plot(this.tvec(iDex), VAFirf_ne12(iDex),'.-r','linewidth',2.5,'markersize',15); hold on;
            plot(this.tvec(iDex), VAFirf_ne21(iDex),'.-g','linewidth',2.5,'markersize',15); hold on;
            plot(this.tvec(iDex), VAFirf_ne22(iDex),'.-k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel(' VAF_{irf}'); ylim([0 100]);
            xlabel('Time (s)'); set(gca,'fontsize',18);
            
            figure;
            subplot(3,1,1); plot(this.tvec(iDex), squeeze(M_q_hat(1,1,iDex)),'.-b',this.tvec(iDex),squeeze(M_q_const(1,1,iDex)),'--b',...
                this.tvec(iDex), squeeze(M_q_hat(1,2,iDex)),'.-r',this.tvec(iDex),squeeze(M_q_const(1,2,iDex)),'--r',...
                this.tvec(iDex), squeeze(M_q_hat(2,1,iDex)),'.-g',this.tvec(iDex),squeeze(M_q_const(2,1,iDex)),'--g',...
                this.tvec(iDex), squeeze(M_q_hat(2,2,iDex)),'.-k',this.tvec(iDex),squeeze(M_q_const(2,2,iDex)),'--k','markersize',20,'linewidth',2.5);
            title('Joint Space'); ylim([0 1]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylabel('Inertia (Kg m^2)'); xlabel('Time (s)'); set(gca,'fontsize',18);
            
            subplot(3,1,2); plot(this.tvec(iDex), squeeze(B_q_hat(1,1,iDex)),'.-b',this.tvec(iDex),squeeze(B_q_const(1,1,iDex)),'--b',...
                this.tvec(iDex), squeeze(B_q_hat(1,2,iDex)),'.-r',this.tvec(iDex),squeeze(B_q_const(1,2,iDex)),'--r',...
                this.tvec(iDex), squeeze(B_q_hat(2,1,iDex)),'.-g',this.tvec(iDex),squeeze(B_q_const(2,1,iDex)),'--g',...
                this.tvec(iDex), squeeze(B_q_hat(2,2,iDex)),'.-k',this.tvec(iDex),squeeze(B_q_const(2,2,iDex)),'--k','markersize',20,'linewidth',2.5);
            ylim([0 10]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylabel('Damping (Nms/rad)'); xlabel('Time (s)'); set(gca,'fontsize',18);
            
            subplot(3,1,3); plot(this.tvec(iDex), squeeze(K_q_hat(1,1,iDex)),'.-b',this.tvec(iDex),squeeze(K_q_const(1,1,iDex)),'--b',...
                this.tvec(iDex), squeeze(K_q_hat(1,2,iDex)),'.-r',this.tvec(iDex),squeeze(K_q_const(1,2,iDex)),'--r',...
                this.tvec(iDex), squeeze(K_q_hat(2,1,iDex)),'.-g',this.tvec(iDex),squeeze(K_q_const(2,1,iDex)),'--g',...
                this.tvec(iDex), squeeze(K_q_hat(2,2,iDex)),'.-k',this.tvec(iDex),squeeze(K_q_const(2,2,iDex)),'--k','markersize',20,'linewidth',2.5);
            ylim([0 200]); xlim([this.tvec(iDex(1)) this.tvec(iDex(end))]); ylabel('Stiffness (Nm/rad)'); xlabel('Time (s)'); set(gca,'fontsize',18);
            
            disp('test');
            
            % Figure with only the expected hand space
            %              figure('Position',[440 274 560 524]);
            %              labelss = {'11','12','21','22'};
            %
            %              ax1 = subplot(3,1,1);
            %              p1=plot(squeeze(M_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            %              p2=plot(this.tvec(iDex), squeeze(M_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
            %              p3=plot(this.tvec(iDex), squeeze(M_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
            %              p4=plot(this.tvec(iDex), squeeze(M_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
            %              grid on;  xlim([0.2 this.tvec(end)]); ylabel('Mass (kg)'); ylim([-1 3]);
            %
            %              ax2 = subplot(3,1,2);
            %              p1=plot(squeeze(B_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            %              p2=plot(this.tvec(iDex), squeeze(B_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
            %              p3=plot(this.tvec(iDex), squeeze(B_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
            %              p4=plot(this.tvec(iDex), squeeze(B_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
            %              grid on;  xlim([0.2 this.tvec(end)]); ylabel('Damping (Ns/m)'); ylim([-50 50]);
            %
            %              ax3 = subplot(3,1,3);
            %              p1=plot(this.tvec(iDex), squeeze(K_ne_const(1,1,iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            %              p2=plot(this.tvec(iDex), squeeze(K_ne_const(1,2,iDex)),'--r','linewidth',2.5,'markersize',15); hold on;
            %              p3=plot(this.tvec(iDex), squeeze(K_ne_const(2,1,iDex)),'--g','linewidth',2.5,'markersize',15); hold on;
            %              p4=plot(this.tvec(iDex), squeeze(K_ne_const(2,2,iDex)),'--k','linewidth',2.5,'markersize',15);
            %              grid on;  xlim([0.2 this.tvec(end)]); xlabel('Time (s)'); ylabel(' Stiffness (N/m)'); ylim([-1000 3000]);
            %              legend([p1(1),p2(1),p3(1),p4(1)],labelss,'location','northwest');legend boxoff;
            %              set(ax1,'fontsize',18);set(ax2,'fontsize',18);set(ax3,'fontsize',18);
            
        end
        
        function [] = get_WAMConstant(this)
            
            pathh = '/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/proficio_working/matlab_Analysis/dataFile_notTracked/';
            %           filee = '20210115aft00.csv';
            %           filee = '20210115aft01.csv';
            %           filee = '20210115aft03.csv';
            %           filee = '20210115aft04.csv';
            %           filee = '20210115aft05.csv'; % Keep normalization, it-max = 8, Kx =1000, Bx = 20
            %           filee = '20210115aft06.csv';
            %           filee = '20210202aft00.csv';
            
            %           filee = '20210203aft00.csv'; % | K_x = 1000, B_x = 20, it_Max = 8, f_pret[1]=0; 2N amplitude |
            %           filee = '20210203aft01.csv'; % | K_x = 1000, B_x = 20, it_Max = 8, f_pret[0]=0; 2N amplitude |
            %           filee = '20210203aft02.csv'; % | K_x = 1000, B_x = 20, it_Max = 8, f_pret[1]=0; 5N amplitude |
            %           filee = '20210203aft03.csv'; % | K_x = 1000, B_x = 20, it_Max = 8, f_pret[0]=0; 5N amplitude |
            
            
            %           filee = '20210203aft04.csv'; % | K_x = 0, B_x = 0, it_Max = 8, pretAmplitude = 0.0 |
            %           filee = '20210203aft05.csv'; % | K_x = 1000, B_x = 20, it_Max = 8, pretAmplitude = 0.0 |
            %           filee = '20210203aft06.csv'; % | K_x = 1000, B_x = 20, it_Max = 8, pretAmplitude =10.0 |
            %           filee = '20210203aft07.csv'; % | K_x = 500,  B_x = 20, it_Max = 8, pretAmplitude =10.0 |
            %           filee = '20210203aft08.csv'; % | K_x = 1000, B_x = 20, it_Max =15, pretAmplitude =10.0 |
            %           filee = '20210203aft09.csv'; % | K_x = 500,  B_x = 20, it_Max =15, pretAmplitude =10.0 |
            %             filee = '20210203aft10.csv';  % | K_x = 500,  B_x = 20, it_Max =15, pretAmplitudex =16.0, pretAmplitudex =8.0 |
            %             filee = '20210203aft11.csv'; % | K_x = 500,  B_x = 20, it_Max =15, pretAmplitudex =16.0, pretAmplitudex =8.0, with I hold |
            %             filee = '20210204aft04.csv'; % K_x=2500, B_x=20, itermax = 8, pretAmplitude = 6
            filee = '20210204aft05.csv';
            
            %             filee = 'noTorquePret.csv';
            %Old trial with bionary random sequency which is either 6 or -6
            % Must be organized to use
            % Import data and convert to X1, X2, Y1, Y2
            
            data = importdata([pathh,filee]);
            
            dexStart = 2e+4;
            dexEnd = size(data,1)-1e+4;
            dexRange = [dexStart:dexEnd];
            t = data(dexRange,1);
            X_m = data(dexRange,10);
            Y_m = data(dexRange,11);
            Pret_X_N = data(dexRange,20);
            Pret_Y_N = data(dexRange,21);
            
            figure;
            subplot(2,1,1); plot(t,X_m); hold on;
            subplot(2,1,2); plot(t,Y_m);hold on;
            
            figure;
            subplot(2,1,1); plot(t,Pret_X_N); hold on;
            subplot(2,1,2); plot(t,Pret_Y_N); hold on;
            
            figure;
            subplot(2,1,1); histogram(Pret_X_N);
            subplot(2,1,2); histogram(Pret_Y_N);
            
            [std(X_m), std(Y_m)]
            
            
            
            % No preturbation [8.26e-2,1.087e-3]
            % 2 N preturbation [1.212e-2, 4.889e-2]
            % 5 N preturbation [7e-4, 2.3e-3]
            
            % figure; plot(1./diff(t),'o');
            
            ensambleLength = 1;
            
            % Define sim parameters
            this.sfrq = 500; % Change later this will change window size
            this.dt = 1/this.sfrq;
            this.tvec = 0:this.dt:ensambleLength+this.dt;
            N = length(this.tvec); % Check N is even
            lagBuffer = 5;
            M1 = 0;
            M2 = 0.3*this.sfrq+lagBuffer;
            L = M2-M1+1;
            %                 R = 50; % number of realizations
            R = floor((length(X_m)-dexStart)/(this.sfrq*ensambleLength)) - 10;
            
            
            %% Filter Postion
            [w, P_x] = this.get_fft(X_m,5);
            [w, P_y] = this.get_fft(Y_m,5);
            
            [w, P_ux] = this.get_fft(Pret_X_N,5);
            [w, P_uy] = this.get_fft(Pret_Y_N,5);
            
            %             cf = 40; % cutoff freqnency
            %             [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            %             X_m = filtfilt(b, a, X_m); % apply fitler
            %             Y_m = filtfilt(b, a, Y_m); % apply fitler
            %
            %             [w, P_x_filt] = this.get_fft(X_m,5);
            %             [w, P_y_filt] = this.get_fft(Y_m,5);
            %             figure;
            %             subplot(2,1,1); plot(w,P_x,'r', w, P_x_filt,'b'); xlim([0 100]); grid on; hold on;
            %             subplot(2,1,2); plot(w,P_y,'r', w, P_y_filt,'b'); xlim([0 100]);
            %             xlabel('Frequency (Hz)');
            %             ylabel('Power/Frequency (dB/Hz)');
            
            figure;
            subplot(2,1,1); plot(w,P_x,w,P_ux,'linewidth',2.5); xlim([0 100]); grid on; hold on; legend('Y','U');
            subplot(2,1,2); plot(w,P_y,w,P_uy,'linewidth',2.5); xlim([0 100]); grid on; hold on;
            xlabel('Frequency (Hz)');
            ylabel('Power/Frequency (dB/Hz)');
            
            figure; plot(X_m,Y_m,'-o'); axis equal;
            xlabel('x');ylabel('y');
            
            
            %% Create Ensamble Data
            
            % Intialize variables
            z1_r = zeros(R,N); % Postion
            z2_r = zeros(R,N); % Postion
            u1_r = zeros(R,N); % Preturbation
            u2_r = zeros(R,N); % Preturbation
            
            % Generate simulated data for ensamble method
            for i = 1:R
                % Interpolate
                dex = (N*i + dexStart):(N*(i+1) + dexStart - 1);
                z1_r(i,:) = interp1(t(dex)-t(dex(1)),detrend(X_m(dex)),this.tvec,'spline');
                z2_r(i,:) = interp1(t(dex)-t(dex(1)),detrend(Y_m(dex)),this.tvec,'spline');
                u1_r(i,:) = interp1(t(dex)-t(dex(1)),Pret_X_N(dex),this.tvec,'nearest');
                u2_r(i,:) = interp1(t(dex)-t(dex(1)),Pret_Y_N(dex),this.tvec,'nearest');
                
            end
            
            u1_r(:,end) = u1_r(:,end-1);
            u2_r(:,end) = u2_r(:,end-1);
            
            
            % Change sampling rate to 500 Hz
            %                 [z_r, u_r, N, L, R, M1, M2] = this.get_lowerSamplingRate(ensambleLength, z_r, u_r, N, L, R, M1, M2);
            
            z1_r_mean = mean(z1_r);
            z2_r_mean = mean(z2_r);
            
            
            for i = 1:R
                z1_r(i,:) = detrend(z1_r(i,:) - z1_r_mean);
                z2_r(i,:) = detrend(z2_r(i,:) - z2_r_mean);
            end
            
            %% Estimate Impedance
            [H_hat_MA, MAwindow] = ensambleSysID_Matrix(this, z1_r, z2_r, u1_r, u2_r, N, L, R, M1, M2);
            
            iDex = floor(linspace(M2+1+40,N+M1-3,10));
            [h_model11, M11, B11, K11, VAFirf11] = this.fitModel(H_hat_MA(:,:,1), N, L, M1, M2,iDex);
            [h_model22, M22, B22, K22, VAFirf22] = this.fitModel(H_hat_MA(:,:,4), N, L, M1, M2,iDex);
            
            
            %% True model
            M11_const = 1;
            M22_const = 1;
            B11_const = 40;
            B22_const = 40;
            K11_const = 2000;
            K22_const = 2000;
            
            s = tf('s');
            h_true_model11 = impulse(inv(M11_const*s^2+B11_const*s+K11_const),0:1/this.sfrq:M2/this.sfrq);
            h_true_model22 = impulse(inv(M22_const*s^2+B22_const*s+K22_const),0:1/this.sfrq:M2/this.sfrq);
            
            %% Plot
            figure;
            ax1 = subplot(2,1,1);
            ax4 = subplot(2,1,2);
            
            for i = iDex(1:end-1)
                axes(ax1);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_hat_MA(:,i,1),'.r'); hold on;
                %                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model11(:),'-g'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model11(:,i),'-b'); hold on;
                
                axes(ax4);
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),H_hat_MA(:,i,4),'.r'); hold on;
                %                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model22(:),'-g'); hold on;
                plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_model22(:,i),'-b'); hold on;
                
            end
            fs = 15;
            axes(ax1); view(35,25); title('11'); xlim([0 M2/this.sfrq]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('x (m)'); grid on; set(gca,'fontsize',fs);%zlim([-0.1 0.1]);
            axes(ax4); view(35,25); title('22'); xlim([0 M2/this.sfrq]); xlabel('lag (s)'); ylabel('Time (s)'); zlabel('y (m)'); grid on; set(gca,'fontsize',fs);%zlim([-0.1 0.1]);
            
            figure;
            labelss = {'11','22'};
            
            ax1 = subplot(3,1,1);
            p1=plot(this.tvec(iDex), M11(iDex),'.-b',this.tvec(iDex), M11_const*ones(size(iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            p2=plot(this.tvec(iDex), M22(iDex),'.-k',this.tvec(iDex), M22_const*ones(size(iDex)),'--k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel('Mass (kg)'); ylim([0 20]);
            legend([p1(1),p2(1)],labelss,'location','northwest');legend boxoff;
            
            ax2 = subplot(3,1,2);
            p1=plot(this.tvec(iDex), B11(iDex),'.-b',this.tvec(iDex), B11_const*ones(size(iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            p4=plot(this.tvec(iDex), B22(iDex),'.-k',this.tvec(iDex), B22_const*ones(size(iDex)),'--k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel(' Damping (N/m)'); ylim([0 150]);
            
            ax3 = subplot(3,1,3);
            p1=plot(this.tvec(iDex), K11(iDex),'.-b',this.tvec(iDex), K11_const*ones(size(iDex)),'--b','linewidth',2.5,'markersize',15); hold on;
            p4=plot(this.tvec(iDex), K22(iDex),'.-k',this.tvec(iDex), K22_const*ones(size(iDex)),'--k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel(' Stiffness (Ns/m)'); ylim([0 9000]);
            
            figure;
            plot(this.tvec(iDex), VAFirf11(iDex),'.-b','linewidth',2.5,'markersize',15); hold on;
            plot(this.tvec(iDex), VAFirf22(iDex),'.-k','linewidth',2.5,'markersize',15);
            grid on;  xlim([0.2 this.tvec(end)]); ylabel('VAF_{irf}'); ylim([0 200]);
            xlabel('Time (s)'); set(gca,'fontsize',18);
            
            disp('test');
            
        end
        
        % Need to remove y_r to runn options other than simData2D
        function [z_r, u_r, N, L, R, M1, M2, y_r] =  importData(this,method,filee)
            %file can also be norm or tan for 2D sim
            
            if(strcmp(method,'simData1D'))
                
                ensambleLength = 2;
                
                % Define sim parameters
                this.sfrq = 2000; % Change later this will change window size
                this.dt = 1/this.sfrq;
                this.tvec = 0:this.dt:ensambleLength+this.dt;
                N = length(this.tvec); % Check N is even
                lagBuffer = 5;
                M1 = 0;
                M2 = 0.2*this.sfrq+lagBuffer;
                L = M2-M1+1;
                R = 400; % number of realizations
                
                % Intialize variables
                y_r = zeros(R,N); % Postion
                z_r = zeros(R,N); % Postion
                u_r = zeros(R,N); % Preturbation
                
                % Generate simulated data for ensamble method
                %                 parfor i = 1:R
                %                     [y_tmp,u_tmp] = this.get_sim();
                %                     y_r(i,:) = y_tmp; % CHECK NORMALIZE
                %                     u_r(i,:) = u_tmp;
                %                     if(mod(i,10)==0)
                %                         disp(i);
                %                     end
                %                 end
                
                % save('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/simValues_Sin_Z.mat','y_r','u_r');
                % load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/simValues_Sin_Z.mat');
                
                
                %%
                figure; plot(this.tvec,u_r(1,:),'linewidth',2.5);
                ylim([-7 7]);xlim([0 0.2]);
                xlabel('Time (s)'); ylabel('U (N)');
                set(gca,'fontsize',18); grid on;
                
                % Add noise
                n_r = 1e-4*randn(size(y_r));
                %                 snr(y_r,n_r)
                z_r = y_r + n_r;
                clear y_r
                
                figure;
                plot(this.tvec,z_r'); hold on; xlim([0 this.tvec(end)]);
                plot(this.tvec,median(z_r),'-k','linewidth',2.5);
                xlabel('Time (s)'); ylabel('x (m)');
                set(gca,'fontsize',18); grid on;
                
                % Change sampling rate to 500 Hz
                %                 [z_r, u_r, N, L, R, M1, M2] = this.get_lowerSamplingRate(ensambleLength, z_r, u_r, N, L, R, M1, M2);
                
                z_r_mean = mean(z_r);
                
                for i = 1:R
                    z_r(i,:) = detrend(z_r(i,:) - z_r_mean);
                end
                
                figure;
                plot(this.tvec,z_r); hold on; xlim([0 this.tvec(end)]);
                xlabel('Time (s)'); ylabel('x (m)');
                set(gca,'fontsize',18); grid on;
                
            elseif (strcmp(method,'fixedAmplitudePRBS'))
                
                %Old trial with bionary random sequency which is either 6 or -6
                % Must be organized to use
                % Import data and convert to X1, X2, Y1, Y2
                load(filee);
                %                 load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/StocasticPoint/Data/UniformNoise_6N_handle_point.mat');
                %                 load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/StocasticPoint/Data/UniformNoise_6N_handle_point_holding.mat');
                
                %                 load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/StocasticPoint/Data/UniformNoise_6N_handle_point_holding.mat');
                %                 load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/ensambleTestData/Gaussian_6N_Handel_15MinCenterPoint_FPGA-17-Nov-2020-16-32.mat');
                %                 load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/ensambleTestData/Uniform_4N_5min_point_FPGA-20-Nov-2020-18-45.mat');
                clear length
                dexStart = 1.424e+4;
                ensambleLength = 1;
                
                % Define sim parameters
                this.sfrq = 2000; % Change later this will change window size
                this.dt = 1/this.sfrq;
                this.tvec = 0:this.dt:ensambleLength+this.dt;
                N = length(this.tvec); % Check N is even
                lagBuffer = 5;
                M1 = 0;
                M2 = 0.3*this.sfrq+lagBuffer;
                L = M2-M1+1;
                %                 R = 50; % number of realizations
                R = floor((length(X_m)-dexStart)/(this.sfrq*ensambleLength)) - 10;
                
                % Intialize variables
                z_r = zeros(R,N); % Postion
                u_r = zeros(R,N); % Preturbation
                
                % Generate simulated data for ensamble method
                for i = 1:R
                    % Interpolate
                    dex = (N*i + dexStart):(N*(i+1) + dexStart - 1);
                    z_r(i,:) = interp1(t(dex)-t(dex(1)),detrend(Y_m(dex)),this.tvec,'spline');
                    u_r(i,:) = interp1(t(dex)-t(dex(1)),Pret_Y_N(dex),this.tvec,'nearest');
                end
                
                u_r(:,end) = u_r(:,end-1);
                
                % Change sampling rate to 500 Hz
                %                 [z_r, u_r, N, L, R, M1, M2] = this.get_lowerSamplingRate(ensambleLength, z_r, u_r, N, L, R, M1, M2);
                
                z_r_median = median(z_r);
                
                for i = 1:R
                    z_r(i,:) = detrend(z_r(i,:) - z_r_median);
                end
                
            elseif(strcmp(method,'gaussianWhite'))
                
                %% Import InMotion data
                %gaussian white noise does not work
                %Import data and convert to X1, X2, Y1, Y2
                load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/StocasticPoint/Data/GaussianNoise_6N_noHandle_point.mat');
                
                %Sanity check
                figure;
                plot(U_X_N); hold on;
                %                 Select points
                [dex_X,dex_Y] = ginput(1);
                dex_X = 15*Hz:(Hz*(15+runningTime));
                dex_X = [round(dex_X):(round(dex_X)+Hz*runningTime)];
                plot([dex_X(1), dex_X(end)],[0,0],'ok','linewidth',2.5);
                
                
                count = 1;
                for i = 1.4e+4:2000:(length(t)-2000)
                    z_r(count,:) = X_m(i:(i+2000));
                    u_r(count,:) = Pret_X_N(i:(i+2000)); % Try with preterbation not U
                    count = count+1;
                end
                R = count-1;
                
            else
                error('Select valid data import method.');
            end
            
        end
        
        function [y_r_1,y_r_2,u_r_1,u_r_2,f_r_1,f_r_2] = get_ensambles_from_simSystem2D(this, R, movDuration,symType,varargin)
            
            if(~isempty(varargin))
                pretScale = varargin;
            else
                pretScale = [];
            end
            
            u_r_1 = [];
            y_r_1 = [];
            f_r_1 = [];
            
            u_r_2 = [];
            y_r_2 = [];
            f_r_2 = [];
            
            pw = PoolWaitbar(R, ['Simulating System 2D ',int2str(R),' times.']);
            
            for i = 1:R
                if(~isempty(pretScale))
                    test = WAM_simSystem2D(movDuration,pretScale);
                else
                    test = WAM_simSystem2D(movDuration);
                end
                test.simulateSystem(symType);
                
                [tmp_y_r_1, tmp_y_r_2, tmp_u_r_1,tmp_u_r_2,tmp_f_r_1, tmp_f_r_2] = test.get_ensambleInput(symType);

                u_r_1 = [u_r_1; tmp_u_r_1];
                y_r_1 = [y_r_1; tmp_y_r_1];
                f_r_1 = [f_r_1; tmp_f_r_1];
                
                u_r_2 = [u_r_2; tmp_u_r_2];
                y_r_2 = [y_r_2; tmp_y_r_2];
                f_r_2 = [f_r_2; tmp_f_r_2];
                
                % disp([int2str(i),'/',int2str(numLoop)]);
                increment(pw); % wait bar display
            end
            
            test = WAM_simSystem2D(movDuration);
            this.sfrq = test.sfrq;
            this.dt = test.dt;
            this.tvec = test.t;
            clear test
             
        end
        
        function [h_hat_MA_filt, MAwindow] = ensambleSysID(this, z_r, u_r, N, L, R, M1, M2)
            
            % Estimates
            h_hat = zeros(L,N);
            
            % MA: moving average with window size MAwindow
            h_hat_MA = zeros(L,N);
            
            Phi_zu_hat = zeros(L,1);
            Phi_uu_hat = zeros(L,L);
            
            tic;
            for i = M2+1:N+M1
                
                for k = M1:M2
                    Phi_zu_hat(k-M1+1,1) = get_PHI(this,z_r,u_r,i,-k,R,M2);
                end
                
                if (i == M2+1)
                    %
                    for k = M1:M2
                        for j = M1:M2
                            Phi_uu_hat(k-M1+1,j-M1+1) = get_PHI(this,u_r,u_r,i-k,k-j,R,M2);
                        end
                    end
                    %
                else
                    Phi_uu_hat(2:L,2:L) = prev_Phi_uu_hat(1:L-1,1:L-1);
                    
                    % Define first column
                    for j = M1:M2
                        k = M1;
                        Phi_uu_hat(k-M1+1,j-M1+1) = get_PHI(this,u_r,u_r,i-k,k-j,R,M2);
                    end
                    
                    % Define first column
                    for k = M1:M2
                        j = M1;
                        Phi_uu_hat(k-M1+1,j-M1+1) = get_PHI(this,u_r,u_r,i-k,k-j,R,M2);
                    end
                    
                end
                %
                prev_Phi_uu_hat = Phi_uu_hat;
                
                %                 figure; plot(Phi_zu_hat);
                %
                %                 figure;
                %                 for i = 1:L
                %                     plot3(1:L,i*ones(1,L),Phi_uu_hat(:,i),'-'); hold on;
                %                 end
                %                 %             plot3(M2,L,X(:,1)','--k','linewidth',2.5);
                %                 xlabel('lag');
                %                 ylabel('i');
                %                 zlabel('magnitude');
                
                % Invert matrix to solve for impulse response function
                
                h_hat(:,i) = this.sfrq*pinv(Phi_uu_hat)*Phi_zu_hat;
                %
                % Use SVD to compute faster?
                %                  [U,S,V]=svd(Phi_uu_hat(:,:));
                %                  h_hat(:,i)=this.sfrq*V*inv(S)*U'*Phi_zu_hat(:,1);
                
                if(mod(i,500)==0)
                    disp([int2str(i),'/',int2str(N+M1)]);
                end
            end
            
            toc;
            
            
            % Moving average with MAstep and MAwindow width
            % h_hat smoothing with moving average window (40 ms)
            MAwindow = 0.04*this.sfrq;  % Moving average window size: 40ms
            MAstep = 0.04*this.sfrq;
            for i = (M2+1+MAwindow/2) : (N+M1-MAwindow/2)
                for j=1:L
                    h_hat_MA(j,i) = mean(h_hat(j,i-MAwindow/2:i+MAwindow/2));
                end
            end
            
            %             figure;plot(h_hat_MA,'linewidth',2);
            %             xlabel('lag (s)');ylabel('magnitude'); set(gca,'fontsize',18);
            %
            %             figure;
            %             for i = M2+1:N+M1
            %                 plot3(this.tvec(1:L),i*ones(1,L),h_hat_MA(:,i),'-','linewidth',2); hold on;
            %             end
            % %             plot3(M2,L,X(:,1)','--k','linewidth',2.5);
            %             xlabel('lag (s)');
            %             ylabel('Time (sec)');
            %             zlabel('magnitude'); set(gca,'fontsize',18);
            
            %             [X,Y] = meshgrid(1:N+M1,1:L);
            %             figure; s = surf(X,Y,h_hat_MA);
            
            cf = 50; % cutoff freqnency
            [b,a] = butter(2,cf/(this.sfrq/2)); % make filter
            h_hat_MA_filt = zeros(L,N);
            dexNonZero = M2-M1+1+MAwindow/2:N+M1-MAwindow/2;
            for i = dexNonZero
                h_hat_MA_filt(:,i) = filtfilt(b,a,h_hat_MA(:,i)); % apply fitler
            end
            
        end
                
        function [H_hat_MA, MAwindow] = ensambleSysID_Matrix(this, z1_r, z2_r, u1_r, u2_r, N, L, R, M1, M2)
            
            
            %% Symbolic santity check
            %             syms Phi_z1u1 Phi_z1u2 Phi_z2u1 Phi_z2u2
            %             syms Phi_u1u1 Phi_u1u2 Phi_u2u1 Phi_u2u2
            %             syms h11 h12 h21 h22
            %             syms z1 z2 u1 u2
            %
            %             [z1;z2] == [h11,h12;h21,h22].'*[u1;u2]
            %
            %
            %             [Phi_z1u1, Phi_z2u1; Phi_z1u2 Phi_z2u2]==[Phi_u1u1, Phi_u2u1; Phi_u1u2, Phi_u2u2]*[h11, h12; h21, h22].'
            
            
            %% Estimates
            h_hat11 = zeros(L,N);
            h_hat12 = zeros(L,N);
            h_hat21 = zeros(L,N);
            h_hat22 = zeros(L,N);
            
            % MA: moving average with window size MAwindow
            h_hat11_MA = zeros(L,N);
            h_hat12_MA = zeros(L,N);
            h_hat21_MA = zeros(L,N);
            h_hat22_MA = zeros(L,N);
            
            Phi_z1u1_hat = zeros(L,1);
            Phi_z1u2_hat = zeros(L,1);
            Phi_z2u1_hat = zeros(L,1);
            Phi_z2u2_hat = zeros(L,1);
            
            Phi_u1u1_hat = zeros(L,L);
            Phi_u1u2_hat = zeros(L,L);
            Phi_u2u1_hat = zeros(L,L);
            Phi_u2u2_hat = zeros(L,L);
            
            tic;
            for i = M2+1:N+M1
                
                for k = M1:M2
                    Phi_z1u1_hat(k-M1+1,1) = get_PHI(this,z1_r,u1_r,i,-k,R,M2);
                    Phi_z1u2_hat(k-M1+1,1) = get_PHI(this,z1_r,u2_r,i,-k,R,M2);
                    Phi_z2u1_hat(k-M1+1,1) = get_PHI(this,z2_r,u1_r,i,-k,R,M2);
                    Phi_z2u2_hat(k-M1+1,1) = get_PHI(this,z2_r,u2_r,i,-k,R,M2);
                end
                
                %                 if (i == M2+1)
                %
                for k = M1:M2
                    for j = M1:M2
                        Phi_u1u1_hat(k-M1+1,j-M1+1) = get_PHI(this,u1_r,u1_r,i-k,k-j,R,M2);
                        Phi_u1u2_hat(k-M1+1,j-M1+1) = get_PHI(this,u1_r,u2_r,i-k,k-j,R,M2);
                        Phi_u2u1_hat(k-M1+1,j-M1+1) = get_PHI(this,u2_r,u1_r,i-k,k-j,R,M2);
                        Phi_u2u2_hat(k-M1+1,j-M1+1) = get_PHI(this,u2_r,u2_r,i-k,k-j,R,M2);
                    end
                end
                %
                %                 else
                %                     Phi_uu_hat(2:L,2:L) = prev_Phi_uu_hat(1:L-1,1:L-1);
                %
                %                     % Define first column
                %                     for j = M1:M2
                %                         k = M1;
                %                         Phi_uu_hat(k-M1+1,j-M1+1) = get_PHI(this,u_r,u_r,i-k,k-j,R,M2);
                %                     end
                %
                %                     % Define first column
                %                     for k = M1:M2
                %                         j = M1;
                %                         Phi_uu_hat(k-M1+1,j-M1+1) = get_PHI(this,u_r,u_r,i-k,k-j,R,M2);
                %                     end
                %
                %                 end
                %
                %                 prev_Phi_uu_hat = Phi_uu_hat;
                
                %                   figure;
                %                   subplot(2,2,1); plot(Phi_z1u1_hat); title(11);
                %                   subplot(2,2,2); plot(Phi_z2u1_hat); title(12);
                %                   subplot(2,2,3); plot(Phi_z1u2_hat); title(21);
                %                   subplot(2,2,4); plot(Phi_z2u2_hat); title(22);
                % %
                % %
                %                 figure;
                %                   ax1 = subplot(2,2,1);
                %                   ax2 = subplot(2,2,2);
                %                   ax3 = subplot(2,2,3);
                %                   ax4 = subplot(2,2,4);
                %                 axes(ax1);
                %                 for i = 1:L
                %                     plot3(1:L,i*ones(1,L),Phi_u1u1_hat(:,i),'-'); hold on;
                %                 end
                %                 axes(ax2);
                %                 for i = 1:L
                %                     plot3(1:L,i*ones(1,L),Phi_u2u1_hat(:,i),'-'); hold on;
                %                 end
                %                 axes(ax3);
                %                 for i = 1:L
                %                     plot3(1:L,i*ones(1,L),Phi_u1u2_hat(:,i),'-'); hold on;
                %                 end
                %                 axes(ax4);
                %                 for i = 1:L
                %                     plot3(1:L,i*ones(1,L),Phi_u2u2_hat(:,i),'-'); hold on;
                %                 end
                %                             plot3(M2,L,X(:,1)','--k','linewidth',2.5);
                %                 xlabel('lag');
                %                 ylabel('i');
                %                 zlabel('magnitude');
                
                % Invert matrix to solve for impulse response function
                
                PHI_uu_hat = [Phi_u1u1_hat',Phi_u2u1_hat';Phi_u1u2_hat',Phi_u2u2_hat'];
                PHI_zu_hat = [Phi_z1u1_hat,Phi_z2u1_hat;Phi_z1u2_hat,Phi_z2u2_hat];
                
                H_hat = this.sfrq*pinv(PHI_uu_hat)*PHI_zu_hat;
                
                h_hat11(:,i) = H_hat(1:L,1);
                h_hat12(:,i) = H_hat(L+1:end,1); % Switch H_hat cross terms to apply transpose
                h_hat21(:,i) = H_hat(1:L,2);
                h_hat22(:,i) = H_hat(L+1:end,2);
                
                % Use SVD to compute faster?
                %                  [U,S,V] = svd(Phi_uu_hat(:,:));
                %                  h_hat(:,i) = this.sfrq*V*inv(S)*U'*Phi_zu_hat(:,1);
                
                if(mod(i,500)==0)
                    disp([int2str(i),'/',int2str(N+M1)]);
                end
            end
            
            toc;
            
            
            % Moving average with MAstep and MAwindow width
            % h_hat smoothing with moving average window (40 ms)
            MAwindow = 0.04*this.sfrq;  % Moving average window size: 40ms
            MAstep = 0.04*this.sfrq;
            for i = (M2+1+MAwindow/2) : (N+M1-MAwindow/2)
                for j=1:L
                    h_hat11_MA(j,i) = mean(h_hat11(j,i-MAwindow/2:i+MAwindow/2));
                    h_hat12_MA(j,i) = mean(h_hat12(j,i-MAwindow/2:i+MAwindow/2));
                    h_hat21_MA(j,i) = mean(h_hat21(j,i-MAwindow/2:i+MAwindow/2));
                    h_hat22_MA(j,i) = mean(h_hat22(j,i-MAwindow/2:i+MAwindow/2));
                end
            end
            
            H_hat_MA(:,:,1) = h_hat11_MA;
            H_hat_MA(:,:,2) = h_hat12_MA;
            H_hat_MA(:,:,3) = h_hat21_MA;
            H_hat_MA(:,:,4) = h_hat22_MA;
            
            %             figure;plot(h_hat_MA,'linewidth',2);
            %             xlabel('lag (s)');ylabel('magnitude'); set(gca,'fontsize',18);
            %
            %             figure;
            %             for i = M2+1:N+M1
            %                 plot3(this.tvec(1:L),i*ones(1,L),h_hat_MA(:,i),'-','linewidth',2); hold on;
            %             end
            % %             plot3(M2,L,X(:,1)','--k','linewidth',2.5);
            %             xlabel('lag (s)');
            %             ylabel('Time (sec)');
            %             zlabel('magnitude'); set(gca,'fontsize',18);
            
            %             [X,Y] = meshgrid(1:N+M1,1:L);
            %             figure; s = surf(X,Y,h_hat_MA);
            
            %             cf = 20; % cutoff freqnency
            %             [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
            %             h_hat11_MA_filt = zeros(L,N);
            %             h_hat12_MA_filt = zeros(L,N);
            %             h_hat21_MA_filt = zeros(L,N);
            %             h_hat22_MA_filt = zeros(L,N);
            %             dexNonZero = M2-M1+1+MAwindow/2:N+M1-MAwindow/2;
            %             for i = dexNonZero
            %                 h_hat11_MA_filt(:,i) = filtfilt(b,a,h_hat11_MA(:,i)); % apply fitler
            %                 h_hat12_MA_filt(:,i) = filtfilt(b,a,h_hat12_MA(:,i)); % apply fitler
            %                 h_hat21_MA_filt(:,i) = filtfilt(b,a,h_hat21_MA(:,i)); % apply fitler
            %                 h_hat22_MA_filt(:,i) = filtfilt(b,a,h_hat22_MA(:,i)); % apply fitler
            %             end
            %
            %             H_hat_MA(:,:,1) = h_hat11_MA_filt;
            %             H_hat_MA(:,:,2) = h_hat12_MA_filt;
            %             H_hat_MA(:,:,3) = h_hat21_MA_filt;
            %             H_hat_MA(:,:,4) = h_hat22_MA_filt;
            
        end
        
        function [h_model, M_hand, B_hand, K_hand, VAFirf] =  fitModel(this, h_hat_MA, N, L, M1, M2, iDex)
            
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
            
            LB = [-20, -200, -10000];
            UB = [20,  200, 10000];
            
            % Bounded Nonlinear Optimization
            h_model = zeros(L,N);
            M_hand = zeros(N,1);
            B_hand = zeros(N,1);
            K_hand = zeros(N,1);
            fval = zeros(N,1);
            exitflags = zeros(N,1);
            buffer = 3;
            lagBuffer = 5;
            %         skip = 300;%floor(50/2);
            %         iDex = M2+1+40:skip:N+M1-3; %-100;
            count = 1;
            for i = iDex
                for j=-buffer:buffer
                    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun',0.01);
                    [x fval(i+j,1) exitflags(i+j,1) outputs] = fminsearchbnd(@(x) this.irfFitting(x,h_hat_MA(:,i+j),M1,M2,lagBuffer),[M0 B0 K0],LB,UB,options);
                    M_opt = x(1);   B_opt = x(2);   K_opt = x(3);
                    
                    % FIX CHECK fit cost function
                    
                    Y_model = tf(1,[M_opt,B_opt,K_opt]);
                    h_model(:,i+j) = impulse(Y_model,0:1/this.sfrq:M2/this.sfrq);
                    
                    M_hand(i+j,1) = M_opt;
                    B_hand(i+j,1) = B_opt;
                    K_hand(i+j,1) = K_opt;
                    
                    M0 = M_opt;     B0 = B_opt;     K0 = K_opt;
                end
                if(mod(count,5)==0)
                    disp([int2str(count),'/',int2str(length(iDex))]);
                end
                count = count+1;
            end
            
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
            %
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
        
        function fval = irfFitting(this,x,h_hat,M1,M2,lagBuffer)
            M2=M2-lagBuffer;
            L=M2-M1+1;
            Y_model = tf(1,[x(1),x(2),x(3)]);
            h_model = impulse(Y_model,this.tvec);
            
            fval=0;
            for j=M1:M2
                %         fval = fval + (h_hat(j-M1+1,1)-h_model(j-M1+1,1))^2;
                
                % Weighted cost function 2
                if(j<M2/2)
                    fval = fval + (h_hat(j-M1+1,1)-h_model(j-M1+1,1))^2;
                else
                    weight = 2*(1-j/M2);
                    fval = fval + (h_hat(j-M1+1,1)-h_model(j-M1+1,1))^2*weight;
                end
            end
            fval = sqrt(fval/L);
        end
        
        function [H_model, M_hat, B_hat, K_hat, VAFirf] =  fitModel_Matrix(this, H_hat_MA, N, L, M1, M2, iDex)
            
            %             M0 = [1.7099, -0.2566; -0.2566, 2.1775];
            %             B0 = [5.2510, -1.0215; -1.0215, 39.0782];
            %             K0 = [105.0196, -20.4292; -20.4292, 781.5645];
            %
            %             XX0 = [M0(1:4),B0(1:4),K0(1:4)]';
            %
            %             LB = [-20*ones(1,4), -200*ones(1,4), -10000*ones(1,4)]';
            %             UB = [ 20*ones(1,4),  200*ones(1,4),  10000*ones(1,4)]';
            
            M0 = [1.7099, -0.2566, 2.1775];
            B0 = [5.2510, -1.0215, 39.0782];
            K0 = [105.0196, -20.4292, 781.5645];
            
            XX0 = [M0(1:3),B0(1:3),K0(1:3)]';
            
            LB = [-20*ones(1,3), -200*ones(1,3), -10000*ones(1,3)]';
            UB = [ 20*ones(1,3),  200*ones(1,3),  10000*ones(1,3)]';
            
            
            s = tf('s');
            H_model = zeros(size(H_hat_MA,1),size(H_hat_MA,2),4);
            
            for i = iDex
                % Change H_hat structure
                H_hat = zeros(size(H_hat_MA,1),2,2);
                H_hat(:,1,1) = H_hat_MA(:,i,1);
                H_hat(:,1,2) = H_hat_MA(:,i,2);
                H_hat(:,2,1) = H_hat_MA(:,i,3);
                H_hat(:,2,2) = H_hat_MA(:,i,4);
                
                %                 [XX,fval] = fmincon(@(x) this.irfFitting_Matrix(x, H_hat, M1, M2), XX0, zeros(9), zeros(9,1),zeros(9), zeros(9,1), LB, UB);
                options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun',0.01);
                [XX fval exitflags outputs] = fminsearchbnd(@(x) this.irfFitting_Matrix(x,H_hat,M1,M2),XX0,LB,UB,options);
                
                
                %                 M_hat(i,:,:) = [XX(1:2)'; XX(3:4)'];
                %                 B_hat(i,:,:) = [XX(5:6)'; XX(7:8)'];
                %                 K_hat(i,:,:) = [XX(9:10)'; XX(11:12)'];
                
                M_hat = [XX(1:2)'; XX(2:3)'];
                B_hat = [XX(4:5)'; XX(5:6)'];
                K_hat = [XX(7:8)'; XX(8:9)'];
                
                Y_model = inv(M_hat*s^2+B_hat*s+K_hat);
                H_model_temp = impulse(Y_model,0:1/this.sfrq:M2/this.sfrq);
                H_model(:,i,1) = H_model_temp(:,1,1);
                H_model(:,i,2) = H_model_temp(:,1,2);
                H_model(:,i,3) = H_model_temp(:,2,1);
                H_model(:,i,4) = H_model_temp(:,2,2);
                
            end
            
        end
        
        function [fval] = irfFitting_Matrix(this,XX,H_hat,M1,M2)
            
            s = tf('s');
            %             M = [XX(1:2)'; XX(3:4)'];
            %             B = [XX(5:6)'; XX(7:8)'];
            %             K = [XX(9:10)'; XX(11:12)'];
            
            M = [XX(1:2)'; XX(2:3)'];
            B = [XX(4:5)'; XX(5:6)'];
            K = [XX(7:8)'; XX(8:9)'];
            
            Y_model = inv(M*s^2+B*s+K);
            H_model = impulse(Y_model,0:1/this.sfrq:M2/this.sfrq);
            
            weight = 2*(1-[1:M2+1]'/M2);
            weight(1:floor(M2/2)) = 1;
            
            %             fval = sum(sum(sum((H_hat - H_model).^2)));
            tmp = (H_hat - H_model).^2;
            fval = sum(weight.*tmp(:,1,1),1) + ...
                sum(weight.*tmp(:,2,1),1) + ...
                sum(weight.*tmp(:,2,2),1);
            
        end
        
        function [x,u] = get_sim(this)
            
            % Define ZFT
            %           [this.X_0] = this.getMinJerkTraj(1,1,0,tvec);
            N = length(this.tvec);
            %             this.X_0 = zeros(1,N);
            this.X_0 = 0.01*sin(1*pi*this.tvec); % variable X_0
            
            %             % Create bandwidthlimited 50 Hz gausian white noise (std before or after filt?)
            fp = 10*rand(1,N);
            fp = fp-mean(fp);
            
            % Filter
            cf = 200; % cutoff freqnency
            [b,a] = butter(2,cf/(this.sfrq/2)); % make filter
            fp = filtfilt(b,a,fp); % apply fitler
            
            %             for i = 1:N
            %                 if(~mod(i-1,10))
            %                     fp(i) = 10*(rand(1)-0.5);
            %                 else
            %                     fp(i) = fp(i-1);
            %                 end
            %             end
            
            % Make bionary amplitude
            fp(find(fp >= 0)) = 5;
            fp(find(fp < 0)) = -5;
            
            % Check Disturbance
            %             figure; plot(this.tvec,fp,'-o');
            %             figure; histogram(fp);
            %             [f,P1] = this.get_fft(fp);
            %             figure; plot(f,P1);
            
            % Create input
            this.U_tmp = [this.X_0;fp];
            
            [tmp,X] = ode45(@(tt,X)this.getStates(tt,X),this.tvec,[0,0]);
            X = X';
            
            x = X(1,:); % Postion
            u = this.U_tmp(2,:); % Preturbation
            
            %             figure;
            %             subplot(2,1,1); plot(this.tvec,u);
            %             subplot(2,1,2); plot(this.tvec,x);
            
            %             figure;
            %             subplot(2,1,1); plot(this.tvec,X(1,:));
            %             subplot(2,1,2); plot(this.tvec,X(2,:));
            
        end
        
        function [X,U] = get_sim_linear(this)
            
            % Define ZFT
            %           [this.X_0] = this.getMinJerkTraj(1,1,0,tvec);
            N = length(this.tvec);
            this.X_0 = zeros(1,N);
            %             this.X_0 = 0.01*sin(1*pi*this.tvec); % variable X_0
            
            %             % Create bandwidthlimited 50 Hz gausian white noise (std before or after filt?)
            fp = 10*rand(2,N);
            fp = fp-mean(fp,2);
            
            % Filter
            cf = 75; % cutoff freqnency
            [b,a] = butter(2,cf/(this.sfrq/2)); % make filter
            fp(1,:) = filtfilt(b,a,fp(1,:)); % apply fitler
            fp(2,:) = filtfilt(b,a,fp(2,:)); % apply fitler
            
            %             for i = 1:N
            %                 if(~mod(i-1,3))
            %                     fp(1,i) = sign(rand(1)-0.5);
            %                     fp(2,i) = sign(rand(1)-0.5);
            %                 else
            %                     fp(1,i) = fp(1,i-1);
            %                     fp(2,i) = fp(2,i-1);
            %                 end
            %             end
            %
            % Make bionary amplitude
            fp(1,find(fp(1,:) >= 0),:) = 5;
            fp(1,find(fp(1,:) < 0),:) = -5;
            
            fp(2,find(fp(2,:) >= 0),:) = 5;
            fp(2,find(fp(2,:) < 0),:) = -5;
            
            
            % Check Disturbance
            %             figure; plot(this.tvec,fp,'-o');
            %             figure; histogram(fp);
            %             [f,P1] = this.get_fft(fp(1,:),1);
            %             figure; plot(f,P1);
            %               figure; plot(this.tvec,fp(1,:),...
            %                   this.tvec(1):1/2000:this.tvec(end),interp1(this.tvec,fp(1,:),this.tvec(1):1/2000:this.tvec(end),'near'));
            %
            % Create input
            this.U_tmp = fp;
            
            [tmp,XX] = ode45(@(tt,X)this.getStates_linear(tt,X),this.tvec,[0;0;0;0]);
            XX = XX';
            X = XX(1:2,:); % Postion
            U = this.U_tmp(1:2,:); % Preturbation
            
            %             figure;
            %             subplot(2,1,1); plot(this.tvec,u);
            %             subplot(2,1,2); plot(this.tvec,x);
            
            %             figure;
            %             subplot(2,1,1); plot(this.tvec,X(1,:));
            %             subplot(2,1,2); plot(this.tvec,X(2,:));
            
        end
        
        function [PHI] = get_PHI(this,a,b,i,k,R,M2)
            
            %             w_a = 0.54 - 0.46*cos(2*pi*(i/M2));
            %             w_b = 0.54 - 0.46*cos(2*pi*(k/M2));
            
            PHI = (1/R)*sum(a(:,i).*b(:,i+k));
            
        end
        
        function [z_r, u_r, R] = importMultiHumanTrials(this,pathh,fileNames)
            
            cycNum = 0;
            Pret_X_N_tot = [];
            Pret_norm_tot = [];
            X_m_tot = [];
            X_radial_tot = [];
            
            for ii = 1:length(fileNames)
                load([pathh,fileNames{ii}]); % Check sign of X!!
                clear length i
                
                % get one cycle
                p = atan2(Y_m, X_m);
                indxs = find(p(:,1)<0);
                p(indxs,1) = p(indxs) + 2*pi;
                
                dexReset = find(abs(diff(p)) > 4);
                % Add back for slow case
                %                 dexReset = dexReset(find(diff(dexReset)>2000*2)); % elliminate multiple samples at same point
                dexReset = dexReset(3:end-2); % Cut first 2 and last cycle
                dexRange = [dexReset(1):dexReset(end)];
                cycNum_current = size(dexReset,2)-1;
                cycNum = cycNum + cycNum_current;
                
                Pret_norm =  Pret_X_N.*cos(p) + Pret_Y_N.*sin(p);
                Pret_tan =  -Pret_X_N.*sin(p) + Pret_Y_N.*cos(p);
                X_radial = sqrt(X_m.^2+Y_m.^2);
                
                Pret_X_N_tot = [Pret_X_N_tot,Pret_X_N(dexRange)];
                X_m_tot = [X_m_tot,X_m(dexRange)];
                X_radial_tot = [X_radial_tot, X_radial(dexRange)];
                Pret_norm_tot = [Pret_norm_tot, Pret_norm(dexRange)];
                
                if (ii == 1)
                    dexReset_tot = [dexReset - dexRange(1) + 1];
                else
                    dexReset_tot = [dexReset_tot,[dexReset(2:end) - dexRange(2) + 1 + dexReset_tot(end)]];
                end
                
            end
            
            R = cycNum-1;
            
            % Intialize variables
            N = length(this.tvec);
            z_r = zeros(R,N); % Postion
            u_r = zeros(R,N); % Preturbation
            
            for i = 1:R
                z_r(i,:) = X_norm_tot(dexReset_tot(i):dexReset_tot(i) + N - 1 );
                u_r(i,:) = Pret_norm_tot(dexReset_tot(i):dexReset_tot(i) + N - 1 );
            end
            
        end
        
        function [z1_r,z2_r, u1_r, u2_r, R] = importMultiHumanTrials_TimeScale(this,pathh,fileNames,ensambleLength)
            
            X_m_tot = [];
            Y_m_tot = [];
            X_norm_tot = [];
            X_tan_tot = [];
            
            Pret_X_N_tot = [];
            Pret_Y_N_tot = [];
            Pret_norm_tot = [];
            Pret_tan_tot = [];
            
            
            for ii = 1:length(fileNames)
                load([pathh,fileNames{ii}]); % Check sign of X!!
                clear length i
                
                % get one cycle
                p = atan2(Y_m, X_m);
                indxs = find(p<0); % CHECK
                p(indxs) = p(indxs) + 2*pi;
                
                dexReset = find(abs(diff(p)) > 4);
                
                if(ensambleLength>10) % In slow case
                    dexReset = dexReset(find(diff(dexReset)>2000*2)); % elliminate multiple samples at same point
                end
                
                dexReset = dexReset(3:end-2); % Cut first 2 and last cycle
                dexRange = [dexReset(1):dexReset(end)];
                cycNum_current = size(dexReset,2)-1;
                
                X_norm = sqrt(X_m.^2+Y_m.^2);
                X_tan = p*0.1;
                
%                 % Temp Test
%                 Pret_X_N = cos(p) - sin(p);
%                 Pret_Y_N = sin(p) + cos(p);
                
                % Preturbations in x y corrinates
%                 Pret_norm = Pret_X_N.*cos(p)    + Pret_Y_N.*(sin(p));
%                 Pret_tan =  Pret_X_N.*(-sin(p)) + Pret_Y_N.*cos(p);

                % Preturbations in normal tangetial coorinates
                Pret_norm = Pret_X_N;
                Pret_tan = Pret_Y_N;                
                
                % Plot verse state
%                 figure;
%                 subplot(2,1,1); plot(p,Pret_X_N,'o',p,Pret_norm,'.'); legend('xy','ne');
%                 subplot(2,1,2); plot(p,Pret_Y_N,'o',p,Pret_tan,'.');
%                 
%                 figure;plot(X_m(1),Y_m(1),'o',[0,0.1*cos(p(1))], [0,0.1*sin(p(1))],'-'); ylim([-0.2 0.2]); xlim([-0.2 0.2]); grid on; axis equal;
% 
%                 figure; plot(X_tan(1:end-1),this.sfrq.*diff(X_tan),'.'); ylim([-0.25 0.25]);
%                 figure; plot(X_tan,X_norm,'.'); ylim([-0.25 0.25]);
                
                 %% Chop up equal Need to write

                 %% Time aline
                for j = 1:cycNum_current
                    dex = [dexReset(j):dexReset(j+1)];
                    t_raw = t(dex)-t(dex(1));
                    N_raw = length(dex);
                    t_resample = linspace(0,N_raw*this.dt-this.dt*2,length(this.tvec));
                    
                    X_m_trial(j,:) =  interp1(t_raw,X_m(dex),t_resample,'spline');
                    Y_m_trial(j,:) =  interp1(t_raw,Y_m(dex),t_resample,'spline');
                    X_norm_trial(j,:) =  interp1(t_raw,X_norm(dex),t_resample,'spline');
                    X_tan_trial(j,:) =  interp1(t_raw,X_tan(dex),t_resample,'spline');
                    
                    Pret_X_N_trial(j,:) =  interp1(t_raw,Pret_X_N(dex),t_resample,'nearest');
                    Pret_Y_N_trial(j,:) =  interp1(t_raw,Pret_Y_N(dex),t_resample,'nearest');
                    Pret_norm_trial(j,:) =  interp1(t_raw,Pret_norm(dex),t_resample,'nearest');
                    Pret_tan_trial(j,:) =  interp1(t_raw,Pret_tan(dex),t_resample,'nearest');
                    
                end

                %% Make Ensembles
                
                X_m_tot = [X_m_tot;X_m_trial];
                Y_m_tot = [Y_m_tot;Y_m_trial];
                X_norm_tot = [X_norm_tot; X_norm_trial];
                X_tan_tot = [X_tan_tot; X_tan_trial];
                
                Pret_X_N_tot = [Pret_X_N_tot; Pret_X_N_trial];
                Pret_Y_N_tot = [Pret_Y_N_tot; Pret_Y_N_trial];
                Pret_norm_tot = [Pret_norm_tot; Pret_norm_trial];
                Pret_tan_tot = [Pret_tan_tot; Pret_tan_trial];
                
                clear X_m_trial Y_m_trial X_norm_trial X_tan_trial
                clear Pret_X_N_trial Pret_Y_N_trial Pret_norm_trial Pret_tan_trial
                
            end
            
            R = size(X_m_tot,1);
            
            z1_r = X_norm_tot;
            z2_r = X_tan_tot;

            u1_r = Pret_norm_tot;
            u2_r = Pret_tan_tot;
            
            u1_r(:,end-1) = u1_r(:,end-2);
            u1_r(:,end) = u1_r(:,end-2);
            u2_r(:,end-1) = u2_r(:,end-2);
            u2_r(:,end) = u2_r(:,end-2);
            
%             figure;plot(z1_r(1,1:1000)','o');
%             figure;plot(z2_r(1,1:1000)','o');
%             for i = 1:R
%                 z1_plot(i,:) = smooth(this.sfrq*diff(z1_r(i,:)));
%                 z2_plot(i,:) = smooth(this.sfrq*diff(z2_r(i,:)));
% 
%             end
%             figure;plot(z2_r(:,1:end-1)',z1_plot','.');ylim([-0.25 0.25]);
%             figure;plot(z2_r(:,1:end-1)',z2_plot','.');ylim([-0.25 0.25]);
            
            %% Check preturbation interpolation
            figure;plot(u1_r(1,:)','o');
            figure;plot(u2_r(1,:)','o');
           
        end
        
        function [K_hat, B_hat, M_hat] = get_parametricLeastSquares(this,iDex, z_r_1, z_r_2, u_r_1, u_r_2)
            
            [R,N] = size(z_r_1);
            
%             cf = 40; % cutoff freqnency
%             [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
%             for r = 1:R
%                 z_r_1(r,:) = filtfilt(b,a,z_r_1(r,:)); % apply fitler
%                 z_r_2(r,:) = filtfilt(b,a,z_r_2(r,:)); % apply fitler
%             end

            K_hat = zeros(N,3);
            B_hat = zeros(N,3);
            M_hat = zeros(N,3);
            
            L = 50;
            for j = iDex(1:end-1) % Time index
                A = [];
                F = [];
                for i = j:j+L % create bin of data to least square fit
                    for r = 1:R % Replication
                        
                        x1 = z_r_1(r,i);
                        x_dot1 = this.sfrq.*diff(z_r_1(r,i:i+1));
                        x_ddot1 = this.sfrq*this.sfrq.*diff(diff(z_r_1(r,i:i+2)));
                        
                        x2 = z_r_2(r,i);
                        x_dot2 = this.sfrq.*diff(z_r_2(r,i:i+1));
                        x_ddot2 = this.sfrq*this.sfrq.*diff(diff(z_r_2(r,i:i+2)));
                        
                        A = [A;...
                            x1, x2, 0, 0, x_dot1, x_dot2, 0, 0, x_ddot1, x_ddot2, 0, 0;...
                            0, 0, x1, x2, 0, 0, x_dot1, x_dot2, 0, 0, x_ddot1, x_ddot2];

                        F = [F;u_r_1(r,i);u_r_2(r,i)];
                    end
                end
                
                p(j,:) = (A'*A)\A'*F;
                
                K_hat(j,:) = p(j,1:3);
%                 B_hat(j,:) = p(j,4:6);
%                 M_hat(j,:) = p(j,7:9);
%                 K_hat(j,:) = p(j,1:4);
%                 B_hat(j,:) = p(j,5:8);
%                 M_hat(j,:) = p(j,9:12);
                
                disp(j);
            end
            
            
%             figure; plot(p(:,1),'o');
%             K_hat_dev = [std(squeeze(p(iDex,1))')',...
%                          std(squeeze(p(iDex,2))')',...
%                          std(squeeze(p(iDex,3))')',...
%                          std(squeeze(p(iDex,4))')'];
% 
%             B_hat_dev = [squeeze(std(p(:,:,5),2)),...
%                          squeeze(std(p(:,:,6),2)),...
%                          squeeze(std(p(:,:,7),2)),...
%                          squeeze(std(p(:,:,8),2))];
%                      
%             figure; 
%             plot(K_hat(:,1),'linewidth',2.5); hold on;
%             plot(K_hat(:,1)+K_hat_dev(:,1),'--');
%             plot(K_hat(:,1)-K_hat_dev(:,1),'--');

% syms x1 x_dot1 x_ddot1 x2 x_dot2 x_ddot2
% syms m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22
% X = [x1;x2];
% X_dot = [x_dot1;x_dot2];
% X_ddot = [x_ddot1;x_ddot2];
% M = [m11,m12;m21,m22];
% B = [b11,b12;b21,b22];
% K = [k11,k12;k21,k22];
% 
% M*X_ddot + B*X_dot + K*X
% 
% % Assume symetry
% h = [m11 m12 m22 b11 b12 b22 k11 k12 k22].';
% z = [x_ddot1 x_ddot2 0 x_dot1 x_dot2 0 x1 x2 0;...
%      0 x_ddot1 x_ddot2 0 x_dot1 x_dot2 0 x1 x2];

            
        end
        
        function [K_hat, B_hat, M_hat] = get_parametricLeastSquares_1D(this,iDex, z_r_1, z_r_2, u_r_1, u_r_2)
            
            [R,N] = size(z_r_1);
            
%             cf = 40; % cutoff freqnency
%             [b,a] = butter(4,cf/(this.sfrq/2)); % make filter
%             for r = 1:R
%                 z_r_1(r,:) = filtfilt(b,a,z_r_1(r,:)); % apply fitler
%                 z_r_2(r,:) = filtfilt(b,a,z_r_2(r,:)); % apply fitler
%             end

            K_hat = zeros(N,1);
            B_hat = zeros(N,1);
            M_hat = zeros(N,1);
            

            L = 50;
            pw = PoolWaitbar(R, ['Least Squares fit ',int2str(iDex(end-1)),' times.']);
            for j = iDex(1:end-1) % Time index
                A = [];
                F = [];
                for i = j:j+L % create bin of data to least square fit
                    for r = 1:R % Replication
                        
                        x1 = z_r_1(r,i);
                        x_dot1 = this.sfrq.*diff(z_r_1(r,i:i+1));
                        x_ddot1 = this.sfrq*this.sfrq.*diff(diff(z_r_1(r,i:i+2)));
                        
                        A = [A;...
                            x1, x_dot1, x_ddot1];

                        F = [F;u_r_1(r,i)];
                    end
                end
                
                p(j,:) = (A'*A)\A'*F;
                
                K_hat(j,:) = p(j,1);
                B_hat(j,:) = p(j,2);
                M_hat(j,:) = p(j,3);
                
                increment(pw);
            end
            
            
%             figure; plot(p(:,1),'o');
%             K_hat_dev = [std(squeeze(p(iDex,1))')',...
%                          std(squeeze(p(iDex,2))')',...
%                          std(squeeze(p(iDex,3))')',...
%                          std(squeeze(p(iDex,4))')'];
% 
%             B_hat_dev = [squeeze(std(p(:,:,5),2)),...
%                          squeeze(std(p(:,:,6),2)),...
%                          squeeze(std(p(:,:,7),2)),...
%                          squeeze(std(p(:,:,8),2))];
%                      
%             figure; 
%             plot(K_hat(:,1),'linewidth',2.5); hold on;
%             plot(K_hat(:,1)+K_hat_dev(:,1),'--');
%             plot(K_hat(:,1)-K_hat_dev(:,1),'--');

% syms x1 x_dot1 x_ddot1 x2 x_dot2 x_ddot2
% syms m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22
% X = [x1;x2];
% X_dot = [x_dot1;x_dot2];
% X_ddot = [x_ddot1;x_ddot2];
% M = [m11,m12;m21,m22];
% B = [b11,b12;b21,b22];
% K = [k11,k12;k21,k22];
% 
% M*X_ddot + B*X_dot + K*X
% 
% % Assume symetry
% h = [m11 m12 m22 b11 b12 b22 k11 k12 k22].';
% z = [x_ddot1 x_ddot2 0 x_dot1 x_dot2 0 x1 x2 0;...
%      0 x_ddot1 x_ddot2 0 x_dot1 x_dot2 0 x1 x2];

            
        end

        function [z1_r, z2_r, u1_r, u2_r, N, L, R, M1, M2] = get_lowerSamplingRate(this,ensambleLength, z1_r, z2_r, u1_r, u2_r, N, L, R, M1, M2)
            
            this.sfrq = 200; % Change later this will change window size
            this.dt = 1/this.sfrq;
            tvec_old = this.tvec;
            this.tvec = 0:this.dt:ensambleLength+this.dt;
            N = length(this.tvec); % Check N is even
            lagBuffer = 5;
            M1 = 0;
            M2 = 0.5*this.sfrq+lagBuffer;
            L = M2-M1+1;
            
            % Intialize variables
            tmp_z1_r = zeros(R,N); % Postion
            tmp_z2_r = zeros(R,N); % Postion
            tmp_u1_r = zeros(R,N); % Preturbation
            tmp_u2_r = zeros(R,N); % Preturbation

            
            for i = 1:R
                tmp_z1_r(i,:) = interp1(tvec_old, z1_r(i,:), this.tvec,'spline');
                tmp_z2_r(i,:) = interp1(tvec_old, z2_r(i,:), this.tvec,'spline');

                tmp_u1_r(i,:) = interp1(tvec_old, u1_r(i,:), this.tvec,'nearest');
                tmp_u2_r(i,:) = interp1(tvec_old, u2_r(i,:), this.tvec,'nearest');

            end
            
            % Quick fix for last value
            tmp_u1_r(:,end) = u1_r(:,end);
            tmp_u2_r(:,end) = u2_r(:,end);

            % Intialize variables
            z1_r = zeros(R,N); % Postion
            z2_r = zeros(R,N); % Postion

            u1_r = zeros(R,N); % Preturbation
            u2_r = zeros(R,N); % Preturbation

            z1_r = tmp_z1_r;
            z2_r = tmp_z2_r;
            u1_r = tmp_u1_r;
            u2_r = tmp_u2_r;
            
        end
        
        function [x,v,a] = getMinJerkTraj(this,D,A,tstart,tt)
            
            % syms c0 c1 c2 c3 c4 c5
            % [c0, c1, c2, c3, c4, c5] = solve([c0 == 0,...
            %     c1 == 0,...
            %     c2 == 0,...
            %     A == c3*D^3 + c4*D^4 + c5*D^5,...
            %     0 == 3*c3*D^2 + 4*c4*D^3 + 5*c5*D^4,...
            %     0 == 6*c3*D + 12*c4*D^2 + 20*c5*D^3], [c0, c1, c2, c3, c4, c5]);
            
            dex = find(tt > tstart & tt < tstart+D);
            if(isempty(dex)) % if single index
                if(tt < tstart)
                    x = 0;
                    v = 0;
                    a = 0;
                elseif(tt > tstart+D)
                    t = D;
                    x = double(A*( (10/D^3)*t.^3 + (-15/D^4)*t.^4 + (6/D^5)*t.^5 ));
                    v = 0;
                    a = 0;
                else
                    t = tt-tstart;
                    x = double(A*( (10/D^3)*t.^3 + (-15/D^4)*t.^4 + (6/D^5)*t.^5 ));
                    v = double(A*( (30/D^3)*t.^2 + (-60/D^4)*t.^3+ (30/D^5)*t.^4 ));
                    a = double(A*( (60/D^3)*t + (-180/D^4)*t.^2+ (120/D^5)*t.^3 ));
                end
            else % is multiple index
                
                x = zeros(size(tt));
                v = zeros(size(tt));
                a = zeros(size(tt));
                
                t = tt(dex)-tstart;
                x(dex) = double(A*( (10/D^3)*t.^3 + (-15/D^4)*t.^4 + (6/D^5)*t.^5 ));
                v(dex) = double(A*( (30/D^3)*t.^2 + (-60/D^4)*t.^3+ (30/D^5)*t.^4 ));
                a(dex) = double(A*( (60/D^3)*t + (-180/D^4)*t.^2+ (120/D^5)*t.^3 ));
                
                x(dex(end):end) = x(dex(end))*ones(size(x(dex(end):end)));
                
            end
            
            % figure;
            % subplot(3,1,1); plot(t,x);
            % title('Problem 1:  Part 2');
            % xlabel('Time (sec)');
            % ylabel('Position (rad)');
            % legend('Minimum Jerk','Cycloidan Motion','Location','southeast');
            %
            % subplot(3,1,2); plot(t,v);
            % xlabel('Time (sec)');
            % ylabel('Velocity (rad/sec)');
            % legend('Minimum Jerk','Cycloidan Motion','Location','southeast');
            %
            % subplot(3,1,3); plot(t,a);
            % xlabel('Time (sec)');
            % ylabel('Acceleration (rad/sec^2)');
            % legend('Minimum Jerk','Location','southeast');
        end
        
        function [X_dot,extra] = getStates(this,tt,X)
            
            % system values
            k = 2500 + 500*sin(2*pi*0.5*tt);
            b = 40 + 10*cos(2*pi*0.5*tt);
            m = 1;
            
            Wn = sqrt(k/m);
            zeta = (b/m)*(1/(2*Wn));
            extra = [Wn,zeta];
            %             ts = 4/(zeta*Wn);
            
            A = [0,1;-k/m,-b/m];
            B = [0,0;k/m,1/m];
            
            % Get time varying input
            U = [interp1(this.tvec,this.U_tmp(1,:),tt,'spline');...
                interp1(this.tvec,this.U_tmp(2,:),tt,'spline')];
            
            X_dot = A*X + B*U;
            
        end
        
        function [XX_dot,extra] = getStates_linear(this,tt,XX)
            
            X = XX(1:2);
            X_dot = XX(3:4);
            
            % system values
            M = [1.7099, -0.2566; -0.2566, 2.1775];
            B = [5.2510, -1.0215; -1.0215, 39.0782];
            K = [105.0196, -20.4292; -20.4292, 781.5645];
            
            
            %             Wn = sqrt(K./M);
            %             zeta = (B./M).*(1./(2*Wn));
            %             extra = [Wn,zeta];
            %             ts = 4./(zeta*Wn);
            
            % Get time varying input
            Fp = [interp1(this.tvec,this.U_tmp(1,:),tt,'near');...
                interp1(this.tvec,this.U_tmp(2,:),tt,'near')];
            
            X_ddot = inv(M)*( K*([0;0] - X) + B*([0;0] - X_dot) + Fp );
            
            XX_dot = [X_dot;X_ddot];
            
        end
        
        function [X_dot,extra] = getStatesStep1(this,tt,X)
            
            % system values
            %             k = 30;
            %             b = 1.1;
            %             m = 0.04;
            
            k = 2500;
            b = 40;
            m = 2;
            
            Wn = sqrt(k/m);
            zeta = b/(2*Wn*m);
            extra = [Wn,zeta];
            ts = 4/(zeta*Wn);
            
            A = [0,1;-k/m,-b/m];
            B = [0,0;k/m,1/m];
            
            %             figure; impulse(ss(A,B,eye(2),zeros(2,2)));
            
            % Get time varying input
            
            if(tt < 0.001)
                pret = 10;
            else
                pret = 0;
            end
            
            U = [0;...
                pret];
            
            X_dot = A*X + B*U;
            
        end
        
        function [X_dot,extra] = getStatesStep2(this,tt,X)
            
            % system values
            k = 15;
            b = 1;
            m = 0.03;
            
            Wn = sqrt(k/m);
            zeta = (b/m)*(1/(2*Wn));
            %             extra = [Wn,zeta];
            %             ts = 4/(zeta*Wn);
            
            A = [0,1;-k/m,-b/m];
            B = [0,0;k/m,1/m];
            
            % Get time varying input
            
            if(tt < 0.001)
                pret = 10;
            else
                pret = 0;
            end
            
            U = [0;...
                pret];
            
            X_dot = A*X + B*U;
            
        end
        
        function [] = get_ideaImpulseResponse(this)
            
            % figure; plot(linspace(0,0.2010, 201),h');hold on;
            [tmp,X1] = ode45(@(tt,X)this.getStatesStep1(tt,X),this.tvec,[0,0]);
            [tmp,X2] = ode45(@(tt,X)this.getStatesStep2(tt,X),this.tvec,[0,0]);
            
            figure; plot(this.tvec,X1(:,1),'--r',...
                this.tvec,X2(:,1),'--b','linewidth',2.5);
            xlim([0 0.2]);
            
        end
        
        function [] = comparePret(this,N)
            
            X1 = rand(1,N);
            X1 = X1 - mean(X1);
            
            X2 = randn(1,N);
            
            X3 = rand(1,N);
            X3 = X3-mean(X3);
            dex1 = find(X3 > 0);
            dex2 = find(X3 < 0);
            
            [f,P1] = this.get_fft(X1);
            [f,P2] = this.get_fft(X1);
            [f,P3] = this.get_fft(X1);
            
            figure;
            subplot(3,1,1); plot(this.tvec,X1);
            subplot(3,1,2); plot(this.tvec,X2);
            subplot(3,1,3); plot(this.tvec,X3);
            
            figure;
            subplot(3,1,1); plot(f,P1);
            subplot(3,1,2); plot(f,P2);
            subplot(3,1,3); plot(f,P3);
            
        end
        
        function [w,P] = get_fft(this,x,windowNumber)
            
            N = (length(x) - mod(length(x),windowNumber))/windowNumber;
            for i = 1:windowNumber
                [pxx(:,i), w] = periodogram(x(N*(i-1)+1:N*i),hamming(N),N,this.sfrq);
            end
            P = 10*log10((1/windowNumber)*(sum(pxx,2)));
            
            %             figure; plot(w,P); grid on; hold on;
            %             xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
            
        end
        
        function [] = testPronyApproach(this)
            
            this.sfrq = 500; % Change later this will change window size
            this.dt = 1/this.sfrq;
            t = 0:this.dt:3;
            
            [C_hat_c11] = get_PronyAnalysis(this,1,1,t);
            [C_hat_c12] = get_PronyAnalysis(this,1,2,t);
            [C_hat_c21] = get_PronyAnalysis(this,2,1,t);
            [C_hat_c22] = get_PronyAnalysis(this,2,2,t);
            
            syms s m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22
            
            M = [m11 m12; m21 m22];
            B = [b11 b12; b21 b22];
            K = [k11 k12; k21 k22];
            
            Z = M*s^2+B*s+K;
            C = inv(Z);
            
            [N,D] = numden(C);
            n(:,1,1) = flip(coeffs(N(1,1),s));
            d(:,1,1) = flip(coeffs(D(1,1),s));
            
            n(:,1,2) = flip(coeffs(N(1,2),s));
            d(:,1,2) = flip(coeffs(D(1,2),s));
            
            n(:,2,1) = flip(coeffs(N(2,1),s));
            d(:,2,1) = flip(coeffs(D(2,1),s));
            
            n(:,2,2) = flip(coeffs(N(2,2),s));
            d(:,2,2) = flip(coeffs(D(2,2),s));
            
            %             eqns = [n(:,1,1) == C_hat_c11.numerator{1}(end-2:end)';...
            %                      d(:,1,1) == C_hat_c11.Denominator{1}(:);...
            %                      n(:,1,2) == C_hat_c12.numerator{1}(end-2:end)';...
            %                      d(:,1,2) == C_hat_c12.Denominator{1}(:);...
            %                      n(:,2,1) == C_hat_c21.numerator{1}(end-2:end)';...
            %                      d(:,2,1) == C_hat_c21.Denominator{1}(:);...
            %                      n(:,2,2) == C_hat_c22.numerator{1}(end-2:end)';...
            %                      d(:,2,2) == C_hat_c22.Denominator{1}(:)];
            
            eqns = [n(:,1,1)./d(1,1,1) == C_hat_c11.numerator{1}(end-2:end)';...
                n(:,1,2)./d(1,1,1) == C_hat_c12.numerator{1}(end-2:end)';...
                n(:,2,1)./d(1,1,1) == C_hat_c21.numerator{1}(end-2:end)';...
                n(:,2,2)./d(1,1,1) == C_hat_c22.numerator{1}(end-2:end)'];
            %
            %             eqns = [d(2:end,1,1)/d(1,1,1) == C_hat_c11.Denominator{1}(2:end)';...
            %                     d(2:end,1,2)/d(1,1,2) == C_hat_c12.Denominator{1}(2:end)';...
            %                     d(2:end,1,2)/d(1,2,2) == C_hat_c22.Denominator{1}(2:end)'];
            
            vars = [m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22];
            
            [m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22] = solve(eqns,vars.');
            
            M = double([m11 m12; m21 m22]);
            B = double([b11 b12; b21 b22]);
            K = double([k11 k12; k21 k22]);
            
        end
        
        function [C_hat_c] = get_PronyAnalysis(this,dex1,dex2,t)
            
            s = tf('s');
            M = [1.7099, -0.2566; -0.2566, 2.1775];
            B = [5.2510, -1.0215; -1.0215, 39.0782];
            K = [105.0196, -20.4292; -20.4292, 781.5645];
            
            Z = M*s^2+B*s+K;
            % Z = M(2,2)*s^2+B(2,2)*s+K(2,2);
            C_total = inv(Z);
            C_c = C_total(dex1,dex2);
            H_c = impulse(C_c,t);
            C_d = c2d(C_c,this.dt,'match');
            H_d = impulse(C_d,t);
            % H_d2 = impz(C_d.numerator{1},C_d.Denominator{1},length(t));
            % figure;plot(t,H_d1,t,H_d2);
            
            %             figure; plot(t, H_c,'-',t,H_d,'o');
            
            [b,a] = prony(H_d,4,4);
            C_hat_d = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
            C_hat_c = d2c(C_hat_d,'matched');
            H_hat_d = impulse(C_hat_d,t);
            H_hat_c = impulse(C_hat_c,t);
            
        end
        
        function [K_ne_hat, h_model_ne11, h_model_ne12, h_model_ne21, h_model_ne22, VAFirf] = get_MIMO_Prony(this,iDex, h_hat_ne11, h_hat_ne12, h_hat_ne21, h_hat_ne22, M2, M2_short)
           
            %             figure;
            
            orderProny = 4;
            for i = iDex
                if(exist('M2_short','var'))
                    [b,a] = prony(h_hat_ne11(1:M2_short,i), orderProny, orderProny);
                else
                   [b,a] = prony(h_hat_ne11(:,i), orderProny, orderProny);
                end
                C_hat_d11 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
                C_hat_c11 = d2c(C_hat_d11,'tustin');
                h_model_ne11(:,i) = impulse(C_hat_c11,0:1/this.sfrq:M2/this.sfrq);
                
                VAFirf(1,i) = this.get_VAF(h_model_ne11(:,i),h_hat_ne11(:,i));
                
%                 subplot(2,2,1); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne11(:,i),0:1/this.sfrq:M2/this.sfrq,h_model_ne11(:,i));
                                
                [b,a] = prony(h_hat_ne12(:,i), orderProny, orderProny);
                C_hat_d12 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
                C_hat_c12 = d2c(C_hat_d12,'tustin');
                h_model_ne12(:,i) = impulse(C_hat_c12,0:1/this.sfrq:M2/this.sfrq);
                
                VAFirf(2,i) = this.get_VAF(h_model_ne12(:,i),h_hat_ne12(:,i));
                
%                 h_model_ne21(:,i) = h_model_ne12(:,i);
%                 VAFirf(3,i) = VAFirf(2,i);
%                 subplot(2,2,2); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne12(:,i),0:1/this.sfrq:M2/this.sfrq,h_model_ne12(:,i));
                
                [b,a] = prony(h_hat_ne21(:,i),orderProny,orderProny);
                C_hat_d21 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
                C_hat_c21 = d2c(C_hat_d21,'tustin');
                H_hat_c21 = impulse(C_hat_c21,0:1/this.sfrq:M2/this.sfrq);
                h_model_ne21(:,i) = impulse(C_hat_c21,0:1/this.sfrq:M2/this.sfrq);
                
                VAFirf(3,i) = this.get_VAF(h_model_ne21(:,i),h_hat_ne21(:,i));

                %            subplot(2,2,3);plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne21(:,i),0:1/this.sfrq:M2/this.sfrq,H_hat_c21);
                
                [b,a] = prony(h_hat_ne22(:,i), orderProny, orderProny);
                C_hat_d22 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
                C_hat_c22 = d2c(C_hat_d22,'tustin');
                h_model_ne22(:,i) = impulse(C_hat_c22,0:1/this.sfrq:M2/this.sfrq);
                
                VAFirf(4,i) = this.get_VAF(h_model_ne22(:,i),h_hat_ne22(:,i));

                
%                 subplot(2,2,4); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne22(:,i),0:1/this.sfrq:M2/this.sfrq,h_model_ne22(:,i));
                
                
                Z_hat = inv([C_hat_c11,C_hat_c12;C_hat_c21,C_hat_c22]);
                
                
%                 figure; pzplot(Z_hat);
%                 figure; bode(Z_hat); hold on;
                %             bode(Z); xlim([10^-1 10^2]);
                
%                 [K_tmp] = bode(Z_hat,1); % Get gain at 1 Hz
                
                [num,den] = tfdata(Z_hat);
                K_ne_hat(:,:,i) = [num{1,1}(end)/den{1,1}(end),num{1,2}(end)/den{1,2}(end);...
                                   num{2,1}(end)/den{2,1}(end),num{2,2}(end)/den{2,2}(end)];
                
            end
            
            function  [dexResponse]  = get_dexSettlingTime(this,h_hat,dex,threshold)
                % Find index when h_hat falls below a threshold and stays
                % there.
%                 for ii = dex-25:dex+25
                    h_hat_abs = abs(h_hat(:,dex));
                    [maxVal,dexMax] = max(h_hat_abs);
                    tmp = find(h_hat_abs > maxVal*threshold,1,'last');
                    if(~isempty(tmp))
                        dexEnd = tmp;
                    else
                        dexEnd = length(h_hat_abs);
                    end
%                 end
                dexResponse = 1:dexEnd;
                
            end
            
%              % Smooth over after fitting
%             for i = iDex
%                 
%                 K_ne_hat(1,1,:) = smooth(K_ne_hat(1,1,:));
%                 K_ne_hat(1,2,:) = smooth(K_ne_hat(1,2,:));
%                 K_ne_hat(2,1,:) = smooth(K_ne_hat(2,1,:));
%                 K_ne_hat(2,2,:) = smooth(K_ne_hat(2,2,:));
% 
%                 VAFirf(1,:) = smooth(VAFirf(1,:))';
%                 VAFirf(2,:) = smooth(VAFirf(2,:))';
%                 VAFirf(3,:) = smooth(VAFirf(3,:))';
%                 VAFirf(4,:) = smooth(VAFirf(4,:))';
%                 
%             end
            
        end
        
        function [z_r_1,z_r_2] = get_noise_subtract_x0(this,y_r_1, y_r_2, u_r_1, u_r_2, f_r_1, f_r_2)
            
            
           %% Add noise 
            n_r_1 = 1e-5*randn(size(y_r_1)); 
%             snr(y_r_ne1,n_r_ne1)
            z_r_1 = y_r_1 + n_r_1;
            
            % Add noise tan
            n_r_2 = 1e-5*randn(size(y_r_2));
%             snr(y_r_ne2,n_r_ne2)
            z_r_2 = y_r_2 + n_r_2;

            %% Remove x0
                        
            R = size(z_r_1,1);         
            
            z_r_1_mean = mean(z_r_1);
            z_r_2_mean = mean(z_r_2);

            figure; title('raw');
            subplot(2,1,1); plot(this.tvec,z_r_1',this.tvec,z_r_1_mean,'-k','linewidth',2.5); xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_1 (m)'); set(gca,'fontsize',18); grid on; set(gca,'fontsize',18); grid on;
            subplot(2,1,2); plot(this.tvec,z_r_2',this.tvec,z_r_2_mean,'-k','linewidth',2.5); xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_2 (m)');
            set(gca,'fontsize',18); grid on;
   
            for i = 1:R
                z_r_1(i,:) = detrend(z_r_1(i,:) - z_r_1_mean);
                z_r_2(i,:) = detrend(z_r_2(i,:) - z_r_2_mean);
            end
                        
            figure; title('mean subtracted');
            subplot(2,1,1); plot(this.tvec,z_r_1,'linewidth',2.5); hold on; xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_1 (m)'); set(gca,'fontsize',18); %ylim([-1e-4 1e-4]);
            subplot(2,1,2); plot(this.tvec,z_r_2,'linewidth',2.5); hold on; xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_2 (m)');
            set(gca,'fontsize',18); grid on; %ylim([-1e-4 1e-4]);
 
            
        end
        
        function [dexFit] = get_dexFit(this,h_hat)
           
            dex = find(abs(h_hat)>max(abs(h_hat))*0.1);
            
%             figure; plot(h_hat); hold on; plot(dex,h_hat(dex),'o');

            dexFit = 1:dex(end);
            
        end
        
        function [K_filt] = get_filtK(this,K_ne_hat,iDex,filterLength)
           
            if(mod(filterLength,2))
                error('The parameter filterLength must be even.');
            end
            
            for i = iDex
                if(i<length(K_ne_hat)-filterLength)
                    K_filt(1,1,i) = mean(K_ne_hat(1,1,i-filterLength/2:i+filterLength/2));
                    K_filt(1,2,i) = mean(K_ne_hat(1,2,i-filterLength/2:i+filterLength/2));
                    K_filt(2,1,i) = mean(K_ne_hat(2,1,i-filterLength/2:i+filterLength/2));
                    K_filt(2,2,i) = mean(K_ne_hat(2,2,i-filterLength/2:i+filterLength/2));
                else
                    K_filt(1,1,i) = K_ne_hat(1,1,i);
                    K_filt(1,2,i) = K_ne_hat(1,2,i);
                    K_filt(2,1,i) = K_ne_hat(2,1,i);
                    K_filt(2,2,i) = K_ne_hat(2,2,i);
                end
            end
                
        end
        
        function [] = get_Z_spectral(this,z_r_ne11, z_r_ne22, u_r_ne11, u_r_ne22)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % SS Dynamic Ankle Impedance Identification (Relaxed)
            % Impedance Identification in principal axis directions (DP and IE)
            
            
            % Load Anklebot Impedance Model: Measured separately without a human subject
            % load('/Users/jhermus/Documents/School/MIT/Research/Limb_Impedance/PreliminaryData/Stocastic/AnklebotImpedanceModel.mat');
            
            
            %% Variables
            nfft=100;
            nFreq = nfft/2+1;
            Hz = this.sfrq;
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
            dex = 1:this.sfrq*2;
            for i = 1:size(u_r_ne11,1)
                X1 = [X1; u_r_ne11(i,dex)'];
                X2 = [X2; u_r_ne22(i,dex)'];
                Y1 = [Y1; z_r_ne11(i,dex)'];
                Y2 = [Y2; z_r_ne22(i,dex)'];
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
                Z11_s_phi(j,1) = 180/pi*this.unwrap2(angle(Z11_s(j,1)),unwrapThreshold,'up');
                Z22_s_phi(j,1) = 180/pi*this.unwrap2(angle(Z22_s(j,1)),unwrapThreshold,'up');
                Z12_s_phi(j,1) = 180/pi*this.unwrap2(angle(Z12_s(j,1)),unwrapThreshold,'up');
                Z21_s_phi(j,1) = 180/pi*this.unwrap2(angle(Z21_s(j,1)),unwrapThreshold,'up');
                
                Z11_l_mag(j,1) = abs(Z11_l(j,1));
                Z22_l_mag(j,1) = abs(Z22_l(j,1));
                Z12_l_mag(j,1) = abs(Z12_l(j,1));
                Z21_l_mag(j,1) = abs(Z21_l(j,1));
                Z11_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z11_l(j,1)),unwrapThreshold,'up');
                Z22_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z22_l(j,1)),unwrapThreshold,'up');
                Z12_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z12_l(j,1)),unwrapThreshold,'up');
                Z21_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z21_l(j,1)),unwrapThreshold,'up');
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
            plot(TF_freq,Z11_l_mag(:,1),'LineWidth',2); grid on; box on;
            axis([xLowerLim xUpperLim yLowerLim11 yUpperLim11]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
            ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
            title('Z11','fontWeight','bold','fontSize',16);
            
            ax2 = subplot(2,2,2,'XScale','log','YScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z22_l_mag(:,1),'LineWidth',2); grid on; box on;
            axis([xLowerLim xUpperLim yLowerLim22 yUpperLim22]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
            ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
            title('Z22','fontWeight','bold','fontSize',16);
            
            ax3 = subplot(2,2,3,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z11_l_phi(:,1),'LineWidth',2); grid on; box on;
            axis([xLowerLim xUpperLim 0 180]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
            
            ax4 = subplot(2,2,4,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,Z22_l_phi(:,1),'LineWidth',2); grid on; box on;
            axis([xLowerLim xUpperLim 0 180]);
            xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
            
            linkaxes([ax1,ax2,ax3,ax4],'x');
            
            %% Partial Coherence Plot
            xLowerLim = 0.5;
            figure(2); hold on;
            set(gcf,'Color',[1,1,1]);
            
            ax1 = subplot(2,2,1,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC11_s(:,1),'LineWidth',2);
            grid on;box on; ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14);
            title('Y11 PC','fontWeight','bold','fontSize',16);
            
            ax2 = subplot(2,2,2,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC12_s(:,1),'LineWidth',2); hold off;
            grid on; box on; ylim([0 1]);
            axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14);
            title('Y12 PC','fontWeight','bold','fontSize',16);
            
            ax3 = subplot(2,2,3,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC21_s(:,1),'LineWidth',2); hold off;
            grid on;box on;ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14);
            title('Y21 PC','fontWeight','bold','fontSize',16);
            
            ax4 = subplot(2,2,4,'XScale','log');
            set(gca,'fontWeight','bold','fontSize',12); hold on;
            plot(TF_freq,PC22_s(:,1),'LineWidth',2); hold off;
            grid on;box on;ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
            xlabel('Hz','fontWeight','bold','fontSize',14);
            title('Y22 PC','fontWeight','bold','fontSize',16);
            
            linkaxes([ax1,ax2,ax3,ax4],'x');
            
%             %% Subtract Robot
%             xLowerLim = 0.5;
%             xUpperLim = 45.0;
%             yLowerLim_sub = 5e+2;
%             yUpperLim_sub = 1e+4;
%             
%             load('/Users/jhermus/Documents/School/MIT/Research/Schwartz_Collaboration/BallisticreleaseAnalysis/matlab/data/model_robot.mat');
%             model_mag = interp1(model_robot(:,1),model_robot(:,2),TF_freq,'linear')';
%             model_phi = interp1(model_robot(:,1),model_robot(:,3),TF_freq,'linear')';
%             
%             figure(3); hold on;
%             set(gcf,'Color',[1,1,1]);
%             
%             ax2 = subplot(2,1,1,'XScale','log','YScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z22_l_mag(:,1)-model_mag,'LineWidth',2,'Color', col_vec);
%             grid on; box on;
%             axis([xLowerLim xUpperLim yLowerLim_sub yUpperLim_sub]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
%             ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
%             title('Z22','fontWeight','bold','fontSize',16);
%             
%             ax4 = subplot(2,1,2,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z22_l_phi(:,1)-model_phi,'LineWidth',2,'Color', col_vec); grid on; box on;
%             axis([xLowerLim xUpperLim -90 45 ]); yticks([-90: 45 :45])
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
%             ylabel('phase (deg)','fontWeight','bold','fontSize',14);
%             
%             linkaxes([ax2,ax4],'x');
            
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
        
        function [TF_freq,Z11_l_mag,Z12_l_mag,Z21_l_mag,Z22_l_mag,Z11_l_phi,Z12_l_phi,Z21_l_phi,Z22_l_phi,PC11_s,PC12_s,PC21_s,PC22_s] = get_Z_spectral_individual(this,binDex, z_r_ne11, z_r_ne22, u_r_ne11, u_r_ne22)

            
            %% Try combining all
            
            nfft=1000;
            nFreq = nfft/2+1;
            Hz = this.sfrq;
            R = size(u_r_ne11,1);

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

            TF_freq = Hz/2*linspace(0,1,nfft/2+1);
            pw = PoolWaitbar(R, ['Simulating System ',int2str(R),' times.']);
            parfor i = 1:R
                X1 = u_r_ne11(i,binDex)';
                X2 = u_r_ne22(i,binDex)';
                Y1 = z_r_ne11(i,binDex)';
                Y2 = z_r_ne22(i,binDex)';
                
                [F,Px1x1(:,i),Px1x2(:,i),Px1y1(:,i),Px1y2(:,i),Px2x1(:,i),Px2x2(:,i),Px2y1(:,i),Px2y2(:,i),Py1x1(:,i),Py1x2(:,i),Py1y1(:,i),Py1y2(:,i),Py2x1(:,i),Py2x2(:,i),Py2y1(:,i),Py2y2(:,i)] = get_specInd(this,nfft,nFreq,Hz,X1,X2,Y1,Y2);
                increment(pw); % wait bar display

            end
            
            Px1x1 = mean(Px1x1,2);
            Px1x2 = mean(Px1x2,2); 
            Px1y1 = mean(Px1y1,2);
            Px1y2 = mean(Px1y2,2);
            
            Px2x1 = mean(Px2x1,2);
            Px2x2 = mean(Px2x2,2);
            Px2y1 = mean(Px2y1,2);
            Px2y2 = mean(Px2y2,2);
            
            Py1x1 = mean(Py1x1,2);
            Py1x2 = mean(Py1x2,2);
            Py1y1 = mean(Py1y1,2);
            Py1y2 = mean(Py1y2,2);
            
            Py2x1 = mean(Py2x1,2);
            Py2x2 = mean(Py2x2,2);
            Py2y1 = mean(Py2y1,2);
            Py2y2 = mean(Py2y2,2);
            
            OCx1x2 = abs(Px1x2).^2./(Px1x1.*Px2x2);
            OCx2x1 = abs(Px2x1).^2./(Px1x1.*Px2x2);
            
            OCx2y1 = abs(Px2y1).^2./(Px2x2.*Py1y1);
            OCx1y2 = abs(Px1y2).^2./(Px1x1.*Py2y2);
            OCx1y1 = abs(Px1y1).^2./(Px1x1.*Py1y1);
            OCx2y2 = abs(Px2y2).^2./(Px2x2.*Py2y2);

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
                Z11_l_mag(j,1) = abs(Z11_l(j,1));
                Z22_l_mag(j,1) = abs(Z22_l(j,1));
                Z12_l_mag(j,1) = abs(Z12_l(j,1));
                Z21_l_mag(j,1) = abs(Z21_l(j,1));
                Z11_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z11_l(j,1)),unwrapThreshold,'up');
                Z22_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z22_l(j,1)),unwrapThreshold,'up');
                Z12_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z12_l(j,1)),unwrapThreshold,'up');
                Z21_l_phi(j,1) = 180/pi*this.unwrap2(angle(Z21_l(j,1)),unwrapThreshold,'up');
            end
            
            %% Impedance plot (Diagonal)
            xLowerLim = 0.5;
            xUpperLim = 50.0;
            yLowerLim11 = 2e+3;
            yUpperLim11 = 1e+6;
            yLowerLim22 = 2e+3;
            yUpperLim22 = 3e+4;
            
%             figure(1); hold on;
%             % Magnitude plot of ankle impedance
%             ax1 = subplot(2,2,1,'XScale','log','YScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_l_mag(:,1),'LineWidth',2); grid on; box on;
% %             axis([xLowerLim xUpperLim yLowerLim11 yUpperLim11]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
%             ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
%             title('Z11','fontWeight','bold','fontSize',16);
%             
%             ax2 = subplot(2,2,2,'XScale','log','YScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z22_l_mag(:,1),'LineWidth',2); grid on; box on;
% %             axis([xLowerLim xUpperLim yLowerLim22 yUpperLim22]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14);
%             ylabel('magnitude (abs)','fontWeight','bold','fontSize',14);
%             title('Z22','fontWeight','bold','fontSize',16);
%             
%             ax3 = subplot(2,2,3,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z11_l_phi(:,1),'LineWidth',2); grid on; box on;
%             axis([xLowerLim xUpperLim 0 180]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
%             
%             ax4 = subplot(2,2,4,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,Z22_l_phi(:,1),'LineWidth',2); grid on; box on;
%             axis([xLowerLim xUpperLim 0 180]);
%             xlabel('frequency(Hz)','fontWeight','bold','fontSize',14); ylabel('phase (deg)','fontWeight','bold','fontSize',14);
%             
%             linkaxes([ax1,ax2,ax3,ax4],'x');
            
            %% Partial Coherence Plot
%             xLowerLim = 0.5;
%             figure(2); hold on;
%             set(gcf,'Color',[1,1,1]);
%             
%             ax1 = subplot(2,2,1,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,PC11_s(:,1),'LineWidth',2);
%             grid on;box on; ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
%             xlabel('Hz','fontWeight','bold','fontSize',14);
%             title('Y11 PC','fontWeight','bold','fontSize',16);
%             
%             ax2 = subplot(2,2,2,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,PC12_s(:,1),'LineWidth',2); hold off;
%             grid on; box on; ylim([0 1]);
%             axis([xLowerLim xUpperLim 0 1]);
%             xlabel('Hz','fontWeight','bold','fontSize',14);
%             title('Y12 PC','fontWeight','bold','fontSize',16);
%             
%             ax3 = subplot(2,2,3,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,PC21_s(:,1),'LineWidth',2); hold off;
%             grid on;box on;ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
%             xlabel('Hz','fontWeight','bold','fontSize',14);
%             title('Y21 PC','fontWeight','bold','fontSize',16);
%             
%             ax4 = subplot(2,2,4,'XScale','log');
%             set(gca,'fontWeight','bold','fontSize',12); hold on;
%             plot(TF_freq,PC22_s(:,1),'LineWidth',2); hold off;
%             grid on;box on;ylim([0 1]); axis([xLowerLim xUpperLim 0 1]);
%             xlabel('Hz','fontWeight','bold','fontSize',14);
%             title('Y22 PC','fontWeight','bold','fontSize',16);
%             
%             linkaxes([ax1,ax2,ax3,ax4],'x');
   
            
            
        end
        
        function [F,Px1x1,Px1x2,Px1y1,Px1y2,Px2x1,Px2x2,Px2y1,Px2y2,Py1x1,Py1x2,Py1y1,Py1y2,Py2x1,Py2x2,Py2y1,Py2y2] = get_specInd(this,nfft,nFreq,Hz,X1,X2,Y1,Y2)
                           
%             figure;
%             subplot(2,1,1); plot(Y1);
%             subplot(2,1,2); plot(Y2);
            
            Y1 = detrend(Y1);
            Y2 = detrend(Y2);
            
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
     
            F = Fx1x1;
        end
        
        function [VAFirf] = get_VAF(this,h_model,h_hat)
           
             % VAF IRF Calculation
                N = length(h_model);
                VAFirf = zeros(N,1);
                VAFirf = 100*(1-std(h_hat-h_model)^2/std(h_hat)^2);
            
        end
        
        function val=unwrap2(this,angle,threshold,mode)
            if(strcmp(mode,'up'))
                if(angle<threshold)
                    val=angle+2*pi;
                else
                    val=angle;
                end
            elseif(strcmp(mode,'down'))
                if(angle>threshold)
                    val=angle-2*pi;
                else
                    val=angle;
                end
            end
        end

    end
end

% %% Kinematic Stiffness Sanity Check
% q_1 = 0;
% q_2 = 0;
% 
% l1a = test.l1a;
% l2a = test.l2a;
% 
% armPosVec_X = [zeros(length(q_1),1), l1a*cos(q_1),...
%     l1a*cos(q_1) + l2a*cos(q_1 + q_2)];
% 
% armPosVec_Y = [zeros(length(q_1),1), l1a*sin(q_1),...
%     l1a*sin(q_1) + l2a*sin(q_1 + q_2)];
%             
% figure; plot(armPosVec_X, armPosVec_Y,'-k',armPosVec_X, armPosVec_Y,'.k','markersize',30,'linewidth',2.5);      
% axis equal; grid on; xlim([-1 1]); ylim([-1 1]); 
% 
% parJ_parq1 = [-( l1a*cos(q_1) + l2a*cos(q_1 + q_2) ),  -l2a*cos(q_1 + q_2);...
%                     -l1a*sin(q_1) - l2a*sin(q_1 + q_2), -l2a*sin(q_1 + q_2)]
%                 
% parJ_parq2 = [- l2a*cos(q_1 + q_2),  -l2a*cos(q_1 + q_2);...
%                     - l2a*sin(q_1 + q_2), -l2a*sin(q_1 + q_2)]
%                 

