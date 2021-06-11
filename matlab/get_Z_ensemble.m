classdef get_Z_ensemble < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sfrq
        dt
        tvec
        K_hat
        time
    end
    
    methods
        function [this] = get_Z_ensemble(Data_pert_ensemble)
            dbstop if error
            
%             L = length(Data_pert_ensemble(1).time);
                        
%             z_r_1 = zeros(length(Data_pert_ensemble),L);
%             z_r_2 = zeros(length(Data_pert_ensemble),L);
%             u_r_1 = zeros(length(Data_pert_ensemble),L);
%             u_r_2 = zeros(length(Data_pert_ensemble),L);
%             
%             for i = 1:length(Data_pert_ensemble)
%               z_r_1(i,1:L) = Data_pert_ensemble(i).tp(1:L,1);
%               z_r_2(i,1:L) = Data_pert_ensemble(i).tp(1:L,2);
%               u_r_1(i,1:L) = Data_pert_ensemble(i).cf(1:L,1);
%               u_r_2(i,1:L) = Data_pert_ensemble(i).cf(1:L,2);
%                 
%             end

            z_r_1 = Data_pert_ensemble.z_r_1;
            z_r_2 = Data_pert_ensemble.z_r_2;
            u_r_1 = Data_pert_ensemble.u_r_1;
            u_r_2 = Data_pert_ensemble.u_r_2;
            
            time = Data_pert_ensemble(1).time_r(1,:);
            this.sfrq = 100;%500;
            this.dt = 1/this.sfrq;
            
            figure;
            subplot(2,1,1); plot(time, z_r_1); ylabel('x_r (m)');
            subplot(2,1,2); plot(time, z_r_2); ylabel('y_r (m)'); xlabel('Time (s)');
            
            z_r_1_mean = mean(z_r_1);
            z_r_2_mean = mean(z_r_2);
            
            z_r_1 = z_r_1 - z_r_1_mean;
            z_r_2 = z_r_2 - z_r_2_mean;
            
            figure;
            subplot(2,1,1); plot(time, z_r_1); ylabel('x_r (m)');
            subplot(2,1,2); plot(time, z_r_2); ylabel('y_r (m)'); xlabel('Time (s)');
            
            %% ID and fit tangential
            
            [R,N] = size(z_r_1);
            lagBuffer = 5;
            M1 = 0;
            M2 = 0.5*this.sfrq+lagBuffer;
            L = M2-M1+1;
            this.tvec = 0:this.dt:this.dt*(L-1);

            [H_hat_MA, MAwindow] = this.ensambleSysID_Matrix(z_r_1, z_r_2, u_r_1, u_r_2, N, L, R, M1, M2);
            h_hat_11 = H_hat_MA(:,:,1);
            h_hat_12 = H_hat_MA(:,:,2);
            h_hat_21 = H_hat_MA(:,:,3);
            h_hat_22 = H_hat_MA(:,:,4);
            
            iDex = M2+1+40:N+M1-40; %floor(linspace(M2+1+40,N+M1-40,25));
         
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
%             [K_hat, h_model_11, h_model_12, h_model_21, h_model_22, VAFirf] = get_MIMO_Prony(this,iDex, h_hat_11_filt, h_hat_12_filt, h_hat_21_filt, h_hat_22_filt, M2);
            
            %% Make fit plots 
%             K_hat(1,1,:) =  K_hat(1,1,:)-2500; % Subtract something
%             B_ne11 = B_ne11-40;
                        
            %% Check Impulse Reponse Function and their fits
            figure('position',[289 231 895 527]);
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,2);
            ax3 = subplot(2,2,3);
            ax4 = subplot(2,2,4);
            
            for i = iDex(1):10:iDex(end)
                axes(ax1);
                plot3(this.tvec,time(i)*ones(1,L),h_hat_11_filt(:,i),'.r'); hold on;
%                 plot3(this.tvec,time(i)*ones(1,L),h_model_11(:,i),'-b'); hold on;
%                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,1,1,i),'-k'); hold on;
                
                axes(ax2);
                plot3(this.tvec,time(i)*ones(1,L),h_hat_12_filt(:,i),'.r'); hold on;
%                 plot3(this.tvec,time(i)*ones(1,L),h_model_12(:,i),'-b'); hold on;
%                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,1,2,i),'-k'); hold on;

                
                axes(ax3);
                plot3(this.tvec,time(i)*ones(1,L),h_hat_21_filt(:,i),'.r'); hold on;
%                 plot3(this.tvec,time(i)*ones(1,L),h_model_21(:,i),'-b'); hold on;
%                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,2,1,i),'-k'); hold on;

                
                axes(ax4);
                plot3(this.tvec,time(i)*ones(1,L),h_hat_22_filt(:,i),'.r'); hold on;
%                 plot3(this.tvec,time(i)*ones(1,L),h_model_22(:,i),'-b'); hold on;
%                 plot3(this.tvec(1:L),this.tvec(i)*ones(1,L),h_true_model(:,2,2,i),'-k'); hold on;

            end
            fs = 15;
            axes(ax1); view(35,25); title(11); xlim([0 M2/this.sfrq]); %zlim([-0.011 0.011]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);
            
            axes(ax2); view(35,25); title(12); xlim([0 M2/this.sfrq]); %zlim([-0.011 0.011]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_t (m)'); grid on; set(gca,'fontsize',fs);
            
            axes(ax3); view(35,25); title(21); xlim([0 M2/this.sfrq]); %zlim([-0.011 0.011]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_n (m)'); grid on; set(gca,'fontsize',fs);
            
            axes(ax4); view(35,25); title(22); xlim([0 M2/this.sfrq]); %zlim([-0.1 0.1]); 
            xlabel('lag (s)'); ylabel('Time (s)'); zlabel('\Delta e_t (m)'); grid on; set(gca,'fontsize',fs);
            
            
%             figure; 
%             plot(time(iDex),VAFirf(1,iDex),'linewidth',2.5); hold on;
%             plot(time(iDex),VAFirf(2,iDex),'linewidth',2.5);
%             plot(time(iDex),VAFirf(3,iDex),'linewidth',2.5);
%             plot(time(iDex),VAFirf(4,iDex),'linewidth',2.5);
%             ylabel('VAF_{irf}'); xlabel('Time (s)'); set(gca, 'fontsize', 16); ylim([0 100]);
%             legend('11','12','21','22');         
            
% %             figure('Position',[440 274 560 524]);
% %             labelss = {'11','12','21','22'};
            
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
% %             p1=plot(time(iDex), squeeze(K_hat(1,1,iDex)),'-b','linewidth',2.5,'markersize',15); hold on;
% %             p2=plot(time(iDex), squeeze(K_hat(1,2,iDex)),'-r','linewidth',2.5,'markersize',15); hold on;
% %             p3=plot(time(iDex), squeeze(K_hat(2,1,iDex)),'-g','linewidth',2.5,'markersize',15); hold on;
% %             p4=plot(time(iDex), squeeze(K_hat(2,2,iDex)),'-k','linewidth',2.5,'markersize',15);
% %             grid on;  xlim([time(1) time(end)]); ylabel(' Stiffness (N/m)'); %ylim([0 3000]);
            
%             figure;
%             plot(this.tvec(iDex), VAFirf_ne11(iDex),'.-b','linewidth',2.5,'markersize',15); hold on;
%             plot(this.tvec(iDex), VAFirf_ne12(iDex),'.-r','linewidth',2.5,'markersize',15); hold on;
%             plot(this.tvec(iDex), VAFirf_ne21(iDex),'.-g','linewidth',2.5,'markersize',15); hold on;
%             plot(this.tvec(iDex), VAFirf_ne22(iDex),'.-k','linewidth',2.5,'markersize',15);
%             grid on;  xlim([0.2 this.tvec(end)]); ylabel(' VAF_{irf}'); ylim([0 100]);
%             xlabel('Time (s)'); set(gca,'fontsize',18);
            
%             this.K_hat = K_hat;
%             this.time = time(1:N+M1-40);
            disp('test');
            
        end
        
        
        % Need to remove y_r to runn options other than simData2D        
        function [z_r_norm,u_r_norm,y_r_norm,z_r_tan,u_r_tan,y_r_tan,f_r_norm,f_r_tan] = get_ensambles_from_simSystem2D(this,numLoop,turnPeriod,symType)
                        
            z_r_norm = [];
            u_r_norm = [];
            y_r_norm = [];
            f_r_norm = [];
            
            z_r_tan = [];
            u_r_tan = [];
            y_r_tan = [];
            f_r_tan = [];
            
            pw = PoolWaitbar(numLoop, ['Simulating System 2D ',int2str(numLoop),' times.']);
            
            parfor i = 1:numLoop
                
                test = simSystem2D(turnPeriod);
                test.simulateSystem(symType);
                
                [tmp_z_r_norm,tmp_u_r_norm,tmp_y_r_norm,tmp_z_r_tan,tmp_u_r_tan,tmp_y_r_tan,tmp_f_r_norm, tmp_f_r_tan] = test.get_ensambleInput(symType);
                
                z_r_norm = [z_r_norm; tmp_z_r_norm];
                u_r_norm = [u_r_norm; tmp_u_r_norm];
                y_r_norm = [y_r_norm; tmp_y_r_norm];
                f_r_norm = [f_r_norm; tmp_f_r_norm];
                
                
                z_r_tan = [z_r_tan; tmp_z_r_tan];
                u_r_tan = [u_r_tan; tmp_u_r_tan];
                y_r_tan = [y_r_tan; tmp_y_r_tan];
                f_r_tan = [f_r_tan; tmp_f_r_tan];
                
                % disp([int2str(i),'/',int2str(numLoop)]);
                increment(pw); % wait bar display
            end
            
            test = simSystem2D(turnPeriod);
            this.sfrq = test.sfrq;
            this.dt = test.dt;
            this.tvec = test.t_cycle;
            clear test
            
            R = size(z_r_norm,1);
            
            figure; title('raw');
            subplot(2,1,1); plot(this.tvec(1:end-2),z_r_norm',this.tvec(1:end-2),mean(z_r_norm),'-k','linewidth',2.5); xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_t norm');
            subplot(2,1,2); plot(this.tvec(1:end-2),z_r_tan',this.tvec(1:end-2),mean(z_r_tan),'-k','linewidth',2.5); xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_t tan');
            set(gca,'fontsize',18); grid on;
            
            z_r_norm_mean = mean(z_r_norm);
            z_r_tan_mean = mean(z_r_tan);
            
%             z_r_norm_tmp = z_r_norm;
%             z_r_tan_tmp = z_r_tan;

%             z_r_norm = z_r_norm_tmp;
%             z_r_tan = z_r_tan_tmp;

            for i = 1:R
                z_r_norm(i,:) = detrend(z_r_norm(i,:) - z_r_norm_mean);
                z_r_tan(i,:) = detrend(z_r_tan(i,:) - z_r_tan_mean);
            end
                        
            figure; title('mean subtracted');
            subplot(2,1,1); plot(this.tvec(1:end-2),z_r_norm,'linewidth',2.5); hold on; xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_t norn'); %ylim([-1e-4 1e-4]);
            subplot(2,1,2); plot(this.tvec(1:end-2),z_r_tan,'linewidth',2.5); hold on; xlim([0 this.tvec(end)]);
            xlabel('Time (s)'); ylabel('z_t tan');
            set(gca,'fontsize',18); grid on; %ylim([-1e-4 1e-4]);
            
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
            MAwindow = 0.06*this.sfrq;  % Moving average window size: 40ms
            MAstep = 0.06*this.sfrq;
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
        
        function [PHI] = get_PHI(this,a,b,i,k,R,M2)
            
            %             w_a = 0.54 - 0.46*cos(2*pi*(i/M2));
            %             w_b = 0.54 - 0.46*cos(2*pi*(k/M2));
            
            PHI = (1/R)*sum(a(:,i).*b(:,i+k));
            
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
                indxs = find(p(:,1)<0);
                p(indxs,1) = p(indxs) + 2*pi;
                
                dexReset = find(abs(diff(p)) > 4);
                
                if(ensambleLength>10) % In slow case
                    dexReset = dexReset(find(diff(dexReset)>2000*2)); % elliminate multiple samples at same point
                end
                
                dexReset = dexReset(3:end-2); % Cut first 2 and last cycle
                dexRange = [dexReset(1):dexReset(end)];
                cycNum_current = size(dexReset,2)-1;
                
                X_norm = sqrt(X_m.^2+Y_m.^2);
                X_tan = p*0.1;
                
                Pret_norm =  Pret_X_N.*cos(p) + Pret_Y_N.*sin(p);
                Pret_tan =  -Pret_X_N.*sin(p) + Pret_Y_N.*cos(p);
                
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
           
        end
        
        function [z1_r, z2_r, u1_r, u2_r, N, L, R, M1, M2] = get_lowerSamplingRate(this,ensambleLength, z1_r, z2_r, u1_r, u2_r, N, L, R, M1, M2)
            
            this.sfrq = 200; % Change later this will change window size
            this.dt = 1/this.sfrq;
            tvec_old = this.tvec;
            this.tvec = 0:this.dt:ensambleLength+this.dt;
            N = length(this.tvec); % Check N is even
            lagBuffer = 5;
            M1 = 0;
            M2 = 0.3*this.sfrq+lagBuffer;
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
        
        function [K_ne_hat, h_model_ne11, h_model_ne12, h_model_ne21, h_model_ne22, VAFirf] = get_MIMO_Prony(this,iDex, h_hat_ne11, h_hat_ne12, h_hat_ne21, h_hat_ne22, M2)
           
            h_model_ne11 = zeros(size(h_hat_ne11));
            h_model_ne12 = zeros(size(h_hat_ne12));
            h_model_ne21 = zeros(size(h_hat_ne21));
            h_model_ne22 = zeros(size(h_hat_ne22));
            VAFirf = zeros(4,length(iDex));
            
            %             figure;
            
            orderProny = 4;
            for i = iDex
                % If second peaks is not half as large as first use window
                % default. If shorter use 0.3 s 
                [PKS,LOCS] = findpeaks(abs(h_hat_ne11(:,i)));
                if(PKS(2) < 0.5*PKS(1))
                    dex = 1:this.sfrq*0.3;
                else
                    dex = 1:length(this.tvec);
                end
                
                [b,a] = prony(h_hat_ne11(dex,i), orderProny, orderProny);
                C_hat_d11 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
%                 z = zero(C_hat_d11);
%                 if(~sum(z<0 & imag(z)==0)>0)
                    C_hat_c11 = d2c(C_hat_d11,'tusting');
                    h_model_ne11(:,i) = impulse(C_hat_c11,0:1/this.sfrq:M2/this.sfrq);
                
                    VAFirf(1,i) = this.get_VAF(h_model_ne11(dex,i),h_hat_ne11(dex,i));
%                 end
%                 subplot(2,2,1); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne11(:,i),0:1/this.sfrq:M2/this.sfrq,h_model_ne11(:,i));
                                
                [b,a] = prony(h_hat_ne12(dex,i), orderProny, orderProny);
                C_hat_d12 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
                 
%                 z = zero(C_hat_d12);
%                 if(~sum(z<0 & imag(z)==0)>0)
                    C_hat_c12 = d2c(C_hat_d12,'tustin');
                    h_model_ne12(:,i) = impulse(C_hat_c12,0:1/this.sfrq:M2/this.sfrq);
                    
                    VAFirf(2,i) = this.get_VAF(h_model_ne12(dex,i),h_hat_ne12(dex,i));
                
                
%                     h_model_ne21(:,i) = h_model_ne12(:,i);
%                     VAFirf(3,i) = VAFirf(2,i);
%                 end
%                 subplot(2,2,2); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne12(:,i),0:1/this.sfrq:M2/this.sfrq,h_model_ne12(:,i));
                
                [b,a] = prony(h_hat_ne21(dex,i),orderProny,orderProny);
                C_hat_d21 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
                C_hat_c21 = d2c(C_hat_d21,'tustin');
                h_model_ne21(:,i) = impulse(C_hat_c21,0:1/this.sfrq:M2/this.sfrq);
                                    
                VAFirf(3,i) = this.get_VAF(h_model_ne21(dex,i),h_hat_ne21(dex,i));

%                            subplot(2,2,3);plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne21_filt(:,i),0:1/this.sfrq:M2/this.sfrq,H_hat_c21);
                
                [b,a] = prony(h_hat_ne22(dex,i), orderProny, orderProny);
                C_hat_d22 = tf(b,this.sfrq*a,this.dt); % Why is this off by a factor of the sampling frequency??
%                 z = zero(C_hat_d22);
%                 if(~sum(z<0 & imag(z)==0)>0)
                    C_hat_c22 = d2c(C_hat_d22,'tustin');
                    h_model_ne22(:,i) = impulse(C_hat_c22,0:1/this.sfrq:M2/this.sfrq);
                    
                    VAFirf(4,i) = this.get_VAF(h_model_ne22(dex,i),h_hat_ne22(dex,i));
%                 end
                
%                 subplot(2,2,4); plot(0:1/this.sfrq:M2/this.sfrq,h_hat_ne22(:,i),0:1/this.sfrq:M2/this.sfrq,h_model_ne22(:,i));
                
                
                Z_hat = inv([C_hat_c11,C_hat_c12;C_hat_c12,C_hat_c22]);
%                 [K,B,M] = get_coeffs(this,C_hat_c11,C_hat_c12,C_hat_c12,C_hat_c22);

%                 K_ne_hat(:,:,i) = K;
%                 B_ne_hat(:,:,i) = B;
%                 M_ne_hat(:,:,i) = M;

%                 figure; pzplot(Z_hat);
%                 figure; bode(Z_hat); hold on;
                %             bode(Z); xlim([10^-1 10^2]);
                
%                 [K_tmp] = bode(Z_hat,1); % Get gain at 1 Hz

                % Take only low frequency approach
                [num,den] = tfdata(Z_hat);
                K_ne_hat(:,:,i) = [num{1,1}(end)/den{1,1}(end),num{1,2}(end)/den{1,2}(end);...
                                   num{2,1}(end)/den{2,1}(end),num{2,2}(end)/den{2,2}(end)];
                
            end
            
        end     
        
        function [K,B,M] = get_coeffs(this,C_hat_c11,C_hat_c12,C_hat_c21,C_hat_c22)
                   
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
%                 d(:,1,1) == C_hat_c11.Denominator{1}(:);...
%                 n(:,1,2) == C_hat_c12.numerator{1}(end-2:end)';...
%                 d(:,1,2) == C_hat_c12.Denominator{1}(:);...
%                 n(:,2,1) == C_hat_c21.numerator{1}(end-2:end)';...
%                 d(:,2,1) == C_hat_c21.Denominator{1}(:);...
%                 n(:,2,2) == C_hat_c22.numerator{1}(end-2:end)';...
%                 d(:,2,2) == C_hat_c22.Denominator{1}(:)];
            
            eqns = [n(:,1,1)./d(1,1,1) == C_hat_c11.numerator{1}(end-2:end)';...
                n(:,1,2)./d(1,1,1) == C_hat_c12.numerator{1}(end-2:end)';...
                n(:,2,1)./d(1,1,1) == C_hat_c21.numerator{1}(end-2:end)';...
                n(:,2,2)./d(1,1,1) == C_hat_c22.numerator{1}(end-2:end)'];
            
%             eqns = [d(2:end,1,1)/d(1,1,1) == C_hat_c11.Denominator{1}(2:end)';...
%                 d(2:end,1,2)/d(1,1,2) == C_hat_c12.Denominator{1}(2:end)';...
%                 d(2:end,1,2)/d(1,2,2) == C_hat_c22.Denominator{1}(2:end)'];
            
            vars = [m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22];
            
            [m11 m12 m21 m22 b11 b12 b21 b22 k11 k12 k21 k22] = solve(eqns,vars.');
            
            M = double([m11 m12; m21 m22]);
            B = double([b11 b12; b21 b22]);
            K = double([k11 k12; k21 k22]);
            
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

