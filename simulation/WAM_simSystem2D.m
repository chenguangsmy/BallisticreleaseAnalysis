
classdef WAM_simSystem2D < handle
    %Data intial processing for InMotion2 ZFT comptutations
    % Filename:	trial.m
    % Author:  James Hermus
    % Date:     24 January 2019
    % Description:	When initalized this class imports subject data from
    %               exported from the FPGA and computes the ZFT.
    
    properties
        
        % Subject Spesific Parameters
        subjectNumber % subject number in study
        l1a % Length of upper arm (m)
        l2a % Length of forearm (from elbow to center of fist) (m)
        ch  % Length from wrist to center of fist (m)
        wt  % Subject weight (kg)
        c1a % Distance from shoulder to uppre arm center of gravity (m)
        c2a % Distance from elbow to forearm center of gravity (m)
        I1a % Interia of upper arm (kg-m^2)
        I2a % Inertia of forearm  (including hand) (kg-m^2)
        m1a % Mass of upper arm (m)
        m2a % Mass of forearm (m)
        
        % Trial Spesific Parameters
        lc    % Crank Radius (m)
        cc    % Distance from crank pivot to crank center of gravity (outboard mass removed) (m)
        Ic    % Inertia of the crank (kg)
        mc    % Mass of crank, measured minus estimated outboard mass (kg)
        d1    % Horiz distance from crank center to shoulder center (m) (right when facing robot is positive)
        d2    % Vertical distance from crank center to shoulder center (m) (postive is toward the robot when facing it)
        rob_d1 % x distance from center of crank to robot base (m) (toward robot is positive)
        rob_d2 % y distance from center of crank to robot base (m) (right when facting robot is postive)
        handSide % boolin which indicates hand used: 0 indicates right, 1 indicates left
        
        % Assumed Human dynamics
        K % Stiffness matrix (N-m/rad)
        B % Damping matrix (N-m-s/rad)
        
        % InMotion2 Parameters
        % Length
        rob_l1
        rob_l2
        rob_l3
        rob_l4
        % Mass for the corresponding links
        rob_m1
        rob_m2
        rob_m3
        rob_m4
        % Distance between center of mass and joint closest to the robot link
        rob_r1
        rob_r2
        rob_r3
        rob_r4
        % Inertia at the end of the link closest to the motor.
        rob_I1
        rob_I2
        rob_I3
        rob_I4
        
        sfrq % Sampling frequency
        dt
        turnPeriod
        N_tot
        movDuration
        
        % Simulation Results
        t % Time (s)
        X % InMotion2 handel position: [Nx2] (m) in cartesian crank coordinates
        V % InMotion2 handel velocity: [Nx2] (m/s) in cartesian crank corrdinates
        A % InMotion2 handel acceleration: [Nx2] (m/s^2) in cartesian corrdinates
        F % InMotion2 handel force: [Nx2] (N) in cartesian corrdinates
        
        % Computed Quantities
        thcp  % Crank angular postion (rad)
        thcv  % Crank angular velocity (rad/s)
        thca  % Crank angular acceleration (rad/s^2)
        F_rot % Force normal and tangental [Nx2] (m)
        q     % shoulder and elbow realitive joint angles [Nx2] (rad)
        q_dot % shoulder and elbow realitive joint velocities [Nx2] (rad)
        q_0   % ZFT postion in joint space [Nx2] (rad)
        q_0_dot
        X_0   % carteasion ZFT (m)
        X_0_dot
        F_pret % Force preturbation [Nx2]
        F_contact
        
        % Plotting params (add comments)
        armPosVec_X
        armPosVec_Y
        crank_x
        crank_y
        constraint_x
        constraint_y
        robPosVec_X
        robPosVec_Y
        pretScale
        
        %         % Ensamble outputs
        %         z_r
        %         u_r
        %         t_cycle
        %         N
        %         L
        %         R
        %         M1
        %         M2
        
    end
    
    methods
        function  [this] = WAM_simSystem2D(movDuration,pretScale)
            this.handSide = 0;
            this.subjectNumber = 1;
            this.lc = 0.1;
            this.sfrq = 100; % Change later this will change window size
            this.dt = 1/this.sfrq;
            this.movDuration = movDuration;
            this.t = 0:this.dt:this.movDuration*4 + this.dt;
            this.N_tot = length(this.t); % Check N is even
            
            if(exist('pretScale','var'))
                this.pretScale = pretScale;
            end
            
            get_subjectParam(this); % Defines subject paramters
            get_InMotionParam(this); % Define InMotion2 paramters
            %             tic;
            %             this.simulateSystem('line');
            %             toc;
            %             this.simulateSystem('point_singlePret');
            %             this.simulateSystem('circle');
            %             this.simulateSystem('circle_singlePret','nn');
            %             this.simulateSystem('line_SISO');
            
            %             [z_r_norm,u_r_norm,y_r_norm,z_r_tan,u_r_tan,y_r_tan,t_cycle,N,R] = this.get_ensambleInput('point');
            %             [z_r_norm,u_r_norm,y_r_norm,z_r_tan,u_r_tan,y_r_tan,t_cycle,N,R] = this.get_ensambleInput('circle'); % Dont have to spesify direction
            %
            %             this.save_zftAnimation('/Users/jhermus/Downloads/Test.mp4');
            %             get_InMotionData(this); % Imports the raw data and filters
            %             get_binnedMeasures(this) % Bins data for plotting
        end
        
        function [] = simulateSystem(this,symType)
            % Filename:	ZFT.m
            % Author:  James Hermus
            % Date:     16 July 2018
            % Description:	Computes the zero force trajectory from subject data
            % collected on the InMotion2with
            % Note: The computation uses realitive damping
            
            % ** For ode45 get everything in shoulder centered coordinates
            
            if(strcmp(symType,'line'))
                
                %% Simulate line move
                
                % Initial State
                % Define the initial state vector
                xx0 = zeros(4,1);
                [xx0(1), xx0(2)] = this.invKino(this.d1,this.d2);
                
                % Define ZFT
                [this.q_0, this.q_0_dot, this.X_0, this.X_0_dot] = this.get_ZFT(this.t);
                
                % Create bandwidthlimited 50 Hz gausian white noise (std before or after filt?)
                %                 this.F_pret = [this.get_instanceOrRandomPreturbation();...
                %                                this.get_instanceOrRandomPreturbation()]'; % (1/8).*
                
                this.F_pret = [this.get_StepPreturbation();...
                    this.get_StepPreturbation()]';
                
                % this.F_pret = this.get_StepPreturbation_2D'; % (1/8).*
                
                % Use ode45 to solve state equations
                [tvec, xvec] = ode45(@(t,y)this.getStates(t,y), this.t, xx0);
                
                % Use SemiImplicitEuler
                %                 [tvec, xvec] = this.SemiImplicitEuler(@(t,y)this.getStates(t,y),[this.t(1) this.t(end)], xx0);
                for i = 1:length(xvec)
                    [tmp,F_contact(:,i)] = this.getStates(tvec(i),xvec(i,:)');
                end
                
            elseif(strcmp(symType,'line_SISO'))
                
                %% Simulate line move SISO
                
                % Initial State
                % Define the initial state vector
                xx0 = zeros(4,1);
                [xx0(1), xx0(2)] = this.invKino(this.d1,this.d2);
                
                % Define ZFT
                [this.q_0, this.q_0_dot, this.X_0, this.X_0_dot] = this.get_ZFT(this.t);
                
                % Create bandwidthlimited 50 Hz gausian white noise (std before or after filt?)
                this.F_pret = [this.get_StepPreturbation();...
                    zeros(1,this.N_tot)]';
                %                 this.F_pret = [this.get_instanceOrRandomPreturbation();...
                %                                zeros(1,this.N_tot)]';
                
                % Use ode45 to solve state equations
                [tvec, xvec] = ode45(@(t,y)this.getStates(t,y), this.t, xx0);
                
                % Use SemiImplicitEuler
                %                 [tvec, xvec] = this.SemiImplicitEuler(@(t,y)this.getStates(t,y),[this.t(1) this.t(end)], xx0);
                for i = 1:length(xvec)
                    [tmp,F_contact(:,i)] = this.getStates(tvec(i),xvec(i,:)');
                end
                
                %                 figure; plot(tvec,this.F_pret(:,1));
                %                 ylim([-6 6]); xlim([0 4]); xlabel('Time (s)'); ylabel('Preturbation (N)');
                %                 set(gca,'fontsize',16);
                
            else
                error('spesify symType as point or circle');
            end
            
            %% Extra plotting things for simulate
            clear tvec
            % Extract states
            this.q = [xvec(:,1),xvec(:,2)];
            [this.X] = this.forwKino( this.q );
            [this.thcp] = get_crankPostion_Vec(this,this.X,'shoulder');
            this.F_contact = F_contact';
            
            % Create data structures for plotting
            [this.armPosVec_X, this.armPosVec_Y] = plotPosArm(this,this.X(:,1),this.X(:,2),'shoulder');
            [this.robPosVec_X, this.robPosVec_Y] = plotPosInMotion(this,this.X(:,1),this.X(:,2),'shoulder');
            this.crank_x = [this.d1*ones(size(this.X(:,1))), this.d1 + this.lc*cos(this.thcp)];
            this.crank_y = [this.d2*ones(size(this.X(:,1))), this.d2 + this.lc*sin(this.thcp)];
            this.constraint_x = this.d1 - this.lc*cos(0:0.01:2*pi)';
            this.constraint_y = this.d2 - this.lc*sin(0:0.01:2*pi)';
            
            % Point Constraint
            %               dex = 1;
            %                         figure;plot(this.armPosVec_X(1,:), this.armPosVec_Y(1,:), '.-b',...
            %                             this.robPosVec_X(1,:), this.robPosVec_Y(1,:), '.-k','linewidth',4,'markersize',40);hold on;
            %                                         this.crank_x(1,:),          this.crank_y(1,:),          '.-k',...
            %                         ylim([0 1.2]); xlim([-0.56 0.56]); grid on; axis equal;
            %                         ylabel('Y Distance (m)','fontsize',14);
            %                         xlabel('X Distance (m)','fontsize',14);
            %                         legend('Arm','InMotion','location','east');legend boxoff;
            %                         set(gca,'fontsize',18);
            
            % Circular Constraint
            %             dex = 3000:length(xvec);
            %             figure;plot(this.constraint_x,     this.constraint_y,     ':k',...
            %                 this.X_0(dex,1),this.X_0(dex,2),'-r',...
            %                 this.armPosVec_X(1,:), this.armPosVec_Y(1,:), '.-b',...
            %                 this.robPosVec_X(1,:), this.robPosVec_Y(1,:), '.-k','linewidth',4,'markersize',40);hold on;
            %             %                 this.crank_x(1,:),          this.crank_y(1,:),          '.-k',...
            %             ylim([0 1.2]); xlim([-0.56 0.56]); grid on; axis equal;
            %             ylabel('Y Distance (m)','fontsize',14);
            %             xlabel('X Distance (m)','fontsize',14);
            %             legend('Constraint','x_0','Arm','InMotion','location','east');legend boxoff;
            %             set(gca,'fontsize',18);
            %
            %             figure; dex = 3000:length(xvec);
            %             plot(this.X(dex,1),this.X(dex,2),'linewidth',2.5); hold on;
            %             plot(this.X_0(dex,1),this.X_0(dex,2),'r','linewidth',2.5); hold on;
            %             plot(this.d1+this.lc*sin(0:0.001:2*pi),this.d2+this.lc*cos(0:0.001:2*pi),'--k','linewidth',2.5);
            %             axis equal;
        end
        
        function get_subjectParam(this)
            
            if(this.subjectNumber == 1)
                
                % This function gets the parameters measured off the subject.
                % Subject Name: Bothabo Ngwenya
                % Age: 22
                % Height: 5'11"
                % Weight: 121 lbs
                % Right handed
                
                %% HUMAN:
                this.l1a = 0.31;         % Length of upper arm
                this.l2a = 0.36;         % Length of forearm (from elbow to center of fist)
                this.ch = 0.085;        % From wrist to center of fist
                this.wt = 121*0.454;    % Weight in kg
                this.d1 = 0;            % Horiz distance from shoulder to crank center
                this.d2 = 0.43; %0.53;         % Vert distance from shoulder to crank center
                
                %% CRANK: (use Joe crank as of now)
                this.cc = (0.7005*(2.55*0.0254) - 0.1418*this.lc)/(0.7005-0.1418);
                % Distance from crank pivot to crank cg. Measured, outboard mass removed
                this.mc = 0.7005 - 0.1418; % mass of crank, measured minus estimated outboard
                
                % HUGE MISTAKE, CORRECTED 22 JUNE 1998
                % Changed back by James Hermus 21 October, 2018 for this application we
                % desire the crank inertia with respect to the pivot not the center of mass
                this.Ic = 0.0221; % Inertia of crank, measured (outboard subtracted).
                
                % Note from thesis Joe's this number comes from (Ic = It - mo*lc^2)
                
                % THE ABOVE INERTIA WAS MEASURED WRT THE CRANK PIVOT!!!!
                % Ic = 3.716e-3 - mc*cc^2; % WRT center of mass
                % Inertia Parameters, from p. 70 of Justin's Thesis:
                
                ka = 0.322*this.l1a;
                mf2 = 0.016*this.wt;
                cf2 = 0.430*(this.l2a-this.ch); kf2 = 0.303*(this.l2a-this.ch);
                mh = 0.006*this.wt + 0.5063; % Hand + (handle + outboard FT mass)
                
                this.c1a = 0.436*this.l1a; % Distance from shoulder to upper arm cg
                this.c2a = (mf2*cf2 + mh*this.l2a)/(mf2+mh); % Distance from elbow to forearm cg
                
                this.m1a = 0.028*this.wt; % mass of upper arm
                this.I1a = this.m1a*ka^2; % Inertia of upper arm
                
                If2 = mf2*kf2^2; % Inertia of forearm (minus hand)
                this.I2a = If2 + mh*this.l2a^2; % Inertia of forearm (including hand)
                this.m2a = mf2+mh; % mass of forearm
                
                % Define joint torques using equlibrium trajectory
                K11 = 29.5;
                K12 = 14.3;
                K22 = 39.3;
                this.K = 2*[K11, K12; K12, K22]; % (Slow stiffness)
                lambda = 0.1; % (Slow Damping)
                this.B = lambda*this.K;
                
            else
                error('Enter valid subject number.');
            end
        end
        
        function get_InMotionParam(this)
            %% InMotion2 Parameters
            % Assume crank center position (CHECK EVERY TIME)
            this.rob_d1 = 0;
            this.rob_d2 = 0.65;
            
            % Length
            this.rob_l1 = 0.4064;
            this.rob_l2 = 0.5144;
            this.rob_l3 = 0.4064;
            this.rob_l4 = 0.1555;
            % Mass for the corresponding links
            this.rob_m1 = 0.756;
            this.rob_m2 = 1.964;
            this.rob_m3 = 0.756;
            this.rob_m4 = 0.378;
            % Distance between center of mass and joint closest to the robot link
            this.rob_r1 = 0.2032;
            this.rob_r2 = 0.3906;
            this.rob_r3 = 0.2032;
            this.rob_r4 = 0.0775;
            % Inertia at the end of the link closest to the motor.
            this.rob_I1 = 0.0416;
            this.rob_I2 = 0.3484;
            this.rob_I3 = 0.0416;
            this.rob_I4 = 0.0030;
        end
        
        function get_InMotionData(this)
            % Filename:	sortDataInMotion.m
            % Author:  James Hermus
            % Date:     18 July 2018
            % Description:
            % -Imports the InMotion2 Data
            % - Computes spline differentiates the postion to get vel and acc
            % - Butterworth filters the force
            % - Ouputs the varibles requried by the ZFT function
            
            % dataTypee spesifies 'raw' or 'wahba' filter (requried to get A --> ZFT)
            
            load(this.dataFile); % imports .mat file
            
            % Convert exported formate to the analysis formate
            t_raw = t; % orignal sampled time for data
            clear t
            x_raw = X_m;
            y_raw = Y_m;
            Fx_raw = X_N;
            Fy_raw = Y_N;
            
            %% Interpolate samping rate to constant
            N = t_raw(end)/(1/this.sfrq);
            
            % account for special case were time is devisible by sfrq exactly
            if(mod(t_raw(end),(1/this.sfrq)) == 0)
                t_resample = 0:1/this.sfrq:t_raw(end); % Define time resampled
            else
                t_resample = 0:1/this.sfrq:(1/this.sfrq)*N;
            end
            
            % Interpolate
            x_resample = interp1(t_raw,x_raw,t_resample,'spline')'; % X postion
            y_resample = interp1(t_raw,y_raw,t_resample,'spline')'; % Y postion
            
            Fx_resample = interp1(t_raw,Fx_raw,t_resample,'spline')'; % X Force
            Fy_resample = interp1(t_raw,Fy_raw,t_resample,'spline')'; % Y Force
            
            % Check data interpolation
            % figure;
            % ax1 = subplot(2,1,1); plot(t_raw, x_raw,'.r', t_raw, x_resample,'or', t_raw, y_raw,'.b', t_resample, y_resample,'ob');
            % ylabel('Position (m)'); title('Resample to Constant Rate');
            % ax2 = subplot(2,1,2); plot(t_raw, Fx_raw,'.r', t_resample, Fx_resample,'or', t_raw, Fy_raw,'.b', t_resample, Fy_resample,'ob');
            % ylabel('Force (N)'); xlabel('Time (s)');
            % linkaxes([ax1,ax2],'x'); xlim([t_raw(round(N/2)), round(t_raw(round(N/2 + 25)))]);
            
            %% Filter Motion using Dohrmann Spline Smoothing Method
            
            % Spline diffrenetiation
            lo_bound = -20;
            hi_bound = 10;
            N_test = 35;
            
            % Look at wahba
            % Bx_min = wahba(x_resample, sfrq, lo_bound, hi_bound, N_test);
            % By_min = wahba(y_resample, sfrq, lo_bound, hi_bound, N_test);
            
            % Choose B
            B_dohrman = 1e-8;
            
            [x,dx,ddx,dddx,Vx] = dohrmann(x_resample, this.sfrq, B_dohrman); % Smooth and differentiate X position
            [y,dy,ddy,dddy,Vy] = dohrmann(y_resample, this.sfrq, B_dohrman); % Smooth and differentiate Y position
            
            % [f_raw,P1_raw] = MakefftPlot(sfrq,x_raw);
            % figure; plot(f_raw,P1_raw)
            
            % Check filter Postion
            %             figure;
            %             ax1 = subplot(3,1,1); hold on; plot(x_raw,'.'); plot(x,'linewidth',1.5); hold off; title('Forces');
            %             xlabel('Sample Number'); ylabel('X-Position (m)'); set(gca, 'fontsize',14);legend('Raw Signal','Filtered Signal');
            %
            %             ax2 = subplot(3,1,2); hold on; plot(diff(smooth(x_raw))*sfrq,'.'); plot(dx,'linewidth',1.5); hold off; title('Forces');
            %             xlabel('Sample Number'); ylabel('X-Velocity (m/s)'); set(gca, 'fontsize',14);legend('Raw Signal','Filtered Signal');
            %
            %             ax3 = subplot(3,1,3); hold on; plot(diff(smooth(diff(smooth(x_raw))))*sfrq^2,'.'); plot(ddx,'linewidth',1.5); hold off; title('Forces');
            %             xlabel('Sample Number'); ylabel('X- Acceleration (m/s^2)'); set(gca, 'fontsize',14); legend('Raw Signal','Filtered Signal');
            %             linkaxes([ax1, ax2, ax3],'x');
            
            %% Fitler Force using Butterworth
            cf = 2;
            [b,a] = butter(3,cf/(this.sfrq/2));
            Fx = filtfilt(b,a,Fx_resample);
            Fy = filtfilt(b,a,Fy_resample);
            
            % Check filter Force
            %             figure;
            %             ax1 = subplot(2,1,1); hold on; plot(Fx_resample); plot(Fx,'linewidth',1.5); hold off; title('Forces'); xlim([1 size(Fx_resample,1)]);
            %             xlabel('Sample Number'); ylabel('X-Force (N)'); set(gca, 'fontsize',14);legend('Raw Signal','Filtered Signal');
            %             ax2 = subplot(2,1,2); hold on; plot(Fy_resample); plot(Fy,'linewidth',1.5); hold off;
            %             xlabel('Sample Number'); ylabel('Y-Force (N)'); set(gca, 'fontsize',14); legend('Raw Signal','Filtered Signal');
            %             linkaxes([ax1, ax2],'x'); xlim([1 size(Fx_resample,1)]);
            
            %% Keep only full crank cycles
            % Compute position as positive angle (assume in crank centered coordinates)
            p = atan2(y, x);
            indxs = find(p(:,1)<0);
            p(indxs,1) = p(indxs) + 2*pi;
            
            dexReset = find(abs(diff(p)) > 4);
            dexRange = [dexReset(1):dexReset(end)];
            
            %% Define outputs
            this.t = t_resample(dexRange)'; % time vector
            this.X = [x(dexRange), y(dexRange)]; % Endpoint postion
            this.V = [dx(dexRange), dy(dexRange)]; % Endpoint velocity
            this.A = [ddx(dexRange), ddy(dexRange)]; % Endpoint acceleration
            this.F = [Fx(dexRange), Fy(dexRange)]; % Endpoint force
            
            % Define Angular Position, Velocity, and Acceleration
            this.thcp = p(dexRange);
            this.thcv = sqrt(sum(this.V.^2,2))/this.lc;
            this.thca = (-this.A(:,1).*sin(this.thcp) + this.A(:,2).*cos(this.thcp))/this.lc;
            
            % Compute Normal and Tangential Forces
            F_norm =  this.F(:,1).*cos(this.thcp) + this.F(:,2).*sin(this.thcp);
            F_tan =  -this.F(:,1).*sin(this.thcp) + this.F(:,2).*cos(this.thcp);
            this.F_rot = [F_norm,F_tan]; % Normal and Tangential Forces
            
            % q1p - Shoulder realative joint angle
            % q1v - Shoulder realative joint angular velocity
            % q1a - Shoulder realative joint angular acceleation
            
            % q2p - Elbow realative joint angle
            % q2v - Elbow realative joint angular velocity
            % q2a - Elbow realative joint angular acceleation
            
            [X_shoulder] = this.convert_cartRefFrame(this.X,0,1); % Convert from crank to shoulder coordinates
            [X_rob] = this.convert_cartRefFrame(this.X,0,2); % Convert from crank to shoulder coordinates
            
            [q1p, q2p, q1v, q2v ] = this.invKino(X_shoulder(:,1),X_shoulder(:,2), this.V(:,1),this.V(:,2));
            this.q = [q1p,q2p];
            this.q_dot = [q1v,q2v];
            
            % Same goes for q1dp but the "d" denotes desired or equlibirum trajectory
            
            % Check output
            % figure;
            % ax1 = subplot(4,1,1); plot(this.t,this.X); ylabel('Pos (m)');set(gca, 'fontsize', 16);
            % ax2 = subplot(4,1,2); plot(this.t,this.V); ylabel('Vel (m/s)');set(gca, 'fontsize', 16);
            % ax3 = subplot(4,1,3); plot(this.t,this.A); ylabel('Acc (m/s^2)');set(gca, 'fontsize', 16);
            % ax4 = subplot(4,1,4); plot(this.t,this.F); ylabel('Force (N)');
            % xlabel('Time (s)'); set(gca, 'fontsize', 16);linkaxes([ax1, ax2, ax3, ax4],'x');
            
            % Check Math with numerical differentiation
            % figure; subplot(2,1,1); plot(this.t, this.thcv, this.t(1:end-1), diff(p)*this.sfrq);
            %         subplot(2,1,2); plot(this.t, this.thca, this.t(1:end-2), diff(diff(p))*this.sfrq^2);
            %         legend('thca','diff'); title('Check compuation of Angular Components');
            
            
        end
        
        function [q_0, q_0_dot, x_0, x_0_dot] = get_ZFT(this,tt)
            
            
            % Line motion
            %             D = 1;
            %             A = 0.1;
            %             tstart = 0;
            %             [x_0,x_0_dot] = getMinJerkTraj(this,D,A,tstart,tt);
            
            % x_0 fixed
            xe = this.d1 + 0.1;
            ye = this.d2;
            x_0 = [xe;ye];
            
            x_0_dot = [zeros(2,1)];
            
            % Right Hand
            alpha1 = atan2(ye,xe);
            alpha2 = acos((xe.^2 + ye.^2 + this.l1a^2 - this.l2a^2)./(2*this.l1a*sqrt(xe.^2 + ye.^2)));
            alpha3 = acos((xe.^2 + ye.^2 + this.l2a^2 - this.l1a^2)./(2*this.l2a*sqrt(xe.^2 + ye.^2)));
            
            % joint space zft postion
            q1p = alpha1 - alpha2;
            q2p = alpha1 + alpha3 - q1p;
            q_0 = [q1p;q2p];
            
            q_0_dot = this.get_jacobian(q_0) \ x_0_dot; % Check outputs incorrect q_0_dot when multiple values are input
            
            %             q_0 = interp1(this.t, this.q_0, tt)';
            %             q_0_dot = interp1(this.t, this.q_0_dot, tt)';
            %             x_0 = interp1(this.t, this.X_0, tt)';
            %             x_0_dot = this.get_jacobian(q_0)*q_0_dot;
            
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
        
        
        function xx_dot = getStates_point(this,tt,xx)
            % Filename:	getStatesZFT.m
            % Author:  James Hermus
            % Date:     16 July 2018
            % Description:	Computes the zero force trajectory from subject data
            % collected on the InMotion2with
            % Note: The computation uses realitive damping
            
            % Extract states
            q=xx(1:2);
            q_dot=xx(3:4);
            
            X = forwKino(this,q);
            
            % Define Jacobian of the arm
            J = this.get_jacobian(q);
            X_dot = J*q_dot; % Compute hand velocity
            
            [thcp,ee,nn] = get_crankPostion(this,X,'shoulder');
            
            M = this.get_massMatrixArm(q);
            
            % Define centifugal and Coriolis forces of the arm
            h = [ -this.m2a*this.l1a*this.c2a*sin(q(2)) * (2*q_dot(1)*q_dot(2) + q_dot(2)^2) ;...
                this.m2a*this.l1a*this.c2a*sin(q(2)) * q_dot(1)^2];
            
            
            % Force from InMotion Constraint
            F_c =  get_ForceFromConstraint_Point(this,X,X_dot,'shoulder');
            
            % Force Preturbation from InMotion
            F_pret = [interp1(this.t, this.F_pret(:,1), tt, 'nearest');...
                interp1(this.t, this.F_pret(:,2), tt, 'nearest')];
            
            % Human Joint toruqe
            q_0 = [0;0];
            [q_0(1), q_0(2)] = this.invKino(this.d1+this.lc,this.d2); % q0 point
            
            tau = this.K*(q_0 - q) + this.B*([0;0]-q_dot); % Point sim
            
            %% Future
            % [theta_r, theta_dr, J_r, J_dr, h_r, M_r, PlotPosVec_X, PlotPosVec_Y] = InMotionPassiveForce(param,x, y, x_dot, y_dot);
            
            % No robot integration
            % Changed preturbation to normal tangential direction
            q_ddot = inv(M) * ( tau + J'*(F_pret+F_c) - h);
            
            
            % With robot constraint integration
            % x34dot = inv( (M + J'*inv(J_r')*M_r*inv(J_r)*J_r)) *(tau - h - J'*inv(J_r')*(M_r*inv(J_r')*(Jd*[q1v;q2v] - J_dr*theta_dr) + J_r'*F + h_r));
            
            % Return the state derivatives to ODE solver
            xx_dot = [q_dot; q_ddot];
        end
        
        function [xx_dot,F_contact] = getStates(this,tt,xx)
            % Filename:	getStatesZFT.m
            % Author:  James Hermus
            % Date:     16 July 2018
            % Description:	Computes the zero force trajectory from subject data
            % collected on the InMotion2with
            % Note: The computation uses realitive damping
            
            % Extract states
            q=xx(1:2);
            q_dot=xx(3:4);
            
            X = forwKino(this,q);
            
            % Define Jacobian of the arm
            J = this.get_jacobian(q);
            X_dot = J*q_dot; % Compute hand velocity
            
            M = this.get_massMatrixArm(q);
            
            % Define centifugal and Coriolis forces of the arm
            h = [ -this.m2a*this.l1a*this.c2a*sin(q(2)) * (2*q_dot(1)*q_dot(2) + q_dot(2)^2) ;...
                this.m2a*this.l1a*this.c2a*sin(q(2)) * q_dot(1)^2];
            
            % Force from InMotion Constraint
            [F_c] = get_ForceFromConstraint_Point(this,X,X_dot,tt,'shoulder');
            
            if(~isempty(this.pretScale))
                
                pretScale_1 = interp1([0;this.pretScale{1}(:,1);2*pi], [this.pretScale{1}(1,2);this.pretScale{1}(:,2);this.pretScale{1}(end,2)], thcp, 'nearest');
                pretScale_2 = interp1([0;this.pretScale{1}(:,1);2*pi], [this.pretScale{1}(1,3);this.pretScale{1}(:,3);this.pretScale{1}(end,3)], thcp, 'nearest');
                
                % Force Preturbation from InMotion
                F_pret = [(10^-4)*(1./pretScale_1)*interp1(this.t, this.F_pret(:,1), tt, 'nearest');...
                    (10^-4)*(1./pretScale_2)*interp1(this.t, this.F_pret(:,2), tt, 'nearest')];
            else
                % Force Preturbation from InMotion
                F_pret = [interp1(this.t, this.F_pret(:,1), tt, 'nearest');...
                    interp1(this.t, this.F_pret(:,2), tt, 'nearest')];
            end
            
            %% Human Joint toruqe
            [q_0,q_0_dot] = get_ZFT(this,tt);
            
            tau = this.K*(q_0 - q) + this.B*(q_0_dot-q_dot);
            
            %% Robot torque
            %           [J_r, J_r_dot, h_r, M_r] = this.InMotionPassiveForce(X, X_dot,'shoulder');
            
            % No robot integration
            % Changed preturbation to normal tangential direction
            q_ddot = inv(M) * ( tau + J'*(F_pret + F_c) - h); % CHECK THIS
            
            % Compute contact force
            F_contact = F_pret + F_c;
            
            % Return the state derivatives to ODE solver
            xx_dot = [q_dot; q_ddot];
            
        end
        
        function [F,F_ne] = get_ForceFromConstraint(this,X,X_dot,coordinatess)
            
            % Check coordinates are in crank space
            X = this.check_crankCoordinates(X,coordinatess);
            
            % Circle
            [theta,ee,nn] = get_crankPostion(this,X,'crank');
            
            delta_x = this.lc - sqrt(X(1)^2+X(2)^2);
            
            rob_K = 2500; % Stiffness
            rob_B = 40; % Uniform damping
            % F = K * [-x_con(dex) + x; -y_con(dex) + y];% + B*[x_dot;y_dot];
            
            %             % Force preturbation
            F = rob_K * delta_x * nn + rob_B*dot(-X_dot,nn)*nn; % Add Damping!!
            F_ne = [rob_K * delta_x + rob_B*dot(-X_dot,nn);0]; % only normal
            
            % Try adding damping in both dof
            %             F = rob_K * delta_x * nn + rob_B*(-X_dot); % Add Damping!!
            %             F_ne = [rob_K * delta_x + rob_B*dot(-X_dot,nn);rob_B*dot(-X_dot,ee)]; % only normal
            
        end
        
        function [F] = get_ForceFromConstraint_Point(this,X,X_dot,tt,coordinatess)
            
            % Check coordinates are in crank space
            X = this.check_crankCoordinates(X,coordinatess);
            
            rob_K = [2500,0;0,2500]; % Stiffness % Added coupling for test
            rob_B = [40,0;0,40]; % Uniform damping
            
            if(tt<=2)
                F = rob_K * ([0;0]-X) + rob_B*([0;0]-X_dot); % Add Damping!!
            else % After release
                F = [0;0];
            end
        end
        
        function [fp] = get_instanceOrRandomPreturbation(this)
            
            fp = 1*rand(1,this.N_tot);
            fp = fp-mean(fp);
            
            %             % Filter
            %             cf = 100; % cutoff freqnency
            %             [b,a] = butter(2,cf/(this.sfrq/2)); % make filter
            %             fp = filtfilt(b,a,fp); % apply fitler
            %
            %             for i = 1:this.N_tot
            %                 if(~mod(i-1,10))
            %                     fp(i) = 1*(rand(1)-0.5);
            %                 else
            %                     fp(i) = fp(i-1);
            %                 end
            %             end
            
            % Make bionary amplitude
            fp(find(fp >= 0)) = 5; % Use 0.1 before
            fp(find(fp < 0)) = -5;
            
        end
        
        function [fp] = get_StepPreturbation(this)
            
            fp = zeros(1,this.N_tot);
            
            Lag = 15;
            % Get random time preturbation
            dex = randi([Lag this.N_tot-Lag],1);
            
            % Get random sign  pretubation
            A = randn(1);
            
            % Make bionary amplitude
            A(find(A >= 0)) = 5; % Use 0.1 before
            A(find(A < 0)) = -5;
            
            if(dex<200)
                A = 4*A;
            end
            
            fp(dex:dex+10) = A;
            
        end
        
        function [fp] = get_linStepPreturbation(this)
            
            fp = zeros(1,this.N_tot);
            
            Lag = 15;
            % Get random time preturbation
            dex = randi([Lag this.N_tot-Lag],1);
            
            % Get random sign  pretubation
            A = randn(1);
            
            % Make bionary amplitude
            A(find(A >= 0)) = 5; % Use 0.1 before
            A(find(A < 0)) = -5;
            
            if(dex<200)
                A = 4*A;
            end
            
            fp(dex:dex+10) = A;
            fp(dex) = A/2;
            fp(dex+10) = A/2;
            
        end
        
        function [fp] = get_continuousStepPreturbation(this)
            
            fp = zeros(1,this.N_tot);
            
            Lag = 15;
            % Get random time preturbation
            dex = randi([Lag this.N_tot-Lag],1);
            
            % Get random sign  pretubation
            A = randn(1);
            
            % Make bionary amplitude
            A(find(A >= 0)) = 5; % Use 0.1 before
            A(find(A < 0)) = -5;
            
            if(dex<200)
                A = 4*A;
            end
            
            w = 24;
            fp(dex-w/2:dex+w/2) = 1;
            
            dexUp = dex-w/4:dex-w/8;
            t_up = linspace(-pi/2, pi/2,length(dexUp));
            fp(dexUp) = (0.5)*(sin(t_up)+1);
            
            dexDown = dex+w/8:dex+w/4;
            t_up = linspace(pi/2, -pi/2,length(dexDown));
            fp(dexDown) = (0.5)*(sin(t_up)+1);
            
            fp = A.*fp;
            
            % Sanity check interp problem
%             clear all
%             close all
%             clc
%             
%             sfrq = 500;
%             dt = 1/sfrq;
%             tmax = 0.5;
%             t = 0:dt:tmax-dt;
%             N = length(t);
%             d = 125;
%             w = 74;
%             for r = 1:5
%                 % Make x1 (box car function)
%                 x1 = zeros(N,1);
%                 x1(d-w/2:d+w/2) = 1;
%                 
%                 % interp x1
%                 sfrq_high = 2000;
%                 dt_high = 1/sfrq_high;
%                 t_high = 0:dt_high:tmax-dt_high;
%                 x1_interp = interp1(t,x1,t_high,'spline');
%                 
%                 % Make x2 (sine edge)
%                 x2 = zeros(N,1);
%                 x2(d-w/2:d+w/2) = 1;
%                 
%                 dexUp = d-w/2:d-w/2+r;
%                 t_up = linspace(-pi/2, pi/2,length(dexUp));
%                 x2(dexUp) = (0.5)*(sin(t_up)+1);
%                 
%                 dexDown = d+w/2-r:d+w/2;
%                 t_up = linspace(pi/2, -pi/2,length(dexDown));
%                 x2(dexDown) = (0.5)*(sin(t_up)+1);
%                 
%                 % interp x2
%                 x2_interp = interp1(t,x2,t_high,'spline');
%                 
%                 % x3 (linear edge)
%                 x3 = x1;
%                 x3(dexUp) = linspace(0,1,length(dexUp));
%                 x3(dexDown) = linspace(1,0,length(dexDown));
%                 
%                 % Actual
%                 t_r = 0:1/5000:tmax;
%                 x_r = zeros(length(t_r),1);
%                 x_r(d*10-w*10/2:d*10+w*10/2) = 1;
%                 
%                 % interp x3
%                 x3_interp = interp1(t,x3,t_high,'spline');
%                 
%                 figure;plot(t,x1,'r-o',t_high,x1_interp,'r-+',...
%                     t,x2,'b-o',t_high,x2_interp,'b-+',...
%                     t,x3,'g-o',t_high,x3_interp,'g-+',...
%                     t_r,x_r,'k--','linewidth',2,'markersize',10);
%                 ylim([-0.25 1.25]); title(['r = ',int2str(r)]);
%             end
            
        end
        
        function [Fp] = get_StepPreturbation_2D(this)
            
            Fp = zeros(2,this.N_tot);
            
            Lag = 15;
            pretLength = 10;
            % Get random time preturbation
            dex = randi([Lag this.N_tot-Lag],1);
            
            % Get random sign  pretubation
            A = randn(1);
            theta = pi*rand();
            
            % Make bionary amplitude
            A(find(A >= 0)) = 5; % Use 0.1 before
            A(find(A < 0)) = -5;
            
            if(dex<200)
                A = 4*A;
            end
            
            Fp(:,dex:dex+pretLength-1) = [A*cos(theta); A*sin(theta)].*ones(2,pretLength);
            
        end
        
        function [q1p, q2p, q1v, q2v ] = invKino(this, xe, ye, xve, yve )
            %Compute joint angles
            %   Computes crank postion and joint angles input params provides the
            %   subject data from the get"kino number"param()
            
            % Requires postion data in shoulder coordinates
            
            if(this.handSide == 1) % If handSide == 1 use left hand
                
                % Left Hand
                alpha1 = atan2(ye,xe);
                alpha2 = acos((xe.^2 + ye.^2 + this.l1a^2 - this.l2a^2)./(2*this.l1a*sqrt(xe.^2 + ye.^2)));
                alpha3 = acos((xe.^2 + ye.^2 + this.l2a^2 - this.l1a^2)./(2*this.l2a*sqrt(xe.^2 + ye.^2)));
                
                % spatial angles
                q1p = alpha1 + alpha2;
                q2p = -(alpha2 + alpha3);
                
                % elbow angle absolute
                thea = q1p + q2p;
                
            elseif(this.handSide == 0) % If handSide 0 use right hand
                
                % Right Hand
                alpha1 = atan2(ye,xe);
                alpha2 = acos((xe.^2 + ye.^2 + this.l1a^2 - this.l2a^2)./(2*this.l1a*sqrt(xe.^2 + ye.^2)));
                alpha3 = acos((xe.^2 + ye.^2 + this.l2a^2 - this.l1a^2)./(2*this.l2a*sqrt(xe.^2 + ye.^2)));
                
                % spatial angles
                q1p = alpha1 - alpha2;
                q2p = alpha1 + alpha3 - q1p;
                
                % elbow angle absolute
                thea = q1p + q2p;
                
            else
                error('No hand side was specified.');
            end
            
            if(nargin > 3 )
                
                for i = 1:length(xe)
                    % Define Jacobian of the arm
                    J = this.get_jacobian([q1p(i),q2p(i)]);
                    qv = inv(J)*[xve(i);yve(i)];
                    q1v(i,1) = qv(1);
                    q2v(i,1) = qv(2);
                    
                end
                
            end
            
        end
        
        function [X] = forwKino(this,q)
            
            [q] = get_VecforArrayFunction(this, q);
            
            q1p = q(:,1);
            q2p = q(:,2);
            % Requires postion data in shoulder coordinates
            % Exports forward kinomatics in Shoulder centered coordinates
            
            x = [this.l1a*cos(q1p) + this.l2a*cos(q1p + q2p)];
            y = [this.l1a*sin(q1p) + this.l2a*sin(q1p + q2p)];
            
            % Recombine for vector and individual
            [X] = recombineVecforArrayFunction(this, x, y);
            
        end
        
        function [p,ee,nn] = get_crankPostion(this,X,coordinatess)
            
            % Check input is for a single sample
            if(size(X,1)==2 && size(X,2)==1)
                % Proceed
            else
                error('Shape of X is incorrect');
            end
            
            % Check coordinates are in crank space
            X = this.check_crankCoordinates(X,coordinatess);
            
            x = X(1);
            y = X(2);
            p = atan2(y, x);
            indxs = find(p(:,1)<0);
            p(indxs,1) = p(indxs) + 2*pi;
            
            % nn % Normal unit vector
            % ee % tangential unit vector
            ee = [ - sin(p); cos(p) ];
            nn = [  cos(p); sin(p) ];
            
        end
        
        function [p] = get_crankPostion_Vec(this,X,coordinatess)
            
            if(size(X,1) >= 2 && size(X,2) >= 2 )
                % Proceed
            else
                error('Shape of X is incorrect');
            end
            
            % Check coordinates are in crank space
            X = this.check_crankCoordinates(X,coordinatess);
            
            x = X(:,1);
            y = X(:,2);
            p = atan2(y, x);
            indxs = find(p(:,1)<0);
            p(indxs,1) = p(indxs) + 2*pi;
            
        end
        
        function [J] = get_jacobian(this,q)
            
            q1p = q(1);
            q2p = q(2);
            
            J11 = -( this.l1a*sin(q1p) + this.l2a*sin(q1p + q2p) );
            J12 = -this.l2a*sin(q1p + q2p);
            J21 = this.l1a*cos(q1p) + this.l2a*cos(q1p + q2p);
            J22 = this.l2a*cos(q1p + q2p);
            J = [J11, J12; J21, J22];
            
        end
        
        function [M] = get_massMatrixArm(this, q)
            
            if(size(q,2)>2)
                error('Shape of X is incorrect');
            end
            
            q1p = q(1);
            q2p = q(2);
            
            % Define inertial matrix of the arm
            M11 = this.m1a*this.c1a^2 + this.m2a*( this.l1a^2 + this.c2a^2 + 2*this.l1a*this.c2a*cos(q2p) ) + this.I1a + this.I2a;
            M12 = this.m2a*( this.c2a^2 + this.l1a*this.c2a*cos(q2p) ) + this.I2a;
            M22 = this.m2a*this.c2a^2 + this.I2a;
            M = [ M11, M12; M12, M22 ];
            
        end
        
        function [q_1, q_2] = robotInvKino(this,x,y)
            % Filename:	RobotInvKino.m
            % Author:  James Hermus
            % Date:		September 13 2017
            % Description:  Takes as inputs the handel x and y postion and computes the
            % correponding joint angles.
            
            % Takes input in InMotion2 coordinates
            
            % IMPORTANT NOTE: These joint angles and masses are spesific to the
            % InMotion two and are diffrent from the global subject joint angles and
            % interial parameters.
            
            % q_1 % InMotion joint angle 1
            % q_2 % InMotion joint angle 2
            
            % Inverse Kinomatics
            % Compute crank postion in carteasion
            ld = sqrt(x.^2 + y.^2);
            
            % Compute elbow angle relative
            q_2_rev = pi+acos((-ld.^2+this.rob_l1^2+this.rob_l2^2)./(2*this.rob_l1*this.rob_l2));
            
            % shoulder angle
            q_1 = atan2(y,x) + acos((ld.^2+this.rob_l1^2-this.rob_l2^2)./(2.*ld.*this.rob_l1));
            
            % elbow angle absolute
            q_2 = q_1 + q_2_rev;
            
        end
        
        function [J_r, J_r_dot, h_r, M_r] = InMotionPassiveForce(this, X, X_dot, coordinatess)
            
            % CHECK THIS FUNCTION COPY AND PASTED
            % Filename:	InMotionPassiveForce.m
            % Author:  James Hermus
            % Date:		September 13 2017
            % Description:  Run parfor decomposition using the LGNB code modified by
            % Hamal on the data imported from hypothesis 1 subject structure.
            % Analysis will exclude the first 1.5 seconds of each trial.
            % Decomposition in sections of 5 seconds.  The fast cases were excluded.
            %Computes the intertial and coriolus torque which results from the
            %InMotion2
            
            % IMPORTANT NOTE: These joint angles and masses are spesific to the
            % InMotion two and are diffrent from the global subject joint angles and
            % interial parameters.
            
            % get individual postions
            xe = X(1);
            ye = X(2);
            
            % get individual velocities
            x_dot = X_dot(1);
            y_dot = X_dot(2);
            
            % Takes in crank coorinates
            if(strcmp(coordinatess,'crank'))
                POS = this.convert_cartRefFrame([xe;ye],0,2);
            elseif(strcmp(coordinatess,'shoulder'))
                POS = this.convert_cartRefFrame([xe;ye],1,2);
            else
                disp('Spesify coordinates');
            end
            
            xe = POS(1);
            ye = POS(2);
            
            [q_1, q_2] = this.robotInvKino(xe, ye);
            
            % Get mass, center of mass distance, and inertia
            [l1, l2, l3, l4, m1, m2, m3, m4, r1, r2, r3, r4, I1, I2, I3, I4] = this.getRobotParam(); % for with in function convenience
            
            % Jacobian end-effector (J_ef)
            J_r = [-l1*sin(q_1), -l2*sin(q_2); l1*cos(q_1), l2*cos(q_2)];
            
            q = J_r\[x_dot;y_dot];
            q_d1 = q(1);
            q_d2 = q(2);
            
            % Jacobian dot
            J_r_dot = [-l1*q_d1*cos(q_1), -l2*q_d2*cos(q_2); -l1*q_d1*sin(q_1), -l2*q_d2*sin(q_2)];
            
            M_r = [ I1 + m2*l1^2 + I3,  (m2*l1*r2 + m3*l4*r3)*cos(q_1 - q_2); ...
                (m2*l1*r2 + m3*l4*r3)*cos(q_1 - q_2),  I2 + m3*l4^2 + I4];
            
            h_r = [ 0, (m2*l2*r2 + m3*l4*r3)*sin(q_1 - q_2)*q_d2;...
                -(m2*l1*r2 + m3*l4*r3)*sin(q_1-q_2)*q_d1, 0 ]*[q_d1;q_d2];
            
            theta_r = [q_1;q_2];
            theta_dr = [q_d1; q_d2];
        end
        
        function [l1, l2, l3, l4, m1, m2, m3, m4, r1, r2, r3, r4, I1, I2, I3, I4] = getRobotParam(this)
            
            % InMotion2 Parameters
            % Length
            l1 = this.rob_l1;
            l2 = this.rob_l2;
            l3 = this.rob_l3;
            l4 = this.rob_l4;
            % Mass for the corresponding links
            m1 = this.rob_m1;
            m2 = this.rob_m2;
            m3 = this.rob_m3;
            m4 = this.rob_m4;
            % Distance between center of mass and joint closest to the robot link
            r1 = this.rob_r1;
            r2 = this.rob_r2;
            r3 = this.rob_r3;
            r4 = this.rob_r4;
            % Inertia at the end of the link closest to the motor.
            I1 = this.rob_I1;
            I2 = this.rob_I2;
            I3 = this.rob_I3;
            I4 = this.rob_I4;
            
        end
        
        function [armPosVec_X, armPosVec_Y] = plotPosArm(this,xe,ye,coordinatess)
            
            % Takes in crank coorinates
            if(strcmp(coordinatess,'crank'))
                POS = this.convert_cartRefFrame([xe,ye],0,1);
            elseif(strcmp(coordinatess,'shoulder'))
                POS = [xe,ye];
            else
                disp('Spesify coordinates');
            end
            
            [q1p, q2p] = this.invKino( POS(:,1),POS(:,2) );
            
            % Returns the arm postion for plotting
            armPosVec_X = [zeros(length(q1p),1), this.l1a*cos(q1p),...
                this.l1a*cos(q1p) + this.l2a*cos(q1p + q2p)];
            
            armPosVec_Y = [zeros(length(q1p),1), this.l1a*sin(q1p),...
                this.l1a*sin(q1p) + this.l2a*sin(q1p + q2p)];
            
            % figure; plot(armPosVec_X, armPosVec_Y,'-ok'); axis equal; grid on;
            
            
        end
        
        function [robPosVec_X, robPosVec_Y] = plotPosInMotion(this,xe,ye,coordinatess)
            
            % Takes in crank coorinates
            if(strcmp(coordinatess,'crank'))
                POS = this.convert_cartRefFrame([xe,ye],0,2);
            elseif(strcmp(coordinatess,'shoulder'))
                POS = this.convert_cartRefFrame([xe,ye],1,2);
            else
                disp('Spesify coordinates');
            end
            
            [q_1,q_2] = this.robotInvKino(POS(:,1),POS(:,2));
            
            % Returns the InMotion postion for plotting in shoulder coordinates
            robPosVec_X = [  this.rob_l4*cos(q_2) + this.rob_l3*cos(q_1), this.rob_l4*cos(q_2),...
                zeros(size(q_1)), ...
                this.rob_l1*cos(q_1),...
                this.rob_l1*cos(q_1) + this.rob_l2*cos(q_2) ]+this.rob_d1 + this.d1;
            robPosVec_Y = [   this.rob_l4*sin(q_2) + this.rob_l3*sin(q_1),...
                this.rob_l4*sin(q_2),...
                zeros(size(q_1)),...
                this.rob_l1*sin(q_1),...
                this.rob_l1*sin(q_1) + this.rob_l2*sin(q_2)]+this.rob_d2 + this.d2;
        end
        
        function [] = get_zftPlot(this)
            % Zero Force Trajectory Plot
            
            figure('position',[700 250 560 420]); hold on;
            plot(this.d1 - this.lc*cos(0:0.01:2*pi),this.d2 - this.lc*sin(0:0.01:2*pi),'--','linewidth',2.5);
            plot(this.X_0(:,1),this.X_0(:,2),'Linewidth',2.5);
            ylabel('Y-Position (m)');
            xlabel('X-Position (m)');
            legend('Crank Path','Zero Force Trajectory','Location','SouthEast');
            ylim([0.2 0.7]);
            hold off; axis equal; grid on;
            set(gca,'FontSize',18);
            
        end
        
        function [q] = get_VecforArrayFunction(this, q)
            
            if(size(q)==[2,1])
                q = q';
            elseif(size(q,1)>2 && size(q,2) == 2)
                % do nothing
            else
                error('Shape of X is incorrect');
            end
            
        end
        
        function [X] = recombineVecforArrayFunction(this, x, y)
            
            if(length(x) == 1)
                X = [x;y];
            elseif(size(x,1) >= 2 && size(y,1) >= 2)
                X = [x,y];
            else
                error('Shape of x or y is incorrect');
            end
            
        end
        
        function [X] = check_crankCoordinates(this,X,coordinatess)
            
            if(strcmp(coordinatess,'crank'))
                % Proceed
            elseif(strcmp(coordinatess, 'shoulder'))
                [X] = convert_cartRefFrame(this,X,1,0);
            else
                error('Spesific crank or shoulder coordinates');
            end
            
        end
        
        function [] = save_zftAnimation(this,videoFileName)
            
            figure;
            v = VideoWriter(videoFileName, 'MPEG-4');
            skipSamples = 5;
            v.FrameRate = this.sfrq/skipSamples;
            open(v);
            
            for i = 1:skipSamples:length(this.t)
                
                
                [thcp,ee,nn] = get_crankPostion(this,this.X(i,:)','shoulder');
                
                ax = plot(this.armPosVec_X(i,:), this.armPosVec_Y(i,:), '-b',...
                    this.armPosVec_X(i,:), this.armPosVec_Y(i,:), '.b',...
                    'markersize',35,'linewidth',3);
                %                 this.crank_x(i,:),     this.crank_y(i,:),     '-r',...
                %                 this.crank_x(i,:),     this.crank_y(i,:),     '.k',...
                %                     this.constraint_x,     this.constraint_y,     ':k',...
                %                     this.X_0(i,1),           this.X_0(i,2),       '.g',...
                %                     this.robPosVec_X(i,:), this.robPosVec_Y(i,:), '-k',...
                %                     this.X(i,1)+0.1*[0,ee(1)],this.X(i,2)+0.1*[0,ee(2)],':b',...
                %                     this.X(i,1)+0.1*[0,nn(1)],this.X(i,2)+0.1*[0,nn(2)],':r',...
                %                     this.robPosVec_X(i,:), this.robPosVec_Y(i,:), '.k',...
                ylim([0 1.2]); xlim([-0.56 0.56]); grid on; axis equal;
                ylabel('Y Distance (m)','fontsize',14);
                xlabel('X Distance (m)','fontsize',14);
                %                 legend([ax(1),ax(2),ax(3),ax(4),ax(5),ax(8)],{'Arm','Constraint Path','Robot','Tangental','Normal','Zero-Force Trajectory'},'location','northwest'); % Add ZFT
                ylabel('Y Distance (m)','fontsize',14);
                xlabel('X Distance (m)','fontsize',14);
                xlim([-0.6 0.6]);
                
                frame = getframe(gcf);
                writeVideo(v,frame);
            end
            close(v);
            
        end
        
        function [y_r_1, y_r_2, u_r_1, u_r_2, f_r_1, f_r_2] = get_ensambleInput(this,symType)
            % Define sim parameters
            N = length(this.t); % Check N is even
            
            % Crop start and end
            % ??? Figure this out soon
            
            %             f_r_1 = zeros(1,N);
            %             f_r_2 = zeros(1,N);
            
            if(strcmp(symType,'line') || strcmp(symType,'line_SISO'))
                % Dimension 1
                y_r_1 = this.X(:,1)';
                u_r_1 = this.F_pret(:,1)';
                f_r_1 = this.F_contact(:,1)';
                
                % Dimension 2
                y_r_2 = this.X(:,2)';
                u_r_2 = this.F_pret(:,2)';
                f_r_2 = this.F_contact(:,2)';
                
            else
                error('symType must be: line.');
            end
            
        end
        
        function [t_list, x_list] = SemiImplicitEuler(this,func,tspan, y0, options)
            
            % Semi-implicit Euler Integrator
            % function @(t,y) e.g. f=@(t,y)(t+y);
            % t0 = initial condition of t
            % y0 = initial condition of y = [q0; Dq0]
            % step size
            
            dt = this.dt;
            t = (tspan(1):dt:tspan(2));
            y = zeros(numel(y0),numel(t));
            y(:,1) = y0;
            
            nq = numel(y0)/2;
            q = zeros(nq, numel(t));
            Dq = zeros(nq, numel(t));
            q(:,1) = y0(1:nq);
            Dq(:,1) = y0(nq+1:end);
            
            % semi-implicit Euler integration
            for i = 1:(numel(t)-1)
                DDy = func(t(i),y(:,i));
                Dq(:,i+1) = Dq(:,i) + dt*DDy(nq+1:end);
                q(:,i+1) = q(:,i) + dt*Dq(:,i+1);
                y(:,i+1) = [q(:,i+1); Dq(:,i+1)];
            end
            
            
            t_list = t';
            x_list = y';
            
        end
        
        function [POS] = convert_cartRefFrame(this,POS,oldFrame,newFrame)
            
            [POS] = get_VecforArrayFunction(this, POS);
            
            % 0: crank coordinates
            % 1: shoulder coordinates
            % 2: inmotion coordinates
            
            % If frame is already correct
            if(oldFrame == newFrame)
                % do nothing
            elseif(oldFrame == 0 && newFrame == 1)
                POS = POS + [this.d1,this.d2];
            elseif(oldFrame == 0 && newFrame == 2)
                POS = POS - [this.rob_d1,this.rob_d2];
            elseif(oldFrame == 1 && newFrame == 0)
                POS = POS - [this.d1,this.d2];
            elseif(oldFrame == 1 && newFrame == 2)
                POS = POS - [this.d1,this.d2] - [this.rob_d1,this.rob_d2];
            elseif(oldFrame == 2 && newFrame == 0)
                POS = POS + [this.rob_d1,this.rob_d2];
            elseif(oldFrame == 2 && newFrame == 1)
                POS = POS + [this.rob_d1,this.rob_d2] + [this.d1,this.d2];
            else
                error('Enter correct tranformation for refrence frame.');
            end
            
            % % Sanity Check Back Up
            % X_shoulder = this.convert_cartRefFrame(this.X,0,2);
            % X_rob = this.convert_cartRefFrame(X_shoulder,2,1);
            %
            % figure;
            % subplot(1,3,1); plot(this.X(:,1),this.X(:,2),0,0,'o'); axis equal; grid on; xlim([-0.3 0.3]); ylim([-1 1]);
            % subplot(1,3,2); plot(X_shoulder(:,1),X_shoulder(:,2),0,0,'o');axis equal; grid on; xlim([-0.3 0.3]); ylim([-1 1]);
            % subplot(1,3,3); plot(X_rob(:,1),X_rob(:,2),0,0,'o');axis equal; grid on; xlim([-0.3 0.3]); ylim([-1 1]);
            
            % If singular change shape back
            if( size(POS,1)== 1 && size(POS,2) == 2)
                POS = POS';
            end
            
        end
        
    end
end

