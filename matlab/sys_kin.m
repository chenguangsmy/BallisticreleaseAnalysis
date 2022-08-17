function results = sys_kin(Data,idx_t,time_t,subj,dir)
% Sys Ident - THIS SECTION OF CODE TAKES a while to run
clc, close all

tic

N = 100;
ang = 0:((2*pi)/N):(2*pi);

gain = 1;
K_joint = gain*[29.5 14.3;14.3 39.3];

xdd = cos(ang);
results.xdd = xdd;
ydd = sin(ang);
results.ydd = ydd;

for f_sel= 1:3
    for d_sel = 1:3

        f_sel
        d_sel
        switch f_sel
            case 1
                Ftg = 15;
            case 2
                Ftg = 20;
            case 3
                Ftg = 25;
            otherwise
                disp('No Force Value')
        end

        switch d_sel
            case 1
                Dtg = 25;
            case 2
                Dtg = 50;
            case 3
                Dtg = 75;
            otherwise
                disp('No Displacement Value')
        end

        % Select Perturbation (1=unperturbed,2-7=perturbed with different starting
        % times)
        unpert = 1;

        trial_l = length(Data(subj,dir,f_sel,d_sel,:,unpert));

        %Create Average Profile of Unperturbed
        time_new = linspace(0,1.25,1.25e3);
        
        %To Compute Average Kinematic Mass at End-Effector across Trials
        Msum_t0 = zeros(2,2);
        Msum_tf = zeros(2,2);
        idx_sum0 = 1;
        idx_sumf = 1;

        Ksum_t0 = zeros(2,2);
        Ksum_tf = zeros(2,2);
        idx_Ksum0 = 1;
        idx_Ksumf = 1;


        for i = 1:trial_l
            if size(Data{subj,dir,f_sel,d_sel,i,unpert}) == [1 1]
                time_up{i,:} = time_t{subj,dir,f_sel,d_sel,i,1};
                force_interp(i,:) = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.f(1,idx_t{subj,dir,f_sel,d_sel,i,1}),time_new); %.f(2,.. for old .f(1,.. for new
                disp_interp(i,:) = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(1,idx_t{subj,dir,f_sel,d_sel,i,1}),time_new); %.x(2,.. for old .x(1,.. for new

                results.hand{f_sel,d_sel,i}.x =  interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(1,idx_t{subj,dir,f_sel,d_sel,i,1},1),time_new);
                results.elbow{f_sel,d_sel,i}.x = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(1,idx_t{subj,dir,f_sel,d_sel,i,1},2),time_new);
                results.shoulder{f_sel,d_sel,i}.x = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(1,idx_t{subj,dir,f_sel,d_sel,i,1},3),time_new);

                results.hand{f_sel,d_sel,i}.y =  interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(2,idx_t{subj,dir,f_sel,d_sel,i,1},1),time_new);
                results.elbow{f_sel,d_sel,i}.y = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(2,idx_t{subj,dir,f_sel,d_sel,i,1},2),time_new);
                results.shoulder{f_sel,d_sel,i}.y = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(2,idx_t{subj,dir,f_sel,d_sel,i,1},3),time_new);

                results.hand{f_sel,d_sel,i}.z =  interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(3,idx_t{subj,dir,f_sel,d_sel,i,1},1),time_new);
                results.elbow{f_sel,d_sel,i}.z = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(3,idx_t{subj,dir,f_sel,d_sel,i,1},2),time_new);
                results.shoulder{f_sel,d_sel,i}.z = interp1(time_up{i,:},Data{subj,dir,f_sel,d_sel,i,1}.ox(3,idx_t{subj,dir,f_sel,d_sel,i,1},3),time_new);


                diff_he = (results.hand{f_sel,d_sel,i}.x-results.elbow{f_sel,d_sel,i}.x).^2 + ...
                    (results.hand{f_sel,d_sel,i}.y-results.elbow{f_sel,d_sel,i}.y).^2 + ...
                    (results.hand{f_sel,d_sel,i}.z-results.elbow{f_sel,d_sel,i}.z).^2;
                diff_es = (results.elbow{f_sel,d_sel,i}.x-results.shoulder{f_sel,d_sel,i}.x).^2 + ...
                    (results.elbow{f_sel,d_sel,i}.y-results.shoulder{f_sel,d_sel,i}.y).^2 + ...
                    (results.elbow{f_sel,d_sel,i}.z-results.shoulder{f_sel,d_sel,i}.z).^2;
                diff_hs = (results.hand{f_sel,d_sel,i}.x-results.shoulder{f_sel,d_sel,i}.x).^2 + ...
                    (results.hand{f_sel,d_sel,i}.y-results.shoulder{f_sel,d_sel,i}.y).^2 + ...
                    (results.hand{f_sel,d_sel,i}.z-results.shoulder{f_sel,d_sel,i}.z).^2;

                %Forearm Length
                results.forearm{f_sel,d_sel,i} = sqrt(diff_he);
                %Arm Length
                results.arm{f_sel,d_sel,i} = sqrt(diff_es);
                %Hand-Shoulder Distance
                results.lhs{f_sel,d_sel,i} = sqrt(diff_hs);

                for tt = length(results.arm{f_sel,d_sel,i}):-1:1
                    if results.arm{f_sel,d_sel,i}(tt) >= 0.4
                        results.hand{f_sel,d_sel,i}.x(tt) = [];
                        results.hand{f_sel,d_sel,i}.y(tt) = [];
                        results.hand{f_sel,d_sel,i}.z(tt) = [];

                        results.elbow{f_sel,d_sel,i}.x(tt) = [];
                        results.elbow{f_sel,d_sel,i}.y(tt) = [];
                        results.elbow{f_sel,d_sel,i}.z(tt) = [];

                        results.shoulder{f_sel,d_sel,i}.x(tt) = [];
                        results.shoulder{f_sel,d_sel,i}.y(tt) = [];
                        results.shoulder{f_sel,d_sel,i}.z(tt) = [];

                        results.forearm{f_sel,d_sel,i}(tt) = [];
                        results.arm{f_sel,d_sel,i}(tt) = [];
                        results.lhs{f_sel,d_sel,i}(tt) = [];
                    end
                end

                lhs_p = abs(results.hand{f_sel,d_sel,i}.x-results.shoulder{f_sel,d_sel,i}.x);
                les_p = abs(results.elbow{f_sel,d_sel,i}.x-results.shoulder{f_sel,d_sel,i}.x);

                %Shoulder Angle
                results.th_s{f_sel,d_sel,i} = acosd(les_p./results.arm{f_sel,d_sel,i});

                alph = acosd(lhs_p./results.lhs{f_sel,d_sel,i});

                th_he = 180 - alph - results.th_s{f_sel,d_sel,i};
                % Elbow Angle
                results.th_e{f_sel,d_sel,i} = asind((results.lhs{f_sel,d_sel,i}/results.forearm{f_sel,d_sel,i}).*sind(th_he));

                %Anthropometric Parameters
                %(based on Winter Biomechanics book 2009)
                M = 75; %[kg] user mass
                m1 = 0.028*M; %[kg] upper arm mass
                m2 = 0.022*M; %[kg] forearm+hand mass
                l1G{f_sel,d_sel,i} = 0.436*results.arm{f_sel,d_sel,i};
                l2G{f_sel,d_sel,i} = 0.682*results.forearm{f_sel,d_sel,i};

                rhoG1{f_sel,d_sel,i} = 0.322*results.arm{f_sel,d_sel,i};
                rhoG2{f_sel,d_sel,i} = 0.468*results.forearm{f_sel,d_sel,i};

                I1{f_sel,d_sel,i} = m1*rhoG1{f_sel,d_sel,i}.^2;
                I2{f_sel,d_sel,i} = m2*rhoG2{f_sel,d_sel,i}.^2;

                
                for tt = 1:length(results.forearm{f_sel,d_sel,i})
                    %Joint Space Mass
                    results.M{f_sel,d_sel,i,tt} = [m1.*l1G{f_sel,d_sel,i}(tt)^2+m2.*(results.arm{f_sel,d_sel,i}(tt)^2+2.*results.arm{f_sel,d_sel,i}(tt)*l2G{f_sel,d_sel,i}(tt)*cosd(180-results.th_e{f_sel,d_sel,i}(tt))+l2G{f_sel,d_sel,i}(tt)^2)+I1{f_sel,d_sel,i}(tt) , m2.*(results.arm{f_sel,d_sel,i}(tt)*l2G{f_sel,d_sel,i}(tt)*cosd(180-results.th_e{f_sel,d_sel,i}(tt))+l2G{f_sel,d_sel,i}(tt)^2); ...
                        m2.*(results.arm{f_sel,d_sel,i}(tt)*l2G{f_sel,d_sel,i}(tt)*cosd(180-results.th_e{f_sel,d_sel,i}(tt))+l2G{f_sel,d_sel,i}(tt)^2) , m2.*l2G{f_sel,d_sel,i}(tt)^2+I2{f_sel,d_sel,i}(tt)];

                    %Jacobian Matrix
                    results.J{f_sel,d_sel,i,tt} = [-results.arm{f_sel,d_sel,i}(tt)*sind(results.th_s{f_sel,d_sel,i}(tt))-results.forearm{f_sel,d_sel,i}(tt)*sind(results.th_s{f_sel,d_sel,i}(tt)+(180-results.th_e{f_sel,d_sel,i}(tt))), -results.forearm{f_sel,d_sel,i}(tt)*sind(results.th_s{f_sel,d_sel,i}(tt)+(180-results.th_e{f_sel,d_sel,i}(tt))); ...
                        +results.arm{f_sel,d_sel,i}(tt)*cosd(results.th_s{f_sel,d_sel,i}(tt))+results.forearm{f_sel,d_sel,i}(tt)*cosd(results.th_s{f_sel,d_sel,i}(tt)+(180-results.th_e{f_sel,d_sel,i}(tt))), +results.forearm{f_sel,d_sel,i}(tt)*cosd(results.th_s{f_sel,d_sel,i}(tt)+(180-results.th_e{f_sel,d_sel,i}(tt)))];

                    %Apparent End-Effector Mass Matrix
                    results.MEE{f_sel,d_sel,i,tt} = inv(results.J{f_sel,d_sel,i,tt}')*results.M{f_sel,d_sel,i,tt}*inv(results.J{f_sel,d_sel,i,tt});

                    %Apparent End-Effector Stiffness
                    results.KEE{f_sel,d_sel,i,tt} = inv(results.J{f_sel,d_sel,i,tt}')*K_joint*inv(results.J{f_sel,d_sel,i,tt});
                    
                    %Computation of Mass Ellipsoids
                    results.Fxx{f_sel,d_sel,i,tt} = results.MEE{f_sel,d_sel,i,tt}(1,:)*[xdd;ydd];
                    results.Fyy{f_sel,d_sel,i,tt} = results.MEE{f_sel,d_sel,i,tt}(2,:)*[xdd;ydd];

                    %Computation of Stiffness Ellipsoids
                    results.Fxx_k{f_sel,d_sel,i,tt} = results.KEE{f_sel,d_sel,i,tt}(1,:)*[xdd;ydd];
                    results.Fyy_k{f_sel,d_sel,i,tt} = results.KEE{f_sel,d_sel,i,tt}(2,:)*[xdd;ydd];
                end
                
                if size(results.MEE{f_sel,d_sel,i,1}) == [2 2]
                    Msum_t0 = Msum_t0+results.MEE{f_sel,d_sel,i,1};
                    idx_sum0 = idx_sum0+1;
                    Ksum_t0 = Ksum_t0+results.KEE{f_sel,d_sel,i,1};
                    idx_Ksum0 = idx_Ksum0+1;
                end
                if size(results.MEE{f_sel,d_sel,i,end}) == [2 2]
                    Msum_tf = Msum_tf+results.MEE{f_sel,d_sel,i,end};
                    idx_sumf = idx_sumf+1;
                    Ksum_tf = Ksum_tf+results.KEE{f_sel,d_sel,i,end};
                    idx_Ksumf = idx_Ksumf+1;
                end
            end
        end
        results.MEE_avg_t0{f_sel,d_sel} = Msum_t0/(idx_sum0-1);
        results.MEE_avg_tf{f_sel,d_sel} = Msum_tf/(idx_sumf-1);
        results.KEE_avg_t0{f_sel,d_sel} = Ksum_t0/(idx_Ksum0-1);
        results.KEE_avg_tf{f_sel,d_sel} = Ksum_tf/(idx_Ksumf-1);
    end
end

results.t = time_new;
toc
end