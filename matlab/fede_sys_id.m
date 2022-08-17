clc, clear, close all

% LOAD THE DATA

% load('ss3938_3949.mat', 'data');  % 6N perturbation on x, 
% load('ss4129_4137.mat', 'data');  % 6N and 12N perturbation on x, only chenguang,
% load('ss4129_4179.mat', 'data');  % 6N and 12N perturbation on x,
% himanshu and chenguang
% load('ss4181_4202.mat', 'data');  % 6N and 12N perturbation on x, himanshu and chenguang
% load('ss4216_4226.mat'); %Multiple Direction Chenguang/Himanshu
% load('ss4216_4239.mat'); %Multiple Direction Chenguang/Himanshu stricter conditions
% load('ss4253_4274.mat') %Single Direction Cheng/Himan with Pulse in force and position hold

%4Subject Multiple Direction
load('ss4310_4356.mat');

% Chenguang adapted code to briefly plot the loaded data and check that
% everythin is fine! ;)


Data = data;
Freq = 500;
t_step = 1/500;
clear axh
fh = figure(2); 

colors = colormap('lines');
r = size(Data, 1); % subj
c = size(Data, 2); % direction
f = size(Data, 3); % force
d = size(Data, 4); % target
l = size(Data, 5); % trials
p = size(Data, 6); % perturbation type
idx_last = 200;
if_subtract = 0;
epoc_type = 2;
plot_type = 1;          % 1 displacement
                        % 2 force
                        % 3 force command
                        % 4 velocity
                        % 5 torque of 3rd joint
                        % 6 displacement (vector length in cartesian space)
                        % 7 force (vector length in cartesian space)
axh = zeros(f,r);


% In this for loop it is defined a variable called idx_t that trims the
% data only over the interested piece (the ballistic release)
% There are two options for it now, one for pure unperturbed ballistic
% release and another for unpert+pulses.
% REMEMBER to switch between the two to correctly analyze data!!!!

for ri = 1:r % subj
    for ci = 1:c
        for fi = 1:f % fce
            for di = 1:d % target distance
                axh(ri, fi) = subplot(d,f,d*(fi-1) + di);grid on;hold on;
                %for di = 3 % target distance
                for li = 1:p % perturbation
                    trial_num = length(Data(ri,ci,fi,di,:,li));
                    for ti = 1:trial_num % each trial
                        if (isempty(Data{ri,ci,fi,di,ti,li}))
                            continue;
                        end

                        switch epoc_type
                            case 1
                                idx = find(Data{ri,ci,fi,di,ti,li}.Fp(2,:)~=0 & ...
                                    Data{ri,ci,fi,di,ti,li}.ts==4);  % pert at y
                                idx = [idx idx(end)+(1:idx_last)]; % may error as the pert not long enough
                                if li == 1
                                    disp('ERROR: should use li == 2!!!');
                                end
                            case 2
                                idx = find(Data{ri,ci,fi,di,ti,li}.ts==5 | Data{ri,ci,fi,di,ti,li}.ts==6);
                                %For Pulses t_in = release - 1.2s
%                                 idx = (idx(1)-600):idx(end);
                                
                                %For Multiple t_in = release - 0.5s
                                idx = (idx(1)-250):idx(end); 
                                
                                %idx = (idx(1)):idx(end);
                        end

                        time = t_step*(idx-idx(1));
                        idx_t{ri,ci,fi,di,ti,li} = idx;
                        time_t{ri,ci,fi,di,ti,li} = t_step*(idx-idx(1));
                        %time = idx-idx(1);
                        switch plot_type
                            case 1
                                dat = Data{ri,ci,fi,di,ti,li}.ox(1,idx);
                                %dat = dat - dat(1);
                                titlestr = 'displacement';
                            case 2
                                dat = Data{ri,ci,fi,di,ti,li}.f(1,idx);
                                titlestr = 'force';
                            case 3
                                dat = Data{ri,ci,fi,di,ti,li}.Fp(1,idx);
                                titlestr = 'Fp';
                            case 4
                                dat = Data{ri,ci,fi,di,ti,li}.v(1,idx);
                                titlestr = 'velocity';
                            case 5
                                dat = Data{ri,ci,fi,di,ti,li}.tq(3,idx);
                                titlestr = 'torque3';
                            case 6
                                dat = Data{ri,ci,fi,di,ti,li}.x(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm displacement';
                            case 7 % the force mode
                                dat = Data{ri,ci,fi,di,ti,li}.f(:,idx);
                                dat_submean = dat - mean(dat(:,1:50),2);
                                dat_norm = sqrt(dat_submean(1,:).^2 + dat_submean(2,:).^2 + dat_submean(3,:).^2);
                                dat = dat_norm .* sign(dat_submean(2,:));
                                titlestr = 'norm force';

                        end
                        if (if_subtract)
                            dat = dat - mean(dat(1:50));
                        end
                        plot(time, dat, 'Color', colors(4*(li-1)+di, :));
                        %                     plot(time, dat, 'Color', [0.7 0.7 0.7]);
                    end
                end
            end
        end
    end
end
%xlim([0 0.7])
linkaxes(axh, 'xy');
%xlim([0 0.5]);
% xlim([0 2])
sgtitle(titlestr);

%% System Identification

clc, close all
%Select Subject
subj = 4;
% Select Direction (1=+x,2=+y,3=-x,4=-y)
dir = 4;

%Here you have the functions to perfrom system identification depending 
%on the kind of test you performed: Unpeturbed+Pulses or UnperturbedOnly

%------Function Computing system Identification Unpert+Pulses (remember to
%change idx_t above)
% results = sys_id_pp(Data,idx_t,time_t,subj);

%------Function Computing system Identification Unpert Only (remember to
%change idx_t above)
results = sys_id_ubr(Data,idx_t,time_t,subj,dir);

%The following function is to compute the kinematic analysis based on the Optotrak
%markers on hand, elbow and shoulder
%(As far as I have seen now, it works well only for direction "+x" e
%subject Chenguamng)

%------Function Computing Arm Kinematic Analysis (Joint Angles, Jacobian, Mass Matrix)
% results = sys_kin(Data,idx_t,time_t,subj,dir);

%% PLOTTING (the rest of the code concerns only results plotting)

% IN the following subsections you can find the different plots for the
% analyzed data, run the specific subsection depending on what you are
% interested in:
% - Time plots;
% - Frequency plots;
% - Bar plots (for Unperturbed only)
% - Kinematic plots
% After the plotting you have also a little piece of code to perform ANOVA
% on the unperturbed ballistic release data

%% TIME PLOTS

clc, close all

%Uncomment for pulses data
% wind_start = 20; %From which time window to start showing

% Axes Limits (POSITIVE directions)
% f_maxx = 30;
% f_minn = -0.5;
% d_minn = -0.01;
% d_maxx = +0.1;

% % Axes Limits (NEGATIVE directions)
f_maxx = +0.5;
f_minn = -30;
d_minn = -0.1;
d_maxx = +0.01;


%Unperturbed Ballistic Release Force [N]
figure(),
set(gcf,'color','w');
sgtitle('Force [N] in Unperturbed Ballistic Release')
idx = 1;

for f_sel = 1:3
    for d_sel = 1:3
        subplot(3,3,idx),
        plot(results.FD_UP{f_sel,d_sel}{1,1},results.FD_UP{f_sel, d_sel}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
        plot(results.avg_FD_UP{f_sel,d_sel}(1,:),results.avg_FD_UP{f_sel,d_sel}(2,:),'r','LineWidth',2), grid on
        ylim([f_minn f_maxx])
        idx = idx+1;
    end
end

%Unperturbed Ballistic Release Displacement [m]
figure(),
set(gcf,'color','w');
sgtitle('Displacement [m] in Unperturbed Ballistic Release')
idx = 1;

for f_sel = 1:3
    for d_sel = 1:3
        subplot(3,3,idx),
         plot(results.FD_UP{f_sel, d_sel}{1, 1},results.FD_UP{f_sel, d_sel}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
         plot(results.avg_FD_UP{f_sel, d_sel}(1,:),results.avg_FD_UP{f_sel, d_sel}(3,:),'r','LineWidth',2), hold on
         plot(results.tt_mod{f_sel, d_sel},results.disp_mod{f_sel, d_sel},'b--','LineWidth',2),
         %Unperturbed 9Trials
         lgd = legend({'','','','','','','','','','',strcat('FIT: ',num2str(round(results.FIT_up_avg(f_sel,d_sel))),'%')},'Location','best');
         %Pulses 20 Trials
         %lgd = legend({'','','','','','','','','','','','','','','','','','','','','',strcat('FIT: ',num2str(round(results.FIT_up_avg(1,1))),'%')},'Location','best');
         lgd.FontSize = 7;
         lgd.Box = 'off';
         grid on
         ylim([d_minn d_maxx])
        idx = idx+1;
    end
end



%% FREQUENCY PLOTS
%Bode Plots
clc, close all

options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
% options.MagVisible = 'off';
options.YLim = {[-80 -40]};

figure(),
set(gcf,'color','w');
sgtitle('Bode Plot - Magnitude - Identified Models')
idx = 1;

for f_sel = 1:3
    for d_sel = 1:3
        subplot(3,3,idx),
        bodemag(results.TF.up{f_sel, d_sel},{0.1*2*pi,10*2*pi},'b',options), hold on,
        grid on
        idx = idx+1;
    end
end

%% 2-D Bar Plots Unperturbed 
clc, close all

ff = [15 20 25];
xx = [25 50 75];
% ff = {'F1' 'F2' 'F3'};
% xx = categorical({'900' '600' '300'});



%Stiffness 

figure(),
subplot(1,3,1)
set(gcf,'color','w');
b = bar(xx,results.K_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.K_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.K_up_avg',results.K_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')



%Damping

subplot(1,3,2)
set(gcf,'color','w');
b = bar(xx,results.B_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.B_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.B_up_avg',results.B_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')


%Inertia 

subplot(1,3,3)
set(gcf,'color','w');
b = bar(xx,results.M_up_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.M_up_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.M_up_avg',results.M_up_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')

% Natural Frequency and Damping Ratio

figure(),
set(gcf,'color','w');

% Natural Frequency

subplot(1,2,1)
set(gcf,'color','w');
b = bar(xx,results.up_omegan_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.up_omegan_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.up_omegan_avg',results.up_omegan_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
% ylim([0 20])
ylabel('Natural Frequency [rad/s]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')

% Damping Ratio

subplot(1,2,2),
set(gcf,'color','w');
b = bar(xx,results.up_dampr_avg'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(results.up_dampr_avg');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',results.up_dampr_avg',results.up_dampr_std','k','linestyle','none');
grid on
xlabel('Displacement [mm]')
% xlabel('Target Stiffness [N/m]')
ylim([0 1])
ylabel('Damping Ratio [-]')
legend('15 N','20 N','25 N')
% legend('S','M','L')
title('Ballistic Release')


%% STATISTICAL ANALYSIS
clc, close all

%ANOVA-test on unperturbed
KKK=[];
BBB=[];
MMM=[];
ooo=[];
ddd=[];
ForceK=[];
ForceB=[];
ForceM=[];
Forceo=[];
Forced=[];
dispK=[];
dispB=[];
dispM=[];
dispo=[];
dispd=[];

for f_sel = 1:3
    for d_sel = 1:3
            KKK = [KKK;results.K_up_tr{f_sel,d_sel,1}'];
            BBB = [BBB;results.B_up_tr{f_sel,d_sel,1}'];
            MMM = [MMM;results.M_up_tr{f_sel,d_sel,1}'];
            ooo = [ooo;results.up_omegan_tr{f_sel,d_sel,1}'];
            ddd = [ddd;results.up_dampr_tr{f_sel,d_sel,1}'];

            if f_sel == 1
                ForceK = [ForceK;15*ones(length(results.K_up_tr{f_sel,d_sel,1}),1)];
                ForceB = [ForceB;15*ones(length(results.B_up_tr{f_sel,d_sel,1}),1)];
                ForceM = [ForceM;15*ones(length(results.M_up_tr{f_sel,d_sel,1}),1)];
                Forceo = [Forceo;15*ones(length(results.up_omegan_tr{f_sel,d_sel,1}),1)];
                Forced = [Forced;15*ones(length(results.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif f_sel == 2
                ForceK = [ForceK;20*ones(length(results.K_up_tr{f_sel,d_sel,1}),1)];
                ForceB = [ForceB;20*ones(length(results.B_up_tr{f_sel,d_sel,1}),1)];
                ForceM = [ForceM;20*ones(length(results.M_up_tr{f_sel,d_sel,1}),1)];
                Forceo = [Forceo;20*ones(length(results.up_omegan_tr{f_sel,d_sel,1}),1)];
                Forced = [Forced;20*ones(length(results.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif f_sel == 3
                ForceK = [ForceK;25*ones(length(results.K_up_tr{f_sel,d_sel,1}),1)];
                ForceB = [ForceB;25*ones(length(results.B_up_tr{f_sel,d_sel,1}),1)];
                ForceM = [ForceM;25*ones(length(results.M_up_tr{f_sel,d_sel,1}),1)];
                Forceo = [Forceo;25*ones(length(results.up_omegan_tr{f_sel,d_sel,1}),1)];
                Forced = [Forced;25*ones(length(results.up_dampr_tr{f_sel,d_sel,1}),1)];
            end

            if d_sel == 1
                dispK = [dispK;25*ones(length(results.K_up_tr{f_sel,d_sel,1}),1)];
                dispB = [dispB;25*ones(length(results.B_up_tr{f_sel,d_sel,1}),1)];
                dispM = [dispM;25*ones(length(results.M_up_tr{f_sel,d_sel,1}),1)];
                dispo = [dispo;25*ones(length(results.up_omegan_tr{f_sel,d_sel,1}),1)];
                dispd = [dispd;25*ones(length(results.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif d_sel == 2
                dispK = [dispK;50*ones(length(results.K_up_tr{f_sel,d_sel,1}),1)];
                dispB = [dispB;50*ones(length(results.B_up_tr{f_sel,d_sel,1}),1)];
                dispM = [dispM;50*ones(length(results.M_up_tr{f_sel,d_sel,1}),1)];
                dispo = [dispo;50*ones(length(results.up_omegan_tr{f_sel,d_sel,1}),1)];
                dispd = [dispd;50*ones(length(results.up_dampr_tr{f_sel,d_sel,1}),1)];
            elseif d_sel == 3
                dispK = [dispK;75*ones(length(results.K_up_tr{f_sel,d_sel,1}),1)];
                dispB = [dispB;75*ones(length(results.B_up_tr{f_sel,d_sel,1}),1)];
                dispM = [dispM;75*ones(length(results.M_up_tr{f_sel,d_sel,1}),1)];
                dispo = [dispo;75*ones(length(results.up_omegan_tr{f_sel,d_sel,1}),1)];
                dispd = [dispd;75*ones(length(results.up_dampr_tr{f_sel,d_sel,1}),1)];
            end

    end
end
      
p_stiff = anovan(KKK,{ForceK dispK},'model','interaction','varnames',{'Force','Displacement'})
p_damp = anovan(BBB,{ForceB dispB},'model','interaction','varnames',{'Force','Displacement'})
p_mass = anovan(MMM,{ForceM dispM},'model','interaction','varnames',{'Force','Displacement'})
p_omega = anovan(ooo,{Forceo dispo},'model','interaction','varnames',{'Force','Displacement'})
p_damp = anovan(ddd,{Forced dispd},'model','interaction','varnames',{'Force','Displacement'})

%%
clearvars KKK MMM BBB ooo ddd ForceK Forced Forceo ForceM ForceB dispM dispB dispK dispo dispd 
clearvars p_stiff p_damp p_mass p_omega p_damp
clearvars results

%% KINEMATIC PLOTS - run this part to plot the results of the kinematic analysis
clc, close all
trial_l = 9;

%Marker Trajectories
figure(),
set(gcf,'color','w');
sgtitle('Trajectories of the Hand-Elbow-Shoulder Joints')
idx = 1;

for f_sel = 1:3
    for d_sel = 1:3
        for i = 1:trial_l
            subplot(3,3,idx),
            plot3(results.hand{f_sel,d_sel,i}.x,results.hand{f_sel,d_sel,i}.y,results.hand{f_sel,d_sel,i}.z,'.k'), hold on
            plot3(results.elbow{f_sel,d_sel,i}.x,results.elbow{f_sel,d_sel,i}.y,results.elbow{f_sel,d_sel,i}.z,'.r'), hold on
            plot3(results.shoulder{f_sel,d_sel,i}.x,results.shoulder{f_sel,d_sel,i}.y,results.shoulder{f_sel,d_sel,i}.z,'.b'), hold on, grid on
            xlabel('X [m]')
            ylabel('Y [m]')
            zlabel('Z [m]')
            xlim([-0.5 0])
            ylim([-0.9 -0.5])
            view(2)
            axis equal
        end
        idx = idx+1;
    end
end

%Mass Ellipsoid
en_f = 8; %Enalrgement factor of Arm Lengths
figure(),
set(gcf,'color','w');
sgtitle('End-Effector Mass Ellipsoid')
idx = 1;
for f_sel = 1:3
    for d_sel = 1:3
        for i = 1:trial_l
            
            subplot(3,3,idx),
            plot(results.Fxx{f_sel,d_sel,i,1},results.Fyy{f_sel,d_sel,i,1},'.k'), hold on
            plot(results.Fxx{f_sel,d_sel,i,end},results.Fyy{f_sel,d_sel,i,end},'.r'), hold on
            plot(en_f*[0,results.forearm{f_sel,d_sel,i}(1).*cosd(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))],en_f*[0,-results.forearm{f_sel,d_sel,i}(1).*sind(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[results.forearm{f_sel,d_sel,i}(1).*cosd(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1)),results.forearm{f_sel,d_sel,i}(1).*cosd(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))-results.arm{f_sel,d_sel,i}(1).*cosd(results.th_s{f_sel,d_sel,i}(1))], ...
                en_f*[-results.forearm{f_sel,d_sel,i}(1).*sind(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1)),-results.forearm{f_sel,d_sel,i}(1).*sind(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))-results.arm{f_sel,d_sel,i}(1).*sind(results.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[0,results.forearm{f_sel,d_sel,i}(end).*cosd(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))],en_f*[0,-results.forearm{f_sel,d_sel,i}(end).*sind(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            plot(en_f*[results.forearm{f_sel,d_sel,i}(end).*cosd(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end)),results.forearm{f_sel,d_sel,i}(end).*cosd(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))-results.arm{f_sel,d_sel,i}(end).*cosd(results.th_s{f_sel,d_sel,i}(end))], ...
                en_f*[-results.forearm{f_sel,d_sel,i}(end).*sind(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end)),-results.forearm{f_sel,d_sel,i}(end).*sind(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))-results.arm{f_sel,d_sel,i}(end).*sind(results.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            grid on
            xlabel('Fx [N]')
            ylabel('Fy [N]')
            [t,s] = title(strcat('M_x(t_0)= ',num2str(results.MEE_avg_t0{f_sel,d_sel}(1,1),'%.1f'),' kg , M_x(t_f)= ',num2str(results.MEE_avg_tf{f_sel,d_sel}(1,1),'%.1f'),' kg'),...
                strcat('M_y(t_0)= ',num2str(results.MEE_avg_t0{f_sel,d_sel}(2,2),'%.1f'),' kg , M_y(t_f)= ',num2str(results.MEE_avg_tf{f_sel,d_sel}(2,2),'%.1f'),' kg'));
            t.FontSize = 9;
            t.FontWeight = 'normal';
            s.FontSize = 9;
            s.FontWeight = 'normal';
            axis equal
        end
        idx = idx+1;
    end
end

%Stiffness Ellipsoid
en_f = 800; %Enalrgement factor of Arm Lengths
figure(),
set(gcf,'color','w');
sgtitle('End-Effector Stiffness Ellipsoid')
idx = 1;
for f_sel = 1:3
    for d_sel = 1:3
        for i = 1:trial_l
            
            subplot(3,3,idx),
            plot(results.Fxx_k{f_sel,d_sel,i,1},results.Fyy_k{f_sel,d_sel,i,1},'.k'), hold on
            plot(results.Fxx_k{f_sel,d_sel,i,end},results.Fyy_k{f_sel,d_sel,i,end},'.r'), hold on
            plot(en_f*[0,results.forearm{f_sel,d_sel,i}(1).*cosd(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))],en_f*[0,-results.forearm{f_sel,d_sel,i}(1).*sind(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[results.forearm{f_sel,d_sel,i}(1).*cosd(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1)),results.forearm{f_sel,d_sel,i}(1).*cosd(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))-results.arm{f_sel,d_sel,i}(1).*cosd(results.th_s{f_sel,d_sel,i}(1))], ...
                en_f*[-results.forearm{f_sel,d_sel,i}(1).*sind(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1)),-results.forearm{f_sel,d_sel,i}(1).*sind(results.th_e{f_sel,d_sel,i}(1)-results.th_s{f_sel,d_sel,i}(1))-results.arm{f_sel,d_sel,i}(1).*sind(results.th_s{f_sel,d_sel,i}(1))],'Color','k','LineWidth',2), hold on
            plot(en_f*[0,results.forearm{f_sel,d_sel,i}(end).*cosd(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))],en_f*[0,-results.forearm{f_sel,d_sel,i}(end).*sind(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            plot(en_f*[results.forearm{f_sel,d_sel,i}(end).*cosd(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end)),results.forearm{f_sel,d_sel,i}(end).*cosd(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))-results.arm{f_sel,d_sel,i}(end).*cosd(results.th_s{f_sel,d_sel,i}(end))], ...
                en_f*[-results.forearm{f_sel,d_sel,i}(end).*sind(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end)),-results.forearm{f_sel,d_sel,i}(end).*sind(results.th_e{f_sel,d_sel,i}(end)-results.th_s{f_sel,d_sel,i}(end))-results.arm{f_sel,d_sel,i}(end).*sind(results.th_s{f_sel,d_sel,i}(end))],'Color','r','LineWidth',2), hold on
            grid on
            xlabel('Fx [N]')
            ylabel('Fy [N]')
            [t,s] = title(strcat('K_x(t_0)= ',num2str(results.KEE_avg_t0{f_sel,d_sel}(1,1),'%.1f'),' N/m , K_x(t_f)= ',num2str(results.KEE_avg_tf{f_sel,d_sel}(1,1),'%.1f'),' N/m'),...
                strcat('K_y(t_0)= ',num2str(results.KEE_avg_t0{f_sel,d_sel}(2,2),'%.1f'),' N/m , K_y(t_f)= ',num2str(results.KEE_avg_tf{f_sel,d_sel}(2,2),'%.1f'),' N/m'));
            t.FontSize = 9;
            t.FontWeight = 'normal';
            s.FontSize = 9;
            s.FontWeight = 'normal';
            axis equal
        end
        idx = idx+1;
    end
end

%Angular Plots
figure(),
set(gcf,'color','w');
sgtitle('Elbow-Shoulder Angles')
idx = 1;
for f_sel = 1:3
    for d_sel = 1:3
        for i = 1:trial_l
            
            subplot(3,3,idx),
            plot(results.th_s{f_sel,d_sel,i},results.th_e{f_sel,d_sel,i}), hold on
            grid on
            xlabel('\theta_s [deg]')
            ylabel('\theta_e [deg]')
%             [t,s] = title(strcat('M_x(t_0)= ',num2str(results.MEE_avg_t0{f_sel,d_sel}(1,1),'%.1f'),' kg , M_x(t_f)= ',num2str(results.MEE_avg_tf{f_sel,d_sel}(1,1),'%.1f'),' kg'),...
%                 strcat('M_y(t_0)= ',num2str(results.MEE_avg_t0{f_sel,d_sel}(2,2),'%.1f'),' kg , M_y(t_f)= ',num2str(results.MEE_avg_tf{f_sel,d_sel}(2,2),'%.1f'),' kg'));
            t.FontSize = 9;
            t.FontWeight = 'normal';
            s.FontSize = 9;
            s.FontWeight = 'normal';
            axis equal
        end
        idx = idx+1;
    end
end

