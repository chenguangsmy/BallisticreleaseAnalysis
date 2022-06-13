% Change ln20, to choose plot things: force, displacement, velocity, etc.
clc, clear, close all

% load('ss3938_3949.mat', 'data');  % 6N perturbation on x, 
% load('ss4129_4137.mat', 'data');  % 6N and 12N perturbation on x, only chenguang,
% load('ss4129_4179.mat', 'data');  % 6N and 12N perturbation on x,
% himanshu and chenguang
% load('ss4181_4202.mat', 'data');  % 6N and 12N perturbation on x, himanshu and chenguang
% load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4216_4220.mat', 'data');  % just release, Chenguang
load('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/processedData/ss4253_4274.mat', 'data'); % release and perturb, Chenguang and Himanshu
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
%                                 idx = (idx(1)-250):idx(end);
                                idx = (idx(1)-600):idx(end);
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
xlim([0 2])
sgtitle(titlestr);

%% Sys Ident - THIS SECTION OF CODE TAKES ABOUT 1000secs to RUN!!!
clc, close all
%Select Subject
subj = 1;

%Function Computing system Identification
results = sys_id_pitt(Data,idx_t,time_t,subj);

%% TO RUN ONLY IF YOU WANT TO CANCEL BAD FITS

for f_sel= 1:3
    for d_sel = 1:3

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

        for pert = 2:4
            if results.FIT_p{pert-1}.p(f_sel,d_sel) < 85
                results.K_p{pert-1}.p(f_sel,d_sel) = NaN;
                results.B_p{pert-1}.p(f_sel,d_sel) = NaN;
                results.M_p{pert-1}.p(f_sel,d_sel)= NaN;
                results.FIT_p{pert-1}.p(f_sel,d_sel) = NaN;
            end
        end
    end
end

%% STATISTICAL ANALYSIS
clc, close all

%T-test
for i = 1:9
    if i == 1
        f_sel = 1;
        d_sel = 1;
    elseif i == 2
        f_sel = 1;
        d_sel = 2;
    elseif i == 3
        f_sel = 1;
        d_sel = 3;
    elseif i == 4
        f_sel = 2;
        d_sel = 1;
    elseif i == 5
        f_sel = 2;
        d_sel = 2;
    elseif i == 6
        f_sel = 2;
        d_sel = 3;
    elseif i == 7
        f_sel = 3;
        d_sel = 1;
    elseif i == 8
        f_sel = 3;
        d_sel = 2;
    elseif i == 9
        f_sel = 3;
        d_sel = 3;
    end

    for j = 1:9

        if j == 1
            f_selt = 1;
            d_selt = 1;
        elseif j == 2
            f_selt = 1;
            d_selt = 2;
        elseif j == 3
            f_selt = 1;
            d_selt = 3;
        elseif j == 4
            f_selt = 2;
            d_selt = 1;
        elseif j == 5
            f_selt = 2;
            d_selt = 2;
        elseif j == 6
            f_selt = 2;
            d_selt = 3;
        elseif j == 7
            f_selt = 3;
            d_selt = 1;
        elseif j == 8
            f_selt = 3;
            d_selt = 2;
        elseif j == 9
            f_selt = 3;
            d_selt = 3;
        end

        Kxx = results.K_up_tr(f_sel,d_sel,:);
        Kyy = results.K_up_tr(f_selt,d_selt,:);
        Kxx = Kxx(:);
        Kyy = Kyy(:);
        [hK,pK] = ttest2(Kxx,Kyy);

        Bxx = results.B_up_tr(f_sel,d_sel,:);
        Byy = results.B_up_tr(f_selt,d_selt,:);
        Bxx = Bxx(:);
        Byy = Byy(:);
        [hB,pB] = ttest2(Bxx,Byy);

        Mxx = results.M_up_tr(f_sel,d_sel,:);
        Myy = results.M_up_tr(f_selt,d_selt,:);
        Mxx = Mxx(:);
        Myy = Myy(:);
        [hM,pM] = ttest2(Mxx,Myy);

        p_val_K(i,j) = pK;
        p_val_B(i,j) = pB;
        p_val_M(i,j) = pM;

        if p_val_K(i,j) > 0.05
            ver_K(i,j) = 1;
        else
            ver_K(i,j) = 0;
        end

        if p_val_B(i,j) > 0.05
            ver_B(i,j) = 1;
        else
            ver_B(i,j)= 0;
        end

        if p_val_M(i,j) > 0.05
            ver_M(i,j) = 1;
        else
            ver_M(i,j) = 0;
        end
    end

end

%ANOVA-test on unperturbed
idx = 1;
for i = 1:3
    for j = 1:3
        for k = 1:length(results.K_up_tr(i,j,:))
            
            KKK(idx,1) = results.K_up_tr(i,j,k);
            BBB(idx,1) = results.B_up_tr(i,j,k);
            MMM(idx,1) = results.M_up_tr(i,j,k);
            
        
            if i == 1
                Force(idx,1) = 15;
            elseif i == 2
                Force(idx,1) = 20;
            elseif i == 3
                Force(idx,1) = 25;
            end

            if j == 1
                Disp(idx,1) = 25;
            elseif j == 2
                Disp(idx,1) = 50;
            elseif j == 3
                Disp(idx,1) = 75;
            end


            idx = idx+1;
        end
    end
end
      
p_stiff = anovan(KKK,{Force Disp},'model','interaction','varnames',{'Force','Displacement'})
p_damp = anovan(BBB,{Force Disp},'model','interaction','varnames',{'Force','Displacement'})
p_mass = anovan(MMM,{Force Disp},'model','interaction','varnames',{'Force','Displacement'})

%ANOVA-test on pulses
% wind_start = 20;
% idx = 1;
% for i = 1:3
%     for j = 1:3
%         for k = wind_start:length(results.wind_v)
%             
%             KKK_fh(idx,1) = results.K_p{1,1}.p(i,j,k);
%             BBB_fh(idx,1) = results.B_p{1,1}.p(i,j,k);
%             MMM_fh(idx,1) = results.M_p{1,1}.p(i,j,k);
% 
%             KKK_m(idx,1) = results.K_p{1,2}.p(i,j,k);
%             BBB_m(idx,1) = results.B_p{1,2}.p(i,j,k);
%             MMM_m(idx,1) = results.M_p{1,2}.p(i,j,k);
% 
%             KKK_ph(idx,1) = results.K_p{1,3}.p(i,j,k);
%             BBB_ph(idx,1) = results.B_p{1,3}.p(i,j,k);
%             MMM_ph(idx,1) = results.M_p{1,3}.p(i,j,k);
% 
% 
%             if i == 1
%                 Forcep(idx,1) = 15;
%             elseif i == 2
%                 Forcep(idx,1) = 20;
%             elseif i == 3
%                 Forcep(idx,1) = 25;
%             end
% 
%             if j == 1
%                 Dispp(idx,1) = 25;
%             elseif j == 2
%                 Dispp(idx,1) = 50;
%             elseif j == 3
%                 Dispp(idx,1) = 75;
%             end
% 
%             idx = idx+1;
%         end
%     end
% end
% 
% 
% p_stiff_fh = anovan(KKK_fh,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% p_damp_fh = anovan(BBB_fh,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% p_mass_fh = anovan(MMM_fh,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% 
% p_stiff_m = anovan(KKK_m,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% p_damp_m = anovan(BBB_m,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% p_mass_m = anovan(MMM_m,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% 
% p_stiff_ph = anovan(KKK_ph,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% p_damp_m = anovan(BBB_ph,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% p_mass_m = anovan(MMM_ph,{Forcep Dispp},'model','interaction','varnames',{'Force','Displacement'})
% 
%% FIGURES
clc
% close all,

tpause = 0.05; %Speed of plotting
wind_start = 20; %From which time window to start showing

ff = [15 20 25];
xx = [25 50 75];

[XX,FF] = meshgrid(xx,ff);

% Interpolate Impedance Matrix to get Surface

fq = 15:1:25;
xq = 25:1:75;

[Xq,Fq] = meshgrid(xq,fq);

Kup_interp = interp2(XX,FF,results.K_up,Xq,Fq);
Bup_interp = interp2(XX,FF,results.B_up,Xq,Fq);
Mup_interp = interp2(XX,FF,results.M_up,Xq,Fq);
FITup_interp = interp2(XX,FF,results.FIT_up,Xq,Fq);


%--------------------------------------------
% CODE TO INSER FOR VIDEOS
% % Initialize video
% myVideo = VideoWriter('myVideoFile'); %open video file
% myVideo.FrameRate = 15;  %can adjust this, 5 - 10 works well for me
% myvideo.Quality = 100;
% open(myVideo)
%     frame = getframe(gcf); %get frame
%     writeVideo(myVideo, frame);
% close(myVideo)


figure(),%'Units','normalized','Position',[0 0 1 1])
set(gcf,'color','w');
plot3(FF,XX,results.K_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.K_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.K_up_avg+results.K_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.K_up_avg-results.K_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
s = surf(Fq,Xq,Kup_interp); hold on
grid on
s.EdgeColor = 'none';
colorbar
colormap summer
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Stiffness [N/m]')
zlim([0 1000])
for widx = wind_start:length(results.wind_v)
% for widx = 10:length(results.wind_v)
    for p = 1%:3
        pplot = plot3(FF,XX,results.K_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
%         elseif p == 2
%             pplot(1).Marker = "o";
%             pplot(2).Marker = "o";
%             pplot(3).Marker = "o";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
%         elseif p == 3
%             pplot(1).Marker = "+";
%             pplot(2).Marker = "+";
%             pplot(3).Marker = "+";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 1000])
    pause(tpause)

end

figure(),
set(gcf,'color','w');
plot3(FF,XX,results.B_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.B_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.B_up_avg+results.B_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.B_up_avg-results.B_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
s = surf(Fq,Xq,Bup_interp); grid on
s.EdgeColor = 'none';
colorbar
colormap cool
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Damping [Ns/m]')
zlim([0 50]) 
for widx = wind_start:length(results.wind_v)
    for p = 1%:3
        pplot = plot3(FF,XX,results.B_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
%         elseif p == 2
%             pplot(1).Marker = "o";
%             pplot(2).Marker = "o";
%             pplot(3).Marker = "o";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
%         elseif p == 3
%             pplot(1).Marker = "+";
%             pplot(2).Marker = "+";
%             pplot(3).Marker = "+";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 50])
    pause(tpause)
end


figure(),
set(gcf,'color','w');
plot3(FF,XX,results.M_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.M_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.M_up_avg+results.M_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.M_up_avg-results.M_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
sM = surf(Fq,Xq,ones(size(Xq))*results.M_avg_tt); hold on
sM.EdgeColor = 'none';
sM.FaceColor = 'b';
sM.FaceAlpha = 0.3;
s = surf(Fq,Xq,Mup_interp); grid on
s.EdgeColor = 'none';
colorbar
colormap hot
caxis([1 5])
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Mass [kg]')
zlim([0 5])
for widx = wind_start:length(results.wind_v)
    for p = 1%:3
        pplot = plot3(FF,XX,results.M_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
%         elseif p == 2
%             pplot(1).Marker = "o";
%             pplot(2).Marker = "o";
%             pplot(3).Marker = "o";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
%         elseif p == 3
%             pplot(1).Marker = "+";
%             pplot(2).Marker = "+";
%             pplot(3).Marker = "+";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 5])
    pause(tpause)
end

figure()
set(gcf,'color','w');
plot3(FF,XX,results.FIT_up,'.k','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.FIT_up_avg,'.r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.FIT_up_avg+results.FIT_up_std,'+r','MarkerSize',10,'LineWidth',2), hold on
plot3(FF,XX,results.FIT_up_avg-results.FIT_up_std,'_r','MarkerSize',10,'LineWidth',2), hold on
s = surf(Fq,Xq,FITup_interp); grid on
s.EdgeColor = 'none';
colorbar
colormap spring
caxis([90 100])
xlabel('Force [N]')
ylabel('Displacement [mm]')
zlabel('Model FIT [%]')
zlim([0 100])
for widx = wind_start:length(results.wind_v)
    for p = 1%:3
        pplot = plot3(FF,XX,results.FIT_p{p}.p(:,:,widx),'k','MarkerSize',10); hold on
        if p == 1
            pplot(1).Marker = "x";
            pplot(2).Marker = "x";
            pplot(3).Marker = "x";
            pplot(1).LineStyle = 'none';
            pplot(2).LineStyle = 'none';
            pplot(3).LineStyle = 'none';
%         elseif p == 2
%             pplot(1).Marker = "o";
%             pplot(2).Marker = "o";
%             pplot(3).Marker = "o";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
%         elseif p == 3
%             pplot(1).Marker = "+";
%             pplot(2).Marker = "+";
%             pplot(3).Marker = "+";
%             pplot(1).LineStyle = 'none';
%             pplot(2).LineStyle = 'none';
%             pplot(3).LineStyle = 'none';
        end
    end
    grid on
    zlim([0 100])
    pause(tpause)
end

%% 2-D Bar Plots Unperturbed vs Pulses over Trials at Constant Window
clc, close all

tpause = 0.05; %Speed of plotting
wind_start = 20; %From which time window to start showing

% ff = [15 20 25];
% xx = [25 50 75];
ff = {'F1' 'F2' 'F3'};
xx = categorical({'900' '600' '300'});

widx =23;

%Stiffness 

figure(),
subplot(1,4,1)
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
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Ballistic Release')

subplot(1,4,2)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.K_p_avg{1}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.K_p_avg{1}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.K_p_avg{1}.p(:,:,widx))',cell2mat(results.K_p_std{1}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Force Hold')

subplot(1,4,3)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.K_p_avg{2}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.K_p_avg{2}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.K_p_avg{2}.p(:,:,widx))',cell2mat(results.K_p_std{2}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Motion')

subplot(1,4,4)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.K_p_avg{3}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.K_p_avg{3}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.K_p_avg{3}.p(:,:,widx))',cell2mat(results.K_p_std{3}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Position Hold')



%Damping

figure(),
subplot(1,4,1)
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
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Ballistic Release')

subplot(1,4,2)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.B_p_avg{1}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.B_p_avg{1}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.B_p_avg{1}.p(:,:,widx))',cell2mat(results.B_p_std{1}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Force Hold')

subplot(1,4,3)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.B_p_avg{2}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.B_p_avg{2}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.B_p_avg{2}.p(:,:,widx))',cell2mat(results.B_p_std{2}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Motion')

subplot(1,4,4)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.B_p_avg{3}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.B_p_avg{3}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.B_p_avg{3}.p(:,:,widx))',cell2mat(results.B_p_std{3}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Position Hold')

%Inertia 

figure(),
subplot(1,4,1)
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
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Ballistic Release')

subplot(1,4,2)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.M_p_avg{1}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.M_p_avg{1}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.M_p_avg{1}.p(:,:,widx))',cell2mat(results.M_p_std{1}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Force Hold')

subplot(1,4,3)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.M_p_avg{2}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.M_p_avg{2}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.M_p_avg{2}.p(:,:,widx))',cell2mat(results.M_p_std{2}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Motion')

subplot(1,4,4)
set(gcf,'color','w');
b = bar(xx,cell2mat(results.M_p_avg{3}.p(:,:,widx))'); hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(cell2mat(results.M_p_avg{3}.p(:,:,widx))');
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',cell2mat(results.M_p_avg{3}.p(:,:,widx))',cell2mat(results.M_p_std{3}.p(:,:,widx))','k','linestyle','none');
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Position Hold')

%% 2-D Bar Plots Unperturbed vs Average Pulses over Windows
clc, close all

tpause = 0.05; %Speed of plotting
wind_start = 20; %From which time window to start showing

% ff = [15 20 25];
% xx = [25 50 75];
ff = {'F1' 'F2' 'F3'};
xx = categorical({'900' '600' '300'});

widx =23;

%Stiffness 

figure(),
subplot(1,4,1)
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
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Ballistic Release')

subplot(1,4,2)
set(gcf,'color','w');
bar(xx,(results.K_p{1}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Force Hold')

subplot(1,4,3)
set(gcf,'color','w');
bar(xx,(results.K_p{2}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Motion')

subplot(1,4,4)
set(gcf,'color','w');
bar(xx,(results.K_p{3}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 800])
ylabel('Stiffness [N/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Position Hold')



%Damping

figure(),
subplot(1,4,1)
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
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Ballistic Release')

subplot(1,4,2)
set(gcf,'color','w');
bar(xx,(results.B_p{1}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Force Hold')

subplot(1,4,3)
set(gcf,'color','w');
bar(xx,(results.B_p{2}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Motion')

subplot(1,4,4)
set(gcf,'color','w');
bar(xx,(results.B_p{3}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 40])
ylabel('Damping [Ns/m]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Position Hold')

%Inertia 

figure(),
subplot(1,4,1)
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
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Ballistic Release')

subplot(1,4,2)
set(gcf,'color','w');
bar(xx,(results.M_p{1}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Force Hold')

subplot(1,4,3)
set(gcf,'color','w');
bar(xx,(results.M_p{2}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Motion')

subplot(1,4,4)
set(gcf,'color','w');
bar(xx,(results.M_p{3}.p(:,:,widx))'); hold on
grid on
% xlabel('Displacement [mm]')
xlabel('Target Stiffness [N/m]')
ylim([0 5])
ylabel('Mass [kg]')
% legend('15 N','20 N','25 N')
legend('S','M','L')
title('Pulse @ Position Hold')
%% TIME PLOTS
clc, close all
wind_start = 20; %From which time window to start showing
% Axes Limits
f_maxx = 55;
f_minn = -0.5;
d_minn = -0.01;
d_maxx = 0.16;

% Unperturbed
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_UP{1, 1}{1, 1},results.FD_UP{1, 1}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,1}(1,:),results.avg_FD_UP{1,1}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,2), plot(results.FD_UP{1, 2}{1, 1},results.FD_UP{1, 2}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,2}(1,:),results.avg_FD_UP{1,2}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,3), plot(results.FD_UP{1, 3}{1, 1},results.FD_UP{1, 3}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,3}(1,:),results.avg_FD_UP{1,3}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,4), plot(results.FD_UP{2, 1}{1, 1},results.FD_UP{2, 1}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,1}(1,:),results.avg_FD_UP{2,1}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,5), plot(results.FD_UP{2, 2}{1, 1},results.FD_UP{2, 2}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,2}(1,:),results.avg_FD_UP{2,2}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,6), plot(results.FD_UP{2, 3}{1, 1},results.FD_UP{2, 3}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,3}(1,:),results.avg_FD_UP{2,3}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,7), plot(results.FD_UP{3, 1}{1, 1},results.FD_UP{3, 1}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,1}(1,:),results.avg_FD_UP{3,1}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,8), plot(results.FD_UP{3, 2}{1, 1},results.FD_UP{3, 2}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,2}(1,:),results.avg_FD_UP{3,2}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,9), plot(results.FD_UP{3, 3}{1, 1},results.FD_UP{3, 3}{2, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,3}(1,:),results.avg_FD_UP{3,3}(2,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
sgtitle('Force [N] in Unperturbed Ballistic Release')


figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_UP{1, 1}{1, 1},results.FD_UP{1, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,1}(1,:),results.avg_FD_UP{1,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 1},results.disp_mod{1, 1},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,2), plot(results.FD_UP{1, 2}{1, 1},results.FD_UP{1, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,2}(1,:),results.avg_FD_UP{1,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 2},results.disp_mod{1, 2},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,3), plot(results.FD_UP{1, 3}{1, 1},results.FD_UP{1, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,3}(1,:),results.avg_FD_UP{1,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 3},results.disp_mod{1, 3},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,4), plot(results.FD_UP{2, 1}{1, 1},results.FD_UP{2, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,1}(1,:),results.avg_FD_UP{2,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 1},results.disp_mod{2, 1},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,5), plot(results.FD_UP{2, 2}{1, 1},results.FD_UP{2, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,2}(1,:),results.avg_FD_UP{2,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 2},results.disp_mod{2, 2},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,6), plot(results.FD_UP{2, 3}{1, 1},results.FD_UP{2, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,3}(1,:),results.avg_FD_UP{2,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 3},results.disp_mod{2, 3},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,7), plot(results.FD_UP{3, 1}{1, 1},results.FD_UP{3, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,1}(1,:),results.avg_FD_UP{3,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 1},results.disp_mod{3, 1},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,8), plot(results.FD_UP{3, 2}{1, 1},results.FD_UP{3, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,2}(1,:),results.avg_FD_UP{3,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 2},results.disp_mod{3, 2},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,9), plot(results.FD_UP{3, 3}{1, 1},results.FD_UP{3, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,3}(1,:),results.avg_FD_UP{3,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 3},results.disp_mod{3, 3},'b--','LineWidth',2), 
                grid on
                ylim([d_minn d_maxx])
sgtitle('Displacement [m] in Unperturbed Ballistic Release')

%Ballistic Response of Pulse Identified Models
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_UP{1, 1}{1, 1},results.FD_UP{1, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,1}(1,:),results.avg_FD_UP{1,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 1},results.disp_mod{1, 1},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{1,1,w},results.disp_modP{1}.p{1,1,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{1,1,w},results.disp_modP{2}.p{1,1,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{1,1,w},results.disp_modP{3}.p{1,1,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,2), plot(results.FD_UP{1, 2}{1, 1},results.FD_UP{1, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on,
                plot(results.avg_FD_UP{1,2}(1,:),results.avg_FD_UP{1,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 2},results.disp_mod{1, 2},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{1,2,w},results.disp_modP{1}.p{1,2,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{1,2,w},results.disp_modP{2}.p{1,2,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{1,2,w},results.disp_modP{3}.p{1,2,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,3), plot(results.FD_UP{1, 3}{1, 1},results.FD_UP{1, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{1,3}(1,:),results.avg_FD_UP{1,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{1, 3},results.disp_mod{1, 3},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{1,3,w},results.disp_modP{1}.p{1,3,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{1,3,w},results.disp_modP{2}.p{1,3,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{1,3,w},results.disp_modP{3}.p{1,3,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,4), plot(results.FD_UP{2, 1}{1, 1},results.FD_UP{2, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,1}(1,:),results.avg_FD_UP{2,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 1},results.disp_mod{2, 1},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{2,1,w},results.disp_modP{1}.p{2,1,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{2,1,w},results.disp_modP{2}.p{2,1,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{2,1,w},results.disp_modP{3}.p{2,1,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,5), plot(results.FD_UP{2, 2}{1, 1},results.FD_UP{2, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,2}(1,:),results.avg_FD_UP{2,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 2},results.disp_mod{2, 2},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{2,2,w},results.disp_modP{1}.p{2,2,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{2,2,w},results.disp_modP{2}.p{2,2,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{2,2,w},results.disp_modP{3}.p{2,2,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,6), plot(results.FD_UP{2, 3}{1, 1},results.FD_UP{2, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{2,3}(1,:),results.avg_FD_UP{2,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{2, 3},results.disp_mod{2, 3},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{2,3,w},results.disp_modP{1}.p{2,3,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{2,3,w},results.disp_modP{2}.p{2,3,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{2,3,w},results.disp_modP{3}.p{2,3,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,7), plot(results.FD_UP{3, 1}{1, 1},results.FD_UP{3, 1}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,1}(1,:),results.avg_FD_UP{3,1}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 1},results.disp_mod{3, 1},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{3,1,w},results.disp_modP{1}.p{3,1,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{3,1,w},results.disp_modP{2}.p{3,1,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{3,1,w},results.disp_modP{3}.p{3,1,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,8), plot(results.FD_UP{3, 2}{1, 1},results.FD_UP{3, 2}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,2}(1,:),results.avg_FD_UP{3,2}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 2},results.disp_mod{3, 2},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{3,2,w},results.disp_modP{1}.p{3,2,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{3,2,w},results.disp_modP{2}.p{3,2,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{3,2,w},results.disp_modP{3}.p{3,2,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
subplot(3,3,9), plot(results.FD_UP{3, 3}{1, 1},results.FD_UP{3, 3}{3, 1},'Color',[0.7 0.7 0.7],'LineWidth',1), hold on, 
                plot(results.avg_FD_UP{3,3}(1,:),results.avg_FD_UP{3,3}(3,:),'r','LineWidth',2), hold on
                plot(results.tt_mod{3, 3},results.disp_mod{3, 3},'b--','LineWidth',2), hold on
                for w = wind_start:length(results.wind_v)
                    plot(results.tt_modP{1}.p{3,3,w},results.disp_modP{1}.p{3,3,w},'g--','LineWidth',1), hold on
                    plot(results.tt_modP{2}.p{3,3,w},results.disp_modP{2}.p{3,3,w},'y--','LineWidth',1), hold on
                    plot(results.tt_modP{3}.p{3,3,w},results.disp_modP{3}.p{3,3,w},'c--','LineWidth',1), hold on
                end
                grid on
                ylim([d_minn d_maxx])
sgtitle('Displacement [m] in Unperturbed Ballistic Release')


% Perturbed on Force Hold
f_minn = -3;
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 1}.p{1, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,1}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,2), plot(results.FD_P{1, 1}.p{1, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,2}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,3), plot(results.FD_P{1, 1}.p{1, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,3}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,4), plot(results.FD_P{1, 1}.p{2, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,1}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,5), plot(results.FD_P{1, 1}.p{2, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,2}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,6), plot(results.FD_P{1, 1}.p{2, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,3}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,7), plot(results.FD_P{1, 1}.p{3, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,1}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,8), plot(results.FD_P{1, 1}.p{3, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,2}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,9), plot(results.FD_P{1, 1}.p{3, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,3}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
sgtitle('Force [N] in Pulse @ Force Hold')


d_minn = -0.5;
d_maxx = -0.35;

figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 1}.p{1, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,2), plot(results.FD_P{1, 1}.p{1, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,3), plot(results.FD_P{1, 1}.p{1, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{1,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,4), plot(results.FD_P{1, 1}.p{2, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,5), plot(results.FD_P{1, 1}.p{2, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,6), plot(results.FD_P{1, 1}.p{2, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{2,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,7), plot(results.FD_P{1, 1}.p{3, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,8), plot(results.FD_P{1, 1}.p{3, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,9), plot(results.FD_P{1, 1}.p{3, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{1}.p{3,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
sgtitle('Displacement [m] in Pulse @ Force Hold')



% Perturbed on Release
f_minn = -3;
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 2}.p{1, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,1}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,2), plot(results.FD_P{1, 2}.p{1, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,2}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,3), plot(results.FD_P{1, 2}.p{1, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,3}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,4), plot(results.FD_P{1, 2}.p{2, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,1}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,5), plot(results.FD_P{1, 2}.p{2, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,2}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,6), plot(results.FD_P{1, 2}.p{2, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,3}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,7), plot(results.FD_P{1, 2}.p{3, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,1}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,8), plot(results.FD_P{1, 2}.p{3, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,2}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
subplot(3,3,9), plot(results.FD_P{1, 2}.p{3, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,3}(1,:),'r','LineWidth',2), grid on
                ylim([f_minn f_maxx])
sgtitle('Force [N] in Pulse @ Relsease')

d_minn = -0.5;
d_maxx = -0.35;

figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 2}.p{1, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,2), plot(results.FD_P{1, 2}.p{1, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,3), plot(results.FD_P{1, 2}.p{1, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{1,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,4), plot(results.FD_P{1, 2}.p{2, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,5), plot(results.FD_P{1, 2}.p{2, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,6), plot(results.FD_P{1, 2}.p{2, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{2,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,7), plot(results.FD_P{1, 2}.p{3, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,8), plot(results.FD_P{1, 2}.p{3, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,9), plot(results.FD_P{1, 2}.p{3, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{2}.p{3,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
sgtitle('Displacement [m] in Pulse @ Release')

% Perturbed on Pos Hold
figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 3}.p{1, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,2), plot(results.FD_P{1, 3}.p{1, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,3), plot(results.FD_P{1, 3}.p{1, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,4), plot(results.FD_P{1, 3}.p{2, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,5), plot(results.FD_P{1, 3}.p{2, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,6), plot(results.FD_P{1, 3}.p{2, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,3}(1,:),'r','LineWidth',2), grid on
subplot(3,3,7), plot(results.FD_P{1, 3}.p{3, 1}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,1}(1,:),'r','LineWidth',2), grid on
subplot(3,3,8), plot(results.FD_P{1, 3}.p{3, 2}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,2}(1,:),'r','LineWidth',2), grid on
subplot(3,3,9), plot(results.FD_P{1, 3}.p{3, 3}{1, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,3}(1,:),'r','LineWidth',2), grid on
sgtitle('Force [N] in Pulse @ Position Hold')

d_minn = -14e-3;
d_maxx = 0.5e-3;

figure(),
set(gcf,'color','w');
subplot(3,3,1), plot(results.FD_P{1, 3}.p{1, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,2), plot(results.FD_P{1, 3}.p{1, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,3), plot(results.FD_P{1, 3}.p{1, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{1,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,4), plot(results.FD_P{1, 3}.p{2, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,5), plot(results.FD_P{1, 3}.p{2, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,6), plot(results.FD_P{1, 3}.p{2, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{2,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,7), plot(results.FD_P{1, 3}.p{3, 1}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,1}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,8), plot(results.FD_P{1, 3}.p{3, 2}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,2}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
subplot(3,3,9), plot(results.FD_P{1, 3}.p{3, 3}{2, 1}','Color',[0.7 0.7 0.7],'LineWidth',1), hold on,  
                plot(results.avg_FD_P{3}.p{3,3}(2,:),'r','LineWidth',2), grid on
                ylim([d_minn d_maxx])
sgtitle('Displacement [m] in Pulse @ Position Hold')

%% FREQUENCY PLOTS
%Bode Plots
options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc.
figure(),
set(gcf,'color','w');
subplot(3,3,1), bodemag(results.TF.up{1, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{1,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{1,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{1,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,2), bodemag(results.TF.up{1, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{1,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{1,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{1,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,3), bodemag(results.TF.up{1, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{1,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{1,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{1,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,4), bodemag(results.TF.up{2, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{2,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{2,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{2,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,5), bodemag(results.TF.up{2, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{2,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{2,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{2,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,6), bodemag(results.TF.up{2, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{2,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{2,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{2,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,7), bodemag(results.TF.up{3, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{3,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{3,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{3,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,8), bodemag(results.TF.up{3, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{3,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{3,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{3,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,9), bodemag(results.TF.up{3, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodemag(results.TF.p{3,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodemag(results.TF.p{3,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodemag(results.TF.p{3,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
sgtitle('Bode Plot - Magnitude - Identified Models')


options.MagVisible = 'off';
figure(),
set(gcf,'color','w');
subplot(3,3,1), bodeplot(results.TF.up{1, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{1,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{1,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{1,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,2), bodeplot(results.TF.up{1, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{1,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{1,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{1,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,3), bodeplot(results.TF.up{1, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{1,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{1,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{1,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,4), bodeplot(results.TF.up{2, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{2,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{2,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{2,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,5), bodeplot(results.TF.up{2, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{2,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{2,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{2,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,6), bodeplot(results.TF.up{2, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{2,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{2,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{2,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,7), bodeplot(results.TF.up{3, 1},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{3,1,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{3,1,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{3,1,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,8), bodeplot(results.TF.up{3, 2},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{3,2,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{3,2,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{3,2,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
subplot(3,3,9), bodeplot(results.TF.up{3, 3},{0.1*2*pi,10*2*pi},'b',options), hold on,
                for w = wind_start:length(results.wind_v)
                    bodeplot(results.TF.p{3,3,1,w},{0.1*2*pi,10*2*pi},options,'g--'), hold on
                    bodeplot(results.TF.p{3,3,2,w},{0.1*2*pi,10*2*pi},options,'y--'), hold on
                    bodeplot(results.TF.p{3,3,3,w},{0.1*2*pi,10*2*pi},options,'c--'), hold on
                end
                grid on
sgtitle('Bode Plot - Phase - Identified Models')



