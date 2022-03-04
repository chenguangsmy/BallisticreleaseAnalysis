% as 1) The stiffness measured from the forceTransducer and OPTOTRAK
% greater than the marked value; 2) The stiffness measured from the
% OPTOTRAK is on the marked value. I'm guessing the forceTransducer read
% the force bigger than the actually marked value.  

% Now I'm testing it use standard mass. The g is 9.80665???

% % % fttmp = SessionScanFT(3957,'KingKongFT03957.forceTest500g.csv');
fttmp.plotForceOrigin()

t_shift = fttmp.elapse - fttmp.elapse(1);
force_org = fttmp.force_origin;
plot(t_shift, force_org);

t_unload = [8 11; 20 23; 32 35; 48 51; 65 68; 82 85];
t_load = [1 6; 15 18; 26 30; 40 46; 55 61; 74 79];

F_measures = zeros(1,size(t_load,1));
for epoc_i = 1:size(t_load,1)
    t_idx_load = t_shift > t_load(epoc_i,1) & t_shift < t_load(epoc_i,2);
    t_idx_unload = t_shift > t_unload(epoc_i,1) & t_shift < t_unload(epoc_i,2);
    Favg_load = mean(force_org(:,t_idx_load),2);
    Favg_unload = mean(force_org(:,t_idx_unload),2);  
    ifplot = 1;
    if (ifplot) 
        clf; hold on;
        plot(t_shift, force_org);
        plot(t_shift(:,t_idx_load), force_org(:,t_idx_load), 's');
        plot(t_shift(:,t_idx_unload), force_org(:,t_idx_unload), 'o');
    end
    
    dF = sqrt(sum(Favg_load - Favg_unload).^2);
    F_measures(epoc_i) = dF;
end
%%
mass_measures = F_measures/9.8 * 1000

%% multiple sessions 
% +y direction, with bearings 
mass_list = [200 500 700 966 1166 1395 1595];
fnames = {  'KingKongFT03958.forceTest200g.csv';
            'KingKongFT03958.forceTest500g.csv';
            'KingKongFT03958.forceTest700g.csv';
            'KingKongFT03958.forceTest966g.csv';
            'KingKongFT03958.forceTest1166g.csv';
            'KingKongFT03958.forceTest1395g.csv';
            'KingKongFT03958.forceTest1595g.csv'};
t_unloadCell = {[ 8 10; 18 22; 32 34; 48 50; 66 70; 82 84;];
                [ 5  9; 17 21; 32 36; 48 52; 64 68; 80 84;];
                [ 6  8; 16 18; 31 34; 47 51; 64 67; 81 85;];
                [16 19; 32 34; 48 50; 64 68; 80 83; 96 99;];
                [ 7  8; 16 19; 32 35; 48 52; 66 68; 81 84;];
                [16 19; 33 35; 48 50; 64 66; 83 85; 96 98;];
                [18 20; 32 34; 51 53; 65 68; 82 86; 98 101;];
            };
t_loadCell = {  [ 2  4; 13 15; 24 28; 38 44; 54 60; 74 78;];
                [ 1  3; 12 14; 26 30; 40 44; 56 60; 72 76;];
                [ 1  3; 12 14; 22 25; 38 44; 57 61; 72 77;];
                [11 14; 22 28; 37 44; 54 62; 71 78; 86 94;];
                [ 1  3; 11 14; 24 29; 38 45; 55 62; 72 78;];
                [12 14; 24 28; 40 44; 54 60; 72 78; 88 94;];
                [12 15; 24 29; 38 45; 56 62; 72 78; 90 94;];
            };
        

F_measuresC = zeros(length(fnames), 6); % 6 epocs
figure('position', [0 0 300 300]);
for file_i = 1:length(fnames)
fttmp = SessionScanFT(3957,fnames{file_i});
fttmp.plotForceOrigin()

t_shift = fttmp.elapse - fttmp.elapse(1);
force_org = fttmp.force_origin;
plot(t_shift, force_org);
close all;
figure('position', [0 800 300 300]);
t_unload = t_unloadCell{file_i};
t_load = t_loadCell{file_i};

F_measures = zeros(1,size(t_load,1));
for epoc_i = 1:size(t_load,1)
    t_idx_load = t_shift > t_load(epoc_i,1) & t_shift < t_load(epoc_i,2);
    t_idx_unload = t_shift > t_unload(epoc_i,1) & t_shift < t_unload(epoc_i,2);
    Favg_load = mean(force_org(:,t_idx_load),2);
    Favg_unload = mean(force_org(:,t_idx_unload),2);
    ifplot = 1;
    
    if (ifplot)
        
        clf; hold on;
        plot(t_shift, force_org);
        plot(t_shift(:,t_idx_load), force_org(:,t_idx_load), 's');
        plot(t_shift(:,t_idx_unload), force_org(:,t_idx_unload), 'o');
    end
    
    dF = sqrt(sum((Favg_load - Favg_unload).^2))
    F_measures(epoc_i) = dF;
%     F_measuresC{file_i} = F_measures/9.8;   % assume g = 9.8
end
    F_measuresC(file_i,:) = F_measures/9.8*1000;
end
% plot 
figure('unit', 'inch', 'position',[0 0 4 4]); hold on;
% for file_i = 1:length(fnames)
%     mass = mass_list(file_i);
%     plot(mass, F_measuresC{file_i}*1000, '.', 'MarkerSize', 10); 
% end
x = ones(6,1)*mass_list;
y = F_measuresC';
lh1 = scatter(x(:), y(:),  10); 
% refline();
lh2 = plot([0 2500], [0 2500], 'Color', [0.8 0.8 0.8], 'LineWidth', 1); 
xlim([0 2500]); ylim([0 2500]);
% yline(157.6 -8.8); yline(157.6 + 8.8);
legend([lh1, lh2], {'measured weight', 'known weight'});
xlabel('known weight (g)');
ylabel('measured weight (g)');
title('force transducer measure known mass X');
        

%% multiple sessions 
% -y direction, with bearings 
mass_list = [200 500 700 966 1166 1395 1595];
fnames = {  'KingKongFT03959.forceTest200g_y.csv';
            'KingKongFT03959.forceTest500g_y.csv';
            'KingKongFT03959.forceTest700g_y.csv';
            'KingKongFT03959.forceTest966g_y.csv';
            'KingKongFT03959.forceTest1166g_y.csv';
            'KingKongFT03959.forceTest1395g_y.csv';
            'KingKongFT03959.forceTest1595g_y.csv'};
t_unloadCell = {[16 20; 32 36; 48 54; 66 70; 82 84; 96 98];
                [16 18; 32 34; 48 50; 64 68; 80 84; 98 100];
                [18 20; 34 38; 50 54; 66 68; 81 84; 97 99];
                [17 19; 32 34; 48 50; 65 68; 82 84; 97 99];
                [17 19; 32 34; 48 50; 65 68; 80 82; 97 99];
                [17 19; 32 34; 48 50; 65 67; 81 84; 97 99];
                [17 19; 33 36; 48 50; 65 67; 80 83; 97 99];
            };
t_loadCell = {  [12 14; 20 30; 40 46; 58 62; 76 78; 90 94];
                [10 14; 25 30; 40 46; 54 60; 72 78; 88 94];
                [10 13; 24 30; 40 46; 57 62; 72 76; 88 95];
                [10 13; 24 30; 40 46; 57 62; 72 76; 88 95];
                [10 13; 24 30; 40 46; 57 62; 72 76; 88 94];
                [10 14; 24 30; 40 46; 57 62; 72 78; 88 94];
                [10 14; 24 30; 40 46; 57 62; 72 78; 88 95];
            };
% % % fttmp = SessionScanFT(3957,fnames{7});
% % fttmp.plotForceOrigin()
% close all; 
% figure('position', [0 800 800 200])
% t_shift = fttmp.elapse - fttmp.elapse(1);
% force_org = fttmp.force_origin;
% plot(t_shift, force_org);
% % %%        

F_measuresC = zeros(length(fnames), 6); % 6 epocs
figure('position', [0 0 300 300]);
for file_i = 1:length(fnames)
fttmp = SessionScanFT(3957,fnames{file_i});
fttmp.plotForceOrigin()

t_shift = fttmp.elapse - fttmp.elapse(1);
force_org = fttmp.force_origin;
plot(t_shift, force_org);
close all;
figure('position', [0 800 300 300]);
t_unload = t_unloadCell{file_i};
t_load = t_loadCell{file_i};

F_measures = zeros(1,size(t_load,1));
for epoc_i = 1:size(t_load,1)
    t_idx_load = t_shift > t_load(epoc_i,1) & t_shift < t_load(epoc_i,2);
    t_idx_unload = t_shift > t_unload(epoc_i,1) & t_shift < t_unload(epoc_i,2);
    Favg_load = mean(force_org(:,t_idx_load),2);
    Favg_unload = mean(force_org(:,t_idx_unload),2);
    ifplot = 1;
    
    if (ifplot)
        
        clf; hold on;
        plot(t_shift, force_org);
        plot(t_shift(:,t_idx_load), force_org(:,t_idx_load), 's');
        plot(t_shift(:,t_idx_unload), force_org(:,t_idx_unload), 'o');
    end
    
    dF = sqrt(sum((Favg_load - Favg_unload).^2))
    F_measures(epoc_i) = dF;
%     F_measuresC{file_i} = F_measures/9.8;   % assume g = 9.8
end
    F_measuresC(file_i,:) = F_measures/9.8*1000;
end
% plot 
figure('unit', 'inch', 'position',[0 0 4 4]); hold on;
% for file_i = 1:length(fnames)
%     mass = mass_list(file_i);
%     plot(mass, F_measuresC{file_i}*1000, '.', 'MarkerSize', 10); 
% end
x = ones(6,1)*mass_list;
y = F_measuresC';
lh1 = scatter(x(:), y(:),  10); 
% refline();
lh2 = plot([0 2500], [0 2500], 'Color', [0.8 0.8 0.8], 'LineWidth', 1); 
xlim([0 2500]); ylim([0 2500]);
% yline(157.6 -8.8); yline(157.6 + 8.8);
legend([lh1, lh2], {'measured weight', 'known weight'});
xlabel('known weight (g)');
ylabel('measured weight (g)');
title('force transducer measure known mass -Y');
        
        
%% multiple sessions 
% -z direction 
mass_list = [200 500 700 966 1166 1395 1595];
fnames = {  'KingKongFT03957.forceTest200g.csv';
            'KingKongFT03957.forceTest500g.csv';
            'KingKongFT03957.forceTest700g.csv';
            'KingKongFT03957.forceTest966g.csv';
            'KingKongFT03957.forceTest1166g.csv';
            'KingKongFT03957.forceTest1395g.csv';
            'KingKongFT03957.forceTest1595g.csv'};
t_unloadCell = {[15 20; 32 36; 48 52; 64 68; 80 85; 98 102];
                [ 8 11; 20 23; 32 35; 48 51; 65 68; 82 85];
                [18 20; 32 34; 48 50; 64 68; 80 84; 98 102];
                [16 18; 32 36; 48 52; 66 70; 82 86; 98 102];
                [ 4  6; 18 22; 32 34; 49 52; 64 67; 78 82];
                [ 4  7; 17 20; 33 36; 50 53; 85 68; 82 85];
                [16 20; 32 36; 50 52; 66 68; 82 84; 98 100];
            };
t_loadCell = {  [10 13; 24 28; 40 44; 55 60; 72 76; 90 94];
                [ 1  6; 15 18; 26 30; 40 46; 55 61; 74 79];
                [12 15; 24 28; 40 44; 54 60; 72 76; 88 94];
                [10 14; 24 28; 40 44; 56 62; 74 78; 92 94];
                [ 1  2; 10 15; 26 30; 40 45; 56 60; 70 74];
                [ 1  3; 11 13; 25 30; 41 46; 57 62; 72 78];
                [12 14; 24 30; 40 46; 56 62; 72 78; 88 94];
            };
%        
        

F_measuresC = zeros(length(fnames), 6); % 6 epocs
figure('position', [0 0 300 300]);
for file_i = 1:length(fnames)
fttmp = SessionScanFT(3957,fnames{file_i});
fttmp.plotForceOrigin()

t_shift = fttmp.elapse - fttmp.elapse(1);
force_org = fttmp.force_origin;
plot(t_shift, force_org);
close all;
figure('position', [0 800 300 300]);
t_unload = t_unloadCell{file_i};
t_load = t_loadCell{file_i};

F_measures = zeros(1,size(t_load,1));
for epoc_i = 1:size(t_load,1)
    t_idx_load = t_shift > t_load(epoc_i,1) & t_shift < t_load(epoc_i,2);
    t_idx_unload = t_shift > t_unload(epoc_i,1) & t_shift < t_unload(epoc_i,2);
    Favg_load = mean(force_org(:,t_idx_load),2);
    Favg_unload = mean(force_org(:,t_idx_unload),2);
    ifplot = 1;
    
    if (ifplot)
        
        clf; hold on;
        plot(t_shift, force_org);
        plot(t_shift(:,t_idx_load), force_org(:,t_idx_load), 's');
        plot(t_shift(:,t_idx_unload), force_org(:,t_idx_unload), 'o');
    end
    
%     dF = sqrt(sum((Favg_load - Favg_unload).^2))
    dF = sqrt(sum(Favg_load - Favg_unload).^2)
    F_measures(epoc_i) = dF;
%     F_measuresC{file_i} = F_measures/9.8;   % assume g = 9.8
end
    F_measuresC(file_i,:) = F_measures/9.8*1000;
end
% plot 
figure('unit', 'inch', 'position',[0 0 4 4]); hold on;
% for file_i = 1:length(fnames)
%     mass = mass_list(file_i);
%     plot(mass, F_measuresC{file_i}*1000, '.', 'MarkerSize', 10); 
% end
x = ones(6,1)*mass_list;
y = F_measuresC';
lh1 = scatter(x(:), y(:),  10); 
% refline();
lh2 = plot([0 2500], [0 2500], 'Color', [0.8 0.8 0.8], 'LineWidth', 1); 
xlim([0 2500]); ylim([0 2500]);
% yline(157.6 -8.8); yline(157.6 + 8.8);
legend([lh1, lh2], {'measured weight', 'known weight'});
xlabel('known weight (g)');
ylabel('measured weight (g)');
title('force transducer measure known mass Y');