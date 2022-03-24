%% optotrak spring stiffness check 
% As the readout from the optotrak and the force transducer is not right,
% this time I began to do the optotrak with the spring. I added mass and
% see the displacement. As mass is known, the optotrak should be able to
% get the exact stiffness of the spring.  

%% optdat = SessionScanOPT(3957, 'KingKongOPT03957.500gmassTestSpring.csv');

optdat = SessionScanOPT(3957, fname{2});
figure();
plot(optdat.datah.t, [optdat.datah.x(1,:);optdat.datah.y(1,:);optdat.datah.z(1,:)]);

t_range = [22, 28];
t_ref_idx = optdat.datah.t > t_range(1,1) & optdat.datah.t<t_range(1,2); 
pos_ref = nanmean([optdat.datah.x(1,t_ref_idx);optdat.datah.y(1,t_ref_idx);optdat.datah.z(1,t_ref_idx)],2);
optdat_dist = sqrt(sum(([optdat.datah.x(1,:);optdat.datah.y(1,:);optdat.datah.z(1,:)] - pos_ref) .^2 ));
F_mass = 4.9; % N, == 500g mass

figure(); 
subplot(2,1,1); title('displacement');
plot(optdat.datah.t, optdat_dist);
xlabel('t'); ylabel('m');
subplot(2,1,2); title('stiffness');
plot(optdat.datah.t, F_mass./optdat_dist);
ylim([0 1000]);
xlabel('t'); ylabel('N/m');



%% measure more data 
%% The files and descriptions 
% 
fname{1} = 'KingKongOPT03957.200gmass.csv'; 
t_zone{1,1} = [22 26; 72 78; 108 112; 138 143]; % uptime
t_zone{1,2} = [10 13; 55 60; 90 95;   122 127]; % downtime 

fname{2} = 'KingKongOPT03957.500gmassTestSpring.csv'; 
t_zone{2,1} = [18 22; 42 46;64 70; 102 110; 150 155]; % uptime
t_zone{2,2} = [8 12; 32 36; 54 60; 86 92;  138 142 ]; % downtime 
% : I only use a 500g mass holding up, and down. 

fname{3} = 'KingKongOPT03957.700gmass.csv'; 
t_zone{3,1} = [20 25; 61 65; 90 95; 125 130]; % uptime
t_zone{3,2} = [10 14; 40 45; 70 75; 110 115]; % downtime 
%: I use a 500g mass and on/off 200g mass

fname{4} = 'KingKongOPT03957.1000gmass.csv'; 
t_zone{4,1} = [25 28; 60 65; 100 105;  135 140; 170 175];
t_zone{4,2} = [10 13; 40 45; 80 85;    115 120; 155 160]; 
% : I use a both 500g mass and the unknown ? 500g mass. (male mass)
fname{5} = 'KingKongOPT03957.1500gmass.csv';
t_zone{5,1} = [17 21;42 48; 68 72; 92 96;120 124];
t_zone{5,2} = [5 10; 30 34; 56 60; 80 84;108 112]; 
% 
% 
fname{6} = 'KingKongOPT03957.1700gmass.csv';
t_zone{6,1} = [22 26;  52 58;  90 95;  135 140; ];
t_zone{6,2} = [6 12;   35 40;  70 75;  110 115; ]; 
% 
F_mass_list = [1.96 4.9, 1.96,4.9, 4.9, 4.9];

%%%%%%%
stiffness_measurements = nan(5,6);

for file_i = 1:6
    optdat = SessionScanOPT(3957, fname{file_i});
    for measure_i = 1:size(t_zone{file_i,1})
        t_range_hold = t_zone{file_i,1}(measure_i,:);
        t_range_measure = t_zone{file_i,2}(measure_i,:);
        t_ref_idx = optdat.datah.t > t_range_hold(1,1) & optdat.datah.t<t_range_hold(1,2);
        t_msu_idx = optdat.datah.t > t_range_measure(1,1) & optdat.datah.t<t_range_measure(1,2);
        pos_ref = nanmean([optdat.datah.x(1,t_ref_idx);optdat.datah.y(1,t_ref_idx);optdat.datah.z(1,t_ref_idx)],2);
        optdat_dist = sqrt(sum(([optdat.datah.x(1,:);optdat.datah.y(1,:);optdat.datah.z(1,:)] - pos_ref) .^2 ));
        F_mass = F_mass_list(file_i); % N, == 500g mass
        stiff_measure = F_mass./optdat_dist;
        
        ifplot = 1;
        if (ifplot)
            clf;
            subplot(2,1,1); title('displacement');
            plot(optdat.datah.t, optdat_dist);
            xlabel('t'); ylabel('m');
            subplot(2,1,2); title('stiffness'); hold on;
            plot(optdat.datah.t, F_mass./optdat_dist);
            plot(optdat.datah.t(t_msu_idx), F_mass./optdat_dist(t_msu_idx), '.', 'MarkerSize', 5);
            ylim([0 1000]);
            xlabel('t'); ylabel('N/m');
        end
        
        stiffness_measurements(measure_i, file_i) = nanmean(stiff_measure(t_msu_idx));
    end
end

% plot a nice figure here 
mass_offset = [200 500 700 1000 1500 1700];
fh = figure('unit', 'inch', 'position',[0 0 3 3]);
plot(mass_offset, stiffness_measurements, '.', 'MarkerSize', 10); 
xlim([0 2500]); ylim([120 320]);
yline(157.6 -8.8); yline(157.6 + 8.8);
xlabel('added mass (g)');
ylabel('measured stiffness (N/m)');
title('measured stiffness at different force');
saveas(fh, '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/sanityCheck_StiffnessMeasurement/stiffnessMeasurement_OPTOTRAK.png')