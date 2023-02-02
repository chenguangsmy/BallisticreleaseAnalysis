% Federico reported the marker 2 and 3 shifted position during subject 1 x
% and y movements. Let me check for that 

sstmp1 = SessionScan(4401); % x
sstmp1 = SessionScan(4427); % x
sstmp2 = SessionScan(4404); % y
sstmp2 = SessionScan(4405); % y

figure(); hold on;
plot(sstmp1.data.t, sstmp1.data.ox(1,:,3)); 
plot(sstmp1.data.t, sstmp1.data.ox(1,:,4)); 
title('marker 3, 4 in x movements ')

figure(); hold on;
plot(sstmp2.data.t, sstmp2.data.ox(1,:,3)); 
plot(sstmp2.data.t, sstmp2.data.ox(1,:,4)); 
title('marker 3, 4 in y movements ')

figure(); hold on;  % x 
plot3(sstmp1.data.ox(1,:,3), sstmp1.data.ox(2,:,3), sstmp1.data.ox(3,:,3)); 
plot3(sstmp1.data.ox(1,:,4), sstmp1.data.ox(2,:,4), sstmp1.data.ox(3,:,4)); 

figure(); hold on; % y
plot3(sstmp2.data.ox(1,:,3), sstmp2.data.ox(2,:,3), sstmp2.data.ox(3,:,3)); 
plot3(sstmp2.data.ox(1,:,4), sstmp2.data.ox(2,:,4), sstmp2.data.ox(3,:,4)); 

%% after export 
% Check data again if there is the marker mistake 
load(['/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/data/' ...
    'processedData/ss4379_4438']);
close all;
for subj_i = 6
    for dir_i = 1:4
        figure(); hold on;
        for fce_i = 1:3
            for dist_i = 1:3
                for trial_i = 1:9 
                    tr = data{subj_i,dir_i,fce_i,dist_i,trial_i,1};
                    plot3(tr.ox(1,:,3), tr.ox(2,:,3), tr.ox(3,:,3), 'b');
                    plot3(tr.ox(1,:,4), tr.ox(2,:,4), tr.ox(3,:,4), 'r');
                end
            end
        end
    end
end
           