%% Conduct the perturbation magnitude checking in the high-stiffness mode,
% showing that the high-stiffness mode cannot yield position differences
% front:

%--------------------------------------------------------------------------
% | KingKong.02256.mat	| KingKongFT02256.csv	| KingKongWAM02256.csv	| randomize step pert, 3N 5cm release. |  
% | KingKong.02268.mat	| KingKongFT02268.csv	| KingKongWAM02268.csv	| randomize step pert, 3N 7.5cm release. | 
% | KingKong.02262.mat	| KingKongFT02262.csv	| KingKongWAM02262.csv	| randomize step pert, 3N 10cm release. |   
% | KingKong.02257.mat	| KingKongFT02257.csv	| KingKongWAM02257.csv	| randomize step pert, 9N 5cm release. |  
% | KingKong.02259.mat	| KingKongFT02259.csv	| KingKongWAM02259.csv	| randomize step pert, 9N 7.5cm release. |  
% | KingKong.02266.mat	| KingKongFT02266.csv	| KingKongWAM02266.csv	| randomize step pert, 9N 10cm release. |  
% | KingKong.02263.mat	| KingKongFT02263.csv	| KingKongWAM02263.csv	| randomize step pert, 15N 7.5cm release. | ??? no pert data! 
% | KingKong.02258.mat	| KingKongFT02258.csv	| KingKongWAM02258.csv	| randomize step pert, 15N 5cm release. |  
% | KingKong.02267.mat	| KingKongFT02267.csv	| KingKongWAM02267.csv	| randomize step pert, 15N 10cm release. |  
% | KingKong.02264.mat	| KingKongFT02264.csv	| KingKongWAM02264.csv	| randomize step pert, 21N 5cm release. |  
% | KingKong.02265.mat	| KingKongFT02265.csv	| KingKongWAM02265.csv	| randomize step pert, 21N 7.5cm release. |  
% | KingKong.02261.mat	| KingKongFT02261.csv	| KingKongWAM02261.csv	| randomize step pert, 21N 10cm release. |  
%--------------------------------------------------------------------------
%   5cm: 2256 2257 2258 2264
% 7.5cm: 2268 2259 2263 2265
%  10cm: 2262 2266 2267 2261

% back:
% | KingKong.2355.mat	| KingKongFT2355.csv	| KingKongWAM2355.csv	| backward ballistic release with pert,   5cm   3N |
% | KingKong.2358.mat	| KingKongFT2358.csv	| KingKongWAM2358.csv	| backward ballistic release with pert, 7.5cm   9N |
% | KingKong.2359.mat	| KingKongFT2359.csv	| KingKongWAM2359.csv	| backward ballistic release with pert,  10cm   15N |
% | KingKong.2360.mat	| KingKongFT2360.csv	| KingKongWAM2360.csv	| backward ballistic release with pert,  10cm   21N |
% | KingKong.2361.mat	| KingKongFT2361.csv	| KingKongWAM2361.csv	| backward ballistic release with pert, 7.5cm   21N |
% | KingKong.2362.mat	| KingKongFT2362.csv	| KingKongWAM2362.csv	| backward ballistic release with pert, 7.5cm   15N |
% | KingKong.2363.mat	| KingKongFT2363.csv	| KingKongWAM2363.csv	| backward ballistic release with pert,  10cm   9N |
% | KingKong.2365.mat	| KingKongFT2365.csv	| KingKongWAM2365.csv	| backward ballistic release with pert, 7.5cm   3N |
% | KingKong.2366.mat	| KingKongFT2366.csv	| KingKongWAM2366.csv	| backward ballistic release with pert,   5cm   9N |
% | KingKong.2367.mat	| KingKongFT2367.csv	| KingKongWAM2367.csv	| backward ballistic release with pert,   5cm   15N |
% | KingKong.2368.mat	| KingKongFT2368.csv	| KingKongWAM2368.csv	| backward ballistic release with pert,   5cm   21N |
% | KingKong.2369.mat	| KingKongFT2369.csv	| KingKongWAM2369.csv	| backward ballistic release with pert,  10cm   3N |
%   5cm: 2355 2366 2367 2368
% 7.5cm: 2365 2358 2362 2361
%  10cm: 2369 2363 2359 2360
sessions_mat_f = [2256 2257 2258 2264; 2268 2259 2263 2265; 2262 2266 2267 2261]; % front 
sessions_mat_b = [2355 2366 2367 2368; 2365 2358 2362 2361; 2369 2363 2359 2360]; % back
sessions_all = [sessions_mat_b(:); sessions_mat_f(:)];
%sessions_all = [2256 2257 2258 2264 2268 2259 2263 2265 2262 2266 2267 2261];
%sessions_mat_b = [2457 2458 2459; 2460 2461 2462; 2455 2463 2456];
for session_i = 1:length(sessions_all)
    sessions_idx = sessions_all(session_i);
    eval(['ss' num2str(sessions_idx) ' = SessionScan(' num2str(sessions_idx) ');']);
end

%%  plot stacked figure;
axh1 = figure();
col_array = colormap(lines);
%force_arr = [3 9 15 21];
force_arr = [3 9 15];
for dist_i = 1:4
%for fce = 1:3
%sessions_allf = sessions_mat_f(:,dist_i)'; % 3, 9, 15, 21N
sessions_allb = sessions_mat_b(:,dist_i)'; % 3, 9, 15, 21N
%sessions_allf = sessions_mat_f(fce,:); % 5 ,7.5, 10cm
%sessions_allb = sessions_mat_b(fce,:); % 5 ,7.5, 10cm
%axh = figure();
axh = subplot(1, 3, fce);
for session_i = 1:length(sessions_allb)
    sessions_idx = sessions_allb(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
ylabel(['force:' num2str(force_arr(fce)) 'N']);
if dist_i == 1
    title('back')
elseif dist_i == 4
    xlabel('time (s)');
end
%axh = figure();
%axh = subplot(4, 2, (dist_i-1)*2+2);
% for session_i = 1:length(sessions_allb)
%     sessions_idx = sessions_allb(session_i);
%     eval(['sstmp = ss' num2str(sessions_idx) ';']);
%     
%     plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
% end
% if dist_i == 1
%     title('back')
% elseif dist_i == 4
%     xlabel('time (s');
% end
end

axh2 = figure();
tar_arr = [5 7.5 10];
for target = 1:3
sessions_allf = sessions_mat_f(target,:); % 5 ,7.5, 10cm
sessions_allb = sessions_mat_b(target,:); % 5 ,7.5, 10cm
%axh = subplot(3, 2, (target-1)*2+1);
axh = subplot(1, 3, target);
for session_i = 1:length(sessions_allf)
    sessions_idx = sessions_allf(session_i);
    eval(['sstmp = ss' num2str(sessions_idx) ';']);
    plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
end
%ylabel(['target:' num2str(tar_arr(target)) 'cm']);
title(['target:' num2str(tar_arr(target)) 'cm']);
% if target == 1
%     title('front')
% elseif target == 3
%     xlabel('time (s');
% end
%axh = figure();
% axh = subplot(3, 2, (target-1)*2+2);
% for session_i = 1:length(sessions_allb)
%     sessions_idx = sessions_allb(session_i);
%     eval(['sstmp = ss' num2str(sessions_idx) ';']);
%     
%     plotStepPertResponse_raw(sstmp, axh, col_array(session_i,:));
% end
% if target == 1
%     title('back')
% elseif target == 3
%     xlabel('time (s');
% end
end
