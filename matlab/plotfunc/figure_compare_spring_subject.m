function figure_compare_spring_subject()
    %% plot spring test when the xs0 and sr0 equilibrium at k0
    addpath('/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab');
% | 3N | 9N | 15N | 21N |
sessions_mat = [2611 2610 2612 2613;        % 160N/m
                2614 2615 2617 2619];       % 320N/m
force_list = [3 9 15 21]; % (N)
sessions_all = sessions_mat(:);
for session_i = 1:length(sessions_all)
     eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
         num2str(sessions_all(session_i)) ');']);
end
%
color_arr = colormap('lines');
figure();
% 1st line
for force_i = 1:size(force_list,2)
    axh = subplot(3,4,force_i); hold on;
    eval(['sstmp = ss' num2str(sessions_mat(1, force_i)) ';']);
    sstmp.plotStepPertResponse_raw(axh,color_arr(255,:)); % how to subtract the mean?
    title(['force ' num2str(force_list(force_i)) 'N']);
    set(gca, 'Ygrid', 'on');
    ylabel('');
    if force_i == 1
        ylabel('K=160 position(m)');
    end
    xlim([-0.1 1]);
    ylim([-0.05 0.01]); % posision
    xlabel('');
end
% 2nd line
for force_i = 1:size(force_list,2)
    axh = subplot(3,4,4+force_i); hold on;
    eval(['sstmp = ss' num2str(sessions_mat(2, force_i)) ';']);
    sstmp.plotStepPertResponse_raw(axh,color_arr(255,:)); % how to subtract the mean?
    title('');
    set(gca, 'Ygrid', 'on');
    ylabel('');
    if force_i == 1
        ylabel('K=320 position(m)');
    end
    xlim([-0.1 1]);
    ylim([-0.05 0.01]); % posision
    xlabel('');
end

% 3nd line
force_all = 3:6:21; 
sessions_mat = [[2468, 2457, 2458, 2459]; [2469, 2460, 2461, 2462]; [2470, 2456, 2463, 2455]; [2471, 2464, 2466, 2467]];
sessions_all = sessions_mat(:);
for session_i = 1:length(sessions_all)
     eval(['ss' num2str(sessions_all(session_i)) '= SessionScan(' ...
         num2str(sessions_all(session_i)) ');']);
end
for fce_i = 1:size(sessions_mat, 1)
    axh = subplot(3,4, fce_i+8); hold on; 
    for dist_i = 1:size(sessions_mat, 1)
        sstmp = eval(['ss' num2str(sessions_mat(fce_i, dist_i))]);
        axh = plotMeantrialPosPert(sstmp, axh, dist_i); % perturbation
    end
    ylim([-0.05 0.01]);
    xlim([-0.1 1]);
    xlabel('time'); 
    title('');
    ylabel('');
    if (fce_i == 1)
        ylabel('subject position(m)');
    end
    set(gca, 'Ygrid', 'on');
    legend off;
    %title([num2str(force_all(fce_i)) 'N']);
    if (fce_i == 4)
        label_handles = [axh.Children(8) axh.Children(6) axh.Children(4) axh.Children(2)];
        label_names = {'2.5cm', '5cm', '7.5cm', '10cm'};
        legend(label_handles, label_names);
    end
        
end

end