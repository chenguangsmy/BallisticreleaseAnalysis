function [nearestTrial,idx_nearestTrial] = getNearestTrial(pertTrial,otherTrials)
%GETNEARESTTRIAL Find the nearest trial of certain condition
%   [nearestTrial,idx_nearestTrial] = getNesrestTrial(pertTrial,otherTrials)
%   Finding nearest trial with certain conditions. 
%   Current condition is have the shortest displacement difference with the
%   unperturbed part. For the data before-release, I will choose the time 
%   limit of [-1, 0]s to be the comparing zone. 
%
%   InpArg: 
%    pertTrial: a 1-by-1 cell only contain the perturbed trial.
%    otherTrials: a n-by-1 cell contain unperturbed trials.
%   OurArg:
%    nearestTrial: a 1-by-1 cell contain the unperturbed trial which is the
%                   cloest to the perturbed trial.
%    idx_narestTrial: a integar indicating the position of the nearest 
%                       trial among the unperturbed trials.
%
%   This function is for the ballistic-release task, estimating stiffness
%   during movement. The place to use this function is when tidying data
%   into formmed data and save the perturbed and un-perturbed trials
%   correspondingly. The major problem this code trying to solve is: avoid
%   the steady-state error before perturbation. 
%   
% Author: Chenguang Zhang
% Date: 2022-03-24

% TODO: 
%   -1. Check the eligibility of pertTrial
    if ~isfield(pertTrial{1}, 'Fp')
        disp('Empty trial, ABORT!');
        nearestTrial = [];
        idx_nearestTrial = nan;
        return
    else
    pertSignal = pertTrial{1}.Fp; 
    if sum(sum(abs(pertSignal)))==0
        disp('No perturbation found in pertTrial, ABORT!');
        nearestTrial = [];
        idx_nearestTrial = nan;
        return
    end
    end
%   0. Specify the comparing needed variables
    t_range = [-1,1];       % comparing time start from -1s after release
    freq = 500;             % intropolate freq
    t_grids = t_range(1):1/freq:t_range(2);
    cpr_var = 'x'; 
    cpr_idx = 2;            % y-axis
    ifplot = 0;             % for debug check
    if (ifplot)
        fh_dbg = figure();
        set(fh_dbg, 'name', 'debug figure');
    else
        fh_dbg = [];
    end
%   1. Interpolate all the data (pertTrial and otherTrials) into one vector
%   and one matrix
    pTrial = pertTrial;
    npTrial= otherTrials;
    nump= 1; % should be 1
    numnp= length(otherTrials); 
    pData = zeros(nump, length(t_grids));
    pData_strobe = zeros(nump, length(t_grids)); % the perturb siangl
    npData = zeros(numnp, length(t_grids));
    
    % interpolate the perturbed trial, both x and Fp
        % deal with time 
    idx_releaset = (pTrial{1}.ts == 5 & diff([1 pTrial{1}.ts]) == 1); 
    t_shift = pTrial{1}.t - pTrial{1}.t(idx_releaset);
        % intropolate to pData
    pData(1,:) = interp1(t_shift, pTrial{1}.(cpr_var)(cpr_idx,:), t_grids, 'linear');
    pData_strobe(1,:) = interp1(t_shift, pTrial{1}.Fp(cpr_idx,:), t_grids, 'linear');
    if (ifplot)
        figure(fh_dbg); clf; hold on;
        plot(t_shift, pTrial{1}.(cpr_var)(cpr_idx,:), ...
            'marker', '.', 'color', 'b'); 
        plot(t_grids, pData, ...
            'marker', '.', 'color', 'r');
        xlim(t_range);
        xlabel('time (s)');
    end
    % interpolate the unperturbed trial
        % deal with time and intropolate to npData
    t_shiftcell = cell(1,numnp);
    for ti = 1:numnp
        idx_releaset = find((npTrial{ti}.ts == 5 & diff([1 npTrial{ti}.ts]) == 1)); 
        t_shift = npTrial{ti}.t - npTrial{ti}.t(idx_releaset);
        npData(ti,:) = interp1(t_shift, npTrial{ti}.(cpr_var)(cpr_idx,:), t_grids, 'linear');
        t_shiftcell{ti} = t_shift;
    end
    if (ifplot)
        figure(fh_dbg); clf; hold on;
        for ti = 1:numnp
            plot(t_shiftcell{ti}, npTrial{ti}.(cpr_var)(cpr_idx,:), ...
                'marker', '.', 'color', 'b');
            plot(t_grids, npData(ti,:), ...
                'marker', '.', 'color', 'r');
        end
        xlim(t_range);
        xlabel('time (s)');
    end
    
    % chop the data depend on the perturb time started 
    pData_strobe(diff([0 abs(pData_strobe)])<0) = 0; % remove the later part
    [~, pData_strobe_idx] = min(abs(abs(pData_strobe) - max(abs(pData_strobe))*0.05));
    if (ifplot)
        figure(fh_dbg); clf; hold on;
        plot(t_grids,pData_strobe); 
        plot(t_grids(pData_strobe_idx), pData_strobe(pData_strobe_idx), 'o');
    end
    
    t_grids_0 = t_grids(1:pData_strobe_idx);
    pData0 = pData(:,1:pData_strobe_idx);
    npData0 = npData(:,1:pData_strobe_idx);
    
    if (ifplot)
        figure(fh_dbg); clf; hold on;
        plot(t_grids_0, pData0, 'r');
        plot(t_grids_0, npData0, 'b');
        legend('pert', 'unpert');
        xlabel('time (s)');
    end
    
%   2. Find the nearest one from them, get the index
    idx_nearestTrial = knnsearch(npData0,pData0); 
    if (ifplot || 0)
        try
            figure(fh_dbg); clf; hold on;
        catch 
            figure(); clf; hold on;
        end
        lnh_u = plot(t_grids_0, npData0, 'b');
        lnh_p = plot(t_grids_0, pData0, 'r', 'lineWidth', 2);
        lnh_n = plot(t_grids_0, npData0(idx_nearestTrial,:), 'm', 'lineWidth', 2);

        legend([lnh_p, lnh_n, lnh_u(1)],{'pert', 'pertNearestNeighbor', 'unpert'});
        xlabel('time (s)');
    end

%   3. From the index, get the 'smallest counterpart'
    nearestTrial = otherTrials(idx_nearestTrial);
    
%   4. Plot the perturbed trial overlay with the unperturb ones
    if (ifplot || 0)
        try
            figure(fh_dbg); clf; hold on;
        catch 
            figure(); clf; hold on;
        end
        lnh_u = plot(t_grids, npData, 'b');
        lnh_p = plot(t_grids, pData, 'r', 'lineWidth', 2);
        lnh_n = plot(t_grids, npData(idx_nearestTrial,:), 'm', 'lineWidth', 2);

        legend([lnh_p, lnh_n, lnh_u(1)],{'pert', 'pertNearestNeighbor', 'unpert'});
        xlabel('time (s)');
    end

end

