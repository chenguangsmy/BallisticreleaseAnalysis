function [psudoK] = getDf_dev_dx(trialtmp, ifplot, if3d, ifpulseInMotion)
%GETDF_DEV_DX get the "psudoStiffness" as a sanityCheck, by using dF/dx
%   the trialtmp is the actual trial thata last for enough length. It has
%   the components all needed for calculating the stiffness 

if ~exist('ifplot', 'var')
    ifplot = 0;
else
    fh_sc = figure();
end

if ~exist('if3d', 'var')
    if3d = 0;
end

if ~exist('ifpulseInMotion', 'var')
    ifpulseInMotion = 0; 
end

psudoK = nan;
 
idx_submean_add = -50:-1; % the -100ms before perturbation happens
% find the index for the perturbation 
idx_pert = trialtmp.Fp(1,:)~=0; % take as granted for x- perturb
idx_pert_stt = find(idx_pert); 
idx_pert_stt = idx_pert_stt(1); 
idx_submean = idx_pert_stt + idx_submean_add; 

% subtract the F and x 
if (if3d==0)
    if (~ifpulseInMotion)
        f_0 = mean(trialtmp.f(1,idx_submean));
        f_sub = trialtmp.f(1,:) - f_0;

        switch length(size(trialtmp.ox))
            case 2
                x_0 = mean(trialtmp.ox(1,idx_submean));
                x_sub = trialtmp.ox(1,:) - x_0;
            case 3
                x_0 = mean(trialtmp.ox(1,idx_submean,1));
                x_sub = trialtmp.ox(1,:,1) - x_0;
        end

        t_shift = trialtmp.t - trialtmp.t(idx_pert_stt);

        % find f_peak and the x_peak, possible range is 0:500ms after perturbation
        % initialtion
        peak_range = 1:500;
        [f_peak, idx_fpeak] = max((f_sub(idx_pert_stt + peak_range)));
        idx_fpeak = peak_range(idx_fpeak) + idx_pert_stt; % idx in this trial
        [x_peak, idx_xpeak] = min((x_sub(idx_pert_stt + peak_range)));
        idx_xpeak = peak_range(idx_xpeak) + idx_pert_stt;


        if (ifplot)
            figure(fh_sc);
            axh(1) = subplot(2,1,1); title('force');
            hold on;
            plot(t_shift, f_sub);
            plot(t_shift(idx_fpeak), f_sub(idx_fpeak), ...
                'o', 'color', [1 0 0]);

            axh(2) = subplot(2,1,2); title('displacement');
            hold on;
            plot(t_shift, x_sub);
            plot(t_shift(idx_xpeak), x_sub(idx_xpeak), ...
                'o', 'color', [1 0 0]);

            linkaxes(axh, 'x');
        end
    else 
        f_sub = trialtmp.df(:);
        x_sub = trialtmp.dx(:);

        t_shift = trialtmp.t - trialtmp.t(idx_pert_stt);

        % find f_peak and the x_peak, possible range is 0:500ms after perturbation
        % initialtion
        peak_range = 1:500;
        [f_peak, idx_fpeak] = max((f_sub(idx_pert_stt + peak_range)));
        idx_fpeak = peak_range(idx_fpeak) + idx_pert_stt; % idx in this trial
        [x_peak, idx_xpeak] = min((x_sub(idx_pert_stt + peak_range)));
        idx_xpeak = peak_range(idx_xpeak) + idx_pert_stt;


        if (ifplot)
            figure(fh_sc);
            axh(1) = subplot(2,1,1); title('force');
            hold on;
            plot(t_shift, f_sub);
            plot(t_shift(idx_fpeak), f_sub(idx_fpeak), ...
                'o', 'color', [1 0 0]);

            axh(2) = subplot(2,1,2); title('displacement');
            hold on;
            plot(t_shift, x_sub);
            plot(t_shift(idx_xpeak), x_sub(idx_xpeak), ...
                'o', 'color', [1 0 0]);

            linkaxes(axh, 'x');
        end
    end

xlim([-0.5 0.5]);
psudoK = -f_peak/x_peak;

else % if3d == 1
    
f_0 = mean(trialtmp.f(1:3,idx_submean),2);
f_sub = trialtmp.f(1:3,:) - f_0;



switch length(size(trialtmp.ox))
    case 2
        x_0 = mean(trialtmp.ox(1:3,idx_submean),2);
        x_sub = trialtmp.ox(1:3,:) - x_0;
    case 3
        x_0 = mean(trialtmp.ox(1:3,idx_submean,1),2);
        x_sub = trialtmp.ox(1:3,:,1) - x_0;
end

f_sub_dist = vecnorm(f_sub);
x_sub_dist = vecnorm(x_sub);

t_shift = trialtmp.t - trialtmp.t(idx_pert_stt);

% find f_peak and the x_peak, possible range is 0:500ms after perturbation
% initialtion 
peak_range = 1:150; % guess the distance will not great than this
[f_peak, idx_fpeak] = max(vecnorm(f_sub(:,idx_pert_stt + peak_range)));
idx_fpeak = peak_range(idx_fpeak) + idx_pert_stt; % idx in this trial
[x_peak, idx_xpeak] = max(vecnorm(x_sub(:,idx_pert_stt + peak_range)));
idx_xpeak = peak_range(idx_xpeak) + idx_pert_stt;


if (ifplot)
    figure(fh_sc); 
    axh(1) = subplot(4,2,1); title('force-x'); 
    hold on; 
    plot(t_shift, f_sub(1,:));
    plot(t_shift(idx_fpeak), f_sub(1,idx_fpeak), ...
        'o', 'color', [1 0 0]);

    axh(3) = subplot(4,2,3); title('force-y'); 
    hold on; 
    plot(t_shift, f_sub(2,:));
    plot(t_shift(idx_fpeak), f_sub(2,idx_fpeak), ...
        'o', 'color', [1 0 0]);

    axh(5) = subplot(4,2,5); title('force-z'); 
    hold on; 
    plot(t_shift, f_sub(3,:));
    plot(t_shift(idx_fpeak), f_sub(3,idx_fpeak), ...
        'o', 'color', [1 0 0]);

    axh(7) = subplot(4,2,7); title('force-norm'); 
    hold on; 
    plot(t_shift, f_sub_dist);
    plot(t_shift(idx_fpeak), f_sub_dist(idx_fpeak), ...
        'o', 'color', [1 0 0]);
    
    %%% displacements 
    axh(2) = subplot(4,2,2); title('displacement-x'); 
    hold on;
    plot(t_shift, x_sub(1,:));
    plot(t_shift(idx_xpeak), x_sub(1,idx_xpeak), ...
        'o', 'color', [1 0 0]);

    axh(4) = subplot(4,2,4); title('displacement-y'); 
    hold on;
    plot(t_shift, x_sub(2,:));
    plot(t_shift(idx_xpeak), x_sub(2,idx_xpeak), ...
        'o', 'color', [1 0 0]);

    axh(6) = subplot(4,2,6); title('displacement-z'); 
    hold on;
    plot(t_shift, x_sub(3,:));
    plot(t_shift(idx_xpeak), x_sub(3,idx_xpeak), ...
        'o', 'color', [1 0 0]);

    axh(8) = subplot(4,2,8); title('displacement-norm'); 
    hold on;
    plot(t_shift, x_sub_dist);
    plot(t_shift(idx_xpeak), x_sub_dist(idx_xpeak), ...
        'o', 'color', [1 0 0]);
    ylim([-0.01, 0.04]);

    linkaxes(axh, 'x');
end

xlim([-0.8 0.5]);
psudoK = f_peak/x_peak;

end




end

