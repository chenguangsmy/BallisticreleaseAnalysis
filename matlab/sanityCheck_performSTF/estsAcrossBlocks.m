%
close all;
clc;
% dir 1: [2 3 4 8 9 12 13 15 17 20]
% dir 2: [2 3 4 8 9 12]
% dir 3: [2 3 4 8 9 12 13 15 17 20]
% dir 4: [2 3 4 8 9 12 15 17 20]
dir_i = 2;
for subj_i = [ 2 3 4 8 9 12 13 14 15 16 17 18 20 21] % subj5, 6, overtrained. 7, not following
    
    fdir = 'test_data'
    fname = sprintf('exampleData_subj%ddir%d.mat', subj_i, dir_i)


    load([fdir, '/' fname], 'headers');


    x_tars = [headers.trialHeader.tarL];
    F_tars = [headers.trialHeader.tarF];
    k_tars = [F_tars./x_tars];
    %
    k_ests = [headers.trialHeader.estK];
    b_ests = [headers.trialHeader.estB];
    m_ests = [headers.trialHeader.estB];
    fit_ests = [headers.trialHeader.estFIT];

    outcome= logical([headers.trialHeader.outcome]);

    dex_qualified = fit_ests > 0 ;

    x_tars = x_tars(dex_qualified);
    F_tars = F_tars(dex_qualified);
    k_tars = k_tars(dex_qualified);
    outcome = outcome(dex_qualified);

    k_ests = k_ests(dex_qualified);
    b_ests = b_ests(dex_qualified);
    m_ests = m_ests(dex_qualified);
    fit_ests = fit_ests(dex_qualified);
    % figure 1 plot k
    figure(); hold on;
    trials_all = 1:length(x_tars); hold on;
    plot(trials_all, k_tars);
    plot(trials_all, k_ests);
    plot(trials_all(outcome), k_ests(outcome), 'g*');
    plot(trials_all(~outcome), k_ests(~outcome), 'r*');
    title(['subject' num2str(subj_i) ' dir' num2str(dir_i)]);
end

%% 
% do a block-trending EMG est 
close all;
clc;

for dir_i = 1:4
    fh(dir_i) = figure('unit', 'inch', 'position', [0 0 6 4]); hold on;
    set(gca, 'fontsize', 15, 'linewidth', 1);
    for b_i = 1:10
        k_ests_poolsubj = [];
        tNo_poolssubj = [];    % the actual trial num in this block
        tNorsp_poolsubj = [];  % resample this to a value [0-9]; 
        outcome_poolsubj = [];

        x_tar = [];
        F_tar = [];
        x_tar_tolerance = [-0.01 0.01];
        F_tar_tolerance = [-2 2];

        for subj_i = [ 2 3 4 8 9 12 13 14 15 16 17 18 20 21] % subj5, 6, overtrained. 7, not following

            fdir = 'test_data';
            fname = sprintf('exampleData_subj%ddir%d.mat', subj_i, dir_i);
            load([fdir, '/' fname], 'headers');

            bNo    = [headers.trialHeader.bNo]; 

            x_tars = [headers.trialHeader.tarL];
            F_tars = [headers.trialHeader.tarF];
            k_tars = [F_tars./x_tars];
            %
            k_ests = [headers.trialHeader.estK];
            b_ests = [headers.trialHeader.estB];
            m_ests = [headers.trialHeader.estB];
            fit_ests = [headers.trialHeader.estFIT];

            outcome= logical([headers.trialHeader.outcome]);

            dex_qualified = fit_ests > 0 & k_ests > 0;
            b_sel  = bNo == b_i;

            x_tars = x_tars(dex_qualified & b_sel);
            F_tars = F_tars(dex_qualified & b_sel);
            k_tars = k_tars(dex_qualified & b_sel);
            outcome = outcome(dex_qualified & b_sel);

            k_ests = k_ests(dex_qualified & b_sel);
            b_ests = b_ests(dex_qualified & b_sel);
            m_ests = m_ests(dex_qualified & b_sel);
            fit_ests = fit_ests(dex_qualified & b_sel);


            k_ests_subj = [k_ests];
            tNo_subj = 1:length(k_ests);    % the actual trial num in this block
            tNorsp_subj = [(tNo_subj-1)/(tNo_subj(end)-1)*9 + 9*(b_i-1)];  % resample this to a value [0-9]; 


            % no need to do this every time 
            x_tar = unique(x_tars); 
            F_tar = unique(F_tars); 
            x_tar_lo = x_tar + x_tar_tolerance(1); 
            x_tar_hi = x_tar + x_tar_tolerance(2); 
            F_tar_lo = F_tar + F_tar_tolerance(1); 
            F_tar_hi = F_tar + F_tar_tolerance(2);

            k_tar = F_tar/x_tar; 
            k_tar_lo = F_tar_lo / x_tar_hi; 
            k_tar_hi = F_tar_hi / x_tar_lo; 

            k_ests_poolsubj = [k_ests_poolsubj k_ests_subj];
            tNo_poolssubj   = [tNo_poolssubj tNo_subj];
            tNorsp_poolsubj = [tNorsp_poolsubj tNorsp_subj];
            outcome_poolsubj = [outcome_poolsubj outcome];

        end
        
            % fixing plots..
            plot(tNorsp_poolsubj(logical(outcome_poolsubj)), k_ests_poolsubj(logical(outcome_poolsubj)), 'g.');
            plot(tNorsp_poolsubj(logical(~outcome_poolsubj)), k_ests_poolsubj(logical(~outcome_poolsubj)), 'r.');
            xs = [min(tNorsp_poolsubj), max(tNorsp_poolsubj)]; 
            ys_med = [k_tar, k_tar];
            ys_min = [k_tar_lo, k_tar_lo];
            ys_max = [k_tar_hi, k_tar_hi];
            line(xs, ys_med, 'linewidth', 2);
            line(xs, ys_min, 'linewidth', 1);
            line(xs, ys_max, 'linewidth', 1);
            title([ ' dir' num2str(dir_i)]);

    end
    ylim([0 1800])
    xlabel('calibrated trial count');
    ylabel('K_{est} (N/m)');
    saveas(fh(dir_i), ['beh_figures/stiffnessEstAcrossTrials.dir' num2str(dir_i) '.png'])
end


%% do a plot compare failed and succeed trial stiffness 
% do a block-trending emg est 
close all;
clc;

for dir_i = 1:4
%     fh(dir_i) = figure('unit', 'inch', 'position', [0 0 6 4]); hold on;
    set(gca, 'fontsize', 15, 'linewidth', 1);
    for b_i = 1:10
        figure();
        k_ests_poolsubj = [];
        tNo_poolssubj = [];    % the actual trial num in this block
        tNorsp_poolsubj = [];  % resample this to a value [0-9]; 
        outcome_poolsubj = [];

        x_tar = [];
        F_tar = [];
        x_tar_tolerance = [-0.01 0.01];
        F_tar_tolerance = [-2 2];

        for subj_i = [ 2 3 4 8 9 12 13 14 15 16 17 18 20 21] % subj5, 6, overtrained. 7, not following

            fdir = 'test_data';
            fname = sprintf('exampleData_subj%ddir%d.mat', subj_i, dir_i);
            load([fdir, '/' fname], 'headers');

            bNo    = [headers.trialHeader.bNo]; 

            x_tars = [headers.trialHeader.tarL];
            F_tars = [headers.trialHeader.tarF];
            k_tars = [F_tars./x_tars];
            %
            k_ests = [headers.trialHeader.estK];
            b_ests = [headers.trialHeader.estB];
            m_ests = [headers.trialHeader.estB];
            fit_ests = [headers.trialHeader.estFIT];

            outcome= logical([headers.trialHeader.outcome]);

            dex_qualified = fit_ests > 0 & k_ests > 0;
            b_sel  = bNo == b_i;

            x_tars = x_tars(dex_qualified & b_sel);
            F_tars = F_tars(dex_qualified & b_sel);
            k_tars = k_tars(dex_qualified & b_sel);
            outcome = outcome(dex_qualified & b_sel);

            k_ests = k_ests(dex_qualified & b_sel);
            b_ests = b_ests(dex_qualified & b_sel);
            m_ests = m_ests(dex_qualified & b_sel);
            fit_ests = fit_ests(dex_qualified & b_sel);


            k_ests_subj = [k_ests];
            tNo_subj = 1:length(k_ests);    % the actual trial num in this block
            tNorsp_subj = [(tNo_subj-1)/(tNo_subj(end)-1)*9 + 9*(b_i-1)];  % resample this to a value [0-9]; 


            % no need to do this every time 
            x_tar = unique(x_tars); 
            F_tar = unique(F_tars); 
            x_tar_lo = x_tar + x_tar_tolerance(1); 
            x_tar_hi = x_tar + x_tar_tolerance(2); 
            F_tar_lo = F_tar + F_tar_tolerance(1); 
            F_tar_hi = F_tar + F_tar_tolerance(2);

            k_tar = F_tar/x_tar; 
            k_tar_lo = F_tar_lo / x_tar_hi; 
            k_tar_hi = F_tar_hi / x_tar_lo; 

            k_ests_poolsubj = [k_ests_poolsubj k_ests_subj];
            tNo_poolssubj   = [tNo_poolssubj tNo_subj];
            tNorsp_poolsubj = [tNorsp_poolsubj tNorsp_subj];
            outcome_poolsubj = [outcome_poolsubj outcome];

        end
        
%             % fixing plots..
%             plot(tNorsp_poolsubj(logical(outcome_poolsubj)), k_ests_poolsubj(logical(outcome_poolsubj)), 'g.');
%             plot(tNorsp_poolsubj(logical(~outcome_poolsubj)), k_ests_poolsubj(logical(~outcome_poolsubj)), 'r.');
%             xs = [min(tNorsp_poolsubj), max(tNorsp_poolsubj)]; 
%             ys_med = [k_tar, k_tar];
%             ys_min = [k_tar_lo, k_tar_lo];
%             ys_max = [k_tar_hi, k_tar_hi];
%             line(xs, ys_med, 'linewidth', 2);
%             line(xs, ys_min, 'linewidth', 1);
%             line(xs, ys_max, 'linewidth', 1);
%             title([ ' dir' num2str(dir_i)]);

        histogram(k_ests_poolsubj(logical(outcome_poolsubj))); hold on;
        histogram(k_ests_poolsubj(logical(~outcome_poolsubj)));
    end
%     ylim([0 1800])
%     xlabel('calibrated trial count');
%     ylabel('K_{est} (N/m)');
%     saveas(fh(dir_i), ['beh_figures/stiffnessEstAcrossTrials.dir' num2str(dir_i) '.png'])
end


