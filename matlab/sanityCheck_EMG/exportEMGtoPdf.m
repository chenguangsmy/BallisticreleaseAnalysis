% export EMG data to a pdf file for easy check with 

blp = ballisticReleaseTaksPlots();
for i = 1:4
    fh_tmp(i) = figure('Visible', 'on');
end
for subj_i = 1:4
    for dir_i = [1 3 2 4]
        for fce_i = 1:3
            for dist_i = 1:3
                display(['dir' num2str(dir_i) 'fce' num2str(fce_i) 'dist' num2str(dist_i)]);
                clf(fh_tmp(1)); 
                clf(fh_tmp(2));
                clf(fh_tmp(3));
                clf(fh_tmp(4));
                fh_tmp = blp.plotEMG_release_specificCond_fce(subj_i,dir_i,fce_i,dist_i, fh_tmp);
                for fig_i = 1:4
%                     set(fh_tmp(fig_i), 'unit', 'inch', 'position', [0 0 8 11]);

                if (subj_i == 1 && dir_i == 1 && fce_i == 1 && dist_i == 1 && fig_i == 1)
%                     export_fig(fh_tmp(fig_i), '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_difftr.pdf', '-preserve_size', '-nocrop');
%                     exportgraphics(fh_tmp(fig_i), '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_difftr.pdf');
                    exportgraphics(fh_tmp(fig_i), '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_difftr_color_btw05.pdf');
                else
%                     export_fig(fh_tmp(fig_i), '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_difftr.pdf', '-preserve_size', '-nocrop', '-append');
%                     exportgraphics(fh_tmp(fig_i), '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_difftr.pdf', 'Append', true);
                    exportgraphics(fh_tmp(fig_i), '/Users/cleave/Documents/projPitt/BallisticreleaseAnalysis/matlab/dataCommuPlots/techrept20220825/subj_EMG_difftr_color_btw05.pdf', 'Append', true);
                end

                end
            end
        end 
    end
end