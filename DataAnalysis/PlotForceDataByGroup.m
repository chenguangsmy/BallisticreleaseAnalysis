%combo number in order: right (x-), left, top(y+ from 0.5), down(- from 0.5),
%x+ point up, - down, y forward/backward, z left/right
%not symmetrical around 0, x is around -0.12, and y is around 0.5
%longer distance -> shorter
%force levels: low, mid high
%array of size 5: 1,3,5,7,9

%Data Structure:
%TrialAverage: 6 rows per condition;
%Time 1 row per condition; 
%AlignedByBlock: total reps (blocks * rep per block) * 6 rows per condition
%In row order: Fx for all reps of a given condition, y, z, torque x, y, z; 
%then repeat for next condition
%RawBeforeAlignment: not aligned at time 0, in order x,y,z,torque x,y,z for
%1 rep of a given condition; then next rep...; after all reps for this
%condition; starts with rep1 of next condition.
%DataIndex: the index used to locate data ror each successful rep
%1 row corresponds to 1 rep, total rows = rep per condition * conditions.

%Group by movement direction: 4 direction: 3x2 traces per plot; 4 plots * 6
%force measurements
%Group by force threshold: 2 distance per plot; 4 directions x 3 plots
%plots
%Group by travel distance: 3 force per plot, 1 direction force per plot; 
% 4 directions*2 distance type * 6 force plots
%Functionize everything

%6 rows x 4 columns, each row: 1 force direction; each plot: this direction
%for all force threshold and movement distance; 
%each column: 1 movement direction
%trial average in order: con1 x,y,z,tx,ty,tz

figure('Name', 'Trial Average Force Grouped By Force Threshold');
tracePerPlot = 6;
cols = totalConditionTypes / tracePerPlot;
indexOrder = zeros(1, 24);
for i = 1 : 4
    for j = 1: 6
      indexOrder(1, (i-1)*6+j) = i + (j-1)*4;
    end
end

style = ['r','y','b','r','y','b'];

for i = 1:4 %col
    for k = 1 :6 %row
        %going horizontal first, in order 1,5,
        subplot(tracePerPlot,cols,indexOrder((k-1)*4+i));
        hold on;
        %to the right x- : 1,5,9,13,17,21, in rows: 1, 5*6+1; 2,6,10,14, in
        %rows: 2*6+1,
        for j = 1:tracePerPlot %one plot, Fx for 3 threshold and 2 distance, con1,5,
            if (j >3)
                plot(AllTime(i,:),AllTrialAverage((indexOrder(j+(i-1)*6)-1)*6 + k, :), style(j),'LineStyle','--');
            else
                plot(AllTime(i,:),AllTrialAverage((indexOrder(j+(i-1)*6)-1)*6 + k, :), style(j));
            end
    %         xlim([-2.5, 1]);
    %         ylim([-15, 10]);
        end
        xline(0, '-.');
        hold off;
    end
end
legend('Low(2)-Far(0.06)', 'Mid(5)-Far','High(8)-Far','Low-Close(0.04)','Mid-Close(0.04)','High-Close(0.04)');

delete(findall(gcf,'type','annotation'));
%top row title
annotation('textbox', [.17,.93,.08,.05],'String','Right (-z)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.38,.93,.08,.05],'String','Left (+z)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.57,.93,.15,.05],'String','Forward (+y)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.77,.93,.15,.05],'String','Backward (-y)', ...
    'EdgeColor', 'none' ,'FontSize', 20);

%column wise title
annotation('textbox', [.02,.85,.095,.06],'String','Fx', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.02,.7,.095,.06],'String','Fy', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.02,.55,.099,.06],'String','Fz', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.015,.41,.12,.06],'String','TorqueX', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.015,.28,.105,.06],'String','TorqueY', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.015,.13,.12,.06],'String','TorqueZ', ...
    'EdgeColor', 'none' ,'FontSize', 12);