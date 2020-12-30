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

%to the right x- : 1,5,9,13,17,21
%6 rows x 4 columns, each column -> 1 movement direction;
%each row: same force threshold, 1->3, force threshold increases;
%bottom 3: closer distance, 4->6, force threshold increases;

% fileName = append('Sonic',sessionNum, 'ForceData.mat');
% load(fileName);

figure('Name', 'Trial Average Force for All Conditions');
totalConditionTypes = 24;
rows = 6;
cols = totalConditionTypes / rows;

for i = 1:totalConditionTypes
    subplot(rows,cols,i);
    hold on;
    for j = 1:rows
        if (j == 1 || j ==2 )
            plot(AllTime(i,:),AllForceTrialAverage((i-1)*forceRowsAverage + j + 6, :));
        else
            plot(AllTime(i,:),AllForceTrialAverage((i-1)*forceRowsAverage + j, :));
        end
        xlim([-2.5, 1]);
        ylim([-15, 10]);
    end
    xline(0, '-.');
    hold off;
end
legend('Fforward','Fup', 'Fz','TorqueX','TorqueY','TorqueZ');

delete(findall(gcf,'type','annotation'));
%top row title
annotation('textbox', [.17,.93,.08,.05],'String','Right (-z)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.38,.93,.08,.05],'String','Left (+z)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.57,.93,.15,.05],'String','Forward', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.77,.93,.15,.05],'String','Backward', ...
    'EdgeColor', 'none' ,'FontSize', 20);

%column wise title
annotation('textbox', [.02,.85,.095,.06],'String','Low Fthreshold (2), Far Distance (0.06)', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.02,.7,.095,.06],'String','Mid Fthreshold (5), Far Distance (0.06)', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.02,.55,.099,.06],'String','High Fthreshold (8), Far Distance (0.06)', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.015,.41,.12,.06],'String','Low Fthreshold (2), Close Distance (0.04)', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.015,.28,.105,.06],'String','Mid Fthreshold (5), Close Distance (0.04)', ...
    'EdgeColor', 'none' ,'FontSize', 12);
annotation('textbox', [.015,.13,.12,.06],'String','High Fthreshold (8), Close Distance (0.04)', ...
    'EdgeColor', 'none' ,'FontSize', 12);

%TODO: some observations: pos not symmetrical around 0, x is around -0.12, and y is around 0.5
titleName = append('Session ', sessionNum,' Trial Average Force');
annotation('textbox', [.43, 0.01,.19,.04],'String',titleName, ...
    'EdgeColor', 'none' ,'FontSize', 15, 'FontWeight', 'bold');