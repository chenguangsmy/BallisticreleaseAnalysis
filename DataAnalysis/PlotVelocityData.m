%Author: Shuqi Liu 
%Date: 2020-01-03 11:17 
%File Description: Plot the velocity data. Plot the velocity reading
%by condition, each row = 1 force and target distance configuration, each
%column = 1 reach direction in order right (-x), left (+x), forward (+y)
%and backward(-y). For position, the z direction is up/down

figure('Name', 'Trial Average Velocity for All Conditions');
totalConditionTypes = 24;
rows = 6;
cols = totalConditionTypes / rows;

for i = 1:totalConditionTypes
    subplot(rows,cols,i);
    plot(AllTime(i,:),AllVelTrialAverage((i-1)*posRows + 1:i*posRows, :));
    %xlim([-2.5, 1]);
    %ylim([-15, 10]);
    xline(0, '-.');
end
legend('Vel-X','Vel-Y','Vel-Z');

delete(findall(gcf,'type','annotation'));
%top row title
annotation('textbox', [.17,.93,.08,.05],'String','Right (-x)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.38,.93,.08,.05],'String','Left (+x)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.57,.93,.15,.05],'String','Forward(+y)', ...
    'EdgeColor', 'none' ,'FontSize', 20);
annotation('textbox', [.77,.93,.15,.05],'String','Backward(-y)', ...
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
titleName = append('Session ', sessionNum,' Trial Average Velocity');
annotation('textbox', [.43, 0.01,.19,.04],'String',titleName, ...
    'EdgeColor', 'none' ,'FontSize', 15, 'FontWeight', 'bold');