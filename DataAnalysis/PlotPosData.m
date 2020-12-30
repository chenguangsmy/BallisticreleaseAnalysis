%Author: Shuqi Liu 
%Date: 2019-12-27 17:06 
%File Description: Plot the position data. Plot the trajectory by
%condition, each row = 1 force and target distance configuration, each
%column = 1 reach direction in order right (-x), left (+x), forward (+y)
%and backward(-y). For position, the z direction is up/down

figure('Name', 'Position Trajectory for All Conditions');
totalConditionTypes = 24;
posRows = 3;
%add the center/home position offset
xoffset = 0.12;
yoffset = -0.5;
zoffset = -0.25;

style=["b","r","g","b--","r--","g--"];

i = 1;
plot3(AllPosTrialAverage((i-1)*posRows + 1, :) +xoffset, ...
    AllPosTrialAverage((i-1)*posRows + 2, :) + yoffset,...
    AllPosTrialAverage((i-1)*posRows + 3, :) + zoffset, style(fix((i-1)/4)+1));

hold on;
%plot data in condition1, 2, 3,... order
for i = 2:totalConditionTypes
    fprintf('Plotting Condition: %d\n',i);

    %TODO?plot 3D doesn't really work here...
%     plot3(AllPosAlignedByBlock(1, :) +xoffset, ...
%         AllPosAlignedByBlock(16, :) + yoffset,...
%         AllPosAlignedByBlock(31, :) + zoffset);
    plot3(AllPosTrialAverage((i-1)*posRows + 1, :) +xoffset, ...
        AllPosTrialAverage((i-1)*posRows + 2, :) + yoffset,...
        AllPosTrialAverage((i-1)*posRows + 3, :) + zoffset, style(fix((i-1)/4)+1));
    %start = (i-1)*blocks*posRows;
%     plot3(AllPosAlignedByBlock(start+1:start+blocks, :) + xoffset,...
%          AllPosAlignedByBlock(start + blocks + 1 : start + 2*blocks, :) + yoffset,...
%          AllPosAlignedByBlock(start+2*blocks+1 : start + 3*blocks, :) + zoffset);

end

xlim([-0.1, 0.1]);
ylim([-0.1, 0.1]);
zlim([-0.1, 0.1])
xlabel('x');
ylabel('y');
zlabel('z');

%adding reference line starting at origin (0,0,0), seems tricky to
%align axis naturally at origin
reference = linspace(-0.2,0.2,100);
plot3(reference,zeros(100), zeros(100),'k');
plot3(zeros(100),reference, zeros(100),'k');
plot3(zeros(100), zeros(100), reference, 'k');

hold off;
%the maxonsetIndex would be the time 0 point.
%xline(0, '-.');
%TODO: figure out the starting and ending target area in real
%coordinates

% delete(findall(gcf,'type','annotation'));
% %top row title
% annotation('textbox', [.17,.93,.08,.05],'String','Right (-z)', ...
%     'EdgeColor', 'none' ,'FontSize', 20);
% annotation('textbox', [.38,.93,.08,.05],'String','Left (+z)', ...
%     'EdgeColor', 'none' ,'FontSize', 20);
% annotation('textbox', [.57,.93,.15,.05],'String','Forward (+y)', ...
%     'EdgeColor', 'none' ,'FontSize', 20);
% annotation('textbox', [.77,.93,.15,.05],'String','Backward (-y)', ...
%     'EdgeColor', 'none' ,'FontSize', 20);
% 
% %column wise title
% annotation('textbox', [.02,.85,.095,.06],'String','Low Fthreshold (2), Far Distance (0.06)', ...
%     'EdgeColor', 'none' ,'FontSize', 12);
% annotation('textbox', [.02,.7,.095,.06],'String','Mid Fthreshold (5), Far Distance (0.06)', ...
%     'EdgeColor', 'none' ,'FontSize', 12);
% annotation('textbox', [.02,.55,.099,.06],'String','High Fthreshold (8), Far Distance (0.06)', ...
%     'EdgeColor', 'none' ,'FontSize', 12);
% annotation('textbox', [.015,.41,.12,.06],'String','Low Fthreshold (2), Close Distance (0.04)', ...
%     'EdgeColor', 'none' ,'FontSize', 12);
% annotation('textbox', [.015,.28,.105,.06],'String','Mid Fthreshold (5), Close Distance (0.04)', ...
%     'EdgeColor', 'none' ,'FontSize', 12);
% annotation('textbox', [.015,.13,.12,.06],'String','High Fthreshold (8), Close Distance (0.04)', ...
%     'EdgeColor', 'none' ,'FontSize', 12);
% 
% %TODO: some observations: pos not symmetrical around 0, x is around -0.12, and y is around 0.5
% titleName = append('Session ', sessionNum,' Position Trajectory');
% annotation('textbox', [.43, 0.01,.19,.04],'String',titleName, ...
%     'EdgeColor', 'none' ,'FontSize', 15, 'FontWeight', 'bold');