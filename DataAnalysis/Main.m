%Author: Shuqi Liu
%Date: 2019-12-27 17:06
%File Description: The main file to run when processing data. User needs to
%update the sessionNum variable to the target session number. The data and
%the current file needs to follow the folder structure '../subjectNameData/', and
%the data file name needs to be subjectName.sessionNum.mat. The processed data is
%then saved to
%'../subjectNameDataProcessed/subjectName.sessionNum.ForceDataAligned.mat'

close all;
clear all;

%TODO: change session number each time
sessionNum = '01317';
subjectName = 'KingKong'; %enter subjec name here
fileName = append('../',subjectName,'Data/',subjectName,'.',sessionNum, '.mat');
fileName = append('../','rawdata/',subjectName,'.',sessionNum, '.mat');

%load data
dataCurr = load(fileName);

%TODO: make the variables clear using global keyword and make things function. 
% global forceRows totalConditionTypes conditionIndex movementOnsetState ...
%     repPerBlock maxCols numberBlocks blocks initPlot plotPos plotForce; 

%set up global variable about if plot will be shown
initPlot = 0;
plotPos = 0;
plotForce = 0;

forceRows= 6;
forceRowsAverage = 8; %include 2 rows each condition for mapped coordinates
posRows = 3;
totalConditionTypes = 24;
maxCols = 1000;
repPerBlock = 5;
movementOnsetState = 4;
%TODO: manually set based on the expeirment, only use complete blocks.
numberBlocks = dataCurr.Data.BlockNo(length(dataCurr.Data.BlockNo)) - 1;

%total number of repititions of one condition. block number * repPerBlock
blocks = numberBlocks * repPerBlock;
%FIXME: for some reason can't use block number read from the array data for the subplot command?

%Current experiment set up: right (-z), left(+z), top(y+), down(-y), 5 rep
%each, left-right x 5, then top-down x5.
%longer distance -> shorter
%force levels: low, mid high
%Hence the trial index of each condition is [1,3,5,7,9; 2,4,6,8,10; 
% 11,13,15,17,19; 12,14,16,18,20, ...]

%set up the condition index array where each row contains the 5 index of
%trials for each condition. 
allConditionIndex = zeros(totalConditionTypes,repPerBlock);
for i = 1:12
    for j = 1:2
        allConditionIndex((i-1)*2+j,:) = (i-1)*10 + j : 2 : i*10;
    end
end

%initialize output variables.
AllForceTrialAverage = zeros(forceRowsAverage * totalConditionTypes, maxCols);
AllForceAlignedByBlock = zeros (forceRows * totalConditionTypes * blocks, maxCols);
AllForceRaw = zeros (forceRows * totalConditionTypes * blocks, maxCols);
AllPosTrialAverage = zeros(posRows * totalConditionTypes, maxCols);
AllPosAlignedByBlock = zeros (posRows * totalConditionTypes * blocks, maxCols);
AllPosRaw = zeros (posRows * totalConditionTypes * blocks, maxCols);
AllVelTrialAverage = zeros(posRows * totalConditionTypes, maxCols);
AllVelAlignedByBlock = zeros (posRows * totalConditionTypes * blocks, maxCols);
AllVelRaw = zeros (posRows * totalConditionTypes * blocks, maxCols);

%initialize time and data index with NaN since real data couldn't take NaN
AllDataIndex = NaN (totalConditionTypes * blocks, maxCols);
AllTime = NaN(totalConditionTypes, maxCols);

%TODO: the max onset time and max movement time if found should also be
%saved in the output.

FindMaxOnsetIndex

for conditionType = 1:totalConditionTypes
    fprintf('Processing Condition: %d\n',conditionType);
    %an array of the trial index for the current condition
    conditionIndex = allConditionIndex(conditionType,:);
    %validate if block number is different than 5, needs to change the code
    FindTrialIndexesForCondition
%     PositionData
    AlignAndPopulateData 
end

outputFileName = append('../',subjectName,'Data/Processed/',subjectName,'.', sessionNum,'.DataAligned.mat');
save(outputFileName, 'AllForceTrialAverage', 'AllForceAlignedByBlock', 'AllForceRaw',... 
    'AllPosTrialAverage', 'AllPosAlignedByBlock', 'AllPosRaw',...
    'AllDataIndex','AllTime', 'numberBlocks', 'maxOnsetIndex', 'blocks');

load(outputFileName);
PlotForceData;
PlotPosData;
PlotVelocityData;
