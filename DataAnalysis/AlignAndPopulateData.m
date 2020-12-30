%Author: Shuqi Liu 
%Date: 2019-12-27 17:06 
%File Description: Pre-process the force data to be aligned at movement
%onset time. Save the processsed data to the output variable.
%Output variable format:
%- AllTime 1 row = the time axis per condition; 
%- AllDataIndex: the index used to locate data ror each successful rep
%1 row corresponds to 1 rep, total rows = rep per condition * conditions.
%condition1 first, then condition 2, etc. for now, each condition contains
%the same number of rows.
%- AllForceAlignedByBlock: total reps (blocks * rep per block) * 6 rows per condition
%In row order: Fx for all reps of a given condition, y, z, torque x, y, z; 
%then repeat for next condition
%- AllForceRaw: not aligned at time 0, in order x,y,z,torque x,y,z for
%1 rep of a given condition; then next rep...; after all reps for this
%condition, go to next condition.
%- AllForceTrialAverage: 6 rows per condition, all force for condition1, then
%next condition, etc.
%- AllPosAlignedByBlock: position for each condition, where x for all reps
%in one condition is listed first, then y,z,; then next condition
%- AllPosRaw: Position for each condition not aligned at time 0.
%- AllPosTrialAverage: position average for condition 1 in x,y,z and then
%next condition.

%show plots as the data is processed
if (plotForce == 1)
    figure('Name','Force'); 
end

%force in 6 directions for each repetition. 
forceByBlock = zeros(forceRows * blocks, maxCols);
%position in 3 direction for each repetition
posByBlock = zeros(posRows * blocks, maxCols);
%velocity in 3 direction for each repetition
velByBlock = zeros(posRows * blocks, maxCols);

%The relative movement onset index (from the beginning of a trial)
onsetIndexByBlock = zeros(1,blocks);
maxMoveTime = 0; %TODO: figure it out in the beginning
%1 block = 1 repetition or trial of a given condition.
for i = 1 : blocks
    IndexPerTrial = trialIndexForCurrentCondition(i, :);
    IndexPerTrial = nonzeros(IndexPerTrial);
    
    if (movementOnsetIndex(1,i) == 0)
        %Current condition didn't run any many repetitions. Should never
        %come here for now since the code setup the blocks to be calculated
        %based on fully completed blocks.
        fprintf('No data for combo no: %d, rep: %d\n\n',conditionIndex, i);
        continue
    end
    
    time0Index = find(IndexPerTrial == movementOnsetIndex(1,i));
    onsetIndexByBlock (1, i) = time0Index;
    moveTime = length(IndexPerTrial) - time0Index;
    %fprintf('Move Time: %d\n',moveTime);
    %TODO: figure out how to calculate max movetime for all in the
    %beginning
    if (moveTime > maxMoveTime)
       maxMoveTime = moveTime;
    end

    force = dataCurr.Data.Force.Sensor(:, IndexPerTrial);
    forceSize = size(force);
    forceByBlock((i-1)*forceRows+1 : i*forceRows, 1 : forceSize(1, 2)) = force;
    
    if (plotForce == 1)
        %construct time point for x axis.
        time = 0 : 0.02 : (length(IndexPerTrial) -1) * 0.02;
        %offset the time array to be 0 at movement onset time.
        time = time - (time0Index-1) * 0.02;
        plotIndex = (i-1)*forceRows+1;
        for k = 1:forceRows
            subplot(blocks,forceRows, plotIndex);
            plot(time, force(k,:));
            plotIndex = plotIndex + 1;
        end
    end
       
    %Original pos data in column x,y,z, each row = 1 time stamp; translates
    %to each row = 1 direction
    pos = dataCurr.Data.Position.Actual(IndexPerTrial, :);
    posSize = size(pos);
    posByBlock((i-1)*posRows+1 : i*posRows, 1 : posSize(1, 1)) = pos';
    
    vel = dataCurr.Data.Velocity.Actual(:,IndexPerTrial);
    velSize = size(vel);
    velByBlock((i-1)*posRows+1 : i*posRows, 1 : velSize(1, 2)) = vel;
    
    %TODO: implement the plotPos; probably optimize this to avoid
    %repetition as well. 
    
end

if (plotForce == 1)
    figure('Name','Force Aligned at time 0'); 
end

%The max number of columns needed for the current condition. 
alignedCols = maxOnsetIndex+maxMoveTime;

%rearrange to be force x for all repetitions first, then y, z, torque x,y,z
alignedForceByBlock = NaN(forceRows*blocks, alignedCols);
%rearrange to be position x for all repetitions first, then y, z
alignedPosByBlock = NaN(posRows*blocks, alignedCols);
%rearrange to be position x for all repetitions first, then y, z
alignedVelByBlock = NaN(posRows*blocks, alignedCols);

%construct time for the x-axis, where t = 0 is movement onset
time = 0 : 0.02 : (alignedCols-1) * 0.02;
time = time - (maxOnsetIndex-1) * 0.02;

for i = 1 : blocks
    time0Index = onsetIndexByBlock(1, i);
    %time stamp padding so that all data will have time 0 at the same
    %column
    offset = maxOnsetIndex - time0Index +1;
    for j = 1 : forceRows
        toFill = forceByBlock((i-1)*forceRows + j,:);
        nonzeroCols = length(nonzeros(toFill));
        alignedForceByBlock(i+(j-1) * blocks, offset : offset+nonzeroCols - 1) = ...
            toFill(1, 1:nonzeroCols);
    end
    
    for j = 1 : posRows
        toFill = posByBlock((i-1)*posRows + j,:);
        nonzeroCols = length(nonzeros(toFill));
        alignedPosByBlock(i+(j-1) * blocks, offset : offset+nonzeroCols - 1) = ...
            toFill(1, 1:nonzeroCols);
    end
    %TODO: implment plotting pos code
    
    for j = 1 : posRows
        toFill = velByBlock((i-1)*posRows + j,:);
        nonzeroCols = length(nonzeros(toFill));
        alignedVelByBlock(i+(j-1) * blocks, offset : offset+nonzeroCols - 1) = ...
            toFill(1, 1:nonzeroCols);
    end
    
    %if plot, plot the force, each row = 1 trial, each column = 1 direction
    if plotForce == 1 
        plotIndex = (i-1)*forceRows+1;
        for k = 1 : forceRows
            %subplot: rows, columns, index from left to right and top to
            %bottom
            subplot(blocks,forceRows, plotIndex);
            plot(time, alignedForceByBlock(i+(k-1)*blocks, :));
            xlim([time(1), time(length(time))]);
            plotIndex = plotIndex + 1;
        end
    end
end

%trial average, where each row = 1 force direction, has 2 more rows than
%other force data, include 2 mapped dimensions
trialAverageForceData = zeros(forceRowsAverage, alignedCols);
trialAveragePosData = zeros(posRows, alignedCols);
trialAverageVelData = zeros(posRows, alignedCols);
if (plotForce == 1)
    figure('Name' ,'Trial Average Force');
end

%populate regular dimensions first
for i = 1 : forceRows
    trialAverageForceData(i, :) = mean(alignedForceByBlock((i-1)*blocks+1:i*blocks, :));
    if (plotForce == 1)
        subplot(forceRows,1,i);
        plot(time,trialAverageForceData(i, :));
    end
end
%add 2 rows of mapped force to horizontal (forward/backward) and vertical
%theta = angle between +x and horizontal plane counterclockwise
theta = pi/4;
%horizontal force = -cosTheta * x + sinTheta * y; assuming forward is
%positive
trialAverageForceData(forceRows+1, :) = - trialAverageForceData(1,:) * cos(theta)...
    + trialAverageForceData(2,:) * sin(theta);
%vertical force = sinTheta * x + cosTheta * y, assuming upward is positive
trialAverageForceData(forceRows+2, :) = trialAverageForceData(1,:) * sin(theta)...
    + trialAverageForceData(2,:) * cos(theta);


for i = 1 : posRows
    trialAveragePosData(i, :) = mean(alignedPosByBlock((i-1)*blocks+1:i*blocks, :));
    trialAverageVelData(i, :) = mean(alignedVelByBlock((i-1)*blocks+1:i*blocks, :));
   
    %TODO: implement plotting position
%     if (plotPos == 1)
%         subplot(forceRows,1,i);
%         plot(time,trialAverageForceData(i, :));
%     end
end

%Assume a variable for all data exists: 
%- AllForceTrialAverage: 6 rows per condition, all force for condition1, then
%next condition, etc.
%- AllTime 1 row = the time axis per condition; 
%- AllForceAlignedByBlock: total reps (blocks * rep per block) * 6 rows per condition
%In row order: Fx for all reps of a given condition, y, z, torque x, y, z; 
%then repeat for next condition
%- AllForceRaw: not aligned at time 0, in order x,y,z,torque x,y,z for
%1 rep of a given condition; then next rep...; after all reps for this
%condition, go to next condition.
%- AllDataIndex: the index used to locate data ror each successful rep
%1 row corresponds to 1 rep, total rows = rep per condition * conditions.
%condition1 first, then condition 2, etc. for now, each condition contains
%the same number of rows.

%starts fill in the data from the first column
AllForceTrialAverage((conditionType-1)*forceRowsAverage+1 : conditionType*forceRowsAverage,...
    1:alignedCols) = trialAverageForceData(:,:);

AllForceAlignedByBlock((conditionType-1)*blocks*forceRows+1 : conditionType*blocks*forceRows,...
    1:alignedCols) = alignedForceByBlock(:,:);

for i = 1 : blocks*forceRows
    toFill = forceByBlock(i,:);
    nonZeroCols = length(nonzeros(toFill));
    AllForceRaw((conditionType-1)*blocks*forceRows+i,1:nonZeroCols) = forceByBlock(i,1:nonZeroCols);
end

AllPosTrialAverage((conditionType-1)*posRows+1 : conditionType*posRows,...
    1:alignedCols) = trialAveragePosData(:,:);

AllPosAlignedByBlock((conditionType-1)*blocks*posRows+1 : conditionType*blocks*posRows,...
    1:alignedCols) = alignedPosByBlock(:,:);
for i = 1 : blocks*posRows
    toFill = posByBlock(i,:);
    nonZeroCols = length(nonzeros(toFill));
    AllPosRaw((conditionType-1)*blocks*posRows+i,1:nonZeroCols) = posByBlock(i,1:nonZeroCols);
end


AllVelTrialAverage((conditionType-1)*posRows+1 : conditionType*posRows,...
    1:alignedCols) = trialAverageVelData(:,:);

AllVelAlignedByBlock((conditionType-1)*blocks*posRows+1 : conditionType*blocks*posRows,...
    1:alignedCols) = alignedVelByBlock(:,:);

for i = 1 : blocks*posRows
    toFill = velByBlock(i,:);
    nonZeroCols = length(nonzeros(toFill));
    AllVelRaw((conditionType-1)*blocks*posRows+i,1:nonZeroCols) = velByBlock(i,1:nonZeroCols);
end

for i = 1 : blocks
    nonZeroCols = length(nonzeros(trialIndexForCurrentCondition(i, :)));
    AllDataIndex((conditionType-1)*blocks+i,1:nonZeroCols) = trialIndexForCurrentCondition(i, 1:nonZeroCols);
end

AllTime(conditionType, 1:length(time)) = time;
