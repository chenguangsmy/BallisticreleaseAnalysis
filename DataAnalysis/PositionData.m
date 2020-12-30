%Author: Shuqi Liu 
%Date: 2019-12-30 11:50 
%File Description: Pre-process the position data to be aligned at movement
%onset time. Save the processsed data to the output variable.
%Output variable format:
%- AllPositionTrialAverage: 6 rows per condition, all force for condition1, then
%next condition, etc.
%- AllPositionAlignedByBlock: total reps (blocks * rep per block) * 6 rows per condition
%In row order: Fx for all reps of a given condition, y, z, torque x, y, z; 
%then repeat for next condition
%- AllPositionRaw: not aligned at time 0, in order x,y,z,torque x,y,z for
%1 rep of a given condition; then next rep...; after all reps for this
%condition, go to next condition.

if (plotPos == 1)
    figure('Name','Position'); 
end

posByBlock = zeros(3*blocks, maxCols);
onsetIndexByBlock = zeros(1,blocks);
maxMoveTime = 0;
for i = 1 : blocks
    %per block
    com1indexPerBlock = com1Index(i, :);
    com1indexPerBlock = nonzeros(com1indexPerBlock);
    time = 0 : 0.02 : (length(com1indexPerBlock) -1) * 0.02;
    if (movementOnsetIndex(1,i) == 0)
        continue
    end
    time0Index = find(com1indexPerBlock == movementOnsetIndex(1,i));
    onsetIndexByBlock (1, i) = time0Index;
    moveTime = length(com1indexPerBlock) - time0Index;
    if (moveTime > maxMoveTime)
        maxMoveTime = moveTime;
    end
    time = time - (time0Index-1) * 0.02;
    pos = dataCurr.Data.Position.Actual(com1indexPerBlock, :);
    posSize = size(pos);
    posByBlock(i*3 - 2 : i*3, 1 : posSize(1, 1)) = pos';
    
    if (plotPos == 1)
        plotIndex = i*3-2;
        subplot(blocks,3, plotIndex);
        plot(time, pos(:,1));
        if i == 1 % temporary, not good
            title('X-Direction');
        end

        subplot(blocks,3,plotIndex+1);
        plot(time, pos(:,2));
        if i == 1
            title('Y-Direction');
        end

        subplot(blocks,3,plotIndex+2);
        plot(time, pos(:,3));
        if i == 1
            title('Z-Direction'); 
        end
    end
end

if (plotPos == 1)
    figure('Name','Position Aligned at time 0'); 
end
maxOnsetIndex = max(onsetIndexByBlock);
alignedCols = maxOnsetIndex+maxMoveTime;
%rearrange to be x1-25, y 26-50, z51-75
alignedPosByBlock = NaN(3*blocks, alignedCols);
for i = 1 : blocks
    time0Index = onsetIndexByBlock(1, i);
    %time stamp padding needed
    offset = maxOnsetIndex - time0Index +1;
    nonzeroCols = nonzeros(posByBlock(i*3-2 : i*3, :));
    nonzeroCols = length(nonzeroCols)/3;
    alignedPosByBlock(i, offset : offset+nonzeroCols - 1) = posByBlock(i*3 - 2, 1:nonzeroCols); 
    alignedPosByBlock(i+blocks, offset : offset+nonzeroCols - 1) = posByBlock(i*3 - 1, 1:nonzeroCols); 
    alignedPosByBlock(i+2*blocks, offset : offset+nonzeroCols - 1) = posByBlock(i*3, 1:nonzeroCols); 

    if plotPos == 1 
        time = 0 : 0.02 : (alignedCols-1) * 0.02;
        time = time - (maxOnsetIndex-1) * 0.02;

        plotIndex = i*3-2;
        subplot(blocks,3, plotIndex);
        plot(time, alignedPosByBlock(i, :));
        xlim([time(1), time(length(time))]);
        if i == 1 % temporary, not good
            title('X-Direction');
        end

        subplot(blocks,3,plotIndex+1);
        plot(time, alignedPosByBlock(i+blocks, :));
        xlim([time(1), time(length(time))]);
        if i == 1
            title('Y-Direction');
        end

        subplot(blocks,3,plotIndex+2);
        plot(time, alignedPosByBlock(i+2*blocks, :));
        xlim([time(1), time(length(time))]);
        if i == 1
            title('Z-Direction'); 
        end
    end
end


trialAverageData = zeros(3, alignedCols);
if (plotPos == 1)
    figure('Name' ,'Trial Average Pos');
end
time = 0 : 0.02 : (alignedCols-1) * 0.02;
time = time - (maxOnsetIndex-1) * 0.02;
for i = 1 : 3
    trialAverageData(i, :) = mean(alignedPosByBlock((i-1)*blocks+1:i*blocks, :));
    if (plotPos == 1)
        subplot(3,1,i)
        plot(time,trialAverageData(i, :));
    end
end

%in format, each row x, y and z with each column corresponds to 1 cell in
%the time variable. 
fieldName = append('trialAverage',num2str(conditionType));
S.(fieldName) = trialAverageData;
fieldName = append('time',num2str(conditionType));
S.(fieldName) = time;
%in format, x 1-25 rows, y 26-50 rows, z 51-75 rows assuming 25 trials per
%condition. Position is aligned at time 0 = movement onset. 
fieldName = append('AlignedPosByBlock', num2str(conditionType));
S.(fieldName) = alignedPosByBlock;
fieldName = append('PosByBlock', num2str(conditionType));
S.(fieldName) = posByBlock;
%The index for all trials in this condition. Each row correspond to 1
%successful trial of the condition. Initialized to 500 columns, but should
%take less and each row might have different number of elements
fieldName = append('DataIndex', num2str(conditionType));
S.(fieldName) = posByBlock;

save('Sonic00144Pos.mat', '-struct', 'S', '-append');
%if first time, don't use append
% save('Sonic00144Pos.mat', '-struct', 'S');

%now plot force



