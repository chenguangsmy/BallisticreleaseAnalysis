%Author: Shuqi Liu 
%Date: 2019-12-27 17:06 
%File Description: Find and set the indexes of all repetitions of the current
%condition, save the index in #trialIndexForCurrentCondition matrix. Find
%and set the movementOnSetIndex (absolute) of all repetitions of the
%current condition.

%function CheckDirectionReading(sessionNumber)

%The indexes of the trials for com1 trials. Each row corresponds to the
%indexes of 1 trial (all times points). number of rows = total number of repetitions of one
%condition. 
trialIndexForCurrentCondition = zeros(blocks,maxCols);
%The index (absolute) of the first movement onset time for each trial.
movementOnsetIndex = zeros(1, blocks);

for i = successTrials
    check = find(conditionIndex == dataCurr.Data.ComboNo(i));
    if (~isempty(check))
        blockNum = dataCurr.Data.BlockNo(i);
        if (blockNum > numberBlocks)
            %skip the incomplete blocks.
            %TODO: not ideal, for now only uses complete blocks which
            %require all conditions to be run the same number of times.
            %should try to use as many reps as possible for each condition.
            continue
        end
        indexToUpdate = (blockNum -1)*repPerBlock + check(1,1);
        %only update the first appearance of state 4 (check for 0 to update
        %it only once).
        %TODO: this can be done in the FindMaxOnsetIndex
        if (movementOnsetIndex(1,indexToUpdate) == 0 &&...
                dataCurr.Data.TaskStateCodes.Values(i) == movementOnsetState)
            movementOnsetIndex(1,indexToUpdate) = i;
        end
        for j = 1:maxCols
            %fill in the next available space. 
            %TODO: this is very inefficient, should be able to just find
            %the start and end index of a trial and copy the interval in 1
            %operation. 
            if (trialIndexForCurrentCondition(indexToUpdate,j) == 0)
                trialIndexForCurrentCondition(indexToUpdate,j) = i;
                break;
            end
        end
    end
end

%Plots to check how many trials are present for this condition
%And if each trial went through the correct 1-7 states
if (initPlot == 1) 
    temp = nonzeros(trialIndexForCurrentCondition);
    temp = sort(temp);
    com1States = dataCurr.Data.TaskStateCodes.Values(temp);
    figure;
    subplot(2,1,1);
    stateP = plot(com1States);
    title('States');

    com1Trials = dataCurr.Data.TrialNo(temp);
    subplot(2,1,2);
    TrialsP = plot(com1Trials);
    title('Trial Number');
end


