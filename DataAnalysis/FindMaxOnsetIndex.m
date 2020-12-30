%Author: Shuqi Liu 
%Date: 2019-12-27 17:06 
%File Description: Find and set the maxOnsetIndex variable. 
%maxOnsetIndex = the max index (samplings) from
%beginning of trial to movement onset. The max is found for all trials in
%one experiment. Used later to line up trials at movement onset time in a
%matrix.

%set up the array of the index success trials, used by other classes later.
successTrials = find(dataCurr.Data.OutcomeMasks.Success == 1);
%TODO: consider only auto success for now, later figure out how to fill in
%the indices of the manual success ones without repetition

%sentinel value. -1 indicates that the loop is looking for a new trial.
%Other positive value indicates the first index of a trial, used as the
%starting point to calclate the maxOnSetIndex.
trialBegin = -1;
%The target variable to set.
maxOnsetIndex = 0;
%TODO: not in use right now. Should capture the longest interval from
%movementonest to end of a trial.
maxMoveTime = 0;
%Not in use for now. maybe used for finding max move time.
trialNo = dataCurr.Data.TrialNo(successTrials(1));
for i = successTrials
%     t = dataCurr.Data.TrialNo(i); if (t~=trialNo)
%         endIndex =  i - 1; moveTime = endIndex - onset + 1; if (moveTime
%         > maxMoveTime)
%             maxMoveTime = moveTime;
%         end
%     end
    
    if (dataCurr.Data.TaskStateCodes.Values(i) == 1 && trialBegin == -1)
        trialBegin = i;
        trialNo = dataCurr.Data.TrialNo(i);
    end
    if (dataCurr.Data.TaskStateCodes.Values(i) == movementOnsetState && trialBegin ~= -1)
        onset = i - trialBegin + 1; %include the first index
        if (onset > maxOnsetIndex)
            maxOnsetIndex = onset;
        end
        trialBegin = -1;
    end
end
fprintf('OnsetIndexMax:%d\n',maxOnsetIndex);
fprintf('MaxMoveTime:%d\n',maxMoveTime);
