function [concatPsthHalf1, concatPsthHalf2] = concat_and_split_response(srData, groupsToSelect)

    global globDat;
    
    allstim = globDat.stim;
    groupIndex = globDat.groupIdx;

    concatPsthHalf1 = [];
    concatPsthHalf2 = [];

    for k = 1:length(groupsToSelect)
        
      gindx = groupsToSelect(k);

      %get stim and response
      ds = srData.datasets{k};
      strflabIndx = find(groupIndex == gindx);
      stim = allstim(strflabIndx, :);      

      %concatenate PSTH haves and predicted PSTH for this trial
      stimLengthMs = (size(stim, 1)/srData.stimSampleRate) * 1e3;
      psthdata = split_psth(ds.resp.rawSpikeTimes, stimLengthMs); 

      concatPsthHalf1 = [concatPsthHalf1; rv(psthdata.psth_half1)];
      concatPsthHalf2 = [concatPsthHalf2; rv(psthdata.psth_half2)];

    end
    