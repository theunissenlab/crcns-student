function predictedResponses = predict_responses(modelParamsTrained)

    global globDat;

    allresp = globDat.resp;
    groupIndex = globDat.groupIdx;
    
    pairCount = length(unique(groupIndex));
    predictedResponses = cell(pairCount, 1);

    for k = 1:pairCount

      %get stim and response
      strflabIndx = find(groupIndex == k);
      resp = allresp(strflabIndx);

      %compute prediction
      [modelParamsTrained, predResp] = strfFwd(modelParamsTrained, strflabIndx);

      %fix any NaNs in response
      predResp(isnan(predResp)) = 0;

      %rectify response
      predResp(predResp < 0) = 0;

      %scale predicted response
      predResp = (predResp / max(predResp)) * max(resp);
      
      predictedResponses{k} = predResp;
      
    end
   