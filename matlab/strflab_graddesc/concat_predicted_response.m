function concatPredResp = concat_predicted_response(predictedResponses, groupsToSelect)

    concatPredResp = [];

    for k = 1:length(groupsToSelect)        
        gk = groupsToSelect(k);
        concatPredResp = [concatPredResp; rv(predictedResponses{gk})];            
    end
    