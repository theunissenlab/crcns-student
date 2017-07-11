function [allstim, allresp, groupIndex] = srdata2strflab(srData, useRaw, preprocOptions)

    if nargin < 3
        preprocOptions = struct;
    end

    pairCount = length(srData.datasets);
    totalStimLength = 0;
    totalRespLength = 0;
    
    %tally up total lengths
    numTrials = 1e10;
    for k = 1:pairCount
        ds = srData.datasets{k};
        numTrials = min(numTrials, length(ds.resp.rawSpikeIndicies));
        totalStimLength = totalStimLength + size(ds.stim.tfrep.spec, 2);
        totalRespLength = totalRespLength + length(ds.resp.psth);
    end
    
    allstim = zeros(totalStimLength, srData.nStimChannels);
    
    if useRaw       
        allresp = cell(numTrials, 1);
    else
        allresp = zeros(1, totalRespLength);
    end
    
    groupIndex = zeros(1, totalRespLength);
    
    %set up matricies    
    currentIndex = 1;
    for k = 1:pairCount
        ds = srData.datasets{k};
        
        stim = ds.stim.tfrep.spec;
        stimLen = size(stim, 2);
        resp = ds.resp.psth;
        
        if size(stim, 2) ~= length(resp)
            error('Stim and response lengths are not the same for dataset %d!\n', k);
        end
        eIndx = currentIndex + length(resp) - 1;
        rng = currentIndex:eIndx;
        
        allstim(rng, :) = stim';
        groupIndex(rng) = k;
        
        if ~useRaw    
            allresp(rng) = resp;
        else            
            for trialNum = 1:numTrials
                spikes = allresp{trialNum};
                st = ds.resp.rawSpikeIndicies{trialNum};
                st(st > stimLen) = [];
                st = st + currentIndex - 1;                
                st = st / srData.stimSampleRate;
                
                allresp{trialNum} = [spikes st(:)'];
            end
        end
        
        currentIndex = eIndx+1;
    end
    
    %do preprocessing if requested
    %% subtract mean from stim    
    if isfield(preprocOptions, 'meanSubtractStim') && preprocOptions.meanSubtractStim
        %subtract off scalar mean
        fprintf('Subtracting off mean stimuli...\n');
        numTimePoints = size(allstim, 1);
        for k = 1:numTimePoints
            allstim(k, :) = allstim(k, :) - srData.stimAvg';
        end
    end

    %% scale input by stdev
    if isfield(preprocOptions, 'scaleStim') && preprocOptions.scaleStim
        fprintf('Scaling input by std deviation...\n');
        allstim = allstim / stdev(allstim(:));
    end

    %% subtract mean from resp
    if isfield(preprocOptions, 'meanSubtractResp') && preprocOptions.meanSubtractResp
        if isfield(preprocOptions, 'tvMean') && preprocOptions.tvMean
            fprintf('Subtracing off time-varying mean rate from response...\n');
            cindx = 1;
            nstims = size(srData.tvRespAvg, 1);
            for k = 1:nstims
                tvresp = srData.tvRespAvg(k, :);
                eindx = cindx + length(tvresp) - 1;
                allresp(cindx:eindx) = allresp(cindx:eindx) - tvresp;
                cindx = eindx + 1;
            end
        else
            fprintf('Subtracing off scalar mean rate from response...\n');
            allresp = allresp - srData.respAvg;
        end
    end
    