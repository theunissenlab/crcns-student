%% Convolve a stimulus with a STRF
%
%   Input:
%       allstim:
%
%       delays:
%
%       strf:
%
%       groupIndex:
%
%   Output:
%
%       modelResponse:
%
%
function modelResponse = conv_strf(allstim, delays, strf, groupIndex)

    nDatasets = length(unique(groupIndex));
    timeLen = size(allstim, 1);
    a = zeros(timeLen, 1);
    
    for k = 1:nDatasets
        
        rng = find(groupIndex == k);
        soff = rng(1) - 1;
        stim = allstim(rng, :);
        for ti = 1:length(delays)
            at = stim * strf(:, ti);

            thisshift = delays(ti);
            if thisshift >= 0
                a(soff+thisshift+1:end) = a(soff+thisshift+1:end) + at(1:end-thisshift);
            else
                offset = mod(thisshift, timeLen);
                a(soff:offset) = a(soff:offset) + at(-thisshift+1:end);
            end
        end
    end
    
    modelResponse = a';
    