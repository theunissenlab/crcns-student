%% Takes a cell array of spike time vectors (one cell for each trial), and
%  converts it to a PSTH. It also splits the trials in half, and creates a
%  PSTH for each half.
%   stimLengthMs: The length of the stimlulus in milliseconds
%
%  Returns a struct psthdata, where:
%       psthdata.psth: PSTH from all trials
%       psthdata.psth_half1: PSTH from even trials
%       psthdata.psth_half2: PSTH from odd trials
%  
function psthdata = split_psth(spikeTrials, stimLengthMs)

    halfSize = round(length(spikeTrials)/2);
    spikeTrials1 = cell(halfSize, 1);
    spikeTrials2 = cell(halfSize, 1);
    for j = 1:length(spikeTrials)
       indx = round(floor(j/2)) + mod(j, 2);
       if mod(j, 2) == 0
           spikeTrials1{indx} = spikeTrials{j};
       else
           spikeTrials2{indx} = spikeTrials{j};
       end
    end

    psth = make_psth(spikeTrials, stimLengthMs, 1);  
    psth1 = make_psth(spikeTrials1, stimLengthMs, 1);
    psth2 = make_psth(spikeTrials2, stimLengthMs, 1);
    
    psthdata.psth = psth;
    psthdata.psth_half1 = psth1;
    psthdata.psth_half2 = psth2;
    