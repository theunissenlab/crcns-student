
%% preliminary stuff: get the directory we're in
% and add the validation subdirectory to the path
cpath = which('coherence_tutorial');
[rootDir, name, ext] = fileparts(cpath);
vpath = fullfile(rootDir, 'validation');
addpath(vpath);
ppath = fullfile(rootDir, 'preprocessing');
addpath(ppath);
dataDir = fullfile(rootDir, '../../', 'data'); %contains stim/response pairs
stimsDir = fullfile(dataDir, 'all_stims'); %contains the .wav files


%% The next three sections allow you to load and visualize single unit data from
% The theunissen lab. Your goals are: 1. Get familiar with this data structure and
% 2. Load you own data in a similar structure. 
%  For the Theunissen data you can specify a directory for three brain
%  regions and three example neurons in each.:
%   'mld' is the avian auditory midbrain
%   'ov'  is the avian auditory thalamus
%   'l2a' is the avian auditory cortex
% each region has a 'good', 'avg', and 'bad' dataset,
% corresponding to the signal to noise ratio,
% quantified by information values.
cellDirName = 'l2a_good';
cellDir = fullfile(dataDir, cellDirName);


%% now we're going to get the stimulus and response
% files from the cell directory using a function that
% was written to deal with this directory structure.
% we'll will pull stim/response files for conspecific
% stimuli.  You should write your own data load function for your data.
datasets = find_datasets(cellDir, stimsDir, 'conspecific');
cellStimDir = datasets{1}.dirname;
stimFiles = datasets{1}.srPairs.stimFiles; %paths to .wav files
respFiles = datasets{1}.srPairs.respFiles; %paths to spike* files


%% now we're going to preprocess the sound stimuli by taking the
% short time fourier transform, and preprocess the raw spike
% times into PSTHs for each stim/response pair
preprocDir = fullfile(cellStimDir, 'preproc'); %cache the preprocessed data here
[s,mess,messid] = mkdir(preprocDir);
preprocOptions = struct; %we'll leave this empty and use default options

srData = preprocess_sound(stimFiles, respFiles, 'ft', struct, preprocDir);
pairCount = length(srData.datasets); %# of stim/response pairs


%% visualize the stimulus/response pairs by setting showSRPairs = 1
showSRPairs = 1;
if showSRPairs
    %go through each pair in the dataset
    for k = 1:pairCount
        ds = srData.datasets{k};
        plot_tf_resp(ds);
    end
end

%% Assignment 1: Calculate the noise for each trial and display it
% Calculate the noise using a signal generated from all trials as well
% as a signal that does not incldue the trial for wich you calculate the
% noise 
% You might not want to plot the noise for all trials (that would be 200
% plots). Plot it for maybe the first trial.  Look at the magnitude of the
% two noise estimates.  Which one is larger?


% These will be the noise and signal traces for all the stim-response pairs
noise_tot = [];          % Noise calculated with a mean signal from all trials
noise_d1_tot = [];       % Noise calculated from the delete one signal
signal_tot = [];

wind1 = hanning(31)/sum(hanning(31));   % 31 ms smoothing for plotting only

for k=1:pairCount
    resp = srData.datasets{k}.resp;
    stim = srData.datasets{k}.stim;
    ntrial = length(resp.rawSpikeTimes);  % Number of trials for this stim-response
    stimdur = stim.stimLength*1000.0;     % Stimulus duration in ms.
    binsize = 1;                          % Sampling rate in for signal and noisems.
    ndur = round(stimdur/binsize)+1;        % Length of signal and noise for this stim-response pair
    trialResp = zeros(ntrial, ndur);
    noise = zeros(ntrial, ndur);
    noise_d1 = zeros(ntrial, ndur);
    
% Transform spike arrival times into a time series at sampling rate binsize
    for itrial=1:ntrial
        
        stimes = resp.rawSpikeTimes{itrial};
        indx = ((stimes >= 0) & (stimes <= stimdur));
        
        stimes = stimes(indx);
        sindxs = round(stimes/binsize) + 1;
        % Prevent exceeding array dimensions (should not happen except
        % because of rounding off).
        for j = 1:length(sindxs)
            if sindxs(j) <= 0
                sindxs(j) = 1;
            end
            if sindxs(j) > ndur
                sindxs(j) = ndur;
            end
            trialResp(itrial, sindxs(j)) = trialResp(itrial, sindxs(j)) + 1;
        end
    end
        % YOUR CODE HERE
        
    
end


%% Assignment 2: Calculate the noise and signal psd obtained by averaging
% with all trials and with all trials minus 1.  Is the noise white?
% Check again the magnitude of the noise.

 % There are many different algorithms for estimating a power spectral
 % density.  The simplest is the periodogram that divides the time series
 % into non-overlapping chunkds of size nfft (in points) and multiplies
 % that segment with the weights given by the window.  If window is null,
 % periodogram uses a rectangular window.
 % Check out both periodogram and pwelch which does overlaping windows
fs = 1000.0;   % The sampling rate of the data is 1 kHz.

%% Assignment 3: Is the noise gaussian?
% Matlab has a nice command called histfit that generates a histogram and
% displays a probability density fit.  The default is a normal
% distribution.  Try it with and wihtou smoothing the noise with the wind1
% window.


%% Assignment 4: Calculate and display the coherence calculated from the signal to noise
% ratio.


%% We are now going to calculate the coherence and channel capacity using 
% Hsu, Borst, Theunissen methodology.  In order to do that we need
% to take the raw spike times, split them into even and odd trials,
% and compute PSTHs for each half. there's already a function to
% do this, so we'll just call it.
[oddPsths, evenPsths] = compute_psth_halves(srData);

%% the next step is to concatenate the PSTHs across all stim/response
% pairs into one big vector. we're going to do the same thing to 
% the psth halves as well. we're also going to take note of the 
% # of spike trials
psthConcat = [];
psthHalf1Concat = [];
psthHalf2Concat = [];
numStimPresentations = -1;
for k = 1:pairCount
  
  ds = srData.datasets{k};  
  numStimPresentations = length(ds.resp.rawSpikeTimes);
  psth = ds.resp.psth;
  psthConcat = [psthConcat psth];
    
  psthHalf1Concat = [psthHalf1Concat oddPsths{k}];
  psthHalf2Concat = [psthHalf2Concat evenPsths{k}]; 
  
end

%% we're going to make a copy of the concatenated PSTH and corrupt
% it with Gaussian noise, pretending it's a PSTH that comes from
% some model
noiseGain = 1e-1; %play with gain to increase or decrease PSTH corruption
noise = randn(size(psthConcat)) * noiseGain; %make some noise!
psthConcatNoisy = psthConcat + noise; %corrupt PSTH
psthConcatNoisy(psthConcatNoisy < 0) = 0; %rectify
psthConcatNoisy(psthConcatNoisy > 1) = 1; %rectify

%% finally, we're going to compute the upper bound of coherence, as
% for the cell itself, as well as the coherence between the noise-
% corrupted PSTH and actual PSTH

infoFreqCutoff = -1; %max frequency in Hz to compute coherence for
infoWindowSize = 0.500; %window size in seconds to compute coherence FFT
[cBound, cModel] = compute_coherence_full(psthConcatNoisy, psthConcat, psthHalf1Concat,...
					  psthHalf2Concat, srData.respSampleRate, numStimPresentations,...
					  infoFreqCutoff, infoWindowSize);

performanceRatio = cModel.info / cBound.info; %how well did our noisy psth do?

%% now we'll make some plots of the coherence values, solid lines
% are the upper bounds, dotted lines are noisy PSTHs
figure; hold on;
plot(cBound.f, cBound.c, 'k-', 'LineWidth', 2);
plot(cBound.f, cBound.cUpper, 'b-', 'LineWidth', 2);
plot(cBound.f, cBound.cLower, 'r-', 'LineWidth', 2);

plot(cModel.f, cModel.c, 'k--', 'LineWidth', 2);
plot(cModel.f, cModel.cUpper, 'b--', 'LineWidth', 2);
plot(cModel.f, cModel.cLower, 'r--', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Coherence');
theTitle = sprintf('Info=%0.0f bits/s out of %0.0f bits/s | Ratio=%0.0f', ...
		   cModel.info, cBound.info, performanceRatio);
title(theTitle);
axis([min(cBound.f), max(cBound.f), 0, 1]);
