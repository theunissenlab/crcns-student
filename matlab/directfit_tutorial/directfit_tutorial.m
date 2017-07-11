
%% preliminary stuff: get the directory we're in
% and add the proper subdirectories to the path
cpath = which('directfit_tutorial');
[rootDir, name, ext] = fileparts(cpath);
spath = fullfile(rootDir, 'strflab');
addpath(genpath(spath));
dfpath = fullfile(rootDir, 'direct_fit');
addpath(dfpath);
vpath = fullfile(rootDir, 'validation');
addpath(vpath);
ppath = fullfile(rootDir, 'preprocessing');
addpath(ppath);
dataDir = fullfile(rootDir, '..', 'data'); %contains stim/response pairs
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
cellDirName = 'l2a_avg';
cellDir = fullfile(dataDir, cellDirName);


%% now we're going to get the stimulus and response
% files from the cell directory using a function that
% was written to deal with this directory structure.
% we'll pull stim/response files for conspecific
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

%% now we're going to set up strflab
nStimChannels = srData.nStimChannels;
  
% initialize linear model
strfLength = 75;
strfDelays = 0:(strfLength-1);
modelParams = linInit(nStimChannels, strfDelays);     

% convert srData into a format that strflab understands
[allstim, allresp, groupIndex] = srdata2strflab(srData, 0);

% put stimulus and response and group assignments into global structure
strfData(allstim, allresp, groupIndex);
 
% create default optimization options structure
optOptions = trnDirectFit();
optOptions.display = 1;

%% run direct fit optimization on all of the data
% The strfOpt is the strf optimizing function.  

datIdx = 1:length(allresp); %the indexes of training data (all of it)
[modelParamsTrained, options] = strfOpt(modelParams, datIdx, optOptions);


%% compute prediction for each stim/response pair and concatenate them,
% also split the real PSTHs in half and concatenate them across pairs
pairCount = length(srData.datasets);
predictionSets = cell(pairCount, 1);
concatPredResp = [];
concatPsthHalf1 = [];
concatPsthHalf2 = [];

displayPredictions = 1;

strfFrequencies = -1;

for k = 1:pairCount
  
  %get stim and response
  ds = srData.datasets{k};
  tfrep = ds.stim.tfrep;
  strfFrequencies = tfrep.f;
  
  strflabIndx = find(groupIndex == k);
  stim = allstim(strflabIndx, :);
  resp = allresp(strflabIndx);
  
  %compute prediction
  [modelParamsTemp, predResp] = strfFwd(modelParamsTrained, strflabIndx);
  clear modelParamsTemp;
  
  %fix any NaNs in response
  predResp(isnan(predResp)) = 0;
  
  %rectify response
  predResp(predResp < 0) = 0;
  
  %scale predicted response
  predResp = (predResp / max(predResp)) * max(resp);
  
  %concatenate PSTH haves and predicted PSTH for this trial
  stimLengthMs = (size(stim, 1)/srData.stimSampleRate) * 1e3;
  psthdata = split_psth(ds.resp.rawSpikeTimes, stimLengthMs); 
  numTrials = length(ds.resp.rawSpikeTimes);
  
  concatPredResp = [concatPredResp; rv(predResp)];
  concatPsthHalf1 = [concatPsthHalf1; rv(psthdata.psth_half1)];
  concatPsthHalf2 = [concatPsthHalf2; rv(psthdata.psth_half2)];
   
  if displayPredictions
    rsint = 1 / srData.respSampleRate;
    tresp = 0:rsint:(length(resp)-1)*rsint;
    tpresp = 0:rsint:(length(predResp)-1)*rsint;
          
    %make plots...
    figure; hold on;    
    
    %plot stimulus
    subplot(2, 1, 1); hold on;
    plot_tfrep(tfrep);
    
    %plot response and prediction
    subplot(2, 1, 2); hold on;
    plot(tresp, resp, 'k-', 'LineWidth', 2);
    plot(tpresp, predResp, 'r-', 'LineWidth', 1);
    legend('Real', 'Model');
    axis tight;    
  end
end

%% display STRF
figure; hold on;
imagesc(strfDelays, strfFrequencies, modelParamsTrained.w1);
axis tight;
xlabel('Delay (ms)');
ylabel('Frequency (Hz)');
title('STRF');
absmax = max(max(abs(modelParamsTrained.w1)));
caxis([-absmax absmax]);
colorbar;


infoFreqCutoff = 90; %90 Hz
infoWindowSize = 0.500; %500ms
%% compute coherence and information values
[cBound, cModel] = compute_coherence_full(concatPredResp, allresp, concatPsthHalf1, concatPsthHalf2, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);

performanceRatio = cModel.info / cBound.info;

%% plot information
figure; hold on;
plot(cBound.f, cBound.c, 'k-', 'LineWidth', 2);
plot(cBound.f, cBound.cUpper, 'b-', 'LineWidth', 2);
plot(cBound.f, cBound.cLower, 'r-', 'LineWidth', 2);

plot(cModel.f, cModel.c, 'k--', 'LineWidth', 2);
plot(cModel.f, cModel.cUpper, 'b--', 'LineWidth', 2);
plot(cModel.f, cModel.cLower, 'r--', 'LineWidth', 2);
theTitle = sprintf('Info=%0.0f out of %0.0f bits | Ratio=%0.2f', ...
		   cModel.info, cBound.info, performanceRatio);
title(theTitle);
axis([min(cBound.f), max(cBound.f), 0, 1]);

%% plot all strfs

% Fist get the file names where you will find all the STRFs saved
strfFiles = cell(length(options.tolerances), 1);
for k = 1:length(options.tolerances)
    fname = [options.outputDir '/strfResult_Tol' num2str(k) '.mat'];
    strfFiles{k} = fname;
end

% Now choose the temporal section of the STRF that is used in the convolution
halfIndx = ceil(max(abs(modelParams.delays))) + 1;   % This is the point corresponding to zero
startIndx = halfIndx + round(min(modelParams.delays));
endIndx = halfIndx + round(max(modelParams.delays));
strfRng = startIndx:endIndx;

figure;
spvals = options.sparsenesses;
maxall = zeros( length(options.tolerances) , length(spvals) );

for k = 1:length(options.tolerances)
    
    % Load filters obtained for that tolerance
    svars = load(strfFiles{k});
    %strfsJN is an MxPxT matrix, where M=# of channels, P=# of STRF
    % delays, T = # of stim/response pairs, each element strfsJN(:, :, k) is a STRF
    % constructed from fitting all but pair k to the data.
    strfsJN = svars.STRFJN_Cell;
    
    %strfsJN_std is also an MxPxT matrix. Element strfsJN_std(:, :, k)
    % is the STRF standard deviation across the set of all jacknifed
    % STRFs excluding the kth STRF
    strfsJN_std = svars.STRFJNstd_Cell;
    
    %strfMean is the mean across all jacknifed STRFs for a given
    % tolerance
    strfMean = svars.STRF_Cell;
    
    %strfStdMean is the average standard deviation across all jacknifed
    % STRFs for a given tolerance
    strfStdMean = squeeze(mean(strfsJN_std,3));
    
    clear svars;
    
    
    for q = 1:length(spvals)
        subplot(length(options.tolerances), length(spvals), (k-1)*length(spvals)+ q );
        smoothedMeanStrf = df_fast_filter_filter(strfMean, strfStdMean, spvals(q));
        imagesc(smoothedMeanStrf(:, strfRng));
        maxall(k,q) = max(max(abs(smoothedMeanStrf(:, strfRng))));
        caxis([-maxall(k,q) maxall(k,q)]);
    end
end

% This is to make all plots have the same axis...
% maxabs = max(max(maxall));
% for k = 1:length(options.tolerances)
%     for q = 1:length(spvals)
%         subplot(length(options.tolerances), length(spvals), (k-1)*length(spvals)+ q );
%         caxis([-maxabs maxabs]);
%     end
% end




