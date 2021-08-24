%% preliminary stuff: get the directory we're in
% and add the proper subdirectories to the path
cpath = which('strflab_graddesc');
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
cellDirName = 'l2a_good';
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
strfLength = 40;
strfDelays = 0:(strfLength-1);
modelParams = linInit(nStimChannels, strfDelays, 'linear');     

% convert srData into a format that strflab understands
[allstim, allresp, groupIndex] = srdata2strflab(srData, 0);

% Normalize the stimulus
[allstimzs, s_stds, s_means] = norm_std_mean(allstim);

% put stimulus and response and group assignments into global structure
strfData(allstimzs, allresp, groupIndex);

%specify training and validation indicies
% The training index will be changed again for early stopping.
% Note that you could do this by Jackknifing by adding an extra loop.
trainSets = 1:18;
holdOutSets = [19 20];
tIndx = zeros(size(allresp));
for k = 1:length(trainSets)    
    tIndx = tIndx | (groupIndex == trainSets(k));    
end
trainingIndex = find(tIndx); 

% Find the mean response to initialize the bias.
meanResp = mean(allresp(trainingIndex));
% This next line depends on the model. 
% modelParams.b1 = log(meanResp);      % This is for exponential 
 modelParams.b1 = meanResp;         % This is for linear

tIndx = zeros(size(allresp));
for k = 1:length(holdOutSets)    
    tIndx = tIndx | (groupIndex == holdOutSets(k));    
end
validationIndex = find(tIndx);

%% Run the direct fit (Ridge Regression)
optOptions = trnDirectFit();
optOptions.display = 1;

% The direct fit actually divides the trainingIndex into STRF fitting and
% ridge parameter fitting by Jackknifing
[modelParamsDirectFit, options] = strfOpt(modelParams, trainingIndex, optOptions);

%% Now calculate second order stimulus

% First get spectrograms at a coarser resolution

preprocDir = fullfile(cellStimDir, 'preproc_sparse'); %cache the preprocessed data here
preprocOptions = struct; 
preprocOptions.fband = 500; 
preprocOptions.nstd = 6; 
preprocOptions.low_freq = 1000; 
preprocOptions.high_freq = 8000; 
preprocOtions.log = 1;
[s,mess,messid] = mkdir(preprocDir);

srData = preprocess_sound(stimFiles, respFiles, 'ft', preprocOptions, preprocDir);
pairCount = length(srData.datasets); %# of stim/response pairs

% convert srData into a format that strflab understands
[allstimsparse, allresp, groupIndex] = srdata2strflab(srData, 0);


%% run gradient descent with a hard limit on the # of iterations
% create default optimization options structure for gradient descent

maxIterations = 1000;   % This will be true for all our runs

optOptions = trnGradDesc();
optOptions.display = 1;
optOptions.maxIter = maxIterations;
optOptions.stepSize = 0.01;           % If the gradient is normalized this is the same size step is taken at each iteration. 
optOptions.gradNorm = 0;              % New flag.  The default is 1 for normalizing. 
modelParams.freqDomain = 0;           % This should be in the optOptions?  It also looks broken

% Run the optimization.  You will want to check the command window to make
% sure that the gradient and error decrease.  If they do not, you should
% try a smaller step size (or you might want to try normalizing the
% gradient in which case the step size will be in units of your parameters
[modelParamsGradDesc, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);


%% run coordinate descent with a hard limit on the # of iterations
% create default optimization options structure for gradient descent
optOptions = trnGradDesc();

optOptions.display = 1;
optOptions.maxIter = maxIterations;
optOptions.stepSize = 1;            % because we are only going in one direction a bigger step works better here
optOptions.coorDesc = 1;
optOptions.gradNorm = 0;              % New flag.  The default is 1 for normalizing. 
modelParams.freqDomain = 0;           % This should be in the optOptions?  It also looks broken


[modelParamsCoorDesc, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);

%% run gradient descent with an early stopping set

%specify training and stopping indicies
trainSets = 1:16;
stopSets = [17 18];
holdOutSets = [19 20];

tIndx = zeros(size(allresp));
for k = 1:length(trainSets)    
    tIndx = tIndx | (groupIndex == trainSets(k));    
end
trainingIndex = find(tIndx);

tIndx = zeros(size(allresp));
for k = 1:length(stopSets)
    tIndx = tIndx | (groupIndex == stopSets(k));
end
stopIndex = find(tIndx);
% Note that the validation set stayed the same.

optOptions = trnGradDesc();
optOptions.display = 1;
optOptions.maxIter = maxIterations;
optOptions.earlyStop = 1;
optOptions.errLastN = 30;  % This is the default
optOptions.errSlope = -1.0000e-003;  % This is also the default
optOptions.stepSize = 0.01;           % If the gradient is normalized this is the same size step is taken at each iteration. 
optOptions.gradNorm = 0;              % New flag.  The default is 1 for normalizing. 
modelParams.freqDomain = 0;           % This should be in the optOptions?  It also looks broken


[modelParamsGradDescES, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, stopIndex);


%% run coordinate descent with an early stopping set
optOptions = trnGradDesc();
optOptions.display = 1;
optOptions.maxIter = maxIterations;
optOptions.stepSize = 1;
optOptions.earlyStop = 1;
optOptions.coorDesc = 1;
optOptions.errLastN = 30;  % This is the default
optOptions.errSlope = -1.0000e-003;  % This is also the default
optOptions.gradNorm = 0;              % New flag.  The default is 1 for normalizing. 
modelParams.freqDomain = 0;           % This should be in the optOptions?  It also looks broken



[modelParamsCoorDescES, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, stopIndex);


%% make plots of STRFs
figure(1); 

subplot(5, 1, 1); 
strf_real = modelParamsDirectFit.w1;
for k=1:nStimChannels
    strf_real(k,:) = strf_real(k,:)./s_stds(k);
end
imagesc(strf_real);
axis tight;
absmax = max(max(abs(strf_real)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Direct Fit: bias=%f', modelParamsDirectFit.b1));

subplot(5, 1, 2);
strf_real = squeeze(modelParamsGradDesc.w1);
for k=1:nStimChannels
    strf_real(k,:) = strf_real(k,:)./s_stds(k);
end
imagesc(strf_real);
axis tight;
absmax = max(max(abs(strf_real)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Grad Desc: bias=%f', modelParamsGradDesc.b1));

subplot(5, 1, 3); 
strf_real = squeeze(modelParamsCoorDesc.w1);
for k=1:nStimChannels
    strf_real(k,:) = strf_real(k,:)./s_stds(k);
end
imagesc(strf_real);
axis tight;
absmax = max(max(abs(strf_real)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Coord Desc: bias=%f', modelParamsCoorDesc.b1));

subplot(5, 1, 4); 
strf_real = squeeze(modelParamsGradDescES.w1);
for k=1:nStimChannels
    strf_real(k,:) = strf_real(k,:)./s_stds(k);
end
imagesc(strf_real);
axis tight;
absmax = max(max(abs(strf_real)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Grad Desc + Early Stopping: bias=%f', modelParamsGradDescES.b1));

subplot(5, 1, 5); 
strf_real = squeeze(modelParamsCoorDescES.w1);
for k=1:nStimChannels
    strf_real(k,:) = strf_real(k,:)./s_stds(k);
end
imagesc(strf_real);
axis tight;
absmax = max(max(abs(strf_real)));
caxis([-absmax absmax]);
colorbar;
title(sprintf('Coord Desc + Early Stopping: bias=%f', modelParamsCoorDescES.b1));


%% compute prediction for each stim/response pair

directFitResps = predict_responses(modelParamsDirectFit);
gradDescResps = predict_responses(modelParamsGradDesc);
coorDescResps = predict_responses(modelParamsCoorDesc);
gradDescESResps = predict_responses(modelParamsGradDescES);
coorDescESResps = predict_responses(modelParamsCoorDescES);

[concatPsthHalf1_train, concatPsthHalf2_train] = concat_and_split_response(srData, trainSets);
[concatPsthHalf1_valid, concatPsthHalf2_valid] = concat_and_split_response(srData, holdOutSets);


infoFreqCutoff = 90; %90 Hz
infoWindowSize = 0.500; %500ms
numTrials = 20;

%% compute coherence and information values for each

concatPredResp_train = concat_predicted_response(directFitResps, trainSets);
concatPredResp_valid = concat_predicted_response(directFitResps, holdOutSets);
[cBoundTrain, cDirectFitTrain] = compute_coherence_full(concatPredResp_train, allresp(trainingIndex), concatPsthHalf1_train, concatPsthHalf2_train, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);
[cBoundValid, cDirectFitValid] = compute_coherence_full(concatPredResp_valid, allresp(validationIndex), concatPsthHalf1_valid, concatPsthHalf2_valid, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);

concatPredResp_train = concat_predicted_response(gradDescResps, trainSets);
concatPredResp_valid = concat_predicted_response(gradDescResps, holdOutSets);
[cBoundTrain, cGradDescTrain] = compute_coherence_full(concatPredResp_train, allresp(trainingIndex), concatPsthHalf1_train, concatPsthHalf2_train, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);
[cBoundValid, cGradDescValid] = compute_coherence_full(concatPredResp_valid, allresp(validationIndex), concatPsthHalf1_valid, concatPsthHalf2_valid, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);

concatPredResp_train = concat_predicted_response(coorDescResps, trainSets);
concatPredResp_valid = concat_predicted_response(coorDescResps, holdOutSets);
[cBoundTrain, cCoorDescTrain] = compute_coherence_full(concatPredResp_train, allresp(trainingIndex), concatPsthHalf1_train, concatPsthHalf2_train, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);
[cBoundValid, cCoorDescValid] = compute_coherence_full(concatPredResp_valid, allresp(validationIndex), concatPsthHalf1_valid, concatPsthHalf2_valid, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);


concatPredResp_train = concat_predicted_response(gradDescESResps, trainSets);
concatPredResp_valid = concat_predicted_response(gradDescESResps, holdOutSets);
[cBoundTrain, cGradDescESTrain] = compute_coherence_full(concatPredResp_train, allresp(trainingIndex), concatPsthHalf1_train, concatPsthHalf2_train, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);
[cBoundValid, cGradDescESValid] = compute_coherence_full(concatPredResp_valid, allresp(validationIndex), concatPsthHalf1_valid, concatPsthHalf2_valid, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);

concatPredResp_train = concat_predicted_response(coorDescESResps, trainSets);
concatPredResp_valid = concat_predicted_response(coorDescESResps, holdOutSets);
[cBoundTrain, cCoorDescESTrain] = compute_coherence_full(concatPredResp_train, allresp(trainingIndex), concatPsthHalf1_train, concatPsthHalf2_train, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);
[cBoundValid, cCoorDescESValid] = compute_coherence_full(concatPredResp_valid, allresp(validationIndex), concatPsthHalf1_valid, concatPsthHalf2_valid, srData.respSampleRate, numTrials, infoFreqCutoff, infoWindowSize);


%% plot training coherences
figure(2); hold on;
plot(cBoundTrain.f, cBoundTrain.c, 'k-', 'LineWidth', 2);
plot(cDirectFitTrain.f, cDirectFitTrain.c, 'b-');
plot(cGradDescTrain.f, cGradDescTrain.c, 'c-');
plot(cGradDescESTrain.f, cGradDescESTrain.c, 'c--');
plot(cCoorDescTrain.f, cCoorDescTrain.c, 'm-');
plot(cCoorDescESTrain.f, cCoorDescESTrain.c, 'm--');
legend('Upper Bound', 'Direct Fit', 'Grad Desc', 'Grad Desc + ES', 'Coord Desc', 'Coord Desc + ES');
title('Training Coherences');
axis([min(cBoundTrain.f), max(cBoundTrain.f), 0, 1]);


%% plot validation coherences
figure(3); hold on;
plot(cBoundValid.f, cBoundValid.c, 'k-', 'LineWidth', 2);
plot(cDirectFitValid.f, cDirectFitValid.c, 'b-');
plot(cGradDescValid.f, cGradDescValid.c, 'c-');
plot(cGradDescESValid.f, cGradDescESValid.c, 'c--');
plot(cCoorDescValid.f, cCoorDescValid.c, 'm-');
plot(cCoorDescESValid.f, cCoorDescESValid.c, 'm--');
legend('Upper Bound', 'Direct Fit', 'Grad Desc', 'Grad Desc + ES', 'Coord Desc', 'Coord Desc + ES');
title('Validation Coherences');
axis([min(cBoundValid.f), max(cBoundValid.f), 0, 1]);

%% Make a bar plot of information values
figure(4);
ybardata = [cDirectFitTrain.info/cBoundTrain.info cGradDescTrain.info/cBoundTrain.info ...
            cGradDescESTrain.info/cBoundTrain.info cCoorDescTrain.info/cBoundTrain.info ...
            cCoorDescESTrain.info/cBoundTrain.info; ...
            cDirectFitValid.info/cBoundValid.info cGradDescValid.info/cBoundValid.info ...
            cGradDescESValid.info/cBoundValid.info cCoorDescValid.info/cBoundValid.info ...
            cCoorDescESValid.info/cBoundValid.info ];
        
 bh = bar(ybardata', 'grouped');
legend('Training', 'Validation');
title('Goodness of Fit');
ylabel('% Total Info');
set( get(bh(1), 'Parent'), 'XTickLabel', ['DF '; 'GD '; 'GDE'; 'CD '; 'CDE']);
        
