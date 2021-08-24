%% set up strflab
% cd('~/Dropbox/Neural Data Analysis/tutorials/directfit_tutorial')

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

load([dataDir '/ecog_data_JH13-1.mat']);
 
% initialize linear model
strfLength = 100;
strfDelays = 0:(strfLength-1);
modelParams = linInit(size(trial_data(1).stim,2), strfDelays);     

% prepare ecog data in strflab format

elec = trial_data(1).chans(1);

allstim = [];
allresp = [];
groupIndex = [];

for i = 1:length(trial_data)
    %upsample to make it work
    %allstim = [allstim; resample(double(trial_data(i).stim),10,1)];
    %allresp  = [allresp resample(double(trial_data(i).resp(:,elec)),10,1)'];
    %groupIndex = [groupIndex i*ones(1,length(trial_data(i).stim)*10)];
    
    allstim = [allstim; trial_data(i).stim];
    allresp  = [allresp trial_data(i).resp(:,elec)'];
    groupIndex = [groupIndex i*ones(1,length(trial_data(i).stim))];
end
%%
% put stimulus and response and group assignments into global structure
strfData(allstim, allresp, groupIndex);
 
% create default optimization options structure
optOptions = trnDirectFit();
optOptions.display = 1;
optOptions.stimSampleRate=trial_data(1).fs;
optOptions.respSampleRate=trial_data(1).fs;
optOptions.infoFreqCutoff=trial_data(1).fs*0.4;

datIdx = 1:length(allresp); %the indexes of training data (all of it)
[modelParamsTrained, options] = strfOpt(modelParams, datIdx, optOptions);
