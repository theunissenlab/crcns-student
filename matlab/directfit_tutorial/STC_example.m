%% Spike triggered-covariance. In this demonstration,
% the spike-triggered covariance methods in the original description and
% as done in strflab (fitting the second order terms) are illustrated
% with actual response from a complex cell recorded in the gallant lab. 
%% Preliminary stuff: get the directory we're in
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

%% load the stimulus and response data and do some pre-processing
load (strcat(dataDir,'/V1_complexcell.mat'));

% reshape the stimulus movie
mov = reshape(mov, [144 24000]);

% Check for NaNs in resp. First few are often NaN to ignore stimulus onset
% effects.
nanidx = find(isnan(resp));

% Take valid part of response and movie and truncate for memory issues
maxlength = 10000;
resp=resp(nanidx(end)+1:maxlength);
mov = mov(:,nanidx(end)+1:maxlength);

% We are going to assume the response comes one frame after the stimulus, 
% so here we align the stimulus and response
mov = mov(:,1:end-1);
resp = resp(2:end);
tlength = length(resp);

% Find time bins with spikes- there are four repetitions and
% time bins can have 0,1,2,3 or 4.
respidx = find(resp>1);

%% Calculate the STA and the first order filter

%Get spike triggered stimuli
STstim = mov(:,respidx);

% The stimulus response cross-correlation or STA
STA = (STstim * resp(respidx))./(tlength*4);
   
% We are using natural stimuli, so you should compute the normalized STA,
% to obtain the linear filter
STS = STstim*STstim'./tlength;

myfilter = STS\STA;

%% Plot the STA and the linear filter.  For this cell, it is noise...
figure(1);
subplot(1,2,1)
imagesc(reshape(STA, [12 12])); axis image; axis off; colormap gray;
subplot(1,2,2)
imagesc(reshape(myfilter, [12 12])); axis image; axis off; colormap gray;

%% Calculate the spike-triggered covariance and plot the results 

%Subtract STA from the spike trigged stimuli
STstim2 = bsxfun(@minus, STstim, myfilter);

% Calculate covariance matrix of spike triggered stimuli
STC = (STstim2*STstim2')./tlength;

% Calculate svd of STC matrix
[u4,s4,v4] = svd(STC);

% Plot showing first order filter, and first 2 eigenvectors of STC matrix
figure(2);
subplot(1,3,1);
imagesc(reshape(myfilter, [12 12])); axis image; axis off; colormap gray;
title('First order filter');
subplot(1,3,2);
imagesc(reshape(u4(:,1), [12 12])); axis image; axis off; colormap gray;
title('Eigen vector of STC');
subplot(1,3,3);
imagesc(reshape(u4(:,2), [12 12])); axis image; axis off; colormap gray;


%% Now let's do the STC as a regression for the second order terms

% Subtract first order filter from stimuli
rawStim2 = bsxfun(@minus, mov, myfilter);

% Make an empty matrix for the second order terms at every time point
clear mov rawStim STstim STstim2
secondStim = zeros([tlength 144 144], 'single');

% Calculate the second order terms for every time point 
for ii = 1:tlength
    secondStim(ii,:,:) = rawStim2(:,ii)*rawStim2(:,ii)';
end

% Reshape secondorder terms of stimuli
secondStim = reshape(secondStim, [size(resp,1) 20736]);

% Calculate the cross-correlation between the second order terms and the response
STR = (secondStim'*resp)./tlength;

%% We can now plot the cross-correlation for the second order terms 
% to show that it is equal to the STC matrix.

STC2 = reshape(nanmean(STR,2), [144 144]);

[u5,s5,v5] = svd(STC2);

figure(3);
subplot(1,3,1);
imagesc(reshape(myfilter, [12 12])); axis image; axis off; colormap gray;
title('First order filter');
subplot(1,3,2);
imagesc(reshape(u5(:,1), [12 12])); axis image; axis off; colormap gray;
title('Eigen vector of 2nd order stim-resp');
subplot(1,3,3);
imagesc(reshape(u5(:,2), [12 12])); axis image; axis off; colormap gray;

%% To do the full regression, you would need to invert the covariance
% of second order terms and this is too large! 
% Calulate covariance of second order terms
% STS = (secondStim*secondStim')./tlength;  % This is too bigg...

% This is simple but wont work.
% STC3 = STS\STR;




%% Now the strflab way
clear all
% load the stimulus and response data and do some pre-processing
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

% load the stimulus and response data 
load (strcat(dataDir,'/V1_complexcell.mat'));

% reshape the stimulus movie
mov = reshape(mov, [144 24000]);

% Check for NaNs in resp. First few are often NaN to ignore stimulus onset
% effects.
nanidx = find(isnan(resp));

% Take valid part of response and movie and truncate for memory issues
maxlength = 10000;
resp=resp(nanidx(end)+1:maxlength);
mov = mov(:,nanidx(end)+1:maxlength);

tlength = length(resp);

% We are going to compute the covariance around the mean pixel values rather
% than the spike triggered mean pixel values. There are arguments to be made for
% either method, but it makes little difference
movmean = mean(mov,2);

rawStim2 = bsxfun(@minus, mov, movmean);

%% Calculate second order terms 
secondStim = zeros([tlength 144 144], 'single');

for ii = 1:tlength
    secondStim(ii,:,:) = rawStim2(:,ii)*rawStim2(:,ii)';
end
% Reshape the second order terms of stimuli. This will be our stimuli for the regression
secondStim = reshape(secondStim, [tlength 144*144]);

%% Initialize strflab 
global globDat
% Initialize a strf with a weight for each covariance value and use 1 delay
strf = linInit(144*144, [1]);

% Set the bias to the mean response
strf.b1 = nanmean(resp);

% Put the stimuli and response in the global variable
strfData(secondStim, resp)

% Set options structure
options = trnSCG;
options.display = 1;

%%  Solve by conjugate gradient descent 
% Try different numbers of iterations to see the effect. What changes? Why?
options.maxIter = 10; 

% Set index to use all data
trainingIdx = 1:globDat.nSample;

% Fit strf model
strfTrained = strfOpt(strf, trainingIdx, options);


%% Display the results and compare to non-normalized results 
STC3 = strfTrained.w1;

STC3 = reshape(STC3, [144 144]);

[u6,s6,v6] = svd(STC3);

figure(4)
subplot(1,3,1);
imagesc(reshape(movmean, [12 12])); axis image; axis off; colormap gray;
title('Mean stimulus');
subplot(1,3,2);
imagesc(reshape(u6(:,1), [12 12])); axis image; axis off; colormap gray;
title('Eigenvector of second order filter');
subplot(1,3,3);
imagesc(reshape(u6(:,2), [12 12])); axis image; axis off; colormap gray;


