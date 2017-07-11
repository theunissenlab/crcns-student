%% Download the latest version of strflab from strflab.berkeley.edu 
% Modify the line below to include the strflab path.
addpath(genpath('C:\Users\Frederic\Documents\MATLAB\strflab_v1.45'));

%% Load voxel data. 
% These are 2322 voxels with 200 time slices with a 2 sec TR.  
% 20 images belonging to 1 of 4 classes are shown in block trials.
% The images are on for .3 seconds and off for .5 seconds. 20*.8 s = 16 second blocks.
% Since TR=2, there are exactly 8 images per block.
% This data has been polynomial detrended and motion corrected.
load ../data/logisticdata.mat
nt = length(trials);                  % Number of time slices = 200
TR = 2;  % TR in seconds

% Make indicies for the different stimulus types: blanks, places, faces, objects
blankidx = find(trials==0);
placesidx = find(trials==1);
facesidx = find(trials==2);
objectidx = find(trials==3);

% Exercise 1.  Plot the trials to understand how the experiment was
% performed.  You might want to try plotting using plot(trials) and
% imagesc(trials').


%% We are going to test for a binary stimulus classification of faces (1) and non-faces (0). 
%  After you run finish the exercise you will try other classifications!
stimclass = zeros(1,nt);
stimclass(facesidx) = 1;

%% We are now going to use strflab to perform the logistic regresssion
% Declare the global Data variable that strflab uses
global globDat;

% Put data in global variable. We are going to be predicting stimulus class from voxel responses, so the 
% voxel responses should go in the "Stimulus" variable and stimulus class goes in the "response" variable
strfData(voxeldata', stimclass');

% Note that delay is negative because voxel responses come after the
% stimulus class!
ndel = 8;      % 8 time slices or 16 seconds - this is also equal to the time the stimulus is on?
strf = linInit(size(globDat.stim,2), -ndel:0, 'logistic');

% Create an options structure
options = trnGradDesc;
options.coorDesc=0;
options.earlyStop=0;
options.stepSize = 1e-2;
options.display = 1;

% We are going to train for the first 150 samples and predict on the last 50 samples
trainIdx = 1:150;
valIdx = 151:200;

%% Fit the strf
strfTrained = strfOpt(strf, trainIdx, options);

%% Use the fit strf to make a prediction on the validation set
[strfTrained, pred] = strfFwd(strfTrained, valIdx);

%% Visualize prediction and actual stimulus class
% Exercise 2.  On figure 2, show the prediction and the actual image
% classes.


%% We can also plot the weight distribution
figure(3);
sw = sum(abs(strfTrained.w1(:,1,:)),3);
plot(sw);

%% And find the voxel with the highest weights to plot the time course

idx = find(sw == max(sw));
figure(4);
t=0:ndel;
t=t*TR;
plot(t,squeeze(strfTrained.w1(idx,1,:)));
xlabel('Time (s)');
ylabel('Weight');
title(sprintf('Weight time course for voxel %d\n', idx));

% Clearly to model is working somewhat, but it is predicting faces where there aren't any.
%% Exercise 3. Repeat the previous using early stopping (see our gradient descent tutorials)
% Remember that you will also have to specify a stopping index. After
% fitting the model, 1) use strfFwd to obtain predictions, 2) plot the
% results (fugure 5), 3) plot the weight distribution as above. (figure
% 6) 4) plot the time course of the weights as above (figure 7);
% 


%% Exercise 4.  Repeat the fit using coordinate descent.  Plot the fit on 
% figure 8 and the weight distribution on figure 9 and the time course of figure 10.  


%% Exercise 5. On figure 11, for the voxel with the largest weight,
% plot the average response (for the 8 time points) 
% in blocks containing faces and compare that to the average response.



