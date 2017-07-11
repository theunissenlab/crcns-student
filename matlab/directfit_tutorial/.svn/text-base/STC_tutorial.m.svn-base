%% For natural stimuli use:
load /auto/k2/share/strflabGOLD/fakedata/mov.mat
rawStim = single(mov(1:10,1:10,1:15000));

%% For white noise use:
% rawStim = single(randn([10 10 15000]));

%% let's create some fake data for a simple cell with a 2D Gabor filter
gparams = [.5 .5 0 3.5 0 .15 .3 0]';
[gabor, gabor90] = make3dgabor([10 10 1], gparams);
gabor = reshape(gabor, [10*10 1]);
gabor90 = reshape(gabor90, [10*10 1]);
rawStim = reshape(rawStim, [10*10 15000]);
resp = dotdelay(gabor, rawStim);
resp90 = dotdelay(gabor90, rawStim);

resp = sqrt(resp.^2 + resp90.^2);

%% Set a spiking threshold
resp(resp<.4*max(resp(:)))=0;
resp(resp>0)=1;

%% Find indicies of spikes
respidx = find(resp==1);

%% Get spike triggered stimuli
STstim = rawStim(:,respidx);

%% Calculate the spike triggered average
STA = mean(STstim,2);

%% If using natural stimuli, you should compute the corrected STA,
%% which has the stimulus correlations removed as in regression
STS = STstim*STstim';
invSTS = inv(STS);

cSTA = invSTS*STA;

%% See that it's garbage
figure(1);
subplot(1,2,1)
imagesc(reshape(STA, [10 10])); axis image; axis off; colormap gray;
subplot(1,2,2)
imagesc(reshape(cSTA, [10 10])); axis image; axis off; colormap gray;

%% Subtract cSTA from the spike trigged stimuli
STstim2 = bsxfun(@minus, STstim, cSTA);

%% Calculate covariance matrix of spike triggered stimuli
STC = STstim2*STstim2';

%% Calculate svd of STC matrix
[u,s,v] = svd(STC);

%% Plot showing STA, Filters, and first 2 eigenvectors of STC matrix
figure(2)
subplot(2,3,2)
imagesc(reshape(gabor, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,3)
imagesc(reshape(gabor90, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,4)
imagesc(reshape(cSTA, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,5)
imagesc(reshape(u(:,1), [10 10])); axis image; axis off; colormap gray;
subplot(2,3,6)
imagesc(reshape(u(:,2), [10 10])); axis image; axis off; colormap gray;


%% Now let's try STC as a regression

%% Subtract STA from stimuli
rawStim2 = bsxfun(@minus, rawStim, cSTA);

%% Make an empty matrix for the stimulus covariance at every frame
covStim = zeros([15000 100 100], 'single');

%% Calculate covariance of every frame
for ii = 1:15000
    covStim(ii,:,:) = rawStim2(:,ii)*rawStim2(:,ii)';
end

%% Reshape covariance of stimuli
covStim = reshape(covStim, [15000 10000]);

%% Calculate stimulus time response
STR = covStim'*resp;

%% Show that the STC matrix is just the uncorrected regression, i.e.
%% stimulus times response
STC2 = reshape(mean(STR,2), [100 100]);

[u2,s2,v2] = svd(STC2);

figure(3)
subplot(2,3,2)
imagesc(reshape(gabor, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,3)
imagesc(reshape(gabor90, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,4)
imagesc(reshape(cSTA, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,5)
imagesc(reshape(u2(:,1), [10 10])); axis image; axis off; colormap gray;
subplot(2,3,6)
imagesc(reshape(u2(:,2), [10 10])); axis image; axis off; colormap gray;

%% To do the full regression you need to invert the covariance
%% of covariance matrix. For white noise, it's approximately the
%% identity matrix so this can be ignored. But for natural stimuli
%% this may be too large for you computer !!!!

%% calulate covariance of covariance terms
% STS = covStim*covStim';

% This is what you would do but look at the dimensions, it's huge!
% [u3,s3,v3] = svd(STS);
% invSTS = v*diag(1/diag(s3))*u3';
% STC3 = invSTS*STR;


clear rawStim rawStim2

%% Now the strflab way

%% Declare global data vairable
global globDat

%% Initialize strf to contain weights for each covariance term and no delay
strf = linInit(10000, [0]);

%% Set bias to be mean response
strf.b1 = mean(resp);

%% Put stimulus covariance and response in global variable
strfData(covStim, resp)

%% Set up options structure for fitting
options = trnSCG;
options.display = 1;
options.maxIter = 100; 

%% Set index to be all the data
trainingIdx = [1:globDat.nSample];

%% Fit the model
strfTrained = strfOpt(strf, trainingIdx, options);


%% look at the STC matrix solved by strflab
STC3 = strfTrained.w1;

STC3 = reshape(STC3, [100 100]);

[u3,s3,v3] = svd(STC3);

figure(4)
subplot(2,3,2)
imagesc(reshape(gabor, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,3)
imagesc(reshape(gabor90, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,4)
imagesc(reshape(cSTA, [10 10])); axis image; axis off; colormap gray;
subplot(2,3,5)
imagesc(reshape(u3(:,1), [10 10])); axis image; axis off; colormap gray;
subplot(2,3,6)
imagesc(reshape(u3(:,2), [10 10])); axis image; axis off; colormap gray;

clear all


%%% Now we are going to try STC with data from a real cell

%% load the stimulus and response data
load /auto/k5/moliver/V1_complexcell.mat

%% reshape the movie
mov = reshape(mov, [144 24000]);

%% Check for NaNs in resp. First few are often NaN to ignore stimulus onset
%% effects.
nanidx = find(isnan(resp));

%% Take valid part of response and movie
resp=resp(nanidx(end)+1:end);
mov = mov(:,nanidx(end)+1:end);

%% We are going to assume the response comes one frame after the stimulus, 
%% so here we allign the stimulus and response
mov = mov(:,1:end-1);
resp = resp(2:end);

respidx = find(resp>1);

%% Get spike triggered stimuli
STstim = mov(:,respidx);

%% Because there are 4 repititions, we weight the spike triggered stimuli by 
%% the response. This is identical to making a copy of the stimulus frame for
%% each spike
STstim = STstim .* resp(respidx);

%% Calculate the spike triggered average and divide by 4 for the number of reps
STA = mean(STstim,2)/4;

%% We are using natural stimuli, so you should compute the corrected STA,
%% which has the stimulus correlations removed as in regression
STS = STstim*STstim';
invSTS = inv(STS);

cSTA = invSTS*STA;

%% See that it's garbage
figure(4);
subplot(1,2,1)
imagesc(reshape(STA, [12 12])); axis image; axis off; colormap gray;
subplot(1,2,2)
imagesc(reshape(cSTA, [12 12])); axis image; axis off; colormap gray;

%% Subtract STA from the spike trigged stimuli
STstim2 = bsxfun(@minus, STstim, cSTA);

%% Calculate covariance matrix of spike triggered stimuli
STC = STstim2*STstim2';

%% Calculate svd of STC matrix
[u4,s4,v4] = svd(STC);

%% Plot showing STA, Filters, and first 2 eigenvectors of STC matrix
figure(5)
subplot(1,3,1)
imagesc(reshape(cSTA, [12 12])); axis image; axis off; colormap gray;
subplot(1,3,2)
imagesc(reshape(u4(:,1), [12 12])); axis image; axis off; colormap gray;
subplot(1,3,3)
imagesc(reshape(u4(:,2), [12 12])); axis image; axis off; colormap gray;


%% Now let's try STC as a regression

%% Subtract STA from stimuli
rawStim2 = bsxfun(@minus, mov, cSTA);

%% Make an empty matrix for the stimulus covariance at every frame
covStim = zeros([size(resp,1) 144 144], 'single');

%% Calculate covariance of every frame
for ii = 1:size(resp,1)
    covStim(ii,:,:) = rawStim2(:,ii)*rawStim2(:,ii)';
end

%% Reshape covariance of stimuli
covStim = reshape(covStim, [size(resp,1) 20736]);

%% Calculate stimulus time response
STR = covStim'*resp;

%% Show that the STC matrix is just the uncorrected regression, i.e.
%% stimulus times response
STC2 = reshape(nanmean(STR,2), [144 144]);

[u5,s5,v5] = svd(STC2);

figure(6)
subplot(1,3,1)
imagesc(reshape(cSTA, [12 12])); axis image; axis off; colormap gray;
subplot(1,3,2)
imagesc(reshape(u5(:,1), [12 12])); axis image; axis off; colormap gray;
subplot(1,3,3)
imagesc(reshape(u5(:,2), [12 12])); axis image; axis off; colormap gray;

%% To do the full regression you need to invert the covariance
%% of covariance matrix. For white noise, it's approximately the
%% identity matrix so this can be ignored. But for natural stimuli
%% this may be too large for you computer !!!!

%% calulate covariance of covariance terms
% STS = covStim*covStim';

% This is what you would do but look at the dimensions, it's huge!
% [u3,s3,v3] = svd(STS);
% invSTS = v*diag(1/diag(s3))*u3';
% STC3 = invSTS*STR;


clear all

%% Now the strflab way

%% load the stimulus and response data
load /auto/k5/moliver/V1_complexcell.mat

%% reshape the movie
mov = reshape(mov, [144 24000]);

%% We are going to compute the covariance around the mean pixel values rather
%% than the spike triggered mean pixel values. There are arguments to be made for
%% either method, but it makes little difference
movmean = mean(mov,2);

rawStim2 = bsxfun(@minus, mov, movmean);

%% Make an empty matrix for the stimulus covariance at every frame
covStim = zeros([24000 144 144], 'single');

%% Calculate covariance of every frame
for ii = 1:size(resp,1)
    covStim(ii,:,:) = rawStim2(:,ii)*rawStim2(:,ii)';
end

%% Reshape covariance of stimuli. This will be our stimuli for the regression
covStim = reshape(covStim, [24000 20736]);

%% Declare global data variable
global globDat

%% Initialize a strf with a weight for each covariance value and use 1 delay
strf = linInit(20736, [1]);

%% Set the bias to the mean response
strf.b1 = nanmean(resp);

%% Put the stimuli and response in the global variable
strfData(covStim, resp)

%% Set options structure
options = trnSCG;
options.display = 1;

%% Try different numbers of iterations to see the effect. What changes? Why?
options.maxIter = 10; 

%% Set index to use all data
trainingIdx = [1:globDat.nSample];

%% Fit strf model
strfTrained = strfOpt(strf, trainingIdx, options);


%% look at the STC matrix solved by strflab
STC3 = strfTrained.w1;

STC3 = reshape(STC3, [144 144]);

[u6,s6,v6] = svd(STC3);

figure(7)
subplot(1,3,1)
imagesc(reshape(movmean, [12 12])); axis image; axis off; colormap gray;
subplot(1,3,2)
imagesc(reshape(u6(:,1), [12 12])); axis image; axis off; colormap gray;
subplot(1,3,3)
imagesc(reshape(u6(:,2), [12 12])); axis image; axis off; colormap gray;


