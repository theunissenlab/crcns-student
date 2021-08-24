%% This exercise uses fake neurons to test the effectivemess of cross-correlation
% measures for assessing functional connectivity.

%% Generate a stimulus

tlength = 15000;            % The stimulus last 15 seconds at 1kHz sampling rate
stim = randn(1,tlength);    % It is gaussian WN


%% Let's create some fake data for two independent cells that respond to our
% stimulus 

% The first cell has an exponential filter with a 20 ms time constant
th=1:100;
h1 = exp(-th./25);
resp1 = conv(stim,h1,'same');

% The second cell has an exponential filter with a 30 ms time constant
h2 = exp(-th./30);
resp2 = conv(stim, h2, 'same');


% Threshold and set stimulus driven rms to 14 spikes/s and
% background to 1 spike/s
resp1(resp1<0) = 0.0;
resp2(resp2<0) = 0.0;
resp1 = resp1.*((0.014)./std(resp1));
resp2 = resp2.*((0.014)./std(resp2));  
resp1 = resp1 + 0.001;   % Background rate set at 1 spike/s
resp2 = resp2 +0.001;

% Exercise 1. Plot the average responses (resp1 and resp2) of the two cells 
% for the first 200 points.




%% Now we're going to generate poisson spikes from these responses.
% We are first going to generate independent spike trials
meanfr = 15;   % poisson_gen assumes that the resp is a profile and will threhold it and adjust it to meanfr
numTrials = 10;
spiketimes1 = poisson_gen_spikes(resp1, meanfr, numTrials);
spiketimes2 = poisson_gen_spikes(resp2, meanfr, numTrials);

% Generate a psth for these spike arrival times and compare to
% resp1 and resp2.
psth1 = zeros(1, length(resp1));
psth2 = zeros(1, length(resp2));

for i=1:numTrials
    trial = zeros(1, length(resp1));
    trial(spiketimes1{i}) = 1;
    psth1 = psth1 + trial;
    trial = zeros(1, length(resp1));
    trial(spiketimes2{i}) = 1;
    psth2 = psth2 + trial;
end

psth1 = psth1./numTrials;
psth2 = psth2./numTrials;

%% Exercise 2.  On the same figure 1, plot the psth obtained from your two
% neurons and compare to theoretical rate.



%% Calculate the cross-correlation.
maxlags = 100;
ccAll = zeros(1, 2*maxlags+1);

for i=1:numTrials
    trial1 = zeros(1, length(resp1));
    trial1(spiketimes1{i}) = 1;

    trial2 = zeros(1, length(resp1));
    trial2(spiketimes2{i}) = 1;
    ccAll = ccAll + xcorr(trial2, trial1, maxlags, 'Unbiased');
end

ccAll = ccAll./numTrials;

%% Exercise 3. Plot the cross-correlation on figure 2. Add labels to the x
% and y axis.  What do you see?


%% Exercise 4. Repeat this calculation after removing the mean firing rate and plot on
% figure 3. What has changed?



%% Exercise 5.  This time calculate the true cross-covariance by removing the
% time varying mean firing rate.  Repeat this calculation after removing the mean firing rate and plot on
% a new figure 4. Use the same y-scale as in figure 3. What do you observe?



%% We are now going to model two connected neurons

% Neuron 1 increases the probability of firing in neuron 2.
th=1:25;
hspike = zeros(1,27);
hspike(3:27) = exp(-th./5);   % 2 ms delay and then exponential with 5 ms decay
hspike = hspike./sum(hspike);
w12 = 1;             % Connectivity weight: 1 is a one to one - one spike causes one spike


psth2 = zeros(1, length(resp2));
clear spiketimes2;

for i=1:numTrials
    trial = zeros(1, length(resp1));
    trial(spiketimes1{i}) = 1;
    resp2_trial = w12.*conv(trial,hspike,'full');
    resp2_tot = resp2 + resp2_trial(1:length(resp2));
    spiketimes2(i) = poisson_gen_spikes(resp2_tot, meanfr, 1);
    trial = zeros(1, length(resp1));
    trial(spiketimes2{i}) = 1;
    psth2 = psth2 + trial;
end
psth2 = psth2./numTrials;

% Let's plot pairs of spike trains in the first 1000 ms.
figure(5);
for i = 1:numTrials
    subplot(numTrials, 1, i);
    hold on;
    nspikes1 = length(spiketimes1{i});
    for is=1:nspikes1
        t1 = spiketimes1{i}(is);
        if t1 > 1000
            break;
        end
        plot([t1 t1], [0 1], 'k');
    end
    nspikes2 = length(spiketimes2{i});
    for is=1:nspikes2
        t2 = spiketimes2{i}(is);
        if t2 > 1000
            break;
        end
        plot([t2 t2], [0 1], 'r');
    end
    axis([0 1000 0 1]);
    axis off;
    hold off;
end

%% Exercise 6. Calculate and plot the cross-correlation and cross-covariance 
% on the same new figure (figure 6)




%% Exercise 7. You are going to normalize the cross-covariance to obtain 
% the cross-cohenrency.  To do so you will divide by the auto-correlation
% of each spike train in the Fourier Domain.
% Plot the auto-covariance on figure 7 and the coherency on figure 8.

% Calculate the auto-correlations (time varying mean subtracted)


%% Now we are going to make neuron 2 a burster by simply adding spikes after
% each of the current spikes in a little gaussian pulse.

th=1:25;
hburst = zeros(1,40);
hburst = exp((th-20).^2./10^2);   % A gaussian pulse
hburst = hburst./sum(hburst);
w22 = 1;  % One will double the firing rate
psth2 = zeros(1, length(resp2));

for i=1:numTrials
    trial = zeros(1, length(resp1));
    trial(spiketimes2{i}) = 1;
    resp2_trial = conv(trial,hburst,'full');
    spiketimes2_added(i) = poisson_gen_spikes(resp2_trial, w22*meanfr, 1);
    trial = zeros(1, length(resp1));
    trial(spiketimes2{i}) = 1;
    trial(spiketimes2_added{i}) = 1;
    psth2 = psth2 + trial;
end
psth2 = psth2./numTrials;

% Let's plot pairs of spike trains in the first 1000 ms.
figure(9);
for i = 1:numTrials
    subplot(numTrials, 1, i);
    
    nspikes1 = length(spiketimes1{i});
    for is=1:nspikes1
        t1 = spiketimes1{i}(is);
        if t1 > 1000
            break;
        end
        plot([t1 t1], [0 1], 'k');
        hold on;
    end
    nspikes2 = length(spiketimes2{i});
    for is=1:nspikes2
        t2 = spiketimes2{i}(is);
        if t2 > 1000
            break;
        end
        plot([t2 t2], [0 1], 'r');
        hold on;
    end
    nspikes2 = length(spiketimes2_added{i});
    for is=1:nspikes2
        t2 = spiketimes2_added{i}(is);
        if t2 > 1000
            break;
        end
        plot([t2 t2], [0 1], 'g');
        hold on;
    end
    axis([0 1000 0 1]);
    axis off;
    hold off;
end

%% Exercise 8.  Repeat exercise 5 and exercise 6 for this new scenario:
% Calculate the cross-correlation, cross-covariance, auto-covariance and 
% cross-coherency. Plot on new figures as in 5 and 6 but give variables
% a different name so that you can also compare the cross-covariance and
% the cross-coherence on a final summary figure



