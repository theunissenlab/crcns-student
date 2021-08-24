
function spiketimes = poisson_gen_spikes(psth, meanfr, numTrials)
% Generates an array of spike arrival times at a mean firing rate given by
% meanfr and a time varying rate given by psth in ms. 
% Written by Alex Huth 11/12/2009
% modified by Mike Schachter
%generates spike times in units of seconds
inbinwidth = 1; %% ms

nt = length(psth);
tmax = nt * inbinwidth; %% ms
times = 1:tmax;
spiketimes = cell(1, numTrials);

%% First generate enough poisson points (with lamstar) to cover our interval
lambda = psth;
lambda = lambda - min(lambda(:)); %% Scale to be min. zero
lambdaMean = mean(lambda(:));
lambda = lambda / lambdaMean * meanfr; %% Scale to be mean 'meanfr'
lamstar = max(lambda(:)) * 1.2; %% Find the maximum value, go a little above it

lamgen = lamstar * tmax / 1000;

for ci = 1:numTrials

   nspikes = round(randn * sqrt(lamgen) + lamgen);
   sptime = floor(rand(1,nspikes) * (tmax-1))+1;
   fastspikes = times(sptime);

   %% Remove each spike with probability lambda / lamstar
   %% This is the 'thinning' method first published by Lewis & Shedler (1979)
   rands = rand(1,numel(fastspikes));
   laminterp = interp1( (1:nt)*inbinwidth, lambda, fastspikes, 'linear');
   throw_out = rands < (1 - (laminterp / lamstar));

   spiketimes{ci} = sort(fastspikes(~throw_out));
end
  