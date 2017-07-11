%% Generate a two dimensional stimulus (e.g two pixels) with correlations and 100 samples (e.g. points in time)
% First pixel data
nsamp = 100;
x1 = randn(1,nsamp);
% Second pixel that is correlated with the first
x2 = .4*x1 + .6*randn(1,nsamp);
% Concatinate into a stimulus matrix
x = [x1; x2];

%% Generate a filter and the corresponding one dimensional response 
% Set weights on each channel
h = [5; 7];
% Make response of system
y = x'*h;


%% set up a range of potential values for each weight (you are solving 
% the equation and are looking for possible correct answers)
[h1,h2] = meshgrid(-1:.2:10, -1:.2:10);
hs = [h1(:)'; h2(:)'];
% get responses from each set of weights
ys = x'*hs;
% calculate error between the response, y, and each of the possible responses, ys.  
err = sum((repmat(y,1,size(ys,2)) - ys).^2,1);
% reshape for plotting
err = reshape(err,size(h1));

%% plot contour of error surface. Note the shape of the surface is angled
% because the two variable are correlated.
figure(1);
contour(h1,h2,err,50);

%% Problem 1.
% Plot the actual solution as a large cross on the contour plot
% What is the value of the error surface at this solution?




%% Problem 2. 
% Generate a new response yr (for y real) that includes noise (SNR ~1). 
% Repeat some of the code above to obtain
% a new error surface as you would with real data.  Plot
% this on figure 2 and compare to figure 1. Calculate the minimum of this
% error and compare to the minimum of the noise free error surface.
% What happens to the error surface if you only have 10 data points (and noise)? 
% Calculate that one and plot in on a third plot.

% The new response:
SNR = 1;
yr = y + (std(y)/sqrt(SNR))*randn(1,nsamp)';

% calculate error between the response and each of the possible responses  


% reshape for plotting


% plot the contour of error surface.


% With fewer data points 
nsampsmall = 10;



%% Problem 3.
% Solve for the Least MSE solution using the analytical solution (calculate the cross and auto correlation and take ratio).
% You migh also check out what the \ does in matlab with the raw x and y
% vectors.
% Plot this solution as well as the correct solutions on your contour plots of figures 2 and 3.  

% your code here

% plot solution on contour
figure(2);
hold on; 

% your code here

hold off;

%% Problem 4. Now we're going to solve using gradient descent to show that you get the same answer as the MSE
%  You are going to solve it using the x and yr data. Use the
%  analytical solution for the gradient.  Fill in the missing code and
%  print out the final solution. Stop the descent when error is less than
%  the noise power.

% set a step size and a fixed number of maximum steps
nsteps = 1000;
vary = var(yr);
hscale = 10;    % This is a guess on the variance of the h parameters
stepsize = (hscale./vary); % This is to get stepsizes with the correct units.

% initialize hhat at origin and allocate space
hhat = zeros(2, nsteps);   %  We will be keeeping track of hhat during our descent

% loop for a certain number of iterations
% and mark when your reach the noise level...
totstep = -1;
for ii=1:nsteps
    
    % calculate the gradient at hhat
    
    
    % update hhat using stepsize and gradient
   
    
    % calculate error
      %err = ...
    
    % set stopping condition when error is below noise power.
    if err < vary/(SNR+1)
       if (totstep == -1)
           totstep = ii;
       end
    end
    
end
if (totstep == -1 ) 
    totstep = ii;
end
fprintf(1, 'Gradient solution of h = [%.2f, %.2f] reached after %d steps\n', hhat(1,totstep), hhat(2,totstep), totstep);


%% Problem 5.  Plot descent path on contour of figure 2.
% Comment on why early stopping could give you a more general solution.  

figure(2);
hold on; 
% your code here...
hold off;

%% Problem 6. 
% Repeat the calculations of Problem 4 and 5 with figure 3. 





%% Problem 6. Solve using coordinate descent
% Fill the code to solve the problem using coordinate descend and plot the 
% trajectories on figure 2 and figure 3 the path with a different
% color line.

nsteps = 1000;
% initialize hhatcd to origin and make space so that we can keep track
% during the descent
hhatcd = zeros(2, nsteps);   %  We will be keeeping track of hhat during our descent

% loop for a certain number of iterations
totstepcd = -1;
for ii=1:nsteps
    
    % calculate the gradient
  
    
    % find the maximum of the gradient
  
    
    % update hhat so that only the dimension with largest gradient changes

    
    % calculate error

    
    % set stopping criteria

    
    % end loop
end


% plot coordinate descent path on contour
figure(2);
hold on; 
% Your code here...
hold off;

% Repeat for the smaller data shown in figure 3. 
% initialize hhatcdsmall to origin and make space so that we can keep track
% during the descent


%% Problem 7.  Show on figures 2 and 3, the path corresponds to ridge regression.
% Remember that you calculated the auto and cross correlation in problem 3.




% hyper parameter lambda as fraction of largest eigenvalue

l=[10 1 0.5 0.1 0.05 0.01 .005 .001 .0005 .0001 .00001 0];
nlambda = length(l);

hr = zeros(2, nlambda);

% Loop over different ridge parameters
for ii=1:nlambda
    
   
    
    
end

% plot ridge path on contour


% Repeat for figure 3

