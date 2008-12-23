function [params, errors] = MarqLev(fitfcn, initparams, xdata, ydata, sigma)

% MarqLev  uses Marquardt-Levenberg routine to fit an arbitrary function
%   to a given set of data.
%
%   [params, errors] = MarqLev(fitfcn, initparams, xdata, ydata, sigma)
%
%   performs the fit, where
%
%    params -- is a column vector of best fit parameters
%    errors -- is a column vector of the estimated errors associated with
%                 each parameter (expressed as a fraction of the parameter
%                 value)
%    fitfcn -- is a string containing the name of the function, which has
%                 format  fitfcn(xdata, params)
%    initparams -- is a column vector specifying the initial guesses for
%                 the parameter values
%    xdata -- is a column vector with the data x values
%    ydata -- is a column vector with the data y values
%    sigma -- (optional) is a column vector with the standard deviations
%                 for the y data values
%
% Final chi squared statistic is returned in the global variable
%   gChiSquared.
%
% Based on discussion in "Data Reduction and Error Analysis for the 
%    Physical Sciences" by Philip R. Bevington.
 
global gChiSquared

if nargin == 4 
  sigma = ones(length(ydata), 1);
end

s = size(ydata);
if s(1) == 1
  ydata = ydata.';
end
s = size(xdata);
if s(1) == 1
  xdata = xdata.';
end
s = size(initparams);
if s(1) == 1
  initparams = initparams.';
end
s = size(sigma);
if s(1) == 1 
  sigma = sigma.';
end


lambda = 0.001;
stopflag = 0;
params = initparams;
iterations = 0;
recalculateflag = 1;
maxiterations = 100;

% calculate chi squared at initial parameters
ycalc = feval(fitfcn, xdata, params);
chisquared = ((ydata - ycalc).^2).' * (1./(sigma.^2))
    
     
for i=1:length(params)
  expandedsigmasquared(:, i) = sigma.^2;
end

% main iteration loop

while (stopflag < 2)

  %  calculate partial derivatives
  
  if (recalculateflag == 1)
    for i=1:length(params)
      delta = zeros(length(params), 1);
      delta(i) = max([.0001, abs(.01*params(i))]);
      derivs(:, i) = (feval(fitfcn, xdata, params + delta) - ...
         feval(fitfcn, xdata, params - delta))/2/sum(delta);
    end
  
  %  calculate alpha matrix
    alpha = (derivs./expandedsigmasquared).'*derivs;
  end
  alphaprime = alpha + lambda * diag(diag(alpha));
  
  % calculate beta vector
  
  beta = derivs.' * ((ydata - feval(fitfcn, xdata, params))./(sigma.^2));

  % calculate change to params
  
  newparams = params + alphaprime\beta
  
  % calculate new chi squared
  
  newchisquared = ((ydata - feval(fitfcn, xdata, newparams)).^2).' * ...
     (1./(sigma.^2));
  
  % check for improvement
  
  deltachisquared = newchisquared - chisquared;
  if (deltachisquared > 0) | (isnan(newchisquared))
    lambda = lambda * 10;
    recalculateflag = 0;
  else
    lambda = lambda/10;
    recalculateflag = 1;
    iterations = iterations + 1;
    if (deltachisquared > -.01) & (-deltachisquared < .001*chisquared)
      stopflag = stopflag + 1;
    else
      stopflag = 0;
    end
    
    if (iterations > maxiterations)
      stopflag = 3;
    end
    
    chisquared = newchisquared;
    params = newparams;
  end
  
  iterations, params, newchisquared
  
end


% Calculate covariance matrix and then errors

  for i=1:length(params)
    delta = zeros(length(params), 1);
    delta(i) = max([.01, abs(.01*params(i))]);
    derivs(:, i) = (feval(fitfcn, xdata, params + delta) - ...
      feval(fitfcn, xdata, params - delta))/2/sum(delta);
  end
    
  alpha = (derivs./expandedsigmasquared).'*derivs;
  covariance = inv(alpha);
  errors = sqrt(diag(covariance))./params;
 
  gChiSquared = newchisquared;