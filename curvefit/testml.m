% testml.m
%
% test out a few things with MarqLev curve fitter

% try out Gaussian noise
%x = (1:200)'; % x axis
%y = exp(-(x-100).^2/(2*20^2)) + randn([length(x),1])/10;

%plot(x, y);

fitdata = load gausswave.txt;
x = fitdata(:,1);
y = fitdata(:,2);
sigma = ones(size(x));
sz = size(fitdata);
if ( sz(2) == 3 ), sigma = fitdata(:,3); end

function yout = testmlfitfunc(x, params)
  A = params(1);
  s = params(2);
  x0 = params(3);
  yout = A*exp(-(x-x0).^2/(2*s^2));
end % function 


A = 0.8;
s = 10;
x0 = 90;
initparams = [A;s;x0];
%plot(x,[y,testmlfitfunc(x,initparams)]);


[params, errors] = MarqLev('testmlfitfunc',initparams,x,y,sigma)

plot(x,[y,testmlfitfunc(x,params)]);

