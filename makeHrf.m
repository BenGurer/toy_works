
function hrf = makeHrf(TR)
% given the TR, return the HRF shape for t = 0 ... 30s
%
% using the equation given in the lecture (simple boynton version)

tau =1.5; % time constant (s)
delta = 3; % time shift  (s)
t = [0:TR:30]; % vector of time points (in steps of TR)
% t = [0:1:30]; % vector of time points (in steps of seconds)

tshift = max(t-delta,0); % shifted, but not < 0
% hrf = (tshift ./ tau) .^2 .* exp(-tshift ./tau )./(2*tau);

hrf = thisGamma(t,1,delta,0,tau,3);
figure; plot(hrf);
hold on
x = 4;
y = 11;
z = 4;
hrf = gampdf(t,x,1)-gampdf(t,y,1)/z;
plot(hrf);

function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);