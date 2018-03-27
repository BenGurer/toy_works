
function [params hrf t] = cal_DoubleGammaHrfModel(params,threshold,sampleDuration,sampleDelay)

% Double gamma HRF model used by GLM V2 plugin

tmax = max(params.y*3, 20); %min length of the hrf model in seconds

if isfield(params, 'tmax')
    tmax = params.tmax;
end

shift = 0;
if isfield(params, 'shift')
    shift = params.shift;
end

dt = 0.05;

t = 0:dt:tmax;
warning('off', 'MATLAB:log:logOfZero');
modelHrf = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
warning('on', 'MATLAB:log:logOfZero');

if shift<0
  modelHrf = [zeros(1, ceil(-shift/dt)), modelHrf];
elseif shift>0
  modelHrf = modelHrf( ceil(shift/dt):end );
end


% if params.includeDerivative
%   % take the derivative
%   modelHrfDerivative = [diff(modelHrf), 0];
%   % orthogonalize
%   modelHrfDerivative = modelHrfDerivative - modelHrf*(modelHrfDerivative/modelHrf);
%   % remove mean
%   modelHrfDerivative = modelHrfDerivative - mean(modelHrfDerivative);
%   % normalize so that its norm equals the Hrf norm
%   modelHrfDerivative = modelHrfDerivative / norm(modelHrfDerivative)*norm(modelHrf);
%   %concatenate
%   modelHrf = [modelHrf; modelHrfDerivative];
% end

%normalise so that integral of sum = 1
modelHrf = modelHrf./sum(modelHrf(:));
    
%downsample with constant integral
dsFactor = round(sampleDuration/dt);
dsDelay = floor(rem(sampleDelay,sampleDuration)/dt)+1;
hrf = mrDownsample(modelHrf', dsFactor, dsDelay);

%remove trailing zeros
%hrf = hrf(1:end-find(flipud(max(abs(hrf),[],2))>threshold,1,'first')+1,:);

%output the max amplitude of the actual model HRF
params.maxModelHrf = sampleDuration/dt * max(modelHrf'); 

% return actual time
t = t(dsDelay:dsFactor:length(t));
% make sure t is same length as hrf
if length(t) < length(hrf)
  t = [t nan(1,length(hrf)-length(t))];
elseif length(t) > length(hrf)
  t = t(1:length(hrf));
end