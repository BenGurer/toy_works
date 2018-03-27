function [ x_doubleGamma, x_Gamma, x_dGamma ] = cal_hrfAverageROI(e,analysisParams)
    %% Centre max response to 0 and average
%     max(a,[],2)
% take mean of all and find max time point and use that
%     hdrMaxTimePoint = 3; % [v i] = max(mean(estimate.hdr(:,:,1),2))
    [v hdrMaxTimePoint] = max(max(mean(estimate.hdr,3)'))
    curve = e.hdr(hdrMaxTimePoint,:,:);
    nOverlays = analysisParams.nhdr;
    resolution = 1; %resolution = 2;
    overlaySize = size(e.hdr);
    overlaySize = nVoxels;
    %find best frequency using Humphries method
    % remove  values that are less than the average activation across conditions or less than 0
    averageActivation = repmat(mean(curve,3),[1 1 nVoxels]);
    curveNaN=curve;
    curveNaN(curve<averageActivation | curve<0)=NaN;
    %compute the weighted average
    indices = permute(repmat(1:nOverlays,[nVoxels 1]),[3 2 1]);
    bestFrequency = round(resolution*nansum( curveNaN .* indices, 2) ./ nansum(curveNaN,2))/resolution;
    
    % recentre tuning curves using estimated best frequency
    bestFrequency = reshape(bestFrequency,prod(overlaySize),1);
    % curve = reshape(curve,prod(overlaySize),nOverlays);
    uniqueTuningCurveIndices = 1-nOverlays:1/resolution:nOverlays-1;
    tuningCurvesIndices=repmat(permute(uniqueTuningCurveIndices,[1 3 4 2]),[overlaySize 1]);
    tuningCurvesIndices=repmat(permute(uniqueTuningCurveIndices,[1 2]),[overlaySize 1]);
    nTuningCurvesIndices=length(uniqueTuningCurveIndices);
    tuningCurves = nan(1, nTuningCurvesIndices,prod(overlaySize));
    HDRS = nan(analysisParams.nHrfComponents,nTuningCurvesIndices,prod(overlaySize));
    uniqueBestFrequencies = unique(bestFrequency);
    uniqueBestFrequencies=uniqueBestFrequencies(~isnan(uniqueBestFrequencies))';
    for i=uniqueBestFrequencies
        bestFrequencyIndices=find(bestFrequency==i);
        %   tuningCurves(bestFrequencyIndices,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))])=curve(bestFrequencyIndices,:);
        tuningCurves(1,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))],bestFrequencyIndices)=curve(:,:,bestFrequencyIndices);
        HDRS(:,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))],bestFrequencyIndices)=e.hdr(:,:,bestFrequencyIndices);
    end
    % tuningCurves = reshape(tuningCurves,[overlaySize nTuningCurvesIndices]);
    unsmoothedTuningCurves = tuningCurves;
    averageTunindCurves = nansum(unsmoothedTuningCurves,3);
    figure; plot([1-nOverlays:1/resolution:nOverlays-1],averageTunindCurves);
    averageHDRTunindCurves = nansum(HDRS,3);
    figure;
    for i = 1:analysisParams.nHrfComponents
        % plot([1:1:glmData.nHrfComponents],[1-nOverlays:1/resolution:nOverlays-1],averageHDRTunindCurves(i,:));
        plot(averageHDRTunindCurves(i,:));
        legendInfo{i} = ['HDR Component ' num2str(i)];
        hold on
    end
    legend(legendInfo)
    
    figure; waterfall([1-nOverlays:1/resolution:nOverlays-1],[1:1:analysisParams.nHrfComponents],averageHDRTunindCurves);
    figure; surf([1-nOverlays:1/resolution:nOverlays-1],[1:1:analysisParams.nHrfComponents],averageHDRTunindCurves);
    a = permute(averageHDRTunindCurves, [2 1 3]);
    figure; surf([1:1:analysisParams.nHrfComponents],[1-nOverlays:1/resolution:nOverlays-1],a);
    
    ydata = a(7,:)/max(a(7,:)) % normalise to allow analysis to scale and offset
    xdata = e.time
    fun = @(x,xdata)x(1).*(((xdata-x(2))/x(4)).^(x(5)-1).*exp(-(xdata-x(2))/x(4)))./(x(4)*factorial(x(5)-1))+x(3); 
    [ result ] = fmribHRF(t, P, p2, p3, p4, p5, scaled)
    x0 = [1 4 3 8 1];
figure;    plot(xdata,fun(x0,xdata))
x = lsqcurvefit(fun,x0,xdata,ydata)


% 'timelag',1,'minmax=[0 inf]','incdec=[-0.5 0.5]','The timelag of the gamma function used to model the HDR. If using gaussian-hdr, this is just the initial value and the actual value will be fit.'};
% {'tau',0.6,'minmax=[0 inf]','incdec=[-0.1 0.1]','The tau (width) of the gamma function used to model the HDR. If using gaussian-hdr, this is just the initial value and the actual value will be fit.'};
% {'exponent',4,'minmax=[0 inf]','incdec=[-1 1]','The exponent of the gamma function used to model the HDR. This is always a fixed param.'};
% 'diffOfGamma',1,'type=checkbox','Set to true if you want the HDR to be a difference of gamma functions - i.e. have a positive and a delayed negative component'};
% paramsInfo{end+1} = {'amplitudeRatio',0.25,'minmax=[0 inf]','incdec=[-0.1 0.1]','Ratio of amplitude of 1st gamma to second gamma','contingent=diffOfGamma'};
% paramsInfo{end+1} = {'timelag2',2,'minmax=[0 inf]','incdec=[-0.5 0.5]','Time lag of 2nd gamma for when you are using a difference of gamma functions','contingent=diffOfGamma'};
% paramsInfo{end+1} = {'tau2',1.2,'minmax=[0 inf]','incdec=[-0.1 0.1]','The tau (width) of the second gamma function.','contingent=diffOfGamma'};
% paramsInfo{end+1} = {'exponent2',11,'minmax=[0 inf]','incdec=[-1 1]','The exponent of the 2nd gamma function.','contingent=diffOfGamma'};
% {'dispHDR',0,'type=pushbutton','buttonString=Display HDR','Display the HDR with the current parameters','callback',@pRFGUIDispHDR,'passParams=1'};


% xdata1 = voxelHRF.time'; %x axis
% ydata1 = voxelHRF.hdr; %y data
scaled = 1;
p_doubleGamma = [6 12 0.9 0.9 0.35 scaled]; %guess
opts = optimset('MaxFunEvals', 500, 'Display', 'off');
[x_doubleGamma, resnorm, ~, exitflag, output] = lsqcurvefit(@myHRF, p_doubleGamma, xdata, ydata, [], [], opts);
myHRF(p_doubleGamma, xdata)
figure
fittedHRF = myHRF(x_doubleGamma,xdata)
plot(xdata,fittedHRF)
hold on
% deconvHRF = a(7,:)
deconvHRF =ydata
plot(deconvHRF)


%% GLM Double Gamma

p_doubleGamma = [6 8 1 1 0]; %guess
lb_DoubleGamma = [0 0 0 0 0];
ub_DoubleGamma = [15 15 15 1 10];
opts = optimset('MaxFunEvals', 500, 'Display', 'off');
[x_doubleGamma, resnorm, ~, exitflag, output] = lsqcurvefit(@get_HRFDoubleGamma, p_doubleGamma, xdata, ydata, lb_DoubleGamma, ub_DoubleGamma, opts);
figure
disp(x_doubleGamma)
fittedHRF = get_fitDoubleGamma(x_doubleGamma,xdata);
plot(xdata,fittedHRF)
hold on
plot(xdata,get_fitDoubleGamma(p_doubleGamma,xdata), '--r')
% deconvHRF = a(7,:)
deconvHRF = ydata;
plot(xdata,deconvHRF)

% time = xdata;
% timelag = x(1);
% tau = x(2);
% exponent = x(3);
% amplitude = x(4);
% offset = x(5);
p_Gamma = [1 0.6 4 1 0]; %guess - Plot?
lb_Gamma = [0 0 0 0 0];
ub_Gamma = [16 inf 16 1 10];
opts = optimset('MaxFunEvals', 500, 'Display', 'off');
[x_Gamma, resnorm, ~, exitflag, output] = lsqcurvefit(@get_HRFGamma, p_Gamma, xdata, ydata, lb_Gamma, ub_Gamma, opts);
figure
disp(x_Gamma)
fittedHRF = get_fitGamma(x_Gamma,xdata);
plot(xdata,fittedHRF)
hold on
plot(xdata,get_fitGamma(p_Gamma,xdata), '--r')
% deconvHRF = a(7,:)
deconvHRF = ydata;
plot(xdata,deconvHRF)

% timelag = x(1);
% tau = x(2);
% exponent = x(3);
% timelag2 = x(4);
% tau2 = x(5);
% exponent2 = x(6);
x_dGamma = [];
p_dGamma = [1 0.6 4 2 1.2 11 0.25]; %guess - Plot?
lb_dGamma = [0 0 0 2 0 0 0];
ub_dGamma = [16 inf 16 16 inf 16 1];
opts = optimset('MaxFunEvals', 500, 'Display', 'off');
[x_dGamma, resnorm, ~, exitflag, output] = lsqcurvefit(@get_HRFDiffOfGamma, p_dGamma, xdata, ydata, lb_dGamma, ub_dGamma, opts);
figure
disp(x_dGamma)
fittedHRF = get_fitGamma(x_dGamma,xdata);
plot(xdata,fittedHRF)
hold on
plot(xdata,get_fitGamma(p_dGamma,xdata), '--r')
plot(xdata,ydata)
% make function to wrap around hrf functions - allow user choose hrf
% function - normalise hrf level and then scale.
% use fitted params for glm and prf analysis - need to use the same hrf
% model!
end
% function doublegammafun = thisDoubleGamma(x,xdata) 
% t = xdata;
% amplitude = x(4);
% offset = x(5);
% % modelHrf = gampdf(t, params.x, 1) - gampdf(t, params.y, 1)/params.z;
% modelHrf = gampdf(t, x(1), 1) - gampdf(t, x(2), 1)/x(3);
% % compare normalising intergral with normalising max
% 
% % create wrapper function to allow normalisation
% 
% %normalise so that integral of sum = 1
% modelHrf = modelHrf./sum(modelHrf(:));
% 
% doublegammafun= (amplitude*modelHrf+offset);
% % end
% function gammafun = thisGamma(x,xdata)
% time = xdata;
% amplitude = x(1);
% timelag = x(2);
% offset = x(3);
% tau = x(4);
% exponent = x(5);
% % exponent = round(exponent);
% % gamma function
% gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));
% 
% % negative values of time are set to zero,
% % so that the function always starts at zero
% gammafun(find((time-timelag) < 0)) = 0;
% 
% % normalize the amplitude
% if (max(gammafun)-min(gammafun))~=0
%     gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
% end
% gammafun = (amplitude*gammafun+offset);
% 
% end