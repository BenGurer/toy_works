% roiTuningWidth.m            
%
%        $Id: dummy.m 2733 2013-05-13 11:47:54Z julien $
%      usage: [  ] = dummyInterrogator(thisView,overlayNum,scanNum,x,y,z,roi)
%         by: julien besle
%       date: 2010-03-09
%     inputs: 
%    outputs: 
%
%    purpose: computes tuning width after re-centering tuning curves on best frequency
%                           and averaging within an ROI. 
%                           plots the averaged tuning width and best gaussian fit
%                           and reports tuning width estimated as a weighted average

function [roiTuningWidth,roiTuningCurve] = roiTuningWidth(thisView,overlayNum,scanNum,x,y,z,roi)
maskingOn = false;
nConditions = 15;

if ieNotDefined('roi')
  roiNums = viewGet(thisView,'visibleROIs');
  cRoi=0;
  for iRoi = roiNums
    cRoi=cRoi+1;
    roi{cRoi}=viewGet(thisView,'roi',iRoi);
  end
end
if ieNotDefined('scanNum')
    scanNum = viewGet(thisView,'curscan');
end


%get overlay data and mask
overlayNum=[2:1+nConditions];
alphaOverlays = zeros(1,length(overlayNum));
for iOverlay = 1:length(overlayNum)
  thisNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay',overlayNum(iOverlay)));
  if ~isempty(thisNum)
    alphaOverlays(iOverlay) = thisNum;
  end
end

[mask,overlayData] = maskOverlay(thisView,overlayNum,scanNum);
overlayData = overlayData{1};

% if maskingOn
%     maskAlpha = maskOverlay(thisView,alphaOverlays,scanNum);
%     
%     mask = mask{1};
%     maskAlpha = maskAlpha{1};
%     %mask the overlay
%     overlayData(~mask) = NaN;
%     %mask with the alpha overlay
%     for iOverlay = 1:length(overlayNum)
%         if alphaOverlays(iOverlay)
%             thisOverlayData = overlayData(:,:,:,iOverlay);
%             thisAlphaMask = maskAlpha(:,:,:,iOverlay);
%             thisOverlayData(~thisAlphaMask)=NaN;
%             overlayData(:,:,:,iOverlay) = thisOverlayData;
%         end
%     end
% end

roiTuningCurve= NaN(13,length(roi));
for iRoi=1:length(roi)
  %get ROI scan coords
  scanCoords = getROICoordinates(thisView,roi{iRoi});
  scanDims = viewGet(thisView,'scandims');
  scanCoordsIndex = sub2ind(scanDims,scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
  nAllVoxels = length(scanCoordsIndex);

  %get the current overlay data
  clear curve
  for i=1:nConditions
%     [mask,overlayData] = maskOverlay(thisView,i+1,scanNum);
%     mask = mask{1};
%     overlayData = overlayData{1};
%     mask the overlay
%     overlayData(~mask) = NaN;
%     see if there`s an alpha overlay and mask with it
%     thisNum = viewGet(thisView,'overlayNum',viewGet(thisView,'alphaOverlay',i+1));
%     if ~isempty(thisNum)
%       maskAlpha = maskOverlay(thisView,thisNum,scanNum);
%       maskAlpha = maskAlpha{1};
%       overlayData(~maskAlpha) = NaN;
%     end
    thisOverlayData = overlayData(:,:,:,i);
    curve(i,:) = thisOverlayData(scanCoordsIndex);
  end

  %remove empty voxels
  curve(:,any(isnan(curve)))=[];
  nVoxels = size(curve,2);
    stdDevCons = mean(std(curve))
    stdDevVox = mean(std(curve,0,2))

  %normalize data by max
  % curve = curve./repmat(max(curve),nConditions,1);

  % %first, set negative values to 0 (?)
  % curve(curve<0)=0;

  %find best frequency using Humphries method
  % remove  values that are less than the average activation across conditions or less than 0
  averageActivation = repmat(mean(curve),[nConditions 1]);
  curveNaN=curve;
  curveNaN(curve<averageActivation | curve<0)=NaN;
  %compute the weighted average
  indices = repmat((1:nConditions)',[1 nVoxels]);
  average = nansum( curveNaN .* indices) ./ nansum(curveNaN);

  resolutions = [1 10];
  nResolutions = length(resolutions);
  figure('name',sprintf('%s: %d/%d voxels',roi{1}.name,nVoxels,nAllVoxels));
  hiresIndices = 1-nConditions:.01:nConditions-1;
  options = optimset('Display','off');
  for i = 1:nResolutions
  %   switch(i) %loop over resolution
  %     case nResolutions
  %       [~,bestFrequency] = max(curve);
  %       resString = 'Max Frequency';
  %     otherwise
        bestFrequency = round(average*resolutions(i))/resolutions(i);
        resString = {'Weighted average',sprintf('Bins=%.3f',1/resolutions(i))};
  %   end

    % recentre tuning curves using estimated best frequency
    indices2 = indices - repmat(bestFrequency,nConditions,1);

    uniqueIndices = unique(indices2);
    uniqueIndices = uniqueIndices(~isnan(uniqueIndices));
    tuningCurve = zeros(length(uniqueIndices),1);
    c=0;
    for j=uniqueIndices'
      c = c+1;
      tuningCurve(c) = mean(curve(indices2 == j));
  %     numberPoints(c) = nnz(indices2 == j);
    end

    nPreProcess=1;
    for j=1:nPreProcess %loop over possible pre-processing
      %first, apply pre-processing
      thisTuningCurve = tuningCurve;
      switch(j)
        case 1
          %no pre-processing
          preProcessString='No pre-processing';
        case 2
          %set negative values to 0
          thisTuningCurve(thisTuningCurve<0)=0;
          preProcessString='Set negative values to 0';
        case 3
          % normalize curve
          thisTuningCurve = (thisTuningCurve-min(thisTuningCurve))/(max(thisTuningCurve)-min(thisTuningCurve));
          preProcessString='Normalize (min/max)';
        case 4
          % normalize curve only if negative value
          if any(thisTuningCurve<0)
            thisTuningCurve = (thisTuningCurve-min(thisTuningCurve))/(max(thisTuningCurve)-min(thisTuningCurve));
          end
          preProcessString='Normalize if negative values';
      end
      %compute weighted mean and standard-deviations
      weightedMean = sum(uniqueIndices.*thisTuningCurve)/sum(thisTuningCurve);
      weightedStdDev = sqrt( sum( thisTuningCurve .* (uniqueIndices - weightedMean ).^2) ./ sum(thisTuningCurve) );

      %fit gaussian
      initialParams=[0 1 1 0];
  %     [params,~,~,exitflag] = lsqcurvefit(@gaussFunction,[0 1 1],double(uniqueIndices),double(thisTuningCurve),[],[],options);
  %     [params,~,~,exitflag] = lsqcurvefit(@gaussFunction2,[0 1 1 0],double(uniqueIndices),double(thisTuningCurve),[-inf 0 0 -inf],[inf inf inf 0],options);
      if min(0,min(thisTuningCurve)) %if minimum value is negative, fit constant parameter with boundaries [min 0]
        [params,~,~,exitflag] = lsqcurvefit(@gaussFunction3,initialParams,double(uniqueIndices),double(thisTuningCurve),[-inf 0 0 min(0,min(thisTuningCurve))],[inf inf inf 0],options);
      else %if minimum value is positive, don't fit constant parameter (set it to 0)
        [params,~,~,exitflag] = lsqcurvefit(@gaussFunction4,initialParams,double(uniqueIndices),double(thisTuningCurve),[],[],options);
      end
      roiTuningWidth(1,iRoi)=params(2);
      subplot(nResolutions,nPreProcess,(i-1)*nPreProcess+j);
      plot(uniqueIndices,thisTuningCurve,'o');
      roiTuningCurve(1:length(thisTuningCurve),iRoi)=thisTuningCurve;
      hold on
  %     plot(hiresIndices,gaussFunction(params,hiresIndices),'g');
      plot(hiresIndices,gaussFunction(params,hiresIndices),'g');
      plot([1-nConditions nConditions-1],[0 0],'--k');
  %     title(sprintf('Weighted std dev = %.2f, fitted std dev = %.2f', weightedStdDev, params(2)));
      title(sprintf('Mean = %.2f, Std-dev = %.2f, scaling = %.2f, constant= %.2f', params(1), params(2), params(3), params(4)));
      xlim([1-nConditions nConditions-1]);
  %     ylim([-0.2 1]);
      if j==1
        ylabel(resString);
      end
      if i==nResolutions
        xlabel(preProcessString);
      end
    end
  end
end

roiTuningWidth = roiTuningWidth';
roiTuningCurve = roiTuningCurve(:);

return;


function kernel = gaussianKernel(FWHM)

sigma_d = FWHM/2.35482;
w = ceil(FWHM); %deals with resolutions that are not integer
%make the gaussian kernel large enough for FWHM
if length(w)==1
  w = [w w w];
  sigma_d = [sigma_d sigma_d sigma_d];
end
kernelDims = 2*w+1;
kernelCenter = ceil(kernelDims/2);
[X,Y,Z] = meshgrid(1:kernelDims(1),1:kernelDims(2),1:kernelDims(3));
kernel = exp(-((X-kernelCenter(1)).^2/(2*sigma_d(1)^2)+(Y-kernelCenter(2)).^2/(2*sigma_d(2)^2)+(Z-kernelCenter(3)).^2/(2*sigma_d(3)^2))); %Gaussian function
kernel = kernel./sum(kernel(:));

function fValues = gaussFunction(params,values)
%params(1) is the mean
%params(2) is the std-deviation
%params(3) is a scaling  factor
%params(4) is a constant term
fValues = params(4) + params(3) * exp(-(values - params(1)).^2/2/params(2)^2);

function fValues = gaussFunction2(params,values)
%params(1) is the mean
%params(2) is the std-deviation
%params(3) is a scaling factor
%params(4) is a constant term and is set to 0
fValues = params(3) * exp(-(values - params(1)).^2/2/params(2)^2);

function fValues = gaussFunction3(params,values)
%params(1) is the mean and is set to 0
%params(2) is the std-deviation
%params(3) is a scaling  factor
%params(4) is a constant term
fValues = params(4) + params(3) * exp(-(values).^2/2/params(2)^2);

function fValues = gaussFunction4(params,values)
%params(1) is the mean and is set to 0
%params(2) is the std-deviation
%params(3) is a scaling  factor
%params(4) is a constant term and is set to 0
fValues = params(3) * exp(-(values).^2/2/params(2)^2);

function fValues = weightedGaussFunction(params,values,weights)
%params(1) is the mean
%params(2) is the std-deviation
%params(3) is a scaling  factor
%params(4) is a constant term
fValues = params(4) + params(3) * exp(-(values - params(1)).^2/2/params(2)^2) .* weights;

function fValues = weightedGaussFunction2(params,values,weights)
%params(1) is the mean
%params(2) is the std-deviation
%params(3) is a scaling  factor
%params(4) is a constant term and is set to 0
fValues = params(3) * exp(-(values - params(1)).^2/2/params(2)^2) .* weights;

function fValues = weightedGaussFunction3(params,values,weights)
%params(1) is the mean and is set to 0
%params(2) is the std-deviation
%params(3) is a scaling  factor
%params(4) is a constant term
fValues = params(4) + params(3) * exp(-(values).^2/2/params(2)^2) .* weights;

function fValues = weightedGaussFunction4(params,values,weights)
%params(1) is the mean and is set to 0
%params(2) is the std-deviation
%params(3) is a scaling  factor
%params(4) is a constant term and is set to 0
fValues = params(3) * exp(-(values).^2/2/params(2)^2) .* weights;