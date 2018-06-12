function [ thisView , ERBdata ] = convertOverlay_GLMjpCF2NERB(thisView,overlayIN,stimInfo)

% convert from stimulus-spacing-space to NERB
% offset all values by the lowests frequencies NERB - this accounts for the
% stimulus spacing
% minus the difference stimuli in NERB to account for be-biasing step -
% this step calucated from 0 to (number of overlays .*1.25)
ERBdata = overlayB.data{1} + (stimInfo.stimNERBs(1)-(stimInfo.stimNERBs(2)-stimInfo.stimNERBs(1)));

scanDims = size(overlayIN.data{1});
% create overlay structure
overlayERB = overlayIN;
dateString = datestr(now);
overlayERB.date = dateString;
overlayERB.name = 'pCF_ERB';

overlayERB.range = [0 max(max(max(overlayIN.data{1})))];
overlayERB.clip = [0 max(max(max(overlayIN.data{1})))];
overlayERB.colorRange = [0 max(max(max(overlayIN.data{1})))];

% difference.colorRange = [0 1];
% difference.range = [0 1];
% % difference.clip = [0 1];
% difference.range = [min(min(min(differenceData))) max(max(max(differenceData)))];
% difference.clip = [min(min(min(differenceData))) max(max(max(differenceData)))];
% difference.colorRange = [min(min(min(differenceData))) max(max(max(differenceData)))];
if exist('brewermap.m', 'file')
    overlayERB.colormap = brewermap(256,'Reds');
else
    overlayERB.colormap = jet(256);
end
% difference.alpha = 1;
% difference.colormapType = 'setRangeToMax';
overlayERB.data{1} = nan(scanDims);
overlayERB.data{1} = ERBdata;

thisView = viewSet(thisView,'newoverlay',overlayERB);

end