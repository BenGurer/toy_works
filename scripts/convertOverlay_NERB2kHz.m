function [ thisView , kHzdata ] = convertOverlay_NERB2kHz(thisView,overlayIN,overlayname)

% convert from NERB to kHzaa
kHzdata = InvNErb(overlayB.data{1});

scanDims = size(overlayIN.data{1});
% create overlay structure
overlaykHz = overlayIN;
dateString = datestr(now);
overlaykHz.date = dateString;
overlaykHz.name = overlayname;

overlaykHz.range = [0 max(max(max(overlayIN.data{1})))];
overlaykHz.clip = [0 max(max(max(overlayIN.data{1})))];
overlaykHz.colorRange = [0 max(max(max(overlayIN.data{1})))];

% difference.colorRange = [0 1];
% difference.range = [0 1];
% % difference.clip = [0 1];
% difference.range = [min(min(min(differenceData))) max(max(max(differenceData)))];
% difference.clip = [min(min(min(differenceData))) max(max(max(differenceData)))];
% difference.colorRange = [min(min(min(differenceData))) max(max(max(differenceData)))];
if exist('brewermap.m', 'file')
    overlaykHz.colormap = brewermap(256,'Reds');
else
    overlaykHz.colormap = jet(256);
end
% difference.alpha = 1;
% difference.colormapType = 'setRangeToMax';
overlaykHz.data{1} = nan(scanDims);
overlaykHz.data{1} = kHzdata;
% difference.data{scanNum}(x(i),y(i),z(i)) = thisr2(i);

thisView = viewSet(thisView,'newoverlay',overlaykHz);
end

function f = InvNErb(nerb)
% ***** lcfInvNErb *****
% Converts ERB number to frequency;
A = 24.7/1000; B = 4.37;
f = 1/B*(exp(A*B*nerb)-1);
end