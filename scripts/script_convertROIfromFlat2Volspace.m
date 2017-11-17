%% need to
% create large ROI for pRF analysis in flatspace
% project through depths
% convert roi back to vol space
% export overlay to from vol space to flat space
% use PAC rois defined in flatmap space for ROI analysis


% set thisView to correct base and analysis
roi = viewGet(thisView,'roi','RightAC');
outputRoi = convertFromFlatVolumeToBase_depth6(roi)
outputRoi.name = [outputRoi.name 'Vol'];
thisView = viewSet(thisView,'newROI',outputRoi);