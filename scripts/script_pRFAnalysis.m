function [thisView, pRFdata] = script_pRFAnalysis(thisView,pRFInfo,iScan)
%
%   usage: script_pRFAnalysis(thisView,pRFInfo)
%      by: Ben Gurer
%    date: 22/11/2017
% purpose: script pRF analysis of fMRI data
%   input: mrView, pRF analysis information
%  output: updated mrView, pRF data
%
% NEED: ROI, save name, fit hdr?, weight stim?
% this function will be run after glm, gradient reversals and roi creation

analysisSaveName = 'pRF - ALL CONS - Scan ';
% for iScan = 1:nScans
% thisView = viewSet(thisView,'curGroup','MotionComp');
% thisView = viewSet(thisView,'curScan', iScan);
if ~isempty(iScan)
    [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
else
    [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
end
pRFParams.saveName = [analysisSaveName mat2str(iScan)];
pRFParams.restrict = ['ROI: ' roi{1, roiNum}.name];
pRFParams.pRFFit.supersampling = 1;
pRFParams.pRFFit.fitHDR = 0;
pRFParams.pRFFit.dispStimScan = iScan;

pRFParams.pRFFit.stimulusWeighting = 'None'; % {'None','SL_level','BOLD','fit'}
pRFParams.pRFFit.SWgradient = 1;
pRFParams.pRFFit.SWoffset = 0;

% [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1',['scanList=' mat2str(iScan)]);
[thisView, pRFParams] = pRF_auditory(thisView,pRFParams,['scanList=' mat2str(iScan)]);
% save analysis
saveAnalysis(thisView,[analysisSaveName mat2str(iScan)])
% analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',[analysisSaveName mat2str(iScan)]));
% end

end