function [thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,roiName,weightStim)
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
if ~weightStim
    for iGroup = 1:length(Info.groupNames)
        thisView = viewSet(thisView,'curGroup',glmInfo.groupNames);
        analysisSaveName = ['pRF' '_' roiName];
        [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
        pRFParams.saveName = [analysisSaveName];
        pRFParams.restrict = ['ROI: ' roiName];
        pRFParams.pRFFit.supersampling = 1;
        pRFParams.pRFFit.fitHDR = 0;
        [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
    end
    
else
    
    %     pRFinfo.analysisNames_Groups{1}{1} = {'pRF'};
    % pRFInfo.stimulusWeighting{1} = {'None'};
    % pRFInfo.stimulusWeighting{2} = {'None','SL_level','BOLD','fit'};
    % for iWeight = 1:length(stimulusWeighting)
    %     pRFinfo.analysisNames_Groups{2}{iWeight} = {['pRF_' stimulusWeighting{2}{iWeight}]};
    % end
    
    
    for iGroup = 1: length(glmInfo.groupNames)
        thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{iGroup});
        stimulusWeighting = pRFInfo.stimulusWeighting{iGroup};
        for iWeight = 1:length(stimulusWeighting)
            analysisSaveName = [ pRFInfo.analysisNames_Groups{iGroup}{iWeight} '_' roiName];
            [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
            pRFParams.saveName = analysisSaveName;
            pRFParams.restrict = ['ROI: ' roiName];
            pRFParams.pRFFit.supersampling = 1;
            pRFParams.pRFFit.fitHDR = 0;
            
            pRFParams.pRFFit.stimulusWeighting = stimulusWeighting{iWeight}; % {'None','SL_level','BOLD','fit'}
            if strcmpi(stimulusWeighting{iWeight},'BOLD')
                pRFParams.pRFFit.SWgradient = glmInfo.m;
                pRFParams.pRFFit.SWoffset = glmInfo.b;
            end
            [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
            
        end
        
    end
    
    
    %     thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{1});
    %     analysisSaveName = 'pRF';
    %     [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
    %     pRFParams.saveName = [analysisSaveName];
    %     pRFParams.restrict = ['ROI: ' roiName];
    %     pRFParams.pRFFit.supersampling = 1;
    %     pRFParams.pRFFit.fitHDR = 0;
    %     [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
    %
    %     thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{2});
    % %     stimulusWeighting = {'None','SL_level','BOLD','fit'};
    %     stimulusWeighting = pRFInfo.stimulusWeighting
    %     for iWeight = 1:length(stimulusWeighting)
    % %         analysisSaveName = 'pRF';
    %         [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1');
    %         pRFParams.saveName = [analysisSaveName '_' stimulusWeighting{iWeight}];
    %         pRFParams.restrict = ['ROI: ' roiName];
    %         pRFParams.pRFFit.supersampling = 1;
    %         pRFParams.pRFFit.fitHDR = 0;
    %
    %         pRFParams.pRFFit.stimulusWeighting = stimulusWeighting{iWeight}; % {'None','SL_level','BOLD','fit'}
    %         if strcmpi(stimulusWeighting{iWeight},'BOLD')
    %             pRFParams.pRFFit.SWgradient = glmInfo.m;
    %             pRFParams.pRFFit.SWoffset = glmInfo.b;
    %         end
    %         [thisView, pRFParams] = pRF_auditory(thisView,pRFParams);
    %
end

% end

%% individual scans
% pRFParams = [];
% thisView = viewSet(thisView,'curGroup','MotionComp');
% for iScan = 1:glmInfo.nScans
%     thisView = viewSet(thisView,'curScan', iScan);
%     [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
%     pRFParams.saveName = ['pRF_Scan_' mat2str(iScan)];
%     pRFParams.restrict = ['ROI: ' roiName];
%     pRFParams.pRFFit.supersampling = 1;
%     % pRFParams.pRFFit.dispStimScan = iScan;
%     pRFParams.pRFFit.fitHDR = 0;
%     [thisView, pRFParams] = pRF_auditory(thisView,pRFParams,['scanList=' mat2str(iScan)]);
% end

end