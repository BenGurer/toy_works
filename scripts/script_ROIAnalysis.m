function data = script_ROIAnalysis(thisView,concatenationGroupNames,analysisNames_group,analysisNames_scans,stimData,conditionSplitIndex)
%% To do
% define analysis names when performing them and save them - input to this function only one you want to analyise _ then only loop over length and index it - use eval to use
% function to create figure based on paper size - poster or paper
% send figure handle to plotting functions - can be used with subplot and plot
% create strucutre for rois or containing rois

% following analysis - save roi results underone structure of subject ID
% plots I want
% ROI average beta weight
% ROI TW conA v B

%% SHOULD SET ANALYSIS NAME ONCE AND PASS IT ON
% loop - create analysis names
% loop length of analysis names
% needs to name based on group - one structure

%% Save data in structures labelled by roi
%% ROI names
% Create ROIs with these names:
% RightAC
% RightPosAC
% RightAntAc
% LeftAC
% LeftPosAC
% LeftAntAC
% AC
ROInames = {'RightAC','RightPosAC','RightAntAC','LeftAC','LeftPosAC','LeftAntAC','AC'};
ROInames = {'RightAC','RightPosAC','RightAntAC','LeftAC','LeftPosAC','LeftAntAC'};

q = char(39);
for iROI = 1:length(ROInames)
    eval([ROInames{iROI} ' = struct;']);
    eval([ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
end

% get data from ROIs
for iROI = 1:length(ROInames)
    for iScan = 1:nScans
        for iStim = 1:length(nStim)
            eval([ROInames{iROI} '.glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '{iScan} = getROIdata_GLM(thisView,GLManalysisData_nCons_' mat2str(nStim(iStim)) '{iScan},' ROInames{iROI} '.roi,iScan);']);
        end
    end
end

%% now pass data to function to plot
% loop over ROIs
% data = roiDataAnalysis(roiData,analName,conditionSplitIndex,restrictIndex,stimFreqs);
%%%% SPLIT UP FUNCTIONS
% new function
% 1 get data from scans
% new function
% get roi data - save in subject strucutre
% new function
% roi analysis

%% NEW FUNCTION %%
% move this outside of function
conditionSplitIndex = {[2,4],[1,3]};
% set what to plot as a function
plotLOGICAL_bin = [1, 0, 1];    % save all data
plotLOGICAL_mv = [1, 0, 1];
plotLOGICAL_all = [1, 0, 1];
% label figures by ROI and what they are
for iROI = 1:length(ROInames)
    for iStim = 1:length(nStim)
        saveName = ['glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim))];
        if nStim(iStim) == 8
            stimFreqs = stimData.stimFreqs_bin;
            plotROIav = 0;
            plotROIpTW = 1;
        elseif nStim(iStim) == 32
            stimFreqs = stimData.stimFreqs_mv;
            plotROIav = 1;
            plotROIpTW = 0;
        end
        eval([ROInames{iROI} '.roiAnalysis_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '= roiSplitGLMDataAnalysis(' ROInames{iROI} ',saveName,conditionSplitIndex,[],stimFreqs,plotROIav,plotROIpTW);'])
    end
end
% loop over ROIs