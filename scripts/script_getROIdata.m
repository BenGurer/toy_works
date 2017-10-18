function ROIdata = script_getROIdata(thisView,analysisData,analysisBaseNames,ROInames,analysisScanNum)
% ,stimInfo,plotInfo,conditionRunIndex)
% load in data
% get roi
% save data within roi to strucutre

ROIdata = struct;
q = char(39);
for iROI = 1:length(ROInames)
    eval(['ROIdata.' ROInames{iROI} ' = struct;']);
    eval(['ROIdata.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
end

% get data from ROIs
for iROI = 1:length(ROInames)
    for iAnal = 1:length(analysisBaseNames)
        eval([['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iAnal}] '{analysisScanNum{iAnal}} = getROIdata_GLM(thisView, analysisData{iAnal}, ROIdata.' ROInames{iROI} '.roi,analysisScanNum{iAnal});'])
    end
end