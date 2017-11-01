function ROIdata = script_getROIdata(thisView,analysisData,analysisBaseNames,ROInames,analysisScanNum,dataType)
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

switch dataType
    case 'GLM'
        for iROI = 1:length(ROInames)
            for iAnal = 1:length(analysisBaseNames)
                eval([['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iAnal}] '{analysisScanNum{iAnal}} = getROIdata_GLM(thisView, analysisData{iAnal}, ROIdata.' ROInames{iROI} '.roi,analysisScanNum{iAnal});'])
            end
        end
        
    case 'overlays'
        
        for iROI = 1:length(ROInames)
                for iAnal = 1:length(analysisBaseNames)
%                     data_flatROI = get_ROIdata(analysisData,ROI);
                    data2get = eval(['analysisData.scans.', analysisBaseNames{iAnal}, '.overlayData{' mat2str(analysisScanNum{iAnal}) '}.data;']);
                    eval(['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iAnal} '{' mat2str(analysisScanNum{iAnal}) '} = get_ROIdata(data2get,ROIdata.' ROInames{iROI} '.roi);']);
                end
        end

        
end