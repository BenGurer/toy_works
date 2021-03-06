function ROIdata = script_getROIdata(thisView,analysisData,analysisBaseNames,ROI,analysisScanNum,dataType)
% ,stimInfo,plotInfo,conditionRunIndex)
% load in data
% get roi
% save data within roi to strucutre


% ROIdata = struct;
% q = char(39);
% for iROI = 1:length(ROInames)
%     eval(['ROIdata.' ROInames{iROI} ' = struct;']);
%     eval(['ROIdata.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
% end

% get data from ROIs


        
        if iscell(analysisScanNum) % see if scan number is defined (i.e. if its individual scans

                for iAnal = 1:length(analysisBaseNames)
                    %                     data_flatROI = get_ROIdata(analysisData,ROI);
%                     data2get = eval(['analysisData.', analysisBaseNames{iAnal}, '.overlayData{' mat2str(analysisScanNum{iAnal}) '}.data;']);

                    data2get = eval(['analysisData.', analysisBaseNames{iAnal}, '.betas{' mat2str(analysisScanNum{iAnal}) '}.data;']);
                    eval(['ROIdata.' analysisBaseNames{iAnal} '{' mat2str(analysisScanNum{iAnal}) '} = get_ROIdata(data2get,ROI);']);
                end

        else
                for iAnal = 1:length(analysisBaseNames)
                    %                     data_flatROI = get_ROIdata(analysisData,ROI);
%                     analysisData.glm_hrfDoubleGamma.overlayData.data
%                     data2get = eval(['analysisData.', analysisBaseNames{iAnal}, '.overlayData.data;']);
                    data2get = analysisData.data;
                    ROIdata = get_ROIdata(data2get,ROI);
%                     eval(['ROIdata.'  analysisBaseNames{iAnal} ' = get_ROIdata(data2get,ROI);']);
                end

            
        end

% ROIdata = struct;
% q = char(39);
% for iROI = 1:length(ROInames)
%     eval(['ROIdata.' ROInames{iROI} ' = struct;']);
%     eval(['ROIdata.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
% end
% 
% % get data from ROIs
% 
% switch dataType
%     case 'GLM'
%         for iROI = 1:length(ROInames)
%             for iAnal = 1:length(analysisBaseNames)
%                 eval([['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iAnal}] '{analysisScanNum{iAnal}} = getROIdata_GLM(thisView, analysisData{iAnal}, ROIdata.' ROInames{iROI} '.roi,analysisScanNum{iAnal});'])
%             end
%         end
%         
%     case 'overlays'
%         
%         if iscell(analysisScanNum) % see if scan number is defined (i.e. if its individual scans
%             for iROI = 1:length(ROInames)
%                 for iAnal = 1:length(analysisBaseNames)
%                     %                     data_flatROI = get_ROIdata(analysisData,ROI);
%                     data2get = eval(['analysisData.scanData.', analysisBaseNames{iAnal}, '.overlayData{' mat2str(analysisScanNum{iAnal}) '}.data;']);
%                     eval(['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iAnal} '{' mat2str(analysisScanNum{iAnal}) '} = get_ROIdata(data2get,ROIdata.' ROInames{iROI} '.roi);']);
%                 end
%             end
%         else
%             for iROI = 1:length(ROInames)
%                 for iAnal = 1:length(analysisBaseNames)
%                     %                     data_flatROI = get_ROIdata(analysisData,ROI);
% %                     analysisData.glm_hrfDoubleGamma.overlayData.data
%                     data2get = eval(['analysisData.', analysisBaseNames{iAnal}, '.overlayData.data;']);
%                     eval(['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iAnal} ' = get_ROIdata(data2get,ROIdata.' ROInames{iROI} '.roi);']);
%                 end
%             end
%             
%         end
        
end