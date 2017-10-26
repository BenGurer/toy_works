function data = get_ROIdata(thisView,analysisData,ROI)
    %
    %   usage: get_ROIdata(thisView,analysisData,ROI,iScan)
    %      by: Ben Gurer
    %    date: 26/10/2017
    % purpose: Get overlay data within ROI
    %   input: mrTools view, analysis data (from overlays), ROI
    %
if iscell(analysisData)
% get roi and then use cords and convert to indices
indices = sub2ind(size(analysisData{1}),ROI.coords(1,:),ROI.coords(2,:),ROI.coords(3,:));
data = cell(1,length(analysisData));
for i = 1:length(analysisData)
data{i} = overlayData.data{i}(indices);
end

else
indices = sub2ind(size(analysisData.data),ROI.coords(1,:),ROI.coords(2,:),ROI.coords(3,:));   
data{1} = overlayData.data(indices);
end

end