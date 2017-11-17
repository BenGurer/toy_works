function pRFdata = cheap_get_pRFoverlayData()
%% NH vs HL
thisView = getMLRView;

roiAC = viewGet(thisView,'roi','RightACVol');
% roiAC = viewGet(thisView,'roi','RIGHT');
roiAC.scanCoords = getROICoordinates(thisView,roiAC);
concatenationGroup = {'ConcatenationNH','ConcatenationHLsim'};
analysisSaveName{1} = {'pRF_auditory'};
analysisSaveName{2} = {'pRF_auditory_notCorrected','pRF_auditory_SLweighting','pRF_auditory_BOLDweighting','pRF_auditory_fitted'};

% define overlay numbers
pCFoverlayNum = 2;
pTWoverlayNum = 3;
%get ROI estimates
for iGroup = 1:length(concatenationGroup)
for i = 1:length(analysisSaveName{iGroup})  
thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
analysisData{i} = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisSaveName{iGroup}{i}));
r2data = analysisData{i}.overlays(1).data{1};
volumeIndices = sub2ind(size(r2data),roiAC.scanCoords(1,:),roiAC.scanCoords(2,:),roiAC.scanCoords(3,:));
r2{i} = r2data(volumeIndices);
pCFdata = analysisData{i}.overlays(pCFoverlayNum).data{1};
pCFest{i} = pCFdata(volumeIndices);
pTWdata = analysisData{i}.overlays(pTWoverlayNum).data{1};
pTWest{i} = pTWdata(volumeIndices);
d{i} = analysisData{i}.d;

pRFdata(iGroup).d{i} = analysisData{i}.d;
pRFdata(iGroup).pCFest{i} = pCFdata(volumeIndices);
pRFdata(iGroup).pTWest{i} = pTWdata(volumeIndices);
pRFdata(iGroup).r2{i} = r2{i};
% plotROIav_pRFEstimates(pCFest{i},pTWest{i},analysisData{i}.d{1}.scale(1,:),r2{i});
end
end

r2Index = pRFdata(1).r2{1} > 0.1;
totalVoxels = sum(r2Index);
% remove nans from condition A
% create figure to plot into
figure('color',[1 1 1])
% fit = zeros(length(analysisSaveName{2}),3);

fit = [];
xMax = max(pRFdata(1).pCFest{1});
xMin = min(pRFdata(1).pCFest{1});
for i = 1:length(analysisSaveName{2})
conA = [];
conB = [];
subplot(round(length(analysisSaveName{2})/2),2,i)
% remove voxels that are nans in condition A & restrict by r2 values from con A
conA = pRFdata(1).pCFest{1}(~isnan(pRFdata(1).pCFest{1}(r2Index)));
conB = pRFdata(2).pCFest{i}(~isnan(pRFdata(1).pCFest{1}(r2Index)));
% remove voxels that are nans in condition B
conA(isnan(conB)) = [];
conB(isnan(conB)) = [];
scatter(conA,conB)
hold on
fit(i,:) = polyfit(conA,conB,1);
plot(xMin:xMax,polyval(fit(i,:),xMin:xMax));
xlim([xMin xMax]);
ylim([xMin xMax]);
title(analysisSaveName{2}{i})

plot(xMin:xMax,xMin:xMax,'--k')

end


figure('color',[1 1 1])
scatter(pRFdata(1).pCFest{1},pRFdata(2).pCFest{2})

figure('color',[1 1 1])
scatter(pRFdata(1).pCFest{1},pRFdata(2).pCFest{3})

figure('color',[1 1 1])
scatter(pRFdata(1).pCFest{1},pRFdata(2).pCFest{4})

for i = 1:length(analysisSaveName{2})
figure
% VoxelStruct(i).conA = data(1).pCFest{1}(nhIndex & r2Index);
% VoxelStruct(i).conB = data(2).pCFest{i}(nhIndex & r2Index);

VoxelStruct(i).conA = pRFdata(1).pCFest{1};
VoxelStruct(i).conB = pRFdata(2).pCFest{i};
g(i) = gramm('x',VoxelStruct(i).conA,'y',VoxelStruct(i).conB);
% g.geom_point('alpha',0.05)
g(i).geom_point()
g(i).stat_glm();
 g(i).set_names('x','Condition A - Frequency (ERB)','y','Condition B - Frequency (ERB)')
% g(i).stat_bin();
g(i).draw()
end

end