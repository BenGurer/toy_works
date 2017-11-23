function pRFdata = cheap_get_pRFoverlayData()
%% NH vs HL
thisView = getMLRView;

roiAC = viewGet(thisView,'roi','RightAC');
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
conA = [];
conB = [];
for i = 1:length(analysisSaveName{2})

subplot(round(length(analysisSaveName{2})/2),2,i)
% remove voxels that are nans in condition A & restrict by r2 values from con A
conA = pRFdata(1).pCFest{1}(~isnan(pRFdata(1).pCFest{1}(r2Index)));
conB{i} = pRFdata(2).pCFest{i}(~isnan(pRFdata(1).pCFest{1}(r2Index)));
% remove voxels that are nans in condition B
conA(isnan(conB{i})) = [];
conB{i}(isnan(conB{i})) = [];
scatter(conA,conB{i})
hold on
fit(i,:) = polyfit(conA,conB{i},1);
plot(xMin:xMax,polyval(fit(i,:),xMin:xMax));
xlim([xMin xMax]);
ylim([xMin xMax]);
title(analysisSaveName{2}{i})

plot(xMin:xMax,xMin:xMax,'--k')

end

figure('color', [1 1 1])
histogram(conA)
hold on

for i = 1:length(analysisSaveName{2})
histogram(conB{i})
end


unpacked.pCF = [conA conB{1} conB{2} conB{3} conB{4}];
% unpacked.con = [repmat(1,1,length(data(1).pCFest{1})) repmat(2,1,length(data(1).pCFest{1})) repmat(3,1,length(data(1).pCFest{1}))];
unpacked.con = [repmat({'Condition A'},1,length(conA)) repmat({'pRF_auditory_notCorrected'},1,sum(r2Index)) repmat({'pRF_auditory_SLweighting'},1,sum(r2Index) repmat({'pRF_auditory_BOLDweighting'},1,sum(r2Index))];
'pRF_auditory_notCorrected','pRF_auditory_SLweighting','pRF_auditory_BOLDweighting','pRF_auditory_fitted'};

unpacked.pCF = conA;
unpacked.con = repmat({'Condition A'},1,length(conA));
for i = 1:length(analysisSaveName{2})
unpacked.pCF = [unpacked.pCF conB{i}];
unpacked.con = [unpacked.con repmat({analysisSaveName{2}{i}},1,length(conB{i}))];
end
    

figure
hist  = gramm('x',unpacked.pCF,'color',unpacked.con);
hist.stat_bin('nbins',8);
hist.set_names('x','Frequency (ERB)','y','Number of Voxels','color','Condition:')
% hist.geom_label()
hist.draw()

hist = [];
figure
hist  = gramm('x',unpacked.pCF,'color',unpacked.con);
hist.set_names('x','Frequency (ERB)','y','Frequency','color','Condition:')
hist.stat_density();
hist.draw()



% g(i) = gramm('x',VoxelStruct(i).conA,'y',VoxelStruct(i).conB);
% g.geom_point('alpha',0.05)

scatData.conA = repmat(conA,1,length(analysisSaveName{2}));
scatData.conB = [];
scatData.pRFmeth = [];
for i = 1:length(analysisSaveName{2})
scatData.conB = [scatData.conB conB{i}];
scatData.pRFmeth = [scatData.pRFmeth repmat({analysisSaveName{2}{i}},1,length(conB{i}))];
end

figure
g  = gramm('x',scatData.conA,'y',scatData.conB,'color',scatData.pRFmeth);
g.set_names('x','Condition B - Frequency (ERB)','y','Condition A - Frequency (ERB)')
g.geom_point()
g.stat_glm();
g.draw()

figure
g  = gramm('x',conA,'y',conB{1});
g.set_names('x','Condition A - Frequency (ERB)','y','Condition B - Frequency (ERB)')
g.geom_point()
g.stat_glm();
g.draw()

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