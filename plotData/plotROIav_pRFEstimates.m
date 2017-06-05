function data = plotROIav_pRFEstimates(ROIpCF,ROIpTW,ROIscaling,ROIR2)
% figure
% add r2 - use it to weight mean
% weight mean by scaling?
% need to indicate number of voxels per bin too
maxStim = 45;
nBins = 9;
binBands = linspace(0,maxStim,nBins);
avROIpCF = zeros(1,nBins);
avROIpTW = zeros(1,nBins);
avROIscaling = zeros(1,nBins);
stdROIpCF = zeros(1,nBins);
stdROIpTW = zeros(1,nBins);
stdROIscaling = zeros(1,nBins);
avrfModel = cell(1,nBins);
x = 0:0.1:maxStim;
c = 0;
for i = 1:nBins-1
    index = ROIpCF>=binBands(i) & ROIpCF<binBands(i+1);
    binnedROIpCF{i} = ROIpCF(index);
    avROIpCF(i) = mean(ROIpCF(index));
    stdROIpCF(i) = std(ROIpCF(index));
    binnedROIpTW{i} = ROIpTW(index);
    avROIpTW(i) = mean(ROIpTW(index));
    stdROIpTW(i) = std(ROIpTW(index));
    binnedROIscaling{i} = ROIscaling(index);
    avROIscaling(i) = mean(ROIscaling(index));
    stdROIscaling(i) = std(ROIscaling(index));
    binnedROIR2{i} = ROIR2(index);
    avROIR2(i) = mean(ROIR2(index));
    stdROIR2(i) = std(ROIR2(index));
    
    avrfModel{i} = exp(-(((x-avROIpCF(i)).^2)/(2*(avROIpTW(i)^2))));
%     avrfModel{i} = avrfModel{i} * avROIscaling(i);
    for ii = 1:length(binnedROIpCF{i})
    binrfModel{i}{ii} = exp(-(((x-binnedROIpCF{i}(ii)).^2)/(2*(binnedROIpTW{i}(ii)^2))));
    
    pRFstruct.bin(c+ii) = i;
    pRFstruct.rfModel{c+ii} = exp(-(((x-binnedROIpCF{i}(ii)).^2)/(2*(binnedROIpTW{i}(ii)^2))));
    pRFstruct.CentredrfModel{c+ii} = exp(-(((x-mean(binnedROIpCF{i}(:))).^2)/(2*(binnedROIpTW{i}(ii)^2))));
    pRFstruct.pCF(c+ii) = binnedROIpCF{i}(ii);
    pRFstruct.pTW(c+ii) = binnedROIpTW{i}(ii);
    pRFstruct.R2(c+ii) = binnedROIR2{i}(ii);
        
    end
    c = c + length(binnedROIpCF{i});
%     plot(x,avrfModel{i})
%     hold on
    
    pRFavStruct.bin(i) = i;
    pRFavStruct.pCF(i) = mean(ROIpCF(index));   
    pRFavStruct.pCFstd(i) = std(ROIpCF(index));
    pRFavStruct.pTW(i) = mean(ROIpTW(index));
    pRFavStruct.pTWstd(i) = std(ROIpTW(index));
    pRFavStruct.rfModel{i} = exp(-(((x-pRFavStruct.pCF(i)).^2)/(2*(pRFavStruct.pTW(i)^2))));
    pRFavStruct.scaling(i) = mean(ROIscaling(index));
    pRFavStruct.scalingstd(i) = std(ROIscaling(index));
    pRFavStruct.hist(i) = sum(index);

end

% figure
% g = gramm('x',x,'y',pRFavStruct.rfModel);
% g.facet_grid([],pRFavStruct.bin)
% g.geom_line()
% g.draw()
% 
% figure
% g_1 = gramm('x',x,'y',pRFstruct.rfModel);
% g_1.facet_grid([],pRFstruct.bin)
% g_1.stat_summary()
% g_1.draw()
% 
% figure
% g_2 = gramm('x',x,'y',pRFstruct.CentredrfModel);
% g_2.facet_grid([],pRFstruct.bin)
% g_2.stat_summary()
% g_2.draw()

figure
g_3 = gramm('x',x,'y',pRFstruct.CentredrfModel,'color',pRFstruct.bin);
g_3.stat_summary()
g_3.draw()