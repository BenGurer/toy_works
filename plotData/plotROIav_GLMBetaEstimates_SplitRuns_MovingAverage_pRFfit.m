function [ROIbetas, ROISte] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage_pRFfit(e1,e2,restrictIndex)

betasAin = squeeze(e1.betas);
betasBin = squeeze(e2.betas);
betasA = betasAin(:,restrictIndex);
betasB = betasBin(:,restrictIndex);

betaSteAin = squeeze(e1.betaSte);
betaSteA = betaSteAin(:,restrictIndex);
betaSteBin = squeeze(e2.betaSte);
betaSteB = betaSteBin(:,restrictIndex);

nBins = 8;
groupSize = size(betasA,1)/nBins;
% loopLength = (length(d.stimNames))/groupSize;
% loopLength = size(betasA,1) - groupSize;
loopLength = size(betasA,1) - groupSize;
betasScanA = zeros(loopLength,size(betasA,2));
betasScanB = zeros(loopLength,size(betasB,2));
betaSteScanA = zeros(loopLength,size(betaSteA,2));
betaSteScanB = zeros(loopLength,size(betaSteB,2));
% 
% [row, col] = find(isnan(ROIbetaSplitSumD));

for i = 1:loopLength
betasScanA(i,:) = mean(betasA(i:i+groupSize-1,:));
betasScanB(i,:) = mean(betasB(i:i+groupSize-1,:));

betaSteScanA(i,:) = mean(betaSteA(i:i+groupSize-1,:));
betaSteScanB(i,:) = mean(betaSteB(i:i+groupSize-1,:));
end

ROIbetaSplitSumC = zeros(size(betasScanA,1));
ROISteSplitSumC = zeros(size(betasScanA,1));
ROIbetaSplitCountC = zeros(size(betasScanA,1),1);
ROIbetaSplitSumD = zeros(size(betasScanA,1));
ROISteSplitSumD = zeros(size(betasScanA,1));
ROIbetaSplitCountD = zeros(size(betasScanA,1),1);

ROIbetaSumA = zeros(size(betasScanA,1));
ROISteSumA = zeros(size(betasScanA,1));
ROIbetaCountA = zeros(size(betasScanA,1),1);
ROIbetaSumB = zeros(size(betasScanA,1));
ROISteSumB = zeros(size(betasScanA,1));
ROIbetaCountB = zeros(size(betasScanA,1),1);

for i = 1:size(betasScanA,2)
    
    [VoxelMax VoxelIndex] = max(betasScanA(:,i));
    ROIbetaSplitSumC(VoxelIndex,:) = ROIbetaSplitSumC(VoxelIndex,:) + betasScanB(:,i)';
    ROISteSplitSumC(VoxelIndex,:) = ROISteSplitSumC(VoxelIndex,:) + betaSteScanB(:,i)';
    ROIbetaSplitCountC(VoxelIndex) = ROIbetaSplitCountC(VoxelIndex) + 1;
    
    [VoxelMax VoxelIndex] = max(betasScanB(:,i));
    ROIbetaSplitSumD(VoxelIndex,:) = ROIbetaSplitSumD(VoxelIndex,:) + betasScanA(:,i)';
    ROISteSplitSumD(VoxelIndex,:) = ROISteSplitSumD(VoxelIndex,:) + betaSteScanA(:,i)';
    ROIbetaSplitCountD(VoxelIndex) = ROIbetaSplitCountD(VoxelIndex) + 1;
    
    [VoxelMax VoxelIndex] = max(betasScanA(:,i));
    ROIbetaSumA(VoxelIndex,:) = ROIbetaSumA(VoxelIndex,:) + betasScanA(:,i)';
    ROISteSumA(VoxelIndex,:) = ROISteSumA(VoxelIndex,:) + betaSteScanA(:,i)';
    ROIbetaCountA(VoxelIndex) = ROIbetaCountA(VoxelIndex) + 1;
    
    [VoxelMax VoxelIndex] = max(betasScanB(:,i));
    ROIbetaSumB(VoxelIndex,:) = ROIbetaSumB(VoxelIndex,:) + betasScanB(:,i)';
    ROISteSumB(VoxelIndex,:) = ROISteSumB(VoxelIndex,:) + betaSteScanB(:,i)';
    ROIbetaCountB(VoxelIndex) = ROIbetaCountB(VoxelIndex) + 1;
    
end
ROIbetasC = ROIbetaSplitSumC./repmat(ROIbetaSplitCountC,1,size(betasScanA,1));
ROISteC = ROISteSplitSumC./repmat(ROIbetaSplitCountC,1,size(betasScanA,1));
ROIbetasC(isnan(ROIbetasC)) = 0;
ROISteC(isnan(ROISteC)) = 0;
ROISteAvC = mean(ROISteC,2);

ROIbetasD = ROIbetaSplitSumD./repmat(ROIbetaSplitCountD,1,size(betasScanA,1));
ROISteD = ROISteSplitSumD./repmat(ROIbetaSplitCountD,1,size(betasScanA,1));
ROIbetasD(isnan(ROIbetasD)) = 0;
ROISteD(isnan(ROISteD)) = 0;
ROISteAvD = mean(ROISteD,2);

ROIbetasA = ROIbetaSumA./repmat(ROIbetaCountA,1,size(betasScanA,1));
ROISteA = ROISteSumA./repmat(ROIbetaCountA,1,size(betasScanA,1));
ROIbetasA(isnan(ROIbetasA)) = 0;
ROISteA(isnan(ROISteA)) = 0;
ROISteAvA = mean(ROISteA,2);

ROIbetasB = ROIbetaSumB./repmat(ROIbetaCountB,1,size(betasScanA,1));
ROISteB = ROISteSumB./repmat(ROIbetaCountB,1,size(betasScanA,1));
ROIbetasB(isnan(ROIbetasB)) = 0;
ROISteB(isnan(ROISteB)) = 0;
ROISteAvB = mean(ROISteB,2);

ROIbetas = (ROIbetasC + ROIbetasD) ./2;
ROISte = (ROISteC + ROISteD) ./2;
ROISteAv = (ROISteAvC + ROISteAvD) ./2;
ROIbetaCountSplit = ROIbetaSplitCountC + ROIbetaSplitCountD;

ROIbetasNorm = ROIbetas ./ repmat(max(ROIbetas')',1,size(ROIbetas,1));
% ROIbetasNorm = ROIbetas ./ repmat(mean(ROIbetas')',1,size(ROIbetas,1));

ROIbetasUncorrected = (ROIbetasA + ROIbetasB) ./2;
ROISteUncorrected = (ROISteA + ROISteB) ./2;
ROISteAvUncorrected = (ROISteAvA + ROISteAvB) ./2;
ROIbetaCountUncorrected= ROIbetaCountA + ROIbetaCountB;

%% make legend lables
ROIbetaslabel = cell(size(ROIbetas,1),1);
for i = 1:size(ROIbetas,1)
    % ROIbetaslabel{i} = ['Condition = ' num2str(i) ' AvSte = ' num2str(ROISteAv(i))];
    
    ROIbetaslabel{i} = ['Condition = ' num2str(i)];
end

%% Condition Estimates - Split data
figure
subplot(2,2,1)
for i = 1:size(ROIbetas,1)
    % errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',2)
    plot(ROIbetas(i,:),'LineWidth',2)
    set(gca,'ColorOrder',jet(size(ROIbetas,1)))
    ylim(abs([min(min(ROIbetas)) max(max(ROIbetas))]));
    xlabel('Condition ID')
    ylabel('Mean Response Level')
    hold on
end
plot([1 size(ROIbetas,1)],[0 0],'--k')
legend(ROIbetaslabel,'Location','southeastoutside')

%% normalised Condition Estimates - Split data
subplot(2,2,2)
for i = 1:size(ROIbetas,1)
    % errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',2)
    plot(ROIbetasNorm(i,:),'LineWidth',2)
    set(gca,'ColorOrder',jet(size(ROIbetas,1)))
    xlim([1 size(ROIbetas,1)]);
    ylim([min(min(ROIbetasNorm)) max(max(ROIbetasNorm))]);
    xlabel('Condition ID')
    ylabel('Mean Response Level')
    hold on
end
plot([1 size(ROIbetas,1)],[0 0],'--k')
% legend(ROIbetaslabel)

%% plot number of voxels responsing max to each Condition
subplot(2,2,3)
% plot(ROIbetaCountSplit,'LineWidth',2)
bar(ROIbetaCountSplit)
hold on
% plot(ROIbetaCountUncorrected,'--','LineWidth',2)
title('Number of Voxels with Max Response to Condition')
xlabel('Condition ID')
ylabel('Voxel Max Count')
% legend('Split Data','Uncorrected')

%% Condition Estimates - uncorrected
subplot(2,2,4)
ROIbetaslabelUncorrected = cell(size(ROIbetas,1),1);
for i = 1:size(ROIbetas,1)
    ROIbetaslabelUncorrected{i} = ['Condition = ' num2str(i) ' AvSte = ' num2str(ROISteAvUncorrected(i))];
end
for i = 1:size(ROIbetasUncorrected,1)
    % errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',2)
    plot(ROIbetasUncorrected(i,:),'LineWidth',2)
    set(gca,'ColorOrder',jet(size(ROIbetasUncorrected,1)))
    xlabel('Condition ID')
    ylabel('Mean Response Level')
    hold on
end
plot([1 size(ROIbetas,1)],[0 0],'--k')
% legend(ROIbetaslabelUncorrected)

%% plot Condition tuning width estimates
figure
subIndex = [size(ROIbetas,1)/(size(ROIbetas,1)/4) size(ROIbetas,1)/4];
fpRF = [];
rfModel = [];
for i = 1:size(ROIbetas,1)
    subplot(subIndex(1),subIndex(2),i)
%     errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',1)
    plot(ROIbetas(i,:),'b','LineWidth',2)
    hold on
    plot(ROIbetasC(i,:),'g--','LineWidth',1)
    plot(ROIbetasD(i,:),'c--','LineWidth',1)
    plot(ROIbetasUncorrected(i,:),'m:','LineWidth',1)
    
    % plot(ROIbetaSum(i,:),'LineWidth',2)
    set(gca,'ColorOrder',jet(size(ROIbetas,1)))
    xlim([1 size(ROIbetas,1)]);
    ylim([min(min(ROIbetas))-min(min(ROISte)) max(max(ROIbetas))+max(max(ROISte))]);
    xlabel(['Condition ID = ' num2str(i)])
    ylabel('Sum Response Level')
    
    plot([1 size(ROIbetas,1)],[0 0],'--k')
    
    % fit pRF
   fun = @(x,xdata) x(4) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(5);
% fun = @(x,xdata) x(4) .* (exp(-(((xdata-x(1)).^2)/(2*(x(2)^2))))).^x(3) + x(5);
% fun = @(x,xdata) max(x(4) .* 1 + x(3).*log10((exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))))) + x(5),0);
xdata = 1:length(ROIbetas);
x0 = [i,1,2,1,1];
[fpRF(i,:),resnorm(i,:),residual(i,:)] = lsqcurvefit(fun,x0,xdata,ROIbetas(i,:));
rfModel{i} = fun(fpRF(i,:),xdata);
plot(rfModel{i},'r--','LineWidth',2)
end
legend('Split Mean','Run A', 'Run b', 'Mean','0','pRF model')

% use average to label - google f-test to see how to do it (degrees of freedom)
% use RSS to compare models
% RSS = resnorm
% n = data points
% p = parameters
% F = (RSSa - RRSb/p2-p1)/(RSSb/n-p2)

figure
xpRFplot = -32:32;
pRFshift = [];
for i = 1:size(ROIbetas,1)
pRFshift(i,:) = fun([0 fpRF(i,2:end)],xpRFplot);
plot(xpRFplot,pRFshift(i,:));
axis tight
hold on
end
plot(xpRFplot,mean(pRFshift),'b','LineWidth',2)
% plot(mean(cell2mat(rfModel')),'b','LineWidth',2)

figure
f = 1:size(betasScanA,1);
for i = 1:size(betasScanA,2)    
    [VoxelMax VoxelIndexA(i)] = max(betasScanA(:,i));
     wA = betasScanA(:,i)';
     wA = (wA-min(wA))./(max(wA)-min(wA));
     VoxelCntrdA(i) = sum(wA.*f)./sum(wA);
     VoxelSprdA(i) = sqrt(sum(wA.*(f-VoxelCntrdA(i)).^2)/sum(wA));
    [VoxelMax VoxelIndexB(i)] = max(betasScanB(:,i));
     wB = betasScanB(:,i)';
     wB = (wB-min(wB))./(max(wB)-min(wB));
     VoxelCntrdB(i) = sum(wB.*f)./(sum(wB));
     VoxelSprdB(i) = sqrt(sum(wB.*(f-VoxelCntrdB(i)).^2)/sum(wB));
     
% %           VoxelCntrdC(i) = sum(f.*wA)/(sum(wA));
%      VoxelSprC(i) = sqrt(sum(wB.*(f-VoxelCntrdA(i)).^2)/sum(wB));
%      VoxelSprD(i) = sqrt(sum(wA.*(f-VoxelCntrdB(i)).^2)/sum(wA));
     VoxelCntrd(i) = mean([VoxelCntrdA(i) VoxelCntrdB(i)]);
     VoxelSprd(i) = mean([VoxelSprdA(i) VoxelSprdB(i)]);
end
% betaEstAvB = hist3([VoxelIndexA',VoxelIndexB'],[8 8]);
% pcolor(betaEstAvB)
% % surf(n)
% colorbar
% xlabel('Condition Estimate Run A')
% ylabel('Condition Estimate Run B')
% 
% figure
subplot(2,2,1)
scatter(VoxelCntrdA,VoxelCntrdB)
title('Centriod')
hold on
fit = polyfit(VoxelCntrdA,VoxelCntrdB,1);
xlim([1 size(betasScanA,1)]);
ylim([1 size(betasScanA,1)]);
plot(polyval(fit,f));
correlation = corrcoef([VoxelCntrdA' VoxelCntrdB']);
text(max(xlim),max(ylim),sprintf('correlation = %.2f \n fit = %.2f %.2f',correlation(2), fit(1), fit(2)))
plot(f,f,'--')


% subplot(2,2,2)
% scatter(VoxelSprdC,VoxelSprdD)
% title('Spread: Split Est')
% hold on
% xlim([1 4]);
% ylim([1 4]);
% fit = polyfit(VoxelSprC,VoxelSprD,1);
% plot(polyval(fit,f));
% correlation = corrcoef([VoxelSprC' VoxelSprD']);
% text(max(xlim),max(ylim),sprintf('correlation = %.2f \n fit = %.2f %.2f',correlation(2), fit(1), fit(2)))

subplot(2,2,3)
scatter(VoxelSprdA,VoxelSprdB)
title('Spread')
hold on
xlim auto
ylim auto
fit = polyfit(VoxelSprdA,VoxelSprdB,1);
plot(polyval(fit,f));
correlation = corrcoef([VoxelSprdA' VoxelSprdB']);
text(max(xlim),max(ylim),sprintf('correlation = %.2f \n fit = %.2f %.2f',correlation(2), fit(1), fit(2)))

subplot(2,2,4)
scatter(VoxelCntrd,VoxelSprd)
title('Centriod vs Spread')
hold on
xlim auto
ylim auto
fit = polyfit(VoxelCntrd,VoxelSprd,1);
plot(polyval(fit,f));

figure
c = 1;
% per group recentre on max and average
ROIbetasSum = zeros(nBins-1,(length(ROIbetas)+groupSize).*2);

subIndexSum = round([size(ROIbetasSum,1)/(size(ROIbetasSum,1)/2) size(ROIbetasSum,1)/2]);
for i = 1:nBins-1
    subplot(subIndexSum(1),subIndexSum(2),i)
%     plot(ROIbetas(c,:),'LineWidth',2)
% ROIbetasSum(i,:) = sum(ROIbetas(c:c + groupSize,:))./groupSize;
c = c + groupSize;
for ii =1:groupSize
    a = c-ii;    
    b = c-ii+length(ROIbetas(ii,:))-1;
%     b = c+groupSize-ii:c-1+groupSize-ii+length(ROIbetas(ii,:));
    ROIbetasSum(i,a:b) = ROIbetasSum(i,a:b) + ROIbetas(a,:);
end
ROIbetasSum(i,:) = ROIbetasSum(i,:)./groupSize;
    plot(ROIbetasSum(i,:),'LineWidth',2)
%     hold on
%     plot(ROIbetasC(c,:),'--','LineWidth',1)
%     plot(ROIbetasD(c,:),'--','LineWidth',1)
%     plot(ROIbetasUncorrected(c,:),':','LineWidth',1)
%     set(gca,'ColorOrder',jet(size(ROIbetas,1)))
%     xlim([1 size(ROIbetas,1)]);
%     ylim([min(min(ROIbetas))-min(min(ROISte)) max(max(ROIbetas))+max(max(ROISte))]);
%     xlabel(['Condition ID = ' num2str((c+groupSize-1)-(groupSize/2))])
%     ylabel('Sum Response Level')    
%     plot([1 size(ROIbetas,1)],[0 0],'--k')
%     c = c + groupSize;
end
legend('Split Mean','Run A', 'Run b', 'Mean')

% combine above and below
% figure out scale/xacis


figure
% c = 1;
for i = 1:nBins-1
fun = @(x,xdata) x(4) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(5);
% fun = @(x,xdata) x(4) .* (exp(-(((xdata-x(1)).^2)/(2*(x(2)^2))))).^x(3) + x(5);
% fun = @(x,xdata) max(x(4) .* 1 + x(3).*log10((exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))))) + x(5),0);
xdata = 1:length(ROIbetasSum);
y = ROIbetasSum(i,:);
[m index] = max(ROIbetasSum(i,:));

x0 = [index,1,2,1,1];
fpRF(i,:) = lsqcurvefit(fun,x0,xdata,y);
rfModel{i} = fun(fpRF(i,:),xdata);

subplot(subIndexSum(1),subIndexSum(2),i)
plot(xdata,y)
hold on
plot(xdata,rfModel{i})
% c = c+1;
end

figure
xpRFplot = -32:32;
for i = 1:nBins-1     
plot(xpRFplot,fun([0 fpRF(i,2:end)],xpRFplot));
axis tight
hold on
end

figure
c = 1;
for i = 1:nBins-1
% fun = @(x,xdata) x(4) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(5);
fun = @(x,xdata) x(4) .* (exp(-(((xdata-x(1)).^2)/(2*(x(2)^2))))).^x(3) + x(5);
% fun = @(x,xdata) max(x(4) .* 1 + x(3).*log10((exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))))) + x(5),0);
xdata = 1:length(ROIbetasSum);
y = ROIbetasSum(i,:);
[m index] = max(ROIbetasSum(i,:));

x0 = [index,1,2,1,1];
fpRF(i,:) = lsqcurvefit(fun,x0,xdata,y);
rfModel{i} = fun(fpRF(i,:),xdata);

subplot(subIndexSum(1),subIndexSum(2),i)
plot(xdata,y)
hold on
plot(xdata,rfModel{i})
c = c+1;
end

figure
xpRFplot = -32:32;
for i = 1:nBins-1     
plot(xpRFplot,fun([0 fpRF(i,2:end)],xpRFplot));
axis tight
hold on
end

figure
c = 1;
for i = 1:nBins-1
% fun = @(x,xdata) x(4) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(5);
% fun = @(x,xdata) x(4) .* (exp(-(((xdata-x(1)).^2)/(2*(x(2)^2))))).^x(3) + x(5);
fun = @(x,xdata) max(x(4) .* 1 + x(3).*log10((exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))))) + x(5),0);
xdata = 1:length(ROIbetasSum);
y = ROIbetasSum(i,:);
[m index] = max(ROIbetasSum(i,:));

x0 = [index,1,2,1,1];
fpRF(i,:) = lsqcurvefit(fun,x0,xdata,y);
rfModel{i} = fun(fpRF(i,:),xdata);

subplot(subIndexSum(1),subIndexSum(2),i)
plot(xdata,y)
hold on
plot(xdata,rfModel{i})
c = c+1;
end

figure
xpRFplot = -32:32;
for i = 1:nBins-1     
plot(xpRFplot,fun([0 fpRF(i,2:end)],xpRFplot));
axis tight
hold on
end
