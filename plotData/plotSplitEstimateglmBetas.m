function plotSplitEstimateglmBetas

dataDir = 'N:\data\CorticalMagnification\03644_012\';

% sparseFiles = {'eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=1PAC.mat','eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=3PAC.mat'};
% contFiles = {'eGLM_GLM Double Gamma 8 bins ContMotionComp_scanNum=2PAC.mat','eGLM_GLM Double Gamma 8 bins ContMotionComp_scanNum=4PAC.mat'};
% %
% fileNameScanA ='eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=1PAC.mat';
% fileNameScanB ='eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=3PAC.mat';
% %
fileNameScanA ={'eGLM_GLM Double Gamma 8 bins ContMotionComp_scanNum=2PAC.mat','eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=1PAC.mat'};
fileNameScanB ={'eGLM_GLM Double Gamma 8 bins ContMotionComp_scanNum=4PAC.mat','eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=3PAC.mat'};
AcqType = {'Cont','Sparse'};
for n = 1:length(AcqType)

dScanA = load([dataDir fileNameScanA{n}]);
dScanB = load([dataDir fileNameScanB{n}]);

betasScanA = squeeze(dScanA.f.betas);
betasScanB = squeeze(dScanB.f.betas);

betaSteScanA = squeeze(dScanA.f.betaSte);
betaSteScanB = squeeze(dScanB.f.betaSte);

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
ROISteAvC = mean(ROISteC,2);

ROIbetasD = ROIbetaSplitSumD./repmat(ROIbetaSplitCountD,1,size(betasScanA,1));
ROISteD = ROISteSplitSumD./repmat(ROIbetaSplitCountD,1,size(betasScanA,1));
ROISteAvD = mean(ROISteD,2);



ROIbetasA = ROIbetaSumA./repmat(ROIbetaCountA,1,size(betasScanA,1));
ROISteA = ROISteSumA./repmat(ROIbetaCountA,1,size(betasScanA,1));
ROISteAvA = mean(ROISteA,2);

ROIbetasB = ROIbetaSumB./repmat(ROIbetaCountB,1,size(betasScanA,1));
ROISteB = ROISteSumB./repmat(ROIbetaCountB,1,size(betasScanA,1));
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
figure('name',AcqType{n});
subplot(2,2,1)
for i = 1:size(ROIbetas,1)
    % errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',2)
    plot(ROIbetas(i,:),'LineWidth',2)
    set(gca,'ColorOrder',jet(size(ROIbetas,1)))
    ylim([min(min(ROIbetas)) max(max(ROIbetas))]);
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
figure('name',AcqType{n});
subIndex = [size(ROIbetas,1)/(size(ROIbetas,1)/2) size(ROIbetas,1)/2];
for i = 1:size(ROIbetas,1)
    subplot(subIndex(1),subIndex(2),i)
%     errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',1)
    plot(ROIbetas(i,:),'LineWidth',2)
    hold on
    plot(ROIbetasC(i,:),'--','LineWidth',1)
    plot(ROIbetasD(i,:),'--','LineWidth',1)
    plot(ROIbetasUncorrected(i,:),':','LineWidth',1)
    
    % plot(ROIbetaSum(i,:),'LineWidth',2)
    set(gca,'ColorOrder',jet(size(ROIbetas,1)))
    xlim([1 size(ROIbetas,1)]);
    ylim([min(min(ROIbetas))-min(min(ROISte)) max(max(ROIbetas))+max(max(ROISte))]);
    xlabel(['Condition ID = ' num2str(i)])
    ylabel('Sum Response Level')
    
    plot([1 size(ROIbetas,1)],[0 0],'--k')
end
legend('Split Mean','Run A', 'Run b', 'Mean')


figure('name',AcqType{n});
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


subplot(2,2,2)
% scatter(VoxelSprC,VoxelSprD)
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
title('Spread: not split')
hold on
xlim([1 4]);
ylim([1 4]);
fit = polyfit(VoxelSprdA,VoxelSprdB,1);
plot(polyval(fit,f));
correlation = corrcoef([VoxelSprdA' VoxelSprdB']);
text(max(xlim),max(ylim),sprintf('correlation = %.2f \n fit = %.2f %.2f',correlation(2), fit(1), fit(2)))

subplot(2,2,4)
scatter(VoxelCntrd,VoxelSprd)
title('Centriod vs Spread')
hold on
xlim([1 8]);
ylim([1 4]);
fit = polyfit(VoxelCntrd,VoxelSprd,1);
plot(polyval(fit,f));

data{n}.ROIbetas = ROIbetas;
data{n}.VoxelbetaMax = mean([VoxelIndexA; VoxelIndexB]);
data{n}.VoxelCntrd = VoxelCntrd;
data{n}.VoxelSprd = VoxelSprd;
data{n}.VoxelCntrdA = VoxelCntrdA;
data{n}.VoxelCntrdB = VoxelCntrdB;
end

figure('name','Cont v Sparse: Beta TW');
lineStyleIndex = {'--',':'};
for n = 1:length(data)
    subIndex = [size(ROIbetas,1)/(size(ROIbetas,1)/2) size(ROIbetas,1)/2];
for i = 1:size(ROIbetas,1)
    subplot(subIndex(1),subIndex(2),i)
%     errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',1)
    plot(data{n}.ROIbetas(i,:),lineStyleIndex{n},'LineWidth',2)
    hold on
%     plot(ROIbetasC(i,:),'--','LineWidth',1)
%     plot(ROIbetasD(i,:),'--','LineWidth',1)
%     plot(ROIbetasUncorrected(i,:),':','LineWidth',1)
    
    % plot(ROIbetaSum(i,:),'LineWidth',2)
    set(gca,'ColorOrder',jet(size(ROIbetas,1)))
    xlim([1 size(ROIbetas,1)]);
    ylim([-2 4]);
    xlabel(['Condition ID = ' num2str(i)])
    ylabel('Sum Response Level')
    plot([1 size(ROIbetas,1)],[0 0],'--k')
end
legend('Cont','Sparse')

end
% for i = 1:size(data{1}.ROIbetas,1)    
%     [VoxelMax VoxelIndexCont(i)] = max(data{1}.ROIbetas(i,:));    
%     [VoxelMax VoxelIndexSparse(i)] = max(data{2}.ROIbetas(i,:)); 
% end
% betaE
figure('name','Cont v Sparse: Voxel Centriod');
scatter(data{1}.VoxelCntrd,data{2}.VoxelCntrd)
hold on
xlim([1 8]);
ylim([1 8]);
fit = polyfit(data{1}.VoxelCntrd,data{2}.VoxelCntrd,1);
plot(polyval(fit,f));
correlation = corrcoef([data{1}.VoxelCntrd' data{2}.VoxelCntrd']);
text(max(xlim),max(ylim),sprintf('correlation = %.2f',correlation(2)))
xlabel('Est Centriod: Cont')
ylabel('Est Centriod: Sparse')
plot(f,f,'--')


figure('name','Cont v Sparse: Voxel Spread');
scatter(data{1}.VoxelSprd,data{2}.VoxelSprd)
hold on
xlim([1 4]);
ylim([1 4]);
fit = polyfit(data{1}.VoxelSprd,data{2}.VoxelSprd,1);
plot(polyval(fit,f));
correlation = corrcoef([data{1}.VoxelSprd' data{2}.VoxelSprd']);
text(max(xlim),max(ylim),sprintf('correlation = %.2f',correlation(2)))
xlabel('Est Spread: Cont')
ylabel('Est Spread: Sparse')

% betaEstAvB = hist3([data{1}.VoxelbetaMax',data{2}.VoxelbetaMax'],[8 8]);
% scatter(data{1}.VoxelbetaMax,data{2}.VoxelbetaMax);
% pcolor(betaEstAvB) 
% colorbar
% xlabel('Condition Estimate Cont')
% ylabel('Condition Estimate Sparse')

