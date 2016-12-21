function glmSplitEstimate

dataDir = 'N:\data\CorticalMagnification\03644_012\';
% fileNameScanA = 'eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=1LeftPAC.mat';
% fileNameScanB = 'eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=3LeftPAC.mat';
% 
% fileNameScanA = 'eGLM_GLM Double Gamma 8 binsMotionComp_scanNum=2LeftPAC.mat';
% fileNameScanB = 'eGLM_GLM Double Gamma 8 binsMotionComp_scanNum=4LeftPAC.mat';


% fileNameScanA ='eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=1PAC.mat';
% fileNameScanB ='eGLM_GLM Double Gamma 8 bins SparseMotionComp_scanNum=3PAC.mat';

% fileNameScanA ='eGLM_GLM Double Gamma 8 binsMotionComp_scanNum=2PAC.mat';
% fileNameScanB ='eGLM_GLM Double Gamma 8 binsMotionComp_scanNum=4PAC.mat';

fileNameScanA ='eGLM_GLM Double Gamma SparseMotionComp_scanNum=1PAC.mat';
fileNameScanB ='eGLM_GLM Double Gamma SparseMotionComp_scanNum=3PAC.mat';

% fileNameScanA ='eGLM_GLM Double Gamma ContMotionComp_scanNum=2PAC.mat';
% fileNameScanB ='eGLM_GLM Double Gamma ContMotionComp_scanNum=4PAC.mat';

dScanA = load([dataDir fileNameScanA]);
dScanB = load([dataDir fileNameScanB]);

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

ROIbetasUncorrected = (ROIbetasA + ROIbetasB) ./2;
ROISteUncorrected = (ROISteA + ROISteB) ./2;
ROISteAvUncorrected = (ROISteAvA + ROISteAvB) ./2;
ROIbetaCountUncorrected= ROIbetaCountA + ROIbetaCountB;

%% make legend lables
ROIbetaslabel = cell(size(ROIbetas,1),1);
for i = 1:size(ROIbetas,1)
% ROIbetaslabel{i} = ['beta = ' num2str(i) ' AvSte = ' num2str(ROISteAv(i))];

ROIbetaslabel{i} = ['beta = ' num2str(i)];
end

%% Beta Estimates - Split data
figure('name',fileNameScanA);
subplot(2,2,1)
for i = 1:size(ROIbetas,1)
% errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',2)
plot(ROIbetas(i,:),'LineWidth',2)
set(gca,'ColorOrder',jet(size(ROIbetas,1)))
ylim([min(min(ROIbetas)) max(max(ROIbetas))]);
xlabel('Beta ID')
ylabel('Mean Response Level')
hold on
end
plot([1 size(ROIbetas,1)],[0 0],'--k')
legend(ROIbetaslabel,'Location','southeastoutside')

%% normalised Beta Estimates - Split data
subplot(2,2,2)
for i = 1:size(ROIbetas,1)
% errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',2)
plot(ROIbetasNorm(i,:),'LineWidth',2)
set(gca,'ColorOrder',jet(size(ROIbetas,1)))
xlim([1 size(ROIbetas,1)]);
ylim([min(min(ROIbetasNorm)) max(max(ROIbetasNorm))]);
xlabel('Beta ID')
ylabel('Mean Response Level')
hold on
end
plot([1 size(ROIbetas,1)],[0 0],'--k')
% legend(ROIbetaslabel)

%% plot number of voxels responsing max to each beta
subplot(2,2,3)
plot(ROIbetaCountSplit,'LineWidth',2)
hold on
plot(ROIbetaCountUncorrected,'--','LineWidth',2)
title('Number of Voxels with Max Response to Beta')
xlabel('Beta ID')
ylabel('Voxel Max Count')
legend('Split Data','Uncorrected')

%% Beta Estimates - uncorrected
subplot(2,2,4)
ROIbetaslabelUncorrected = cell(size(ROIbetas,1),1);
for i = 1:size(ROIbetas,1)
ROIbetaslabelUncorrected{i} = ['beta = ' num2str(i) ' AvSte = ' num2str(ROISteAvUncorrected(i))];
end
for i = 1:size(ROIbetasUncorrected,1)
% errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',2)
plot(ROIbetasUncorrected(i,:),'LineWidth',2)
set(gca,'ColorOrder',jet(size(ROIbetasUncorrected,1)))
xlabel('Beta ID')
ylabel('Mean Response Level')
hold on
end
plot([1 size(ROIbetas,1)],[0 0],'--k')
% legend(ROIbetaslabelUncorrected)

%% plot beta tuning width estimates
figure('name',fileNameScanA);
subIndex = [size(ROIbetas,1)/(size(ROIbetas,1)/2) size(ROIbetas,1)/2];
for i = 1:size(ROIbetas,1)
subplot(subIndex(1),subIndex(2),i)
errorbar(ROIbetas(i,:),ROISte(i,:),'LineWidth',1)
% plot(ROIbetaSum(i,:),'LineWidth',2)
set(gca,'ColorOrder',jet(size(ROIbetas,1)))
xlim([1 size(ROIbetas,1)]);
ylim([min(min(ROIbetas))-min(min(ROISte)) max(max(ROIbetas))+max(max(ROISte))]);
xlabel(['Beta ID = ' num2str(i)])
ylabel('Sum Response Level')
hold on
plot([1 size(ROIbetas,1)],[0 0],'--k')
end

% legend(label)