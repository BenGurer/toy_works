function [condition_splitMean, voxel_Mean_Peak, Voxel_Cntrd, Voxel_Sprd] = cal_splitMean(betas_A,betas_B)

% function [condition_splitMean, voxel_MeanPeak, voxel_Mean, condition_splitMean_max, condition_splitMean_peak, VoxelIndex_A, VoxelIndex_B, VoxelMax_A, VoxelMax_B] = cal_splitMean(betas_A,betas_B)

% add error bars?
% variation in betas estiamted and indiviusdal voxel beta estimation

% output
% ROI centriod and spread, pRF properties - fit to bins
% voxel centriod and spread, pRF properties - fit to voxels


[VoxelMax_A, VoxelIndex_A] = max(betas_A);
[VoxelMax_B, VoxelIndex_B] = max(betas_B);

for i = 1:size(betas_A,1)
    splitA(i,:) = sum(betas_A(:,VoxelIndex_B==i),2)/sum(VoxelIndex_B==i);
    splitB(i,:) = sum(betas_B(:,VoxelIndex_A==i),2)/sum(VoxelIndex_A==i);
end

condition_splitMean = (splitA + splitB) / 2;
% [condition_splitMean_max, condition_splitMean_peak] = max(condition_splitMean);

voxel_Mean_Index = (VoxelIndex_A + VoxelIndex_B) / 2;
voxel_Mean = (betas_A + betas_B) / 2;
[voxel_Mean_Max voxel_Mean_Peak] = max(voxel_Mean); 
% 
% roi_av = mean([mean(betas_A,2),mean(betas_A,2)],2);

% calculate weighted mean and spread
f = repmat((1:size(voxel_Mean,1))',1,size(voxel_Mean,2));
voxel_wA = (voxel_Mean-min(voxel_Mean))./(max(voxel_Mean)-min(voxel_Mean));
Voxel_Cntrd = sum(voxel_wA.*f)./sum(voxel_wA);
Voxel_Sprd = sqrt(sum(voxel_wA.*(f-Voxel_Cntrd).^2)/sum(voxel_wA));


% fit gaussian pRF

% output params of both

% voxel by voxel and condition by condition
% weight mean
%      wA = betasScanA(:,i)';
%      wA = (wA-min(wA))./(max(wA)-min(wA));
%      VoxelCntrdA(i) = sum(wA.*f)./sum(wA);
%      VoxelSprdA(i) = sqrt(sum(wA.*(f-VoxelCntrdA(i)).^2)/sum(wA));
%           wB = betasScanB(:,i)';
%      wB = (wB-min(wB))./(max(wB)-min(wB));
%      VoxelCntrdB(i) = sum(wB.*f)./(sum(wB));
%      VoxelSprdB(i) = sqrt(sum(wB.*(f-VoxelCntrdB(i)).^2)/sum(wB));
%      VoxelCntrd(i) = mean([VoxelCntrdA(i) VoxelCntrdB(i)]);
%      VoxelSprd(i) = mean([VoxelSprdA(i) VoxelSprdB(i)]);
% 
% % pRF
% for i = 1:nBins-1
% fun = @(x,xdata) x(4) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(5);
% % fun = @(x,xdata) x(4) .* (exp(-(((xdata-x(1)).^2)/(2*(x(2)^2))))).^x(3) + x(5);
% % fun = @(x,xdata) max(x(4) .* 1 + x(3).*log10((exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))))) + x(5),0);
% xdata = 1:length(ROIbetasSum);
% y = ROIbetasSum(i,:);
% [m index] = max(ROIbetasSum(i,:));
% 
% x0 = [index,1,2,1,1];
% fpRF(i,:) = lsqcurvefit(fun,x0,xdata,y);
% rfModel{i} = fun(fpRF(i,:),xdata);
% 
% subplot(subIndexSum(1),subIndexSum(2),i)
% plot(xdata,y)
% hold on
% plot(xdata,rfModel{i})
% % c = c+1;
% end
% 
% figure
% xpRFplot = -32:32;
% for i = 1:nBins-1     
% plot(xpRFplot,fun([0 fpRF(i,2:end)],xpRFplot));
% axis tight
% hold on
% end

end