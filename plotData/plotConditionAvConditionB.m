function [ROIbetas, ROISte] = plotConditionAvConditionB(e1,e2,restrictIndex)

betasAin = squeeze(e1.betas);
betasBin = squeeze(e2.betas);
betasA = betasAin(:,restrictIndex);
betasB = betasBin(:,restrictIndex);

betaSteAin = squeeze(e1.betaSte);
betaSteA = betaSteAin(:,restrictIndex);
betaSteBin = squeeze(e2.betaSte);
betaSteB = betaSteBin(:,restrictIndex);

% Compare activity for voxels between conditions
% ConA - categorise voxels
% Index conB by conA results
% Compare differences

% loop over voxels to compare conditions

% weighted mean

% compare estiamtes between cons

%% weighted mean
nBins = 8;
groupSize = size(betasA,1)/nBins;
loopLength = size(betasA,1) - groupSize;
betas_MA_A = zeros(loopLength,size(betasA,2));
betas_MA_B = zeros(loopLength,size(betasB,2));
betaSte_MA_A = zeros(loopLength,size(betaSteA,2));
betaSte_MA_B = zeros(loopLength,size(betaSteB,2));

for i = 1:loopLength
betas_MA_A(i,:) = mean(betasA(i:i+groupSize-1,:));
betas_MA_B(i,:) = mean(betasB(i:i+groupSize-1,:));

betaSte_MA_A(i,:) = mean(betaSteA(i:i+groupSize-1,:));
betaSte_MA_B(i,:) = mean(betaSteB(i:i+groupSize-1,:));
end

% loop over voxels
for i = 1:size(betas_MA_A,2)
    % new variables - mean and splitmean - size = betasScan
    betas_SplitMean(:,i) = 
    
    % save voxel index and use previous method
    
end


end