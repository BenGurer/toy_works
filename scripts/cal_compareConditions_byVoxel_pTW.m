function [ConA, ConB] = cal_compareConditions_byVoxel_pTW(betas_A,betas_B)

[VoxelMax_A, VoxelIndex_A] = max(betas_A);

for i = 1:size(betas_A,1)
    ConA(i,:) = sum(betas_A(:,VoxelIndex_A==i),2)/sum(VoxelIndex_A==i);
    ConB(i,:) = sum(betas_B(:,VoxelIndex_A==i),2)/sum(VoxelIndex_A==i);
end

