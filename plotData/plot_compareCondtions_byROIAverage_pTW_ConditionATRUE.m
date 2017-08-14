function data = plot_compareCondtions_byROIAverage_pTW_ConditionATRUE(splitMean_ConA,voxel_MeanPeak_ConA,runA,runB)
for i = 1:size(betas_A,1)
    splitA(i,:) = sum(betas_A(:,VoxelIndex_B==i),2)/sum(VoxelIndex_B==i);
    splitB(i,:) = sum(betas_B(:,VoxelIndex_A==i),2)/sum(VoxelIndex_A==i);
end

end