function plotEstimatedBetas(filenames)
% filenames={'SparseROIbetas.mat','ContROIbetas.mat'};
filenames={'SparseUnbinnedROIbetas.mat','ContUnbinnedROIbetas.mat'};
% betas= zeros(length(e.betas),length(filenames));
% Ste= zeros(length(e.betas),length(filenames));
% normalisedBetas= zeros(length(e.betas),length(filenames));
% normalisedSte= zeros(length(e.betas),length(filenames));
betas= [];
Ste= [];
normalisedBetas= [];
normalisedSte= [];


figure;
for i = 1:length(filenames)
load(filenames{i})
betas(:,i)= e.betas;
Ste(:,i) = e.betaSte(:,:,1);
normalisedBetas(:,i) = betas(:,i)./(mean(betas(:,i)));
normalisedSte(:,i) = Ste(:,i)./(mean(betas(:,i)));


subplot(length(filenames),2,i)
errorbar(1:length(betas(:,i)),betas(:,i),Ste(:,i));
title([filenames{i} '_mean Ste = ' num2str(mean(Ste(:,i)))])
% ylim ([-0.5 3.5])
ylim ([-2 6])
xlim ([0 length(normalisedBetas(:,i))+1])
hold on
plot([0 length(betas(:,i))+1],[0 0],'--k')

subplot(length(filenames),2,i+length(filenames))
errorbar(1:length(normalisedBetas(:,i)),normalisedBetas(:,i),normalisedSte(:,i));
title(['mean Ste = ' num2str(mean(normalisedSte(:,i)))])
% ylim ([-1 2.5])
ylim ([-1.5 3])
xlim ([0 length(normalisedBetas(:,i))+1])
hold on
plot([0 length(betas(:,i))+1],[0 0],'--k')
end
