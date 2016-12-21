function plotSplitEstimatepRF
dataDir = 'N:\data\CorticalMagnification\03644_012\MotionComp\pRFAnal\';
labels = {'Sparse','Cont'};
figure
f = 1:45;
for i = 1:2
fileNameA = ['pRF_auditory_SCAN' num2str(i)];
fileNameB = ['pRF_auditory_SCAN' num2str(i+2)];

A = load([dataDir fileNameA  '.mat']);
B = load([dataDir fileNameB  '.mat']);

A = struct2cell(A);
B = struct2cell(B);
runParamsA = A{1}.d{1,i}.params;
runParamsB = B{1}.d{1,i+2}.params;

subplot(2,2,i)
scatter(runParamsA(1,:),runParamsB(1,:))

xlabel('Est pCF: Run A')
ylabel('Est pCF: Run B')
hold on
fit = polyfit(runParamsA(1,:),runParamsB(1,:),1);
plot(polyval(fit,f));
pCFr = corrcoef([runParamsA(1,:)',runParamsB(1,:)']);
title([labels{i} sprintf(' pCF\n correlation = %.2f \n fit = %.2f %.2f',pCFr(2), fit(1), fit(2))])
plot(f,f,'--')

subplot(2,2,i+2)
scatter(runParamsA(3,:),runParamsB(3,:))
xlabel('Est pTW: Run A')
ylabel('Est pTW: Run B')
hold on
fit = polyfit(runParamsA(3,:),runParamsB(3,:),1);
plot(polyval(fit,f));
pTWr = corrcoef([runParamsA(3,:)',runParamsB(3,:)']);
title([labels{i} sprintf(' pTW\n correlation = %.2f \n fit = %.2f %.2f',pTWr(2), fit(1), fit(2))])
plot(f,f,'--')

pCF{i} = mean([runParamsA(1,:)',runParamsB(1,:)'],2);
pTW{i} = mean([runParamsA(3,:)',runParamsB(3,:)'],2);
end

figure
scatter(pCF{1},pCF{2})
figure
scatter(pTW{1},pTW{2})