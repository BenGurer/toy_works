function data = plot_compareConditions_pCF(ConA,ConB,conMax)

figure('color',[1 1 1])

scatter(ConA,ConB)
xlim([0 conMax])
ylim([0 conMax])
hold on
f = 0:conMax;
fit = polyfit(ConA,ConB,1);
plot(polyval(fit,f));
correlation = corrcoef([ConA' ConB']);
text(conMax.*0.75,conMax.*0.25,sprintf('Correlation = %.2f \n y = %.2fx + %.2f',correlation(2), fit(1), fit(2)))
plot(f,f,'k--')
xlabel('Condition A - Estimated pCF')
ylabel('Condition B - Estimated pCF')
title('Voxel pCF estimates')
legend('Voxel pCF estimate','Line of best fit','Unity line','Location','best')

figure('color',[1 1 1])
histogram(ConA,conMax)
hold on
histogram(ConB,conMax)
xlim([0 conMax])
xlabel('Estimated pCF')
ylabel('Frequency')
title('Distribution Voxel pCF estimates')
legend('Condition A','Condition B')

end