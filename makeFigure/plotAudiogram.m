function plotAudiogram(tFrequencies,left,right)
left = [5 -5 5 5 5 0 0 -5 5 10 15];
right = [0 0 5 5 0 5 0 -5 5 10 10];
tFrequencies = [125, 250, 500,750, 1000, 1500, 2000, 3000, 4000, 6000, 8000];
f = tFrequencies/1000;

% figure('Color',[1 1 1]);
% plot(f,left)
% hold on
% plot(f,right)
% xlim([min(f) max(f)])
% xlabel('Frequency (kHz)')
% ylabel('dB HL')

figure('Color',[1 1 1]);
semilogx(f,left,'-<','MarkerSize',5,'LineWidth',2)
hold on
semilogx(f,right,'->','MarkerSize',5,'LineWidth',2)
% xlim([min(f) max(f)])
ylim([-10 70])
a = gca;
axis(a,'ij');
% ydirection('ij')
xlabel('Frequency (kHz)')
ylabel('dB HL')
xticks(f)

end