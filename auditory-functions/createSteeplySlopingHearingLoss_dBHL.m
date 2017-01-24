function Threshold_dBHL = createSteeplySlopingHearingLoss_dBHL(f)

Slope = 36.5; % dB/oct 
F_HearingLoss_lower = 3;
F_HearingLoss_upper = 8;

Threshold_HL = max(Slope*log2(f/F_HearingLoss_lower),0);
Threshold_HL_upper = Threshold_HL(find(f>=F_HearingLoss_upper,1,'first'));

Threshold_dBHL = min(Threshold_HL,Threshold_HL_upper);

figure
plot(f,Threshold_dBHL)